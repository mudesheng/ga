package constructcf

import (
	"bufio"
	"compress/gzip"
	"encoding/binary"
	"fmt"
	"io"
	"log"
	"os"
	"runtime"
	"runtime/pprof"
	"strconv"
	"strings"
	"time"

	"github.com/jwaldrip/odin/cli"
	"github.com/mudesheng/ga/bnt"
	"github.com/mudesheng/ga/cuckoofilter"
	"github.com/mudesheng/ga/utils"
)

const (
	AllState    = 1
	ScaffState  = 2
	GapState    = 3
	EndString   = "end"
	ReadSeqSize = 100
)

type LibInfo struct {
	Name          string // name of library
	NumberReads   int64  // the number of reads
	Diverse       uint8
	Paired        uint8
	AsmFlag       uint8 // denote which assembly phase used, note 1 used for all step of assembly pipeline, note 2 used for scaffold phase only, 3 used for filling gap only
	SeqProfile    uint8 // denote the data origin
	QualBenchmark uint8 // the benchmark of quality score, strength encourage use phred+33
	TotalBasesNum int64
	ReadLen       int
	InsertSize    int // paired read insert size
	InsertSD      int // Standard Deviation
	//	fnNum         int      // the number of files
	FnName []string // the files name slice
}

type CfgInfo struct {
	MaxRdLen int // maximum read length
	MinRdLen int // minimum read length
	//	n        int // the number of library
	Libs []LibInfo
}

type ReadBnt struct {
	Seq    []byte
	Length int // length of sequence
}

type ReadSeqBucket struct {
	ReadBuf [ReadSeqSize]ReadBnt
	Count   int
}

/*type ChanStruct struct {
	readseq chan []string
	state   chan int // note if goroutinue can return
} */

func (s1 ReadBnt) BiggerThan(s2 ReadBnt) bool {
	if s1.Length < s2.Length {
		return false
	} else if s1.Length > s2.Length {
		return true
	} else {
		for i := 0; i < len(s1.Seq); i++ {
			if s1.Seq[i] < s2.Seq[i] {
				return false
			} else if s1.Seq[i] > s2.Seq[i] {
				return true
			}
		}
	}

	return false
}

/*func (s1 ReadBnt) Equal(s2 ReadBnt) bool {
	if s1.Length != s2.Length {
		return false
	} else {
		for i := 0; i < len(s1.Seq); i++ {
			if s1.Seq[i] != s2.Seq[i] {
				return false
			}
		}
	}

	return true
} */

func ParseCfg(fn string) (cfgInfo CfgInfo, e error) {
	var inFile *os.File
	var err error
	if inFile, err = os.Open(fn); err != nil {
		log.Fatal(err)
	}
	var libInfo LibInfo
	defer inFile.Close()
	reader := bufio.NewReader(inFile)
	eof := false
	for !eof {
		var line string
		line, err = reader.ReadString('\n')
		if err == io.EOF {
			err = nil
			eof = true
		} else if err != nil {
			log.Fatal(err)
		}
		fields := strings.Fields(line)
		if len(fields) == 0 {
			continue
		}
		var v int
		switch fields[0] {
		case "[global_setting]":
		case "[LIB]":
			if libInfo.Name != "" {
				cfgInfo.Libs = append(cfgInfo.Libs, libInfo)
				var nli LibInfo
				libInfo = nli
			}
		case "max_rd_len":
			v, err = strconv.Atoi(fields[2])
			cfgInfo.MaxRdLen = v
		case "min_rd_len":
			v, err = strconv.Atoi(fields[2])
			cfgInfo.MinRdLen = v
		case "name":
			libInfo.Name = fields[2]
		case "avg_insert_len":
			v, err = strconv.Atoi(fields[2])
			libInfo.InsertSize = v
		case "insert_SD":
			v, err = strconv.Atoi(fields[2])
			libInfo.InsertSD = v
		case "diverse_rd_len":
			v, err = strconv.Atoi(fields[2])
			libInfo.Diverse = uint8(v)
		case "asm_flag":
			v, err = strconv.Atoi(fields[2])
			libInfo.AsmFlag = uint8(v)
		case "seq_profile":
			v, err = strconv.Atoi(fields[2])
			libInfo.SeqProfile = uint8(v)
		case "qual_benchmark":
			v, err = strconv.Atoi(fields[2])
			libInfo.QualBenchmark = uint8(v)
		case "f1", "f2":
			libInfo.FnName = append(libInfo.FnName, fields[2])
		default:
			if fields[0][0] != '#' && fields[0][0] != ';' {
				log.Fatalf("noknown line: %s", line)
			}
		}
		if err != nil {
			e = err
			return
		}
	}
	if libInfo.Name != "" {
		cfgInfo.Libs = append(cfgInfo.Libs, libInfo)
	}

	return
}

func ExtendReadBnt2Byte(rb ReadBnt) (extRB ReadBnt) {
	extRB.Length = rb.Length
	extRB.Seq = make([]byte, extRB.Length)
	var crb ReadBnt
	crb.Length = rb.Length
	crb.Seq = make([]byte, len(rb.Seq))
	copy(crb.Seq, rb.Seq)
	// fmt.Printf("[ExtendReadBnt2Byte] rb: %v, crb: %v\n", rb, crb)

	for i := crb.Length - 1; i >= 0; i-- {
		base := crb.Seq[i/bnt.NumBaseInByte] & bnt.BaseMask
		extRB.Seq[i] = base
		crb.Seq[i/bnt.NumBaseInByte] >>= bnt.NumBitsInBase
	}

	return extRB
}

func GetReadBntKmer(extRBnt ReadBnt, startPos int, kmerlen int) (rbk ReadBnt) {
	rbk.Length = kmerlen
	rbk.Seq = make([]byte, (rbk.Length+bnt.NumBaseInByte-1)/bnt.NumBaseInByte)

	for i := 0; i < kmerlen; i++ {
		rbk.Seq[i/bnt.NumBaseInByte] <<= bnt.NumBitsInBase
		rbk.Seq[i/bnt.NumBaseInByte] |= extRBnt.Seq[i+startPos]
	}

	return rbk
}

func ReverseComplet(ks ReadBnt) (rs ReadBnt) {
	var tmp ReadBnt
	tmp.Length = ks.Length
	tmp.Seq = make([]byte, len(ks.Seq))
	copy(tmp.Seq, ks.Seq)
	rs.Length = ks.Length
	rs.Seq = make([]byte, len(ks.Seq))
	for i := tmp.Length - 1; i >= 0; i-- {
		base := tmp.Seq[i/bnt.NumBaseInByte] & bnt.BaseMask
		tmp.Seq[i/bnt.NumBaseInByte] >>= bnt.NumBitsInBase
		rs.Seq[(rs.Length-i-1)/bnt.NumBaseInByte] <<= bnt.NumBitsInBase
		rs.Seq[(rs.Length-i-1)/bnt.NumBaseInByte] |= bnt.BntRev[base]
	}
	return rs
}

func ParaConstructCF(cf cuckoofilter.CuckooFilter, cs chan ReadSeqBucket, wc chan ReadSeqBucket) {

	var wrsb ReadSeqBucket
	for {
		rsb := <-cs
		if rsb.Count == 0 {
			wc <- wrsb
			wc <- rsb
			break
		} else {
			if rsb.Count < ReadSeqSize {
				fmt.Printf("rsb.ReadBuf length is :%d\n", len(rsb.ReadBuf))
			}
			for i := 0; i < rsb.Count; i++ {
				extRBnt := ExtendReadBnt2Byte(rsb.ReadBuf[i])
				for j := 0; j < extRBnt.Length-cf.Kmerlen+1; j++ {
					ks := GetReadBntKmer(extRBnt, j, cf.Kmerlen)
					rs := ReverseComplet(ks)
					if ks.BiggerThan(rs) {
						ks, rs = rs, ks
					}
					//fmt.Printf("ks: %v, rs: %v\n", ks, rs)
					count, suc := cf.Insert(ks.Seq)
					if suc == false {
						log.Fatal("Insert to the CuckooFilter false")
					}
					//fmt.Printf("retrun count : %d\n", count)
					//fmt.Printf("count set: %d\n", cf.GetCount(ks.Seq))
					if count == 1 {
						if wrsb.Count >= ReadSeqSize {
							wc <- wrsb
							var nrsb ReadSeqBucket
							wrsb = nrsb
						}
						//fmt.Printf("[paraConstructCF] wrsb.count: %d\n", wrsb.count)
						wrsb.ReadBuf[wrsb.Count] = ks
						wrsb.Count++
					}
				}
			}
		}
	}
}

func WriteKmer(wrfn string, wc <-chan ReadSeqBucket, Kmerlen, numCPU int) {
	outfp, err := os.Create(wrfn)
	if err != nil {
		log.Fatal(err)
	}
	gzwriter := gzip.NewWriter(outfp)
	// bufwriter := bufio.NewWriter(gzwriter)
	// defer outfp.Close()
	defer gzwriter.Close()

	// KBntByteNum := (Kmerlen + bnt.NumBaseInByte - 1) / bnt.NumBaseInByte
	endFlagCount := 0
	writeKmerCount := 0
	for {
		rsb := <-wc
		if rsb.Count == 0 {
			endFlagCount++
			if endFlagCount == numCPU {
				break
			} else {
				continue
			}
		}
		// fmt.Printf("[writeKmer] rsb.count: %d\n", rsb.count)
		for i := 0; i < rsb.Count; i++ {
			err := binary.Write(gzwriter, binary.LittleEndian, rsb.ReadBuf[i].Seq)
			if err != nil {
				log.Fatalf("[writeKmer] write kmer seq err: %v\n", err)
			}
			//fmt.Printf("[writeKmer] Seq: %v\n", rsb.ReadBuf[i].Seq)
			/* if n != KBntByteNum {
				log.Fatalf("[writeKmer] n(%d) != KBntByteNum(%d)\n", n, KBntByteNum)
			} */
			writeKmerCount++
		}
	}
	gzwriter.Flush()
	fmt.Printf("[writeKmer] total write kmer number is : %d\n", writeKmerCount)
}

func Trans2Byte(s string) (rb ReadBnt) {
	rb.Length = len(s)
	rb.Seq = make([]byte, (rb.Length+bnt.NumBaseInByte-1)/bnt.NumBaseInByte)
	for i := 0; i < rb.Length; i++ {
		b := bnt.Base2Bnt[s[i]]
		if b > 3 {
			fmt.Printf("[Trans2Byte]found input sequence base '%c' not belong 'ACTG/actg', please check\n", s[i])
			log.Fatal("error found")
		}
		rb.Seq[i/bnt.NumBaseInByte] <<= bnt.NumBitsInBase
		rb.Seq[i/bnt.NumBaseInByte] |= b
	}

	return rb
}

type Options struct {
	utils.ArgsOpt
	CFSize int64
}

func checkArgs(c cli.Command) (opt Options, suc bool) {
	tmp, err := strconv.Atoi(c.Flag("S").String())
	if err != nil {
		log.Fatalf("[checkArgs] argument 'S': %v set error: %v\n", c.Flag("S"), err)
	}
	if tmp < 1024*1024 {
		log.Fatalf("the argument 'S': %v must bigger than 1024 * 1024\n", c.Flag("S"))
	}
	opt.CFSize = int64(tmp)
	suc = true
	return opt, suc
}

func GetReadSeqBucket(libs []LibInfo, cs chan<- ReadSeqBucket) {
	var processNumReads int
	var rsb ReadSeqBucket
	//var bucketCount int
	// iteration read cfgInfo files
	for _, lib := range libs {
		// seqProfile == 1 note Illumina
		if lib.AsmFlag != AllState || lib.SeqProfile != 1 {
			continue
		}
		for _, fn := range lib.FnName {
			sfn := strings.Split(fn, ".")
			tmp := sfn[len(sfn)-2]
			var format bool // fa format == true, fq format == false
			if tmp == "fa" || tmp == "fasta" {
				format = true
			} else if tmp == "fq" || tmp == "fastq" {
				format = false
			} else {
				log.Fatalf("[GetReadSeqBucket] reads file: %v need suffix end with '*.fa.gz | *.fasta.gz | *.fq.gz | *.fastq.gz'\n", fn)
			}
			infile, err := os.Open(fn)
			if err != nil {
				log.Fatal(err)
			}
			gzreader, err := gzip.NewReader(infile)
			if err != nil {
				log.Fatalf("[GetReadSeqBucket] file: %v maybe not gzip compressed file, err: %v\n", fn, err)
			}
			defer gzreader.Close()
			defer infile.Close()
			bufin := bufio.NewReader(gzreader)

			eof := false
			count := 0
			for !eof {
				var line string
				line, err = bufin.ReadString('\n')
				if err != nil {
					if err == io.EOF {
						eof = true
						continue
					} else {
						log.Fatal(err)
					}
				}

				if (format && count%2 == 1) || (!format && count%4 == 1) {
					s := strings.TrimSpace(line)
					bs := Trans2Byte(s)
					if rsb.Count >= ReadSeqSize {
						cs <- rsb
						var nsb ReadSeqBucket
						rsb = nsb
					}
					rsb.ReadBuf[rsb.Count] = bs
					rsb.Count++
					//fmt.Printf("%s\n", s)
					processNumReads++
					if processNumReads%(1024*1024) == 0 {
						fmt.Printf("processed reads number: %d\n", processNumReads)
					}
				}
				count++
			}
		}
	}
	if rsb.Count > 0 {
		cs <- rsb
	}
	// send read finish signal
	close(cs)
}

func CCF(c cli.Command) {
	//fmt.Println(c.Flags(), c.Parent().Flags())
	//argsCheck(c)
	gOpt, suc := utils.CheckGlobalArgs(c.Parent())
	if suc == false {
		log.Fatalf("[CCF] check global Arguments error, opt: %v\n", gOpt)
	}
	opt := Options{gOpt, 0}
	tmp, suc := checkArgs(c)
	if suc == false {
		log.Fatalf("[CCF] check Arguments error, opt: %v\n", tmp)
	}
	opt.CFSize = tmp.CFSize
	fmt.Printf("[CCF] opt: %v\n", opt)
	cpuprofilefp, err := os.Create(opt.Cpuprofile)
	if err != nil {
		log.Fatalf("[CCF] open cpuprofile file: %v failed\n", opt.Cpuprofile)
	}
	pprof.StartCPUProfile(cpuprofilefp)
	defer pprof.StopCPUProfile()
	cfgInfo, err := ParseCfg(opt.CfgFn)
	if err != nil {
		log.Fatalf("[CCF] ParseCfg 'C': %v err :%v\n", opt.CfgFn, err)
	}
	fmt.Println(cfgInfo)

	// make CuckooFilter

	t0 := time.Now()
	cf := cuckoofilter.MakeCuckooFilter(uint64(opt.CFSize), opt.Kmer)
	numCPU := opt.NumCPU
	runtime.GOMAXPROCS(numCPU + 2)
	bufsize := 10
	cs := make(chan ReadSeqBucket, bufsize)
	wc := make(chan ReadSeqBucket, bufsize)
	we := make(chan int)
	defer close(we)
	defer close(wc)
	defer close(cs)
	go GetReadSeqBucket(cfgInfo.Libs, cs)
	for i := 0; i < numCPU; i++ {
		go ParaConstructCF(cf, cs, wc)
	}
	// write goroutinue
	wrfn := opt.Prefix + ".uniqkmerseq.gz"
	WriteKmer(wrfn, wc, opt.Kmer, numCPU)

	// end signal from write goroutinue
	// prefix := c.Parent().Flag("p").String()
	// cfinfofn := prefix + ".cfInfo"
	// if err := cf.WriteCuckooFilterInfo(cfinfofn); err != nil {
	// 	log.Fatal(err)
	// }

	cfmmapfn := opt.Prefix + ".cfmmap"
	err = cf.MmapWriter(cfmmapfn)
	if err != nil {
		log.Fatal(err)
	}
	// output stat
	cf.GetStat()
	t1 := time.Now()
	fmt.Printf("construct CuckooFilter took %v to run\n", t1.Sub(t0))

	if err != nil {
		log.Fatal("error found in the CCF")
	}
}
