package constructcf

import (
	"bufio"
	"compress/gzip"
	"encoding/binary"
	"fmt"
	"ga/bnt"
	"ga/cuckoofilter"
	"github.com/jwaldrip/odin/cli"
	"io"
	"log"
	"os"
	"runtime"
	"strconv"
	"strings"
	"time"
)

const (
	AllState    = 1
	ScaffState  = 2
	GapState    = 3
	EndString   = "end"
	ReadSeqSize = 100
)

type LibInfo struct {
	name          string // name of library
	numberReads   int64  // the number of reads
	diverse       uint8
	paired        uint8
	asmFlag       uint8 // denote which assembly phase used, note 1 used for all step of assembly pipeline, note 2 used for scaffold phase only, 3 used for filling gap only
	seqProfile    uint8 // denote the data origin
	qualBenchmark uint8 // the benchmark of quality score, strength encourage use phred+33
	totalBasesNum int64
	readLen       int
	insertSize    int // paired read insert size
	insertSD      int // Standard Deviation
	//	fnNum         int      // the number of files
	fnName []string // the files name slice
}

type CfgInfo struct {
	maxRdLen int // maximum read length
	minRdLen int // minimum read length
	//	n        int // the number of library
	libs []LibInfo
}

type ReadBnt struct {
	Seq    []byte
	Length int // length of sequence
}

type ReadSeqBucket struct {
	End     bool
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

func (s1 ReadBnt) Equal(s2 ReadBnt) bool {
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
}

func parseCfg(fn string) (cfgInfo CfgInfo, e error) {
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
			if libInfo.name != "" {
				cfgInfo.libs = append(cfgInfo.libs, libInfo)
				var nli LibInfo
				libInfo = nli
			}
		case "max_rd_len":
			v, err = strconv.Atoi(fields[2])
			cfgInfo.maxRdLen = v
		case "min_rd_len":
			v, err = strconv.Atoi(fields[2])
			cfgInfo.minRdLen = v
		case "name":
			libInfo.name = fields[2]
		case "avg_insert_len":
			v, err = strconv.Atoi(fields[2])
			libInfo.insertSize = v
		case "insert_SD":
			v, err = strconv.Atoi(fields[2])
			libInfo.insertSD = v
		case "diverse_rd_len":
			v, err = strconv.Atoi(fields[2])
			libInfo.diverse = uint8(v)
		case "asm_flag":
			v, err = strconv.Atoi(fields[2])
			libInfo.asmFlag = uint8(v)
		case "seq_profile":
			v, err = strconv.Atoi(fields[2])
			libInfo.seqProfile = uint8(v)
		case "qual_benchmark":
			v, err = strconv.Atoi(fields[2])
			libInfo.qualBenchmark = uint8(v)
		case "f1", "f2":
			libInfo.fnName = append(libInfo.fnName, fields[2])
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
	if libInfo.name != "" {
		cfgInfo.libs = append(cfgInfo.libs, libInfo)
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
	tmp := ks
	tmp.Seq = make([]byte, len(ks.Seq))
	copy(tmp.Seq, ks.Seq)
	rs.Length = tmp.Length
	rs.Seq = make([]byte, len(tmp.Seq))
	for i := tmp.Length - 1; i >= 0; i-- {
		base := tmp.Seq[i/bnt.NumBaseInByte] & bnt.BaseMask
		rs.Seq[(rs.Length-i-1)/bnt.NumBaseInByte] <<= bnt.NumBitsInBase
		rs.Seq[(rs.Length-i-1)/bnt.NumBaseInByte] |= (^base & bnt.BaseMask)
		tmp.Seq[i/bnt.NumBaseInByte] >>= bnt.NumBitsInBase
	}
	return rs
}

func paraConstructCF(cf cuckoofilter.CuckooFilter, cs chan ReadSeqBucket, wc chan ReadSeqBucket) {

	var wrsb ReadSeqBucket
	for {
		rsb := <-cs
		if rsb.End == true {
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

func writeKmer(wrfn string, we chan int, wc chan ReadSeqBucket, Kmerlen, numCPU int) {
	outfp, err := os.Create(wrfn)
	if err != nil {
		log.Fatal(err)
	}
	defer outfp.Close()
	gzwriter := gzip.NewWriter(outfp)
	if err != nil {
		log.Fatal(err)
	}
	defer gzwriter.Close()
	// bufwriter := bufio.NewWriter(gzwriter)

	// KBntByteNum := (Kmerlen + bnt.NumBaseInByte - 1) / bnt.NumBaseInByte
	endFlagCount := 0
	writeKmerCount := 0
	for {
		rsb := <-wc
		if rsb.End == true {
			endFlagCount++
			if endFlagCount == numCPU {
				we <- 1
				break
			} else {
				continue
			}
		}
		//fmt.Printf("[writeKmer] rsb.count: %d\n", rsb.count)
		for i := 0; i < rsb.Count; i++ {
			err := binary.Write(gzwriter, binary.LittleEndian, rsb.ReadBuf[i].Seq)
			if err != nil {
				log.Fatalf("[writeKmer] write kmer seq err: %v\n", err)
			}
			// if n != KBntByteNum {
			// 	log.Fatalf("[writeKmer] n(%d) != KBntByteNum(%d)\n", n, KBntByteNum)
			// }
			writeKmerCount++
			// gzwriter.Write([]byte("\n"))
		}
	}
	// bufwriter.Flush()
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

func CCF(c cli.Command) {
	fmt.Println(c.Flags(), c.Parent().Flags())

	//argsCheck(c)
	fnName := c.Parent().Flag("C").String()
	cfgInfo, err := parseCfg(string(fnName))
	if err != nil {
		log.Fatal("[parseCfg] found err")
	}
	fmt.Println(cfgInfo)

	// make CuckooFilter
	var size, kmerlen int
	if size, err = strconv.Atoi(c.Flag("S").String()); err != nil {
		log.Fatal("flag 'S' set error")
	}
	if kmerlen, err = strconv.Atoi(c.Parent().Flag("K").String()); err != nil {
		log.Fatal("flag 'K' set error")
	}
	t0 := time.Now()
	cf := cuckoofilter.MakeCuckooFilter(uint64(size), kmerlen)
	numCPU, err := strconv.Atoi(c.Parent().Flag("t").String())
	runtime.GOMAXPROCS(numCPU + 2)
	bufsize := 10
	cs := make(chan ReadSeqBucket, bufsize)
	wc := make(chan ReadSeqBucket, bufsize)
	we := make(chan int)
	defer close(we)
	defer close(wc)
	defer close(cs)
	for i := 0; i < numCPU; i++ {
		go paraConstructCF(cf, cs, wc)
	}
	// write goroutinue
	wrfn := c.Parent().Flag("p").String() + ".uniqkmerseq.gz"
	go writeKmer(wrfn, we, wc, kmerlen, numCPU)

	var processNumReads int
	var rsb ReadSeqBucket
	//var bucketCount int
	// iteration read cfgInfo files
	for _, lib := range cfgInfo.libs {
		if lib.asmFlag == AllState {
			for _, fn := range lib.fnName {
				infile, err := os.Open(fn)
				if err != nil {
					log.Fatal(err)
				}
				defer infile.Close()
				gzreader, err := gzip.NewReader(infile)
				bufin := bufio.NewReader(gzreader)
				defer func() {
					if err != nil {
						log.Fatal(err)
					}
					gzreader.Close()
				}()

				eof := false
				count := 0
				for !eof {
					var line string
					line, err = bufin.ReadString('\n')
					if err != nil {
						if err == io.EOF {
							err = nil
							eof = true
							continue
						} else {
							log.Fatal(err)
						}
					}

					if count%4 == 1 {
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
	}
	cs <- rsb
	// send read finish signal
	for i := 0; i < numCPU; i++ {
		var nrsb ReadSeqBucket
		nrsb.End = true
		cs <- nrsb
	}

	<-we // end signal from wirte goroutinue
	prefix := c.Parent().Flag("p").String()
	// cfinfofn := prefix + ".cfInfo"
	// if err := cf.WriteCuckooFilterInfo(cfinfofn); err != nil {
	// 	log.Fatal(err)
	// }

	cfmmapfn := prefix + ".cfmmap"
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
