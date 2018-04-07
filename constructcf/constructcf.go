package constructcf

import (
	"bufio"
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
	"github.com/mudesheng/ga/cbrotli"
	"github.com/mudesheng/ga/cuckoofilter"
	"github.com/mudesheng/ga/utils"
)

const (
	AllState    = 1
	ScaffState  = 2
	GapState    = 3
	EndString   = "end"
	ReadSeqSize = 1000
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

/*type ReadBnt struct {
	Seq    []byte
	Length int // length of sequence
}*/

// KmerBnt compact kmer base sequence to the uint64 Array
type KmerBnt struct {
	Seq []uint64
	Len int // length of Sequence
}

type KmerBntBucket struct {
	KmerBntBuf [ReadSeqSize]KmerBnt
	Count      int
}

type ReadSeqBucket struct {
	ReadBuf [ReadSeqSize][]byte
	Count   int
}

/*type ChanStruct struct {
	readseq chan []string
	state   chan int // note if goroutinue can return
} */

func (k1 KmerBnt) BiggerThan(k2 KmerBnt) bool {
	if k1.Len < k2.Len {
		return false
	} else if k1.Len > k2.Len {
		return true
	} else {
		for i := 0; i < len(k1.Seq); i++ {
			if k1.Seq[i] < k2.Seq[i] {
				return false
			} else if k1.Seq[i] > k2.Seq[i] {
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

func ExtendKmerBnt2Byte(rb KmerBnt) (extRB []byte) {
	if rb.Len <= 0 || len(rb.Seq) == 0 {
		log.Fatalf("[ExtendKmerBnt2Byte] rb: %v\n", rb)
	}
	extRB = make([]byte, rb.Len)
	var crb KmerBnt
	crb.Len = rb.Len
	crb.Seq = make([]uint64, len(rb.Seq))
	copy(crb.Seq, rb.Seq)
	// fmt.Printf("[ExtendReadBnt2Byte] rb: %v, crb: %v\n", rb, crb)

	for i := crb.Len - 1; i >= 0; i-- {
		base := crb.Seq[i/bnt.NumBaseInUint64] & bnt.BaseMask
		extRB[i] = uint8(base)
		crb.Seq[i/bnt.NumBaseInUint64] >>= bnt.NumBitsInBase
	}

	return extRB
}

func GetReadBntKmer(seq []byte, startPos, kmerlen int) (kBnt KmerBnt) {
	kBnt.Len = kmerlen
	kBnt.Seq = make([]uint64, (kBnt.Len+bnt.NumBaseInUint64-1)/bnt.NumBaseInUint64)

	for i := 0; i < kmerlen; i++ {
		kBnt.Seq[i/bnt.NumBaseInUint64] <<= bnt.NumBitsInBase
		kBnt.Seq[i/bnt.NumBaseInUint64] |= uint64(seq[i+startPos])
	}

	return
}

func ReverseComplet(ks KmerBnt) (rs KmerBnt) {
	rs.Len = ks.Len
	rs.Seq = make([]uint64, len(ks.Seq))
	var tmp KmerBnt
	tmp.Len = ks.Len
	tmp.Seq = make([]uint64, len(ks.Seq))
	copy(tmp.Seq, ks.Seq)
	for i := tmp.Len - 1; i >= 0; i-- {
		base := tmp.Seq[i/bnt.NumBaseInUint64] & bnt.BaseMask
		tmp.Seq[i/bnt.NumBaseInUint64] >>= bnt.NumBitsInBase
		rs.Seq[(rs.Len-i-1)/bnt.NumBaseInUint64] <<= bnt.NumBitsInBase
		rs.Seq[(rs.Len-i-1)/bnt.NumBaseInUint64] |= uint64(bnt.BntRev[base])
	}
	return rs
}

func DeleteLastBaseKmer(kb KmerBnt) (dk KmerBnt, base uint64) {
	nLen := (kb.Len + bnt.NumBaseInUint64 - 1) / bnt.NumBaseInUint64
	dk.Len = kb.Len - 1
	base = kb.Seq[nLen-1] & bnt.BaseMask
	if nLen > (kb.Len-1+bnt.NumBaseInUint64-1)/bnt.NumBaseInUint64 {
		dk.Seq = make([]uint64, nLen-1)
		copy(dk.Seq, kb.Seq)
	} else {
		dk.Seq = make([]uint64, nLen)
		copy(dk.Seq, kb.Seq)
		dk.Seq[nLen-1] >>= bnt.NumBitsInBase
	}
	return
}

func DeleteFirstBaseKmer(kb KmerBnt) (dk KmerBnt, base uint64) {
	nLen := (kb.Len + bnt.NumBaseInUint64 - 1) / bnt.NumBaseInUint64
	dk.Len = kb.Len - 1
	seqLen := (dk.Len + bnt.NumBaseInUint64 - 1) / bnt.NumBaseInUint64
	dk.Seq = make([]uint64, seqLen)
	for i := nLen - 1; i >= 0; i-- {
		if i == nLen-1 {
			offset := (uint64(kb.Len)%bnt.NumBaseInUint64 - 1) * bnt.NumBitsInBase
			base = (kb.Seq[i] >> offset) & bnt.BaseMask
			mask := uint64(1<<offset) - 1
			if seqLen == nLen {
				dk.Seq[i] = kb.Seq[i] & mask
			}
		} else {
			tb := (kb.Seq[i] >> 62) & bnt.BaseMask
			dk.Seq[i] = (kb.Seq[i] << bnt.NumBitsInBase) | base
			base = tb
		}
	}

	return
}

func GetNextKmer(kb KmerBnt, base uint64, kmerlen int) (next KmerBnt) {
	nLen := (kmerlen + bnt.NumBaseInUint64 - 1) / bnt.NumBaseInUint64
	next.Len = kmerlen
	next.Seq = make([]uint64, nLen)
	if kb.Len == kmerlen-1 {
		copy(next.Seq, kb.Seq)
		next.Seq[nLen-1] <<= bnt.NumBitsInBase
		next.Seq[nLen-1] |= base
		return
	} else if kb.Len != kmerlen {
		log.Fatalf("[GetNextKmer] length of kb set error, kmerlen:%d,   kb length: %v\n", kmerlen, kb.Len)
	}
	for i := nLen - 1; i >= 0; i-- {
		//fmt.Printf("[GetNextKmer] base: %v, kb.Seq[%d]: %b\n", base, i, kb.Seq[i])
		if i == nLen-1 {
			offset := (uint64(kmerlen)%bnt.NumBaseInUint64 - 1) * bnt.NumBitsInBase
			tb := (kb.Seq[i] >> offset) & bnt.BaseMask
			mask := uint64(1<<offset) - 1
			next.Seq[i] = kb.Seq[i] & mask
			next.Seq[i] = (next.Seq[i] << bnt.NumBitsInBase) | base
			base = tb
		} else {
			tb := (kb.Seq[i] >> 62) & bnt.BaseMask
			next.Seq[i] = (kb.Seq[i] << bnt.NumBitsInBase) | base
			base = tb
		}
		//fmt.Printf("[GetNextKmer] base: %v, next.Seq[%d]: %b\n", base, i, next.Seq[i])
	}
	return
}

func GetPreviousKmer(kb KmerBnt, base uint64, kmerlen int) (pkb KmerBnt) {
	nLen := (kmerlen + bnt.NumBaseInUint64 - 1) / bnt.NumBaseInUint64
	pkb.Len = kmerlen
	pkb.Seq = make([]uint64, nLen)
	if kb.Len != kmerlen-1 && kb.Len != kmerlen {
		log.Fatalf("[GetPreviousKmer] length of kb set error, kmerlen:%d,   kb length: %v\n", kmerlen, kb.Len)
	}
	for i := 0; i < nLen; i++ {
		if i == nLen-1 {
			offset := (uint64(kmerlen)%bnt.NumBaseInUint64 - 1) * bnt.NumBitsInBase
			if kb.Len == kmerlen {
				pkb.Seq[i] = kb.Seq[i] >> bnt.NumBitsInBase
			} else {
				pkb.Seq[i] = kb.Seq[i]
			}
			pkb.Seq[i] |= (base << offset)
		} else {
			tb := kb.Seq[i] & bnt.BaseMask
			pkb.Seq[i] = kb.Seq[i] >> bnt.NumBitsInBase
			pkb.Seq[i] |= (base << 62)
			base = tb
		}
	}
	return
}

func ParaConstructCF(cf cuckoofilter.CuckooFilter, cs chan ReadSeqBucket, wc chan KmerBntBucket) {

	var wrsb KmerBntBucket
	for {
		rsb, ok := <-cs
		if !ok {
			if wrsb.Count > 0 {
				wc <- wrsb
			}
			var tmp KmerBntBucket
			wc <- tmp
			break
		}
		if rsb.Count < ReadSeqSize {
			fmt.Printf("[ParaConstructCF]rsb.ReadBuf length is :%d\n", len(rsb.ReadBuf))
		}
		for i := 0; i < rsb.Count; i++ {
			rBntSeq := rsb.ReadBuf[i]
			//fmt.Printf("[ParaConstructCF]rBntSeq :%v\n", rBntSeq)
			lenS := len(rBntSeq)
			ks := GetReadBntKmer(rBntSeq, 0, cf.Kmerlen-1)
			rs := ReverseComplet(ks)
			for j := cf.Kmerlen - 1; j < lenS; j++ {
				ks = GetNextKmer(ks, uint64(rBntSeq[j]), cf.Kmerlen)
				rs = GetPreviousKmer(rs, uint64(bnt.BntRev[rBntSeq[j]]), cf.Kmerlen)
				//extRBntks := ExtendKmerBnt2Byte(ks)
				//extRBntrs := ExtendKmerBnt2Byte(rs)
				//fmt.Printf("[ParaConstructCF] j: %d, kmerSeq: %v\n\tks: %v\n\trs: %v\n", j, rBntSeq[j+1-cf.Kmerlen:j+1], extRBntks, extRBntrs)
				min := ks
				if ks.BiggerThan(rs) {
					min = rs
				}
				//fmt.Printf("ks: %v, rs: %v\n", ks, rs)
				count, suc := cf.Insert(min.Seq)
				if suc == false {
					log.Fatal("[ParaConstructCF] Insert to the CuckooFilter false")
				}
				//fmt.Printf("retrun count : %d\n", count)
				//fmt.Printf("count set: %d\n", cf.GetCount(ks.Seq))
				if count == 2 {
					if wrsb.Count >= ReadSeqSize {
						wc <- wrsb
						var nrsb KmerBntBucket
						wrsb = nrsb
					}
					//fmt.Printf("[paraConstructCF] wrsb.count: %d\n", wrsb.count)
					wrsb.KmerBntBuf[wrsb.Count] = ks
					wrsb.Count++
				}
			}
		}
	}
}

func WriteKmer(wrfn string, wc <-chan KmerBntBucket, Kmerlen, numCPU int) {
	outfp, err := os.Create(wrfn)
	if err != nil {
		log.Fatal(err)
	}
	defer outfp.Close()
	brfp := cbrotli.NewWriter(outfp, cbrotli.WriterOptions{Quality: 1})
	// bufwriter := bufio.NewWriter(gzwriter)
	// defer outfp.Close()
	buffp := bufio.NewWriterSize(brfp, 1<<25)
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
			err := binary.Write(buffp, binary.LittleEndian, rsb.KmerBntBuf[i].Seq)
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
	if err := buffp.Flush(); err != nil {
		log.Fatalf("[writeKmer] write kmer seq err: %v\n", err)
	}
	if err := brfp.Flush(); err != nil {
		log.Fatalf("[writeKmer] write kmer seq err: %v\n", err)
	}

	fmt.Printf("[writeKmer] total write kmer number is : %d\n", writeKmerCount)
}

/*func Trans2Byte(s string) (rb ReadBnt) {
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
}*/

func Transform2BntByte(ks []byte) []byte {
	ls := make([]byte, len(ks))
	for i, b := range ks {
		//ls = append(ls, alphabet.Letter(bnt.BitNtCharUp[b]))
		ls[i] = bnt.Base2Bnt[b]
	}

	return ls
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

func GetReadsFileFormat(fn string) (format string) {
	sfn := strings.Split(fn, ".")
	if len(sfn) < 3 {
		log.Fatalf("[GetReadsFileFormat] reads file: %v need suffix end with '*.fa.br | *.fasta.br | *.fq.br | *.fastq.br'\n", fn)
	}
	tmp := sfn[len(sfn)-2]
	if tmp == "fa" || tmp == "fasta" {
		format = "fa"
	} else if tmp == "fq" || tmp == "fastq" {
		format = "fq"
	} else {
		log.Fatalf("[GetReadsFileFormat] reads file: %v need suffix end with '*.fa.gz | *.fasta.gz | *.fq.gz | *.fastq.gz'\n", fn)
	}

	return format
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
			format := GetReadsFileFormat(fn)
			fp, err := os.Open(fn)
			if err != nil {
				log.Fatalf("[GetReadSeqBucket] file: %v open failed, err: %v\n", fn, err)
			}
			defer fp.Close()
			brfp := cbrotli.NewReaderSize(fp, 1<<25)
			defer brfp.Close()
			buffp := bufio.NewReader(brfp)
			var blockLineNum int
			if format == "fa" {
				blockLineNum = 2
			} else {
				blockLineNum = 4
			}

			//var count int
			var EOF error
			for EOF != io.EOF {
				var b [][]byte
				b = make([][]byte, blockLineNum)
				var x int
				var err1 error
				for ; x < blockLineNum; x++ {
					b[x], err1 = buffp.ReadBytes('\n')
					if err1 != nil {
						break
					}
				}
				if err1 != nil {
					if err1 == io.EOF {
						if x == 0 {
							break
						}
						if x != blockLineNum-1 {
							log.Fatalf("[GetReadSeqBucket] x: %d, file: %s found not unbroken record\n", x, fn)
						}
						EOF = io.EOF
					} else {
						log.Fatalf("[GetReadSeqBucket] file : %v encounter err: %v\n", fn, err1)
					}
				}
				/*id, err2 := strconv.Atoi(string(b[0][1 : len(b[0])-1]))
				if err2 != nil {
					log.Fatalf("[LoadNGSReads] load fn: '%v' file, read ID: %v not digits, please convert to digits...\n", fn, b[0][:len(b[0])-1])
				}*/
				bs := Transform2BntByte(b[1][:len(b[1])-1])
				if rsb.Count >= ReadSeqSize {
					cs <- rsb
					var nsb ReadSeqBucket
					rsb = nsb
				}
				rsb.ReadBuf[rsb.Count] = bs
				rsb.Count++
				processNumReads++
			}
			if rsb.Count > 0 {
				cs <- rsb
			}
		}
	}
	// send read finish signal
	fmt.Printf("[GetReadSeqBucket] processed reads number is: %d\n", processNumReads)
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
	profileFn := opt.Prefix + ".CCF.prof"
	cpuprofilefp, err := os.Create(profileFn)
	if err != nil {
		log.Fatalf("[CCF] open cpuprofile file: %v failed\n", profileFn)
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
	bufsize := 20
	cs := make(chan ReadSeqBucket, bufsize)
	wc := make(chan KmerBntBucket, bufsize)
	//we := make(chan int)
	//defer close(we)
	defer close(wc)
	go GetReadSeqBucket(cfgInfo.Libs, cs)
	for i := 0; i < numCPU; i++ {
		go ParaConstructCF(cf, cs, wc)
	}
	// write goroutinue
	wrfn := opt.Prefix + ".uniqkmerseq.br"
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
