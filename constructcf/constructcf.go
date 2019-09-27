package constructcf

import (
	"bufio"
	"encoding/binary"
	"fmt"
	"io"
	"log"
	"os"
	"runtime"
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
	Seq []byte
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

type ReadInfo struct {
	ID        int64
	Seq       []byte
	Qual      []byte
	Anotition string
}

/*type ChanStruct struct {
	readseq chan []string
	state   chan int // note if goroutinue can return
} */

/*func (k1 KmerBnt) BiggerThan(k2 KmerBnt) bool {
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
}*/

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

func ParseCfg(fn string, correct bool) (cfgInfo CfgInfo, e error) {
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
			if correct == true {
				libInfo.FnName = append(libInfo.FnName, fields[2])
			} else { // correct == false, need read corrected file
				if fields[0] == "f1" {
					idx := strings.LastIndex(fields[2], "1")
					if idx < 0 || (fields[2][idx+1:] != ".fa.br" && fields[2][idx+1:] != ".fasta.br" && fields[2][idx+1:] != ".fq.br" && fields[2][idx+1:] != ".fastq.br") {
						log.Fatalf("[paraProcessReadsFile] fn: %v, must used suffix *[1|2].[fasta|fa|fq|fastq].br\n", fields[2])
					}
					brfn := fields[2][:idx] + ".Correct.fa.br"
					libInfo.FnName = append(libInfo.FnName, brfn)
				}
			}
		case "f":
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

func ExtendKmerBnt2Byte(kb []uint64, kmerlen int) (extRB []byte) {
	if kmerlen <= 0 || len(kb) == 0 {
		log.Fatalf("[ExtendKmerBnt2Byte] rb: %v\n", kb)
	}
	extRB = make([]byte, kmerlen)

	tmp := make([]uint64, kmerlen)
	copy(tmp, kb)
	// fmt.Printf("[ExtendReadBnt2Byte] rb: %v, crb: %v\n", rb, crb)

	for i := kmerlen - 1; i >= 0; i-- {
		base := tmp[i/bnt.NumBaseInUint64] & bnt.BaseMask
		extRB[i] = uint8(base)
		tmp[i/bnt.NumBaseInUint64] >>= bnt.NumBitsInBase
	}

	return extRB
}

func GetReadBntKmer(seq []byte, startPos, kmerlen int) (sb []uint64) {
	//kBnt.Len = kmerlen
	sb = make([]uint64, (kmerlen+bnt.NumBaseInUint64-1)/bnt.NumBaseInUint64)
	//fmt.Printf("[GetReadBntKmer] len(seq): %d, seq: %v\n\tstartPos: %d, kmerlen: %d\n", len(seq), seq, startPos, kmerlen)

	for i := 0; i < kmerlen; i++ {
		sb[i/bnt.NumBaseInUint64] <<= bnt.NumBitsInBase
		sb[i/bnt.NumBaseInUint64] |= uint64(seq[i+startPos])
	}

	return
}

func NoAllocGetReadBntKmer(seq []byte, startPos, kmerlen int, kb KmerBnt) KmerBnt {
	kb.Len = kmerlen
	//kBnt.Seq = make([]uint64, (kBnt.Len+bnt.NumBaseInUint64-1)/bnt.NumBaseInUint64)
	//fmt.Printf("[GetReadBntKmer] len(seq): %d, seq: %v\n\tstartPos: %d, kmerlen: %d\n", len(seq), seq, startPos, kmerlen)
	copy(kb.Seq, seq[startPos:startPos+kmerlen])

	return kb
}

func ReverseCompletBnt(kb []uint64, klen int) (rb []uint64) {
	rb = make([]uint64, len(kb))
	tmp := make([]uint64, len(kb))
	copy(tmp, kb)
	for i := klen - 1; i >= 0; i-- {
		base := tmp[i/bnt.NumBaseInUint64] & bnt.BaseMask
		tmp[i/bnt.NumBaseInUint64] >>= bnt.NumBitsInBase
		rb[(klen-i-1)/bnt.NumBaseInUint64] <<= bnt.NumBitsInBase
		rb[(klen-i-1)/bnt.NumBaseInUint64] |= uint64(bnt.BntRev[base])
	}
	return
}

func GetReverseComplet(kb KmerBnt) KmerBnt {
	var rb KmerBnt
	rb.Len = kb.Len
	rb.Seq = make([]byte, rb.Len)
	for i := 0; i < kb.Len; i++ {
		rb.Seq[kb.Len-1-i] = bnt.BntRev[kb.Seq[i]]
	}
	return rb
}

/*func DeleteLastBaseKmer(kb KmerBnt) (dk KmerBnt, base uint64) {
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
}*/

/*func DeleteFirstBaseKmer(kb KmerBnt) (dk KmerBnt, base uint64) {
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
}*/

/*func GetNextKmer(kb KmerBnt, base uint64, kmerlen int) (next KmerBnt) {
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
func NoAllocGetNextKmer(kb1, kb2 KmerBnt, base uint64, kmerlen int) KmerBnt {
	nLen := (kmerlen + bnt.NumBaseInUint64 - 1) / bnt.NumBaseInUint64
	kb2.Len = kmerlen
	if kb1.Len == kmerlen-1 {
		//kb2.Len = kmerlen
		copy(kb2.Seq, kb1.Seq)
		kb2.Seq[nLen-1] <<= bnt.NumBitsInBase
		kb2.Seq[nLen-1] |= base
		return kb2
	} else if kb1.Len != kmerlen {
		log.Fatalf("[GetNextKmer] length of kb set error, kmerlen:%d,   kb1 length: %v\n", kmerlen, kb1.Len)
	}
	for i := nLen - 1; i >= 0; i-- {
		//fmt.Printf("[GetNextKmer] base: %v, kb.Seq[%d]: %b\n", base, i, kb.Seq[i])
		if i == nLen-1 {
			offset := (uint64(kmerlen)%bnt.NumBaseInUint64 - 1) * bnt.NumBitsInBase
			tb := (kb1.Seq[i] >> offset) & bnt.BaseMask
			mask := uint64(1<<offset) - 1
			kb2.Seq[i] = kb1.Seq[i] & mask
			kb2.Seq[i] = (kb2.Seq[i] << bnt.NumBitsInBase) | base
			base = tb
		} else {
			tb := (kb1.Seq[i] >> 62) & bnt.BaseMask
			kb2.Seq[i] = (kb1.Seq[i] << bnt.NumBitsInBase) | base
			base = tb
		}
		//fmt.Printf("[GetNextKmer] base: %v, next.Seq[%d]: %b\n", base, i, next.Seq[i])
	}
	return kb2
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

func NoAllocGetPreviousKmer(rb1, rb2 KmerBnt, base uint64, kmerlen int) KmerBnt {
	nLen := (kmerlen + bnt.NumBaseInUint64 - 1) / bnt.NumBaseInUint64
	rb2.Len = kmerlen
	if rb1.Len != kmerlen-1 && rb1.Len != kmerlen {
		log.Fatalf("[GetPreviousKmer] length of rb set error, kmerlen:%d,   rb1 length: %v\n", kmerlen, rb1.Len)
	}
	for i := 0; i < nLen; i++ {
		if i == nLen-1 {
			offset := (uint64(kmerlen)%bnt.NumBaseInUint64 - 1) * bnt.NumBitsInBase
			if rb1.Len == kmerlen {
				rb2.Seq[i] = rb1.Seq[i] >> bnt.NumBitsInBase
			} else {
				rb2.Seq[i] = rb1.Seq[i]
			}
			rb2.Seq[i] |= (base << offset)
		} else {
			tb := rb1.Seq[i] & bnt.BaseMask
			rb2.Seq[i] = rb1.Seq[i] >> bnt.NumBitsInBase
			rb2.Seq[i] |= (base << 62)
			base = tb
		}
	}
	return rb2
}*/

func BiggerThan(kb, rb []byte) bool {
	if len(kb) != len(rb) {
		log.Fatalf("[BiggerThan]len(kb):%v != len(rb):%v\n", len(kb), len(rb))
	}
	for i := 0; i < len(kb); i++ {
		if kb[i] > rb[i] {
			return true
		} else if kb[i] < rb[i] {
			return false
		}
	}
	return false
}

func BiggerThanBnt(kb, rb []uint64) bool {
	if len(kb) != len(rb) {
		log.Fatalf("[BiggerThan]len(kb):%v != len(rb):%v\n", len(kb), len(rb))
	}
	for i := 0; i < len(kb); i++ {
		if kb[i] > rb[i] {
			return true
		} else if kb[i] < rb[i] {
			return false
		}
	}
	return false
}

func ParaConstructCF(cf cuckoofilter.CuckooFilter, cs <-chan ReadSeqBucket, wc chan<- KmerBntBucket) {
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
		/*if rsb.Count < ReadSeqSize {
			fmt.Printf("[ParaConstructCF]rsb.ReadBuf length is :%d\n", rsb.Count))
		}*/
		for i := 0; i < rsb.Count; i++ {
			rBntSeq := rsb.ReadBuf[i]
			//fmt.Printf("[ParaConstructCF]rBntSeq :%v\n", rBntSeq)
			lenS := len(rBntSeq)
			var Bnt KmerBnt
			Bnt.Len = len(rBntSeq)
			Bnt.Seq = rBntSeq
			//kb1 = NoAllocGetReadBntKmer(rBntSeq, 0, cf.Kmerlen-1, kb1)
			cBntSeq := GetReverseComplet(Bnt)
			/*ks := GetReadBntKmer(rBntSeq, 0, cf.Kmerlen-1)
			rs := ReverseComplet(ks)
			if !reflect.DeepEqual(kb1.Seq, ks.Seq) {
				log.Fatalf("[ParaConstructCF]kb1.Seq: %v != ks.Seq: %v\n", kb1, ks)
			}*/
			for j := 0; j < lenS-(cf.Kmerlen-1); j++ {
				kb := rBntSeq[j : j+cf.Kmerlen]
				rb := cBntSeq.Seq[cBntSeq.Len-cf.Kmerlen-j : cBntSeq.Len-j]
				//kb2 = NoAllocGetNextKmer(kb1, kb2, uint64(rBntSeq[j]), cf.Kmerlen)
				//rb2 = NoAllocGetPreviousKmer(rb1, rb2, uint64(bnt.BntRev[rBntSeq[j]]), cf.Kmerlen)
				/*ks = GetNextKmer(ks, uint64(rBntSeq[j]), cf.Kmerlen)
					rs = GetPreviousKmer(rs, uint64(bnt.BntRev[rBntSeq[j]]), cf.Kmerlen)
					if !reflect.DeepEqual(kb2.Seq, ks.Seq) {
						log.Fatalf("[ParaConstructCF]kb2.Seq: %v != ks.Seq: %v\n", kb2, ks)
					}
					if !reflect.DeepEqual(rb2.Seq, rs.Seq) {
						log.Fatalf("[ParaConstructCF]rb2.Seq: %v != rs.Seq: %v\n", rb2, rs)
					}
				  extRBntkb1 := ExtendKmerBnt2Byte(kb1)
				  extRBntkb2 := ExtendKmerBnt2Byte(kb2)
				  extRBntrb1 := ExtendKmerBnt2Byte(rb1)
				  extRBntrb2 := ExtendKmerBnt2Byte(rb2)
				  fmt.Printf("[ParaConstructCF] j: %d, rBntSeq[j]: %v\n\tkb1: %v\n\tkb2: %v\n\trb1: %v\n\trb2: %v\n", j, rBntSeq[j], extRBntkb1, extRBntkb2, extRBntrb1, extRBntrb2)
				*/
				min := kb
				if BiggerThan(kb, rb) {
					min = rb
				}
				//fmt.Printf("ks: %v, rs: %v\n", ks, rs)
				count, suc := cf.Insert(min)
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
					//fmt.Printf("[ParaConstructCF] wrsb.count: %d\n", wrsb.Count)
					var nb KmerBnt
					nb.Len = cf.Kmerlen
					nb.Seq = make([]byte, nb.Len)
					copy(nb.Seq, min)
					wrsb.KmerBntBuf[wrsb.Count] = nb
					wrsb.Count++
				}
			}
		}
	}
}

func ConcurrentConstructCF(fn string, cf cuckoofilter.CuckooFilter, wc chan<- KmerBntBucket, concurrentNum int, kmerlen int, processT chan int) {
	bufSize := 30
	cs := make(chan ReadSeqBucket, bufSize)
	for i := 0; i < concurrentNum; i++ {
		go ParaConstructCF(cf, cs, wc)
	}
	GetReadSeqBucket(fn, cs, kmerlen)
	processT <- 1
}

func WriteKmer(wrfn string, wc <-chan KmerBntBucket, Kmerlen, totalThreadsNum int) {
	outfp, err := os.Create(wrfn)
	if err != nil {
		log.Fatal(err)
	}
	defer outfp.Close()
	brfp := cbrotli.NewWriter(outfp, cbrotli.WriterOptions{Quality: 1})
	defer brfp.Close()
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
			if endFlagCount == totalThreadsNum {
				break
			} else {
				continue
			}
		}
		//fmt.Printf("[writeKmer] rsb.Count: %d\n", rsb.Count)
		for i := 0; i < rsb.Count; i++ {
			err := binary.Write(buffp, binary.LittleEndian, rsb.KmerBntBuf[i].Seq)
			if err != nil {
				log.Fatalf("[writeKmer] write kmer seq err: %v\n", err)
			}
			//fmt.Printf("[writeKmer] Seq: %v\n", rsb.KmerBntBuf[i].Seq)
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
		if b == 'N' {
			return ls[:i]
		}
		//ls = append(ls, alphabet.Letter(bnt.BitNtCharUp[b]))
		ls[i] = bnt.Base2Bnt[b]
	}

	return ls
}

type Options struct {
	utils.ArgsOpt
	CFSize  int64
	Correct bool
}

func checkArgs(c cli.Command) (opt Options, suc bool) {
	tmp, err := strconv.Atoi(c.Flag("S").String())
	if err != nil {
		log.Fatalf("[checkArgs] argument 'S': %v set error: %v\n", c.Flag("S"), err)
	}
	if tmp < 1024*1024 {
		log.Fatalf("[checkArgs]the argument 'S': %v must bigger than 1024 * 1024\n", c.Flag("S"))
	}
	opt.CFSize = int64(tmp)
	cor := c.Flag("Correct").Get().(bool)
	if cor == true {
		opt.Correct = true
	} else if cor == false {
		opt.Correct = false
	} else {
		log.Fatalf("[checkArgs] argument 'Correct': %v set error, must set true|false\n", c.Flag("Correct"))
	}
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

func GetReadFileRecord(buffp *bufio.Reader, format string, annotion bool) (ri ReadInfo, err error) {
	var blockLineNum int
	if format == "fa" {
		blockLineNum = 2
	} else { // format == "fq"
		blockLineNum = 4
	}
	b := make([][]byte, blockLineNum)
	i := 0
	for ; i < blockLineNum; i++ {
		b[i], err = buffp.ReadBytes('\n')
		if err != nil {
			break
		}
	}
	if err != nil {
		if err == io.EOF {
			if i == 0 {
				return
			}
			if i != blockLineNum-1 {
				log.Fatalf("[GetReadSeqBucket] x: %d, fp: %v found not unbroken record\n", i, buffp)
			}
		} else {
			log.Fatalf("[GetReadSeqBucket] fp : %v encounter err: %v\n", buffp, err)
		}
	}
	flist := strings.Fields(string(b[0]))
	id, err2 := strconv.Atoi(flist[0][1:])
	if err2 != nil {
		log.Fatalf("[LoadNGSReads] load fp: '%v' file, read ID: %v not digits, please convert to digits...\n", buffp, flist[0][1:])
	}
	ri.ID = int64(id)
	if annotion {
		ri.Anotition = strings.Join(flist[1:], "\t")
	}
	ri.Seq = Transform2BntByte(b[1][:len(b[1])-1])
	if format == "fq" {
		ri.Qual = b[3][:len(ri.Seq)]
	}

	return
}

func GetReadSeqBucket(fn string, cs chan<- ReadSeqBucket, kmerlen int) {
	var processNumReads int
	var rsb ReadSeqBucket
	//var bucketCount int
	format := GetReadsFileFormat(fn)
	fp, err := os.Open(fn)
	if err != nil {
		log.Fatalf("[GetReadSeqBucket] file: %v open failed, err: %v\n", fn, err)
	}
	defer fp.Close()
	fmt.Printf("[GetReadSeqBucket] processed reads in file: %v\n", fn)
	brfp := cbrotli.NewReaderSize(fp, 1<<25)
	defer brfp.Close()
	buffp := bufio.NewReader(brfp)

	//var count int
	var err1 error
	for err1 != io.EOF {
		ri, err1 := GetReadFileRecord(buffp, format, false)
		if ri.ID == 0 {
			if err1 == io.EOF {
				break
			} else {
				log.Fatalf("[GetReadSeqBucket] file: %s encounter err: %v\n", fn, err1)
			}
		}
		if rsb.Count >= ReadSeqSize {
			cs <- rsb
			var nsb ReadSeqBucket
			rsb = nsb
		}
		if len(ri.Seq) < kmerlen+10 {
			continue
		}
		rsb.ReadBuf[rsb.Count] = ri.Seq
		rsb.Count++
		processNumReads++
	}
	if rsb.Count > 0 {
		cs <- rsb
	}
	// send read finish signal
	fmt.Printf("[GetReadSeqBucket] processed reads number is: %d, finished processed file: %v\n", processNumReads, fn)
	close(cs)
}

func CCF(c cli.Command) {
	CCFFunc(c)
}

func CCFFunc(c cli.Command) cuckoofilter.CuckooFilter {
	fmt.Println(c.Flags(), c.Parent().Flags())
	//argsCheck(c)
	gOpt, suc := utils.CheckGlobalArgs(c.Parent())
	if suc == false {
		log.Fatalf("[CCF] check global Arguments error, opt: %v\n", gOpt)
	}
	opt := Options{gOpt, 0, false}
	tmp, suc := checkArgs(c)
	if suc == false {
		log.Fatalf("[CCF] check Arguments error, opt: %v\n", tmp)
	}
	opt.CFSize = tmp.CFSize
	opt.Correct = tmp.Correct
	fmt.Printf("[CCF] opt: %v\n", opt)
	/*profileFn := opt.Prefix + ".CCF.prof"
	cpuprofilefp, err := os.Create(profileFn)
	if err != nil {
		log.Fatalf("[CCF] open cpuprofile file: %v failed\n", profileFn)
	}
	pprof.StartCPUProfile(cpuprofilefp)
	defer pprof.StopCPUProfile()*/
	cfgInfo, err := ParseCfg(opt.CfgFn, opt.Correct)
	if err != nil {
		log.Fatalf("[CCF] ParseCfg 'C': %v err :%v\n", opt.CfgFn, err)
	}
	fmt.Println(cfgInfo)

	// make CuckooFilter

	t0 := time.Now()
	cf := cuckoofilter.MakeCuckooFilter(uint64(opt.CFSize), opt.Kmer)
	numCPU := opt.NumCPU
	runtime.GOMAXPROCS(numCPU + 2)
	bufsize := 60000
	wc := make(chan KmerBntBucket, bufsize)
	defer close(wc)
	//we := make(chan int)
	//defer close(we)
	var fnArr []string
	for _, lib := range cfgInfo.Libs {
		if lib.AsmFlag != AllState && lib.SeqProfile != 1 {
			continue
		}
		fnArr = append(fnArr, lib.FnName...)
	}

	//fmt.Printf("[CCF] fileName array: %v\n", fnArr)
	totalFileNum := len(fnArr)
	//concurrentNum := 6
	//totalNumT := opt.NumCPU/(concurrentNum+1) + 1
	concurrentNum := 6
	/*if opt.Kmer < 128 {
		concurrentNum = 6
	} else {
		concurrentNum = 12
	}*/
	totalNumT := opt.NumCPU/(concurrentNum+1) + 1
	if totalNumT > totalFileNum {
		totalNumT = totalFileNum
	}

	processT := make(chan int, totalNumT)
	defer close(processT)
	for i := 0; i < totalNumT; i++ {
		processT <- 1
	}

	// write goroutinue
	wrfn := opt.Prefix + ".uniqkmerseq.br"
	go WriteKmer(wrfn, wc, opt.Kmer, totalFileNum*concurrentNum)

	for _, fn := range fnArr {
		<-processT
		//fmt.Printf("[CCF] processing file: %v\n", lib.FnName[i])
		go ConcurrentConstructCF(fn, cf, wc, concurrentNum, opt.Kmer, processT)
	}

	for i := 0; i < totalNumT; i++ {
		<-processT
	}
	time.Sleep(time.Second * 3)
	fmt.Printf("[CCF] construct CuckooFilter took %v to run\n", time.Now().Sub(t0))

	// end signal from write goroutinue
	// prefix := c.Parent().Flag("p").String()
	// cfinfofn := prefix + ".cfInfo"
	// if err := cf.WriteCuckooFilterInfo(cfinfofn); err != nil {
	// 	log.Fatal(err)
	// }
	t0 = time.Now()
	if !opt.Correct {
		cfinfofn := opt.Prefix + ".cf.Info"
		err = cf.WriteCuckooFilterInfo(cfinfofn)
		if err != nil {
			log.Fatalf("[CCF]WriteCuckooFilterInfo file: %v error: %v\n", cfinfofn, err)
		}
		cffn := opt.Prefix + ".cf.Hash.br"
		err = cf.HashWriter(cffn)
		if err != nil {
			log.Fatalf("[CCF]HashWriter file: %v error: %v\n", cffn, err)
		}
		fmt.Printf("[CCF] write CuckooFilter hash used:  %v to run\n", time.Now().Sub(t0))
	}

	// output stat
	cf.GetStat()
	if err != nil {
		log.Fatal("error found in the CCF")
	}
	return cf
}
