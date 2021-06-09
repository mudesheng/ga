package constructcf

import (
	"bufio"
	"bytes"
	"compress/gzip"
	"fmt"
	"io"
	"log"
	"os"
	"runtime"
	"runtime/pprof"
	"strconv"
	"strings"
	"sync"
	"time"

	"github.com/jwaldrip/odin/cli"
	"github.com/mudesheng/ga/bnt"
	"github.com/mudesheng/ga/cbrotli"
	//"github.com/google/brotli/go/cbrotli"
	"github.com/mudesheng/ga/cuckoofilter"
	"github.com/mudesheng/ga/utils"
)

const (
	AllState          = 1
	ScaffState        = 2
	GapState          = 3
	EndString         = "end"
	ReadSeqSize       = 1000
	UniqueKmerFileNum = 10
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
	InsertSize    int   // paired read insert size
	InsertSD      int   // Standard Deviation
	Merged        uint8 // set library merge flag, 0 note not merged, 1 need merged, 2 has merged
	//	fnNum         int      // the number of files
	FnName     []string // the files name slice
	UnMergedFn []string
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

type ReadInfoS struct {
	ID        []byte
	Seq       []byte
	Qual      []byte
	Anotition []byte
}

type ReadInfo struct {
	ID        uint32
	Seq       []byte
	Qual      []byte
	Anotition []byte
}

type ReadInfoSBuf struct {
	Arr []ReadInfo
	Buf []byte
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

func ParseCfg(fn string, merge, correct bool) (cfgInfo CfgInfo, e error) {
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
		case "merge_flag": // = 0 note not merged, =1 need merge, =2 has been merged
			v, err = strconv.Atoi(fields[2])
			libInfo.Merged = uint8(v)
		case "p":
			fs := strings.Split(fields[2], ".")
			if len(fs) < 2 || fs[len(fs)-1] != "br" || (fs[len(fs)-2] != "fa" && fs[len(fs)-2] != "fasta" && fs[len(fs)-2] != "fq" && fs[len(fs)-2] != "fastq") {
				log.Fatalf("[ParseCfg] fn: %v, must used suffix *.[fasta|fa|fq|fastq].br\n", fields[2])
			}
			if correct == true {
				libInfo.FnName = append(libInfo.FnName, fields[2])
			} else if merge == true {
				if libInfo.Merged < 2 {
					brfn := strings.Join(fs[:len(fs)-2], ".") + ".Correct." + fs[len(fs)-2] + "." + fs[len(fs)-1]
					libInfo.FnName = append(libInfo.FnName, brfn)
				} else {
					brfn := strings.Join(fs[:len(fs)-2], ".") + ".Merged." + fs[len(fs)-2] + "." + fs[len(fs)-1]
					libInfo.FnName = append(libInfo.FnName, brfn)
					brfn = strings.Join(fs[:len(fs)-2], ".") + ".UnMerged." + fs[len(fs)-2] + "." + fs[len(fs)-1]
					libInfo.UnMergedFn = append(libInfo.UnMergedFn, brfn)
				}
			} else {
				fs := strings.Split(fields[2], ".")
				if len(fs) < 2 || fs[len(fs)-1] != "br" || (fs[len(fs)-2] != "fa" && fs[len(fs)-2] != "fasta" && fs[len(fs)-2] != "fq" && fs[len(fs)-2] != "fastq") {
					log.Fatalf("[ParseCfg] fn: %v, must used suffix *.[fasta|fa|fq|fastq].br\n", fields[2])
				}
				brfn := strings.Join(fs[:len(fs)-2], ".") + ".Merged." + fs[len(fs)-2] + "." + fs[len(fs)-1]
				libInfo.FnName = append(libInfo.FnName, brfn)
				brfn = strings.Join(fs[:len(fs)-2], ".") + ".UnMerged." + fs[len(fs)-2] + "." + fs[len(fs)-1]
				libInfo.UnMergedFn = append(libInfo.UnMergedFn, brfn)
			}
		case "LR":
			libInfo.FnName = append(libInfo.FnName, fields[2])
		default:
			if fields[0][0] != '#' && fields[0][0] != ';' {
				log.Fatalf("[ParseCfg]noknown line: %s", line)
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

func GetReverseCompletUint64(bs uint64, slen uint32) uint64 {
	bs = ^bs
	bs = (bs&0x3333333333333333)<<2 | (bs&0xCCCCCCCCCCCCCCCC)>>2
	bs = (bs&0x0F0F0F0F0F0F0F0F)<<4 | (bs&0xF0F0F0F0F0F0F0F0)>>4
	bs = (bs&0x00FF00FF00FF00FF)<<8 | (bs&0xFF00FF00FF00FF00)>>8
	bs = (bs&0x0000FFFF0000FFFF)<<16 | (bs&0xFFFF0000FFFF0000)>>16
	bs = (bs&0x00000000FFFFFFFF)<<32 | (bs&0xFFFFFFFF00000000)>>32
	bs >>= (64 - (slen << 1))
	return bs
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

func GetReverseComplet2(kb, rb KmerBnt) KmerBnt {
	rb.Len = kb.Len
	if len(rb.Seq) < rb.Len {
		log.Fatalf("[GetReverseComplet2] len(rb.Seq):%d < rb.Len:%d\n", len(rb.Seq), rb.Len)
	}
	rb.Seq = rb.Seq[:rb.Len]
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

func ParaConstructCF(cf cuckoofilter.CuckooFilter, cs <-chan ReadSeqBucket, wc chan<- KmerBntBucket, MinKmerFreq int) {
	var wrsb KmerBntBucket
	for {
		rsb, _ := <-cs
		if rsb.Count == 0 {
			if wrsb.Count > 0 {
				wc <- wrsb
			}
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
				//cpbMin := CompressToBnt(min)
				//fmt.Printf("ks: %v, rs: %v\n", ks, rs)
				count, suc := cf.Insert(min)
				if suc == false {
					log.Fatal("[ParaConstructCF] Insert to the CuckooFilter false")
				}
				//fmt.Printf("retrun count : %d\n", count)
				//fmt.Printf("count set: %d\n", cf.GetCount(ks.Seq))
				if count == MinKmerFreq-1 {
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

func ConcurrentConstructCF(fn string, cf cuckoofilter.CuckooFilter, wc chan<- KmerBntBucket, concurrentNum int, kmerlen int, processT chan int, MinKmerFreq int) {
	bufSize := 5*concurrentNum + 5
	cs := make(chan ReadSeqBucket, bufSize)
	for i := 0; i < concurrentNum; i++ {
		go ParaConstructCF(cf, cs, wc, MinKmerFreq)
	}
	kmerNum := GetReadSeqBucket(fn, cs, kmerlen)
	for len(cs) > 0 {
		time.Sleep(time.Second)
	}
	time.Sleep(time.Second)
	processT <- kmerNum
}

func WriteKmer(prefix string, wc <-chan KmerBntBucket, Kmerlen int) {
	kmerNumC := make(chan int, 1)
	for i := 0; i < UniqueKmerFileNum; i++ {
		wrfn := prefix + "." + strconv.Itoa(i) + ".uniqkmerseq.br"
		go ParaWriteKmer(wrfn, wc, Kmerlen, kmerNumC)
	}

	kmerNum := 0
	for i := 0; i < UniqueKmerFileNum; i++ {
		kmerNum += <-kmerNumC
	}
	close(kmerNumC)
	fmt.Printf("[writeKmer] total write kmer number is : %d\n", kmerNum)
	return
}

// concurrent write kmer
func ParaWriteKmer(wrfn string, wc <-chan KmerBntBucket, Kmerlen int, kmerNumC chan<- int) {
	outfp, err := os.Create(wrfn)
	if err != nil {
		log.Fatal(err)
	}
	defer outfp.Close()
	brfp := cbrotli.NewWriter(outfp, cbrotli.WriterOptions{Quality: 1})
	//brfp := cbrotli.NewWriter(outfp, cbrotli.WriterOptions{LGWin: 1 << 10})
	defer brfp.Close()
	// bufwriter := bufio.NewWriter(gzwriter)
	// defer outfp.Close()
	buffp := bufio.NewWriterSize(brfp, 1<<20)
	// KBntByteNum := (Kmerlen + bnt.NumBaseInByte - 1) / bnt.NumBaseInByte
	writeKmerCount := 0
	for {
		rsb, _ := <-wc
		if rsb.Count == 0 {
			break
		}
		for i := 0; i < rsb.Count; i++ {
			buffp.Write(rsb.KmerBntBuf[i].Seq)
			writeKmerCount++
		}
	}
	if err := buffp.Flush(); err != nil {
		log.Fatalf("[ParaWriteKmer] write kmer seq err: %v\n", err)
	}
	if err := brfp.Flush(); err != nil {
		log.Fatalf("[ParaWriteKmer] write kmer seq err: %v\n", err)
	}
	kmerNumC <- writeKmerCount
	return
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

func Transform2BntByte2(ks []byte) []byte {
	//ls := make([]byte, len(ks))
	for i, b := range ks {
		if b == 'N' {
			return ks[:i]
		}
		//ls = append(ls, alphabet.Letter(bnt.BitNtCharUp[b]))
		ks[i] = bnt.Base2Bnt[b]
	}

	return ks
}

type Options struct {
	utils.ArgsOpt
	CFSize      int64
	Correct     bool
	Merge       bool
	MinKmerFreq int
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
	opt.Merge = c.Flag("Merge").Get().(bool)

	minKmerFreq := c.Flag("MinKmerFreq").Get().(int)
	if minKmerFreq < 1 || minKmerFreq > cuckoofilter.MAX_C {
		log.Fatalf("[checkArgs]the argument 'MinKmerFreq': %v must [%v ~ %v]\n", minKmerFreq, 1, cuckoofilter.MAX_C)
	}
	opt.MinKmerFreq = minKmerFreq
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

func WriteBr(fn string, wc <-chan []byte, wbytePool *sync.Pool, wfinish chan<- bool) {
	outfp, err1 := os.Create(fn)
	if err1 != nil {
		log.Fatal(err1)
	}
	defer outfp.Close()
	brfp := cbrotli.NewWriter(outfp, cbrotli.WriterOptions{Quality: 1, LGWin: 21})
	defer brfp.Close()
	//outbuffp := bufio.NewWriterSize(brfp, 1<<20)
	for {
		wb, ok := <-wc
		if !ok {
			break
		}
		brfp.Write(wb)
		wbytePool.Put(wb)
	}

	if err := brfp.Flush(); err != nil {
		log.Fatalf("[WriteBrFa] write read seq err: %v\n", err)
	}
	wfinish <- true
}

/*func WriteSplitBrFa(prefix, suffix string, wc <-chan ReadInfo, splitRecordNum int) {
	fileNum := 1
	fn := prefix + "_" + strconv.Itoa(fileNum) + suffix
	fileNum++
	outfp, err1 := os.Create(fn)
	if err1 != nil {
		log.Fatal(err1)
	}
	recordCount := 0
	//defer
	brfp := cbrotli.NewWriter(outfp, cbrotli.WriterOptions{Quality: 1})
	brfp.Write(p)
	//defer
	outbuffp := bufio.NewWriterSize(brfp, 1<<20)
	for {
		ri, ok := <-wc
		if !ok {
			break
		}
		id := strconv.Itoa(int(ri.ID))
		s1 := ">" + id + "\n" + string(ri.Seq) + "\n"
		outbuffp.WriteString(s1)
		recordCount++
		if recordCount > 0 && recordCount%splitRecordNum == 0 {
			if err := outbuffp.Flush(); err != nil {
				log.Fatalf("[WriteBrFa] write read seq err: %v\n", err)
			}
			if err := brfp.Flush(); err != nil {
				log.Fatalf("[WriteBrFa] write read seq err: %v\n", err)
			}
			brfp.Close()
			outfp.Close()
			fn = prefix + "_" + strconv.Itoa(fileNum) + "_" + suffix
			fileNum++
			outfp, err1 = os.Create(fn)
			if err1 != nil {
				log.Fatal(err1)
			}
			recordCount = 0
			brfp = cbrotli.NewWriter(outfp, cbrotli.WriterOptions{Quality: 1})
			outbuffp = bufio.NewWriterSize(brfp, 1<<20)
		}
	}

	if err := outbuffp.Flush(); err != nil {
		log.Fatalf("[WriteBrFa] write read seq err: %v\n", err)
	}
	if err := brfp.Flush(); err != nil {
		log.Fatalf("[WriteBrFa] write read seq err: %v\n", err)
	}
	brfp.Close()
	outfp.Close()
}*/
func ReadBrFile(brfn string, fncs chan<- []byte, bufSize int) {

}

func ReadStdin(fncs chan<- []byte, bytesPool *sync.Pool) {
	bp := bytesPool.Get().([]byte)
	bp = bp[:cap(bp)]
	bplen := 0
	size := cap(bp)
	minSize := size * 9 / 10
	buffp := bufio.NewReaderSize(os.Stdin, size)
	//fmt.Printf("[ReadStdin] buffp size:%d\n", size)
	var err error
	var n int
	for err != io.EOF {
		n, err = buffp.Read(bp[bplen:])
		if n > 0 {
			if bplen+n > minSize {
				fncs <- bp[:bplen+n]
				//fmt.Printf("[ReadStdin]read %d bytes, err:%v\n", bplen+n, err)
				bp = bytesPool.Get().([]byte)
				bp = bp[:cap(bp)]
				bplen = 0
			} else {
				bplen += n
			}
		}

		if err != nil {
			if err != io.EOF {
				log.Fatalf("[ReadStdin] found err: %v\n", err)
			}
			break
		}
	}
	fncs <- bp[:bplen]
	close(fncs)
	return
}

func ReadBrFile2(brfn string, bytesPool *sync.Pool, cs chan<- []byte) {
	fp, err1 := os.Open(brfn)
	if err1 != nil {
		log.Fatalf("[ReadBrFile] open file: %v failed..., err: %v\n", brfn, err1)
	}
	defer fp.Close()
	//size := (1 << 26)
	//brfp1 := cbrotli.NewReader(fp1)
	brfp := cbrotli.NewReader(fp)
	defer brfp.Close()
	//buffp := brfp
	//buffp := bufio.NewReaderSize(brfp, 1<<20)
	var err error
	var num int
	for err != io.EOF {
		//buf := make([]byte, size)
		buf := bytesPool.Get().([]byte)
		buf = buf[:cap(buf)]
		//num, err = buffp.Read(buf)
		num, err = brfp.Read(buf)
		if err != nil {
			if err != io.EOF {
				log.Fatalf("[ReadBrFile] read file: %s, found err: %v\n", brfn, err)
			}
		}

		//fmt.Printf("[ReadBrFile]read %d bytes\n", num)
		if num > 0 {
			cs <- buf[:num]
		}
		//time.Sleep(time.Second)
	}

	close(cs)
	//time.Sleep(time.Second * 10)

	return
}

func ReadGzFile(gzfn string, cs chan<- []byte, size int) {
	fp, err1 := os.Open(gzfn)
	if err1 != nil {
		log.Fatalf("[ReadGzFile] open file: %v failed..., err: %v\n", gzfn, err1)
	}
	defer fp.Close()
	//size := (1 << 26)
	//brfp1 := cbrotli.NewReader(fp1)
	gzfp, err2 := gzip.NewReader(fp)
	if err2 != nil {
		log.Fatalf("[ReadGzFile] open file: %v failed..., err: %v\n", gzfn, err2)
	}
	defer gzfp.Close()
	buffp := bufio.NewReaderSize(gzfp, 1<<20)
	var err error
	for err != io.EOF {
		buf := make([]byte, size)
		var num int
		num, err = buffp.Read(buf)
		if err != nil {
			if err != io.EOF {
				log.Fatalf("[ReadGzFile] read file: %v, found err: %v\n", gzfn, err)
			}
		}
		if num > 0 {
			cs <- buf[:num]
		}
	}
	//time.Sleep(time.Second * 10)
	close(cs)

	return
}

func GetNextLine(buf []byte, pos int) ([]byte, int, bool) {
	ok := true
	i := pos
	var bl []byte
	//var bl []byte
	for ; i < len(buf); i++ {
		if buf[i] == '\n' {
			break
		}
	}
	if i >= len(buf) {
		ok = false
	} else {
		/*if len(bl) < i-pos {
			bl = make([]byte, i-pos)
		}
		bl = bl[:i-pos]
		copy(bl, buf[pos:i])*/
		bl = buf[pos:i]
		pos = i + 1
	}
	return bl, pos, ok
}

func GetReadFileRecord(fncs <-chan []byte, recordChan chan<- ReadInfo, format string, bufSize int, bnt bool) int {
	var blockLineNum int
	if format == "fa" {
		blockLineNum = 2
	} else { // format == "fq"
		blockLineNum = 4
	}
	var readNum int
	buf := make([]byte, 0, bufSize*2)
	lastPos := 0
	nb, ok := <-fncs
	if !ok {
		log.Fatalf("[GetReadFileRecord] not found any data\n")
	}
	buf = append(buf, nb...)
	var suc bool
	b := make([][]byte, blockLineNum)
	for {
		//m := make([]byte, 20+252)
		//b[0] = m[:20:20]
		//b[1] = m[20:]
		finished := false
		i := 0
		pos := lastPos
		for i < blockLineNum {
			b[i], pos, suc = GetNextLine(buf, pos)
			//fmt.Printf("[GetReadFileRecord]i:%d pos:%d, b[i]:%s\n", i, pos, string(b[i]))
			if !suc {
				nb, okfncs := <-fncs
				if !okfncs {
					finished = true
					break
				}
				if lastPos < len(buf) {
					//fmt.Printf("[GetReadFileRecord]buf[%d:]:%v\n", lastPos, buf[lastPos:])
					copy(buf[0:], buf[lastPos:len(buf)])
				}
				buf = buf[:len(buf)-lastPos]
				lastPos = 0
				buf = append(buf, nb...)
				//fmt.Printf("[GetReadFileRecord]buf:%v\n", buf)
				break
			}
			i++
		}
		if finished {
			break
		}
		if i != blockLineNum {
			continue
			//log.Fatalf("[GetReadFileRecord] i: %d, found not unbroken record\n", i)
		}
		lastPos = pos
		var ri ReadInfo
		idx := bytes.IndexByte(b[0], '\t')
		var ba []byte
		if idx < 0 {
			ba = b[0]
		} else {
			ba = b[0][:idx]
			if idx < len(b[0])-1 {
				ri.Anotition = b[0][idx+1:]
			}
		}
		//flist := strings.Fields(string(b[0]))
		id, err2 := strconv.Atoi(string(ba[1:]))
		if err2 != nil {
			log.Fatalf("[GetReadFileRecord]  read ID: %v not digits, please convert to digits..., readNum:%d, b[0]:%s\n", ba[1:], readNum, string(b[0]))
		} else {
			ri.ID = uint32(id)
		}
		//for x := 1; x < len(flist); x++ {
		//	ri.Anotition += flist[x] + "\t"
		//}
		//ri.Anotition = strings.Join(flist[1:], "\t")
		ri.Seq = make([]byte, len(b[1]))
		copy(ri.Seq, b[1])
		if bnt {
			ri.Seq = Transform2BntByte2(ri.Seq)
		}

		if format == "fq" {
			ri.Qual = make([]byte, len(b[3]))
			for j, q := range b[3] {
				ri.Qual[j] = q - 33
			}
			//copy(ri.Qual, b[3])
			if len(ri.Seq) != len(ri.Qual) {
				log.Fatalf("[GetReadFileRecord]eID:%d len(ri.Seq):%d != len(ri.Qual):%d\n", ri.ID, len(ri.Seq), len(ri.Qual))
			}
		}
		//fmt.Printf("[GetReadFileRecord]ri: %v\n\tlen(Seq):%d\n", ri, len(ri.Seq))
		readNum++
		recordChan <- ri
	}
	close(recordChan)
	fmt.Printf("[GetReadFileRecord] read Record num is: %v\n", readNum)
	//time.Sleep(time.Second * 60)
	return readNum
}

func GetRecordS(nb []byte, format string) (ri ReadInfoS) {
	lastIdx := 0
	if format == "fa" {
		sz := bytes.IndexAny(nb[lastIdx:], " \t\n")
		if nb[lastIdx+sz-2] == '/' {
			ri.ID = nb[lastIdx+1 : lastIdx+sz-2]
		} else {
			ri.ID = nb[lastIdx+1 : lastIdx+sz]
		}
		lastIdx += sz + 1
		if nb[lastIdx-1] != '\n' {
			sz = bytes.IndexByte(nb[lastIdx:], '\n')
			ri.Anotition = nb[lastIdx : lastIdx+sz]
			lastIdx += sz + 1
		}
		sz = bytes.IndexByte(nb[lastIdx:], '\n')
		ri.Seq = nb[lastIdx : lastIdx+sz]
	} else {
		sz := bytes.IndexAny(nb[lastIdx:], " \t\n")
		if nb[lastIdx+sz-2] == '/' {
			ri.ID = nb[lastIdx+1 : lastIdx+sz-2]
		} else {
			ri.ID = nb[lastIdx+1 : lastIdx+sz]
		}
		lastIdx += sz + 1
		if nb[lastIdx-1] != '\n' {
			sz = bytes.IndexByte(nb[lastIdx:], '\n')
			ri.Anotition = nb[lastIdx : lastIdx+sz]
			lastIdx += sz + 1
		}
		sz = bytes.IndexByte(nb[lastIdx:], '\n')
		ri.Seq = nb[lastIdx : lastIdx+sz]
		lastIdx += sz + 1
		sz = bytes.IndexByte(nb[lastIdx:], '\n')
		lastIdx += sz + 1
		sz = bytes.IndexByte(nb[lastIdx:], '\n')
		ri.Qual = nb[lastIdx : lastIdx+sz]
	}
	//fmt.Printf("[GetRecordS]ri.ID:%s\nri.Seq:%s\nri.Qual:%s\n", string(ri.ID), string(ri.Seq), string(ri.Qual))
	return
}

func GetReadFileRecordS(nb []byte, notFoundNewRecordBuf []byte, riBuf []ReadInfoS, format string, bnt bool) ([]ReadInfoS, int) {
	var blockLineNum int
	var n int
	if format == "fa" {
		blockLineNum = 2
	} else { // format == "fq"
		blockLineNum = 4
	}
	riBuf = riBuf[:0]

	m := bytes.Count(notFoundNewRecordBuf, []byte("\n"))
	idx := 0
	i := 0
	for ; i < blockLineNum-m; i++ {
		if idx >= len(nb) {
			break
			//log.Fatalf("[GetReadFileRecordS] not found proper Record, nb:%s\n", string(nb))
		}
		sz := bytes.IndexByte(nb[idx:], '\n')
		if sz < 0 {
			break
			//log.Fatalf("[GetReadFileRecordS] not found proper Record, nb:%s\n", string(nb))
		}
		idx += sz + 1
	}
	if i != blockLineNum-m {
		return riBuf, n
	}
	notFoundNewRecordBuf = append(notFoundNewRecordBuf, nb[:idx]...)
	//fmt.Printf("[GetReadFileRecordS]notFoundNewRecordBuf:%slen(nb):%d\n", notFoundNewRecordBuf, len(nb))
	riBuf = append(riBuf, GetRecordS(notFoundNewRecordBuf, format))
	lastIdx := idx
	for {
		//m := make([]byte, 20+252)
		//b[0] = m[:20:20]
		//b[1] = m[20:]
		i := 0
		for ; i < blockLineNum; i++ {
			if idx >= len(nb) {
				break
			}
			sz := bytes.IndexByte(nb[idx:], '\n')
			if sz < 0 {
				break
			}
			idx += sz + 1
			//fmt.Printf("[GetReadFileRecord]i:%d, idx:%d\n", i, idx)
		}
		if i != blockLineNum {
			break
		}
		//fmt.Printf("[GetReadFileRecord]i:%d, nb[%d:%d]: %s\n", i, lastIdx, idx, string(nb[lastIdx:idx]))
		riBuf = append(riBuf, GetRecordS(nb[lastIdx:], format))
		lastIdx = idx
	}
	n = lastIdx
	//fmt.Printf("[GetReadFileRecord] read Record num is: %v\n", len(riBuf))
	//time.Sleep(time.Second * 60)
	return riBuf, n
}

func GetReadSeqBucket(fn string, cs chan<- ReadSeqBucket, kmerlen int) (kmerNum int) {
	var processNumReads int
	var rsb ReadSeqBucket
	//var bucketCount int
	format := GetReadsFileFormat(fn)
	bufSize := (1 << 20)
	size := 8000
	fncs1 := make(chan []byte, 10)
	recordChan1 := make(chan ReadInfo, size)
	fmt.Printf("[GetReadSeqBucket] begin processe file: %v...\n", fn)
	//go ReadBrFile(fn, fncs1, bufSize)
	go GetReadFileRecord(fncs1, recordChan1, format, bufSize, true)

	//var count int
	for {
		ri, _ := <-recordChan1
		if ri.ID == 0 {
			break
		}
		if rsb.Count >= ReadSeqSize {
			cs <- rsb
			var nb ReadSeqBucket
			rsb = nb
		}
		if len(ri.Seq) < kmerlen {
			continue
		}
		rsb.ReadBuf[rsb.Count] = ri.Seq
		rsb.Count++
		processNumReads++
		kmerNum += (len(ri.Seq) - (kmerlen - 1))
	}
	if rsb.Count > 0 {
		cs <- rsb
	}
	// send read finish signal
	fmt.Printf("[GetReadSeqBucket] processed reads number is: %d, avg kmer num: %d, finished processed file: %v\n", processNumReads, kmerNum/(processNumReads+1), fn)
	close(cs)
	//time.Sleep(time.Second * 60)
	return
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
	opt := Options{gOpt, 0, false, false, 3}
	tmp, suc := checkArgs(c)
	if suc == false {
		log.Fatalf("[CCF] check Arguments error, opt: %v\n", tmp)
	}
	opt.CFSize = tmp.CFSize
	opt.Correct = tmp.Correct
	opt.Merge = tmp.Merge
	opt.MinKmerFreq = tmp.MinKmerFreq
	fmt.Printf("[CCF] opt: %v\n", opt)
	/*profileFn := opt.Prefix + ".CCF.prof"
	cpuprofilefp, err := os.Create(profileFn)
	if err != nil {
		log.Fatalf("[CCF] open cpuprofile file: %v failed\n", profileFn)
	}
	pprof.StartCPUProfile(cpuprofilefp)
	defer pprof.StopCPUProfile()*/
	cfgInfo, err := ParseCfg(opt.CfgFn, opt.Merge, opt.Correct)
	if err != nil {
		log.Fatalf("[CCF] ParseCfg 'C': %v err :%v\n", opt.CfgFn, err)
	}
	fmt.Println(cfgInfo)

	// make CuckooFilter

	t0 := time.Now()
	cf := cuckoofilter.MakeCuckooFilter(uint64(opt.CFSize), opt.Kmer)
	numCPU := opt.NumCPU
	runtime.GOMAXPROCS(numCPU + 2)
	//bufsize := 100
	wc := make(chan []byte)
	//we := make(chan int)
	//defer close(we)
	var fnArr []string
	for _, lib := range cfgInfo.Libs {
		if lib.AsmFlag != AllState && lib.SeqProfile != 1 {
			continue
		}
		fnArr = append(fnArr, lib.FnName...)
		fnArr = append(fnArr, lib.SingleFn...)
	}

	//fmt.Printf("[CCF] fileName array: %v\n", fnArr)
	totalFileNum := len(fnArr)
	//concurrentNum := 6
	//totalNumT := opt.NumCPU/(concurrentNum+1) + 1
	concurrentNum := opt.Kmer / 15
	if opt.Correct {
		if concurrentNum < 7 {
			concurrentNum = 7
		}
	} else {
		if concurrentNum < 12 {
			concurrentNum = 12
		}
	}
	/*if opt.Kmer < 128 {
		concurrentNum = 6
	} else {
		concurrentNum = 12
	}*/
	totalNumT := opt.NumCPU / concurrentNum
	if totalNumT > totalFileNum {
		totalNumT = totalFileNum
	} else if totalNumT < 1 {
		totalNumT = 1
	}

	processT := make(chan int, totalNumT)
	defer close(processT)
	for i := 0; i < totalNumT; i++ {
		processT <- 0
	}

	// write goroutinue
	go WriteKmer(opt.Prefix, wc, opt.Kmer)
	totalKmerNum := 0
	for _, fn := range fnArr {
		totalKmerNum += <-processT
		//fmt.Printf("[CCF] processing file: %v\n", lib.FnName[i])
		go ConcurrentConstructCF(fn, cf, wc, concurrentNum, opt.Kmer, processT, opt.MinKmerFreq)
	}

	for i := 0; i < totalNumT; i++ {
		totalKmerNum += <-processT
	}
	close(wc)
	for len(wc) > 0 {
		time.Sleep(time.Second)
	}
	time.Sleep(time.Second * 10)
	fmt.Printf("[CCF] totalKmerNum: %v \n", totalKmerNum)
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

func CopyReadInfoS(dst, src *ReadInfoS) {
	if dst.ID == nil {
		dst.ID = append(dst.ID, src.ID...)
		dst.Seq = append(dst.Seq, src.Seq...)
		dst.Qual = append(dst.Qual, src.Qual...)
		dst.Anotition = append(dst.Anotition, src.Anotition...)
	} else {
		dst.ID = append(dst.ID[:0], src.ID...)
		dst.Seq = append(dst.Seq[:0], src.Seq...)
		dst.Qual = append(dst.Qual[:0], src.Qual...)
		dst.Anotition = append(dst.Anotition[:0], src.Anotition...)
	}
}

func ExtractPairEnd(c cli.Command) {
	var err error
	input := c.Flag("input").String()
	IDFile := c.Flag("IDFile").String()
	prefix := c.Flag("prefix").String()
	format := c.Flag("format").String()
	startID := c.Flag("startID").Get().(int)
	//splitRecordNum := c.Flag("SplitRecordNum").Get().(int)
	//splitNum := c.Flag("splitNum").Get().(int)
	//splitNum := 1
	var inputbuffp *bufio.Reader
	/*if inputfa == "" {
		inputbuffp = bufio.NewReader(os.Stdin)
	} else {
		fp, err1 := os.Open(inputfa)
		if err1 != nil {
			log.Fatalf("[GetReadSeqBucket] file: %v open failed, err: %v\n", inputfa, err1)
		}
		defer fp.Close()
		inputgzfp, err1 := gzip.NewReader(fp)
		if err1 != nil {
			log.Fatalf("[GetReadSeqBucket] file: %v not gzip format or other error, err: %v\n", inputfa, err1)
		}
		defer inputgzfp.Close()
		inputbuffp = bufio.NewReaderSize(inputgzfp, 1<<25)
	} */

	var outbuffp1, outbuffp2, outbufSingle *bufio.Writer
	//var brSingle *cbrotli.Writer
	/*var gzfp1, gzfp2 *gzip.Writer
	if format == "gz" {
		outf1, outf2 := prefix+"_1.fa.gz", prefix+"_2.fa.gz"
		outfp1, err1 := os.Create(outf1)
		if err1 != nil {
			log.Fatal(err1)
		}
		defer outfp1.Close()
		gzfp1 = gzip.NewWriter(outfp1)
		defer gzfp1.Close()
		outbuffp1 = bufio.NewWriterSize(gzfp1, 1<<20)

		outfp2, err2 := os.Create(outf2)
		if err2 != nil {
			log.Fatal(err2)
		}
		defer outfp2.Close()
		gzfp2 = gzip.NewWriter(outfp2)
		defer gzfp2.Close()
		outbuffp2 = bufio.NewWriterSize(gzfp2, 1<<20)
	}*/

	var pairNum, singleNum, writePairNum int
	if IDFile != "" {
		// read IDFile
		idfp, err1 := os.Open(IDFile)
		if err1 != nil {
			log.Fatalf("[GetReadSeqBucket] file: %v open failed, err: %v\n", IDFile, err1)
		}
		defer idfp.Close()
		idbuffp := bufio.NewReader(idfp)
		idmap := make(map[int]int)
		err = nil
		for err != io.EOF {
			var t string
			t, err = idbuffp.ReadString('\n')
			if err != nil {
				if err != io.EOF {
					log.Fatalf("[ExtractPairEnd] read file: %v error\n", IDFile)
				}
				break
			}
			id, err2 := strconv.Atoi(t[:len(t)-1])
			if err2 != nil {
				log.Fatalf("[ExtractPairEnd] id: %v set error\n", t)
			}
			idmap[id]++
		}
		paired := 8
		hasoutput := 4
		var errorNum int
		for k, v := range idmap {
			if v == 2 {
				idmap[k] = paired
				pairNum += 2
			} else if v == 1 {
				singleNum++
			} else {
				errorNum++
			}
		}
		fmt.Printf("[ExtractPairEnd] pairNum: %v, singleNum: %v, errorNum: %v\n", pairNum, singleNum, errorNum)
		if errorNum > 0 {
			log.Fatalf("[ExtractPairEnd]  errorNum: %v > 0\n", errorNum)
		}
		//var buffp1, buffp2 *bufio.Writer
		//var nbrfp1, nbrfp2 *cbrotli.Writer
		//buffp1, buffp2 = outbuffp1, outbuffp2
		outputNum := 0
		rd := make([]string, 2)
		err = nil
		for err != io.EOF {
			for i := 0; i < 2; i++ {
				rd[i], err = inputbuffp.ReadString('\n')
				if err != nil {
					break
				}
			}
			if err != nil && err != io.EOF {
				log.Fatalf("[ExtractPairEnd] err: %v set error\n", err)
			}
			if err == io.EOF {
				break
			}
			//fmt.Printf("[ExtractPairEnd] id: %v", rd[0])
			fd := strings.Fields(rd[0])
			id, err2 := strconv.Atoi(fd[0][1:])
			if err2 != nil {
				log.Fatalf("[ExtractPairEnd] id: %v set error\n", fd[0][1:])
			}
			v, ok := idmap[id]
			if ok {
				rd[0] = rd[0][:1] + strconv.Itoa(id+startID) + "\n"
				if v == paired {
					for i := 0; i < 2; i++ {
						outbuffp1.WriteString(rd[i])
					}
					idmap[id] = hasoutput
					outputNum++
				} else if v == hasoutput {
					for i := 0; i < 2; i++ {
						outbuffp2.WriteString(rd[i])
					}
					idmap[id] = 0
					writePairNum += 2
					outputNum++
				} else {
					for i := 0; i < 2; i++ {
						outbufSingle.WriteString(rd[i])
					}
					idmap[id] = 0
					outputNum++
				}
			} else {
				log.Fatalf("[ExtractPairEnd] read id: %v not include id file: %v\n", id, IDFile)
			}
		}
		// check code
		if writePairNum != pairNum {
			log.Fatalf("[ExtractPairEnd] writePairNum: %v != pairNum: %v\n", writePairNum, pairNum)
		}
	} else {
		bufSize := (1 << 19) * 5
		size := 5000
		nextSize := 600
		//recordChan := make(chan ReadInfo, size)

		if false {
			profileFn := prefix + ".ExtractPairEnd.prof"
			cpuprofilefp, err := os.Create(profileFn)
			if err != nil {
				log.Fatalf("[ExtractPairEnd] open cpuprofile file: %v failed\n", profileFn)
			}
			pprof.StartCPUProfile(cpuprofilefp)
			defer pprof.StopCPUProfile()
		}

		//fArr := strings.Split(inputfa, ".")
		//if len(fArr) < 2 {
		//	log.Fatalf("[ExtractPairEnd]input file:%v can identify file compress format:%v\n", inputfa, fArr)
		//}
		rbytesPool := sync.Pool{New: func() interface{} {
			bytes := make([]byte, bufSize)
			return bytes
		}}

		wbytesPool := sync.Pool{New: func() interface{} {
			bytes := make([]byte, bufSize+nextSize)
			return bytes
		}}
		//riPool := sync.Pool{New: func() interface{} {
		//	readInfos := make([]ReadInfo, size)
		//	return readInfos
		//}}

		fncs := make(chan []byte, 1)
		wc := make(chan []byte, 1)
		var f string
		if input == "/dev/stdin" {
			f = "fq"
			go ReadStdin(fncs, &rbytesPool)
		} else {
			f = GetReadsFileFormat(input)
			go ReadBrFile2(input, &rbytesPool, fncs)
		}
		riBuf := make([]ReadInfoS, size)
		notFoundNewLineBuf := make([]byte, 0, nextSize)
		var lastRI ReadInfoS
		var wfn string
		var n int
		if format == "fa" {
			wfn = prefix + ".fa.br"
		} else {
			wfn = prefix + ".fq.br"
		}
		var readBlockNum, writeBlockNum int
		wfinish := make(chan bool)
		go WriteBr(wfn, wc, &wbytesPool, wfinish)
		for {
			if len(fncs) == 0 {
				readBlockNum++
			}
			cs, okcs := <-fncs
			if !okcs {
				break
			}
			//cs = append(cs, notFoundNewLineBuf...)
			riBuf, n = GetReadFileRecordS(cs, notFoundNewLineBuf, riBuf, f, false)
			if n == 0 {
				notFoundNewLineBuf = append(notFoundNewLineBuf, cs...)
				continue
			}
			var last *ReadInfoS
			idx := 0
			//fmt.Printf("[ExtractPairEnd] lastRI.ID:%s, riBuf[0].ID:%s\n", string(lastRI.ID), string(riBuf[0].ID))
			if len(lastRI.ID) == 0 {
				last = &riBuf[0]
				idx = 1
			} else {
				last = &lastRI
			}
			//fmt.Printf("[ExtractPairEnd]last.ID:%s\n", string(last.ID))
			wbp := wbytesPool.Get().([]byte)
			wbp = wbp[:0]
			for j := idx; j < len(riBuf); j++ {
				if last == nil {
					last = &riBuf[j]
				} else {
					//fmt.Printf("[ExtractPairEnd]last.ID:%s\triBuf[%d]:%s\n", string(last.ID), j, string(riBuf[j].ID))
					if utils.BytesEqual(last.ID, riBuf[j].ID) {
						if format == "fq" {
							wbp = append(wbp, '@')
							wbp = append(wbp, strconv.Itoa(startID)...)
							wbp = append(wbp, '\n')
							wbp = append(wbp, last.Seq...)
							wbp = append(wbp, "\n+\n"...)
							wbp = append(wbp, last.Qual...)

							wbp = append(wbp, "\n@"...)
							wbp = append(wbp, strconv.Itoa(startID)...)
							wbp = append(wbp, '\n')
							wbp = append(wbp, riBuf[j].Seq...)
							wbp = append(wbp, "\n+\n"...)
							wbp = append(wbp, riBuf[j].Qual...)
							wbp = append(wbp, '\n')
						} else {
							wbp = append(wbp, '>')
							wbp = append(wbp, strconv.Itoa(startID)...)
							wbp = append(wbp, '\n')
							wbp = append(wbp, last.Seq...)

							wbp = append(wbp, "\n>"...)
							wbp = append(wbp, strconv.Itoa(startID)...)
							wbp = append(wbp, '\n')
							wbp = append(wbp, riBuf[j].Seq...)
							wbp = append(wbp, '\n')
						}
						startID++
						pairNum++
						last = nil
					} else {
						// write to single file
						//fmt.Printf("[ExtractPairEnd]last.ID:%s\triBuf[%d]:%s\n", string(last.ID), j, string(riBuf[j].ID))
						singleNum++
						last = &riBuf[j]
					}
				}
			}
			if last != nil {
				CopyReadInfoS(&lastRI, last)
			} else {
				lastRI.ID = lastRI.ID[:0]
			}
			notFoundNewLineBuf = append(notFoundNewLineBuf[:0], cs[n:]...)
			rbytesPool.Put(cs)
			if len(wc) == 1 {
				writeBlockNum++
			}
			wc <- wbp
		}
		close(wc)
		<-wfinish
		fmt.Printf("[ExtractPairEnd] pairNum:%d, single:%d, readBlockNum:%d, writeBlockNum:%d\n", pairNum, singleNum, readBlockNum, writeBlockNum)
	}
	//readMap := make(map[int64]ReadInfo)
	//outf1, outf2 := prefix+"_1.fa.br", prefix+"_2.fa.br"
	//go WriteSplitBrFa(prefix, ".single.fa.br", wcS, splitRecordNum)
	//go WriteBrFa(outf1, wc1)
	//go WriteBrFa(outf2, wc2)
	// read input fa

	//if format == "gz" {
	//	if err := gzfp1.Flush(); err != nil {
	//		log.Fatalf("[ExtractPairEnd] write read1 seq err: %v\n", err)
	//	}
	//	if err := gzfp2.Flush(); err != nil {
	//		log.Fatalf("[ExtractPairEnd] write read2 seq err: %v\n", err)
	//	}
	//}
	return
}

func GetMeanQual(ri ReadInfo) (mq int) {
	tq := 0
	//qArr = make([]int, (len(ri.Qual)+99)/100)
	//wq := 0
	for _, q := range ri.Qual {
		tq += int(q) - 33
	}
	if len(ri.Qual) > 0 {
		mq = tq / len(ri.Qual)
	}
	return
}

func FilterLong(c cli.Command) {
	inputfq := c.Flag("input").String()
	//output := c.Flag("output").String()
	startID := c.Flag("startID").Get().(int)
	minLen := c.Flag("minLen").Get().(int)
	minMeanQuailty := c.Flag("minMeanQuality").Get().(int)
	bufSize := (1 << 20)
	size := 1000
	fncs := make(chan []byte, 10)
	recordChan := make(chan ReadInfo, size)
	f := GetReadsFileFormat(inputfq)

	go ReadGzFile(inputfq, fncs, bufSize)
	go GetReadFileRecord(fncs, recordChan, f, bufSize, false)

	wc := make(chan ReadInfo, size)
	//outf1, outf2 := prefix+"_1.fa.br", prefix+"_2.fa.br"
	//go WriteBrFa(output, wc)
	shortFilterNum := 0
	lowQualFilterNum := 0
	filterSeqLen := 0
	passNum := 0
	maxQual := 50
	qualArr := make([]int, maxQual)
	// read input fa
	for {
		ri, ok := <-recordChan
		if !ok {
			break
		}

		if len(ri.Seq) < minLen {
			shortFilterNum++
			filterSeqLen += len(ri.Seq)
			continue
		}
		mq := GetMeanQual(ri)
		if mq > maxQual-1 {
			mq = maxQual - 1
		}
		qualArr[mq]++
		if mq < minMeanQuailty {
			lowQualFilterNum++
			filterSeqLen += len(ri.Seq)
			continue
		}
		ri.ID = uint32(startID + passNum)
		//ri.Anotition = "MQ:" + strconv.Itoa(mq)
		/*for _, q := range arr {
			ri.Anotition += strconv.Itoa(q) + ":"
		}*/
		wc <- ri
		passNum++
	}
	//time.Sleep(time.Second * 10)
	close(wc)
	time.Sleep(time.Second * 3)

	qs, rn := 0, 0
	for i, num := range qualArr {
		fmt.Printf("[FilterLong]qual: %d readNum:%d\n", i, num)
		if i >= minMeanQuailty {
			qs += i * num
			rn += num
		}
	}
	fmt.Printf("[FilterLong] shortFilterNum:%d lowQualFilterNum:%d filterSeqLen:%d passNum:%d avg qual: %d\n", shortFilterNum, lowQualFilterNum, filterSeqLen, passNum, qs/(rn+1))

	return
}
