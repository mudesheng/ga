package main

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

	"github.com/google/brotli/go/cbrotli"
	"github.com/jwaldrip/odin/cli"
	"github.com/klauspost/compress/zstd"
)

const (
	AllState          = 1
	ScaffState        = 2
	GapState          = 3
	EndString         = "end"
	ReadSeqSize       = 1000
	UniqueKmerFileNum = 2
	ReadLen           = 250
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

type IDSeq struct {
	ID  []byte
	Seq []byte
}

type ReadInfoS struct {
	ID        []byte
	Seq       []byte
	Qual      []byte
	Anotition []byte
}

type ReadInfo struct {
	ID        int
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
			if len(fs) < 2 || fs[len(fs)-1] != "zst" || (fs[len(fs)-2] != "fa" && fs[len(fs)-2] != "fasta" && fs[len(fs)-2] != "fq" && fs[len(fs)-2] != "fastq") {
				log.Fatalf("[ParseCfg] fn: %v, must used suffix *.[fasta|fa|fq|fastq].zst\n", fields[2])
			}
			if correct {
				libInfo.FnName = append(libInfo.FnName, fields[2])
			} else if merge {
				if libInfo.Merged < 2 {
					zstfn := strings.Join(fs[:len(fs)-2], ".") + ".Correct." + "fa." + fs[len(fs)-1]
					libInfo.FnName = append(libInfo.FnName, zstfn)
				} else {
					zstfn := strings.Join(fs[:len(fs)-2], ".") + ".Merged." + fs[len(fs)-2] + "." + fs[len(fs)-1]
					libInfo.FnName = append(libInfo.FnName, zstfn)
				}
			} else {
				fs := strings.Split(fields[2], ".")
				if len(fs) < 2 || fs[len(fs)-1] != "zst" || (fs[len(fs)-2] != "fa" && fs[len(fs)-2] != "fasta" && fs[len(fs)-2] != "fq" && fs[len(fs)-2] != "fastq") {
					log.Fatalf("[ParseCfg] fn: %v, must used suffix *.[fasta|fa|fq|fastq].zst\n", fields[2])
				}
				zstfn := strings.Join(fs[:len(fs)-2], ".") + ".Merged.fa." + fs[len(fs)-1]
				libInfo.FnName = append(libInfo.FnName, zstfn)
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
		base := tmp[i/NumBaseInUint64] & BaseMask
		extRB[i] = uint8(base)
		tmp[i/NumBaseInUint64] >>= NumBitsInBase
	}

	return extRB
}

func ExtendKmerBnt2Byte2(kb []uint64, seq []byte) {
	if len(kb) == 0 || len(seq) > len(kb)*32 || len(seq) == 0 {
		log.Fatalf("[ExtendKmerBnt2Byte]kb:%v seq:%v\n", kb, seq)
	}
	//fmt.Printf("[ExtendReadBnt2Byte]kb:%v seq:%v\n", kb, seq)
	i := len(seq) - 1
	for ; i >= 0; i-- {
		//fmt.Printf("[ExtendReadBnt2Byte]i:%d\n", i)
		base := kb[i/NumBaseInUint64] & BaseMask
		seq[i] = uint8(base)
		kb[i/NumBaseInUint64] >>= NumBitsInBase
	}
}

func GetReadBntKmer(seq []byte, startPos, kmerlen int) (sb []uint64) {
	//kLen = kmerlen
	sb = make([]uint64, (kmerlen+NumBaseInUint64-1)/NumBaseInUint64)
	//fmt.Printf("[GetReadBntKmer] len(seq): %d, seq: %v\n\tstartPos: %d, kmerlen: %d\n", len(seq), seq, startPos, kmerlen)

	for i := 0; i < kmerlen; i++ {
		sb[i/NumBaseInUint64] <<= NumBitsInBase
		sb[i/NumBaseInUint64] |= uint64(seq[i+startPos])
	}

	return
}

func NoAllocGetReadBntKmer(seq []byte, startPos, kmerlen int, kb KmerBnt) KmerBnt {
	kb.Len = kmerlen
	//kSeq = make([]uint64, (kLen+NumBaseInUint64-1)/NumBaseInUint64)
	//fmt.Printf("[GetReadBntKmer] len(seq): %d, seq: %v\n\tstartPos: %d, kmerlen: %d\n", len(seq), seq, startPos, kmerlen)
	copy(kb.Seq, seq[startPos:startPos+kmerlen])

	return kb
}

func ReverseCompletBnt(kb []uint64, klen int) (rb []uint64) {
	rb = make([]uint64, len(kb))
	tmp := make([]uint64, len(kb))
	copy(tmp, kb)
	for i := klen - 1; i >= 0; i-- {
		base := tmp[i/NumBaseInUint64] & BaseMask
		tmp[i/NumBaseInUint64] >>= NumBitsInBase
		rb[(klen-i-1)/NumBaseInUint64] <<= NumBitsInBase
		rb[(klen-i-1)/NumBaseInUint64] |= uint64(BntRev[base])
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

func GetReverseCompletUint16(bs uint16, slen int) (rs uint16) {
	for i := 0; i < slen; i++ {
		rs <<= NumBitsInBase
		rs |= uint16(BntRev[bs&0x3])
		bs >>= NumBitsInBase
	}
	return rs
}

func GetReverseComplet(kb KmerBnt) KmerBnt {
	var rb KmerBnt
	rb.Len = kb.Len
	rb.Seq = make([]byte, rb.Len)
	for i := 0; i < kb.Len; i++ {
		rb.Seq[kb.Len-1-i] = BntRev[kb.Seq[i]]
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
		rb.Seq[kb.Len-1-i] = BntRev[kb.Seq[i]]
	}
	return rb
}

func GetReverseCompletBytes(kb, rb []byte) []byte {
	kl := len(kb)
	if kl > cap(rb) {
		log.Fatalf("[GetReverseCompletBytes] cap(rb):%d < len(kb):%d\n", cap(rb), kl)
	}
	rb = rb[:kl]
	for i := 0; i < kl; i++ {
		rb[kl-1-i] = BntRev[kb[i]]
	}
	return rb
}

func ReverseCompletIncoming(nd *DBGNode) {
	nd.EdgeIDIncoming[0], nd.EdgeIDOutcoming[3] = nd.EdgeIDOutcoming[3], nd.EdgeIDIncoming[0]
	nd.EdgeIDIncoming[1], nd.EdgeIDOutcoming[2] = nd.EdgeIDOutcoming[2], nd.EdgeIDIncoming[1]
	nd.EdgeIDIncoming[2], nd.EdgeIDOutcoming[1] = nd.EdgeIDOutcoming[1], nd.EdgeIDIncoming[2]
	nd.EdgeIDIncoming[3], nd.EdgeIDOutcoming[0] = nd.EdgeIDOutcoming[0], nd.EdgeIDIncoming[3]
}

/*func DeleteLastBaseKmer(kb KmerBnt) (dk KmerBnt, base uint64) {
	nLen := (kb.Len + NumBaseInUint64 - 1) / NumBaseInUint64
	dk.Len = kb.Len - 1
	base = kb.Seq[nLen-1] & BaseMask
	if nLen > (kb.Len-1+NumBaseInUint64-1)/NumBaseInUint64 {
		dk.Seq = make([]uint64, nLen-1)
		copy(dk.Seq, kb.Seq)
	} else {
		dk.Seq = make([]uint64, nLen)
		copy(dk.Seq, kb.Seq)
		dk.Seq[nLen-1] >>= NumBitsInBase
	}
	return
}*/

/*func DeleteFirstBaseKmer(kb KmerBnt) (dk KmerBnt, base uint64) {
	nLen := (kb.Len + NumBaseInUint64 - 1) / NumBaseInUint64
	dk.Len = kb.Len - 1
	seqLen := (dk.Len + NumBaseInUint64 - 1) / NumBaseInUint64
	dk.Seq = make([]uint64, seqLen)
	for i := nLen - 1; i >= 0; i-- {
		if i == nLen-1 {
			offset := (uint64(kb.Len)%NumBaseInUint64 - 1) * NumBitsInBase
			base = (kb.Seq[i] >> offset) & BaseMask
			mask := uint64(1<<offset) - 1
			if seqLen == nLen {
				dk.Seq[i] = kb.Seq[i] & mask
			}
		} else {
			tb := (kb.Seq[i] >> 62) & BaseMask
			dk.Seq[i] = (kb.Seq[i] << NumBitsInBase) | base
			base = tb
		}
	}

	return
}*/

/*func GetNextKmer(kb KmerBnt, base uint64, kmerlen int) (next KmerBnt) {
	nLen := (kmerlen + NumBaseInUint64 - 1) / NumBaseInUint64
	next.Len = kmerlen
	next.Seq = make([]uint64, nLen)
	if kb.Len == kmerlen-1 {
		copy(next.Seq, kb.Seq)
		next.Seq[nLen-1] <<= NumBitsInBase
		next.Seq[nLen-1] |= base
		return
	} else if kb.Len != kmerlen {
		log.Fatalf("[GetNextKmer] length of kb set error, kmerlen:%d,   kb length: %v\n", kmerlen, kb.Len)
	}
	for i := nLen - 1; i >= 0; i-- {
		//fmt.Printf("[GetNextKmer] base: %v, kb.Seq[%d]: %b\n", base, i, kb.Seq[i])
		if i == nLen-1 {
			offset := (uint64(kmerlen)%NumBaseInUint64 - 1) * NumBitsInBase
			tb := (kb.Seq[i] >> offset) & BaseMask
			mask := uint64(1<<offset) - 1
			next.Seq[i] = kb.Seq[i] & mask
			next.Seq[i] = (next.Seq[i] << NumBitsInBase) | base
			base = tb
		} else {
			tb := (kb.Seq[i] >> 62) & BaseMask
			next.Seq[i] = (kb.Seq[i] << NumBitsInBase) | base
			base = tb
		}
		//fmt.Printf("[GetNextKmer] base: %v, next.Seq[%d]: %b\n", base, i, next.Seq[i])
	}
	return
}
func NoAllocGetNextKmer(kb1, kb2 KmerBnt, base uint64, kmerlen int) KmerBnt {
	nLen := (kmerlen + NumBaseInUint64 - 1) / NumBaseInUint64
	kb2.Len = kmerlen
	if kb1.Len == kmerlen-1 {
		//kb2.Len = kmerlen
		copy(kb2.Seq, kb1.Seq)
		kb2.Seq[nLen-1] <<= NumBitsInBase
		kb2.Seq[nLen-1] |= base
		return kb2
	} else if kb1.Len != kmerlen {
		log.Fatalf("[GetNextKmer] length of kb set error, kmerlen:%d,   kb1 length: %v\n", kmerlen, kb1.Len)
	}
	for i := nLen - 1; i >= 0; i-- {
		//fmt.Printf("[GetNextKmer] base: %v, kb.Seq[%d]: %b\n", base, i, kb.Seq[i])
		if i == nLen-1 {
			offset := (uint64(kmerlen)%NumBaseInUint64 - 1) * NumBitsInBase
			tb := (kb1.Seq[i] >> offset) & BaseMask
			mask := uint64(1<<offset) - 1
			kb2.Seq[i] = kb1.Seq[i] & mask
			kb2.Seq[i] = (kb2.Seq[i] << NumBitsInBase) | base
			base = tb
		} else {
			tb := (kb1.Seq[i] >> 62) & BaseMask
			kb2.Seq[i] = (kb1.Seq[i] << NumBitsInBase) | base
			base = tb
		}
		//fmt.Printf("[GetNextKmer] base: %v, next.Seq[%d]: %b\n", base, i, next.Seq[i])
	}
	return kb2
}

func GetPreviousKmer(kb KmerBnt, base uint64, kmerlen int) (pkb KmerBnt) {
	nLen := (kmerlen + NumBaseInUint64 - 1) / NumBaseInUint64
	pkb.Len = kmerlen
	pkb.Seq = make([]uint64, nLen)
	if kb.Len != kmerlen-1 && kb.Len != kmerlen {
		log.Fatalf("[GetPreviousKmer] length of kb set error, kmerlen:%d,   kb length: %v\n", kmerlen, kb.Len)
	}
	for i := 0; i < nLen; i++ {
		if i == nLen-1 {
			offset := (uint64(kmerlen)%NumBaseInUint64 - 1) * NumBitsInBase
			if kb.Len == kmerlen {
				pkb.Seq[i] = kb.Seq[i] >> NumBitsInBase
			} else {
				pkb.Seq[i] = kb.Seq[i]
			}
			pkb.Seq[i] |= (base << offset)
		} else {
			tb := kb.Seq[i] & BaseMask
			pkb.Seq[i] = kb.Seq[i] >> NumBitsInBase
			pkb.Seq[i] |= (base << 62)
			base = tb
		}
	}
	return
}

func NoAllocGetPreviousKmer(rb1, rb2 KmerBnt, base uint64, kmerlen int) KmerBnt {
	nLen := (kmerlen + NumBaseInUint64 - 1) / NumBaseInUint64
	rb2.Len = kmerlen
	if rb1.Len != kmerlen-1 && rb1.Len != kmerlen {
		log.Fatalf("[GetPreviousKmer] length of rb set error, kmerlen:%d,   rb1 length: %v\n", kmerlen, rb1.Len)
	}
	for i := 0; i < nLen; i++ {
		if i == nLen-1 {
			offset := (uint64(kmerlen)%NumBaseInUint64 - 1) * NumBitsInBase
			if rb1.Len == kmerlen {
				rb2.Seq[i] = rb1.Seq[i] >> NumBitsInBase
			} else {
				rb2.Seq[i] = rb1.Seq[i]
			}
			rb2.Seq[i] |= (base << offset)
		} else {
			tb := rb1.Seq[i] & BaseMask
			rb2.Seq[i] = rb1.Seq[i] >> NumBitsInBase
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

func ParaConstructCF(cf *CuckooFilter, rbytesPool *sync.Pool, seqArrPool *sync.Pool, seqPoolChan <-chan SeqPool, wc chan<- []byte, wbytesPool *sync.Pool, MinKmerFreq int, kmerNumC chan<- int) {
	ws := wbytesPool.Get().([]byte)
	ws = ws[:0]
	var kmerNum int
	rSeq := make([]byte, 800)
	var seq, kb, rb, min []byte
	var lenS int

	for {
		sp, ok := <-seqPoolChan
		if !ok {
			if len(ws) > 0 {
				wc <- ws
			}
			break
		}

		for i := 0; i < len(sp.SeqArr); i++ {
			seq = sp.SeqArr[i]
			lenS = len(seq)
			if lenS < cf.Kmerlen+1 {
				continue
			}
			//kb1 = NoAllocGetReadBntKmer(rBntSeq, 0, cf.Kmerlen-1, kb1)
			rSeq = GetReverseCompletBytes(seq, rSeq)
			/*ks := GetReadBntKmer(rBntSeq, 0, cf.Kmerlen-1)
			rs := ReverseComplet(ks)
			if !reflect.DeepEqual(kb1.Seq, ks.Seq) {
				log.Fatalf("[ParaConstructCF]kb1.Seq: %v != ks.Seq: %v\n", kb1, ks)
			}*/
			kmerNum += lenS - (cf.Kmerlen - 1)
			for j := 0; j < lenS-(cf.Kmerlen-1); j++ {
				kb = seq[j : j+cf.Kmerlen]
				rb = rSeq[lenS-cf.Kmerlen-j : lenS-j]
				//kb2 = NoAllocGetNextKmer(kb1, kb2, uint64(rBntSeq[j]), cf.Kmerlen)
				//rb2 = NoAllocGetPreviousKmer(rb1, rb2, uint64(BntRev[rBntSeq[j]]), cf.Kmerlen)
				/*ks = GetNextKmer(ks, uint64(rBntSeq[j]), cf.Kmerlen)
					rs = GetPreviousKmer(rs, uint64(BntRev[rBntSeq[j]]), cf.Kmerlen)
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
				min = kb
				if BiggerThan(min, rb) {
					min = rb
				}
				//cpbMin := CompressToBnt(min)
				count, suc := cf.Insert(min)
				if !suc {
					log.Fatal("[ParaConstructCF] Insert to the CuckooFilter false")
				}
				//fmt.Printf("retrun count : %d\n", count)
				//fmt.Printf("count set: %d\n", cf.GetCount(ks.Seq))
				if count == uint16(MinKmerFreq)-1 {
					if len(ws)+cf.Kmerlen > cap(ws) {
						wc <- ws
						ws = wbytesPool.Get().([]byte)
						ws = ws[:0]
					}
					//fmt.Printf("[ParaConstructCF] wrsb.count: %d\n", wrsb.Count)
					ws = append(ws, min...)
				}
			}
		}
		// put pool
		rbytesPool.Put(sp.Cs)
		//notFoundNewRecordPool.Put(sp.NotFoundRecordBytes)
		seqArrPool.Put(sp.SeqArr)
	}
	kmerNumC <- kmerNum
}

type SeqPool struct {
	SeqArr              [][]byte
	Cs                  []byte
	NotFoundRecordBytes []byte
}

type IDSeqPool struct {
	IDSeqArr            []IDSeq
	Cs                  []byte
	NotFoundRecordBytes []byte
}

type RIPool struct {
	RIArr               []ReadInfo
	Cs                  []byte
	NotFoundRecordBytes []byte
	MergeSeqArr         [][]byte
}

type RISPool struct {
	RIArr               []ReadInfoS
	Cs                  []byte
	NotFoundRecordBytes []byte
	MergeSeqArr         [][]byte
}

func ConcurrentConstructCF(fn string, cf *CuckooFilter, wc chan<- []byte, writeBytesPool *sync.Pool, concurrentNum int, kmerlen int, MinKmerFreq int) (kmerNum int) {
	cs := make(chan []byte, concurrentNum)
	kmerNumC := make(chan int)
	defer close(kmerNumC)
	bufSize := WindowSize

	rbytesPool := sync.Pool{New: func() interface{} {
		bytes := make([]byte, bufSize)
		return bytes
	}}

	seqArrPool := sync.Pool{New: func() interface{} {
		arr := make([][]byte, bufSize/ReadLen)
		return arr
	}}

	seqPoolChan := make(chan SeqPool, concurrentNum)
	//defer close(seqPoolChan)

	//go ReadBrFile2(fn, &rbytesPool, cs)
	go ReadZstdFile(fn, &rbytesPool, cs)
	for i := 0; i < concurrentNum; i++ {
		go ParaConstructCF(cf, &rbytesPool, &seqArrPool, seqPoolChan, wc, writeBytesPool, MinKmerFreq, kmerNumC)
	}

	GetReadSeqBucket(fn, cs, &seqArrPool, seqPoolChan)

	for i := 0; i < concurrentNum; i++ {
		kmerNum += <-kmerNumC
	}
	return
}

// write kmer
func WriteZstd(wrfn string, wc <-chan []byte, writeBytesPool *sync.Pool, wn chan<- int) {
	outfp, err := os.Create(wrfn)
	if err != nil {
		log.Fatal(err)
	}
	defer outfp.Close()

	zw, err1 := zstd.NewWriter(outfp, zstd.WithEncoderCRC(false), zstd.WithEncoderConcurrency(1), zstd.WithEncoderLevel(1))
	if err1 != nil {
		log.Fatalf("[WriteZstd]open write file:%s err: %v\n", wrfn, err)
	}
	fmt.Printf("zw:%v\n", zw)
	var emptyCount, fullCount int
	//brfp := cbrotli.NewWriter(outfp, cbrotli.WriterOptions{Quality: 1, LGWin: 15})
	//brfp := cbrotli.NewWriter(outfp, cbrotli.WriterOptions{LGWin: 1 << 10})
	defer zw.Close()
	// bufwriter := bufio.NewWriter(gzwriter)
	// defer outfp.Close()
	// KBntByteNum := (Kmerlen + NumBaseInByte - 1) / NumBaseInByte
	writeCount := 0
	for {
		if len(wc) == 0 {
			emptyCount++
		} else if len(wc) == cap(wc) {
			fullCount++
		}

		bs, ok := <-wc
		if !ok {
			break
		}
		_, err := zw.Write(bs)
		if err != nil {
			log.Fatalf("[WriteZstd]open write file:%s err: %v\n", wrfn, err)
		}
		writeCount += len(bs)
		writeBytesPool.Put(bs)
	}
	//if err = zw.Flush(); err != nil {
	//	log.Fatalf("[WriteZstd]flush write err: %v\n", err)
	//}
	wn <- writeCount
	fmt.Printf("[WriteZstd] file:%s write total size is %d\n", wrfn, writeCount)
	fmt.Printf("[WriteZstd]emptyCount:%d fullCount:%d\n", emptyCount, fullCount)
	return
}

func WriteUniqKmerseqInfo(kmerInfofn string, writeKmerNum int) (err error) {
	outfp, err := os.Create(kmerInfofn)
	if err != nil {
		log.Fatal(err)
	}
	defer outfp.Close()
	_, err = outfp.WriteString(fmt.Sprintf("kmerNum:%d\n", writeKmerNum))
	if err != nil {
		return err
	}
	return
}

func ReadUniqKmerseqInfo(kmerInfofn string) (writeKmerNum int, err error) {
	fp, err := os.Open(kmerInfofn)
	if err != nil {
		log.Fatal(err)
	}
	defer fp.Close()
	//var bs []byte
	//bs, err = ioutil.ReadAll(fp)
	_, err = fmt.Fscanf(fp, "kmerNum:%d\n", &writeKmerNum)
	if err != nil {
		log.Fatalf("[ReadUniqKmerseqInfo] scan kmerNum error:%v\n", err)
	}
	//fs := strings.Split(string(bs), ":\n")
	//writeKmerNum, err = strconv.Atoi(fs[1])
	return
}

/*func Trans2Byte(s string) (rb ReadBnt) {
	rb.Length = len(s)
	rb.Seq = make([]byte, (rb.Length+NumBaseInByte-1)/NumBaseInByte)
	for i := 0; i < rb.Length; i++ {
		b := Base2Bnt[s[i]]
		if b > 3 {
			fmt.Printf("[Trans2Byte]found input sequence base '%c' not belong 'ACTG/actg', please check\n", s[i])
			log.Fatal("error found")
		}
		rb.Seq[i/NumBaseInByte] <<= NumBitsInBase
		rb.Seq[i/NumBaseInByte] |= b
	}

	return rb
}*/

func Transform2BntByte2(ks []byte) []byte {
	//ls := make([]byte, len(ks))
	for i, b := range ks {
		if b == 'N' {
			return ks[:i]
		}
		//ls = append(ls, alphabet.Letter(BitNtCharUp[b]))
		ks[i] = Base2Bnt[b]
	}

	return ks
}

type optionsCCF struct {
	ArgsOpt
	CFSize      int
	Correct     bool
	Merge       bool
	MinKmerFreq int
}

func checkArgsCCF(c cli.Command) (opt optionsCCF, suc bool) {
	gOpt, suc := CheckGlobalArgs(c.Parent())
	if suc == false {
		log.Fatalf("[checkArgsCCF] check global Arguments error, opt: %v\n", gOpt)
	}
	opt.ArgsOpt = gOpt
	opt.CFSize = c.Flag("S").Get().(int)
	if opt.CFSize < 1024*1024 {
		log.Fatalf("[checkArgsCCF]the argument 'S':%d must bigger than 1024 * 1024\n", opt.CFSize)
	}
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
	if minKmerFreq < 1 || minKmerFreq > MaxC {
		log.Fatalf("[checkArgs]the argument 'MinKmerFreq': %v must [%v ~ %v]\n", minKmerFreq, 1, MaxC)
	}
	opt.MinKmerFreq = minKmerFreq
	suc = true
	return opt, suc
}

func GetReadsFileFormat(fn string) (format string) {
	sfn := strings.Split(fn, ".")
	if len(sfn) < 3 {
		log.Fatalf("[GetReadsFileFormat] reads file: %v need suffix end with '*.fa.zst | *.fasta.zst | *.fq.zst | *.fastq.zst'\n", fn)
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

func WriteGz(fn string, wc <-chan []byte, wbytePool *sync.Pool, wfinish chan<- bool) {
	outfp, err1 := os.Create(fn)
	if err1 != nil {
		log.Fatal(err1)
	}
	defer outfp.Close()
	gzfp := gzip.NewWriter(outfp)
	defer gzfp.Close()
	//outbuffp := bufio.NewWriterSize(brfp, 1<<20)
	for {
		wb, ok := <-wc
		if !ok {
			break
		}
		gzfp.Write(wb)
		wbytePool.Put(wb)
	}

	if err := gzfp.Flush(); err != nil {
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
				log.Fatalf("[ReadBrFile2] read file: %s, found err: %v\n", brfn, err)
			}
		}

		fmt.Printf("[ReadBrFile2]read %d bytes\n", num)
		if num > 0 {
			cs <- buf[:num]
		}
		//time.Sleep(time.Second)
	}

	close(cs)
	//time.Sleep(time.Second * 10)

	return
}

func ReadGzFile2(gzfn string, bytesPool *sync.Pool, cs chan<- []byte) {
	fp, err1 := os.Open(gzfn)
	if err1 != nil {
		log.Fatalf("[ReadBrFile] open file: %v failed..., err: %v\n", gzfn, err1)
	}
	defer fp.Close()
	//size := (1 << 26)
	//brfp1 := cbrotli.NewReader(fp1)
	gzfp, _ := gzip.NewReader(fp)
	defer gzfp.Close()
	//buffp := brfp
	//buffp := bufio.NewReaderSize(brfp, 1<<20)
	var err error
	var num int
	for err != io.EOF {
		//buf := make([]byte, size)
		buf := bytesPool.Get().([]byte)
		buf = buf[:cap(buf)]
		//num, err = buffp.Read(buf)
		num, err = gzfp.Read(buf)
		if err != nil {
			if err != io.EOF {
				log.Fatalf("[ReadGzFile2] read file: %s, found err: %v\n", gzfn, err)
			}
		}

		fmt.Printf("[ReadGzFile2]read %d bytes\n", num)
		if num > 0 {
			cs <- buf[:num]
		}
		//time.Sleep(time.Second)
	}

	close(cs)
	//time.Sleep(time.Second * 10)
	return
}

const WindowSize = 1 << 16

func ReadZstdFile(zstdfn string, bytesPool *sync.Pool, cs chan<- []byte) {
	fp, err1 := os.Open(zstdfn)
	if err1 != nil {
		log.Fatalf("[ReadZstdFile] open file: %v failed..., err: %v\n", zstdfn, err1)
	}
	defer fp.Close()
	//size := (1 << 26)
	//brfp1 := cbrotli.NewReader(fp1)
	zr, err2 := zstd.NewReader(fp, zstd.WithDecoderConcurrency(1))
	if err2 != nil {
		log.Fatalf("[ReadZstdFile] zstd open file: %v failed..., err: %v\n", zstdfn, err1)
	}
	defer zr.Close()
	//buffp := brfp
	//buffp := bufio.NewReaderSize(brfp, 1<<20)
	var fullCount, emptyCount int
	var err error
	var num int
	buf := bytesPool.Get().([]byte)
	buf = buf[:cap(buf)]
	count := 0
	csSum := 0
	for err != io.EOF {
		//buf := make([]byte, size)
		//num, err = buffp.Read(buf)
		num, err = zr.Read(buf[count:])
		if err != nil {
			if err != io.EOF {
				log.Fatalf("[ReadZstdFile] read file: %s, found err: %v\n", zstdfn, err)
			}
		}
		count += num
		if err == nil && count < len(buf)*8/10 {
			continue
		}
		if len(cs) == cap(cs) {
			fullCount++
		} else if len(cs) == 0 {
			emptyCount++
		}
		cs <- buf[:count]
		csSum++
		buf = bytesPool.Get().([]byte)
		buf = buf[:cap(buf)]
		count = 0
		//time.Sleep(time.Second)
	}
	close(cs)
	fmt.Printf("[ReadZstdFile]emptyCount:%d fullCount:%d csSum:%d\n", emptyCount, fullCount, csSum)
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

/*func GetReadFileRecord(fncs <-chan []byte, recordChan chan<- ReadInfo, format string, bufSize int, bnt bool) int {
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
}*/

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

func GetIDSeq(nb []byte, bnt bool) (idSeq IDSeq) {
	idx := 0
	sz := bytes.IndexAny(nb[idx:], " \t\n")
	if nb[idx+sz-2] == '/' {
		idSeq.ID = nb[idx+1 : idx+sz-2]
	} else {
		idSeq.ID = nb[idx+1 : idx+sz]
	}
	idx += sz
	if nb[idx] != '\n' {
		sz = bytes.IndexByte(nb[idx:], '\n')
		idx += sz
	}
	idx++
	sz = bytes.IndexByte(nb[idx:], '\n')
	idSeq.Seq = nb[idx : idx+sz]
	if bnt {
		idSeq.Seq = Transform2BntByte2(idSeq.Seq)
	}
	return idSeq

}

func GetRecord(nb []byte, format string, bnt bool) (ri ReadInfo) {
	lastIdx := 0
	if format == "fa" {
		sz := bytes.IndexAny(nb[lastIdx:], " \t\n")
		var bs []byte
		if nb[lastIdx+sz-2] == '/' {
			bs = nb[lastIdx+1 : lastIdx+sz-2]
		} else {
			bs = nb[lastIdx+1 : lastIdx+sz]
		}
		id, _ := strconv.Atoi(string(bs))
		ri.ID = id
		lastIdx += sz + 1
		if nb[lastIdx-1] != '\n' {
			sz = bytes.IndexByte(nb[lastIdx:], '\n')
			ri.Anotition = nb[lastIdx : lastIdx+sz]
			lastIdx += sz + 1
		}
		sz = bytes.IndexByte(nb[lastIdx:], '\n')
		ri.Seq = nb[lastIdx : lastIdx+sz]
		if bnt {
			ri.Seq = Transform2BntByte2(ri.Seq)
		}
	} else {
		sz := bytes.IndexAny(nb[lastIdx:], " \t\n")
		var bs []byte
		if nb[lastIdx+sz-2] == '/' {
			bs = nb[lastIdx+1 : lastIdx+sz-2]
		} else {
			bs = nb[lastIdx+1 : lastIdx+sz]
		}
		id, _ := strconv.Atoi(string(bs))
		ri.ID = id
		lastIdx += sz + 1
		if nb[lastIdx-1] != '\n' {
			sz = bytes.IndexByte(nb[lastIdx:], '\n')
			ri.Anotition = nb[lastIdx : lastIdx+sz]
			lastIdx += sz + 1
		}
		sz = bytes.IndexByte(nb[lastIdx:], '\n')
		ri.Seq = nb[lastIdx : lastIdx+sz]
		if bnt {
			ri.Seq = Transform2BntByte2(ri.Seq)
		}
		lastIdx += sz + 1
		sz = bytes.IndexByte(nb[lastIdx:], '\n')
		lastIdx += sz + 1
		sz = bytes.IndexByte(nb[lastIdx:], '\n')
		ri.Qual = nb[lastIdx : lastIdx+sz]
	}
	//fmt.Printf("[GetRecordS]ri.ID:%s\nri.Seq:%s\nri.Qual:%s\n", string(ri.ID), string(ri.Seq), string(ri.Qual))
	return
}

func GetRecordSeq(nb []byte, format string, bnt bool) (ri ReadInfo) {
	lastIdx := 0
	if format == "fa" {
		sz := bytes.IndexByte(nb[lastIdx:], '\n')
		lastIdx += sz + 1
		sz = bytes.IndexByte(nb[lastIdx:], '\n')
		ri.Seq = nb[lastIdx : lastIdx+sz]
	} else {
		sz := bytes.IndexByte(nb[lastIdx:], '\n')
		lastIdx += sz + 1
		sz = bytes.IndexByte(nb[lastIdx:], '\n')
		ri.Seq = nb[lastIdx : lastIdx+sz]
		if bnt {
			ri.Seq = Transform2BntByte2(ri.Seq)
		}
	}
	//fmt.Printf("[GetRecordS]ri.ID:%s\nri.Seq:%s\nri.Qual:%s\n", string(ri.ID), string(ri.Seq), string(ri.Qual))
	return
}

/*func GetReadFileRecord(nb []byte, notFoundNewRecordBuf []byte, riBuf []ReadInfo, format string, bnt bool) ([]ReadInfo, int) {
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
	riBuf = append(riBuf, GetRecord(notFoundNewRecordBuf, format, true))
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
		riBuf = append(riBuf, GetRecord(nb[lastIdx:], format, bnt))
		lastIdx = idx
	}
	n = lastIdx
	//fmt.Printf("[GetReadFileRecord] read Record num is: %v\n", len(riBuf))
	//time.Sleep(time.Second * 60)
	return riBuf, n
}*/

/*func GetReadFileRecordS(nb []byte, notFoundNewRecordBuf []byte, riBuf []ReadInfoS, format string, bnt bool) ([]ReadInfoS, int) {
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
}*/

func GetNextReadSeq(bs []byte, format string) (seq []byte, idx int) {
	bl := 2
	if format == "fq" {
		bl = 4
	}
	//idx := 0
	for i := 0; i < bl; i++ {
		sz := bytes.IndexByte(bs[idx:], '\n')
		if sz < 0 {
			seq = nil
			idx = -1
			break
			//log.Fatalf("[GetNextReadSeq] not found proper Record, bs:%s\n", string(bs))
		}
		if i == 1 {
			seq = bs[idx : idx+sz]
		}
		idx += sz + 1
	}
	return seq, idx
}

func GetNextReadIDSeq(bs []byte, format string) (idSeq IDSeq, idx int) {
	bl := 2
	if format == "fq" {
		bl = 4
	}
	//idx := 0
	for i := 0; i < bl; i++ {
		sz := bytes.IndexByte(bs[idx:], '\n')
		if sz < 0 {
			idSeq.Seq = nil
			idx = -1
			break
			//log.Fatalf("[GetNextReadSeq] not found proper Record, bs:%s\n", string(bs))
		}
		if i == 0 {
			c := bytes.IndexAny(bs[idx:], " \t\n")
			if idx+c-2 > 0 && bs[idx+c-2] == '/' {
				idSeq.ID = bs[idx+1 : idx+c-2]
			} else {
				idSeq.ID = bs[idx+1 : idx+c]
			}
		}
		if i == 1 {
			idSeq.Seq = bs[idx : idx+sz]
		}
		idx += sz + 1
	}
	return idSeq, idx
}

func GetNextReadInfo(bs []byte, format string) (ri ReadInfo, idx int) {
	bl := 2
	if format == "fq" {
		bl = 4
	}
	//idx := 0
	for i := 0; i < bl; i++ {
		sz := bytes.IndexByte(bs[idx:], '\n')
		if sz < 0 {
			ri.Seq = nil
			idx = -1
			break
			//log.Fatalf("[GetNextReadSeq] not found proper Record, bs:%s\n", string(bs))
		}
		if i == 0 {
			c := bytes.IndexAny(bs[idx:], " \t\n")
			var id []byte
			if idx+c-2 > 0 && bs[idx+c-2] == '/' {
				id = bs[idx+1 : idx+c-2]
			} else {
				id = bs[idx+1 : idx+c]
			}
			uid, _ := strconv.Atoi(string(id))
			ri.ID = uid
			if bs[idx+c] != '\n' {
				ri.Anotition = bs[idx+c+1 : idx+sz]
			}
		} else if i == 1 {
			ri.Seq = bs[idx : idx+sz]
		} else if format == "fq" && i == 3 {
			ri.Qual = bs[idx : idx+sz]
		}
		idx += sz + 1
	}
	return ri, idx
}

func GetNextReadInfoS(bs []byte, format string) (ri ReadInfoS, idx int) {
	bl := 2
	if format == "fq" {
		bl = 4
	}
	//idx := 0
	for i := 0; i < bl; i++ {
		sz := bytes.IndexByte(bs[idx:], '\n')
		if sz < 0 {
			ri.Seq = nil
			idx = -1
			break
			//log.Fatalf("[GetNextReadSeq] not found proper Record, bs:%s\n", string(bs))
		}
		if i == 0 {
			c := bytes.IndexAny(bs[idx:], " \t\n")
			if c > 2 && bs[idx+c-2] == '/' {
				c -= 2
			}
			ri.ID = bs[idx+1 : idx+c]
			if bs[idx+c] != '\n' {
				ri.Anotition = bs[idx+c+1 : idx+sz]
			}
		} else if i == 1 {
			ri.Seq = bs[idx : idx+sz]
		} else if format == "fq" && i == 3 {
			ri.Qual = bs[idx : idx+sz]
		}
		idx += sz + 1
	}
	return ri, idx
}

func GetReadFileSeq(sp *SeqPool, format string, bnt bool) (idx int) {

	if len(sp.NotFoundRecordBytes) > 0 {
		var blockLineNum int
		if format == "fa" {
			blockLineNum = 2
		} else { // format == "fq"
			blockLineNum = 4
		}
		m := bytes.Count(sp.NotFoundRecordBytes, []byte("\n"))
		if m >= blockLineNum {
			log.Fatalf("m:%d > blockLineNum:%d\n", m, blockLineNum)
		}
		i := 0
		for ; i < blockLineNum-m; i++ {
			if idx >= len(sp.Cs) {
				break
			}
			sz := bytes.IndexByte(sp.Cs[idx:], '\n')
			if sz < 0 {
				break
			}
			idx += sz + 1
		}
		if i != blockLineNum-m {
			log.Fatalf("[GetReadFileSeq] not found proper Record, nb:%s\n", string(sp.Cs))
		}
		sp.NotFoundRecordBytes = append(sp.NotFoundRecordBytes, sp.Cs[:idx]...)
		seq, _ := GetNextReadSeq(sp.NotFoundRecordBytes, format)
		if bnt {
			seq = Transform2BntByte2(seq)
		}
		sp.SeqArr = append(sp.SeqArr, seq)
	}

	for {
		sq, nextIdx := GetNextReadSeq(sp.Cs[idx:], format)
		if nextIdx < 0 {
			break
		}
		if bnt {
			sq = Transform2BntByte2(sq)
		}
		sp.SeqArr = append(sp.SeqArr, sq)
		idx += nextIdx
	}
	return
}

func GetReadFileIDSeq(sp *IDSeqPool, format string, bnt bool) (idx int) {

	if len(sp.NotFoundRecordBytes) > 0 {
		var blockLineNum int
		if format == "fa" {
			blockLineNum = 2
		} else { // format == "fq"
			blockLineNum = 4
		}
		m := bytes.Count(sp.NotFoundRecordBytes, []byte("\n"))
		if m >= blockLineNum {
			log.Fatalf("m:%d > blockLineNum:%d\n", m, blockLineNum)
		}
		i := 0
		for ; i < blockLineNum-m; i++ {
			if idx >= len(sp.Cs) {
				break
			}
			sz := bytes.IndexByte(sp.Cs[idx:], '\n')
			if sz < 0 {
				break
			}
			idx += sz + 1
		}
		if i != blockLineNum-m {
			log.Fatalf("[GetReadFileIDSeq] not found proper Record, nb:%s\n", string(sp.Cs))
		}
		sp.NotFoundRecordBytes = append(sp.NotFoundRecordBytes, sp.Cs[:idx]...)
		idSeq, _ := GetNextReadIDSeq(sp.NotFoundRecordBytes, format)
		if bnt {
			Transform2BntByte2(idSeq.Seq)
		}
		sp.IDSeqArr = append(sp.IDSeqArr, idSeq)
	}

	for {
		sq, nextIdx := GetNextReadIDSeq(sp.Cs[idx:], format)
		if nextIdx < 0 {
			break
		}
		if bnt {
			Transform2BntByte2(sq.Seq)
		}
		sp.IDSeqArr = append(sp.IDSeqArr, sq)
		idx += nextIdx
	}

	return
}

func GetReadFileReadInfo(rp *RIPool, format string, bnt bool, pair bool) (idx int) {
	if len(rp.NotFoundRecordBytes) > 0 {
		var blockLineNum int
		if format == "fa" {
			blockLineNum = 2
		} else { // format == "fq"
			blockLineNum = 4
		}
		if pair {
			blockLineNum *= 2
		}
		m := bytes.Count(rp.NotFoundRecordBytes, []byte("\n"))
		if m >= blockLineNum {
			log.Fatalf("m:%d > blockLineNum:%d\n", m, blockLineNum)
		}
		i := 0
		for ; i < blockLineNum-m; i++ {
			if idx >= len(rp.Cs) {
				break
			}
			sz := bytes.IndexByte(rp.Cs[idx:], '\n')
			if sz < 0 {
				break
			}
			idx += sz + 1
		}
		if i != blockLineNum-m {
			//log.Fatalf("[GetReadFileReadInfo] not found proper Record, nb:%s\n", string(rp.Cs))
			//rp.NotFoundRecordBytes = append(rp.NotFoundRecordBytes, rp.Cs...)
			idx = -1
			return
		}
		rp.NotFoundRecordBytes = append(rp.NotFoundRecordBytes, rp.Cs[:idx]...)
		ri, p := GetNextReadInfo(rp.NotFoundRecordBytes, format)
		if p < 0 {
			return
		}
		if bnt {
			Transform2BntByte2(ri.Seq)
		}
		rp.RIArr = append(rp.RIArr, ri)

		if pair {
			ri, p = GetNextReadInfo(rp.NotFoundRecordBytes[p:], format)
			if p < 0 {
				return
			}
			if bnt {
				Transform2BntByte2(ri.Seq)
			}
			rp.RIArr = append(rp.RIArr, ri)
		}
	}
	if pair {
		for {
			ri1, idx1 := GetNextReadInfo(rp.Cs[idx:], format)
			if idx1 < 0 {
				break
			}
			ri2, idx2 := GetNextReadInfo(rp.Cs[idx+idx1:], format)
			if idx2 < 0 {
				break
			}
			if bnt {
				Transform2BntByte2(ri1.Seq)
				Transform2BntByte2(ri2.Seq)
			}
			rp.RIArr = append(rp.RIArr, ri1, ri2)
			idx += idx1 + idx2
		}
	} else {
		for {
			ri1, idx1 := GetNextReadInfo(rp.Cs[idx:], format)
			if idx1 < 0 {
				break
			}
			if bnt {
				Transform2BntByte2(ri1.Seq)
			}
			rp.RIArr = append(rp.RIArr, ri1)
			idx += idx1
		}
	}

	return
}

func GetReadFileReadInfoS(rp *RISPool, format string, bnt bool, pair bool) (idx int) {

	if len(rp.NotFoundRecordBytes) > 0 {
		var blockLineNum int
		if format == "fa" {
			blockLineNum = 2
		} else { // format == "fq"
			blockLineNum = 4
		}
		if pair {
			blockLineNum *= 2
		}
		m := bytes.Count(rp.NotFoundRecordBytes, []byte("\n"))
		if m >= blockLineNum {
			log.Fatalf("m:%d > blockLineNum:%d\n", m, blockLineNum)
		}
		i := 0
		for ; i < blockLineNum-m; i++ {
			if idx >= len(rp.Cs) {
				break
			}
			sz := bytes.IndexByte(rp.Cs[idx:], '\n')
			if sz < 0 {
				break
			}
			idx += sz + 1
		}
		if i != blockLineNum-m {
			log.Fatalf("[GetReadFileReadInfo] not found proper Record, nb:%s\n", string(rp.Cs))
		}
		rp.NotFoundRecordBytes = append(rp.NotFoundRecordBytes, rp.Cs[:idx]...)
		ri, p := GetNextReadInfoS(rp.NotFoundRecordBytes, format)
		if bnt {
			Transform2BntByte2(ri.Seq)
		}
		rp.RIArr = append(rp.RIArr, ri)

		if pair {
			ri, _ = GetNextReadInfoS(rp.NotFoundRecordBytes[p:], format)
			if bnt {
				Transform2BntByte2(ri.Seq)
			}
			rp.RIArr = append(rp.RIArr, ri)
		}
	}

	if pair {
		for {
			ri1, idx1 := GetNextReadInfoS(rp.Cs[idx:], format)
			if idx1 < 0 {
				break
			}
			ri2, idx2 := GetNextReadInfoS(rp.Cs[idx+idx1:], format)
			if idx2 < 0 {
				break
			}
			if bnt {
				Transform2BntByte2(ri1.Seq)
				Transform2BntByte2(ri2.Seq)
			}
			rp.RIArr = append(rp.RIArr, ri1, ri2)
			idx += idx1 + idx2
		}
	} else {
		for {
			ri1, idx1 := GetNextReadInfoS(rp.Cs[idx:], format)
			if idx1 < 0 {
				break
			}
			if bnt {
				Transform2BntByte2(ri1.Seq)
			}
			rp.RIArr = append(rp.RIArr, ri1)
			idx += idx1
		}
	}

	return
}

/*func GetReadFileReadInfo(cs, notFoundNewRecord []byte, recordChan chan<- ReadInfo, format string, bnt bool) (idx int) {

	if len(notFoundNewRecord) > 0 {
		var blockLineNum int
		if format == "fa" {
			blockLineNum = 2
		} else { // format == "fq"
			blockLineNum = 4
		}
		m := bytes.Count(notFoundNewRecord, []byte("\n"))
		if m >= blockLineNum {
			log.Fatalf("m:%d > blockLineNum:%d\n", m, blockLineNum)
		}
		i := 0
		for ; i < blockLineNum-m; i++ {
			if idx >= len(cs) {
				break
			}
			sz := bytes.IndexByte(cs[idx:], '\n')
			if sz < 0 {
				break
			}
			idx += sz + 1
		}
		if i != blockLineNum-m {
			log.Fatalf("[GetReadFileIDSeq] not found proper Record, nb:%s\n", string(cs))
		}
		notFoundNewRecord = append(notFoundNewRecord, cs[:idx]...)
		ri, _ := GetNextReadInfo(notFoundNewRecord, format)
		if bnt {
			Transform2BntByte2(ri.Seq)
		}
		recordChan <- ri
	}

	for {
		ri, nextIdx := GetNextReadInfo(cs[idx:], format)
		if nextIdx < 0 {
			break
		}
		if bnt {
			Transform2BntByte2(ri.Seq)
		}
		recordChan <- ri
		idx += nextIdx
	}

	return
}*/

func GetReadSeqBucket(fn string, cs <-chan []byte, seqArrPool *sync.Pool, seqPoolChan chan<- SeqPool) (readNum int) {
	//var processNumReads int
	//var bucketCount int
	format := GetReadsFileFormat(fn)
	//notFoundNewRecord := notFoundNewRecordPool.Get().([]byte)
	var notFoundNewRecord []byte
	//bufSize := (1 << 20)
	//size := 8000
	//fncs1 := make(chan []byte, 10)
	//recordChan1 := make(chan ReadInfo, size)
	fmt.Printf("[GetReadSeqBucket] begin processe file: %v...\n", fn)
	//go ReadBrFile(fn, fncs1, bufSize)
	//go GetReadFileRecord(fncs1, recordChan1, format, bufSize, true)

	//var count int
	for {
		bs, ok := <-cs
		if !ok {
			if len(notFoundNewRecord) != 0 {
				log.Fatalf("[GetReadSeqBucket] to the file end, but len(notFoundNewRecord) = %d\n", len(notFoundNewRecord))
			}
			break
		}

		var sp SeqPool
		sp.SeqArr = seqArrPool.Get().([][]byte)
		sp.SeqArr = sp.SeqArr[:0]
		sp.Cs = bs
		sp.NotFoundRecordBytes = notFoundNewRecord
		idx := GetReadFileSeq(&sp, format, true)
		readNum += len(sp.SeqArr)
		var nf []byte
		if idx < len(sp.Cs) {
			nf = append(nf, sp.Cs[idx:]...)
			//sp.Cs = sp.Cs[:idx]
		}
		notFoundNewRecord = nf
		//fmt.Printf("[GetReadSeqBucket] processed reads number is: %d,finished processed file: %v\n", len(sp.Cs),len(notFoundNewRecord), len(sp.SeqArr))
		seqPoolChan <- sp
	}
	// send read finish signal
	close(seqPoolChan)
	fmt.Printf("[GetReadSeqBucket] processed reads number is: %d,finished processed file: %v\n", readNum, fn)
	return
}

func GetReadIDSeqBucket(fn string, cs <-chan []byte, idSeqArrPool *sync.Pool, idSeqPoolChan chan<- IDSeqPool) (readNum int) {
	//var processNumReads int
	//var bucketCount int
	format := GetReadsFileFormat(fn)
	var notFoundNewRecord []byte
	//bufSize := (1 << 20)
	//size := 8000
	//fncs1 := make(chan []byte, 10)
	//recordChan1 := make(chan ReadInfo, size)
	fmt.Printf("[GetReadIDSeqBucket] begin processe file: %v...\n", fn)
	//go ReadBrFile(fn, fncs1, bufSize)
	//go GetReadFileRecord(fncs1, recordChan1, format, bufSize, true)

	//var count int
	for {
		bs, ok := <-cs
		if !ok {
			if len(notFoundNewRecord) != 0 {
				log.Fatalf("[GetReadIDSeqBucket] to the file end, but len(notFoundNewRecord) = %d\n", len(notFoundNewRecord))
			}
			break
		}

		var sp IDSeqPool
		sp.IDSeqArr = idSeqArrPool.Get().([]IDSeq)
		sp.IDSeqArr = sp.IDSeqArr[:0]
		sp.Cs = bs
		sp.NotFoundRecordBytes = notFoundNewRecord
		idx := GetReadFileIDSeq(&sp, format, true)
		readNum += len(sp.IDSeqArr)
		var nf []byte
		if idx < len(sp.Cs) {
			nf = append(nf, sp.Cs[idx:]...)
		}
		notFoundNewRecord = nf
		idSeqPoolChan <- sp
	}
	// send read finish signal
	close(idSeqPoolChan)
	fmt.Printf("[GetReadIDSeqBucket] processed reads number is: %d,finished processed file: %v\n", readNum, fn)
	return
}

func GetReadInfoBucket(fn string, cs <-chan []byte, riArrPool *sync.Pool, riPoolChan chan<- RIPool, pair bool) (readNum int) {
	//var processNumReads int
	//var bucketCount int
	format := GetReadsFileFormat(fn)
	var notFoundNewRecord []byte
	//bufSize := (1 << 20)
	//size := 8000
	//fncs1 := make(chan []byte, 10)
	//recordChan1 := make(chan ReadInfo, size)
	fmt.Printf("[GetReadInfoBucket] begin processe file: %v...\n", fn)
	//go ReadBrFile(fn, fncs1, bufSize)
	//go GetReadFileRecord(fncs1, recordChan1, format, bufSize, true)

	//var count int
	for {
		bs, ok := <-cs
		if !ok {
			if len(notFoundNewRecord) != 0 {
				log.Fatalf("[GetReadInfoBucket] to the file end, but len(notFoundNewRecord) = %d\n", len(notFoundNewRecord))
			}
			break
		}

		var rp RIPool
		rp.RIArr = riArrPool.Get().([]ReadInfo)
		rp.RIArr = rp.RIArr[:0]
		rp.Cs = bs
		rp.NotFoundRecordBytes = notFoundNewRecord
		idx := GetReadFileReadInfo(&rp, format, true, pair)
		readNum += len(rp.RIArr)
		var nf []byte
		if idx < len(rp.Cs) {
			nf = append(nf, rp.Cs[idx:]...)
		}
		notFoundNewRecord = nf
		riPoolChan <- rp
	}
	// send read finish signal
	close(riPoolChan)
	fmt.Printf("[GetReadInfoBucket] processed reads number is: %d,finished processed file: %v\n", readNum, fn)
	return
}

func GetReadInfoSBucket(fn string, cs <-chan []byte, riArrPool *sync.Pool, riPoolChan chan<- RISPool, pair bool) (readNum int) {
	//var processNumReads int
	//var bucketCount int
	format := GetReadsFileFormat(fn)
	var notFoundNewRecord []byte
	//bufSize := (1 << 20)
	//size := 8000
	//fncs1 := make(chan []byte, 10)
	//recordChan1 := make(chan ReadInfo, size)
	fmt.Printf("[GetReadInfoBucket] begin processe file: %v...\n", fn)
	//go ReadBrFile(fn, fncs1, bufSize)
	//go GetReadFileRecord(fncs1, recordChan1, format, bufSize, true)

	//var count int
	for {
		bs, ok := <-cs
		if !ok {
			if len(notFoundNewRecord) != 0 {
				log.Fatalf("[GetReadInfoBucket] to the file end, but len(notFoundNewRecord) = %d\n", len(notFoundNewRecord))
			}
			break
		}

		var rp RISPool
		rp.RIArr = riArrPool.Get().([]ReadInfoS)
		rp.RIArr = rp.RIArr[:0]
		rp.Cs = bs
		rp.NotFoundRecordBytes = notFoundNewRecord
		idx := GetReadFileReadInfoS(&rp, format, true, pair)
		readNum += len(rp.RIArr)
		var nf []byte
		if idx < len(rp.Cs) {
			nf = append(nf, rp.Cs[idx:]...)
		}
		notFoundNewRecord = nf
		riPoolChan <- rp
	}
	// send read finish signal
	close(riPoolChan)
	fmt.Printf("[GetReadInfoBucket] processed reads number is: %d,finished processed file: %v\n", readNum, fn)
	return
}

/*func GetReadFileRecord(fn string, cs <-chan []byte, recordChan chan<- ReadInfo) {
	//var processNumReads int
	//var bucketCount int
	format := GetReadsFileFormat(fn)
	var notFoundNewRecord []byte
	//bufSize := (1 << 20)
	//size := 8000
	//fncs1 := make(chan []byte, 10)
	//recordChan1 := make(chan ReadInfo, size)
	fmt.Printf("[GetReadFileRecord] begin processe file: %v...\n", fn)
	//go ReadBrFile(fn, fncs1, bufSize)
	//go GetReadFileRecord(fncs1, recordChan1, format, bufSize, true)

	//var count int
	for {
		bs, ok := <-cs
		if !ok {
			if len(notFoundNewRecord) != 0 {
				log.Fatalf("[GetReadFileRecord] to the file end, but len(notFoundNewRecord) = %d\n", len(notFoundNewRecord))
			}
			break
		}
		idx := GetReadFileReadInfo(bs, notFoundNewRecord, recordChan, format, true)
		var nr []byte
		if idx < len(bs) {
			nr = append(nr, bs[idx:]...)
		}
		notFoundNewRecord = nr
	}
	// send read finish signal
	close(recordChan)
	//fmt.Printf("[GetReadFileRecord] processed reads number is: %d,finished processed file: %v\n", readNum, fn)
	return
}*/

func CCF(c cli.Command) {
	//fmt.Println(c.Flags(), c.Parent().Flags())
	opt, suc := checkArgsCCF(c)
	if suc == false {
		log.Fatalf("[CCF] check Arguments error, opt:%v\n", opt)
	}
	fmt.Printf("[CCF]opt:%v\n", opt)
	/*profileFn := opt.Prefix + ".CCF.prof"
	cpuprofilefp, err := os.Create(profileFn)
	if err != nil {
		log.Fatalf("[CCF] open cpuprofile file: %v failed\n", profileFn)
	}
	pprof.StartCPUProfile(cpuprofilefp)
	defer pprof.StopCPUProfile()*/

	//defer profile.Start().Stop()
	cfgInfo, err := ParseCfg(opt.CfgFn, opt.Merge, opt.Correct)
	if err != nil {
		log.Fatalf("[CCF] ParseCfg 'C': %v err :%v\n", opt.CfgFn, err)
	}
	fmt.Println(cfgInfo)

	// make CuckooFilter

	t0 := time.Now()
	cf := MakeCuckooFilter(uint64(opt.CFSize), opt.Kmer)
	numCPU := opt.NumCPU
	runtime.GOMAXPROCS(numCPU)
	//bufsize := 100
	//we := make(chan int)
	//defer close(we)
	var fnArr []string
	for _, lib := range cfgInfo.Libs {
		if lib.AsmFlag != AllState && lib.SeqProfile != 1 {
			continue
		}
		fnArr = append(fnArr, lib.FnName...)
		fnArr = append(fnArr, lib.UnMergedFn...)
	}

	//fmt.Printf("[CCF] fileName array: %v\n", fnArr)
	//concurrentNum := 6
	//totalNumT := opt.NumCPU/(concurrentNum+1) + 1
	//concurrentNum := opt.Kmer / 20
	concurrentNum := opt.NumCPU - 1

	// write goroutinue
	wbytesPool := sync.Pool{New: func() interface{} {
		bytes := make([]byte, WindowSize)
		return bytes
	}}
	wc := make(chan []byte, concurrentNum)
	wrfn := opt.Prefix + ".uniqkmerseq.zst"
	cc := make(chan int)
	defer close(cc)
	go WriteZstd(wrfn, wc, &wbytesPool, cc)

	totalKmerNum := 0
	for _, fn := range fnArr {
		//fmt.Printf("[CCF] processing file: %v\n", lib.FnName[i])
		totalKmerNum += ConcurrentConstructCF(fn, &cf, wc, &wbytesPool, concurrentNum, opt.Kmer, opt.MinKmerFreq)
	}
	close(wc)
	wn := <-cc
	if wn%opt.Kmer != 0 {
		log.Fatalf("wn:%d wn//pt.Kmer:%d != 0\n", wn, wn%opt.Kmer)
	}
	writeKmerNum := wn / opt.Kmer
	kmerInfofn := opt.Prefix + ".uniqkmerseq.info"
	WriteUniqKmerseqInfo(kmerInfofn, writeKmerNum)
	fmt.Printf("[CCF] total write unique kmer number is:%d cf.Count:%d\n", writeKmerNum, cf.Count)
	//time.Sleep(time.Second * 10)
	fmt.Printf("[CCF] totalKmerNum:%d\n", totalKmerNum)
	fmt.Printf("[CCF] construct CuckooFilter took %v to run\n", time.Now().Sub(t0))

	// end signal from write goroutinue
	// prefix := c.Parent().Flag("p").String()
	// cfinfofn := prefix + ".cfInfo"
	// if err := cf.WriteCuckooFilterInfo(cfinfofn); err != nil {
	// 	log.Fatal(err)
	// }
	t0 = time.Now()
	/*if !opt.Correct {
		cffn := opt.Prefix + ".cf.Hash"
		err = cf.HashWriter(cffn)
		if err != nil {
			log.Fatalf("[CCF]HashWriter file: %v error: %v\n", cffn, err)
		}
		fmt.Printf("[CCF] write CuckooFilter hash used:  %v to run\n", time.Now().Sub(t0))
	}*/

	// output stat
	cf.GetStat()
	if err != nil {
		log.Fatal("error found in the CCF")
	}
	return
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
	input1 := c.Flag("read1").String()
	input2 := c.Flag("read2").String()
	IDFile := c.Flag("IDFile").String()
	prefix := c.Flag("prefix").String()
	//format := c.Flag("format").String()
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
		//size := 5000
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

		seqArrPool := sync.Pool{New: func() interface{} {
			arr := make([]ReadInfoS, bufSize/(ReadLen*2))
			return arr
		}}

		wbytesPool := sync.Pool{New: func() interface{} {
			bytes := make([]byte, bufSize+nextSize)
			return bytes
		}}
		//riPool := sync.Pool{New: func() interface{} {
		//	readInfos := make([]ReadInfo, size)
		//	return readInfos
		//}}

		//var f string
		//if input == "/dev/stdin" {
		//f = "fq"
		//go ReadStdin(fncs, &rbytesPool)
		//f = GetReadsFileFormat(input)
		idMap := make(map[string]uint8)
		{
			fncs := make(chan []byte)
			seqPoolChan := make(chan RISPool, 2)
			go ReadGzFile2(input1, &rbytesPool, fncs)
			go GetReadInfoSBucket(input1, fncs, &seqArrPool, seqPoolChan, false)
			for {
				sp, ok := <-seqPoolChan
				if !ok {
					break
				}
				for _, seq := range sp.RIArr {
					idMap[string(seq.ID)]++
				}
				rbytesPool.Put(sp.Cs)
				seqArrPool.Put(sp.RIArr)
			}
		}

		{
			fncs := make(chan []byte)
			seqPoolChan := make(chan RISPool, 2)
			go ReadGzFile2(input2, &rbytesPool, fncs)
			go GetReadInfoSBucket(input2, fncs, &seqArrPool, seqPoolChan, false)
			for {
				sp, ok := <-seqPoolChan
				if !ok {
					break
				}
				for _, seq := range sp.RIArr {
					idMap[string(seq.ID)]++
				}
				rbytesPool.Put(sp.Cs)
				seqArrPool.Put(sp.RIArr)
			}
		}

		{
			pairNum, singleNum = 0, 0
			wfn := prefix + "1.fq.gz"
			wc := make(chan []byte, 2)
			wfinish := make(chan bool, 1)
			go WriteGz(wfn, wc, &wbytesPool, wfinish)
			fncs := make(chan []byte)
			seqPoolChan := make(chan RISPool, 2)
			go ReadGzFile2(input1, &rbytesPool, fncs)
			go GetReadInfoSBucket(input1, fncs, &seqArrPool, seqPoolChan, false)
			for {
				sp, ok := <-seqPoolChan
				if !ok {
					break
				}
				wbp := wbytesPool.Get().([]byte)
				wbp = wbp[:0]
				for _, seq := range sp.RIArr {
					if idMap[string(seq.ID)] == 2 {
						wbp = append(wbp, '@')
						wbp = append(wbp, seq.ID...)
						wbp = append(wbp, '\t')
						wbp = append(wbp, seq.Anotition...)
						wbp = append(wbp, '\n')
						Transform2Char2(seq.Seq)
						wbp = append(wbp, seq.Seq...)
						wbp = append(wbp, "\n+\n"...)
						wbp = append(wbp, seq.Qual...)
						wbp = append(wbp, '\n')
						pairNum++
					} else {
						singleNum++
					}
				}
				wc <- wbp
				rbytesPool.Put(sp.Cs)
				seqArrPool.Put(sp.RIArr)
			}
			close(wc)
			<-wfinish
			fmt.Printf("[ExtractPairEnd] pairNum:%d single:%d\n", pairNum, singleNum)
		}

		{
			pairNum, singleNum = 0, 0
			wfn := prefix + "2.fq.gz"
			wc := make(chan []byte, 2)
			wfinish := make(chan bool, 1)
			go WriteGz(wfn, wc, &wbytesPool, wfinish)
			fncs := make(chan []byte)
			seqPoolChan := make(chan RISPool, 2)
			go ReadGzFile2(input2, &rbytesPool, fncs)
			go GetReadInfoSBucket(input2, fncs, &seqArrPool, seqPoolChan, false)
			for {
				sp, ok := <-seqPoolChan
				if !ok {
					break
				}
				wbp := wbytesPool.Get().([]byte)
				wbp = wbp[:0]
				for _, seq := range sp.RIArr {
					if idMap[string(seq.ID)] == 2 {
						wbp = append(wbp, '@')
						wbp = append(wbp, seq.ID...)
						wbp = append(wbp, '\t')
						wbp = append(wbp, seq.Anotition...)
						wbp = append(wbp, '\n')
						Transform2Char2(seq.Seq)
						wbp = append(wbp, seq.Seq...)
						wbp = append(wbp, "\n+\n"...)
						wbp = append(wbp, seq.Qual...)
						wbp = append(wbp, '\n')
						pairNum++
					} else {
						singleNum++
					}
				}
				wc <- wbp
				rbytesPool.Put(sp.Cs)
				seqArrPool.Put(sp.RIArr)
			}
			close(wc)
			<-wfinish
			fmt.Printf("[ExtractPairEnd] pairNum:%d single:%d\n", pairNum, singleNum)
		}
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
	//f := GetReadsFileFormat(inputfq)

	go ReadGzFile(inputfq, fncs, bufSize)
	//go GetReadFileRecordS(fncs, recordChan, f, bufSize, false)

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
		ri.ID = startID + passNum
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
