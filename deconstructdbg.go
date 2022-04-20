package main

import (
	"fmt"
	"io"
	"log"
	"math"
	"os"
	"runtime"
	"runtime/pprof"
	"sort"
	"strconv"
	"sync"
	"time"

	"github.com/jwaldrip/odin/cli"
	"github.com/klauspost/compress/zstd"
)

const BitNum = 2
const SeedLen = 21 // SeedLen allow <=23
const SeedMUSK = (uint64(1) << (uint64(SeedLen) * NumBitsInBase)) - 1
const SeedBITNUM = (uint64(SeedLen) - 1) * NumBitsInBase

const SyncmerMinimizerLen = 5 // must <= 8
const SyncmerMinimizerMUSK = (uint16(1) << (uint16(SyncmerMinimizerLen) * NumBitsInBase)) - 1
const SyncmerMinimizerBITNUM = (uint64(SyncmerMinimizerLen) - 1) * NumBitsInBase
const SyncmerPos1 = 3
const SyncmerPos2 = SeedLen - SyncmerPos1 - SyncmerMinimizerLen

const MaxAllowSeqLen = (1 << seqPositionBitNum) - 1
const MinEdgeLenDiff = 20

//"Debug const set Debug model if true"
var DebugModel bool = false
var BuckReadsNum int = 100
var BuckSize = BuckReadsNum * 20000
var EdgeMapKmerMin int = 60 // minimum kmer number of edge used map Long Reads
var MinMinimerNum int = 5
var MinEdgeLen int = 300 // kmerlen*2 - MinEdgeLenDiff

const mapLenBitNum = 32 - 1 - seqPositionBitNum
const mapingKmerInfoMaxMapLen = (1 << mapLenBitNum) - 1
const infoReSetMapLen = (1 << (seqPositionBitNum + 1)) - 1
const ePositionMask = (1 << seqPositionBitNum) - 1
const infoReSetPosition = (((1 << mapLenBitNum) - 1) << (seqPositionBitNum + 1)) | 0x1

type MapingKmerInfo struct {
	SeedInfo // kmerBase|Position|Strand [SeedLen*NumBitsInBase:seqPositionBitNum:1] just store 2-bits base, allow max 32 kmer base
	EID      uint32
	Info     uint32 // MapLen|Position|Strand [32-1-seqPositionBitNum:seqPositionBitNum:1]
	//Len          uint16 // mapping base length
}

func (m *MapingKmerInfo) GetMapStrand() bool {
	if (uint64(m.SeedInfo) & 0x1) == uint64(m.Info&0x1) {
		return true
	} else {
		return false
	}
}

func (m *MapingKmerInfo) GetMapLen() uint32 {
	return (m.Info >> (seqPositionBitNum + 1))
}

func (m *MapingKmerInfo) SetMapLen(p uint32) {
	if p > mapingKmerInfoMaxMapLen {
		p = mapingKmerInfoMaxMapLen
	}
	m.Info = (p << (seqPositionBitNum + 1)) | (m.Info & infoReSetMapLen)
}

func (m *MapingKmerInfo) GetEPosition() uint32 {
	return (m.Info >> 1) & positionMask
}
func (m *MapingKmerInfo) SetEPosition(p uint32) {
	m.Info = (p << 1) | (m.Info & infoReSetPosition)
}

func (m *MapingKmerInfo) GetEStrand() bool {
	if (m.Info & 0x1) > 0 {
		return true
	} else {
		return false
	}
}

func (m *MapingKmerInfo) SetEStrand(s bool) {
	if s == PLUS {
		m.Info |= 0x1
	} else {
		m.Info >>= 1
		m.Info <<= 1
	}
}

type MapPos struct {
	Info uint32 // Position|Strand:31|1
}

const StrandMask = ((1 << 31) - 1) << 1

func (m *MapPos) GetPosition() uint32 {
	return (m.Info >> 1)
}
func (m *MapPos) SetPosition(p uint32) {
	m.Info = (p << 1) | (m.Info & 0x1)
}

func (m *MapPos) GetStrand() bool {
	if (m.Info & 0x1) > 0 {
		return true
	} else {
		return false
	}
}

func (m *MapPos) SetStrand(s bool) {
	if s == PLUS {
		m.Info |= 0x1
	} else {
		m.Info >>= 1
		m.Info <<= 1
	}
}

type optionsDDBG struct {
	ArgsOpt
	//MaxMapEdgeLen int
	SeedLen                      int
	MinDepth                     int
	AvgDepth                     int
	MinUniqueEdgeLen             int
	AvgReadLen                   int
	WinSize                      int
	MaxNGSReadLen                int
	MinMapFreq                   int
	ExtLen                       int
	MinChainScoreIdentityPercent float32
	ONTFn                        string
	Correct                      bool
	Haplotype                    bool
	Debug                        bool
}

func checkArgsDDBG(c cli.Command) (opt optionsDDBG, succ bool) {
	gOpt, suc := CheckGlobalArgs(c.Parent())
	if suc == false {
		log.Fatalf("[checkArgsDDBG] check global Arguments error, opt: %v\n", gOpt)
	}
	opt.ArgsOpt = gOpt
	var ok bool
	opt.MinUniqueEdgeLen, ok = c.Flag("MinUniqueEdgeLen").Get().(int)
	if !ok {
		log.Fatalf("[checkArgsDDBG] argument 'MinUniqueEdgeLen': %v set error\n ", c.Flag("MinUniqueEdgeLen").String())
	}
	if opt.MinUniqueEdgeLen < 200 || opt.MinUniqueEdgeLen > 10000 {
		log.Fatalf("[checkArgsDDBG] argument 'MinUniqueEdgeLen': %v must between [200~10k]\n", c.Flag("MinUniqueEdgeLen"))
	}
	opt.AvgReadLen, ok = c.Flag("AvgReadLen").Get().(int)
	if !ok {
		log.Fatalf("[checkArgsDDBG] argument 'AvgReadLen': %v set error\n ", c.Flag("AvgReadLen").String())
	}
	if opt.AvgReadLen < 1000 || opt.AvgReadLen > 1000000 {
		log.Fatalf("[checkArgsDDBG] argument 'AvgReadLen': %v must between [1k~1m]\n", c.Flag("AvgReadLen"))
	}

	opt.WinSize, ok = c.Flag("WinSize").Get().(int)
	if !ok {
		log.Fatalf("[checkArgsDDBG] argument 'WinSize': %v set error\n ", c.Flag("WinSize").String())
	}
	if opt.WinSize < 1 || opt.WinSize > 100 {
		log.Fatalf("[checkArgsDDBG] argument 'WinSize': %v must between 1~100\n", c.Flag("WinSize"))
	}
	opt.MaxNGSReadLen, ok = c.Flag("MaxNGSReadLen").Get().(int)
	if !ok {
		log.Fatalf("[checkArgsDDBG] argument 'MaxNGSReadLen': %v set error\n ", c.Flag("MaxNGSReadLen").String())
	}
	if opt.MaxNGSReadLen < opt.Kmer+50 {
		log.Fatalf("[checkArgsDDBG] argument 'MaxNGSReadLen': %v must bigger than K+50\n", c.Flag("MaxNGSReadLen").String())
	}

	opt.MinMapFreq, ok = c.Flag("MinMapFreq").Get().(int)
	if !ok {
		log.Fatalf("[checkArgsDDBG] argument 'MinMapFreq': %v set error\n ", c.Flag("MinMapFreq").String())
	}
	if opt.MinMapFreq < 5 && opt.MinMapFreq >= 20 {
		log.Fatalf("[checkArgsDDBG] argument 'MinMapFreq': %v must 5 <= MinMapFreq < 20\n", c.Flag("MinMapFreq").String())
	}

	opt.ExtLen, ok = c.Flag("ExtLen").Get().(int)
	if !ok {
		log.Fatalf("[checkArgsDDBG] argument 'ExtLen': %v set error\n ", c.Flag("ExtLen").String())
	}
	if opt.ExtLen < 1000 {
		log.Fatalf("[checkArgsDDBG] argument 'ExtLen': %v must bigger than 1000\n", c.Flag("ExtLen").String())
	}
	var minPerc float64
	minPerc, ok = c.Flag("MinChainScoreIdentityPercent").Get().(float64)
	if !ok {
		log.Fatalf("[checkArgsDDBG] argument 'MinChainScoreIdentityPercent': %v set error\n ", c.Flag("MinChainScoreIdentityPercent").String())
	}
	if minPerc <= 0 || minPerc >= 1 {
		log.Fatalf("[checkArgsDDBG] argument 'MinChainScoreIdentityPercent': %v must between [1~100]\n", c.Flag("MinChainScoreIdentityPercent").String())
	}
	opt.MinChainScoreIdentityPercent = float32(minPerc)

	opt.Correct, ok = c.Flag("Correct").Get().(bool)
	if !ok {
		log.Fatalf("[checkArgsDDBG] argument 'Correct': %v set error\n ", c.Flag("Correct").String())
	}
	opt.Haplotype, ok = c.Flag("Haplotype").Get().(bool)
	if !ok {
		log.Fatalf("[checkArgsDDBG] argument 'Haplotype': %v set error\n ", c.Flag("Haplotype").String())
	}
	opt.Debug, ok = c.Flag("Debug").Get().(bool)
	if !ok {
		log.Fatalf("[checkArgsDDBG] argument 'Debug': %v set error\n ", c.Flag("Debug").String())
	}

	opt.ONTFn = c.Flag("LongReadFile").String()

	succ = true
	return opt, succ
}

func GetConstructDBGSampleSize(edgesArr []DBGEdge, winSize, seedLen, kmerlen int) (size int) {
	var e *DBGEdge
	var edgeNum int
	for i := range edgesArr {
		e = &edgesArr[i]
		if e.GetDeleteFlag() > 0 || e.ID < 2 || e.GetBubbleRepeatFlag() > 0 {
			continue
		}
		if IsUniqueEdge(e) || len(e.Ks) > kmerlen*2-MinEdgeLenDiff {

		} else {
			continue
		}
		edgeNum++
		//fmt.Printf("[GetCuckoofilterDBGSampleSize] e: %v\n", e)
		size += ((len(e.Ks) - seedLen + winSize - 1) / winSize) + 1
		/*if el < 2*maxNGSReadLen-kmerLen { // get whole length sliding windowns
			cfSize += ((el - kmerLen + winSize - 1) / winSize) + 1
		} else { // just sample two ends
			cfSize += ((maxNGSReadLen-kmerLen+winSize-1)/winSize + 1) * 2
		}*/
	}
	fmt.Printf("[GetConstructDBGSampleSize]need Construct DBG edge number:%d\n", edgeNum)
	return
}

type ReadsBucketInfo struct {
	Cs                               []byte
	ReadsArr                         []ReadInfo
	KmerSortArr                      []KI
	RbytesPool, RIArrPool, KIArrPool *sync.Pool
}

func CompressToUint64(ks []byte) (cb uint64) {
	for _, b := range ks {
		cb <<= NumBitsInBase
		cb |= uint64(b)
	}
	return
}

func CompressToUint16(ks []byte) (cb uint16) {
	for _, b := range ks {
		cb <<= NumBitsInBase
		cb |= uint16(b)
	}
	return
}

func GetKmerInfo(cs uint64, ID uint32, pos uint64, strand bool) (ki KI) {
	bs := cs
	bs = bs<<seqPositionBitNum | pos
	bs <<= 1
	if strand == PLUS {
		bs |= 0x1
	}
	ki.SeedInfo = SeedInfo(bs)
	ki.ID = ID
	//fmt.Printf("[GetKmerInfo]ki.SeedInfo:%b\n", ki.GetKmer())
	return
}

/*func FindMinKmerID(ka []KmerID) (min KmerID, idx int) {
	min = ka[0]
	idx = 0
	for i, ki := range ka[1:] {
		if min.GetKmer() > ki.GetKmer() {
			min = ki
			idx = i + 1
		}
	}
	return
}*/

func FindMinKI(ka []KI) (min KI, idx int) {
	min = ka[0]
	idx = 0
	for i := 1; i < len(ka); i++ {
		if min.SeedInfo.GetKmer() > ka[i].SeedInfo.GetKmer() {
			min = ka[i]
			idx = i
		}
	}
	return
}

type KIArr []KI

func (arr KIArr) Len() int {
	return len(arr)
}

func (arr KIArr) Less(i, j int) bool {
	return arr[i].GetKmer() < arr[j].GetKmer()
}
func (arr KIArr) Swap(i, j int) {
	arr[i], arr[j] = arr[j], arr[i]
}

func GetSeqMiniKmerInfoArr(ks []byte, buf, kmerInfoArr []KI, ID uint32, SeedLen, WinSize int) []KI {
	if len(ks) < SeedLen+WinSize {
		var t []KI
		return t
	}
	buf = buf[:cap(buf)]
	kmerInfoArr = kmerInfoArr[:0]
	//kmerInfoArr := make([]KI, 0, len(ks)/WinSize*2)
	//rs = GetReverseCompByteArr2(ks, rs)
	//fmt.Printf("[GetSeqMiniKmerInfoArr]eID:%d seq:%v\n\trseq:%v\n", ID, ks, rs)
	sl := len(ks)
	//buf := make([]KI, bufSize)
	idx := 0
	last := -1
	bs := CompressToUint64(ks[:SeedLen])
	rs := GetReverseCompletUint64(bs, uint32(SeedLen))
	//fmt.Printf("[GetSeqMiniKmerInfoArr]bs:%X\n", bs)
	//fmt.Printf("[GetSeqMiniKmerInfoArr]rs:%X\n", rs)
	ki := GetKmerInfo(bs, ID, 0, PLUS)
	rki := GetKmerInfo(rs, ID, 0, MINUS)
	if bs < rs {
		buf[idx] = ki
	} else {
		buf[idx] = rki
	}
	idx++
	//fmt.Printf("[GetSeqMiniKmerInfos]i:%d ks:%v\n\tki.SeedInfo:%b rki.SeedInfo:%b\n", i, ks[i:i+SeedLen], ki.GetKmer(), rki.GetKmer())
	//ki := GetKmerInfo(ks[WinSize-1:WinSize-1+SeedLen], ID, uint32(WinSize-1), true)
	//rki := GetKmerInfo(rs[sl-(WinSize-1)-SeedLen:sl-(WinSize-1)], ID, uint32(WinSize-1), false)
	BITNUM := (uint64(SeedLen) - 1) * NumBitsInBase
	MUSK := (uint64(1) << (uint64(SeedLen) * NumBitsInBase)) - 1
	for i := SeedLen; i < sl; i++ {
		bs <<= NumBitsInBase
		bs |= uint64(ks[i])
		bs &= MUSK //rki.GetKmer() &= MUSK
		ki = GetKmerInfo(bs, ID, uint64(i-SeedLen-1), PLUS)

		rs >>= NumBitsInBase
		//bt = uint64(rs[sl-i-SeedLen])
		//fmt.Printf("[GetSeqMiniKmerInfos]bt:%b\n", bt<<BITNUM)
		rs |= (uint64(BntRev[ks[i]]) << BITNUM)
		rki = GetKmerInfo(rs, ID, uint64(i-(SeedLen-1)), MINUS)
		if bs < rs {
			buf[idx] = ki
		} else {
			buf[idx] = rki
		}
		//fmt.Printf("[GetSeqMiniKmerInfoArr]bs:%X\n", bs)
		//fmt.Printf("[GetSeqMiniKmerInfoArr]rs:%X\n", rs)
		idx++

		//fmt.Printf("[GetSeqMiniKmerInfos]i:%d rb:%b ki.SeedInfo:%b rki.SeedInfo:%b\n", i, rs[sl-i-SeedLen], ki.GetKmer(), rki.GetKmer())
		if idx < WinSize {
			continue
		} else if last >= idx-WinSize {
			if buf[idx-1].GetKmer() < buf[last].GetKmer() {
				kmerInfoArr = append(kmerInfoArr, buf[idx-1])
				last = idx - 1
			}
		} else {
			min, x := FindMinKI(buf[idx-WinSize : idx])
			j := idx - WinSize + x
			kmerInfoArr = append(kmerInfoArr, min)
			last = j
		}

		//fmt.Printf("[GetSeqMiniKmerInfos]last: %v, j: %v, min: %v\n", last, j, min)

		// adjust buf
		if idx == len(buf) {
			copy(buf[0:], buf[last:])
			idx -= last
			last = 0
		}
	}
	return kmerInfoArr
}

func GetMinSyncmerMinimizerIdx(buf []uint16, start, end int) int {
	min := buf[start]
	idx := start
	for i := start + 1; i < end; i++ {
		if buf[i] < min {
			min = buf[i]
			idx = i
		}
	}
	return idx
}

func GetSeqSyncmerArr(ks []byte, buf []uint16, kmerInfoArr []KI, ID uint32, SeedLen int) []KI {
	kmerInfoArr = kmerInfoArr[:0]
	if len(ks) < SeedLen {
		return kmerInfoArr
	}
	buf = buf[:cap(buf)]
	if len(buf) < SeedLen*2 {
		log.Fatalf("[GetSeqSyncmerArr] len(buf):%d < SeedLen*2:%d\n", len(buf), SeedLen*2)
	}
	//kmerInfoArr := make([]KI, 0, len(ks)/WinSize*2)
	//rs = GetReverseCompByteArr2(ks, rs)
	//fmt.Printf("[GetSeqSyncmerArr]eID:%d seq:%v\n", ID, ks)
	//buf := make([]KI, bufSize)
	idx := 0
	bufLeastIdx := 0
	sm := CompressToUint16(ks[:SyncmerMinimizerLen])
	rsm := GetReverseCompletUint16(sm, SyncmerMinimizerLen)
	if sm < rsm {
		buf[idx] = sm
	} else {
		buf[idx] = rsm
	}
	minsm := buf[idx]
	minIdx := idx
	idx++
	for i := SyncmerMinimizerLen; i < SeedLen; i++ {
		sm <<= NumBitsInBase
		sm &= SyncmerMinimizerMUSK
		sm |= uint16(ks[i])
		rsm >>= NumBitsInBase
		rsm |= (uint16(BntRev[ks[i]]) << uint16(SyncmerMinimizerBITNUM))

		if sm < rsm {
			buf[idx] = sm
		} else {
			buf[idx] = rsm
		}
		if buf[idx] < minsm {
			minsm = buf[idx]
			minIdx = idx
		}
		//fmt.Printf("[GetSeqSyncmerArr]buf[%d]:%b minIdx:%d\n", idx, buf[idx], minIdx)
		idx++
	}
	bs := CompressToUint64(ks[:SeedLen])
	rbs := GetReverseCompletUint64(bs, uint32(SeedLen))
	if minIdx == SyncmerPos1 || minIdx == SyncmerPos2 {
		if bs < rbs {
			kmerInfoArr = append(kmerInfoArr, GetKmerInfo(bs, ID, 0, PLUS))
		} else {
			kmerInfoArr = append(kmerInfoArr, GetKmerInfo(rbs, ID, 0, MINUS))
		}
	}
	//fmt.Printf("[GetSeqMiniKmerInfoArr]rs:%X\n", rs)
	//fmt.Printf("[GetSeqMiniKmerInfos]i:%d ks:%v\n\tki.SeedInfo:%b rki.SeedInfo:%b\n", i, ks[i:i+SeedLen], ki.GetKmer(), rki.GetKmer())
	//ki := GetKmerInfo(ks[WinSize-1:WinSize-1+SeedLen], ID, uint32(WinSize-1), true)
	//rki := GetKmerInfo(rs[sl-(WinSize-1)-SeedLen:sl-(WinSize-1)], ID, uint32(WinSize-1), false)
	sl := len(ks)
	for i := SeedLen; i < sl; i++ {
		bs <<= NumBitsInBase
		bs &= SeedMUSK //rki.GetKmer() &= MUSK
		bs |= uint64(ks[i])
		rbs >>= NumBitsInBase
		rbs |= (uint64(BntRev[ks[i]]) << SeedBITNUM)

		sm <<= NumBitsInBase
		sm &= SyncmerMinimizerMUSK
		sm |= uint16(ks[i])
		rsm >>= NumBitsInBase
		rsm |= (uint16(BntRev[ks[i]]) << uint16(SyncmerMinimizerBITNUM))

		if sm < rsm {
			buf[idx] = sm
		} else {
			buf[idx] = rsm
		}
		if buf[idx] < minsm {
			minsm = buf[idx]
			minIdx = idx
		} else if bufLeastIdx+minIdx < i+1-SeedLen {
			minIdx = GetMinSyncmerMinimizerIdx(buf, idx+SyncmerMinimizerLen-SeedLen, idx+1)
			minsm = buf[minIdx]
		}
		//fmt.Printf("[GetSeqSyncmerArr]buf[%d]:%b minIdx:%d\n", idx, buf[idx], minIdx)
		idx++
		//rki = GetKmerInfo(rbs, ID, uint64(i-(SeedLen-1)), MINUS)
		if bufLeastIdx+minIdx == i-(SeedLen-1)+SyncmerPos1 || bufLeastIdx+minIdx == i-(SeedLen-1)+SyncmerPos2 {
			if bs < rbs {
				kmerInfoArr = append(kmerInfoArr, GetKmerInfo(bs, ID, uint64(i-(SeedLen-1)), PLUS))
			} else {
				kmerInfoArr = append(kmerInfoArr, GetKmerInfo(rbs, ID, uint64(i-(SeedLen-1)), MINUS))
			}
			//fmt.Printf("[GetSeqSyncmerArr]add KI pos:%d minIdx:%d\n", i-(SeedLen-1), bufLeastIdx+minIdx)
		}
		// adjust buf
		if idx == len(buf) {
			bufLeastIdx += minIdx
			copy(buf[0:], buf[minIdx:])
			idx = len(buf) - minIdx
			minIdx = 0
		}
	}
	return kmerInfoArr
}

const MinmerLen = 7
const MinmerMUSK = (uint16(1) << (uint16(MinmerLen) * NumBitsInBase)) - 1
const MinmerBITLen = uint64(MinmerLen) * NumBitsInBase
const MinmerBITNUM = (uint64(MinmerLen) - 1) * NumBitsInBase

func GetSeqMinmerArr(ks []byte, buf []uint16, kmerInfoArr []KI, ID uint32, SeedLen int) []KI {
	if len(ks) < SeedLen {
		return kmerInfoArr[:0]
	}
	buf = buf[:cap(buf)]
	if len(buf) < SeedLen*2 {
		log.Fatalf("[GetSeqMinmerArr] len(buf):%d < SeedLen*2:%d\n", len(buf), SeedLen*2)
	}
	kmerInfoArr = kmerInfoArr[:0]
	//kmerInfoArr := make([]KI, 0, len(ks)/WinSize*2)
	//rs = GetReverseCompByteArr2(ks, rs)
	//fmt.Printf("[GetSeqSyncmerArr]eID:%d seq:%v\n", ID, ks)
	//buf := make([]KI, bufSize)
	idx := 0
	bufLeastIdx := 0
	sm := CompressToUint16(ks[:MinmerLen])
	rsm := GetReverseCompletUint16(sm, MinmerLen)
	if sm < rsm {
		buf[idx] = sm
	} else {
		buf[idx] = rsm
	}
	minsm := buf[idx]
	minIdx := idx
	idx++
	for i := MinmerLen; i < SeedLen; i++ {
		sm <<= NumBitsInBase
		sm &= MinmerMUSK
		sm |= uint16(ks[i])
		rsm >>= NumBitsInBase
		rsm |= (uint16(BntRev[ks[i]]) << uint16(MinmerBITNUM))

		if sm < rsm {
			buf[idx] = sm
		} else {
			buf[idx] = rsm
		}
		if buf[idx] < minsm {
			minsm = buf[idx]
			minIdx = idx
		}
		//fmt.Printf("[GetSeqSyncmerArr]buf[%d]:%b minIdx:%d\n", idx, buf[idx], minIdx)
		idx++
	}
	var bs uint64
	lastMinsm := minsm
	lastIdx := minIdx
	lastIdx2, lastIdx3, lastIdx4 := -1, -1, -1
	//fmt.Printf("[GetSeqMinmerArr]minIdx:%d minsm:%b\n", bufLeastIdx+minIdx, minsm)
	//fmt.Printf("[GetSeqMiniKmerInfos]i:%d ks:%v\n\tki.SeedInfo:%b rki.SeedInfo:%b\n", i, ks[i:i+SeedLen], ki.GetKmer(), rki.GetKmer())
	//ki := GetKmerInfo(ks[WinSize-1:WinSize-1+SeedLen], ID, uint32(WinSize-1), true)
	//rki := GetKmerInfo(rs[sl-(WinSize-1)-SeedLen:sl-(WinSize-1)], ID, uint32(WinSize-1), false)
	sl := len(ks)
	for i := SeedLen; i < sl; i++ {
		sm <<= NumBitsInBase
		sm &= MinmerMUSK
		sm |= uint16(ks[i])
		rsm >>= NumBitsInBase
		rsm |= (uint16(BntRev[ks[i]]) << uint16(MinmerBITNUM))

		if sm < rsm {
			buf[idx] = sm
		} else {
			buf[idx] = rsm
		}
		var added bool
		if bufLeastIdx+minIdx < i+1-SeedLen {
			minIdx = GetMinSyncmerMinimizerIdx(buf, idx+MinmerLen-SeedLen, idx+1)
			minsm = buf[minIdx]
			added = true
		} else if buf[idx] < minsm {
			minsm = buf[idx]
			minIdx = idx
			added = true
		}
		//fmt.Printf("[GetSeqSyncmerArr]buf[%d]:%b minIdx:%d\n", idx, buf[idx], minIdx)
		idx++
		//rki = GetKmerInfo(rbs, ID, uint64(i-(SeedLen-1)), MINUS)
		if added {
			//fmt.Printf("[GetSeqMinmerArr]minIdx:%d minsm:%b\n", bufLeastIdx+minIdx, minsm)
			if lastIdx+3 >= bufLeastIdx+minIdx {
				lastMinsm = minsm
				lastIdx = bufLeastIdx + minIdx
			} else {
				bs <<= MinmerBITLen
				bs |= uint64(lastMinsm)
				lastMinsm = minsm
				lastIdx4 = lastIdx3
				lastIdx3 = lastIdx2
				lastIdx2 = lastIdx
				lastIdx = bufLeastIdx + minIdx
				if lastIdx4 >= 0 {
					bs &= SeedMUSK
					kmerInfoArr = append(kmerInfoArr, GetKmerInfo(bs, ID, uint64(lastIdx4), PLUS))
					//fmt.Printf("[GetSeqMinmerArr]kaIdx:%d bs:%b\n", lastIdx4, bs)
				}
			}
			//fmt.Printf("[GetSeqSyncmerArr]add KI pos:%d minIdx:%d\n", i-(SeedLen-1), bufLeastIdx+minIdx)
		}
		// adjust buf
		if idx == len(buf) {
			bufLeastIdx += minIdx
			copy(buf[0:], buf[minIdx:])
			idx = len(buf) - minIdx
			minIdx = 0
		}
	}
	if lastIdx3 >= 0 {
		bs <<= MinmerBITLen
		bs |= uint64(lastMinsm)
		bs &= SeedMUSK
		kmerInfoArr = append(kmerInfoArr, GetKmerInfo(bs, ID, uint64(lastIdx3), PLUS))
		//fmt.Printf("[GetSeqMinmerArr]kaIdx:%d bs:%b\n", lastIdx3, bs)
	}
	return kmerInfoArr
}

func RadixSortKIArr3Bit(kmerArr []KI, bulkArr [][]KI, SeedLen, BitNum, NumBulk uint64) ([]KI, [][]KI) {
	if len(kmerArr) < 2 {
		return kmerArr, bulkArr
	}
	if len(kmerArr) < 2000 {
		sort.Sort(KIArr(kmerArr))
		return kmerArr, bulkArr
	}
	MUSK := uint64((1 << BitNum) - 1)
	overflowNum := 0
	countArr := make([]int, NumBulk)
	sz := (SeedLen*NumBitsInBase + BitNum - 1) / BitNum
	bulkArr[0] = kmerArr
	for i := uint64(0); i < sz; i++ {
		// clean countArr
		for j := range countArr {
			countArr[j] = 0
		}
		shiftNum := i * BitNum
		for _, ki := range kmerArr {
			idx := (ki.GetKmer() >> shiftNum) & MUSK
			if countArr[idx] >= len(bulkArr[idx]) {
				fmt.Printf("[RadixSortKIArr3Bit] overflow idx:%d\n", idx)
				overflowNum++
				extSize := countArr[idx] * 6 / 5
				if countArr[idx] < 100 {
					extSize = 2 * countArr[idx]
				}
				na := make([]KI, extSize)
				copy(na, bulkArr[idx])
				bulkArr[idx] = na
			}
			bulkArr[idx][countArr[idx]] = ki
			countArr[idx]++
		}
		// collect
		ct := countArr[0]
		//fmt.Printf("[RadixSortKIArr3Bit]i:%d countArr[0]:%d\n", i, countArr[0])
		for j := uint64(1); j < NumBulk; j++ {
			//fmt.Printf("[RadixSortKIArr3Bit]i:%d countArr[%d]:%d\n", i, j, countArr[j])
			for x := 0; x < countArr[j]; x++ {
				kmerArr[ct] = bulkArr[j][x]
				ct++
			}
		}
		if ct != len(kmerArr) {
			log.Fatalf("[RadixSortKIArr3Bit] ct:%d != len(kmerArr):%d\n", ct, len(kmerArr))
		}
	}
	fmt.Printf("[RadixSortKIArr3Bit] overflowNum:%d\n", overflowNum)
	return kmerArr, bulkArr
}

func SortLongReadsKmer(rc <-chan RIPool, wc chan<- ReadsBucketInfo, SeedLen, WinSize, RISeqSize int, rbytesPool, RIArrPool *sync.Pool) {
	//var addReadKmerLen int64
	//rs := make([]byte, 10000)
	KIArrPool := sync.Pool{New: func() interface{} {
		arr := make([]KI, RISeqSize*2/WinSize)
		return arr
	}}
	bufSize := SeedLen * 2
	buf := make([]uint16, bufSize)
	NumBulk := uint64(1) << BitNum
	bulkArr := make([][]KI, NumBulk)
	bukSize := RISeqSize / int(NumBulk)
	for i := 1; i < len(bulkArr); i++ {
		if i < len(bulkArr) {
			bulkArr[i] = make([]KI, bukSize*9/5)
		} else {
			bulkArr[i] = make([]KI, bukSize*8/5)
		}
	}
	ka := make([]KI, 10000)
	for {
		ra, ok := <-rc
		if !ok {
			break
		}
		//tl := GetReadArrSeqLen(ra)
		//var kmerArr []KI
		var sl int
		kmerArr := KIArrPool.Get().([]KI)
		kmerArr = kmerArr[:0]
		for _, ri := range ra.RIArr {
			if len(ri.Seq) > MaxAllowSeqLen {
				log.Fatalf("[SortLongReadsKmer]read ID:%d seqLen:%d > MaxAllowSeqLen:%d\n", ri.ID, len(ri.Seq), MaxAllowSeqLen)
			}
			sl += len(ri.Seq)
			//ka = GetSeqMiniKmerInfoArr(ri.Seq, buf, ka, uint32(ri.ID), SeedLen, WinSize)
			ka = GetSeqSyncmerArr(ri.Seq, buf, ka, uint32(ri.ID), SeedLen)
			//ka = GetSeqMinmerArr(ri.Seq, buf, ka, uint32(ri.ID), SeedLen)
			kmerArr = append(kmerArr, ka...)
		}
		kmerArr, bulkArr = RadixSortKIArr3Bit(kmerArr, bulkArr, uint64(SeedLen), BitNum, NumBulk)
		//sort.Sort(KIArr(kmerArr))
		var rbi ReadsBucketInfo
		rbi.Cs = ra.Cs
		rbi.RbytesPool = rbytesPool
		rbi.RIArrPool = RIArrPool
		rbi.KIArrPool = &KIArrPool
		rbi.ReadsArr = ra.RIArr
		rbi.KmerSortArr = kmerArr
		wc <- rbi
		fmt.Printf("[SortLongReadsKmer]finished sort reads num:%d seqLen:%d miniKmer num:%d cap(kmerArr):%d\n", len(ra.RIArr), sl, len(kmerArr), cap(kmerArr))
	}
}

func LoadLongReadsAndSort(fn string, wc chan<- ReadsBucketInfo, SeedLen, WinSize int, buckSize int) {
	avgLongReadLen := 10000
	rbytesPool := sync.Pool{New: func() interface{} {
		bytes := make([]byte, buckSize)
		return bytes
	}}

	size := buckSize / avgLongReadLen
	RIArrPool := sync.Pool{New: func() interface{} {
		arr := make([]ReadInfo, size)
		return arr
	}}
	cs := make(chan []byte, 1)
	RIPoolChan := make(chan RIPool, 1)
	go ReadZstdFile(fn, &rbytesPool, cs)
	go GetReadInfoBucket(fn, cs, &RIArrPool, RIPoolChan, false)
	SortLongReadsKmer(RIPoolChan, wc, SeedLen, WinSize, buckSize, &rbytesPool, &RIArrPool)
}

func GetLongReadsFile(cfgInfo CfgInfo) (fnArr []string) {
	for _, lib := range cfgInfo.Libs {
		if lib.SeqProfile == 3 || lib.SeqProfile == 4 {
			fnArr = append(fnArr, lib.FnName...)
		}
	}
	return
}

func ParaLoadLongReadsAndSort(cfgInfo CfgInfo, cs chan<- ReadsBucketInfo, SeedLen, WinSize, BuckSize int) {
	fnArr := GetLongReadsFile(cfgInfo)
	for _, fn := range fnArr {
		LoadLongReadsAndSort(fn, cs, SeedLen, WinSize, BuckSize)
	}
}

func IsEdgeNoDB(e *DBGEdge) bool {
	if e.GetBubbleRepeatFlag() > 0 {
		return true
	}
	if IsUniqueEdge(e) || len(e.Ks) > MinEdgeLen {
		return false
	} else {
		return true
	}
}

type IDInfo struct {
	ID   uint32 // edgeID or reads ID
	Info uint32 // the first high 31-bit for Position,and last bit for strand info
}

func (k *IDInfo) GetPosition() uint32 {
	return (k.Info >> 1)
}
func (k *IDInfo) SetPosition(p uint32) {
	k.Info = (p << 1) | (k.Info & 0x1)
}

func (k *IDInfo) GetStrand() bool {
	if (k.Info & 0x1) > 0 {
		return true
	} else {
		return false
	}
}

func (k *IDInfo) SetStrand(s bool) {
	if s == PLUS {
		k.Info |= 0x1
	} else {
		k.Info = k.Info & (((1 << 31) - 1) << 1)
	}
}

type IDInfoArr []IDInfo

func (arr IDInfoArr) Len() int {
	return len(arr)
}

func (arr IDInfoArr) Less(i, j int) bool {
	return arr[i].ID < arr[j].ID
}
func (arr IDInfoArr) Swap(i, j int) {
	arr[i], arr[j] = arr[j], arr[i]
}

type DBGKmerInfoPac struct {
	KmerinfoArr           []KI                // KI arr, if found high freq kmer, just store one kmer that ID == math.MaxUint32
	HighFreqMap           map[uint64][]IDInfo // key is KI.Kmer, value is map[uint32]uint32, key is KI.ID, value is KI.Info arr
	LowFreqPercentEdgeMap map[uint64][]uint32 // low normal kmer coverage edge info map, key is KI.Kmer, value is arr of eID
}

func SortDBGEdgesKmer(dbgCS chan<- DBGKmerInfoPac, edgesArr []DBGEdge, SeedLen int, opt optionsDDBG, HighFreqKmerMin, HighFreqKmerMax, MinimerKmerNum int) {
	var dbgKmerInfoPac DBGKmerInfoPac
	var deleteKmerNum int
	var edgeNum, edgeLenSum int
	kmerlen := opt.Kmer
	MinEdgeLen = kmerlen*2 - MinEdgeLenDiff
	kmerArr := make([]KI, 0, MinimerKmerNum)
	bufSize := SeedLen * 2
	buf := make([]uint16, bufSize)
	ka := make([]KI, 1000)
	//rs := make([]byte, 10000)
	for i := 0; i < len(edgesArr); i++ {
		e := &edgesArr[i]
		if e.ID < 2 || e.GetDeleteFlag() > 0 {
			continue
		}
		if IsEdgeNoDB(e) {
			continue
		}
		edgeNum++
		edgeLenSum += e.GetSeqLen()
		if e.GetSeqLen() > MaxAllowSeqLen {
			log.Fatalf("[SortDBGEdgesKmer]eID:%d seqLen:%d > MaxAllowSeqLen:%d\n", e.ID, e.GetSeqLen(), MaxAllowSeqLen)
		}
		//ka = GetSeqMiniKmerInfoArr(e.Ks, buf, ka, e.ID, SeedLen, opt.WinSize)
		ka = GetSeqSyncmerArr(e.Ks, buf, ka, e.ID, SeedLen)
		//ka = GetSeqMinmerArr(e.Ks, buf, ka, e.ID, SeedLen)
		//fmt.Printf("[SortDBGEdgesKmer]eID:%d el:%d len(ka):%d\n", e.ID, len(e.Ks), len(ka))
		//sort.Sort(KIArr(ka))
		//e.SeedInfoArr = make([]SeedInfo, len(ka))
		//for x, ki := range ka {
		//	e.SeedInfoArr[x] = ki.SeedInfo
		//}
		//edgesArr[e.ID].SeedInfoArr = e.SeedInfoArr
		//var ka []KI
		//PrintKIArr(ka)
		kmerArr = append(kmerArr, ka...)
		//freqMat[i] = make([]FreqInfo, len(ka))
		//for x, ki := range ka {
		//	freqMat[i][x].Pos = uint32(ki.GetPosition())
		//}
		//for x, ki := range ka {
		//	fmt.Printf("[SortDBGEdgesKmer]%d\t Pos:%d Plus:%t kmer:%b\n", x, ki.GetPosition(), ki.GetStrand(), ki.GetKmer())
		//}
	}
	fmt.Printf("[SortDBGEdgesKmer]cap(kmerArr):%d len(kmerArr):%d\n", MinimerKmerNum, len(kmerArr))

	NumBulk := uint64(1) << BitNum
	bulkArr := make([][]KI, NumBulk)
	//bulkArr[0] = make([]KI, len(kmerArr)/int(NumBulk)*2)
	fmt.Printf("[SortDBGEdgesKmer]buk Size:%d\n", len(kmerArr)/int(NumBulk))
	bukSize := len(kmerArr) / int(NumBulk)
	for i := 1; i < len(bulkArr); i++ {
		if i < len(bulkArr)/2 {
			bulkArr[i] = make([]KI, bukSize*9/5)
		} else {
			bulkArr[i] = make([]KI, bukSize*8/5)
		}
	}
	kmerArr, _ = RadixSortKIArr3Bit(kmerArr, bulkArr, uint64(SeedLen), BitNum, NumBulk)
	totalKALen := len(kmerArr)
	//sort.Sort(KIArr(kmerArr))
	deleteCount := 0
	UniqueNum := 0
	highFreqKmerMap := make(map[uint64][]IDInfo, 2000)

	//sortKmerArr := make([]KI, count)
	//lowPercentMinimersEdgesMap := make(map[uint64][]uint32, 2000) // key is kmer, value eID array
	//lowFreqNumArr := make([]uint8, len(edgesArr))
	count := 0
	skaIdx := 0
	//highFreqInDBNum := 0
	for i := 0; i < len(kmerArr); {
		k1 := kmerArr[i].GetKmer()
		j := i + 1
		for ; j < len(kmerArr); j++ {
			k2 := kmerArr[j].GetKmer()
			if k1 != k2 {
				break
			}
		}
		if j-i <= 3 {
			UniqueNum += j - i
		}
		if j-i < HighFreqKmerMin {
			for x := i; x < j; x++ {
				kmerArr[skaIdx] = kmerArr[x]
				skaIdx++
			}
			count += j - i
		} else if j-i < HighFreqKmerMax {
			arr := make([]IDInfo, j-i)
			var idi IDInfo
			for x := i; x < j; x++ {
				idi.ID = kmerArr[x].ID
				idi.SetPosition(uint32(kmerArr[x].GetPosition()))
				idi.SetStrand(kmerArr[x].GetStrand())
				arr[x-i] = idi
			}
			sort.Sort(IDInfoArr(arr))
			highFreqKmerMap[kmerArr[j-1].GetKmer()] = arr
			kmerArr[j-1].ID = math.MaxUint32
			kmerArr[skaIdx] = kmerArr[j-1]
			skaIdx++
			count++
		} else {
			deleteKmerNum += (j - i)
			deleteCount++
		}
		i = j
	}

	if skaIdx != count {
		fmt.Printf("[SortDBGEdgesKmer] skaIdx:%v != count:%d\n", skaIdx, count)
	}
	kmerArr = kmerArr[:skaIdx]

	miniNumArr := make([]uint16, len(edgesArr))
	for i := range kmerArr {
		eID := kmerArr[i].ID
		if eID == math.MaxUint32 {
			continue
		}
		if miniNumArr[eID] < math.MaxUint16 {
			miniNumArr[eID]++
		}
	}

	// set edgesArr miniNum
	var LitterMinMinimerEdgeNum int
	for i := 0; i < len(edgesArr); i++ {
		e := &edgesArr[i]
		if e.ID < 2 || e.GetDeleteFlag() > 0 {
			continue
		}
		if IsEdgeNoDB(e) {
			continue
		}
		edgesArr[i].CovD = miniNumArr[i]
		if int(miniNumArr[i]) < MinMinimerNum {
			LitterMinMinimerEdgeNum++
		}
		//if int(miniNumArr[i]) >= MinMinimerNum {
		//	edgesArr[i].SeedInfoArr = nil
		//}

		//if int(miniNumArr[i]) < MinMinimerNum {
		//fmt.Printf("[SortDBGEdgesKmer]eID:%d eLen:%d miniNum:%d repeat:%t\n", e.ID, e.GetSeqLen(), miniNumArr[i], IsRepeatEdge(e))
		//}

	}

	totalHighFreq := 0
	for _, v := range highFreqKmerMap {
		totalHighFreq += len(v)
	}

	//fmt.Printf("[SortDBGEdgesKmer]highFreqInDBNum:%d percent:%.3f\n", highFreqInDBNum, float64(highFreqInDBNum)/float64(totalKALen))
	fmt.Printf("[SortDBGEdgesKmer]edgeNum:%d edgeLenSum:%d Seed Number <%d edgeNum:%d\n", edgeNum, edgeLenSum, MinMinimerNum, LitterMinMinimerEdgeNum)
	fmt.Printf("[SortDBGEdgesKmer]deleteCount:%d deleteKmerNum:%d percent:%.4f\n", deleteCount, deleteKmerNum, float64(deleteKmerNum)/float64(totalKALen))
	fmt.Printf("[SortDBGEdgesKmer]HighFreqKmerMin:%d MinimerKmerNum:%d kmerArr len:%d count:%d \n", HighFreqKmerMin, MinimerKmerNum, len(kmerArr), count)
	fmt.Printf("[SortDBGEdgesKmer]UniqueNum:%d Unique percent:%v high Freq kmer num:%d avg high freq:%d high freq percent:%v\n", UniqueNum, float32(UniqueNum)/float32(totalKALen), len(highFreqKmerMap), totalHighFreq/(len(highFreqKmerMap)+1), float32(totalHighFreq)/float32(totalKALen))

	//dbgKmerInfoPac.KmerinfoArr = make([]KI, len(kmerArr))
	//copy(dbgKmerInfoPac.KmerinfoArr, kmerArr)
	dbgKmerInfoPac.KmerinfoArr = kmerArr
	dbgKmerInfoPac.HighFreqMap = highFreqKmerMap
	//dbgKmerInfoPac.LowFreqPercentEdgeMap = lowPercentMinimersEdgesMap
	dbgCS <- dbgKmerInfoPac
	close(dbgCS)
	//time.Sleep(time.Second)
}

func WriteLongPathToFile(wc chan string, pathfn, ONTMapInfofn string, numCPU int) {
	pathfp, err := os.Create(pathfn)
	if err != nil {
		log.Fatalf("[WriteLongPathToFile] file %s create error, err: %v\n", pathfn, err)
	}
	defer pathfp.Close()
	zfp, err1 := zstd.NewWriter(pathfp, zstd.WithEncoderCRC(false), zstd.WithEncoderConcurrency(1), zstd.WithEncoderLevel(1))
	defer zfp.Close()
	if err1 != nil {
		log.Fatalf("[WriteLongPathToFile]open write file:%s err:%v\n", pathfn, err1)
	}
	var finishNum int
	mapRecordNum := 0
	for {
		ms := <-wc
		if len(ms) == 0 {
			finishNum++
			if finishNum == numCPU {
				break
			}
			continue
		}

		fmt.Fprintf(zfp, "%s\n", ms)
		mapRecordNum++
	}
	mapfp, err := os.Create(ONTMapInfofn)
	if err != nil {
		log.Fatalf("[WriteLongPathToFile] file %s create error, err: %v\n", ONTMapInfofn, err)
	}
	defer mapfp.Close()
	fmt.Fprintf(mapfp, "MapRecordSize:%d\n", mapRecordNum)
}

func LoadPathFromFile(pathfn string, cs chan<- MapingInfo) {
	bufSize := WindowSize
	rbytesPool := sync.Pool{New: func() interface{} {
		bytes := make([]byte, bufSize)
		return bytes
	}}
	fncs := make(chan []byte, 2)
	go ReadZstdFile(pathfn, &rbytesPool, fncs)
	buf := make([]byte, bufSize+(1<<12))
	buf = buf[:0]
	//lb := make([]byte, 100)
	var lb []byte
	pos := 0
	var suc bool
	for {
		lb, pos, suc = GetNextLine(buf, pos)
		if !suc {
			nb, ok := <-fncs
			if !ok {
				break
			}
			if pos < len(buf) {
				copy(buf[0:], buf[pos:])
			}
			buf = buf[:len(buf)-pos]
			pos = 0
			buf = append(buf, nb...)
			rbytesPool.Put(nb)
			continue
		}
		var mi MapingInfo
		colList := SplitEIDArr(lb, '\t')
		id, err := strconv.Atoi(string(colList[0]))
		if err != nil {
			log.Fatalf("[LoadPathFromFile]%s transform read id error\n", colList[0])
		}
		mi.ID = uint32(id)
		mi.Anotition = string(colList[6])
		var pathMat [][][]uint32
		var path []uint32
		//flist := strings.Fields(string(lb))
		//eIDArr := make([]uint32, len(colList))
		list := SplitEIDArr(colList[5], '>')
		for _, ele := range list {
			bubbleArr := SplitEIDArr(ele, '|')
			if len(bubbleArr) == 1 {
				bs := bubbleArr[0]
				eID, err := strconv.Atoi(string(bs))
				if err != nil {
					log.Fatalf("[LoadPathFromFile]%s transform eID error\n", bs)
				}
				path = append(path, uint32(eID))
			} else {
				pathMat = append(pathMat, [][]uint32{path})
				var tmp []uint32
				path = tmp
				pathArr := make([][]uint32, len(bubbleArr))
				for x, path := range bubbleArr {
					arr := SplitEIDArr(path, '-')
					for _, bs := range arr {
						//eID, err := ByteArrInt(id)
						eID, err := strconv.Atoi(string(bs))
						if err != nil {
							log.Fatalf("[LoadPathFromFile]%s transform eID error\n", bs)
						}
						pathArr[x] = append(pathArr[x], uint32(eID))
					}
				}
				pathMat = append(pathMat, pathArr)
			}
		}
		//fmt.Printf("[LoadPathFromFile]eIDArr:%v\n", eIDArr)
		mi.PathMat = pathMat
		cs <- mi
	}
	close(cs)
	return
}

func LoadONTMapInfoFile(mapDBGInfoFn string) (mapRecordSize int) {
	mapfp, err := os.Open(mapDBGInfoFn)
	if err != nil {
		log.Fatalf("[LoadONTMapInfoFile] file %s open error, err: %v\n", mapDBGInfoFn, err)
	}
	defer mapfp.Close()
	if _, err = fmt.Fscanf(mapfp, "MapRecordSize:%d\n", &mapRecordSize); err != nil {
		log.Fatalf("[LoadONTMapInfoFile] file %s scanf error, err: %v\n", mapDBGInfoFn, err)
	}
	return
}

type SeedNumInfo struct {
	//EdgeInfo           uint32 // high 12-bits for EdgeLen, middle 10-bits for highest pos of edge, low 10-bits for lowest pos of edge
	SeedNumP, SeedNumM uint8 // SeedNumP note PLUS strand, SeedNumM note MINUS strand
	EdgeLowFreqSeedNum uint8
}

func GetEdgeLowFreqSeedNumArr(edgesArr []DBGEdge, LongEdgeMinLen int) []SeedNumInfo {
	snArr := make([]SeedNumInfo, len(edgesArr))
	for i, e := range edgesArr {
		sn := int(e.CovD)
		if sn > math.MaxUint8 {
			snArr[i].EdgeLowFreqSeedNum = math.MaxUint8
		} else {
			snArr[i].EdgeLowFreqSeedNum = uint8(sn)
		}
	}
	return snArr
}

func GetInterSectionKmerInfoB(dbgKmerInfoPac DBGKmerInfoPac, kmerInfoArrB []KI, StartID uint32, SeedLen uint16, bucketInterSecKmerInfoArr []MapingKmerInfo) []MapingKmerInfo {
	kmerInfoArrA := dbgKmerInfoPac.KmerinfoArr
	bucketInterSecKmerInfoArr = bucketInterSecKmerInfoArr[:0]
	i, j := 0, 0
	rk := kmerInfoArrB[j]
	for ; i < len(kmerInfoArrA) && j < len(kmerInfoArrB); i++ {
		//fmt.Printf("[GetInterSectionKmerInfo]kmerInfoArrB[%v]: %v, kmerInfoArrA[%v]: %v\n", i, kmerInfoArrB[i], j, kmerInfoArrA[j])
		ek := kmerInfoArrA[i]
		if ek.GetKmer() < rk.GetKmer() {
			continue
		} else if ek.GetKmer() > rk.GetKmer() {
			for j = j + 1; j < len(kmerInfoArrB); j++ {
				rk = kmerInfoArrB[j]
				if ek.GetKmer() <= rk.GetKmer() {
					break
				}
			}
		}

		if ek.GetKmer() < rk.GetKmer() {
			continue
		} else if ek.GetKmer() > rk.GetKmer() {
			break
		}
		// rk.GetKmer() == ek.GetKmer()
		var mk MapingKmerInfo
		mk.EID = ek.ID
		//if mk.EID > 748030 {
		//	fmt.Fprintf(os.Stderr, "[GetInterSectionKmerInfoB]mk.EID:%d\n", mk.EID)
		//}
		mk.SetEPosition(uint32(ek.GetPosition()))
		mk.SetEStrand(ek.GetStrand())
		if rk.ID-StartID > mapingKmerInfoMaxMapLen {
			log.Fatalf("[GetInterSectionKmerInfoB]id:%d > %d\n", rk.ID-StartID, mapingKmerInfoMaxMapLen)
		}
		mk.SetMapLen(rk.ID - StartID)
		for m := j; m < len(kmerInfoArrB); m++ {
			tk := kmerInfoArrB[m]
			if tk.GetKmer() > ek.GetKmer() {
				break
			}
			mk.SeedInfo = tk.SeedInfo
			//fmt.Printf("[GetInterSectionKmerInfo]rk: %v, tk: %v\n", rk, tk)
			//fmt.Printf("[GetInterSectionKmerInfo]edgeID: %v: ReadPos: %v, EdgePos: %v.SeedInfo: %v\n", rk.ID, tk.GetPos(), rk.GetPos(), tk.GetKmer())
			bucketInterSecKmerInfoArr = append(bucketInterSecKmerInfoArr, mk)
		}
	}

	return bucketInterSecKmerInfoArr
}

type MapingKmerInfoReadIDArr []MapingKmerInfo

func (arr MapingKmerInfoReadIDArr) Len() int {
	return len(arr)
}

func (arr MapingKmerInfoReadIDArr) Less(i, j int) bool {
	if arr[i].GetMapLen() == arr[j].GetMapLen() {
		if arr[i].GetKmer() == arr[j].GetKmer() {
			if arr[i].EID == arr[j].EID {
				return arr[i].GetEPosition() < arr[j].GetEPosition()
			} else {
				return arr[i].EID < arr[j].EID
			}
		} else {
			return arr[i].GetKmer() < arr[j].GetKmer()
		}
	} else {
		return arr[i].GetMapLen() < arr[j].GetMapLen()
	}
}

func (arr MapingKmerInfoReadIDArr) Swap(i, j int) {
	arr[i], arr[j] = arr[j], arr[i]
}

type EdgeStrand struct {
	EdgeID    uint32
	SNPercent float32
	Strand    bool
}

type EdgeStrandArr []EdgeStrand

func (arr EdgeStrandArr) Len() int {
	return len(arr)
}

func (arr EdgeStrandArr) Less(i, j int) bool {
	return arr[i].SNPercent > arr[j].SNPercent
}

func (arr EdgeStrandArr) Swap(i, j int) {
	arr[i], arr[j] = arr[j], arr[i]
}

func DeleteLowSeedNumMinimersEdge(SeedNumCountArr []SeedNumInfo, LongEdgeMinLen uint32, WinSize int, kmerlen, MQ int) ([]EdgeStrand, float32) {
	SeedNumSum := float32(0)
	//EdgeMinLen := uint32((kmerlen-1)*2) - 10
	//minNum := uint8(EdgeMinLen * 2 / WinSize * 4 / 5)
	//dif := uint8(20)
	if MQ < 10 {
		MQ = 10
	}
	MinPerCent := float32(MQ) / 100 * 0.5
	//edgeCountMap := make(map[uint32]uint32, 200)
	esArr := make([]EdgeStrand, 0, 200)
	//t := int(LongEdgeMinLen) * 2 / WinSize * 4 / 5
	//if t > math.MaxUint8 {
	//	t = math.MaxUint8
	//}
	//LongEdgeSeedNumMin := uint8(t) // LongEdgeMinLen/WinSize * 2
	//min := (kmerlen+EdgeMapKmerMin)/dif + 1
	for i, sni := range SeedNumCountArr {
		numP, numM := sni.SeedNumP, sni.SeedNumM
		if numP == 0 && numM == 0 {
			continue
		}
		num := MaxUint8(numP, numM)
		SeedNumCountArr[i].SeedNumP = 0
		SeedNumCountArr[i].SeedNumM = 0
		//num1, num2 := GetKmerNum(ct, PLUS), GetKmerNum(ct, MINUS)
		//el := sni.GetEdgeLen()
		edgeLowFreqSeedNum := sni.EdgeLowFreqSeedNum

		/*if num > 10000 {
			fmt.Printf("[DeleteLowSeedNumMinimersEdge]num1:%d num2:%d el:%d\n", num1, num2, el)
		}*/
		if float32(num) > float32(edgeLowFreqSeedNum)*MinPerCent && num > 1 {
			var es EdgeStrand
			es.EdgeID = uint32(i)
			if num == numP {
				es.Strand = PLUS
			} else {
				es.Strand = MINUS
			}
			es.SNPercent = float32(num) / float32(edgeLowFreqSeedNum)
			esArr = append(esArr, es)
			SeedNumSum += es.SNPercent
			//fmt.Printf("[DeleteLowSeedNumMinimersEdge]el:%d SeedNum:%d edgeLowFreqSeedNum:%d\n", el, num, edgeLowFreqSeedNum)
		}
		/*if el >= EdgeMinLen {
			minRegLen := uint32(edgeLowFreqSeedNum) * el / (el * 2 / uint32(WinSize)) * 2 / 3
			if num >= edgeLowFreqSeedNum/dif+1 && sni.GetHighPos()-sni.GetLowPos() > minRegLen {
				SeedNumCountArr[i].SeedNum = num
				var esn EdgeSeedNumInfo
				esn.EdgeID, esn.SeedNInfo.SeedNum, esn.SeedNInfo.EdgeLowFreqSeedNum = uint32(i), sni.SeedNum, sni.EdgeLowFreqSeedNum
				EdgeSeedNumCountArr = append(EdgeSeedNumCountArr, esn)
				b++
				longEdgeNum++
				longEdgeSeedNumSum += int(edgeLowFreqSeedNum) / int(num)
				totalNum++
				SeedNumSum += int(edgeLowFreqSeedNum) / int(num)
			} else {
				SeedNumCountArr[i].SeedNum = 0
			}
		} else {
			minRegLen := uint32(edgeLowFreqSeedNum) * el / (el * 2 / uint32(WinSize)) * 3 / 4
			if num >= edgeLowFreqSeedNum/dif+1 && uint32(num) < el && sni.GetHighPos()-sni.GetLowPos() > minRegLen {
				SeedNumCountArr[i].SeedNum = num
				var esn EdgeSeedNumInfo
				esn.EdgeID, esn.SeedNInfo.SeedNum, esn.SeedNInfo.EdgeLowFreqSeedNum = uint32(i), sni.SeedNum, sni.EdgeLowFreqSeedNum
				EdgeSeedNumCountArr = append(EdgeSeedNumCountArr, esn)
				b++
				totalNum++
				SeedNumSum += int(edgeLowFreqSeedNum) / int(num)
			} else {
				SeedNumCountArr[i].SeedNum = 0
			}
		}*/
	}
	if DebugModel {
		fmt.Printf("[DeleteLowSeedNumMinimersEdge]MQ:%d totalNum:%d\n", MQ, len(esArr))
	}
	sp := SeedNumSum / float32(len(esArr))
	MaxNum := 200
	if len(esArr) > MaxNum {
		sort.Sort(EdgeStrandArr(esArr))
		esArr = esArr[:MaxNum]
	}
	//longEdgePercent := uint16(float32(1) / (float32(longEdgeSeedNumSum) / float32(longEdgeNum)) * 0.9 * 100)
	// filter by long edge percent
	//lesnm := uint16(LongEdgeSeedNumMin)
	/*if longEdgePercent > 5 {
		count := 0
		for _, esnc := range EdgeSeedNumCountArr {
			num := uint16(esnc.SeedNInfo.SeedNum)
			edgeLowFreqSeedNum := uint16(esnc.SeedNInfo.EdgeLowFreqSeedNum)
			el := esnc.SeedNInfo.GetEdgeLen()
			if el >= LongEdgeMinLen {
				if num >= edgeLowFreqSeedNum*longEdgePercent/100 {
					EdgeSeedNumCountArr[count] = esnc
					count++
				}
			} else {
				if num >= edgeLowFreqSeedNum*longEdgePercent/100+1 {
					EdgeSeedNumCountArr[count] = esnc
					count++
				}
			}
		}
		EdgeSeedNumCountArr = EdgeSeedNumCountArr[:count]
	} else {
		longEdgePercent = 5
	}*/

	return esArr, sp
}

const Step = 10

func IsInEdgeStrandArr(esArr []EdgeStrand, eID uint32) (ok bool, strand bool) {
	i := len(esArr) * int(eID) / int(esArr[len(esArr)-1].EdgeID)
	if i >= len(esArr) {
		i = len(esArr) - 1
	}
	count := 0
	if esArr[i].EdgeID == eID {
		ok = true
		strand = esArr[i].Strand
	} else if esArr[i].EdgeID < eID {
		j := i + 1
		for ; j < len(esArr); j += Step {
			count++
			if esArr[j].EdgeID == eID {
				ok = true
				strand = esArr[j].Strand
				return
			} else if esArr[j].EdgeID > eID {
				break
			}
		}
		i = j - Step
		if i < 0 {
			i = 0
		}
		j--
		if j >= len(esArr) {
			j = len(esArr) - 1
		}
		idx := (i + j) / 2
		if esArr[idx].EdgeID == eID {
			ok = true
			strand = esArr[idx].Strand
			return
		} else if esArr[idx].EdgeID < eID {
			i = idx + 1
		} else {
			j = idx - 1
		}
		for x := i; x <= j; x++ {
			if esArr[x].EdgeID == eID {
				ok = true
				strand = esArr[x].Strand
				break
			}
		}
	} else {
		j := i - 1
		for ; j >= 0; j -= Step {
			count--
			if esArr[j].EdgeID == eID {
				ok = true
				strand = esArr[j].Strand
				return
			} else if esArr[j].EdgeID < eID {
				break
			}
		}
		i = j + Step
		if i >= len(esArr) {
			i = len(esArr) - 1
		}
		j++
		if j < 0 {
			j = 0
		}
		idx := (i + j) / 2
		if esArr[idx].EdgeID == eID {
			ok = true
			strand = esArr[idx].Strand
			return
		} else if esArr[idx].EdgeID < eID {
			j = idx + 1
		} else {
			i = idx - 1
		}
		for x := j; x <= i; x++ {
			if esArr[x].EdgeID == eID {
				ok = true
				strand = esArr[x].Strand
				break
			}
		}
	}

	//fmt.Fprintf(os.Stderr, "[IsInEdgeStrandArr]count:%d\n", count)
	return
}

func AddHighFreqKmer2(ka []MapingKmerInfo, i int, dbgKmerInfoPac DBGKmerInfoPac, SeedLen uint16, LongEdgeMinLen int, SeedNumCountArr []SeedNumInfo, WinSize, kmerlen, MQ int) (int, int, []EdgeStrand, float32) {
	//MinChainScoreIdentityPercent *= 2
	var freqRepeatNum int
	highFreqMap := make(map[uint64][]MapPos)
	var lastMKI MapingKmerInfo
	id := ka[i].GetMapLen()
	i0 := i
	for ; i < len(ka); i++ {
		mki := ka[i]
		if mki.GetMapLen() != id {
			break
		}
		if mki.EID == math.MaxUint32 {
			kmer := mki.GetKmer()
			v, ok := highFreqMap[kmer]
			var mp MapPos
			mp.SetPosition(uint32(mki.GetPosition()))
			mp.SetStrand(mki.GetStrand())
			if ok {
				v = append(v, mp)
				highFreqMap[kmer] = v
				continue
			} else {
				highFreqMap[kmer] = []MapPos{mp}
			}

			/*eInfoArr, ok := dbgKmerInfoPac.LowFreqPercentEdgeMap[mki.GetKmer()]
			if !ok {
				log.Fatalf("[AddHighFreqKmer] mki:%v not include LowFreqPercentEdgeMap\n", PrintMKI(mki))
			}
			for _, eID := range eInfoArr {
				if SeedNumCountArr[eID].SeedNum < math.MaxUint8 {
					SeedNumCountArr[eID].SeedNum++
				}
				//v := SetKmerNum(SeedNumCountArr[eID], PLUS)
			}*/
		} else {
			if mki.EID == lastMKI.EID {
				if mki.GetEPosition() == lastMKI.GetEPosition() || mki.GetPosition() == lastMKI.GetPosition() {
					freqRepeatNum++
					continue
				}
			}
			strand := mki.GetStrand() == mki.GetEStrand()
			if strand {
				if SeedNumCountArr[mki.EID].SeedNumP < math.MaxUint8 {
					SeedNumCountArr[mki.EID].SeedNumP++
				}
			} else {
				if SeedNumCountArr[mki.EID].SeedNumM < math.MaxUint8 {
					SeedNumCountArr[mki.EID].SeedNumM++
				}
			}
			lastMKI = mki
		}
	}
	//beforeDN := len(edgeMiniKmerNumMap)
	//edgeMiniKmerNumMap = DeleteLowNumMinimersEdge(edgesArr, edgeMiniKmerNumMap, LongEdgeMinLen, MinChainScoreIdentityPercent)
	//fmt.Printf("[AddHighFreqKmer]highFreqNum:%d beforeDN:%d len(edgeMiniKmerNumMap):%d\n", highFreqNum, beforeDN, len(edgeMiniKmerNumMap))
	esArr, seedNumPercent := DeleteLowSeedNumMinimersEdge(SeedNumCountArr, uint32(LongEdgeMinLen), WinSize, kmerlen, MQ)
	if DebugModel {
		fmt.Printf("[AddHighFreqKmer2]len(ka):%d len(highFreqMap):%d freqRepeatNum:%d  after delete:%d seedNumPercent:%.2f\n", i-i0, len(highFreqMap), freqRepeatNum, len(esArr), seedNumPercent)
	}
	//fmt.Printf("[AddHighFreqKmer]esArr:%v\n", esArr)
	//aka = aka[:0]
	// EdgeSeedNumMap
	if len(esArr) == 0 {
		return i, i, esArr, seedNumPercent
	}

	idxKa := i0
	for j := i0; j < i; j++ {
		mki := ka[j]
		if mki.EID != math.MaxUint32 {
			ok, strand := IsInEdgeStrandArr(esArr, mki.EID)
			if ok && strand == mki.GetMapStrand() {
				ka[idxKa] = mki
				idxKa++
			}
		}
	}
	// add high freq mki
	addHighFreqNum := 0
	for kmer, RinfoArr := range highFreqMap {
		arr, ok := dbgKmerInfoPac.HighFreqMap[kmer]
		if !ok {
			log.Fatalf("[AddHighFreqKmer2] kmer:%v not in the HighFreqMap\n", kmer)
		}
		j := 0
		for _, es := range esArr {
			for ; j < len(arr); j++ {
				if arr[j].ID >= es.EdgeID {
					break
				}
			}
			for ; j < len(arr) && arr[j].ID == es.EdgeID; j++ {
				for _, RInfo := range RinfoArr {
					strand := RInfo.GetStrand() == arr[j].GetStrand()
					if strand != es.Strand {
						continue
					}
					var nmki MapingKmerInfo
					//nmki.SetKmer(kmer)
					nmki.SetPosition(uint64(RInfo.GetPosition()))
					nmki.SetStrand(RInfo.GetStrand())
					nmki.SetEPosition(arr[j].GetPosition())
					nmki.SetEStrand(arr[j].GetStrand())
					nmki.SetMapLen(id)
					//fmt.Printf("[AddHighFreqKmer]EPos:%d\n", arr[j].GetPosition())
					//nmki.SeedInfo = uint64(RInfo.Info)
					//nmki.Info = arr[j].Info;
					if idxKa >= i {
						log.Fatalf("[AddHighFreqKmer2]idxKa:%d >= i:%d\n", idxKa, i)
					}
					nmki.EID = es.EdgeID
					ka[idxKa] = nmki
					idxKa++
					addHighFreqNum++
				}
			}
		}
	}
	if DebugModel {
		fmt.Printf("[AddHighFreqKmer2]addHighFreqNum:%d\n", addHighFreqNum)
	}
	return i, idxKa, esArr, seedNumPercent
}

type MapingKmerInfoEPosArr []MapingKmerInfo

func (arr MapingKmerInfoEPosArr) Len() int {
	return len(arr)
}

func (arr MapingKmerInfoEPosArr) Less(i, j int) bool {
	if arr[i].EID == arr[j].EID {
		if arr[i].GetEPosition() == arr[j].GetEPosition() {
			return arr[i].GetMapLen() > arr[j].GetMapLen()
		} else {
			return arr[i].GetEPosition() < arr[j].GetEPosition()
		}
	} else {
		return arr[i].EID < arr[j].EID
	}
}

func (arr MapingKmerInfoEPosArr) Swap(i, j int) {
	arr[i], arr[j] = arr[j], arr[i]
}

func GetChainBlocks(ka []MapingKmerInfo) []MapingKmerInfo {
	for i, mk := range ka {
		if mk.GetMapLen() == 0 {
			continue
		}
		//fmt.Printf("[GetChainBlocks] %v\n", PrintMKI(mk))
		last := mk
		mE := last.GetEPosition()
		j := i + 1
		for ; j < len(ka); j++ {
			nk := ka[j]
			if nk.GetMapLen() == 0 {
				continue
			}
			nR := uint32(nk.GetPosition())
			nE := nk.GetEPosition()
			strand := nk.GetMapStrand()
			if nk.EID != mk.EID || nE > last.GetEPosition()+last.GetMapLen() {
				break
			}
			if strand != last.GetMapStrand() {
				continue
			}
			//fmt.Printf("[GetChainBlocks] %v\n", PrintMKI(nk))
			if strand == PLUS {
				if nR-uint32(last.GetPosition()) == nE-last.GetEPosition() {
					mk.SetMapLen(nE + nk.GetMapLen() - mE)
					ka[j].SetMapLen(0)
					last = nk
				}
			} else {
				if nE-last.GetEPosition() == (uint32(last.GetPosition())+last.GetMapLen())-(nR+nk.GetMapLen()) {
					//fmt.Printf("[GetChainBlocks] EID: %v, mk RPos: %v, EPos: %v,Len: %v, nk RPos: %v, EPos: %v, Len:%v\n", mk.EID, mk.GetPosition(), mk.GetEPosition(), mk.Len, nk.GetPosition(), nk.GetEPosition(), nk.GetMapLen())
					mk.SetPosition(uint64(nR))
					mk.SetMapLen(nE + nk.GetMapLen() - mE)
					//mk.SetEPosition(nE)
					//mk.GetMapLen() += uint16(((nk.GetPosition() + uint32(nk.GetMapLen())) - (mk.GetPosition() + uint32(mk.GetMapLen()))))
					ka[j].SetMapLen(0)
					last = nk
					//fmt.Printf("[GetChainBlocks] mk RPos: %v, EPos: %v,Len: %v\n", mk.GetPosition(), mk.GetEPosition(), mk.GetMapLen())
				}
			}
		}
		ka[i] = mk
	}

	// clean ka
	idx := 0
	for _, mk := range ka {
		if mk.GetMapLen() > 0 {
			ka[idx] = mk
			idx++
		}
	}
	fmt.Printf("[GetChainBlocks]len(ka):%d after ChainBlocks:%d\n", len(ka), idx)
	return ka[:idx]
}

func PrintMKIArr(mkiArr []MapingKmerInfo) {
	for i, mki := range mkiArr {
		fmt.Printf("[PrintMKIArr]arr[%d] EID:%v EPos:%d strand:%t RPos:%d Len:%d\n", i, mki.EID, mki.GetEPosition(), mki.GetStrand() == mki.GetEStrand(), mki.GetPosition(), mki.GetMapLen())
	}
}

func PrintKIArr(arr []KI) {
	for i, ki := range arr {
		fmt.Printf("[PrintKIArr]arr[%d] ID:%v Pos:%d strand:%t\n", i, ki.ID, ki.GetPosition(), ki.GetStrand())
	}
}

const ArrSize = 16

type MaxPathInfo struct {
	EID    uint32
	Base   uint32
	Score  float32
	Arr    [ArrSize]uint8 // Arr just store base index of Kb, if need mkiArr[i] = Kb[Base+Arr[i]]
	LenArr uint8
	Kb     []MapingKmerInfo // source mkiArr
	//Flag  uint8
	//Score2             uint16
	//FirstHalfArrScore  uint16 // if LenArr == ArrSize
	//SecondHalfArrScore uint16 // if LenArr == ArrSize
	//Path               []uint32
	//DistE int
}

func CountBit(covArr []uint8) (sum int) {
	for _, c := range covArr {
		for j := 0; j < 8; j++ {
			sum += int(c & 0x1)
			c >>= 1
		}
	}
	return
}

func GetSumMappingKIArrLen(ka []MapingKmerInfo, el int) (int, int, int, int) {
	eID := ka[0].EID
	num := (el + 8 - 1) / 8
	covArr := make([]uint8, num)
	a, b := uint32(math.MaxInt32), uint32(0)
	lp := ka[0].GetEPosition()
	llen := ka[0].GetMapLen()
	j := 0
	for ; j < len(ka); j++ {
		mki := ka[j]
		if mki.EID != eID {
			break
		}
		l := mki.GetMapLen()
		start, end := mki.GetEPosition(), mki.GetEPosition()+l
		if j > 0 && start == lp && llen == l {
			continue
		}
		for j := start; j < end; j++ {
			covArr[j>>3] |= (0x1 << (j & 0x7))
		}
		if start < a {
			a = start
		}
		if b < end {
			b = end
		}
		lp = start
		llen = l
	}
	sum := CountBit(covArr)
	//sum *= BaseNumInBit
	return sum, int(a), int(b), j
}

func CleanUint8Arr(arr []uint8) {
	for i, c := range arr {
		if c != math.MaxUint8 {
			arr[i] = 0
		}
	}
	return
}

func GetMKIDistance2(cen, mk MapingKmerInfo) (dX, dY int) {
	cenS := cen.GetStrand() == cen.GetEStrand()
	mkS := mk.GetStrand() == mk.GetEStrand()
	if cenS != mkS {
		dX, dY = -1, -1
		return
	}

	cR, cE, cLen := int(cen.GetPosition()), int(cen.GetEPosition()), int(cen.GetMapLen())
	mR, mE, mLen := int(mk.GetPosition()), int(mk.GetEPosition()), int(mk.GetMapLen())
	if cenS == PLUS {
		if cE < mE {
			dX = mE - (cE + cLen)
			dY = mR - (cR + cLen)
		} else {
			dX = cE - (mE + mLen)
			dY = cR - (mR + mLen)
		}
	} else {
		if cE < mE {
			dX = mE - (cE + cLen)
			dY = cR - (mR + mLen)
		} else {
			dX = cE - (mE + mLen)
			dY = mR - (cR + cLen)
		}
	}
	return
}

func GetTwoBlockScore2(cen, mk MapingKmerInfo, GapCost int) (sc float32) {
	var gapLen int
	dX, dY := GetMKIDistance2(cen, mk)

	//fmt.Printf("[GetTwoBlockScore2]dX:%d dY:%d\n", dX, dY)
	/*if distanceE != 0 {
		fmt.Printf("[GetTwoBlockScore]cen EID:%v, mk EID: %v, dX: %v,dY: %v\n", cen.EID, mk.EID, dX, dY)
	}*/
	if dX < 0 || dY < 0 {
		//log.Fatalf("[GetTwoBlockScore] dX: %v dY: %v, cen: %v, mk: %v\n", dX, dY, PrintMKI(cen), PrintMKI(mk))
		//fmt.Fprintf(os.Stderr, "[GetTwoBlockScore] dX: %v dY: %v, cen: %v, mk: %v\n", dX, dY, PrintMKI(cen), PrintMKI(mk))
		sc = 0
		return
	}
	/*if dX > 2000 || dY > 2000 {
		sc = 0
		return
	}*/
	gapLen = AbsInt(dX - dY)
	min := Min(dX, dY)
	if gapLen > min/8+10 {
		sc = 0
		return
	}
	a := float32(mk.GetMapLen())
	// a gap length of gapLen costs
	var b float32
	if gapLen > 1 {
		//b = float64(0.01)*float64(mk.GetMapLen())*float64(gapLen) + float64(0.5)*math.Log2(float64(gapLen))
		//b = float32(0.01)*float32(mk.GetMapLen())*float32(gapLen) + float32(0.5)*float32(gapLen/10)
		//b = float32(0.01) * float32(mk.GetMapLen()) * float32(gapLen) + 0.1 * min
		b = 0.1*float32(min) + float32(gapLen)
	} else {
		if dX == 0 && dY == 0 {
			b = 0
		} else {
			b = 1
		}
	}
	sc = a - b
	//sc = int(mk.GetMapLen()) - gapLen*GapCost
	//fmt.Printf("[GetTwoBlockScore]dX: %v, dY: %v, sc: %v, a: %v, b: %v\n", dX, dY, sc, a, b)
	//fmt.Printf("[GetTwoBlockScore2]sc:%d\n", sc)
	return
}

func BlockToNearBlocksScoreArr(kb []MapingKmerInfo, Idx int, GapCost int, bsA []uint8, SeedLen uint16) {
	cen := kb[Idx]
	strand := cen.GetEStrand() == cen.GetStrand()
	//maxExtendNum := 6
	//count := 0
	extendSize := 400
	cenR, cenE, cenLen := int(cen.GetPosition()), int(cen.GetEPosition()), int(cen.GetMapLen())
	//lastEPos := cenE
	//fmt.Printf("[BlockToNearBlocksScoreArr]Idx:%d cen.RPos:%d cen.EPos:%d\n", Idx, cen.GetPosition(), cen.GetEPosition())
	// BACKWARD
	for i := Idx - 1; i >= 0; i-- {
		bs := bsA[i]
		if bs == math.MaxUint8-1 {
			break
		} else if bs == math.MaxUint8 {
			continue
		}
		mk := kb[i]
		mkR, mkE, mkLen := int(mk.GetPosition()), int(mk.GetEPosition()), int(mk.GetMapLen())
		if mkE+mkLen < cenE-extendSize {
			break
		}

		if strand == PLUS {
			if mkR+mkLen > cenR || mkE+mkLen > cenE {
				dX, dY := cenR-(mkR+mkLen), cenE-(mkE+mkLen)
				if dX > dY {
					dX, dY = dY, dX
				}
				if (dX < 0 || dY < 0) && dY-dX < 8 && mkLen+dX > 10 {
					mk.SetMapLen(uint32(mkLen + dX))
					kb[i] = mk
				} else {
					bsA[i] = math.MaxUint8 - 1
					continue
				}
			}
			sc := GetTwoBlockScore2(cen, mk, GapCost)
			//fmt.Printf("[BlockToNearBlocksScoreArr]sc:%v\n", sc)
			if int(sc) > int(bs) {
				if sc >= math.MaxUint8-1 {
					sc = math.MaxUint8 - 2
				}
				bsA[i] = uint8(sc)
			}
		} else { // strand == MINUS
			if mkR < cenR+cenLen || mkE+mkLen > cenE {
				dX, dY := mkR-(cenR+cenLen), cenE-(mkE+mkLen)
				if dX > dY {
					dX, dY = dY, dX
				}
				if (dX < 0 || dY < 0) && dY-dX < 8 && mkLen+dX > 10 {
					mk.SetMapLen(uint32(mkLen + dX))
					mk.SetPosition(uint64(mkR - dX))
					kb[i] = mk
				} else {
					bsA[i] = math.MaxUint8 - 1
					continue
				}
			}
			sc := GetTwoBlockScore2(cen, mk, GapCost)
			if int(sc) > int(bs) {
				if sc >= math.MaxUint8-1 {
					sc = math.MaxUint8 - 2
				}
				bsA[i] = uint8(sc)
			}
		}
	}

	// FORWARD
	for i := Idx + 1; i < len(kb); i++ {
		bs := bsA[i]
		if bs == math.MaxUint8-1 {
			break
		} else if bs == math.MaxUint8 {
			continue
		}
		mk := kb[i]
		mkR, mkE, mkLen := int(mk.GetPosition()), int(mk.GetEPosition()), int(mk.GetMapLen())
		if mkE > cenE+cenLen+extendSize {
			break
		}
		var sc float32
		if strand == PLUS {
			if cenR+cenLen > mkR || cenE+cenLen > mkE {
				dX, dY := mkR-(cenR+cenLen), mkE-(cenE+cenLen)
				if dX > dY {
					dX, dY = dY, dX
				}
				if (dX < 0 || dY < 0) && dY-dX < 8 && mkLen+dX > 10 {
					mk.SetMapLen(uint32(mkLen + dX))
					mk.SetPosition(uint64(mkR - dX))
					mk.SetEPosition(uint32(mkE - dX))
					kb[i] = mk
				} else {
					bsA[i] = math.MaxUint8 - 1
					continue
				}
			}
			sc = GetTwoBlockScore2(cen, mk, GapCost)
			if int(sc) > int(bs) {
				if sc >= math.MaxUint8-1 {
					sc = math.MaxUint8 - 2
				}
				bsA[i] = uint8(sc)
			}
		} else { //strand == MINUS
			if cenE+cenLen > mkE || mkR+mkLen > cenR {
				dX, dY := mkE-(cenE+cenLen), cenR-(mkR+mkLen)
				if dX > dY {
					dX, dY = dY, dX
				}
				if (dX < 0 || dY < 0) && dY-dX < 8 && mkLen+dX > 10 {
					mk.SetMapLen(uint32(mkLen + dX))
					mk.SetEPosition(uint32(mkE - dX))
					kb[i] = mk
				} else {
					bsA[i] = math.MaxUint8 - 1
					continue
				}
			}

			sc = GetTwoBlockScore2(cen, mk, GapCost)
			if int(sc) > int(bs) {
				if sc >= math.MaxUint8-1 {
					sc = math.MaxUint8 - 2
				}
				bsA[i] = uint8(sc)
			}
		}
	}

	//sort.Sort(ScoreArr(bsA))
	//return idxScArr
}

type IdxScore struct {
	Idx uint16
	Sc  uint16
}

func GetMaxScoreFromBaSA(blockToAnchorsScoreArr []uint8) IdxScore {
	var idxSc IdxScore
	for i, sc := range blockToAnchorsScoreArr {
		if sc >= math.MaxUint8-1 {
			continue
		}
		if uint16(sc) > idxSc.Sc {
			idxSc.Sc = uint16(sc)
			idxSc.Idx = uint16(i)
		}
	}
	return idxSc
}

type Uint16Slice []uint16

func (p Uint16Slice) Len() int           { return len(p) }
func (p Uint16Slice) Less(i, j int) bool { return p[i] < p[j] }
func (p Uint16Slice) Swap(i, j int)      { p[i], p[j] = p[j], p[i] }

func GetMappingKIArrScore(mkArr []MapingKmerInfo, indexArr []uint16, GapCost int) (score float32) {
	if len(indexArr) < 1 {
		return
	}
	last := mkArr[indexArr[0]]
	score = float32(last.GetMapLen())
	for i := 1; i < len(indexArr); i++ {
		nk := mkArr[indexArr[i]]
		sc := GetTwoBlockScore2(last, nk, GapCost)
		//fmt.Printf("[GetMappingKIArrScore] nk EdgePos: %v, ReadPos: %v, Len: %v, score: %v\n", nk.GetEPosition(), nk.GetPosition(), nk.GetMapLen(), sc)
		score += sc
		last = nk
	}
	return
}

type MaxPathInfoArr []MaxPathInfo

func (arr MaxPathInfoArr) Len() int {
	return len(arr)
}

func (arr MaxPathInfoArr) Less(i, j int) bool {
	return arr[i].Score < arr[j].Score
}
func (arr MaxPathInfoArr) Swap(i, j int) {
	arr[i], arr[j] = arr[j], arr[i]
}

func GetMaxScoreArr(ka []MapingKmerInfo, SeedLen, GapCost, MaxGapLen int, SeqType uint8, el int, bsA []uint8, mkiAllFlag bool) (mp MaxPathInfo) {
	//var blockToBlocksScorePat [][]Score
	//sort.Sort(MapingKmerInfoArr(ki.MkInfoArr))
	/*for i, mki := range ki.MkInfoArr {
		fmt.Printf("[GetMaxScoreArr] ki.MkInfoArr[%v] EdgePos: %v, ReadPos : %v, Len: %v\n", i, mki.GetEPosition(), mki.GetPosition(), mki.GetMapLen())
	}*/
	//anchorMKI, idx := GetAnchorMKI(ki)
	//PrintMKIArr(ka)
	if len(ka) == 1 {
		//mp.Arr = make([]uint16, 1)
		mp.Arr[0] = 0
		mp.LenArr = 1
		mp.Score = float32(ka[0].GetMapLen())
		return
	}
	kaLen := len(ka)
	if kaLen > math.MaxUint16 {
		fmt.Fprintf(os.Stderr, "[GetMaxScoreArr] kaLen:%d > %d, el:%d, return\n", kaLen, math.MaxUint16, el)
		return
	}

	maxLoop := kaLen/10 + 1
	if maxLoop > 10 {
		maxLoop = 10
	}
	usedArr := make([]int, 0, maxLoop)
	bCLen := el / 40
	if bCLen > kaLen {
		bCLen = kaLen
	}
	var maxScore float32
	bCArr := make([]uint16, 0, bCLen)
	maxArr := make([]uint16, bCLen)
	maxIndexArr := maxArr[:0]
	//fmt.Printf("[GetMaxScoreArr]kaLen:%d maxLoop:%d\n", kaLen, maxLoop)
	for i := 0; i < maxLoop; i++ {
		// found anchor
		maxLen := uint32(0)
		idx := -1
		//pos := time.Now().Nanosecond() % kaLen
		for j, cb := range ka {
			if maxLen < uint32(cb.GetMapLen()) && !IsInIntArr(usedArr, j) {
				maxLen = uint32(cb.GetMapLen())
				idx = j
			}
		}
		if idx < 0 {
			break
		}
		usedArr = append(usedArr, idx)
		CleanUint8Arr(bsA)
		//fmt.Printf("[GetMaxScoreArr] anchor[%v] edgeID: %v, mki ReadPos : %v,EdgePos: %v, Len: %v\n", idx, anchorMKI.EID, anchorMKI.GetPosition(), anchorMKI.GetEPosition(), anchorMKI.GetMapLen())
		//anchorMKI := ka[idx]
		bCArr = bCArr[:0]
		bCArr = append(bCArr, uint16(idx))
		bsA[idx] = math.MaxUint8 - 1
		BlockToNearBlocksScoreArr(ka, idx, GapCost, bsA, uint16(SeedLen))
		//blockToBlocksScorePat = append(blockToBlocksScorePat, baSA)
		//var maxSc IdxScore
		for {
			//fmt.Printf("[GetMaxScoreArr]bsA: %v\n", bsA)
			maxSc := GetMaxScoreFromBaSA(bsA)
			if maxSc.Sc <= 0 {
				break
			}
			/*tk, sc := CheckChainBoundary(maxA, ki.MkInfoArr[maxSc.Idx], maxSc)
			if sc.Sc < maxSc.Sc {
				continue
			}*/
			// add to Max Chain
			//visitArr[maxSc.Idx] = true
			bsA[maxSc.Idx] = math.MaxUint8 - 1
			bCArr = append(bCArr, uint16(maxSc.Idx))
			//fmt.Printf("[GetMaxScoreArr]added mki idx:%d bsA:%v\n", maxSc.Idx, bsA)
			BlockToNearBlocksScoreArr(ka, int(maxSc.Idx), GapCost, bsA, uint16(SeedLen))
			//blockToBlocksScorePat = append(blockToBlocksScorePat, baSA)
		}
		//sort.Sort(MapingKmerInfoArr(bChainArr))
		//fmt.Printf("[GetMaxScoreArr]bsA:   %v\n", bsA)
		//fmt.Printf("[GetMaxScoreArr]bCarr: %v\n", bCArr)
		//fmt.Printf("[GetMaxScoreArr]bCArr:%v ", bCArr)
		sort.Sort(Uint16Slice(bCArr))
		arrScore := GetMappingKIArrScore(ka, bCArr, GapCost)
		//fmt.Printf(" arrScore:%.1f\n", arrScore)
		if arrScore-maxScore > 0.5 {
			if len(bCArr) > len(maxArr) {
				maxArr = make([]uint16, len(bCArr))
			}
			copy(maxArr, bCArr)
			maxIndexArr = maxArr[:len(bCArr)]
			maxScore = arrScore
			//i = -1
		}
	}
	//mp.Arr = GetMKIArrByIndex(ka, maxIndexArr)
	if len(maxIndexArr) < ArrSize {
		for x, idx := range maxIndexArr {
			mp.Arr[x] = uint8(idx)
		}
		mp.LenArr = uint8(len(maxIndexArr))
	} else {
		mp.Arr[0] = uint8(maxIndexArr[0])
		mp.Arr[ArrSize-1] = uint8(maxIndexArr[len(maxIndexArr)-1])
		mp.LenArr = ArrSize
	}
	mp.Score = maxScore
	if maxScore == 0 {
		fmt.Printf("[GetMaxScoreArr]maxScore:%.1f\n", maxScore)
		PrintMKIArr(ka)
	}
	return mp
}

type MaxPathInfoTmpArr []MaxPathInfo

func (arr MaxPathInfoTmpArr) Len() int {
	return len(arr)
}

func (arr MaxPathInfoTmpArr) Less(i, j int) bool {
	ai, bi := arr[i].Kb[arr[i].Base+uint32(arr[i].Arr[0])].GetPosition(), arr[i].Kb[arr[i].Base+uint32(arr[i].Arr[arr[i].LenArr-1])].GetPosition()
	aj, bj := arr[j].Kb[arr[j].Base+uint32(arr[j].Arr[0])].GetPosition(), arr[j].Kb[arr[j].Base+uint32(arr[j].Arr[arr[j].LenArr-1])].GetPosition()
	if ai > bi {
		ai, bi = bi, ai
	}
	if aj > bj {
		aj, bj = bj, aj
	}
	if ai < aj {
		return true
	} else if ai > aj {
		return false
	} else {
		return arr[i].Score > arr[j].Score
	}
}

func (arr MaxPathInfoTmpArr) Swap(i, j int) {
	arr[i], arr[j] = arr[j], arr[i]
}

func GetMaxChainBlockArr(kb []MapingKmerInfo, maxPathInfoArr []MaxPathInfo, edgesArr []DBGEdge, nodesArr []DBGNode, SeedLen, GapCost, MaxGapLen, kmerlen int, SeqType uint8, seedNumPercent float32, rl int, logfp io.Writer) []MaxPathInfo {
	//edgesKIArr := GetEdgesKIArr(kb, edgesArr, MinKmerScore, LongEdgeMinLen)
	//sort.Sort(EdgeKIArr(edgesKIArr))
	//var maxScore int
	//var maxScArr []int
	//var loopNum int
	//var cBMRIArr []ChainBlockMapLongReadInfo
	//var eIDArr []uint32
	//MinChainScoreIdentityPercent = seedNumPercent * 3 / 2
	if seedNumPercent > 0.9 {
		seedNumPercent = 0.9
	}
	MinChainScoreIdentityPercent := seedNumPercent * 2 / 3
	maxPathInfoArr = maxPathInfoArr[:0]
	//EdgeMinLen := (kmerlen-1)*2 - 10
	bsShare := make([]uint8, 30) // visit just used lowest-bit of bsA element &0x1
	//PrintMKIArr(kb)
	//var tmpArr []MaxPathInfo
	// search long read local best mapping to the edge
	totalNum, GetMaxScoreNum := 0, 0
	for i := 0; i < len(kb); {
		eID := kb[i].EID
		strand := kb[i].GetStrand() == kb[i].GetEStrand()
		totalNum++
		//fmt.Printf("[GetMaxChainBlockArr]i:%d j:%d\n", i, j)
		//num := MaxInt(num1, num2)
		//e := edgesArr[ki.MkInfoArr[0].EID]
		/*if e.StartNID == e.EndNID || e.GetTwoEdgesCycleFlag() > 0 {
			continue
		}
		if len(e.Ks) < kmerlen*2-3 {
			continue
		}*/
		//PrintAddedMpArr(ki.MkInfoArr)
		el := edgesArr[eID].GetSeqLen()
		sc, a, b, j := GetSumMappingKIArrLen(kb[i:], el)
		j += i
		if DebugModel {
			fmt.Fprintf(logfp, "eID:%d el:%d edgeMinNum:%d len(ka):%d mapping region:%d sc:%d\n", eID, el, edgesArr[eID].CovD, j-i, b-a, sc)
		}
		/*if sc <= SeedLen {
			fmt.Fprintf(os.Stderr, "[GetMaxChainBlockArr]eID:%d el:%d edgeMinNum:%d len(ka):%d mapping region:%d sc:%d\n", eID, el, edgesArr[eID].CovD, j-i, b-a, sc)
			for _, mki := range kb[i:j] {
				fmt.Fprintf(os.Stderr, "%s\n", PrintMKI(mki))
			}
		}*/
		if float32(sc) < float32(el)*MinChainScoreIdentityPercent {
			i = j
			continue
		}
		PrintMKIArr(kb[i:j])
		if j-i > cap(bsShare) {
			bsShare = make([]uint8, (j - i))
		}
		bsA := bsShare[:j-i]
		//bsA = InitStrandInfo(kb[i:j], bsA, strand)
		//PrintMKIArr(kb[i:j])
		//PrintMKI(mki)
		maxPI := GetMaxScoreArr(kb[i:j], SeedLen, GapCost, MaxGapLen, SeqType, el, bsA, false)
		fmt.Printf("[GetMaxChainBlockArr]maxPI.Arr:%v score:%.1f\n", maxPI.Arr[:maxPI.LenArr], maxPI.Score)
		if maxPI.LenArr < 1 || maxPI.Score <= float32(SeedLen) {
			i = j
			continue
		}
		GetMaxScoreNum++

		mka, mkb := kb[i+int(maxPI.Arr[0])], kb[i+int(maxPI.Arr[maxPI.LenArr-1])]
		ml := int(mkb.GetMapLen())
		startE, endE := int(mka.GetEPosition()), int(mkb.GetEPosition())+ml
		var startR, endR int
		if strand == MINUS {
			startR, endR = int(mkb.GetPosition()), int(mka.GetPosition())+int(mka.GetMapLen())
		} else {
			startR, endR = int(mka.GetPosition()), int(mkb.GetPosition())+ml
		}
		//fmt.Printf("[GetMaxChainBlockArr]maxPI.Score:%f start:%d end:%d startE:%d endE:%d el:%d\n", maxPI.Score, startR, endR, startE, endE, el)
		readRegLen := endR - startR
		edgeRegLen := endE - startE
		min := MinInt(readRegLen, edgeRegLen)
		added := false
		if float32(maxPI.Score) > float32(el)*MinChainScoreIdentityPercent && AbsInt(readRegLen-edgeRegLen) < min/8 {
			added = true
		}
		if DebugModel {
			fmt.Fprintf(logfp, "added:%v mpi LenArr:%d score:%.1f\n", added, maxPI.LenArr, maxPI.Score)
		}
		if added {
			maxPI.Base = uint32(i)
			maxPI.EID = eID
			maxPI.Kb = kb
			maxPathInfoArr = append(maxPathInfoArr, maxPI)
		}
		i = j
		//var cb ChainBlockMapLongReadInfo
		//lastMK := maxPI.Arr[len(maxPI.Arr)-1]
		//cb.Start, cb.End, cb.Score = int(maxPI.Arr[0].GetPosition()), int(lastMK.GetPosition())+int(lastMK.GetMapLen()), maxPI.Score
		//tmpArr = append(tmpArr, maxPI)
		//if AnchorFilter(cBMRIArr, cb) == false {
		//eIDArr = append(eIDArr, e.ID)
		//}
	}
	fmt.Fprintf(logfp, "[GetMaxChainBlockArr]totalNum:%d GetMaxScoreNum:%d len(maxPathInfoArr):%d\n", totalNum, GetMaxScoreNum, len(maxPathInfoArr))
	sort.Sort(MaxPathInfoTmpArr(maxPathInfoArr))
	//fmt.Printf("[GetMaxChainBlockArr] ReadLen: %v\n", ReadLen)
	/*for i, mpi := range maxPathInfoArr {
		start, end := mpi.Arr[0].GetPosition(), mpi.Arr[len(mpi.Arr)-1].GetPosition()+uint32(mpi.Arr[len(mpi.Arr)-1].GetMapLen())
		fmt.Printf("[GetMaxChainBlockArr] maxPathInfoArr[%v]: Edge ID: %v, EdgeLen: %v, Read region[%v---%v] Len:%v, Edge region[%v---%v], score: %v, Percent: %v\n", i, mpi.Arr[0].EID, len(edgesArr[mpi.Arr[0].EID].Ks), start, end, end-start, mpi.Arr[0].GetEPosition(), mpi.Arr[len(mpi.Arr)-1].GetEPosition(), mpi.Score, mpi.Score*100/(int(end)-int(start)))
	}*/
	/*for y, mki := range ele.Arr {
		var cb ChainBlockMapLongReadInfo
		lastMK := mpi.Arr[len(mpi.Arr)-1]
		cb.Start, cb.End, cb.Score = int(mpi.Arr[0].GetPosition()), int(lastMK.GetPosition())+int(lastMK.GetMapLen()), mpi.Score
		fmt.Printf("[GetMaxChainBlockArr][%v] EID: %v, start: %v, end: %v, score: %v\n", i, lastMK.EID, cb.Start, cb.End, cb.Score)
	}*/

	return maxPathInfoArr
}

func GetMPIInfo(mpi *MaxPathInfo) (start, end int, strand bool, startE, endE int) {
	//fmt.Printf("[GetMPIInfo]Base:%d LenArr:%d len(Arr):%d\n", mpi.Base, mpi.LenArr, len(mpi.Arr))
	if mpi.LenArr == 0 {
		return
	}
	mka, mkb := mpi.Kb[mpi.Base+uint32(mpi.Arr[0])], mpi.Kb[mpi.Base+uint32(mpi.Arr[mpi.LenArr-1])]
	startE, endE = int(mka.GetEPosition()), int(mkb.GetEPosition()+mkb.GetMapLen())
	strand = mka.GetMapStrand()
	if strand == PLUS {
		start, end = int(mka.GetPosition()), int(mkb.GetPosition())+int(mkb.GetMapLen())
	} else {
		start, end = int(mkb.GetPosition()), int(mka.GetPosition())+int(mka.GetMapLen())
	}
	return
}

func GetMPIRegionScore(mpi *MaxPathInfo, start, end uint32) (sc float32) {
	var lastMKI MapingKmerInfo
	mk0 := mpi.Kb[mpi.Base+uint32(mpi.Arr[0])]
	strand := mk0.GetStrand() == mk0.GetEStrand()
	if strand == PLUS {
		for _, idx := range mpi.Arr {
			mki := mpi.Kb[mpi.Base+uint32(idx)]
			s := uint32(mki.GetPosition())
			e := s + mki.GetMapLen()
			if e <= start || s >= end {
				continue
			}
			maxS := MaxUint32(start, s)
			minE := MinUint32(end, e)
			if minE <= maxS {
				log.Fatalf("[GetMPIRegionScore] minE:%d <= maxS:%d\n", minE, maxS)
			}
			mkLen := minE - maxS
			if lastMKI.EID > 0 {
				sc += GetTwoBlockScore2(lastMKI, mki, 1)
				if mki.GetMapLen() > mkLen {
					sc -= float32(mki.GetMapLen() - mkLen)
				}
			} else {
				sc += float32(mkLen)
			}
			lastMKI = mki
		}
	} else {
		for i := len(mpi.Arr) - 1; i >= 0; i-- {
			idx := mpi.Arr[i]
			mki := mpi.Kb[mpi.Base+uint32(idx)]
			s := uint32(mki.GetPosition())
			e := s + mki.GetMapLen()
			if e <= start || s >= end {
				continue
			}
			maxS := MaxUint32(start, s)
			minE := MinUint32(end, e)
			if minE <= maxS {
				log.Fatalf("[GetMPIRegionScore] minE:%d <= maxS:%d\n", minE, maxS)
			}
			mkLen := minE - maxS
			if lastMKI.EID > 0 {
				sc += GetTwoBlockScore2(lastMKI, mki, 1)
				if mki.GetMapLen() > mkLen {
					sc -= float32(mki.GetMapLen() - mkLen)
				}
			} else {
				sc += float32(mkLen)
			}
			lastMKI = mki
		}
	}

	return
}

func DiffFloat32(f1, f2 float32) (diff float32) {
	if f1 > f2 {
		diff = f1 - f2
	} else {
		diff = f2 - f1
	}
	return diff
}

func PrintMPI(mpi *MaxPathInfo, edgesArr []DBGEdge) string {
	var s string
	if mpi.LenArr == 0 {
		return s
	}
	start, end, strand, startE, endE := GetMPIInfo(mpi)
	if startE == endE {
		endE++
	}
	s = fmt.Sprintf("EdgeID:%d R[%d--%d] E[%d--%d] std:%t eLen:%d Len:%d score:%0.1f Percent:%0.3f", mpi.EID, start, end, startE, endE, strand, edgesArr[mpi.EID].GetSeqLen(), endE-startE, mpi.Score, mpi.Score/float32(edgesArr[mpi.EID].GetSeqLen()))
	return s
}

func FilterMPIArr(maxPathInfoArr []MaxPathInfo, kmerlen, ReadLen int, edgesArr []DBGEdge, logfp io.Writer) []MaxPathInfo {
	// filter low mapping quality edge and filter only maping middle region of edge
	boundLen := kmerlen
	diff := 10
	filterFlag := make([]bool, len(maxPathInfoArr))
	//fmt.Printf("[FilterMPIArr]filterFlag:%v\n", filterFlag)
	for i := range maxPathInfoArr {
		if filterFlag[i] {
			continue
		}
		mpi := &maxPathInfoArr[i]
		if mpi.Score < SeedLen*2 {
			filterFlag[i] = true
			continue
		}
		//fmt.Printf("[FilterMPIArr]mpiArr[%d]:%v\n", i, *mpi)
		start, end, strand, startE, endE := GetMPIInfo(mpi)
		//Start, End := int(mka.GetPosition()), int(mkb.GetPosition())+int(mkb.GetMapLen())
		e := &edgesArr[mpi.EID]
		el := e.GetSeqLen()
		if el > 6*kmerlen {
			boundLen = 2 * kmerlen
		}
		if end-start < kmerlen/2 {
			filterFlag[i] = true
			continue
		}
		//eS, eE := int(mka.GetEPosition()), int(mkb.GetEPosition())+int(mkb.GetMapLen())
		if strand {
			if start > boundLen && startE > boundLen {
				filterFlag[i] = true
				continue
			}
			if end < ReadLen-boundLen && endE < el-boundLen {
				filterFlag[i] = true
				continue
			}
		} else {
			if start > boundLen && endE < el-boundLen {
				filterFlag[i] = true
				continue
			}
			if end < ReadLen-boundLen && startE > boundLen {
				filterFlag[i] = true
				continue
			}
		}

		for j := i + 1; j < len(maxPathInfoArr); j++ {
			if filterFlag[j] {
				continue
			}
			hmpi := &maxPathInfoArr[j]
			if hmpi.Score < SeedLen*2 {
				filterFlag[i] = true
				continue
			}
			hStart, hEnd, _, _, _ := GetMPIInfo(hmpi)
			if end <= hStart {
				break
			}
			//fmt.Printf("[FilterMPIArr]i:%d j:%d start:%d end:%d hStart:%d hEnd:%d mpi.Score:%d hmpi.Score:%d\n", i, j, start, end, hStart, hEnd, mpi.Score, hmpi.Score)
			if hStart-diff <= start && end <= hEnd+diff && end-start <= hEnd-hStart && mpi.Score-hmpi.Score > 0.5 {
				filterFlag[j] = true
			} else if start-diff <= hStart && hEnd <= end+diff && hEnd-hStart <= end-start && hmpi.Score-mpi.Score > 0.5 {
				filterFlag[i] = true
				break
			} else if hStart-diff <= start && end <= hEnd+diff && end-start <= hEnd-hStart && hmpi.Score/float32(hEnd-hStart)-mpi.Score/float32(end-start) > 0.005 {
				filterFlag[i] = true
				break
			} else if start-diff <= hStart && hEnd <= end+diff && hEnd-hStart <= end-start && mpi.Score/float32(end-start)-hmpi.Score/float32(hEnd-hStart) > 0.005 {
				filterFlag[j] = true
			} /*else if hStart <= Start && End <= hEnd && ((hEnd-hStart)-(End-Start)) < (End-Start)/10 && int(mpi.Score)*100/(End-Start) > (int(hmpi.Score)*100/(hEnd-hStart)) {
				filterFlag[i] = true
				filterFlag[j] = true
				break
			} else if Start <= hStart && hEnd <= End && ((End-Start)-(hEnd-hStart)) < (hEnd-hStart)/10 && int(mpi.Score)*100/(End-Start) < (int(hmpi.Score)*100/(hEnd-hStart)) {
				filterFlag[i] = true
				filterFlag[j] = true
				break
			}*/
		}
	}

	idx := 0
	for i, mpi := range maxPathInfoArr {
		if !filterFlag[i] {
			maxPathInfoArr[idx] = mpi
			idx++
		}
	}

	fmt.Fprintf(logfp, "[FilterMPIArr]deleted filter maxPathInfoArr num:%d\n", len(maxPathInfoArr)-idx)
	maxPathInfoArr = maxPathInfoArr[:idx]

	// filter contained mpi
	filterFlag = filterFlag[:len(maxPathInfoArr)]
	for i := range filterFlag {
		filterFlag[i] = false
	}
	for i := range maxPathInfoArr {
		if filterFlag[i] {
			continue
		}
		mpi := &maxPathInfoArr[i]
		start, end, _, _, _ := GetMPIInfo(mpi)
		//start, end := uint32(mka.GetPosition()), uint32(mkb.GetPosition())+mkb.GetMapLen()
		//e := edgesArr[mpi.EID]
		//el := e.GetSeqLen()
		for j := i + 1; j < len(maxPathInfoArr); j++ {
			if filterFlag[j] {
				continue
			}
			nm := &maxPathInfoArr[j]
			nstart, nend, _, _, _ := GetMPIInfo(nm)
			//nstart, nend := uint32(a.GetPosition()), uint32(b.GetPosition())+b.GetMapLen()
			if nend > end || (nstart == start && nend == end) {
				break
			}
			//e2 := edgesArr[nm.EID]
			//el2 := e2.GetSeqLen()
			if start <= nstart && nend <= end && mpi.LenArr < ArrSize {
				//sz := len(mpi.Arr)
				//dr := FORWARD
				sc := GetMPIRegionScore(mpi, uint32(nstart), uint32(nend))
				//fmt.Printf("[FilterMPIArr]i:%d j:%d sc:%.1f nm.Score:%.1f\n", i, j, sc, nm.Score)
				/*for j, idx := range mpi.Arr[:mpi.LenArr] {
					fmt.Printf("mpi[%d]%s\n", j, PrintMKI(mpi.Kb[mpi.Base+uint32(idx)]))
				}
				for j, idx := range nm.Arr[:nm.LenArr] {
					fmt.Printf("nm[%d]%s\n", j, PrintMKI(nm.Kb[nm.Base+uint32(idx)]))
				}*/
				if sc-nm.Score > 0.5 {
					filterFlag[j] = true
				} else if DiffFloat32(sc, nm.Score) < 0.5 {
					//AlignmentBlocks()
				}
			}
		}
	}

	count := 0
	for i := range maxPathInfoArr {
		if filterFlag[i] {
			continue
		}
		maxPathInfoArr[count] = maxPathInfoArr[i]
		count++
	}
	if DebugModel {
		for i := range maxPathInfoArr[:count] {
			fmt.Fprintf(logfp, "mpiArr[%d]:%s\n", i, PrintMPI(&maxPathInfoArr[i], edgesArr))
		}
	}
	fmt.Fprintf(logfp, "[FilterMPIArr]len(maxPathInfoArr):%d after filter contained mpi count:%d\n", len(maxPathInfoArr), count)
	maxPathInfoArr = maxPathInfoArr[:count]
	return maxPathInfoArr
}

func IsBubble(e1, e2 *DBGEdge, nodesArr []DBGNode) bool {
	if (e1.StartNID == e1.EndNID) || (e2.StartNID == e2.EndNID) || (e1.StartNID < 2) || (e1.EndNID < 2) || (e2.StartNID < 2) || (e2.EndNID < 2) {
		return false
	}
	// test startNID
	if IsInComing(nodesArr[e1.StartNID].EdgeIDIncoming, e1.ID) {
		if !IsInComing(nodesArr[e1.StartNID].EdgeIDIncoming, e2.ID) {
			return false
		}
	} else {
		if !IsInComing(nodesArr[e1.StartNID].EdgeIDOutcoming, e2.ID) {
			return false
		}
	}
	// test endNID
	if IsInComing(nodesArr[e1.EndNID].EdgeIDIncoming, e1.ID) {
		if !IsInComing(nodesArr[e1.EndNID].EdgeIDIncoming, e2.ID) {
			return false
		}
	} else {
		if !IsInComing(nodesArr[e1.EndNID].EdgeIDOutcoming, e2.ID) {
			return false
		}
	}

	return true
}

func GetAnchorMPIIdx(maxPathInfoArr []MaxPathInfo, priorIdxArr []uint16, edgesArr []DBGEdge) int {
	anchorIdx := -1
	maxScore := float32(0)
	for _, idx := range priorIdxArr {
		mpi := maxPathInfoArr[idx]
		if mpi.Score <= maxScore {
			continue
		}
		e := edgesArr[mpi.EID]
		if e.GetUniqueFlag() > 0 || e.GetBubbleFlag() > 0 || e.GetBubbleRepeatFlag() > 0 {
			maxScore = mpi.Score
			anchorIdx = int(idx)
		}
	}
	return anchorIdx
}

func GetMPIArrPrior(highQMPIArr []MaxPathInfo, kmerlen int, edgesArr []DBGEdge, nodesArr []DBGNode) ([]MaxPathInfo, int, []uint16) {
	// prior choose unique edge or bubble and bubbleRepeat edge
	/*if priorIdxArrLen < 5 {
		priorIdxArrLen = 5
	}*/
	priorIdxArr := make([]uint16, len(highQMPIArr))
	//var score float32
	//var repeatOK bool
	//MinRegLen := uint32(800)
	for i, mpi := range highQMPIArr {
		if priorIdxArr[i] == math.MaxUint16 {
			continue
		}
		//mka := mpi.Kb[mpi.Base+uint32(mpi.Arr[0])]
		mka, mkb := mpi.Kb[mpi.Base+uint32(mpi.Arr[0])], mpi.Kb[mpi.Base+uint32(mpi.Arr[mpi.LenArr-1])]
		var start, end int
		strand := mka.GetMapStrand()
		if strand == PLUS {
			start, end = int(mka.GetPosition()), int(mkb.GetPosition())+int(mkb.GetMapLen())
		} else {
			start, end = int(mkb.GetPosition()), int(mka.GetPosition())+int(mka.GetMapLen())
		}

		e := &edgesArr[mpi.EID]
		//bubbleIdx := -1
		//bubbleCount := 0
		//var bubblea, bubbleb uint32
		for j := i + 1; j < len(highQMPIArr); j++ {
			if priorIdxArr[j] == math.MaxUint16 {
				continue
			}
			nm := highQMPIArr[j]
			a, b := nm.Kb[nm.Base+uint32(nm.Arr[0])], nm.Kb[nm.Base+uint32(nm.Arr[nm.LenArr-1])]
			var nstart, nend int
			nstrand := a.GetMapStrand()
			if nstrand == PLUS {
				nstart, nend = int(a.GetPosition()), int(b.GetPosition())+int(b.GetMapLen())
			} else {
				nstart, nend = int(b.GetPosition()), int(a.GetPosition())+int(a.GetMapLen())
			}
			if nstart >= end {
				break
			}
			e2 := &edgesArr[nm.EID]
			//fmt.Printf("[GetMPIArrPrior]start:%d end:%d nstart:%d nend:%d\n", start, end, nstart, nend)
			if nstart == start && nend == end {
				if !IsBubble(e, e2, nodesArr) {
					priorIdxArr[i] = math.MaxUint16
					priorIdxArr[j] = math.MaxUint16
				}
			} else if nend <= end {
				priorIdxArr[i] = math.MaxUint16
				priorIdxArr[j] = math.MaxUint16
			}
		}
	}

	//compact priorIdxArr
	count := 0
	for i, idx := range priorIdxArr {
		if idx != math.MaxUint16 {
			priorIdxArr[count] = uint16(i)
			count++
		}
	}
	priorIdxArr = priorIdxArr[:count]

	// set anchorIdx
	anchorIdx := GetAnchorMPIIdx(highQMPIArr, priorIdxArr, edgesArr)

	return highQMPIArr, anchorIdx, priorIdxArr
}

func MapingONTReadsChunk(dbgKmerInfoPac DBGKmerInfoPac, cs <-chan ReadsBucketInfo, SeedLen int, opt optionsDDBG, GapCost, MaxGapLen int, wc chan<- string, edgesArr []DBGEdge, nodesArr []DBGNode, SeqType uint8, MaxPathLen, MaxBubbleSeqLen, BuckSize, processID int) {
	winSize := opt.WinSize
	kmerlen := opt.Kmer
	//MinChainScoreIdentityPercent := opt.MinChainScoreIdentityPercent
	logfn := opt.Prefix + "." + strconv.Itoa(processID) + ".log"
	fp, err := os.Create(logfn)
	if err != nil {
		log.Fatalf("[MapingONTReadsChunk] file %s create error, err: %v\n", logfn, err)
	}
	defer fp.Close()
	logfp := os.Stdout
	//logfp := bufio.NewWriter(fp)
	//defer logfp.Flush()
	//if err := logfp.Flush(); err != nil {
	//	log.Fatalf("[MapingONTReadsChunk] failed to flush file: %s, err: %v\n", logfn, err)
	//}
	//extendLen := 3 * kmerlen
	//TwoEdgeCyleMaxLen := 800
	//MaxEnlongationSeqLen := uint32(3 * kmerlen)
	//MaxFlankLen := 2000
	//MaxStackAllow := 5000
	LongEdgeMinLen := 600
	SeedNumCountArr := GetEdgeLowFreqSeedNumArr(edgesArr, LongEdgeMinLen)
	highRegPercentArr := make([]int, 10)
	var totalReadNum, notFoundSeedEdgeNum, mapOneEdgeNum, ExtendNum int
	bucketInterSecKmerInfoArr := make([]MapingKmerInfo, BuckSize*10)
	//buf := make([]MapingKmerInfo, 0, 10000)
	//EdgeSeedNumCountArr := make([]EdgeSeedNumInfo, 0, 1200)
	maxPathInfoArr := make([]MaxPathInfo, 0, 1000)
	//highQMPIArr := make([]MaxPathInfo, 0, 3000)
	//stk := make([]ExtendPathInfo, 0, 1000)
	//StkCleanSize := 100
	//maxScoreArr := make([]ScorePos, 0, 1000)
	//ScoreWinSize := 50
	//highQualEPathInfoArr := make([]ExtendPathInfo, 0, 10)
	for {
		rb, ok := <-cs
		if !ok {
			var ts string
			wc <- ts
			break
		}
		totalReadNum += len(rb.ReadsArr)
		startReadID := rb.ReadsArr[0].ID
		bucketInterSecKmerInfoArr = GetInterSectionKmerInfoB(dbgKmerInfoPac, rb.KmerSortArr, uint32(startReadID), uint16(SeedLen), bucketInterSecKmerInfoArr)
		sort.Sort(MapingKmerInfoReadIDArr(bucketInterSecKmerInfoArr))
		//loop process every read
		for i := 0; i < len(bucketInterSecKmerInfoArr); {
			id := bucketInterSecKmerInfoArr[i].GetMapLen()
			rID := startReadID + int(id)
			ri := rb.ReadsArr[id]
			rl := len(ri.Seq)
			fmt.Fprintf(logfp, "read[%d] length:%d %s\n", rID, rl, ri.Anotition)
			MQ := 20
			if len(ri.Anotition) > 3 {
				x, err := strconv.Atoi(string(ri.Anotition[3:]))
				if err != nil {
					log.Fatalf("[MapingONTReadsChunk]read[%d].Anotition:%s error!!!!\n", ri.ID)
				} else {
					MQ = x
				}
			}
			i1, idx, esArr, seedNumPercent := AddHighFreqKmer2(bucketInterSecKmerInfoArr, i, dbgKmerInfoPac, uint16(SeedLen), LongEdgeMinLen, SeedNumCountArr, winSize, kmerlen, MQ)
			if idx-i < 10 || len(esArr) == 0 {
				fmt.Fprintf(os.Stderr, "[MapingONTReadsChunk]read[%d] len(ka):%d < 10 or len(esArr):%d!!!!\n", ri.ID, idx-i, len(esArr))
				//sort.Sort(MapingKmerInfoEPosArr(bucketInterSecKmerInfoArr[i:idx]))
				//PrintMKIArr(bucketInterSecKmerInfoArr[i:idx])
				i = i1
				continue
			}
			if len(esArr) > 500 {
				fmt.Fprintf(os.Stderr, "[MapingONTReadsChunk]read[%d] len(enArr):%d > 500\n", ri.ID, len(esArr))
			}
			fmt.Fprintf(logfp, "seedNumPercent:%.2f\n", seedNumPercent)
			// set MapLen
			for j := i; j < idx; j++ {
				bucketInterSecKmerInfoArr[j].SetMapLen(uint32(SeedLen))
			}
			sort.Sort(MapingKmerInfoEPosArr(bucketInterSecKmerInfoArr[i:idx]))
			ka := GetChainBlocks(bucketInterSecKmerInfoArr[i:idx])
			//PrintMKIArr(ka)
			maxPathInfoArr = GetMaxChainBlockArr(ka, maxPathInfoArr, edgesArr, nodesArr, SeedLen, GapCost, MaxGapLen, kmerlen, SeqType, seedNumPercent, rl, logfp)
			fmt.Fprintf(logfp, "len(ka):%d len(maxPathInfoArr):%d\n", len(ka), len(maxPathInfoArr))
			if DebugModel {
				for x := range maxPathInfoArr {
					fmt.Fprintf(logfp, "mpiArr[%d]:%s\n", x, PrintMPI(&maxPathInfoArr[x], edgesArr))
				}
			}
			var anchorIdx int
			var priorIdxArr []uint16
			maxPathInfoArr = FilterMPIArr(maxPathInfoArr, kmerlen, rl, edgesArr, logfp)
			maxPathInfoArr, anchorIdx, priorIdxArr = GetMPIArrPrior(maxPathInfoArr, kmerlen, edgesArr, nodesArr)

			fmt.Fprintf(logfp, "priorIdxArr:%v  anchorIdx:%d\n", priorIdxArr, anchorIdx)
			if len(maxPathInfoArr) == 0 || anchorIdx < 0 {
				notFoundSeedEdgeNum++
				fmt.Fprintf(logfp, "readID[%d] not found high quality mapping region\n", rID)
				i = i1
				continue
			} else if len(maxPathInfoArr) == 1 {
				mpi := &maxPathInfoArr[0]
				start, end, _, _, _ := GetMPIInfo(mpi)
				if int(start) < 2*kmerlen && int(end) > rl-2*kmerlen {
					mapOneEdgeNum++
					fmt.Fprintf(logfp, "readID[%d] only map one long edge[%d]\n", rID, mpi.EID)
				}
				i = i1
				continue
			}
			anchorMPI := &maxPathInfoArr[anchorIdx]
			fmt.Fprintf(logfp, "len(mpiArr):%d anchorMPI:%s\n", len(maxPathInfoArr), PrintMPI(anchorMPI, edgesArr))

			/*// find DBG path
			{
				mapArr := FindMappingPath(seedInfoArr, maxPathInfoArr, anchorIdx, priorIdxArr, edgesArr, nodesArr, kmerlen, rl, buf, MinChainScoreIdentityPercent, winSize, logfp)
				if len(mapArr) > 0 {
					for _, pathMat := range mapArr {
						wc <- TransformPathMat2String(pathMat, maxPathInfoArr, int(ri.ID), rl, string(ri.Anotition))
					}
				}
			}*/
			i = i1
		}
	}
	fmt.Printf("[MapingONTReadsChunk] total reads num: %v, not found seed edge num: %v, map one edge num: %v, extend read num: %v\n", totalReadNum, notFoundSeedEdgeNum, mapOneEdgeNum, ExtendNum)
	for i, num := range highRegPercentArr {
		fmt.Printf("[MapingONTReadsChunk]highRegPercentArr[%d]:%d\n", i, num)
	}

	fmt.Printf("[MapingONTReadsChunk]cap(bucketInterSecKmerInfoArr):%d\n", cap(bucketInterSecKmerInfoArr))
	fmt.Printf("[MapingONTReadsChunk]maxPathInfoArr cap:%d\n", cap(maxPathInfoArr))
}

func DeconstructDBG(c cli.Command) {
	// check arguments
	opt, suc := checkArgsDDBG(c)
	if suc == false {
		log.Fatalf("[DeconstructDBG] check Arguments error, opt: %v\n", opt)
	}
	fmt.Printf("Arguments: %+v\n", opt)

	DebugModel = opt.Debug
	//BuckReadsNum = 1000
	EdgeMapKmerMin = 40
	GapCost := 1
	MaxGapLen := 2000
	SeqType := uint8(1)
	MaxPathLen := 10
	MaxBubbleSeqLen := 3 * opt.Kmer
	HighFreqKmerMin := 1000
	HighFreqKmerMax := 100000
	//Kmerlen = opt.Kmer

	if DebugModel {
		profileFn := opt.Prefix + ".decdbg.prof"
		cpuprofilefp, err := os.Create(profileFn)
		if err != nil {
			log.Fatalf("[DeconstructDBG] open cpuprofile file: %v failed\n", profileFn)
		}
		pprof.StartCPUProfile(cpuprofilefp)
		defer pprof.StopCPUProfile()
	}

	// read files and construt DBG
	DBGInfofn := opt.Prefix + ".MapDBG.DBGInfo"
	eSize, nSize := DBGInfoReader(DBGInfofn)
	nodesfn := opt.Prefix + ".nodes.MapDBG.Arr"
	nodesArr := make([]DBGNode, nSize)
	fc := make(chan int, 1)
	go NodesArrReader(nodesfn, nodesArr, opt.Kmer, fc)
	//if len(nodesArr) != nSize {
	//	log.Fatalf("[DeconstructDBG] len(nodesArr): %v != nodesArr Size: %v in file: %v\n", len(nodesArr), nSize, DBGInfofn)
	//}
	//edgesArr := make([]DBGEdge, eSize)
	edgesfn := opt.Prefix + ".edges.MapDBG.fa.zst"
	edgesArr := ReadEdgesFromFile(edgesfn, uint32(eSize), opt.Kmer)
	<-fc

	go CheckInterConnectivity(edgesArr, nodesArr)

	AddNodeInfo2DBGEdgeArr(edgesArr, nodesArr)
	uniqueNum, bubbleRepeatNum, twoCycleNum, selfCycle, selfCycleSameComingNum, bubbleEdgeNum := SetDBGEdgesUniqueFlag(edgesArr, nodesArr, opt.Kmer)
	//nodesArr = nil
	fmt.Printf("[DeconstructDBG] unique edge number is: %v, bubbleRepeatNum: %v, twoCycleNum:%v, selfCycle:%v, selfCycleSameComingNum: %v, bubbleEdgeNum: %v\n", uniqueNum, bubbleRepeatNum, twoCycleNum, selfCycle, selfCycleSameComingNum, bubbleEdgeNum)
	runtime.GC()
	t0 := time.Now()

	PrintTmpNodesArr(nodesArr, opt.Prefix)

	pathFn := opt.Prefix + ".Path.zst"
	mapDBGInfoFn := opt.Prefix + ".ONTMapDBGInfo"
	// no any other align tools for seed edge by ONT reads, use same as daligner,
	// first sort DBG edges by seed kmer, and split ONT reads by bucket, every bucket sort by seed kmer, and found share kmers from edges sort kmersArr
	if _, err := os.Stat(pathFn); err != nil {
		//BucketSize := (1 << 28) // ==2**28
		cfgInfo, err := ParseCfg(opt.CfgFn, false, false)
		if err != nil {
			log.Fatalf("[DeconstructDBG] ParseCfg 'C': %v err :%v\n", opt.CfgFn, err)
		}
		fmt.Printf("[DeconstructDBG]cfgInfo:%+v\n", cfgInfo)

		minimerSize := GetConstructDBGSampleSize(edgesArr, opt.WinSize, SeedLen, opt.Kmer)
		//fmt.Printf("[DeconstructDBG]total found edges minimizers number is: %v\n", len(edgesKmerSortArr))
		//bufSize := int64(minimerSize) / 50
		//bufSize := int64(1000000)

		cs := make(chan ReadsBucketInfo, 2)
		go ParaLoadLongReadsAndSort(cfgInfo, cs, SeedLen, opt.WinSize, BuckSize)
		dbgCS := make(chan DBGKmerInfoPac, 1)
		SortDBGEdgesKmer(dbgCS, edgesArr, SeedLen, opt, HighFreqKmerMin, HighFreqKmerMax, minimerSize)
		dbgKmerInfoPac := <-dbgCS
		runtime.GC()
		fmt.Printf("[DeconstructDBG] Sort DBGEdge used: %v\n", time.Now().Sub(t0))
		//go PrintDBGEdges(edgesArr, opt.Prefix)
		t0 = time.Now()
		wc := make(chan string, 1000)
		for i := 0; i < opt.NumCPU-1; i++ {
			go MapingONTReadsChunk(dbgKmerInfoPac, cs, SeedLen, opt, GapCost, MaxGapLen, wc, edgesArr, nodesArr, SeqType, MaxPathLen, MaxBubbleSeqLen, BuckSize, i)
		}
		//func FindONTReadsPath(edgesKmerSortArr []KI, cs <-chan ReadsBucketInfo, SeedLen, winSize, GapCost, MaxGapLen, kmerlen int, wc chan<- LongReadMappingInfo, edgesArr []DBGEdge, nodesArr []DBGNode, SeqType uint8) {

		//pathFn := opt.Prefix + ".Path.protobuf"
		WriteLongPathToFile(wc, pathFn, mapDBGInfoFn, opt.NumCPU-1)
		/*if !DebugModel {
			profileFn := opt.Prefix + ".decdbg.Mem.prof"
			f, err := os.Create(profileFn)
			if err != nil {
				log.Fatal(err)
			}
			pprof.WriteHeapProfile(f)
			f.Close()
		}*/
		fmt.Printf("[DeconstructDBG] Maping long reads to DBG used: %v\n", time.Now().Sub(t0))
	}
	// remap NGS reads to the new samplify DBG
	/*var copt optionsDDBG
	copt.SeedInfo = opt.Kmer
	copt.NumCPU = opt.NumCPU
	copt.WinSize = opt.WinSize
	copt.MaxNGSReadLen = opt.MaxNGSReadLen
	copt.CfgFn = opt.CfgFn
	wrFn := opt.Prefix + ".decdbg.NGSAlignment"
	MapNGS2DBG(copt, nodesArr, edgesArr, wrFn)
	AddPathToDBGEdge(edgesArr, wrFn)
	MergePathMat(edgesArr, nodesArr, opt.MinMapFreq)*/

	// get ont reads mapping info by minimap2
	//paffn := opt.Prefix + ".paf"

	// get ont Long reads Mapping info by minimap2, use smfy edges for reference, ONT reads for query
	//bufSize := 5000
	/*if _, err := os.Stat(pathfn); err != nil {
		paffn := opt.Prefix + ".paf"
		rc := make(chan []PAFInfo, opt.NumCPU*bufSize)
		wc := make(chan LongReadMappingInfo, opt.NumCPU*bufSize)

		go GetPAFRecord(paffn, opt.ONTFn, rc, opt.NumCPU)

		for i := 0; i < opt.NumCPU; i++ {
			go paraFindLongReadsMappingPath(rc, wc, edgesArr, nodesArr, opt)
		}

		WriteLongPathToFile(wc, pathfn, opt.NumCPU)
		fmt.Printf("[DeconstructDBG] mapping long reads to DBG edges used: %v\n", time.Now().Sub(t0))
	}*/

	/*t1 := time.Now()
	{
		cs := make(chan MapingInfo, 100)
		go LoadPathFromFile(pathFn, cs)
		// Simplify using Long Reads Mapping info
		mapRecordSize := LoadONTMapInfoFile(mapDBGInfoFn)
		//joinPathArr := SimplifyByLongReadsPath(edgesArr, nodesArr, cs, opt, MaxPathLen, MaxBubbleSeqLen, mapRecordSize)

		//graphfn := opt.Prefix + ".afterLR.dot"
		//GraphvizDBGArr(nodesArr, edgesArr, graphfn)
		DcDBGEdgesfn := opt.Prefix + ".edges.DcDBG.fa"
		ExtractSeq(edgesArr, nodesArr, joinPathArr, DcDBGEdgesfn, opt.Kmer)
		fmt.Printf("[DeconstructDBG] find max edge path and produce sequence used: %v\n", time.Now().Sub(t1))
	}*/
	fmt.Printf("[DeconstructDBG] total used: %v\n", time.Now().Sub(t0))
	//StoreEdgesToFn(DcDBGEdgesfn, edgesArr)
}
