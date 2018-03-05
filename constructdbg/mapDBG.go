package constructdbg

import (
	// "bufio"
	// "compress/gzip"
	// "encoding/binary"
	"fmt"
	"reflect"

	"github.com/mudesheng/GA/bnt"
	"github.com/mudesheng/GA/constructcf"
	//"github.com/mudesheng/GA/constructdbg"
	// "ga/cuckoofilter"
	"io"
	"log"

	"github.com/biogo/biogo/alphabet"
	"github.com/biogo/biogo/io/seqio/fasta"
	"github.com/biogo/biogo/seq/linear"
	"github.com/jwaldrip/odin/cli"
	// "math"
	"os"
	"runtime"
	"sort"
	"strconv"
	// "strings"
	// "time"
)

const (
	AllState         = 1
	ScaffState       = 2
	GapState         = 3
	EndString        = "end"
	ReadSeqSize      = 100
	ErrorRate        = 0.15
	MIN_KMER_MAP_LEN = 60
)

//var Kmerlen = -1

type SeedKmerInfo struct {
	Kmer uint32
	Info uint64 // compact of eID, start position of index and strand(0: plus, 1: minus)
}

type Seq struct {
	ID int
	S  []byte
}

type MergeKmerInfo struct {
	QID, RID     uint32
	Qinfo, Rinfo uint32 // info compact start position of index and strand(0:plus, 1:minus)
}

func GetEdge(edgesArr []DBGEdge, rc chan Seq, numCPU int) {
	for _, e := range edgesArr {
		if e.ID == 0 {
			continue
		}
		var seq Seq
		seq.ID = int(e.ID)
		seq.S = e.Utg.Ks
		rc <- seq
	}
	// seed terminal signals
	for i := 0; i < numCPU; i++ {
		var ns Seq
		rc <- ns
	}
}

type SkSlice []SeedKmerInfo

func CollectSeedKmerInfo(wc chan SeedKmerInfo, numCPU int) (skSlice []SeedKmerInfo) {
	terminalNum := 0
	for {
		skInfo := <-wc
		if skInfo.Info == 0 {
			terminalNum++
			if terminalNum == numCPU {
				break
			} else {
				continue
			}
		}
		// fmt.Printf("[CollectSeedKmerInfo] ID: %d\n", skInfo.Info>>32)

		skSlice = append(skSlice, skInfo)

	}

	return skSlice
}

func GetKmer(seq []byte, len int) (kmer uint32) {
	for i := 0; i < len; i++ {
		kmer <<= bnt.NumBitsInBase
		kmer |= uint32(seq[i])
	}

	return kmer
}

func GetRCSeed(k uint32, size int) (rk uint32) {
	k = ^k
	for i := 0; i < size; i++ {
		rk <<= bnt.NumBitsInBase
		rk |= (k & 0x3)
		k >>= bnt.NumBitsInBase
	}

	return rk
}

func GetMinSeed(k uint32, seed int) (uint32, uint64) {
	rk := GetRCSeed(k, seed)
	if rk < k {
		return rk, 1
	} else {
		return k, 0
	}
}

// compact by 32,31,1
func CompactInfo(ID, pos int, strand uint64) (ci uint64) {
	ci = uint64(ID) << 32
	ci |= (uint64(pos) << 1)
	ci |= strand

	return ci
}

func GetMinKmerFromWidth(minArr []SeedKmerInfo, width int) (min SeedKmerInfo, p int) {
	p = 0
	min = minArr[p]
	for i := 1; i < len(minArr); i++ {
		if minArr[i].Kmer < min.Kmer {
			min = minArr[i]
			p = i
		}
	}

	return min, p
}

func paraCollectMinKmer(rc chan Seq, wc chan SeedKmerInfo, seed int, width int) {

	for {
		data := <-rc
		seq := data.S
		if len(seq) == 0 {
			var ski SeedKmerInfo
			wc <- ski
			break
		}
		var minArr []SeedKmerInfo
		w := width
		for i := 0; i < len(seq)-seed+1; i++ {
			k := GetKmer(seq[i:], seed)
			mk, strand := GetMinSeed(k, seed)
			var min SeedKmerInfo
			min.Kmer = mk
			min.Info = CompactInfo(data.ID, i, strand)
			// fmt.Printf("[paraCollectMinKmer] ID: %d\n", min.Info>>32)
			minArr = append(minArr, min)
			if i > 0 && i%128 == 0 {
				for j := 0; j < width; j++ {
					minArr[j] = minArr[len(minArr)-width+j]
				}
				minArr = minArr[:width]
			}
			if len(minArr) >= width {
				if w < width {
					if min.Kmer < minArr[len(minArr)-1-w].Kmer {
						wc <- min
						w = 0
					}
				} else {
					m, p := GetMinKmerFromWidth(minArr[len(minArr)-width:], width)
					wc <- m
					w = width - 1 - p
				}
				w++
			}

		}
	}
}

const B = 8 // test B get better performance
const SIZE = (1 << 8)
const MASK = SIZE - 1

func SumBucket(src []SeedKmerInfo, bucket [][SIZE]int, B, part int) {
	num := len(bucket)
	// get part summary
	for i := 0; i < num; i++ {
		end := (i + 1) * part
		if i == num-1 {
			end = len(src)
		}
		for j := i * part; j < end; j++ {
			c := src[j].Kmer
			// b := c >> begin
			// if b&MASK == 0 {
			// 	fmt.Printf("[RadixSort]b: %v\tsrc[%d]:%v\n", b, j, src[j])
			// }
			bucket[i][c&MASK]++
		}
	}

	// accumulative total sum
	sum := 0
	for i := 0; i < SIZE; i++ {
		for j := 0; j < num; j++ {
			tmp := bucket[j][i]
			bucket[j][i] = sum
			sum += tmp
		}
	}
	fmt.Printf("[SumBucket] sum: %v\tlen(src): %v\n", sum, len(src))
}

func SumNext(next [][][SIZE]int, bucket [][SIZE]int) {
	num := len(bucket)
	// sum threads
	for i := 0; i < num; i++ {
		for j := 0; j < SIZE; j++ {
			bucket[i][j] = 0
			for z := 0; z < num; z++ {
				bucket[i][j] += next[z][i][j]
			}
		}
	}

	// accumulative total
	sum := 0
	for i := 0; i < SIZE; i++ {
		for j := 0; j < num; j++ {
			tmp := bucket[j][i]
			bucket[j][i] = sum
			sum += tmp
		}
	}
	fmt.Printf("[SumNext] sum: %v\n", sum)

}

func InitNext(next [][][SIZE]int) {
	for i := 0; i < len(next); i++ {
		for j := 0; j < len(next[i]); j++ {
			for z := 0; z < SIZE; z++ {
				next[i][j][z] = 0
			}
		}
	}
}

func paraBucketSort(src, trg []SeedKmerInfo, bucket *[SIZE]int, next [][SIZE]int, begin uint64, part int, wc chan int) {
	for _, c := range src {
		b := c.Kmer >> begin
		// fmt.Printf("[paraBucketSort] b: %v\n", b)
		s := b & MASK
		x := bucket[s]
		bucket[s]++
		trg[x] = c
		// fmt.Printf("[paraBucketSort] %d/%d: %v\n", x, part, x/part)
		next[x/part][(b>>B)&MASK]++
	}

	wc <- 1
}

func RadixSort(src []SeedKmerInfo, hbits, numCPU int) []SeedKmerInfo {
	trg := make([]SeedKmerInfo, len(src))
	part := (len(src)-1)/numCPU + 1
	bucket := make([][SIZE]int, numCPU)
	next := make([][][SIZE]int, numCPU)
	for i := 0; i < len(next); i++ {
		next[i] = make([][SIZE]int, numCPU)
	}

	for i := 0; i < hbits; i += B {
		if i == 0 {
			SumBucket(src, bucket, B, part)
		} else {
			SumNext(next, bucket)
		}
		// fmt.Printf("[RadixSort] part: %v\nbucket: %v\n", part, bucket)
		InitNext(next)
		wc := make(chan int)
		for j := 0; j < numCPU; j++ {
			if j == numCPU-1 {
				go paraBucketSort(src[j*part:len(src)], trg, &bucket[j], next[j], uint64(i), part, wc)
			} else {
				go paraBucketSort(src[j*part:(j+1)*part], trg, &bucket[j], next[j], uint64(i), part, wc)
			}
		}

		// wait for sort threads finished
		for j := 0; j < numCPU; j++ {
			<-wc
		}
		// fmt.Printf("[RadixSort] next: %v\n", next)

		src, trg = trg, src
	}

	return src
}

func GetRawReads(rc chan []Seq, cfgInfo constructcf.CfgInfo, blockSize int) {
	var slice []Seq
	ID := 0
	num := 0
	// fmt.Printf("[GetRawReads] cfgInfo: %v\n", cfgInfo)
	for _, lib := range cfgInfo.Libs {
		// seqProfile == 2 note Pacbio
		if lib.SeqProfile == 2 {
			for _, fn := range lib.FnName {
				infile, err := os.Open(fn)
				if err != nil {
					log.Fatal(err)
				}
				defer infile.Close()
				// gzreader, err := gzip.NewReader(infile)
				// bufin := bufio.NewReader(gzreader)
				fafp := fasta.NewReader(infile, linear.NewSeq("", nil, alphabet.DNA))
				for {
					if s, err := fafp.Read(); err != nil {
						if err == io.EOF {
							break
						} else {
							log.Fatalf("[GetRawReads] read file: %s error: %v\n", fn, err)
						}
					} else {
						l := s.(*linear.Seq)
						// fmt.Printf("[GetRawReads] l: %v\n", l)
						var seq Seq
						seq.ID = ID
						ID++
						seq.S = make([]byte, len(l.Seq))
						for j, v := range l.Seq {
							seq.S[j] = bnt.Base2Bnt[v]
						}
						slice = append(slice, seq)
						num += len(seq.S)
						if num > blockSize {
							rc <- slice
							var ns []Seq
							slice = ns
							num = 0
						}

					}
				}
			}
		}
	}

	rc <- slice
	var ns []Seq
	rc <- ns
}

func GetSeq(seqSlice []Seq, rc chan Seq, numCPU int) {
	for _, s := range seqSlice {
		rc <- s
		// fmt.Printf("[GetSeq]ID: %d\n", s.ID)
	}

	for i := 0; i < numCPU; i++ {
		var s Seq
		rc <- s
	}
}

func FindRefStart(skArr []SeedKmerInfo, initPos int, kmer uint32) int {
	if skArr[initPos].Kmer < kmer {
		for i := initPos + 1; i < len(skArr); i++ {
			if skArr[i].Kmer >= kmer {
				return i
			}
		}
	} else { // skArr[initPos].Kmer <= kmer
		for i := initPos; i >= 0; i-- {
			if skArr[i].Kmer == kmer {
				for j := i - 1; j >= 0; j-- {
					if skArr[i].Kmer < kmer {
						return j + 1
					}
				}
			} else if skArr[i].Kmer < kmer {
				return i + 1
			}
		}
	}

	// unreachable state
	log.Fatalf("[FindRefStart] not found start\n")

	return -1
}

func FindRefEnd(skArr []SeedKmerInfo, start int) (end int) {
	initK := skArr[start].Kmer
	end = len(skArr)
	for i := start + 1; i < len(skArr); i++ {
		if skArr[i].Kmer != initK {
			end = i
			break
		}
	}

	return end
}

const MAX_UINT32 = (1 << 32) - 1

func SplitFromSeedKmer(ski SeedKmerInfo) (id, info uint32) {
	info = uint32(ski.Info & MAX_UINT32)
	id = uint32(ski.Info >> 32)

	return id, info
}

func paraMergeRefAndQuy(refSlice, MinQuery []SeedKmerInfo, begin, end int, wc chan MergeKmerInfo) {
	refStart := 0
	if begin > 0 {
		refStart = FindRefStart(refSlice, (begin*len(refSlice))/len(MinQuery), MinQuery[begin].Kmer)
	}

	rs, re := refStart, refStart
	for i := begin; i < end; i++ {
		if MinQuery[i].Kmer < refSlice[rs].Kmer {
			continue
		} else if MinQuery[i].Kmer > refSlice[rs].Kmer {
			for j := re; j < len(refSlice); j++ {
				if MinQuery[i].Kmer <= refSlice[j].Kmer {
					rs = j
					break
				}
			}
		}
		if MinQuery[i].Kmer == refSlice[rs].Kmer {
			if re <= rs {
				re = FindRefEnd(refSlice, rs)
			}
			QID, Qinfo := SplitFromSeedKmer(MinQuery[i])
			for j := rs; j < re; j++ {
				// if MinQuery[i].Kmer != refSlice[j].Kmer {
				// 	log.Fatalf("[MergeRadixSort] MinQuery != refSlice\n")
				// }
				var mk MergeKmerInfo
				mk.QID, mk.Qinfo = QID, Qinfo
				mk.RID, mk.Rinfo = SplitFromSeedKmer(refSlice[j])
				wc <- mk
			}
		}
	}

	var n MergeKmerInfo
	wc <- n
}

func CollectMergeKmerInfo(wc chan MergeKmerInfo, numCPU int) (arr []MergeKmerInfo) {
	terminalNum := 0
	for {
		mki := <-wc
		if mki.RID == 0 {
			terminalNum++
			if terminalNum == numCPU {
				break
			} else {
				continue
			}
		}
		arr = append(arr, mki)
	}

	return arr
}

func SumMergedBucket(src []MergeKmerInfo, bucket [][SIZE]int, B, part int) {
	num := len(bucket)
	// get part summary
	for i := 0; i < num; i++ {
		end := (i + 1) * part
		if i == num-1 {
			end = len(src)
		}
		for j := i * part; j < end; j++ {
			c := src[j].Qinfo
			// b := c >> begin
			// if b&MASK == 0 {
			// 	fmt.Printf("[RadixSort]b: %v\tsrc[%d]:%v\n", b, j, src[j])
			// }
			bucket[i][c&MASK]++
		}
	}

	// accumulative total sum
	sum := 0
	for i := 0; i < SIZE; i++ {
		for j := 0; j < num; j++ {
			tmp := bucket[j][i]
			bucket[j][i] = sum
			sum += tmp
		}
	}
	fmt.Printf("[SumMergedBucket] sum: %v\tlen(src): %v\n", sum, len(src))

}

func paraQinfoBucketSort(src, trg []MergeKmerInfo, bucket *[SIZE]int, next [][SIZE]int, begin uint64, part int, wc chan int) {
	if begin < 24 {
		for _, c := range src {
			b := c.Qinfo >> begin
			// fmt.Printf("[paraBucketSort] b: %v\n", b)
			s := b & MASK
			x := bucket[s]
			bucket[s]++
			trg[x] = c
			// fmt.Printf("[paraBucketSort] %d/%d: %v\n", x, part, x/part)
			next[x/part][(b>>B)&MASK]++
		}

	} else {
		for _, c := range src {
			b := c.Qinfo >> begin
			// fmt.Printf("[paraBucketSort] b: %v\n", b)
			s := b & MASK
			x := bucket[s]
			bucket[s]++
			trg[x] = c
			// fmt.Printf("[paraBucketSort] %d/%d: %v\n", x, part, x/part)
			next[x/part][c.RID&MASK]++
		}

	}

	wc <- 1

}

func paraRIDBucketSort(src, trg []MergeKmerInfo, bucket *[SIZE]int, next [][SIZE]int, begin uint64, part int, wc chan int) {
	if begin < 24 {
		for _, c := range src {
			b := c.RID >> begin
			// fmt.Printf("[paraBucketSort] b: %v\n", b)
			s := b & MASK
			x := bucket[s]
			bucket[s]++
			trg[x] = c
			// fmt.Printf("[paraBucketSort] %d/%d: %v\n", x, part, x/part)
			next[x/part][(b>>B)&MASK]++
		}

	} else {
		for _, c := range src {
			b := c.RID >> begin
			// fmt.Printf("[paraBucketSort] b: %v\n", b)
			s := b & MASK
			x := bucket[s]
			bucket[s]++
			trg[x] = c
			// fmt.Printf("[paraBucketSort] %d/%d: %v\n", x, part, x/part)
			next[x/part][c.QID&MASK]++
		}

	}

	wc <- 1

}

func paraQIDBucketSort(src, trg []MergeKmerInfo, bucket *[SIZE]int, next [][SIZE]int, begin uint64, part int, wc chan int) {
	if begin < 24 {
		for _, c := range src {
			b := c.QID >> begin
			// fmt.Printf("[paraBucketSort] b: %v\n", b)
			s := b & MASK
			x := bucket[s]
			bucket[s]++
			trg[x] = c
			// fmt.Printf("[paraBucketSort] %d/%d: %v\n", x, part, x/part)
			next[x/part][(b>>B)&MASK]++
		}

	} else {
		for _, c := range src {
			b := c.QID >> begin
			// fmt.Printf("[paraBucketSort] b: %v\n", b)
			s := b & MASK
			x := bucket[s]
			bucket[s]++
			trg[x] = c
			// fmt.Printf("[paraBucketSort] %d/%d: %v\n", x, part, x/part)
			// next[x/part][c.RID&MASK]++
		}

	}

	wc <- 1

}

func MergeRadixSort(src []MergeKmerInfo, numCPU int) {
	trg := make([]MergeKmerInfo, len(src))
	part := (len(src)-1)/numCPU + 1
	bucket := make([][SIZE]int, numCPU)
	next := make([][][SIZE]int, numCPU)
	for i := 0; i < len(next); i++ {
		next[i] = make([][SIZE]int, numCPU)
	}

	// sort by Qinfo
	hbits := 32
	for i := 0; i < hbits; i += B {
		if i == 0 {
			SumMergedBucket(src, bucket, B, part)
		} else {
			SumNext(next, bucket)
		}
		// fmt.Printf("[MergeRadixSort] part: %v\nbucket: %v\n", part, bucket)
		InitNext(next)
		wc := make(chan int)
		for j := 0; j < numCPU; j++ {
			if j == numCPU-1 {
				go paraQinfoBucketSort(src[j*part:len(src)], trg, &bucket[j], next[j], uint64(i), part, wc)
			} else {
				go paraQinfoBucketSort(src[j*part:(j+1)*part], trg, &bucket[j], next[j], uint64(i), part, wc)
			}
		}

		// wait for sort threads finished
		for j := 0; j < numCPU; j++ {
			<-wc
		}
		// fmt.Printf("[RadixSort] next: %v\n", next)

		src, trg = trg, src
	}

	// sort by RID
	for i := 0; i < hbits; i += B {
		SumNext(next, bucket)
		InitNext(next)
		wc := make(chan int)
		for j := 0; j < numCPU; j++ {
			if j == numCPU-1 {
				go paraRIDBucketSort(src[j*part:len(src)], trg, &bucket[j], next[j], uint64(i), part, wc)
			} else {
				go paraRIDBucketSort(src[j*part:(j+1)*part], trg, &bucket[j], next[j], uint64(i), part, wc)
			}
		}

		// wait for sort threads finished
		for j := 0; j < numCPU; j++ {
			<-wc
		}
		// fmt.Printf("[RadixSort] next: %v\n", next)

		src, trg = trg, src
	}

	// sort by QID
	for i := 0; i < hbits; i += B {
		SumNext(next, bucket)
		InitNext(next)
		wc := make(chan int)
		for j := 0; j < numCPU; j++ {
			if j == numCPU-1 {
				go paraQIDBucketSort(src[j*part:len(src)], trg, &bucket[j], next[j], uint64(i), part, wc)
			} else {
				go paraQIDBucketSort(src[j*part:(j+1)*part], trg, &bucket[j], next[j], uint64(i), part, wc)
			}
		}

		// wait for sort threads finished
		for j := 0; j < numCPU; j++ {
			<-wc
		}
		// fmt.Printf("[RadixSort] next: %v\n", next)

		src, trg = trg, src
	}
}

func GetMki(mkiArr []MergeKmerInfo, rc chan []MergeKmerInfo, numCPU int) {
	// fmt.Printf("[GetMki]len(mkiArr): %v\n", len(mkiArr))
	var slc []MergeKmerInfo
	slc = append(slc, mkiArr[0])
	smkiArr := mkiArr[1:]
	for _, m := range smkiArr {
		// fmt.Printf("[GetMki]m: %v\n", m)
		if m.QID != slc[len(slc)-1].QID {
			if len(slc) > 10 {
				rc <- slc
				// fmt.Printf("[GetMki]slc: %v\n", slc)
			}
			var n []MergeKmerInfo
			slc = n
		}
		slc = append(slc, m)
	}
	if len(slc) > 10 {
		rc <- slc
	}

	// send terminal signal
	for i := 0; i < numCPU; i++ {
		var n []MergeKmerInfo
		rc <- n
	}
}

type MapEdgeKmers struct {
	EID                    uint32
	MapBN                  int32 // the total mapping base number
	Strand                 uint8 // the Ref strand ^ Query strand
	Index, Max, Start, End int32 // Index and Max note start and end index of region, Start and End note QPos
}

type MEKArr []MapEdgeKmers

func (arr MEKArr) Len() int {
	return len(arr)
}

func (arr MEKArr) Less(i, j int) bool {
	return arr[i].MapBN < arr[j].MapBN
}

func (arr MEKArr) Swap(i, j int) {
	arr[i], arr[j] = arr[j], arr[i]
}

func SummaryRegionMapBaseNum(noContainMk MapEdgeKmers, start, end, seed int, mks []MergeKmerInfo) (mapBaseN int) {
	lasts := 0 - 2*seed
	for i := noContainMk.Index; ; i++ {
		mk := mks[i]
		if mk.RID == 0 {
			continue
		} else if mk.RID != noContainMk.EID {
			break
		}
		s := int(mk.Qinfo >> 1)
		if start <= s && s <= end-seed {
			if s-lasts >= seed {
				mapBaseN += seed
			} else {
				mapBaseN += (s - lasts)
			}
		} else if s > end-seed {
			break
		}
		lasts = s
	}

	return mapBaseN
}

func FoundMaxUniqueRegion(depth []uint8) (start, end int) {
	s, e := -1, -1
	for i, c := range depth {
		if c != 1 {
			if s >= 0 && e > s && e+1-s > end-start {
				start = s
				end = e + 1
			}
			s, e = -1, -1
			continue
		}
		if s < 0 {
			s = i
		}
		e = i

	}

	return start, end
}

func GetSeedEdge(NoContainMKArr []MapEdgeKmers, start, end int32) (index, overlapLen int) {
	for i, m := range NoContainMKArr {

		maxs, mine := start, end
		if m.Start > maxs {
			maxs = m.Start
		}
		if m.End < mine {
			mine = m.End
		}
		if int(mine-maxs) > overlapLen {
			index = i
			overlapLen = int(mine - maxs)
		}
	}

	return index, overlapLen
}

func GetNextEArr(eID DBG_MAX_INT, node DBGNode) (eIDArr []DBG_MAX_INT) {
	var direction uint8
	for i := 0; i < bnt.BaseTypeNum; i++ {
		if node.EdgeIDIncoming[i] == eID {
			direction = BACKWARD
			break
		}
		if node.EdgeIDOutcoming[i] == eID {
			direction = FORWARD
			break
		}
	}
	if direction != BACKWARD && direction != FORWARD {
		log.Fatalf("[GetNextEArr] direction not set, eID:%d, node:%v\n", eID, node)
	}
	for i := 0; i < bnt.BaseTypeNum; i++ {
		if direction == FORWARD {
			if node.EdgeIDIncoming[i] > 0 {
				eIDArr = append(eIDArr, node.EdgeIDIncoming[i])
			}
		} else {
			if node.EdgeIDOutcoming[i] > 0 {
				eIDArr = append(eIDArr, node.EdgeIDOutcoming[i])
			}
		}
	}
	return eIDArr
}

func FindMEKArr(NoContainMKArr []MapEdgeKmers, eIDArr []DBG_MAX_INT, qPos int32) (proArr []MapEdgeKmers) {
	for _, eID := range eIDArr {
		var arr []MapEdgeKmers
		for _, m := range NoContainMKArr {
			if DBG_MAX_INT(m.EID) == eID {
				arr = append(arr, m)
			}
		}
		if len(arr) > 1 {
			min := int32(1000000)
			idx := -1
			for x, mk := range arr {
				if qPos < mk.Start {
					tm := mk.Start - qPos
					if tm < min {
						min = tm
						idx = x
					}
				} else if qPos <= mk.End {
					if mk.End-qPos < qPos-mk.Start && mk.End-qPos < min {
						min = mk.End - qPos
						idx = x
					} else if mk.End-qPos > qPos-mk.Start && qPos-mk.Start < min {
						min = qPos - mk.Start
						idx = x
					}
				} else {
					if mk.End-qPos < min {
						min = mk.End - qPos
						idx = x
					}
				}
			}

			if idx >= 0 {
				proArr = append(proArr, arr[idx])
			}
		} else if len(arr) == 1 {
			proArr = append(proArr, arr[0])
		}
	}
	return proArr
}

func GetRefFlankIdx(marr []int32, q int32, direction uint8) (r int32) {
	r = -1
	if direction == FORWARD {
		for i, m := range marr {
			if m == q {
				r = int32(i)
				break
			}
		}
	} else {
		for i := len(marr) - 1; i >= 0; i-- {
			if marr[i] == q {
				r = int32(i)
				break
			}
		}
	}

	if r < 0 {
		log.Fatalf("[GetRefFlankIdx] r < 0\n")
	}

	return r
}

func GetNoSeedRegionEnd(marr []int32, rs int, QEnd int32) (re, qe int) {
	for i := rs; i < len(marr); i++ {
		if marr[i] > 0 {
			re = i
			qe = int(marr[i])
			return re, qe
		}
	}

	re = len(marr)
	qe = int(QEnd)
	return re, qe
}

func GetNoSeedRegionStart(marr []int32, re int, QS int32) (rs, qs int) {
	for i := re; i >= 0; i-- {
		if marr[i] > 0 {
			rs = i
			qs = int(marr[i])
			return rs, qs
		}
	}

	rs = 0
	qs = int(QS)
	return rs, qs
}

func DynamicAlign(cSeq, readSeq []byte) (matchNum, inst, del, mis int32) {
	if len(cSeq) > 255 && len(readSeq) > 255 {
		log.Fatalf("[DynamicAlign] just allow align length <256, len(cSeq): %d, len(readSeq): %d\n", len(cSeq), len(readSeq))
	}

	fmt.Printf("len(cSeq): %d, len(readSeq): %d\n", len(cSeq), len(readSeq))
	maxScore, x, y := uint8(0), 0, 0
	scoreArr := make([][]uint8, len(cSeq)+1)
	for i, _ := range scoreArr {
		scoreArr[i] = make([]uint8, len(readSeq)+1)
	}

	for i := 1; i < len(cSeq)+1; i++ {
		offDiagonal := int(ErrorRate*float32(i)) + 5
		min, max := 1, len(readSeq)+1
		if i-offDiagonal > min {
			min = i - offDiagonal
		}
		if i+offDiagonal < max {
			max = i + offDiagonal
		}
		// fmt.Printf("min: %d, max: %d\n", min, max)
		for j := min; j < max; j++ {
			sc := scoreArr[i-1][j]
			if sc < scoreArr[i][j-1] {
				sc = scoreArr[i][j-1]
			}
			if cSeq[i-1] == readSeq[j-1] && scoreArr[i-1][j-1]+1 > sc {
				sc = scoreArr[i-1][j-1] + 1
				if sc > maxScore {
					maxScore = sc
					x, y = i, j
				}
			}
			scoreArr[i][j] = sc
		}
	}
	fmt.Printf("    %v\n", readSeq)
	for i := 0; i <= len(cSeq); i++ {
		if i > 0 {
			fmt.Printf("%v|", cSeq[i-1])
		} else {
			fmt.Printf("  ")
		}
		fmt.Printf("%v\n", scoreArr[i])

	}
	fmt.Printf("x: %d y: %d\n", x, y)

	return matchNum, inst, del, mis
}

func MatchSeq(readSeq []byte, mk MapEdgeKmers, mks []MergeKmerInfo, edgeSeq []byte, seed uint32, direction uint8) (matchNum int32) {
	cSeq := edgeSeq
	if mk.Strand == 1 {
		tmp := make([]byte, len(edgeSeq))
		copy(tmp, edgeSeq)
		ReverseCompByteArr(tmp)
		cSeq = tmp
	}
	fmt.Printf("[MatchSeq]mk: %v\n", mk)
	marr := make([]int32, len(edgeSeq))
	for i := mk.Index; i <= mk.Max; i++ {
		if mks[i].RID == 0 {
			continue
		}
		var rs uint32
		if mk.Strand == 0 {
			rs = mks[i].Rinfo >> 0x1
		} else {
			rs = uint32(len(cSeq)) - 1 - (mks[i].Rinfo >> 0x1) - (seed - 1)
		}
		// fmt.Printf("[MatchSeq]mks[%d]: %v\trs: %v\tlen(marr): %v\n", i, mks[i], rs, len(marr))
		marr[rs] = int32(mks[i].Qinfo >> 0x1)
		marr[rs+seed-1] = marr[rs] + int32(seed) - 1
		for j := rs + 1; j < rs+seed-1; j++ {
			if marr[j] == 0 {
				marr[j]++
			}
		}
	}

	max_gap_len := int32(3)
	ExtendSeed := int32(4)
	inst, del, mis := int32(0), int32(0), int32(0)
	if direction == FORWARD {
		rs := GetRefFlankIdx(marr, mk.Start, direction)
		for i, j := rs, mk.Start; i < int32(len(marr)) && j < mk.End; {
			if marr[i] > 0 || cSeq[i] == readSeq[j] {
				matchNum++
				i++
				j++
				continue
			}
			found := false
			if i < int32(len(marr))-max_gap_len-ExtendSeed &&
				j < mk.End-max_gap_len-ExtendSeed {
				for z := int32(1); z <= max_gap_len; z++ {
					y := 2*(z-1) + 1
					if reflect.DeepEqual(cSeq[i:i+ExtendSeed], readSeq[j+y:j+y+ExtendSeed]) {
						inst += y
						matchNum += ExtendSeed
						i += ExtendSeed
						j += ExtendSeed + y
						found = true
						break
					} else if reflect.DeepEqual(cSeq[i:i+ExtendSeed], readSeq[j+y+1:j+y+1+ExtendSeed]) {
						inst += y + 1
						matchNum += ExtendSeed
						i += ExtendSeed
						j += ExtendSeed + y + 1
						found = true
						break
					} else if reflect.DeepEqual(cSeq[i+z:i+z+ExtendSeed], readSeq[j:j+ExtendSeed]) {
						del += z
						i += z + ExtendSeed
						j += ExtendSeed
						found = true
						break
					} else if reflect.DeepEqual(cSeq[i+z:i+z+ExtendSeed], readSeq[j+z:j+z+ExtendSeed]) {
						mis += z
						i += z + ExtendSeed
						j += z + ExtendSeed
						found = true
						break
					}
				}
			}
			if found == false {
				re, qe := GetNoSeedRegionEnd(marr, int(i), mk.End)
				ma, in, de, mi := DynamicAlign(cSeq[i:re], readSeq[j:qe])
				matchNum += ma
				inst += in
				del += de
				mis += mi
				i = int32(re)
				j = int32(qe)
				// log.Fatalf("[MatchSeq] cSeq:%v\nreadSeq: %v\n", cSeq[i:i+20], readSeq[j:j+20])
			}

		}

	} else { // BACKWARD
		re := GetRefFlankIdx(marr, mk.End-1, direction)
		fmt.Printf("[MatchSeq]i: %d, j: %d\n", re, mk.End-1)
		for i, j := re, mk.End-1; i >= 0 && j >= mk.Start; {
			if marr[i] > 0 || cSeq[i] == readSeq[j] {
				matchNum++
				i--
				j--
				continue
			}
			found := false
			if i > max_gap_len+ExtendSeed && j > max_gap_len+ExtendSeed {
				for z := int32(1); z <= max_gap_len; z++ {
					y := 2*(z-1) + 1
					if reflect.DeepEqual(cSeq[i-ExtendSeed+1:i+1], readSeq[j-y-ExtendSeed+1:j-y+1]) {
						inst += y
						matchNum += ExtendSeed
						i -= ExtendSeed
						j -= (ExtendSeed + y)
						found = true
						break
					} else if reflect.DeepEqual(cSeq[i-ExtendSeed+1:i+1], readSeq[j-y-1-ExtendSeed+1:j-y-1+1]) {
						inst += y + 1
						matchNum += ExtendSeed
						i -= ExtendSeed
						j -= (ExtendSeed + y + 1)
						found = true
						break
					} else if reflect.DeepEqual(cSeq[i-z-ExtendSeed+1:i-z+1], readSeq[j-ExtendSeed+1:j+1]) {
						del += z
						i -= (z + ExtendSeed)
						j -= ExtendSeed
						found = true
						break
					} else if reflect.DeepEqual(cSeq[i-z-ExtendSeed+1:i-z+1], readSeq[j-z-ExtendSeed+1:j-z+1]) {
						mis += z
						i -= (z + ExtendSeed)
						j -= (z + ExtendSeed)
						found = true
						break
					}
				}
			}

			if found == false {
				rs, qs := GetNoSeedRegionStart(marr, int(i), mk.Start)
				ma, in, de, mi := DynamicAlign(cSeq[rs:i+1], readSeq[qs:j+1])
				matchNum += ma
				inst += in
				del += de
				mis += mi
				i = int32(rs) - 1
				j = int32(qs) - 1
				// log.Fatalf("[MatchSeq] cSeq:%v\nreadSeq: %v\n", cSeq[i:i+20], readSeq[j:j+20])
			}

		}

	}
	// fmt.Printf("[MatchSeq]marr: %v\nQ:%v\nR:%v\n", marr[rs:], readSeq[mk.Start:mk.End], cSeq[rs:])

	return matchNum
}

func PrintMatch(readSeq, edgeSeq []byte) {
	rS := make([]byte, len(readSeq))
	eS := make([]byte, len(edgeSeq))
	for i, c := range readSeq {
		rS[i] = bnt.BitNtCharUp[c]
	}
	for i, c := range edgeSeq {
		eS[i] = bnt.BitNtCharUp[c]
	}
	fmt.Fprintf(os.Stdout, ">readSeq\n%s\n", rS)
	fmt.Fprintf(os.Stdout, ">edgeSeq\n%s\n", eS)
}

func ExtendSeedEdge(seedEdgeIndex int, NoContainMKArr []MapEdgeKmers, edgesArr []DBGEdge, nodesArr []DBGNode, mks []MergeKmerInfo, seed int, ReadLen int32, readSeq []byte) (readPath []DBG_MAX_INT) {
	// extend left flank
	{
		smk := NoContainMKArr[seedEdgeIndex]
		readPath = append(readPath, DBG_MAX_INT(smk.EID))
		var nID DBG_MAX_INT
		var qPos int32
		direction := BACKWARD
		if smk.Strand == 0 {
			nID = edgesArr[smk.EID].StartNID
		} else {
			nID = edgesArr[smk.EID].EndNID
		}
		qPos = smk.Start
		for nID > 0 && qPos > 0 {
			eIDArr := GetNextEArr(readPath[len(readPath)-1], nodesArr[nID])
			if len(eIDArr) == 0 {
				break
			}
			fmt.Printf("[ExtendSeedEdge:left] eIDArr: %v\tqPos: %d\n", eIDArr, qPos)
			proArr := FindMEKArr(NoContainMKArr, eIDArr, qPos)
			fmt.Printf("[ExtendSeedEdge:left] proArr: %v\n", proArr)
			if len(proArr) == 0 {
				if len(eIDArr) == 1 {
					readPath = append(readPath, eIDArr[0])
					qPos -= int32(len(edgesArr[eIDArr[0]].Utg.Ks) - Kmerlen)
				} else {
					fmt.Printf("[ExtendSeedEdge] not found alignment!\n")
					break
				}
			} else if len(proArr) == 1 {
				readPath = append(readPath, DBG_MAX_INT(proArr[0].EID))
				qPos = proArr[0].Start
			} else if len(proArr) == 2 {
				if proArr[1].End == proArr[1].End {
					if proArr[0].Start != proArr[1].Start {
						if proArr[0].Start > proArr[1].Start {
							proArr[0], proArr[1] = proArr[1], proArr[0]
						}
						proArr[0].Start = proArr[1].Start
					}
					// need alignment
					var matchBase [2]int32
					matchBase[0] = MatchSeq(readSeq, proArr[0], mks, edgesArr[proArr[0].EID].Utg.Ks, uint32(seed), direction)
					//fmt.Printf("[ExtendSeedEdge] len(readSeq): %d, len(edgeSeq)matchBase[0]: %d\n", matchBase[0])
					PrintMatch(readSeq[proArr[0].Start:proArr[0].End], edgesArr[proArr[0].EID].Utg.Ks)
					matchBase[1] = MatchSeq(readSeq, proArr[1], mks, edgesArr[proArr[1].EID].Utg.Ks, uint32(seed), direction)
					if matchBase[1] > matchBase[0] {
						readPath = append(readPath, DBG_MAX_INT(proArr[1].EID))
						qPos = proArr[1].Start
					} else if matchBase[1] < matchBase[0] {
						readPath = append(readPath, DBG_MAX_INT(proArr[0].EID))
						qPos = proArr[0].Start
					} else {
						fmt.Printf("[ExtendSeedEdge] matchBase: %v\nproArr: %v\n", matchBase, proArr)
						break
					}

				} else {
					fmt.Printf("[ExtendSeedEdge] proArr: %v\n", proArr)
					break
				}
			} else {
				log.Fatalf("[ExtendSeedEdge] not processed case len(proArr) > 2")
			}

			eID := readPath[len(readPath)-1]
			if nID == edgesArr[eID].StartNID {
				nID = edgesArr[eID].EndNID
			} else {
				nID = edgesArr[eID].StartNID
			}

		}
	}

	readPath = ReverseDBG_MAX_INTArr(readPath)

	// extend rigth flank
	{
		smk := NoContainMKArr[seedEdgeIndex]
		var nID DBG_MAX_INT
		var qPos int32
		direction := FORWARD
		if smk.Strand == 0 {
			nID = edgesArr[smk.EID].EndNID
		} else {
			nID = edgesArr[smk.EID].StartNID
		}
		qPos = smk.End
		for nID > 0 && qPos < ReadLen {
			eIDArr := GetNextEArr(readPath[len(readPath)-1], nodesArr[nID])
			if len(eIDArr) == 0 {
				break
			}
			fmt.Printf("[ExtendSeedEdge:rigth] eIDArr: %v\tqPos: %d\n", eIDArr, qPos)
			proArr := FindMEKArr(NoContainMKArr, eIDArr, qPos)
			fmt.Printf("[ExtendSeedEdge:right] proArr: %v\n", proArr)
			if len(proArr) == 0 {
				if len(eIDArr) == 1 {
					readPath = append(readPath, eIDArr[0])
					qPos += int32(len(edgesArr[eIDArr[0]].Utg.Ks) - Kmerlen)
				} else {
					fmt.Printf("[ExtendSeedEdge] not found alignment!\n")
					break
				}
			} else if len(proArr) == 1 {
				readPath = append(readPath, DBG_MAX_INT(proArr[0].EID))
				qPos = proArr[0].End
			} else if len(proArr) == 2 {
				if proArr[0].Start == proArr[1].Start {
					if proArr[1].End != proArr[1].End {
						if proArr[0].End > proArr[1].End {
							proArr[0], proArr[1] = proArr[1], proArr[0]
						}
						proArr[1].End = proArr[0].End
					}
					// need alignment
					var matchBase [2]int32
					matchBase[0] = MatchSeq(readSeq, proArr[0], mks, edgesArr[proArr[0].EID].Utg.Ks, uint32(seed), direction)
					matchBase[1] = MatchSeq(readSeq, proArr[1], mks, edgesArr[proArr[1].EID].Utg.Ks, uint32(seed), direction)
					if matchBase[1] > matchBase[0] {
						readPath = append(readPath, DBG_MAX_INT(proArr[1].EID))
						qPos = proArr[1].Start
					} else if matchBase[1] < matchBase[0] {
						readPath = append(readPath, DBG_MAX_INT(proArr[0].EID))
						qPos = proArr[0].Start
					} else {
						fmt.Printf("[ExtendSeedEdge] matchBase: %v\nproArr: %v\n", matchBase, proArr)
						break
					}

				} else {
					fmt.Printf("[ExtendSeedEdge] proArr: %v\n", proArr)
					break
				}
			} else {
				log.Fatalf("[ExtendSeedEdge] not processed case len(proArr) > 2")
			}

			eID := readPath[len(readPath)-1]
			if nID == edgesArr[eID].StartNID {
				nID = edgesArr[eID].EndNID
			} else {
				nID = edgesArr[eID].StartNID
			}

		}
	}

	return readPath
}

func paraFoundPath(rc chan []MergeKmerInfo, wc chan []DBG_MAX_INT, edgesArr []DBGEdge, nodesArr []DBGNode, seed int32, seqSlice []Seq) {

	minKmerNum := 3
	for {
		mks := <-rc
		fmt.Printf("[paraFoundPath] len(mks): %v\n", len(mks))
		if len(mks) == 0 {
			break
		}
		var mEdgeKArr []MapEdgeKmers
		var ReadLen int32
		for i := 0; i < len(mks); {
			// end := i
			// fmt.Printf("[paraFoundPath] mks[%d]: %v\n", i, mks[i])
			var mek MapEdgeKmers
			j := i + 1
			var ss [2]int
			for ; j < len(mks); j++ {
				if mks[i].RID != mks[j].RID {
					break
				}
				x := (mks[j].Qinfo & 0x1) ^ (mks[j].Rinfo & 0x1)
				ss[x]++
			}
			if j-i < minKmerNum || ss[0] == ss[1] {
				i = j
				continue
			}
			mek.EID = mks[i].RID
			if ss[0] < ss[1] {
				mek.Strand = 1
			} else {
				mek.Strand = 0
			}

			// last := 0 - 2*seed
			lastI := int32(-1)
			// fmt.Printf("[paraFoundPath]i: %d, j:%d\n", i, j)
			for z := i; z < j; z++ {
				x := uint8((mks[z].Qinfo & 0x1) ^ (mks[z].Rinfo & 0x1))
				if x != mek.Strand {
					mks[z].RID = 0
					continue
				}
				if int(mks[z].Rinfo>>1) > len(edgesArr[mks[z].RID].Utg.Ks) {
					log.Fatalf("[paraFoundPath]: len(edge): %v\n", len(edgesArr[mks[z].RID].Utg.Ks))
				}

				idx := int32(mks[z].Qinfo >> 1)
				if lastI < 0 {
					lastI = int32(z)
				} else {
					lastRidx := int32(mks[lastI].Rinfo >> 1)
					lastQidx := int32(mks[lastI].Qinfo >> 1)
					ridx := int32(mks[z].Rinfo >> 1)
					if idx > lastQidx && idx > ReadLen {
						ReadLen = idx + seed
					} else if lastQidx > idx && lastQidx > ReadLen {
						ReadLen = lastQidx + seed
					}
					if (x == 0 && float32((idx-lastQidx)-(ridx-lastRidx)) < ErrorRate*float32(idx-lastQidx)+5) ||
						(x == 1 && float32((idx-lastQidx)-(lastRidx-ridx)) < ErrorRate*float32(idx-lastQidx)+5) {
						if mek.Start == 0 && mek.End == 0 {
							mek.Index = lastI
							mek.Start = lastQidx
							mek.MapBN += seed
						}
						if idx-lastQidx >= seed {
							mek.MapBN += seed
						} else {
							mek.MapBN += (idx - lastQidx)
						}
						// last = idx
						mek.End = idx + seed
						mek.Max = int32(z)
					} else {
						if mek.End == 0 {
							mks[lastI].RID = 0
						}
						mks[z].RID = 0
					}
					lastI = int32(z)
				}
			}
			if mek.MapBN > MIN_KMER_MAP_LEN && mek.End-mek.Start > int32(Kmerlen)*2/3 {
				mEdgeKArr = append(mEdgeKArr, mek)
				// fmt.Printf("[paraFoundPath] mek: %v\n", mek)
			}
			i = j
		}

		// sort the mEdgeKArr
		if len(mEdgeKArr) < 2 {
			log.Fatalf("[paraFoundPath] mEdgeKArr: %v\n", mEdgeKArr)
		}
		sort.Sort(sort.Reverse(MEKArr(mEdgeKArr)))

		// fmt.Printf("[paraFoundPath] mEdgeKArr: %v\nReadLen: %d\n", mEdgeKArr, ReadLen)
		fmt.Printf("[paraFoundPath]ReadLen: %d\n", ReadLen)
		zoom := int32(100)
		depth := make([]uint8, (ReadLen+zoom-1)/zoom)
		var NoContainMKArr []MapEdgeKmers
		for i := 0; i < len(mEdgeKArr); i++ {
			unique := true
			for j := 0; j < len(NoContainMKArr); j++ {
				if mEdgeKArr[i].Start == NoContainMKArr[j].Start && mEdgeKArr[i].End == NoContainMKArr[j].End {
					if mEdgeKArr[i].MapBN == NoContainMKArr[j].MapBN {
						break
					} else {
						unique = false
						break
					}

				}
				if mEdgeKArr[i].Start >= NoContainMKArr[j].Start &&
					mEdgeKArr[i].End <= NoContainMKArr[j].End {
					mapBaseN := SummaryRegionMapBaseNum(NoContainMKArr[j], int(mEdgeKArr[i].Start), int(mEdgeKArr[i].End), int(seed), mks)
					// if int32(mapBaseN) < mEdgeKArr[i].MapBN {
					// 	log.Fatalf("[paraFoundPath] mEdgeKArr[%d]: %v\n NoContainMKArr[%d]: %v, mapBaseN: %v\n", i, mEdgeKArr[i], j, NoContainMKArr[j], mapBaseN)
					if int32(mapBaseN) >= mEdgeKArr[i].MapBN {
						unique = false
						break
					}
				}
			}
			if unique {
				NoContainMKArr = append(NoContainMKArr, mEdgeKArr[i])
				for j := mEdgeKArr[i].Start; j < mEdgeKArr[i].End; j += zoom {
					if int(j/zoom) >= len(depth) {
						log.Fatalf("[paraFoundPath] j/zoom: %d\tlen(depth): %v\n", j/zoom, len(depth))
					}
					depth[j/zoom] += 1
				}
			}
		}
		fmt.Printf("[paraFoundPath] NoContainMKArr: %v\ndepth: %v\n", NoContainMKArr, depth)

		// found read path
		// found max unique region
		// min_union_num := 2
		start, end := FoundMaxUniqueRegion(depth)
		fmt.Printf("[paraFoundPath]start, end = %d, %d\n", start, end)
		seedEdgeIndex, overlapLen := GetSeedEdge(NoContainMKArr, int32(start)*zoom, int32(end)*zoom)
		if int32(overlapLen) < 2*zoom {
			log.Fatalf("[paraFoundPath] ovlelapLen: %v\n", overlapLen)
		}
		fmt.Printf("[paraFoundPath]SeedMK : %v\n", NoContainMKArr[seedEdgeIndex])
		readPath := ExtendSeedEdge(seedEdgeIndex, NoContainMKArr, edgesArr, nodesArr, mks, int(seed), ReadLen, seqSlice[mks[0].QID].S)
		wc <- readPath

	}

	var tmp []DBG_MAX_INT
	wc <- tmp
}

func WriteReadPath(wc chan []DBG_MAX_INT, numCPU int) {
	terminalNum := 0
	for {
		rp := <-wc
		if len(rp) == 0 {
			terminalNum++
			if terminalNum == numCPU {
				break
			}
		}

		fmt.Printf("[WriteReadPath]ReadPath: %v\n", rp)
	}
}

func FindPath(mkiArr []MergeKmerInfo, numCPU, seed int, edgesArr []DBGEdge, nodesArr []DBGNode, seqSlice []Seq) {
	rc := make(chan []MergeKmerInfo, numCPU)
	wc := make(chan []DBG_MAX_INT, numCPU)

	fmt.Printf("[FindPath]len(mkiArr): %v\n", len(mkiArr))
	go GetMki(mkiArr, rc, numCPU)

	for i := 0; i < numCPU; i++ {
		go paraFoundPath(rc, wc, edgesArr, nodesArr, int32(seed), seqSlice)
	}

	WriteReadPath(wc, numCPU)

}

func AlignQuery(qrc chan []Seq, refSlice []SeedKmerInfo, seed, width, numCPU int, edgesArr []DBGEdge, nodesArr []DBGNode) {
	for {
		seqSlice := <-qrc
		if len(seqSlice) == 0 {
			break
		}
		rc := make(chan Seq, numCPU)
		wc := make(chan SeedKmerInfo, numCPU)
		go GetSeq(seqSlice, rc, numCPU)

		for i := 0; i < numCPU; i++ {
			go paraCollectMinKmer(rc, wc, seed, width)
		}
		MinQuery := CollectSeedKmerInfo(wc, numCPU)
		fmt.Printf("[AlignQuery]len(MinQuery): %d\n", len(MinQuery))
		// seqSlice = nil
		// radix sort
		RadixSort(MinQuery, seed*bnt.NumBitsInBase, numCPU)

		// merge ref and query MinQuery
		mwc := make(chan MergeKmerInfo, numCPU*2)
		part := (len(MinQuery)-1)/numCPU + 1
		for i := 0; i < numCPU; i++ {
			if i == numCPU-1 {
				go paraMergeRefAndQuy(refSlice, MinQuery, i*part, len(MinQuery), mwc)
			} else {
				go paraMergeRefAndQuy(refSlice, MinQuery, i*part, (i+1)*part, mwc)
			}
		}

		mergeKInfoArr := CollectMergeKmerInfo(mwc, numCPU)
		MinQuery = nil
		fmt.Printf("[AlignQuery] len(mergeKInfoArr): %d\n", len(mergeKInfoArr))

		// radix sort for mergeedKmerInfo
		MergeRadixSort(mergeKInfoArr, numCPU)
		fmt.Printf("[AlignQuery]QID\tRID\tQStart\tQstrand\tRStart\tRstrand\n")
		for i := 0; i < len(mergeKInfoArr); i++ {
			mk := mergeKInfoArr[i]
			fmt.Printf("[AlignQuery] %d\t%d\t%d\t%d\t%d\t%d\n", mk.QID, mk.RID, mk.Qinfo>>1, mk.Qinfo&0x1, mk.Rinfo>>1, mk.Rinfo&0x1)
		}

		// found most probality path in the edges
		FindPath(mergeKInfoArr, numCPU, seed, edgesArr, nodesArr, seqSlice)

	}
}

func MapDBG(c cli.Command) {
	fmt.Println(c.Flags(), c.Parent().Flags())

	//argsCheck(c)
	fnName := c.Parent().Flag("C").String()
	cfgInfo, err := constructcf.ParseCfg(string(fnName))
	if err != nil {
		log.Fatal("[parseCfg] found parseCfg err")
	}
	fmt.Println(cfgInfo)

	if Kmerlen, err = strconv.Atoi(c.Parent().Flag("K").String()); err != nil {
		log.Fatal("flag 'K' set error")
	}
	prefix := c.Parent().Flag("p").String()
	numCPU, err := strconv.Atoi(c.Parent().Flag("t").String())
	if err != nil {
		log.Fatal("'t' set error")
	}
	seed, err := strconv.Atoi(c.Flag("Seed").String())
	if err != nil {
		log.Fatal("Seed length set error")
	}
	width, err := strconv.Atoi(c.Flag("Width").String())
	if err != nil {
		log.Fatal("window of bandwidth set error")
	}

	runtime.GOMAXPROCS(numCPU)

	// read nodes file and transform to array mode, for more quickly access
	smfyNodesfn := prefix + ".nodes.smfy.mmap"
	nodeMap := NodeMapMmapReader(smfyNodesfn)
	DBGStatfn := prefix + ".DBG.stat"
	nodesSize, edgesSize := DBGStatReader(DBGStatfn)
	nodesArr := make([]DBGNode, nodesSize)
	NodeMap2NodeArr(nodeMap, nodesArr)

	// Restore edges info
	//edgesStatfn := prefix + ".edges.stat"
	//edgesSize := EdgesStatReader(edgesStatfn)
	edgesArr := make([]DBGEdge, edgesSize)
	edgesfn := prefix + ".edges.smfy.fq"
	LoadEdgesfqFromFn(edgesfn, edgesArr, false)

	//construct target sequence index
	rc := make(chan Seq, numCPU)
	wc := make(chan SeedKmerInfo, numCPU)
	go GetEdge(edgesArr, rc, numCPU)

	for i := 0; i < numCPU; i++ {
		go paraCollectMinKmer(rc, wc, seed, width)
	}
	skSlice := CollectSeedKmerInfo(wc, numCPU)

	// radix sort
	{
		fmt.Printf("[MapDBG] len(skSlice): %v\n", len(skSlice))
		RadixSort(skSlice, seed*bnt.NumBitsInBase, numCPU)
		// fmt.Printf("[MapDBG] len(skSlice): %v\n", len(skSlice))
		// for _, s := range skSlice {
		// 	fmt.Printf("[MapDBG]%v\n", s)
		// }
	}

	// Get min kmer of query
	blockSize := len(skSlice) * width
	qrc := make(chan []Seq, 1)
	go GetRawReads(qrc, cfgInfo, blockSize)
	AlignQuery(qrc, skSlice, seed, width, numCPU, edgesArr, nodesArr)
}
