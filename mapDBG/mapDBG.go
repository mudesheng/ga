package mapDBG

import (
	// "bufio"
	// "compress/gzip"
	// "encoding/binary"
	"fmt"
	"ga/bnt"
	"ga/constructcf"
	"ga/constructdbg"
	// "ga/cuckoofilter"
	"github.com/biogo/biogo/alphabet"
	"github.com/biogo/biogo/io/seqio/fasta"
	"github.com/biogo/biogo/seq/linear"
	"github.com/jwaldrip/odin/cli"
	"io"
	"log"
	// "math"
	"os"
	"runtime"
	"sort"
	"strconv"
	// "strings"
	// "time"
)

const (
	AllState    = 1
	ScaffState  = 2
	GapState    = 3
	EndString   = "end"
	ReadSeqSize = 100
)

var Kmerlen = -1

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

func GetEdge(edgesArr []constructdbg.DBGEdge, rc chan Seq, numCPU int) {
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
	var slc []MergeKmerInfo
	slc = append(slc, mkiArr[0])
	for _, m := range mkiArr[1:] {
		if m.QID != slc[len(slc)-1].QID {
			if len(slc) > 10 {
				rc <- slc
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
	EID               uint32
	MapBN             int32 // the total mapping base number
	Strand            uint8 // the Ref strand ^ Query strand
	Index, Start, End int32
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

func paraFoundPath(rc chan []MergeKmerInfo, edgesArr []constructdbg.DBGEdge, nodesArr []constructdbg.DBGNode, seed int32, seqSlice []Seq) {

	minKmerNum := 3
	for {
		mks := <-rc
		if len(mks) == 0 {
			break
		}
		var mEdgeKArr []MapEdgeKmers
		for i := 0; i < len(mks); {
			end := i
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
			if j-i < 3 || ss[0] == ss[1] {
				i = j
				continue
			}
			mek.EID = mks[i].RID
			if ss[0] < ss[1] {
				mek.Strand = 1
			} else {
				mek.Strand = 0
			}
			// clear MergeKmerInfo by reference start point
			// forOB z := i; z < j-1; z++ {

			// }

			last := 0 - 2*seed
			lastI := -1
			for z := i; z < j; z++ {
				x := uint8((mks[z].Qinfo & 0x1) ^ (mks[z].Rinfo & 0x1))
				if x != mek.Strand {
					continue
				}

				idx := int32(mks[z].Qinfo >> 1)
				if lastI == -1 && z+1 < j {
					x2 := uint8((mks[z+1].Qinfo & 0x1) ^ (mks[z+1].Rinfo & 0x1))
					if x == x2 {
						idx2 := int32(mks[z+1].Qinfo >> 1)
						ridx := int32(mks[z].Rinfo >> 1)
						ridx2 := int32(mks[z+1].Rinfo >> 1)
						if x == 0 {
							if (idx2-idx)-(ridx2-ridx) < ErrorRate*(idx2-idx)+5 {
								mek.Index = int32(z)
								mek.Start = idx
								lastI = idx
							}
						} else {
							if (idx2-idx)-(ridx-ridx2) < Errorr*(idx2-idx)+5 {
								mek.Index = int32(z)
								mek.Start = idx
								lastI = idx
							}
						}
					}
					mek.Index = int32(z)
					mek.Start = idx
				}
				if idx-last >= seed {
					mek.MapBN += seed
				} else {
					mek.MapBN += (idx - last)
				}
				last = idx
				mek.End = idx + seed
			}
			mEdgeKArr = append(mEdgeKArr, mek)
			i = j
		}

		// sort the mEdgeKArr
		sort.Sort(sort.Reverse(MEKArr(mEdgeKArr)))

	}
}

func FindPath(mkiArr []MergeKmerInfo, numCPU, seed int, edgesArr []constructdbg.DBGEdge, nodesArr []constructdbg.DBGNode, seqSlice []Seq) {
	rc := make(chan []MergeKmerInfo, numCPU)

	go GetMki(mkiArr, rc, numCPU)

	for i := 0; i < numCPU; i++ {
		go paraFoundPath(rc, edgesArr, nodesArr, seed, seqSlice)
	}

}

func AlignQuery(qrc chan []Seq, refSlice []SeedKmerInfo, seed, width, numCPU int, edgesArr []constructdbg.DBGEdge, nodesArr []constructdbg.DBGNode) {
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
	nodeMap := constructdbg.NodeMapMmapReader(smfyNodesfn)
	nodesStatfn := prefix + ".nodes.stat"
	nodesSize := constructdbg.NodesStatReader(nodesStatfn)
	nodesArr := make([]constructdbg.DBGNode, nodesSize)
	constructdbg.NodeMap2NodeArr(nodeMap, nodesArr)

	// Restore edges info
	edgesStatfn := prefix + ".edges.stat"
	edgesSize := constructdbg.EdgesStatReader(edgesStatfn)
	edgesArr := make([]constructdbg.DBGEdge, edgesSize)
	edgesfn := prefix + ".edges.smfy.fq"
	constructdbg.LoadEdgesfqFromFn(edgesfn, edgesArr, false)

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
