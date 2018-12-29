package deconstructdbg

import (
	"bufio"
	"fmt"
	"io"
	"log"
	"math"
	"os"
	"reflect"
	"sort"
	"strconv"
	"strings"

	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/sam"
	"github.com/mudesheng/ga/bnt"
	"github.com/mudesheng/ga/cbrotli"
	"github.com/mudesheng/ga/constructcf"
	"github.com/mudesheng/ga/constructdbg"
)

const (
	OutComing = true
	InComing  = false
)

func GetSamRecord(bamfn string, rc chan []sam.Record, numCPU int) {
	fp, err := os.Open(bamfn)
	if err != nil {
		log.Fatalf("[GetSamRecord] open file: %s failed, err: %v\n", bamfn, err)
	}
	bamfp, err := bam.NewReader(fp, numCPU/5+1)
	if err != nil {
		log.Fatalf("[GetSamRecord] create bam.NewReader err: %v\n", err)
	}
	defer bamfp.Close()
	defer fp.Close()
	var rArr []sam.Record
	// var cigar sam.Cigar
	//var NM = sam.Tag{'N', 'M'}
	//var AS = sam.Tag{'A', 'S'}
	for {
		r, err := bamfp.Read()
		if err != nil {
			break
		}
		if r.Flags&sam.Unmapped != 0 {
			continue
		}
		/*v := r.AuxFields.Get(NM).Value()
		vi := 0
		switch v.(type) {
		case uint8:
			vi = int(v.(uint8))
		case uint16:
			vi = int(v.(uint16))
		}
		if vi != 0 {
			continue
		}
		as := r.AuxFields.Get(AS).Value()
		asi := 0
		switch as.(type) {
		case uint8:
			asi = int(as.(uint8))
		case uint16:
			asi = int(as.(uint16))
		case uint32:
			asi = int(as.(uint32))
		default:
			log.Fatalf("[GetSamRecord] as type unknown\n")
		}*/
		// fmt.Printf("[GetSamRecord]AS:%d, cigar: %v\n", asi, r.Cigar.String())
		/*if asi < Kmerlen {
			continue
		}*/
		// Debug
		// if r.Cigar.String() == "400M" {
		// 	continue
		// }
		//Debug
		if len(rArr) > 0 && rArr[0].RefID() != r.RefID() {
			rc <- rArr
			rArr = nil
		}
		rArr = append(rArr, *r)
	}
	if len(rArr) > 0 {
		rc <- rArr
	}
	// send terminal signals
	for i := 0; i < numCPU; i++ {
		var nilArr []sam.Record
		rc <- nilArr
	}
}

type LRRecord struct {
	Ins, Del, Mis, Mch int // Insertion, Deletion and Mismatch base number
	RefID, ID          constructdbg.DBG_MAX_INT
	//Pair                        uint8 // 1 note Ref Name RefID/1, 2 note RefID/2
	RefStart, RefEnd, RefLen int
	Start, End, Len          int
	//RefSeq, Seq                 string
	MapNum, GapMapNum int  // same as minimap2 column 10 and 11, number mapping base, add gap mapping length
	Strand            bool // = constructdbg.PLUS or constructdbg.MINUS
	ReadSeqBnt        []byte
}

type MapInfo struct {
	ReadPos int // note read sequence mapping start position
	EdgePos int // note edge sequence mapping start position
	EIDIdx  int // mapping edge ID in the extArr[0]
	//NID     constructdbg.DBG_MAX_INT // note edge to next node ID
	Strand bool // note PLUS or MINUS strand mapping to the readSeq
}

func (R LRRecord) GetNM() int {
	return R.Mis + R.Ins + R.Del
}

func (R LRRecord) GetRefCon() int {
	return R.RefEnd - R.RefStart
}

func AccumulateCigar(cigar sam.Cigar, strand bool) (Mnum, Inum, Dnum, Pos int) {
	for i, co := range cigar {
		switch co.Type() {
		case sam.CigarMatch:
			Mnum += co.Len()
		case sam.CigarDeletion:
			Dnum += co.Len()
		case sam.CigarInsertion:
			Inum += co.Len()
		case sam.CigarHardClipped:
			if strand && i == len(cigar)-1 {
				Pos = co.Len()
			} else if !strand && i == 0 {
				Pos = co.Len()
			}
		case sam.CigarSoftClipped:
			if strand && i == len(cigar)-1 {
				Pos = co.Len()
			} else if !strand && i == 0 {
				Pos = co.Len()
			}
		}
	}

	return
}

func GetAuxUint(v interface{}) int {
	var nmv int
	switch v.(type) {
	case uint8:
		nmv = int(v.(uint8))
	case uint16:
		nmv = int(v.(uint16))
	case uint32:
		nmv = int(v.(uint32))
	default:
		log.Fatalf("[GetAuxNM] unknown type of nm: %v\n", v)
	}
	return nmv
}

type LRRecordArr []LRRecord

func (a LRRecordArr) Len() int {
	return len(a)
}

func (a LRRecordArr) Less(i, j int) bool {
	return a[i].Start < a[j].Start
}

func (a LRRecordArr) Swap(i, j int) {
	a[i], a[j] = a[j], a[i]
}

type Diff struct {
	Idx        int
	Unidentity float64
}
type DiffArr []Diff

func (a DiffArr) Len() int {
	return len(a)
}

func (a DiffArr) Less(i, j int) bool {
	return a[i].Unidentity < a[j].Unidentity
}

func (a DiffArr) Swap(i, j int) {
	a[i], a[j] = a[j], a[i]
}

type LRRecordScoreArr []LRRecord

func (a LRRecordScoreArr) Len() int {
	return len(a)
}

func (a LRRecordScoreArr) Less(i, j int) bool {
	return a[i].MapNum > a[j].MapNum
}

func (a LRRecordScoreArr) Swap(i, j int) {
	a[i], a[j] = a[j], a[i]
}

/*func GetIDArr(arr []LRRecord) (IDArr []constructdbg.DBG_MAX_INT) {
	for _, R := range arr {
		if constructdbg.IsInDBG_MAX_INTArr(refIDArr, R.RefID) == false {
			refIDArr = append(refIDArr, R.RefID)
		}
	}
	return refIDArr
} */

func PairIndex(arr []LRRecord, idx int, item LRRecord) int {
	nidx := -1
	for i := idx; i < len(arr); i++ {
		if arr[i].RefID == item.RefID {
			nidx = i
			break
		}
	}

	return nidx
}

func findSeedEID(arr []LRRecord, edgesArr []constructdbg.DBGEdge, flankAllow int) int {
	pos := -1
	var maxMapNum int
	for i, item := range arr {
		e := edgesArr[item.RefID]
		if e.GetUniqueFlag() > 0 && maxMapNum < item.MapNum {
			if item.Start < flankAllow {
				if item.Strand == constructdbg.PLUS {
					if item.RefLen-item.RefEnd > flankAllow {
						continue
					}
				} else {
					if item.RefStart > flankAllow {
						continue
					}
				}
			} else if item.Len-item.End < flankAllow {
				if item.Strand == constructdbg.PLUS {
					if item.RefStart > flankAllow {
						continue
					}
				} else {
					if item.RefLen-item.RefEnd > flankAllow {
						continue
					}
				}
			} else {
				continue
			}
			pos = i
			maxMapNum = item.MapNum
		}
	}
	// check read only mapping one edge
	if pos >= 0 && arr[pos].Start < flankAllow && arr[pos].Len-arr[pos].End < flankAllow {
		pos = -1
	}

	return pos
}

/*func findMapStrand(RefID constructdbg.DBG_MAX_INT, pos int, arr []LRRecord) (strand bool) {
	var max int
	for i := pos; i < len(arr); i++ {
		if arr[i].RefID == RefID {
			if max == 0 {
				max = arr[i].Mch
				strand = arr[i].Strand
			} else {
				if max < arr[i].Mch {
					max = arr[i].Mch
					strand = arr[i].Strand
				}
			}
		} else {
			if max > 0 {
				break
			}
		}
	}

	if max == 0 {
		log.Fatalf("[findMapStrand] not found RefID: %v in the arr[%v:]\n\tarr: %v\n", RefID, pos, arr)
	}
	return strand
}*/

func ReverseLRRcordArr(arr []LRRecord) []LRRecord {
	for i := 0; i < len(arr)/2; i++ {
		arr[i], arr[len(arr)-1-i] = arr[len(arr)-1-i], arr[i]
	}

	return arr
}

func GetExtendPath(IDArr []constructdbg.DBG_MAX_INT, eID constructdbg.DBG_MAX_INT, direction uint8) (arr []constructdbg.DBG_MAX_INT) {
	idx := constructdbg.IndexEID(IDArr, eID)
	if idx < 0 {
		log.Fatalf("[GetExtendPath] not found eID: %v in the array: %v\n", eID, IDArr)
	}

	if direction == constructdbg.FORWARD {
		if idx < len(IDArr)-1 {
			arr = make([]constructdbg.DBG_MAX_INT, len(IDArr[idx+1:]))
			copy(arr, IDArr[idx+1:])
		}
	} else {
		if idx > 0 {
			arr = make([]constructdbg.DBG_MAX_INT, len(IDArr[:idx]))
			copy(arr, IDArr[:idx])
			constructdbg.ReverseDBG_MAX_INTArr(arr)
		}
	}

	return arr
}

func SearchLRRecordArr(arr []LRRecord, aidx int, direction uint8, ea []constructdbg.DBG_MAX_INT, readPos, searchLen int) (isa []int) {
	count := 0
	if direction == constructdbg.FORWARD {
		// search before partition
		for i := aidx - 1; i >= 0; i-- {
			if arr[i].Start >= readPos {
				for _, eID := range ea {
					if eID == arr[i].RefID {
						isa = append(isa, i)
						count++
						break
					}
				}
			} else {
				break
			}
		}
		max := readPos + searchLen
		for i := aidx + 1; i < len(arr); i++ {
			if arr[i].Start > max {
				break
			}
			for _, eID := range ea {
				if eID == arr[i].RefID {
					isa = append(isa, i)
					count++
					break
				}
			}
			if count == len(ea) {
				break
			}
		}
	} else {
		// search after partition
		/*for i := aidx + 1; i < len(arr); i++ {
			if arr[i].End <= readPos {
				for _, eID := range ea {
					if eID == arr[i].RefID {
						isa = append(isa, i)
						count++
						break
					}
				}
			} else {
				break
			}
		}*/
		min := readPos - searchLen
		for i := aidx - 1; i >= 0; i-- {
			if arr[i].End < min {
				continue
			}
			for _, eID := range ea {
				if eID == arr[i].RefID {
					isa = append(isa, i)
					count++
					break
				}
			}
			if count == len(ea) {
				break
			}
		}
	}

	return
}

func GetMinLenEID(edgesArr []constructdbg.DBGEdge, eIDArr []constructdbg.DBG_MAX_INT) int {
	if len(eIDArr) <= 1 {
		log.Fatalf("[GetMinLenEID] len(eIDArr): %v <= 1 in the array: %v\n", len(eIDArr), eIDArr)
	}
	minLen := edgesArr[eIDArr[0]].GetSeqLen()
	for _, eID := range eIDArr[1:] {
		if edgesArr[eID].GetSeqLen() < minLen {
			minLen = edgesArr[eID].GetSeqLen()
		}
	}
	return minLen
}

// score bigger than  suboptimal edge
/*func GetMaxMchIdx(arr []LRRecord) (idx, subidx int, diffScore int) {
	diffScore = -1
	idx, subidx = 0, -1
	score := arr[0].Mch - arr[0].Mis - arr[0].Ins - arr[0].Del
	subscore := -1
	for i := 1; i < len(arr); i++ {
		R := arr[i]
		sc := R.Mch - R.Mis - R.Ins - R.Del
		if sc > score {
			subidx = idx

			idx = 1
			score = sc
			diffScore = sc - score
		} else if sc < score {
			subidx
			diffScore = score - sc
		}
	}

	return idx, diffScore
}*/

func FindShareNID(eID1, eID2 constructdbg.DBG_MAX_INT, edgesArr []constructdbg.DBGEdge) constructdbg.DBG_MAX_INT {
	e1, e2 := edgesArr[eID1], edgesArr[eID2]
	if e1.StartNID == e2.StartNID || e1.StartNID == e2.EndNID {
		return e1.StartNID
	} else if e1.EndNID == e2.StartNID || e1.EndNID == e2.EndNID {
		return e1.EndNID
	}
	log.Fatalf("[FindShareNID] not found share node ID, e1: %v\n\te2: %v\n", e1, e2)
	return 0
}

func GetFlankLen(R LRRecord, diretion uint8) int {
	if diretion == constructdbg.FORWARD {
		if R.Strand == constructdbg.MINUS {
			return R.RefStart
		} else {
			return R.RefLen - R.RefEnd
		}
	} else { // BACKWARD
		if R.Strand == constructdbg.MINUS {
			return R.RefLen - R.RefEnd
		} else {
			return R.RefStart
		}
	}
}

/*func RelocateLRRecordArr(arr []LRRecord, si int, direction uint8, newArr []constructdbg.DBG_MAX_INT, searchLen int, edgesArr []constructdbg.DBGEdge, opt Options) (int, int) {
	nl := searchLen
	for _, eID := range newArr {
		nl += edgesArr[eID].GetSeqLen() - (opt.Kmer - 1)
		idx, ok := SearchLRRecordArr(arr, si, direction, eID, nl*6/5)
		if ok && arr[idx].RefID == eID {
			si = idx
			nl = GetFlankLen(arr[idx], direction)
		}
	}

	return si, nl
}*/

type ExtPathInfo struct {
	Path        []constructdbg.DBG_MAX_INT
	LRRecordArr []LRRecord
	ExtLen      int
	Score       int
	Nd          constructdbg.DBGNode
}

type ExtPathInfoArr []ExtPathInfo

func (a ExtPathInfoArr) Len() int {
	return len(a)
}

func (a ExtPathInfoArr) Less(i, j int) bool {
	return a[i].Score > a[j].Score
}

func (a ExtPathInfoArr) Swap(i, j int) {
	a[i], a[j] = a[j], a[i]
}

/*func GetNM(R LRRecord) (mch, mis, ins, del int) {
	for i, t := range R.Seq {
		c := byte(t)
		if c == R.RefSeq[i] {
			mch++
		} else if c == '-' {
			del++
		} else if R.RefSeq[i] == '-' {
			ins++
		} else {
			mis++
		}
	}

	return
}*/

/*func Getscore(R LRRecord, begin, end int) (score int) {
	if begin == R.Start && end == R.Start+R.ConLen {
		score = R.Mch - R.Mis - R.Ins - R.Del
		return score
	}
	ins, del, mis, mch := 0, 0, 0, 0
	pos := 0
	for i, t := range R.Seq {
		c := byte(t)
		if begin <= pos && pos < end {
			if c == R.RefSeq[i] {
				mch++
			} else if c == '-' {
				del++
			} else if R.RefSeq[i] == '-' {
				ins++
			} else {
				mis++
			}
		}
		if c != '-' {
			pos++
		}
	}

	score = mch - mis - ins - del
	fmt.Printf("[Getscore] LongRead Start: %v, End : %v, mch: %v, ins: %v, del: %v, mis: %v, score: %v\n", begin, end, mch, ins, del, mis, score)
	return score
}*/
func priorIndex(arr []LRRecord, idx int) int {
	pi := -1
	for i := idx - 1; i >= 0; i-- {
		if arr[i].RefID != 0 {
			pi = i
			break
		}
	}

	return pi
}

func nextIndex(arr []LRRecord, idx int) int {
	ni := -1
	for i := idx + 1; i < len(arr); i++ {
		if arr[i].RefID != 0 {
			ni = i
			break
		}
	}

	return ni
}

/*// region special by Refence(Long Reads)
func GetScoreFixedLen(edgesArr []constructdbg.DBGEdge, epi ExtPathInfo, flankPos int, direction uint8) int {
	var score int
	if direction == constructdbg.FORWARD {
		for i, _ := range epi.Path {
			if epi.LRRecordArr[i].RefID != 0 {
				begin, end := epi.LRRecordArr[i].Start, epi.LRRecordArr[i].Start+epi.LRRecordArr[i].ConLen
				pi := priorIndex(epi.LRRecordArr, i)
				ni := nextIndex(epi.LRRecordArr, i)
				if pi >= 0 && epi.LRRecordArr[pi].Start+epi.LRRecordArr[pi].ConLen > begin {
					begin += (epi.LRRecordArr[pi].Start + epi.LRRecordArr[pi].ConLen - begin) / 2
				}
				if ni > i && epi.LRRecordArr[ni].Start < end {
					end -= (end - epi.LRRecordArr[ni].Start) / 2
				}
				if begin > flankPos {
					break
				}
				if end > flankPos {
					end = flankPos
				}
				score += Getscore(epi.LRRecordArr[i], begin, end)
			}
		}
	} else { // BACKWARD
		for i, _ := range epi.Path {
			if epi.LRRecordArr[i].RefID != 0 {
				begin, end := epi.LRRecordArr[i].Start, epi.LRRecordArr[i].Start+epi.LRRecordArr[i].ConLen
				pi := priorIndex(epi.LRRecordArr, i)
				ni := nextIndex(epi.LRRecordArr, i)
				if pi >= 0 && epi.LRRecordArr[pi].Start < end {
					end -= (end - epi.LRRecordArr[pi].Start) / 2
				}
				if ni > i && epi.LRRecordArr[ni].Start+epi.LRRecordArr[ni].ConLen > begin {
					begin += (epi.LRRecordArr[ni].Start + epi.LRRecordArr[ni].ConLen - begin) / 2
				}
				if end < flankPos {
					break
				}
				if begin < flankPos {
					begin = flankPos
				}
				score += Getscore(epi.LRRecordArr[i], begin, end)
			}
		}
	}

	return score
}

func IsInRefRegionContain(a1, a2 LRRecord) (constructdbg.DBG_MAX_INT, bool) {
	a1E := a1.Start + a1.ConLen
	a2E := a2.Start + a2.ConLen
	if (a1.Start <= a2.Start && a2E <= a1E) || (a2.Start < a1.Start && a1E < a2E) {
		s1, s2 := a1.Score, a2.Score
		if s1 > s2 {
			return a1.RefID, true
		} else {
			return a2.RefID, true
		}
	}
	return 0, false
}

func extendPathDistinguish(edgesArr []constructdbg.DBGEdge, nd constructdbg.DBGNode, nodesArr []constructdbg.DBGNode, eIDArr []constructdbg.DBG_MAX_INT, lrArr []LRRecord, si, searchLen int, direction uint8, opt Options) (extArr []constructdbg.DBG_MAX_INT) {
	var epiArr []ExtPathInfo
	st := list.New()
	for _, eID := range eIDArr {
		var epi ExtPathInfo
		epi.Path = append(epi.Path, eID)
		epi.ExtLen = edgesArr[eID].GetSeqLen()
		if edgesArr[eID].StartNID == nd.ID {
			epi.Nd = nodesArr[edgesArr[eID].EndNID]
		} else {
			epi.Nd = nodesArr[edgesArr[eID].StartNID]
		}
		if epi.ExtLen > opt.ExtLen {
			epiArr = append(epiArr, epi)
		} else {
			st.PushBack(epi)
		}
	}
	for st.Len() > 0 {
		t := st.Back()
		st.Remove(t)
		epi := t.Value.(ExtPathInfo)
		ea := constructdbg.GetNearEdgeIDArr(epi.Nd, epi.Path[len(epi.Path)-1])
		if len(ea) == 0 {
			epiArr = append(epiArr, epi)
			continue
		}
		for _, id := range ea {
			var nepi ExtPathInfo
			nepi.Path = append(nepi.Path, epi.Path...)
			nepi.Path = append(nepi.Path, id)
			nepi.ExtLen = epi.ExtLen + edgesArr[id].GetSeqLen() - (opt.Kmer - 1)
			if edgesArr[id].StartNID == epi.Nd.ID {
				nepi.Nd = nodesArr[edgesArr[id].EndNID]
			} else {
				nepi.Nd = nodesArr[edgesArr[id].StartNID]
			}
			if nepi.ExtLen > opt.ExtLen {
				epiArr = append(epiArr, nepi)
			} else {
				st.PushBack(nepi)
			}
		}
	}
	//fmt.Printf("[extendPathDistinguish] epiArr: %v\n", epiArr)
	// compute total score each epi in epiArr
	for i, epi := range epiArr {
		lrArrIdx := si
		sl := searchLen
		epi.LRRecordArr = make([]LRRecord, len(epi.Path))
		//var mapArr []LRRecord
		fmt.Printf("[extendPathDistinguish] epi: %v\n", epi)
		for j, eID := range epi.Path {
			sl += edgesArr[eID].GetSeqLen() - (opt.Kmer - 1)
			idx, ok := SearchLRRecordArr(lrArr, lrArrIdx, direction, eID, sl*6/5)
			if ok && lrArr[idx].RefID == eID {
				fmt.Printf("[extendPathDistinguish] lrArr[%v]: %v\n", idx, lrArr[idx])
				lrArrIdx = idx
				sl = GetFlankLen(lrArr[idx], direction)
				if j > 0 && epi.LRRecordArr[j-1].RefID != 0 { // if has a contain edge mapped, ignore...
					if id, ok := IsInRefRegionContain(epi.LRRecordArr[j-1], lrArr[idx]); ok {
						if id == lrArr[idx].RefID {
							epi.LRRecordArr[j-1].RefID = 0
							epi.LRRecordArr[j] = lrArr[idx]
						}
					} else {
						epi.LRRecordArr[j] = lrArr[idx]
					}
				} else {
					epi.LRRecordArr[j] = lrArr[idx]
				}
			}
		}
		var flankPos int
		if direction == constructdbg.FORWARD {
			flankPos = lrArr[si].Start + lrArr[si].ConLen + opt.ExtLen
		} else {
			flankPos = lrArr[si].Start - opt.ExtLen
		}
		epi.Score = GetScoreFixedLen(edgesArr, epi, flankPos, direction) // region special by Refence(Long Reads)
		epiArr[i] = epi
	}

	sort.Sort(ExtPathInfoArr(epiArr))
	for i, item := range epiArr {
		fmt.Printf("[extendPathDistinguish] epiArr[%v]: %v\n", i, item)
	}
	if epiArr[0].Score > epiArr[1].Score {
		extArr = append(extArr, epiArr[0].Path...)
	}

	return extArr
}*/

func IndexArrPos(subA []int, arr []LRRecord, eID constructdbg.DBG_MAX_INT) int {
	for _, idx := range subA {
		if arr[idx].RefID == eID {
			return idx
		}
	}
	return -1
}

func GetConsisPath(pathMat []constructdbg.Path, eID, neID constructdbg.DBG_MAX_INT) (eIDArr []constructdbg.DBG_MAX_INT) {
	var pm [][]constructdbg.DBG_MAX_INT
	maxLen := 0
	for _, p := range pathMat {
		j := constructdbg.IndexEID(p.IDArr, eID)
		if j < 0 {
			log.Fatalf("[GetConsisPath] eID: %v not in p: %v\n", eID, p)
		}
		var na []constructdbg.DBG_MAX_INT
		if j > 1 && p.IDArr[j-1] == neID {
			na = constructdbg.GetReverseDBG_MAX_INTArr(p.IDArr[:j+1])
		} else if j < len(p.IDArr)-1 && p.IDArr[j+1] == neID {
			na = p.IDArr[j:]
		}
		if len(na) > 2 {
			pm = append(pm, na)
			if len(na) > maxLen {
				maxLen = len(na)
			}
		}
	}
	// get consis path
	consis := true
	for i := 0; i < maxLen; i++ {
		var id constructdbg.DBG_MAX_INT
		for _, p := range pm {
			if len(p) <= i {
				continue
			}
			if id == 0 {
				id = p[i]
			} else if id != p[i] {
				consis = false
				break
			}
		}
		if consis {
			eIDArr = append(eIDArr, id)
		} else {
			break
		}
	}
	return
}

func GetNGSPath(extArr [2][]constructdbg.DBG_MAX_INT, edgesArr []constructdbg.DBGEdge) (eIDArr []constructdbg.DBG_MAX_INT) {
	if len(extArr[0]) < 2 {
		return
	}
	for i := len(extArr[0]) - 1; i >= 0; i-- {
		e := edgesArr[extArr[0][i]]
		if e.GetUniqueFlag() > 0 && e.GetTwoEdgesCycleFlag() == 0 && i < len(extArr[0])-1 {
			//fmt.Printf("[GetNGSPath] eID: %v, pathMat: %v\n", e.ID, e.PathMat)
			arr := GetConsisPath(e.PathMat, e.ID, extArr[0][i+1])
			if len(arr) > len(extArr[0])-i && reflect.DeepEqual(extArr[0][i:], arr[:len(extArr[0])-i]) {
				eIDArr = arr[len(extArr[0])-i:]
				return
			}
		}
	}
	return
}

func GetMostProbableEdge(edgesArr []constructdbg.DBGEdge, extArr [2][]constructdbg.DBG_MAX_INT, subA []int, arr, na []LRRecord) int {
	j := -1
	if float32(na[0].MapNum)/float32(na[0].GapMapNum) > float32(na[1].MapNum)/float32(na[1].GapMapNum) {
		j = IndexArrPos(subA, arr, na[0].RefID)
	} else if na[0].RefLen < na[1].RefLen {
		j = IndexArrPos(subA, arr, na[0].RefID)
	} else {
		eIDArr := GetNGSPath(extArr, edgesArr)
		if int(eIDArr[0]) > 1 {
			j = IndexArrPos(subA, arr, eIDArr[0])
		}
	}
	return j
}

/*func findMostProbablePath(arr []LRRecord, edgesArr []constructdbg.DBGEdge, nodesArr []constructdbg.DBGNode, opt Options) (extArr [2][]constructdbg.DBG_MAX_INT) {
	pos := findSeedEID(arr, edgesArr, opt.Kmer)
	if pos < 0 || arr[pos].GapMapNum < opt.Kmer {
		return
	}
	fmt.Printf("[findMostProbablePath] found seed Record: %v\n", arr[pos])

	maxFlank := opt.Kmer * 4 / 5
	searchLen := opt.Kmer * 6 / 5
	extArr[0] = append(extArr[0], arr[pos].RefID)
	extArr[1] = append(extArr[1], 0)
	edge := edgesArr[arr[pos].RefID]
	for i := 0; i < 2; i++ { // i == 0  note BACKWARD, i == 1 note FORWARD
		j := pos
		var nID constructdbg.DBG_MAX_INT
		var direction uint8
		if i == 0 {
			if arr[j].Strand == constructdbg.PLUS {
				nID = edge.StartNID
			} else {
				nID = edge.EndNID
			}
			direction = constructdbg.BACKWARD
		} else {
			if arr[j].Strand == constructdbg.PLUS {
				nID = edge.EndNID
			} else {
				nID = edge.StartNID
			}
			direction = constructdbg.FORWARD
		}
		if nID == 0 {
			continue
		}
		e := edge
		for {
			var alter constructdbg.DBG_MAX_INT // alternative edge, maybe cause by haplotype
			var flankLen int
			var readPos int
			if direction == constructdbg.BACKWARD {
				if arr[j].Strand == constructdbg.PLUS {
					flankLen = arr[j].RefStart
				} else {
					flankLen = arr[j].RefLen - arr[j].RefEnd
				}
				readPos = arr[j].Start - flankLen + searchLen
			} else {
				if arr[j].Strand == constructdbg.PLUS {
					flankLen = arr[j].RefLen - arr[j].RefEnd
				} else {
					flankLen = arr[j].RefStart
				}
				readPos = arr[j].End + flankLen - searchLen
			}
			if flankLen >= maxFlank {
				break
			}
			ea := constructdbg.GetNearEdgeIDArr(nodesArr[nID], e.ID)
			if len(ea) == 0 {
				break
			}
			fmt.Printf("[findMostProbablePath]j: %v, direction: %v, readPos: %v, searchLen: %v\n", j, direction, readPos, searchLen)
			subA := SearchLRRecordArr(arr, j, direction, ea, readPos, searchLen)
			if len(subA) == 0 {
				for z, r := range arr {
					fmt.Printf("[findMostProbablePath]arr[%v]: %v\n", z, r)
				}
				log.Fatalf("[findMostProbablePath]not found mapping edge, ea: %v, extArr[0]: %v\n", ea, extArr[0])
			} else if len(subA) == 1 {
				j = subA[0]
			} else { // len(subA) > 1
				na := make([]LRRecord, len(subA))
				for i, idx := range subA {
					na[i] = arr[idx]
				}
				sort.Sort(LRRecordScoreArr(na))
				fmt.Printf("[findMostProbablePath]na: %v\n", na)
				j = GetMostProbableEdge(edgesArr, extArr, subA, arr, na)
				if j < 0 {
					if len(na) == 2 && constructdbg.IsBubble(na[0].RefID, na[1].RefID, edgesArr) {
						j = IndexArrPos(subA, arr, na[0].RefID)
						alter = na[1].RefID
					} else {
						for z, r := range arr {
							fmt.Printf("[findMostProbablePath]arr[%v]: %v\n", z, r)
						}
						log.Fatalf("[findMostProbablePath] encounter unknow case, na: %v\n\textArr: %v\n", na, extArr)
					}
				}
			}

			extArr[0] = append(extArr[0], arr[j].RefID)
			extArr[1] = append(extArr[1], alter)
			e = edgesArr[arr[j].RefID]
			if e.StartNID == nID {
				nID = e.EndNID
			} else {
				nID = e.StartNID
			}
		}

		// extArr process
		if direction == constructdbg.BACKWARD {
			constructdbg.ReverseDBG_MAX_INTArr(extArr[0])
			constructdbg.ReverseDBG_MAX_INTArr(extArr[1])
		}
	}
	return
}*/

func ComputingNextEdgeInfo(e, ed constructdbg.DBGEdge, strand bool, nd constructdbg.DBGNode) (bool, constructdbg.DBG_MAX_INT) {
	var nID constructdbg.DBG_MAX_INT
	if e.StartNID == nd.ID {
		if ed.StartNID == nd.ID {
			strand = !strand
			nID = ed.EndNID
		} else {
			nID = ed.StartNID
		}
	} else {
		if ed.EndNID == nd.ID {
			strand = !strand
			nID = ed.StartNID
		} else {
			nID = ed.EndNID
		}
	}
	return strand, nID
}

/*func GetMappingEdgePathSeq(mi MapInfo, extArr []constructdbg.DBG_MAX_INT, edgesArr []constructdbg.DBGEdge, nodesArr []constructdbg.DBGNode, direction uint8, extendEID constructdbg.DBG_MAX_INT, extendMinLen, kmerLen int) (edgeSeq []byte, strand bool, minLen int) {
	fmt.Printf("[GetMappingEdgePathSeq]mi: %v, direction: %v, extendEID: %v, extendMinLen: %v, kmerLen: %v\n", mi, direction, extendEID, extendMinLen, kmerLen)
	// add first edge sequence
	e := edgesArr[extArr[mi.EIDIdx]]
	strand = mi.Strand
	//var strand bool
	if direction == constructdbg.BACKWARD {
		if strand == constructdbg.PLUS {
			edgeSeq = append(edgeSeq, constructdbg.GetReverseByteArr(e.Utg.Ks[:mi.EdgePos])...)
		} else {
			edgeSeq = append(edgeSeq, constructdbg.GetCompByteArr(e.Utg.Ks[mi.EdgePos:])...)
		}
	} else {
		if strand == constructdbg.PLUS {
			edgeSeq = append(edgeSeq, e.Utg.Ks[mi.EdgePos:]...)
		} else {
			edgeSeq = append(edgeSeq, constructdbg.GetReverseCompByteArr(e.Utg.Ks[:mi.EdgePos])...)
		}
	}
	//fmt.Printf("[GetMappingEdgePathSeq]after add first edge len(edgeSeq): %v\n", len(edgeSeq))
	// add middle path edges sequence
	nID := mi.NID
	for _, eID := range extArr[mi.EIDIdx+1:] {
		nd := nodesArr[nID]
		ed := edgesArr[eID]
		strand, nID = ComputingNextEdgeInfo(e, ed, strand, nd)
		if direction == constructdbg.BACKWARD {
			if strand == constructdbg.PLUS {
				edgeSeq = append(edgeSeq, constructdbg.GetReverseByteArr(ed.Utg.Ks[:len(ed.Utg.Ks)-kmerLen+1])...)
			} else {
				edgeSeq = append(edgeSeq, constructdbg.GetCompByteArr(ed.Utg.Ks[kmerLen-1:])...)
			}
		} else {
			if strand == constructdbg.PLUS {
				edgeSeq = append(edgeSeq, ed.Utg.Ks[kmerLen-1:]...)
			} else {
				edgeSeq = append(edgeSeq, constructdbg.GetReverseCompByteArr(ed.Utg.Ks[:len(ed.Utg.Ks)-kmerLen+1])...)
			}
		}
		extendMinLen -= len(ed.Utg.Ks) - (kmerLen - 1)
		e = ed
	}
	//fmt.Printf("[GetMappingEdgePathSeq]after add middle edge len(edgeSeq): %v\n", len(edgeSeq))
	// add last edge sequence
	ed := edgesArr[extendEID]
	strand, _ = ComputingNextEdgeInfo(e, ed, strand, nodesArr[nID])
	extendMinLen += kmerLen - 1
	//fmt.Printf("[GetMappingEdgePathSeq]after add middle edge extendMinLen: %v\n", extendMinLen)
	if direction == constructdbg.BACKWARD {
		if strand == constructdbg.PLUS {
			edgeSeq = append(edgeSeq, constructdbg.GetReverseByteArr(ed.Utg.Ks[len(ed.Utg.Ks)-extendMinLen:len(ed.Utg.Ks)-kmerLen+1])...)
		} else {
			edgeSeq = append(edgeSeq, constructdbg.GetCompByteArr(ed.Utg.Ks[kmerLen-1:extendMinLen])...)
		}
	} else {
		if strand == constructdbg.PLUS {
			edgeSeq = append(edgeSeq, ed.Utg.Ks[kmerLen-1:extendMinLen]...)
		} else {
			edgeSeq = append(edgeSeq, constructdbg.GetReverseCompByteArr(ed.Utg.Ks[len(ed.Utg.Ks)-extendMinLen:len(ed.Utg.Ks)-kmerLen+1])...)
		}
	}
	//fmt.Printf("[GetMappingEdgePathSeq]after add last edge len(edgeSeq): %v\n", len(edgeSeq))
	minLen = extendMinLen
	return
}*/

func GetMappingReadSeq(direction uint8, mi MapInfo, readSeqBnt []byte, extendMinLen, kmerLen int, edgeSeq []byte) ([]byte, []byte, int, int) {
	var readSeq []byte
	var readLen int
	var nextReadPos, notMappingLen int
	if direction == constructdbg.BACKWARD {
		readLen = mi.ReadPos
	} else {
		readLen = len(readSeqBnt) - mi.ReadPos
	}
	flankLen := 20
	if readLen < len(edgeSeq) {
		if readLen > len(edgeSeq)-(extendMinLen-(kmerLen-1))+flankLen {
			notMappingLen = len(edgeSeq) - readLen
			edgeSeq = edgeSeq[:readLen]
		} else {
			fmt.Printf("[DPAlignEdgePath] readLen: %v < len(edgeSeq): %v\n", readLen, len(edgeSeq))
			return readSeq, edgeSeq, nextReadPos, notMappingLen
		}
	} else {
		if readLen > len(edgeSeq)+len(edgeSeq)/20 {
			readLen = len(edgeSeq) + len(edgeSeq)/20
		}
	}
	if direction == constructdbg.BACKWARD {
		readSeq = constructdbg.GetReverseByteArr(readSeqBnt[mi.ReadPos-readLen : mi.ReadPos])
		nextReadPos = mi.ReadPos - readLen
	} else {
		readSeq = readSeqBnt[mi.ReadPos : mi.ReadPos+readLen]
		nextReadPos = mi.ReadPos + readLen
	}
	return readSeq, edgeSeq, nextReadPos, notMappingLen
}

const PosWidth = 16
const MaxSeedKmerLen = (64 - PosWidth) / bnt.NumBitsInBase
const MUSKPos = (1 << PosWidth) - 1

type SeedKmer struct {
	Seed uint64 // the most left (64 - PosWidth) with kmer sequence, the right partition for kmer position in the sequence
}

type SKArr []SeedKmer

func (arr SKArr) Len() int {
	return len(arr)
}

func (arr SKArr) Less(i, j int) bool {
	return arr[i].Seed < arr[j].Seed
}

func (arr SKArr) Swap(i, j int) {
	arr[i], arr[j] = arr[j], arr[i]
}

type Chain struct {
	X, Y uint32 // note edgesSeq and readSeq seedKmer start position
	Len  uint32 // the length of seed
}

type ChainArr []Chain

func (arr ChainArr) Len() int {
	return len(arr)
}

func (arr ChainArr) Less(i, j int) bool {
	if arr[i].X < arr[j].X {
		return true
	} else if arr[i].X > arr[j].X {
		return false
	} else { // arr[i].X == arr[j].X
		return arr[i].Y <= arr[j].Y
	}
}

func (arr ChainArr) Swap(i, j int) {
	arr[i], arr[j] = arr[j], arr[i]
}

type Score struct {
	Sc  int
	Idx int
}

type ScoreArr []Score

func (arr ScoreArr) Len() int {
	return len(arr)
}

func (arr ScoreArr) Less(i, j int) bool {
	return arr[i].Sc < arr[j].Sc
}

func (arr ScoreArr) Swap(i, j int) {
	arr[i], arr[j] = arr[j], arr[i]
}

func GetKmerList(seq []byte, skLen uint32) (list []SeedKmer) {
	if uint32(len(seq)) < skLen {
		log.Fatalf("[GetKmerList] len(seq): %v < skLen: %v\n", len(seq), skLen)
	}
	if skLen > MaxSeedKmerLen {
		log.Fatalf("[GetKmerList] seedKmer length: %v just allow <=%v, please check\n", skLen, MaxSeedKmerLen)
	}
	if len(seq) >= (1 << PosWidth) {
		log.Fatalf("[GetKmerList] len(seq): %v just allow %v\n", len(seq), 1<<PosWidth)
	}
	list = make([]SeedKmer, len(seq)-int(skLen)+1)
	MUSKSK := (uint64(1)<<(uint64(skLen*bnt.NumBitsInBase)) - 1)
	var sk SeedKmer
	for i := 0; i < int(skLen)-1; i++ {
		sk.Seed <<= bnt.NumBitsInBase
		sk.Seed |= uint64(seq[i])
	}
	for i := int(skLen) - 1; i < len(seq); i++ {
		sk.Seed <<= bnt.NumBitsInBase
		sk.Seed |= uint64(seq[i])
		sk.Seed &= MUSKSK
		j := uint64(i - (int(skLen) - 1))
		list[j].Seed = (sk.Seed << PosWidth) | j
	}
	return
}

func GetInterSection(listA, listB []SeedKmer, skLen uint32) []Chain {
	var interSecArr []Chain
	i, j := 0, 0
	sb := listB[j].Seed >> PosWidth
	for ; i < len(listA) && j < len(listB); i++ {
		sa := listA[i].Seed >> PosWidth
		if sa < sb {
			continue
		} else if sa > sb {
			for ; j < len(listB); j++ {
				sb = listB[j].Seed >> PosWidth
				if sa <= sb {
					break
				}
			}
		}

		if sa < sb {
			continue
		} else if sa > sb {
			break
		}
		// sa == sb
		//fmt.Printf("[GetInterSection] i: %v, j: %v\n", i, j)
		for x := j; x < len(listB); x++ {
			tmp := listB[x].Seed >> PosWidth
			//fmt.Printf("[GetInterSection] tmp: %v,sa: %v,  x: %v\n", tmp, sa, x)
			if tmp > sa {
				break
			}
			var ch Chain
			ch.X, ch.Y = uint32(listA[i].Seed&MUSKPos), uint32(listB[x].Seed&MUSKPos)
			ch.Len = skLen
			//fmt.Printf("[GetInterSection] ch: %v\n", ch)
			interSecArr = append(interSecArr, ch)
		}
	}

	return interSecArr
}

func GetChainArr(interSecArr []Chain, seedKmerLen uint32) (chainA []Chain) {
	for i := 0; i < len(interSecArr); {
		ch := interSecArr[i]
		d := int(ch.X) - int(ch.Y)
		if d < 0 {
			d = int(ch.Y) - int(ch.X)
		}
		j := i + 1
		for ; j < len(interSecArr); j++ {
			if ch.X != interSecArr[j].X {
				break
			}
		}
		for i++; i < j; i++ {
			nch := interSecArr[i]
			td := int(nch.X) - int(nch.Y)
			if td < 0 {
				td = int(nch.Y) - int(nch.X)
			}
			//fmt.Printf("[GetChainArr] interSecArr[%v]: %v, ch: %v\n", i, nch, ch)
			if td < d {
				ch = nch
				d = td
			}
		}

		if d > 50 && d > int(ch.X)/5 {
			continue
		}

		if len(chainA) > 0 {
			lch := chainA[len(chainA)-1]
			//fmt.Printf("[GetChainArr] lch: %v\n", lch)
			if int(ch.Y)-int(lch.Y) < 0 {
				continue
			}
			if int(ch.X)-int(lch.X) < int(lch.Len) || int(ch.Y)-int(lch.Y) < int(lch.Len) {
				if ch.X-lch.X == ch.Y-lch.Y {
					lch.Len += (ch.X + ch.Len - (lch.X + lch.Len))
					chainA[len(chainA)-1] = lch
				}
				continue
			}
		}
		chainA = append(chainA, ch)
	}
	return
}

func AbsInt(a int) int {
	if a >= 0 {
		return a
	} else {
		return 0 - a
	}
}

func FindBestNearChain(chainBlocks []Chain, lastChain Chain, direction uint8, MaxGapLen, GapCost int) (max int, pos int) {
	max = math.MinInt32
	for i := 0; i < len(chainBlocks); i++ {
		ch := chainBlocks[i]
		if direction == constructdbg.BACKWARD {
			if (ch.X+ch.Len > lastChain.X+lastChain.Len) || (ch.Y+ch.Len > lastChain.Y+lastChain.Len) {
				continue
			}
		} else {

		}
		sc := GetNearChainScore(chainBlocks[i], lastChain, direction, MaxGapLen, GapCost)
		if sc > max {
			max = sc
			pos = i
		} else if sc == max {
			if direction == constructdbg.BACKWARD {
				if chainBlocks[pos].X < chainBlocks[i].X {
					pos = i
				}
			} else {
				if chainBlocks[pos].X > chainBlocks[i].X {
					pos = i
				}
			}
		}
	}
	return
}
func GetNearChainScore(nc, lastc Chain, direction uint8, MaxGapLen, GapCost int) (sc int) {
	var gapLen int
	if direction == constructdbg.BACKWARD {
		gapLen = AbsInt((int(lastc.X) - int(nc.X)) - (int(lastc.Y) - int(nc.Y)))
	} else {
		gapLen = AbsInt((int(nc.X) - int(lastc.X)) - (int(nc.Y) - int(lastc.Y)))
	}
	sc = int(nc.Len) - gapLen*GapCost
	/*if gapLen > MaxGapLen {
		sc = math.MinInt32
	} else {
	}*/

	return
}
func GetTwoAnchorScore(lastc, nc Chain, GapCost int) (sc int) {
	var gapLen, distance int
	if lastc.X > nc.X {
		distance = int(lastc.X) - int(nc.X)
		gapLen = AbsInt((int(lastc.X) - int(nc.X)) - (int(lastc.Y) - int(nc.Y)))
	} else {
		distance = int(nc.X) - int(lastc.X)
		gapLen = AbsInt((int(nc.X) - int(lastc.X)) - (int(nc.Y) - int(lastc.Y)))
	}
	if gapLen > distance/5 {
		sc = int(nc.Len) - (gapLen-distance/5)*GapCost
	} else {
		sc = int(nc.Len)
	}

	return
}

func FindInterSecBestChain(chainBlocks []Chain, idx int, cb, lastc Chain, direction uint8, MaxGapLen, GapCost int) (maxY, pY int) {
	maxY = math.MinInt32
	k := idx
	if direction == constructdbg.BACKWARD {
		for ; k >= 0; k-- {
			nc := chainBlocks[k]
			if int(cb.X)-(int(nc.X)+int(nc.Len)) > MaxGapLen {
				break
			}
			if nc.Y+nc.Len >= cb.Y && nc.Y+nc.Len < lastc.Y {
				sc := GetNearChainScore(nc, lastc, constructdbg.BACKWARD, MaxGapLen, GapCost)
				if sc > maxY {
					maxY = sc
					pY = k
				}
			}
		}
	} else { // FORWARD
		for ; k < len(chainBlocks); k++ {
			nc := chainBlocks[k]
			if int(nc.X)-(int(cb.X)+int(cb.Len)) > MaxGapLen {
				break
			}
			if nc.Y <= cb.Y+cb.Len && nc.Y > lastc.Y+lastc.Len {
				sc := GetNearChainScore(nc, lastc, constructdbg.FORWARD, MaxGapLen, GapCost)
				if sc > maxY {
					maxY = sc
					pY = k
				}
			}
		}
	}

	return
}

func GetChainScore(arr []Chain, MaxGapLen, GapCost int) (score int) {
	score = int(arr[0].Len)
	for i := 1; i < len(arr); i++ {
		cb := arr[i-1]
		nc := arr[i]
		score += GetNearChainScore(nc, cb, constructdbg.FORWARD, MaxGapLen, GapCost)
	}
	return score
}

func ReverseChainArr(arr []Chain) {
	for i := 0; i < len(arr)/2; i++ {
		arr[i], arr[len(arr)-1-i] = arr[len(arr)-1-i], arr[i]
	}
}

func IsInIntArr(IdxArr []int, idx int) bool {
	for _, a := range IdxArr {
		if a == idx {
			return true
		}
	}
	return false
}

func SearchNearCb(idxArr []int, chainBlocks []Chain, Idx int, direction uint8) (idx int) {
	if Idx < 0 || Idx >= len(chainBlocks) {
		log.Fatalf("[SearchNearCb] Idx: %d must between [0, %v]\n", Idx, len(chainBlocks))
	}
	if direction == constructdbg.BACKWARD {
		idx = -1
		for i := Idx - 1; i >= 0; i-- {
			//fmt.Printf("[CheckBoundary] Idx: %v, i: %v\n", Idx, i)
			if !IsInIntArr(idxArr, i) {
				continue
			}
			idx = i
			break
		}
	} else {
		idx = len(chainBlocks)
		for i := Idx + 1; i < len(chainBlocks); i++ {
			if !IsInIntArr(idxArr, i) {
				continue
			}
			idx = i
			break
		}
	}
	return
}

func CheckBoundary(chainBlocks []Chain, maxSc Score, visit []bool, arr []int) (cb Chain, sc Score) {
	cb = chainBlocks[maxSc.Idx]
	sc = maxSc
	m := SearchNearCb(arr, chainBlocks, maxSc.Idx, constructdbg.BACKWARD)
	n := SearchNearCb(arr, chainBlocks, maxSc.Idx, constructdbg.FORWARD)
	//fmt.Printf("[CheckBoundary] m: %v, n: %v\n", m, n)
	if m >= 0 {
		nb := chainBlocks[m]
		//fmt.Printf("[CheckBoundary] nb: %v, cb: %v\n", nb, cb)
		if nb.X+nb.Len > cb.X {
			if cb.X+cb.Len <= nb.X+nb.Len {
				sc.Sc = -1
				return
			}
			maxSc.Sc -= int(nb.X + nb.Len - cb.X)
			cb.Len -= (nb.X + nb.Len - cb.X)
			sc.Sc -= int(nb.X + nb.Len - cb.X)
			cb.Y += (nb.X + nb.Len - cb.X)
			cb.X = nb.X + nb.Len
		}
		if cb.Y < nb.Y || cb.Y+cb.Len <= nb.Y+nb.Len {
			sc.Sc = -1
			return
		} else if cb.Y < nb.Y+nb.Len && cb.Y+cb.Len > nb.Y+nb.Len {
			maxSc.Sc -= int(nb.Y + nb.Len - cb.Y)
			cb.Len -= (nb.Y + nb.Len - cb.Y)
			sc.Sc -= int(nb.Y + nb.Len - cb.Y)
			cb.X += (nb.Y + nb.Len - cb.Y)
			cb.Y = nb.Y + nb.Len
		}
	}

	if n < len(chainBlocks) {
		nb := chainBlocks[n]
		//fmt.Printf("[CheckBoundary] nb: %v, cb: %v\n", nb, cb)
		if cb.X+cb.Len > nb.X {
			if cb.X >= nb.X || cb.X+cb.Len >= nb.X+nb.Len {
				sc.Sc = -1
				return
			}
			maxSc.Sc -= int(cb.X + cb.Len - nb.X)
			cb.Len -= (cb.X + cb.Len - nb.X)
			sc.Sc -= int(cb.X + cb.Len - nb.X)
			//cb.X = nb.X - cb.Len
		}
		if cb.Y >= nb.Y || cb.Y+cb.Len >= nb.Y+nb.Len {
			sc.Sc = -1
			return
		} else if cb.Y+cb.Len > nb.Y {
			maxSc.Sc -= int(cb.Y + cb.Len - nb.Y)
			cb.Len -= (cb.Y + cb.Len - nb.Y)
			sc.Sc -= int(cb.Y + cb.Len - nb.Y)
			//cb.Y = nb.Y - cb.Len
		}
	}
	return
}

func BlockToAnchorsScoreArr(chainBlocks []Chain, idx int, visit []bool, GapCost, MaxGapLen int) (baSA []Score) {
	cen := chainBlocks[idx]
	baSA = make([]Score, 0, 15)
	// BACKWARD
	for i := idx - 1; i >= 0; i-- {
		if visit[i] == true {
			continue
		}
		cb := chainBlocks[i]
		sc := GetTwoAnchorScore(cen, cb, GapCost)
		if sc > 0 {
			var s Score
			s.Sc = sc
			s.Idx = i
			baSA = append(baSA, s)
		}
		if int(cen.X)-int(cb.X) > MaxGapLen {
			break
		}
	}
	//FORWARD
	for i := idx + 1; i < len(chainBlocks); i++ {
		if visit[i] == true {
			continue
		}
		cb := chainBlocks[i]
		sc := GetTwoAnchorScore(cen, cb, GapCost)
		if sc > 0 {
			var s Score
			s.Sc = sc
			s.Idx = i
			baSA = append(baSA, s)
		}
		if int(cb.X+cb.Len)-int(cen.X) > MaxGapLen {
			break
		}
	}
	sort.Sort(ScoreArr(baSA))
	return
}

func GetMaxScoreFromBTASP(blockToAnchorsScorePat [][]Score, visit []bool) (max Score) {
	for i, sa := range blockToAnchorsScorePat {
		if len(sa) <= 0 {
			continue
		}
		for j := len(sa) - 1; j >= 0; j-- {
			sc := sa[j]
			if visit[sc.Idx] == true {
				blockToAnchorsScorePat[i] = sa[:j]
				continue
			}

			if max.Sc < sc.Sc {
				max = sc
			}
			break
		}
	}
	return
}

func GetMaxChainArr(interSecSortArr []Chain, strand bool, seedKmerLen uint32) (maxChainArr []Chain, maxScore int) {
	var chainBlocks []Chain
	//var stack []int // restore the need retraverse index of interSecSortArr
	flagArr := make([]bool, len(interSecSortArr)) // denote corresponding chain has been used
	loopNum := 0
	if strand == constructdbg.MINUS {
		log.Fatalf("[GetMaxChainArr] mapping by '-' strand unprocessed\n")
	}
	for { // until to the no flag false chain
		loopNum++
		chl := len(chainBlocks)
		for i := 0; i < len(flagArr); {
			if flagArr[i] == true {
				i++
				continue
			}
			cb := interSecSortArr[i]
			flagArr[i] = true
			j := i + 1
			for ; j < len(interSecSortArr); j++ {
				if flagArr[j] == true {
					continue
				}
				nc := interSecSortArr[j]
				nextX := cb.X + cb.Len
				nextY := cb.Y + cb.Len
				if nc.X > nextX {
					break
				}
				// nc.X == nextX
				if nc.Y < cb.Y {
					continue
				}
				if nc.Y > nextY {
					continue
				}
				//
				if nc.X-cb.X == nc.Y-cb.Y {
					cb.Len += (nc.X + nc.Len) - (cb.X + cb.Len)
					flagArr[j] = true
				} else {
					break
				}
			}
			chainBlocks = append(chainBlocks, cb)
			i = j
		}

		// no new added chainBlock
		if chl == len(chainBlocks) {
			break
		}
	}

	// check code
	/*fmt.Printf("[GetMaxChainArr] loop number is : %d\n", loopNum)
	for i, fg := range flagArr {
		if fg == false {
			log.Fatalf("[GetMaxChainArr] interSecSortArr have not been finished processed, flagArr[%d]\n", i)
		}
	}*/

	/*for i := 0; i < len(chainBlocks)-1; i++ {
		cb := chainBlocks[i]
		nc := chainBlocks[i+1]
		if cb.X+cb.Len > nc.X || cb.Y+cb.Len > nc.Y {
			x := int(cb.X) + int(cb.Len) - int(nc.X)
			y := int(cb.Y) + int(cb.Len) - int(nc.Y)
			max := x
			if x < y {
				max = y
			}
			cb.Len -= uint32(max)
			chainBlocks[i] = cb
		}
	}*/

	sort.Sort(ChainArr(chainBlocks))
	scArr := make([]Score, len(chainBlocks))
	for i, cb := range chainBlocks {
		scArr[i].Sc = int(cb.Len)
		scArr[i].Idx = i
	}
	//fmt.Printf("[GetMaxChainArr] chainBlocks: %v\n", chainBlocks)
	sort.Sort(ScoreArr(scArr))
	// find max score
	GapCost := 1       // one base Gap cost
	MaxGapLen := 1000  // allow max Gap length
	MaxRepeatNum := 10 // max score have been repeat display times
	maxScore1 := math.MinInt32
	//var maxScore2 int
	var maxNum int // note maxScore1 has been repeat times
	var maxArr1 []Chain
	for i := len(scArr) - 1; i >= 0; i-- {
		idx := scArr[i].Idx
		var arr []int
		var chainArr []Chain
		var blockToAnchorsScorePat [][]Score
		visit := make([]bool, len(chainBlocks))
		visit[idx] = true
		arr = append(arr, idx)
		chainArr = append(chainArr, chainBlocks[idx])
		baSA := BlockToAnchorsScoreArr(chainBlocks, idx, visit, GapCost, MaxGapLen)
		blockToAnchorsScorePat = append(blockToAnchorsScorePat, baSA)
		for {
			maxSc := GetMaxScoreFromBTASP(blockToAnchorsScorePat, visit)
			//fmt.Printf("[GetMaxChainArr] blockToAnchorsScorePat: %v\n\tarr: %v\n", blockToAnchorsScorePat, arr)
			visit[maxSc.Idx] = true
			if maxSc.Sc <= 0 {
				break
			}
			//fmt.Printf("[GetMaxChainArr] maxSc: %v\n", maxSc)
			// check boundary
			//cb := chainBlocks[maxSc.Idx]
			cb, sc := CheckBoundary(chainBlocks, maxSc, visit, arr)
			//fmt.Printf("[GetMaxChainArr] cb:%v, sc: %v\n", cb, sc)
			if sc.Sc <= 0 {
				continue
			}
			if cb.Len < seedKmerLen {
				continue
			}
			// add to Max Chain
			arr = append(arr, sc.Idx)
			chainArr = append(chainArr, cb)
			baSA := BlockToAnchorsScoreArr(chainBlocks, sc.Idx, visit, MaxGapLen, GapCost)
			blockToAnchorsScorePat = append(blockToAnchorsScorePat, baSA)
		}
		sort.Sort(ChainArr(chainArr))
		arrScore := GetChainScore(chainArr, MaxGapLen, GapCost)
		//fmt.Printf("[GetMaxChainArr]score: %v\tchainArr: %v\n", arrScore, chainArr)
		if arrScore > maxScore1 {
			maxScore1 = arrScore
			maxArr1 = chainArr
			maxNum = 0
		} else if arrScore == maxScore1 {
			if len(chainArr) == len(maxArr1) && reflect.DeepEqual(chainArr, maxArr1) {
				maxNum++
			} else {
				//maxScore2 = arrScore
				//maxArr2 = chainArr
				//fmt.Printf("[GetMaxChainArr]display two same max score\n\tarrScore2:%v, maxArr2: %v\n\tarrScoce1: %v, maxArr1: %v\n", maxScore2, maxArr2, maxScore1, maxArr1)
			}
		} else {
			maxNum++
		}

		if maxNum >= MaxRepeatNum {
			break
		}
	}

	maxChainArr = maxArr1
	maxScore = maxScore1
	return
	/*for j := idx - 1; j >= 0; {
		k := j - 1
		lastc := arr[len(arr)-1]
		for ; k >= 0; k-- {
			if chainBlocks[k].X+chainBlocks[k].Len < chainBlocks[j].X {
				break
			}
		}
		max, p := FindBestNearChain(chainBlocks[k+1:j+1], lastc, constructdbg.BACKWARD, MaxGapLen, GapCost)
		// found Y side have intersection with X
		p += k + 1
		cb := chainBlocks[p]
		maxY, pY := FindInterSecBestChain(chainBlocks, p-1, cb, lastc, constructdbg.BACKWARD, MaxGapLen, GapCost)
		if maxY > max {
			p = pY
		}
		arr = append(arr, chainBlocks[p])
		j = p - 1
	}
	ReverseChainArr(arr)

	// FORWARD search
	for j := idx + 1; j < len(chainBlocks); {
		k := j + 1
		lastc := arr[len(arr)-1]
		for ; k < len(chainBlocks); k++ {
			if chainBlocks[k].X > chainBlocks[j].X+chainBlocks[j].Len {
				break
			}
		}
		max, p := FindBestNearChain(chainBlocks[j:k], lastc, constructdbg.FORWARD, MaxGapLen, GapCost)
		// found Y side have intersection with X
		p += j
		cb := chainBlocks[p]
		maxY, pY := FindInterSecBestChain(chainBlocks, p+1, cb, lastc, constructdbg.FORWARD, MaxGapLen, GapCost)
		if maxY > max {
			p = pY
		}
		arr = append(arr, chainBlocks[p])
		j = p + 1
	}*/
}

type CIGAR struct {
	Ins, Del, Mis, Mch uint32
}

func Max(a, b, c int) int {
	max := a
	if b > max {
		max = b
	}
	if c > max {
		max = c
	}

	return max
}

func GetKArr(cycle, D uint, low, high int) []int {
	kArr := make([]int, high-low+1)
	i := 0
	if cycle&0x1 == 0 {
		for j := low; j <= high; j++ {
			if uint(int(D)-j)&0x1 == 0 {
				kArr[i] = j
				i++
			}
		}
		for j := low; j <= high; j++ {
			if uint(int(D)-j)&0x1 == 1 {
				kArr[i] = j
				i++
			}
		}
	} else {
		for j := low; j <= high; j++ {
			if uint(int(D)-j)&0x1 == 1 {
				kArr[i] = j
				i++
			}
		}
		for j := low; j <= high; j++ {
			if uint(int(D)-j)&0x1 == 0 {
				kArr[i] = j
				i++
			}
		}
	}
	return kArr
}

func PrintBvec(b uint64, l int) []byte {
	if l > 64 {
		l = 64
	}
	bi := make([]byte, l)
	for i := 0; i < l; i++ {
		if (b & 0x1) > 0 {
			bi[i] = '1'
		} else {
			bi[i] = '0'
		}
		b >>= 1
	}
	bi = []byte(Reverse(string(bi)))
	return bi
}

func GlobalAlignment(seqa, seqb []byte, localFirst bool) (cg CIGAR, lastY int) {
	//var lastY int
	//C := uint(64)
	LAG := 20
	//MASKC := 1 << (C - 1)
	var max, min int
	if len(seqa) >= len(seqb) {
		max = len(seqa)
		min = len(seqb)
	} else {
		max = len(seqb)
		min = len(seqa)
	}
	var D uint32
	if max < 200 {
		D = uint32(max) + 2
	} else {
		D = uint32(max) / 2
	}
	if D < 10 {
		D = 10
	}
	if max-min > int(D)-min/2 {
		D = uint32(max - min + min/2)
	}
	if D&0x1 == 1 { // D must even
		D++
	}
	mid := int(D)
	width := 2*D + 1
	B := make([]uint64, width) // diagonal value is D
	M := make([]int, width)    // store match score
	W := make([]int, width)    // store max x+y == (y+k) + y == 2* y + k
	for i, _ := range W {
		W[i] = -1
	}
	//W[mid] = -1
	low, high := 1, -1
	var finished bool
	//fmt.Printf("[GlobalAlignment]mid: %v\n", mid)
	maxBound, minBound := 0, 0
	i := uint32(0)
	for ; i <= D; i++ {
		//first process the digits propery(even or odd) same as i
		//kArr := GetKArr(i, D, low-1, high+1)
		//for _, k := range kArr
		//fmt.Printf("kArr: %v\n", kArr)
		/*start := low - 1
		if ((D - (uint(low) - 1)) & 0x1) != (i & 0x1) {
			start++
		}*/
		var lastB uint64
		var lastM int
		var ap, ac, am int
		high++
		low--
		if high <= maxBound {
			W[mid+high] = -1
			M[mid+high] = 0
			B[mid+high] = 0
		}
		am = W[mid+high]
		lastM, lastB = 0, 0
		if low >= minBound {
			M[mid+low], M[mid+low-1] = 0, 0
			W[mid+low], W[mid+low-1] = -1, -1
			B[mid+low], B[mid+low-1] = 0, 0
		}
		ac = -1
		for k := high; k >= low; k-- { // k = i - j
			ap = ac
			ac = am
			am = W[mid+k-1]
			max := Max(ap, ac, am)
			/*if (i & 0x1) != (uint(int(D)-k) & 0x1) { // D and k one is odd ,another is even
				if max != W[mid+k]+1 {
					continue
				}
			}*/
			var a, m int
			var b uint64
			//fmt.Printf("k: %v,ap: %v, m: %v,  ac: %v, m: %v, am: %v,m: %v, max: %v\n", k, ap, M[mid+k+1], ac, M[mid+k], am, M[mid+k-1], max)
			// check have been encounter boundary
			var tx, ty int
			if max >= 0 && max >= k {
				/*if max < k {
					log.Fatalf("[GlobalAlignment] max:%v must bigger than k: %v\n", max, k)
				}*/
				ty = (max - k) >> 1
				tx = ty + k
			} else {
				tx, ty = -1, -1
			}
			if tx >= len(seqa) {
				a, m, b = max+1, M[mid+k-1], B[mid+k-1]<<1
			} else if ty >= len(seqb) {
				a, m, b = max+1, lastM, lastB
			} else {
				if max == ap {
					a, m, b = max+1, lastM, lastB
				} else if max == am {
					a, m, b = max+1, M[mid+k-1], B[mid+k-1]<<1
				} else {
					a, m, b = max+2, M[mid+k], B[mid+k]<<1
				}
			}
			y := -1
			if a >= k {
				y = (a - k) >> 1
			}
			if y < 0 || y+k < 0 {
				log.Fatalf("[GlobalAlignment] y: %v, y+k: %v must >= 0\n", y, y+k)
			}
			//fmt.Printf("[GlobalAlignment]y: %v, y+k: %v, a: %v, m: %v\n", y, y+k, a, m)
			/*if y < 0 || y+k < 0 {
				continue
			}*/
			for ; y < len(seqb) && y+k < len(seqa) && seqb[y] == seqa[y+k]; y++ {
				b = (b << 1) | 1
				m++
				a += 2
				//fmt.Printf("[GlobalAlignment]seqb[y] == seqa[y+k], y: %v, y+k: %v\n", y, y+k)
			}

			//fmt.Printf("[GlobalAlignment]a: %v, y: %v, y+k: %v\n", a, y, y+k)
			/*if y > len(seqb) || y+k > len(seqa) {
				continue
			}*/

			if localFirst && y < len(seqb) && y+k == len(seqa) {
				cg.Mch = uint32(m)
				cg.Mis = i
				lastY = y
				finished = true
				//fmt.Printf("locaFirst, y: %v, y+k: %v, cg: %v\n", y, y+k, cg)
				break
			} else if y == len(seqb) && y+k == len(seqa) {
				cg.Mch = uint32(m)
				cg.Mis = i
				lastY = y
				finished = true
				//fmt.Printf("right finished, y: %v, y+k: %v, cg: %v\n", y, y+k, cg)
				//fmt.Println("B:" + string(PrintBvec(b, len(seqa))))
				break
			}
			lastB, lastM = B[mid+k], M[mid+k]
			W[mid+k], M[mid+k], B[mid+k] = a, m, b
		}

		if finished {
			break
		}

		// trimming furthest reaching points
		var besta, besty int
		for k := low; k <= high; k++ {
			if besta < W[mid+k] { // i + j == y+k + y
				besta = W[mid+k]
				besty = (W[mid+k] - k) >> 1 // W[mid+k] = 2*y + k
			}
		}

		if LAG > 10 && LAG < besty/5 {
			LAG = besty / 5
		}
		for ; low <= high && W[mid+low] <= besta-LAG; low++ {
		}
		for ; high >= low && W[mid+high] <= besta-LAG; high-- {
		}

		if low-2 <= 0-int(D) || high+1 >= int(D) || low > high {
			break
			//log.Fatalf("[GlobalAlignment] low:%v, high: %v must between (-%v ~ %v)\n", low, high, D, D)
		}
		if minBound > low {
			minBound = low
		}
		if maxBound < high {
			maxBound = high
		}

		if localFirst && i > 50 && besty < 2*int(i) {
			break
		}
		//fmt.Printf("i: %v, low: %v, high: %v\nW: %v\n", i, low, high, W)
		//fmt.Printf("M: %v\n", M)
	}
	if !finished {
		maxM, maxa := 0, 0
		//var b uint64
		for j, mch := range M {
			if mch > 0 && mch >= maxM {
				if mch == maxM && W[j] >= maxa {
					continue
				}
				k := j - mid
				if W[j] > k {
					lastY = (W[j] - k) >> 1
				} else {
					lastY = -1
				}
				maxM = mch
				maxa = W[j]
				//b = B[j]
				//cg.Mis = uint32(len(seqa) - (lastY + j - mid))
				//fmt.Printf("!!!!!unfinished\tmaxM: %v, y: %v, y+k: %v\n", maxM, lastY, lastY+j-mid)
			}
		}
		cg.Mch = uint32(maxM)
		cg.Mis = i
		fmt.Printf("!!!!!unfinished\tmaxM: %v, cg: %v\n", maxM, cg)
		//fmt.Println("B:" + string(PrintBvec(b, len(seqa))))
		if !localFirst {
			log.Fatalf("[GlobalAlignment] not found proper alignment, i: %v,  \nseqa: %v\nseqb: %v\n", i, seqa, seqb)
		}
	}
	// get score
	// i note the alignment encounter total number base of  MisMatch+Insertion+Deletion, any Mis or InDel base penlize score -1
	//cg.Mis += uint32(i)
	return
}

func DPLocalAlign(edgesSeq, readSeq []byte, chainA []Chain) (cg CIGAR) {
	//enlongation := uint32(200) // same as minimap2 mapping step '-g' argument
	for i := 0; i < len(chainA); i++ {
		var pch Chain
		if i == 0 {
			pch.X, pch.Y, pch.Len = 0, 0, 0
		} else {
			pch = chainA[i-1]
		}
		ch := chainA[i]
		pXEnd := pch.X + pch.Len
		pYEnd := pch.Y + pch.Len
		xlen, ylen := int(ch.X)-int(pXEnd), int(ch.Y)-int(pYEnd)
		//fmt.Printf("[DPLocalAlign] len(edgesSeq): %v, len(readSeq): %v, pch: %v, ch: %v\n", xlen, ylen, pch, ch)
		if xlen < 0 || ylen < 0 {
			if xlen < ylen {
				cg.Ins += uint32(ylen - xlen)
			} else {
				cg.Del += uint32(xlen - ylen)
			}
			max := math.MinInt32
			if xlen < 0 {
				max = AbsInt(xlen)
			}
			if ylen < 0 && max < AbsInt(ylen) {
				max = AbsInt(ylen)
			}
			if ch.Len <= uint32(max) {
				log.Fatalf("[DPLocalAlign] ch.Len: %v <= max: %v, pch: %v, ch: %v\n", ch.Len, max, pch, ch)
			}
			cg.Mch += (ch.Len - uint32(max))
			continue
		}

		cg.Mch += ch.Len
		if pXEnd == ch.X {
			cg.Ins += ch.Y - pYEnd
			continue
		} else if pYEnd == ch.Y {
			cg.Del += ch.X - pXEnd
			continue
		} else if xlen != ylen {
			var finished bool // test if a substring in another
			if xlen < ylen {
				for j := 1; j <= ylen-xlen; j++ {
					if reflect.DeepEqual(edgesSeq[pXEnd:ch.X], readSeq[int(pYEnd)+j:int(pYEnd)+j+xlen]) {
						cg.Ins += uint32(ylen - xlen)
						cg.Mch += uint32(xlen)
						finished = true
						break
					}
				}
			} else {
				for j := 1; j <= xlen-ylen; j++ {
					if reflect.DeepEqual(readSeq[pYEnd:ch.Y], edgesSeq[int(pXEnd)+j:int(pXEnd)+j+ylen]) {
						cg.Del += uint32(xlen - ylen)
						cg.Mch += uint32(ylen)
						finished = true
						break
					}
				}
			}
			if finished {
				continue
			}
		} else { // ch.X - pXEnd == ch.Y - pYEnd
			if xlen == 1 {
				if edgesSeq[pXEnd] != readSeq[pYEnd] {
					cg.Mis++
					continue
				} else {
					cg.Mch++
					//log.Fatalf("[DPLocalAlign] found  between two chain region have 100 percent identity sequence, %v == %v\n", edgesSeq[pXEnd], readSeq[pYEnd])
				}
			} else if xlen == 2 {
				if edgesSeq[pXEnd] == readSeq[pYEnd+1] || edgesSeq[pXEnd+1] == readSeq[pYEnd] {
					cg.Mis++
					cg.Mch++
				} else {
					cg.Mis += 2
				}
				continue
			} else if i > 0 && (int(ch.X)-1)-(int(pXEnd)+1) > 0 && (int(ch.X)-1)-(int(pXEnd)+1) == (int(ch.Y)-1)-(int(pYEnd)+1) && reflect.DeepEqual(edgesSeq[pXEnd+1:ch.X-1], readSeq[pYEnd+1:ch.Y-1]) { // just the first and last diffence
				cg.Mis += 2
				cg.Mch += (ch.X - 1 - (pXEnd + 1))
				continue
			}
		}
		//fmt.Printf("[DPLocalAlign] cg: %v,ch: %v\n", cg, ch)
		//fmt.Printf("[DPLocalAlign] edge part seq: %v\n", edgesSeq[pXEnd:ch.X])
		//fmt.Printf("[DPLocalAlign] read part seq: %v\n", readSeq[pYEnd:ch.Y])
		// need DP alignment
		localFirst := true
		cigar, _ := GlobalAlignment(edgesSeq[pXEnd:ch.X], readSeq[pYEnd:ch.Y], localFirst)
		//fmt.Printf("[DPLocalAlign]ch: %v, cigar: %v, cg: %v, y: %v\n", ch, cigar, cg, y)
		cg.Ins += cigar.Ins
		cg.Del += cigar.Del
		cg.Mis += cigar.Mis
		cg.Mch += cigar.Mch
	}

	return
}

/*func GetChainScore(chainA []Chain) (chainScore int) {
	for _, ch := range chainA {
		chainScore += int(ch.Len)
	}
	return
}*/

func DPAlignEdgePathExtend(pathSeq, readSeq []byte, edgesArr []constructdbg.DBGEdge, nodesArr []constructdbg.DBGNode) (cg CIGAR, chainScore int, alignEdgePos, alignReadPos int) {
	//debug code
	//readSeq, edgesSeq = readSeq[:651], edgesSeq[:653]
	//fmt.Printf("[DPAlignEdgePath] len(edgesSeq): %v, len(readSeq): %v\n", len(edgesSeq), len(readSeq))
	//fmt.Printf("[DPAlignEdgePath]edgesSeq string: %v\n", bnt.TransformSeq(pathSeq))
	//fmt.Printf("[DPAlignEdgePath]readSeq string: %v\n", bnt.TransformSeq(readSeq))
	//fmt.Printf("[DPAlignEdgePath] reverse edgesSeq string: %v\n\t\t reverse readSeq string: %v\n", bnt.TransformSeq(constructdbg.GetReverseByteArr(edgesSeq)), bnt.TransformSeq(constructdbg.GetReverseByteArr(readSeq)))

	// found two sequences seed chaining
	seedKmerLen := uint32(15)
	listA := GetKmerList(pathSeq, seedKmerLen)
	listB := GetKmerList(readSeq, seedKmerLen)
	/*i := 0
	for ; i < len(listA) && i < len(listB); i++ {
		fmt.Printf("[GetChainArr] listA[%v]: %v:%v\t listB[%v]: %v:%v\n", i, listA[i].Seed>>PosWidth, listA[i].Seed&MUSKPos, i, listB[i].Seed>>PosWidth, listB[i].Seed&MUSKPos)
	}
	if i < len(listA) {
		fmt.Printf("[GetChainArr] listA[%v:]: %v\n", i, listA[i:])
	} else if i < len(listB) {
		fmt.Printf("[GetChainArr]\t\t\tlistB[%v:]: %v\n", i, listB[i:])
	}*/
	sort.Sort(SKArr(listA))
	sort.Sort(SKArr(listB))
	interSecArr := GetInterSection(listA, listB, seedKmerLen)
	sort.Sort(ChainArr(interSecArr))
	/*for i, ch := range interSecArr {
		fmt.Printf("[GetChainArr] interSecArr[%v]: %v\n", i, ch)
	}*/
	maxChain, maxScore := GetMaxChainArr(interSecArr, constructdbg.PLUS, seedKmerLen)
	//chainA := GetChainArr(interSecArr, seedKmerLen)
	//chainScore = GetChainScore(chainA)
	chainScore = maxScore
	//fmt.Printf("[DPAlignEdgePath] MaxScore: %v\n", maxScore)
	//fmt.Printf("[DPAlignEdgePath] edgesSeq: %v\n\t\t\treadSeq: %v\nchainScore: %v\n", edgesSeq, readSeq, chainScore)
	//fmt.Printf("[DPAlignEdgePath] maxChain: %v\n", maxChain)
	AlignmentMinLen := constructdbg.Min(len(pathSeq), len(readSeq))
	if int(maxScore) < AlignmentMinLen/15 {
		fmt.Printf("[DPAlignEdgePath] chainScore: %v < min(len(edgesSeq),len(readSeq))/15: %v\n", chainScore, AlignmentMinLen/15)
		return
	}
	lch := maxChain[len(maxChain)-1]
	alignEdgePos = int(lch.X + lch.Len)
	alignReadPos = int(lch.Y + lch.Len)
	//notChainEdgeLen = notMappinglen + (len(edgesSeq) - int(lch.X+lch.Len))

	// Rapid Local Alignment Discovery from seeds chain
	cg = DPLocalAlign(pathSeq, readSeq, maxChain)

	//fmt.Printf("[DPAlignEdgePath] cg: %v\n", cg)
	//score = int(cg.Mch) - int(cg.Ins+cg.Del+cg.Mis)
	return
}

func GetMinLen(ea []constructdbg.DBG_MAX_INT, edgesArr []constructdbg.DBGEdge, kmerLen int) (minLen int) {
	minLen = math.MaxInt64
	for _, id := range ea {
		if len(edgesArr[id].Utg.Ks)-(kmerLen-1) < minLen {
			minLen = len(edgesArr[id].Utg.Ks) - (kmerLen - 1)
		}
	}
	return
}

func IndexCigarArr(arr []CIGAR, score int) int {
	for i, v := range arr {
		if score == int(v.Mch) {
			return i
		}
	}
	return -1
}

/*func GetMappingExtendPathArr(extArr, ea []constructdbg.DBG_MAX_INT, nID constructdbg.DBG_MAX_INT, edgesArr []constructdbg.DBGEdge, nodesArr []constructdbg.DBGNode, MinExtendLen, kmerlen int) (mappingPathArr [][]constructdbg.DBG_MAX_INT, lastEdgeArr []constructdbg.DBG_MAX_INT, newMinLen int) {
	type pathExtend struct {
		Path      []constructdbg.DBG_MAX_INT
		NID       constructdbg.DBG_MAX_INT
		ExtendLen int
	}
	newMinLen = math.MaxInt64
	var pathQueue []pathExtend
	for _, eID := range ea {
		var pe pathExtend
		pe.Path = append(pe.Path, extArr...)
		pe.Path = append(pe.Path, eID)
		te := edgesArr[eID]
		if te.StartNID == nID {
			pe.NID = te.EndNID
		} else {
			pe.NID = te.StartNID
		}
		pe.ExtendLen = len(te.Utg.Ks) - (kmerlen - 1)
		pathQueue = append(pathQueue, pe)
	}
	for len(pathQueue) > 0 {
		pe := pathQueue[len(pathQueue)-1]
		pathQueue = pathQueue[:len(pathQueue)-1]
		if pe.ExtendLen > MinExtendLen {
			mappingPathArr = append(mappingPathArr, pe.Path[:len(pe.Path)-1])
			leID := pe.Path[len(pe.Path)-1]
			lastEdgeArr = append(lastEdgeArr, leID)
			if pe.ExtendLen < newMinLen {
				newMinLen = pe.ExtendLen
			}
			continue
		}
		// add next edge
		tea := constructdbg.GetNearEdgeIDArr(nodesArr[pe.NID], pe.Path[len(pe.Path)-1])
		for _, eID := range tea {
			ne := edgesArr[eID]
			var npe pathExtend
			npe.Path = append(npe.Path, pe.Path...)
			npe.Path = append(npe.Path, eID)
			if ne.StartNID == pe.NID {
				npe.NID = ne.EndNID
			} else {
				npe.NID = ne.StartNID
			}
			npe.ExtendLen = pe.ExtendLen + len(ne.Utg.Ks) - (kmerlen - 1)
			pathQueue = append(pathQueue, npe)
		}
	}
	return
}*/

/*func BackWardDeduceStrand(mappingPath []constructdbg.DBG_MAX_INT, consisPathLen int, edgesArr []constructdbg.DBGEdge, nID constructdbg.DBG_MAX_INT, strand bool) bool {
	for y := len(mappingPath) - 1; y >= consisPathLen; y-- {
		ae := edgesArr[mappingPath[y]]
		be := edgesArr[mappingPath[y-1]]
		if ae.StartNID == nID {
			if be.StartNID == nID {
				strand = !strand
				nID = be.EndNID
			} else {
				nID = be.StartNID
			}
		} else {
			if be.EndNID == nID {
				strand = !strand
				nID = be.StartNID
			} else {
				nID = be.EndNID
			}
		}
	}
	return strand
}*/

func GetLongConsisPath(mappingPathArr [][]constructdbg.DBG_MAX_INT, idxArr []int) (consisPath []constructdbg.DBG_MAX_INT) {
	for i := 0; ; i++ {
		var eID constructdbg.DBG_MAX_INT
		for _, idx := range idxArr {
			if len(mappingPathArr[idx]) <= i {
				eID = 0
				break
			}
			if eID == 0 {
				eID = mappingPathArr[idx][i]
			} else if eID != mappingPathArr[idx][i] {
				eID = 0
				break
			}
		}
		if eID > 0 {
			consisPath = append(consisPath, eID)
		} else {
			break
		}
	}

	return
}

func GetExtendLen(extArr []constructdbg.DBG_MAX_INT, idx int, edgesArr []constructdbg.DBGEdge, kmerLen int) (extendLen int) {
	for x := idx; x < len(extArr); x++ {
		te := edgesArr[extArr[x]]
		extendLen += len(te.Utg.Ks) - (kmerLen - 1)
	}
	return
}

/*func GetEdgePos(e constructdbg.DBGEdge, direction uint8, strand bool, kmerLen, minLen, flankLen int) int {
	if direction == constructdbg.BACKWARD {
		if strand == constructdbg.PLUS {
			return len(e.Utg.Ks) - (kmerLen - 1) - minLen + flankLen
		} else {
			return (kmerLen - 1) + minLen - flankLen
		}
	} else {
		if strand == constructdbg.PLUS {
			return (kmerLen - 1) + minLen - flankLen
		} else {
			return len(e.Utg.Ks) - (kmerLen - 1) - minLen + flankLen
		}
	}
}*/

type Item struct {
	EIDArr []constructdbg.DBG_MAX_INT
	//NID       constructdbg.DBG_MAX_INT
	ExtendLen int
	Strand    bool
}

func GetComingEdgeArr(edgeIDComing [bnt.BaseTypeNum]constructdbg.DBG_MAX_INT) (ea []constructdbg.DBG_MAX_INT) {
	for _, id := range edgeIDComing {
		if id > 1 {
			ea = append(ea, id)
		}
	}
	return
}

func GetExtendMappingPathArr(ea []constructdbg.DBG_MAX_INT, strArr []bool, direction uint8, extMinLen int, nodesArr []constructdbg.DBGNode, edgesArr []constructdbg.DBGEdge, kmerlen int) (extPathArr [][]constructdbg.DBG_MAX_INT) {
	var stack []Item
	for i, id := range ea {
		var t Item
		e := edgesArr[id]
		t.Strand = strArr[i]
		t.EIDArr = append(t.EIDArr, id)
		t.ExtendLen += len(e.Utg.Ks) - (kmerlen - 1)
		stack = append(stack, t)
	}

	// process stack
	for len(stack) > 0 {
		t := stack[len(stack)-1]
		stack = stack[:len(stack)-1]
		if t.ExtendLen >= extMinLen {
			extPathArr = append(extPathArr, t.EIDArr)
			continue
		}
		e := edgesArr[t.EIDArr[len(t.EIDArr)-1]]
		_, ea, strArr := GetNextEdgeInfo(direction, e, t.Strand, nodesArr, edgesArr)
		for j, id := range ea {
			var nt Item
			nt.EIDArr = make([]constructdbg.DBG_MAX_INT, len(t.EIDArr)+1)
			copy(nt.EIDArr[:len(t.EIDArr)], t.EIDArr)
			nt.EIDArr[len(nt.EIDArr)-1] = id
			nt.Strand = strArr[j]
			e := edgesArr[id]
			nt.ExtendLen = t.ExtendLen + (len(e.Utg.Ks) - (kmerlen - 1))
			stack = append(stack, nt)
		}
	}

	return
}
func GetPathSeqArr(mi MapInfo, initPath []constructdbg.DBG_MAX_INT, direction uint8, extendPathArr [][]constructdbg.DBG_MAX_INT, edgesArr []constructdbg.DBGEdge, nodesArr []constructdbg.DBGNode, kmerlen int, extLen int) ([][]byte, []int, []bool, int) {
	// first Get initArr Sequence
	var sharePathLen int
	seqArr := make([][]byte, len(extendPathArr))
	edgeEndPosArr := make([]int, len(extendPathArr))
	strandArr := make([]bool, len(extendPathArr))
	var initSeq []byte
	var e constructdbg.DBGEdge
	strand := mi.Strand
	//fmt.Printf("[GetPathSeqArr]mi: %v, initPath: %v, direction:%v\n", mi, initPath, direction)
	if direction == constructdbg.FORWARD {
		for i, eID := range initPath {
			e = edgesArr[eID]
			if i == 0 {
				if strand {
					initSeq = append(initSeq, e.Utg.Ks[mi.EdgePos:]...)
				} else {
					initSeq = append(initSeq, constructdbg.GetReverseCompByteArr(e.Utg.Ks[:mi.EdgePos])...)
				}
			} else {
				if strand {
					initSeq = append(initSeq, e.Utg.Ks[kmerlen-1:]...)
				} else {
					initSeq = append(initSeq, constructdbg.GetReverseCompByteArr(e.Utg.Ks[:len(e.Utg.Ks)-kmerlen+1])...)
				}
			}

			// change strand
			if i < len(initPath)-1 {
				ne := edgesArr[initPath[i+1]]
				strand = GetNextEdgeStrand(e, ne, nodesArr, strand)
			}
		}
	} else { // direction == constructdbg.BACKWARD
		for i, eID := range initPath {
			e = edgesArr[eID]
			if i == 0 {
				if strand {
					initSeq = append(initSeq, constructdbg.GetReverseByteArr(e.Utg.Ks[:mi.EdgePos])...)
				} else {
					initSeq = append(initSeq, constructdbg.GetCompByteArr(e.Utg.Ks[mi.EdgePos:])...)
				}
			} else {
				if strand {
					initSeq = append(initSeq, constructdbg.GetReverseByteArr(e.Utg.Ks[:len(e.Utg.Ks)-(kmerlen-1)])...)
				} else {
					initSeq = append(initSeq, constructdbg.GetCompByteArr(e.Utg.Ks[kmerlen-1:])...)
				}
			}

			// change strand
			if i < len(initPath)-1 {
				ne := edgesArr[initPath[i+1]]
				strand = GetNextEdgeStrand(e, ne, nodesArr, strand)
			}
		}
	}

	if e.ID == 0 {
		log.Fatalf("[GetPathSeqArr] must used start edge info, e.ID: %v\n", e.ID)
	}
	sharePathLen = len(initSeq)
	// Get extendPathArr seq
	for i, path := range extendPathArr {
		strd := strand
		le := e
		pSeq := make([]byte, len(initSeq))
		copy(pSeq, initSeq)
		var edgeEndPos int
		var seqLen int
		//fmt.Printf("[GetPathSeqArr]path: %v,e.ID: %v, direction:%v\n", path, e.ID, direction)
		for j, eID := range path {
			ne := edgesArr[eID]
			// set strd
			strd = GetNextEdgeStrand(le, ne, nodesArr, strd)
			le = ne
			//_, ea, strArr := GetNextEdgeInfo(direction, e, strd, nodesArr, edgesArr)
			// add seq
			//fmt.Printf("[GetPathSeqArr]ne: %v, strand:%v, len(ne.Utg.Ks): %v, extLen: %v, seqLen: %v\n", ne.ID, strd, len(ne.Utg.Ks), extLen, seqLen)
			if direction == constructdbg.FORWARD {
				if j == len(path)-1 {
					addLen := extLen - seqLen
					if strd {
						edgeEndPos = kmerlen - 1 + addLen
						pSeq = append(pSeq, ne.Utg.Ks[kmerlen-1:edgeEndPos]...)
					} else {
						edgeEndPos = len(ne.Utg.Ks) - (kmerlen - 1) - addLen
						pSeq = append(pSeq, constructdbg.GetReverseCompByteArr(ne.Utg.Ks[edgeEndPos:len(ne.Utg.Ks)-(kmerlen-1)])...)
					}
				} else {
					if strd {
						pSeq = append(pSeq, ne.Utg.Ks[kmerlen-1:]...)
					} else {
						pSeq = append(pSeq, constructdbg.GetReverseCompByteArr(ne.Utg.Ks[:len(ne.Utg.Ks)-kmerlen+1])...)
					}
					seqLen += len(ne.Utg.Ks) - (kmerlen - 1)
				}
			} else { // direction == constructdbg.BACKWARD
				if j == len(path)-1 {
					addLen := extLen - seqLen
					if strd {
						edgeEndPos = len(ne.Utg.Ks) - (kmerlen - 1) - addLen
						pSeq = append(pSeq, constructdbg.GetReverseByteArr(ne.Utg.Ks[edgeEndPos:len(ne.Utg.Ks)-(kmerlen-1)])...)
					} else {
						edgeEndPos = (kmerlen - 1) + addLen
						pSeq = append(pSeq, constructdbg.GetCompByteArr(ne.Utg.Ks[kmerlen-1:edgeEndPos])...)
					}
				} else {
					if strd {
						pSeq = append(pSeq, constructdbg.GetReverseByteArr(ne.Utg.Ks[:len(ne.Utg.Ks)-(kmerlen-1)])...)
					} else {
						pSeq = append(pSeq, constructdbg.GetCompByteArr(ne.Utg.Ks[kmerlen-1:])...)
					}
					seqLen += len(ne.Utg.Ks) - (kmerlen - 1)
				}
			}
		}
		// add to array
		seqArr[i] = pSeq
		edgeEndPosArr[i] = edgeEndPos
		strandArr[i] = strd
	}

	return seqArr, edgeEndPosArr, strandArr, sharePathLen
}

func TraceMappingEdgePos(pa []constructdbg.DBG_MAX_INT, strand bool, edgeEndPos int, direction uint8, notMappingEdgeLen int, edgesArr []constructdbg.DBGEdge, nodesArr []constructdbg.DBGNode, kmerlen int) ([]constructdbg.DBG_MAX_INT, int, bool, int) {
	var edgePos, backLen int
	strd := strand
	//fmt.Printf("[TraceMappingEdgePos]pa: %v, strand: %v, edgeEndPos: %v, direction: %v, notMappingEdgeLen: %v\n", pa, strand, edgeEndPos, direction, notMappingEdgeLen)
	y := len(pa) - 1
	if direction == constructdbg.FORWARD {
		for ; y >= 0; y-- {
			e := edgesArr[pa[y]]
			if y == len(pa)-1 {
				if strd {
					if edgeEndPos-(kmerlen-1) > notMappingEdgeLen {
						edgePos = edgeEndPos - notMappingEdgeLen
						notMappingEdgeLen = 0
						break
					} else {
						notMappingEdgeLen -= (edgeEndPos - (kmerlen - 1))
					}
				} else {
					if len(e.Utg.Ks)-(kmerlen-1)-edgeEndPos > notMappingEdgeLen {
						edgePos = edgeEndPos + notMappingEdgeLen
						notMappingEdgeLen = 0
						break
					} else {
						notMappingEdgeLen -= (len(e.Utg.Ks) - (kmerlen - 1) - edgeEndPos)
					}
				}
			} else {
				if strd {
					if len(e.Utg.Ks)-(kmerlen-1) > notMappingEdgeLen {
						edgePos = len(e.Utg.Ks) - notMappingEdgeLen
						notMappingEdgeLen = 0
						break
					} else {
						notMappingEdgeLen -= (len(e.Utg.Ks) - (kmerlen - 1))
					}
				} else {
					if len(e.Utg.Ks)-(kmerlen-1) > notMappingEdgeLen {
						edgePos = (kmerlen - 1) + notMappingEdgeLen
						notMappingEdgeLen = 0
						break
					} else {
						notMappingEdgeLen -= (len(e.Utg.Ks) - (kmerlen - 1))
					}
				}
			}

			// change strd
			if y > 0 {
				ne := edgesArr[pa[y-1]]
				strd = GetNextEdgeStrand(e, ne, nodesArr, strd)
			}
		}
		if notMappingEdgeLen > 0 {
			log.Fatalf("[TraceMappingEdgePos] trace over pa: %v, but notMappingEdgeLen: %v, please check if pa is correct", pa, notMappingEdgeLen)
		}
		// if  last edge mapping length too short will be produce error path
		e := edgesArr[pa[y]]
		if strd {
			if edgePos-(kmerlen-1) < constructdbg.Min(kmerlen, (len(e.Utg.Ks)-(kmerlen-1)+10)*4/5) {
				if y <= 0 {
					pa = nil
					return pa, -1, false, backLen
				}
				backLen = edgePos - (kmerlen - 1)
				y--
				ne := edgesArr[pa[y]]
				strd = GetNextEdgeStrand(e, ne, nodesArr, strd)
				if strd {
					edgePos = len(ne.Utg.Ks) - 1
				} else {
					edgePos = kmerlen
				}
			}
		} else {
			if len(e.Utg.Ks)-(kmerlen-1)-edgePos < constructdbg.Min(kmerlen, (len(e.Utg.Ks)-(kmerlen-1)+10)*4/5) {
				if y <= 0 {
					pa = nil
					return pa, -1, false, backLen
				}
				backLen = len(e.Utg.Ks) - (kmerlen - 1) - edgePos
				y--
				ne := edgesArr[pa[y]]
				strd = GetNextEdgeStrand(e, ne, nodesArr, strd)
				if strd {
					edgePos = len(ne.Utg.Ks) - 1
				} else {
					edgePos = kmerlen
				}
			}
		}
		pa = pa[:y+1]
	} else { // direction == constructdbg.BACKWARD
		for ; y >= 0; y-- {
			e := edgesArr[pa[y]]
			if y == len(pa)-1 {
				if strd {
					if len(e.Utg.Ks)-(kmerlen-1)-edgeEndPos > notMappingEdgeLen {
						edgePos = edgeEndPos + notMappingEdgeLen
						notMappingEdgeLen = 0
						break
					} else {
						notMappingEdgeLen -= (len(e.Utg.Ks) - (kmerlen - 1) - edgeEndPos)
					}
				} else {
					if edgeEndPos-(kmerlen-1) > notMappingEdgeLen {
						edgePos = edgeEndPos - notMappingEdgeLen
						notMappingEdgeLen = 0
						break
					} else {
						notMappingEdgeLen -= (edgeEndPos - (kmerlen - 1))
					}
				}
			} else {
				if strd {
					if len(e.Utg.Ks)-(kmerlen-1) > notMappingEdgeLen {
						edgePos = (kmerlen - 1) + notMappingEdgeLen
						notMappingEdgeLen = 0
						break
					} else {
						notMappingEdgeLen -= (len(e.Utg.Ks) - (kmerlen - 1))
					}
				} else {
					if len(e.Utg.Ks)-(kmerlen-1) > notMappingEdgeLen {
						edgePos = len(e.Utg.Ks) - (kmerlen - 1) - notMappingEdgeLen
						notMappingEdgeLen = 0
						break
					} else {
						notMappingEdgeLen -= (len(e.Utg.Ks) - (kmerlen - 1))
					}
				}
			}

			// change strd
			if y > 0 {
				ne := edgesArr[pa[y-1]]
				strd = GetNextEdgeStrand(e, ne, nodesArr, strd)
			}
		}
		if notMappingEdgeLen > 0 {
			log.Fatalf("[TraceMappingEdgePos] trace over pa: %v, but notMappingEdgeLen: %v, please check if pa is correct", pa, notMappingEdgeLen)
		}
		// if  last edge mapping length too short will be produce error path
		e := edgesArr[pa[y]]
		if strd {
			if len(e.Utg.Ks)-(kmerlen-1)-edgePos < constructdbg.Min(kmerlen, (len(e.Utg.Ks)-(kmerlen-1)+10)*4/5) {
				if y <= 0 {
					pa = nil
					return pa, -1, false, backLen
				}
				backLen = len(e.Utg.Ks) - (kmerlen - 1) - edgePos
				y--
				ne := edgesArr[pa[y]]
				strd = GetNextEdgeStrand(e, ne, nodesArr, strd)
				if strd {
					edgePos = kmerlen
				} else {
					edgePos = len(ne.Utg.Ks) - 1
				}
			}
		} else {
			if edgePos-(kmerlen-1) < constructdbg.Min(kmerlen, (len(e.Utg.Ks)-(kmerlen-1)+10)*4/5) {
				if y <= 0 {
					pa = nil
					return pa, -1, false, backLen
				}
				backLen = edgePos - (kmerlen - 1)
				y--
				ne := edgesArr[pa[y]]
				strd = GetNextEdgeStrand(e, ne, nodesArr, strd)
				if strd {
					edgePos = kmerlen
				} else {
					edgePos = len(ne.Utg.Ks) - 1
				}
			}
		}
		pa = pa[:y+1]
	}
	//fmt.Printf("[TraceMappingEdgePos]y: %v, strd: %v, edgePos: %v, notMappingEdgeLen: %v\n", y, strd, edgePos, notMappingEdgeLen)

	return pa, edgePos, strd, backLen
}

func GetNextEdgeStrand(e, ne constructdbg.DBGEdge, nodesArr []constructdbg.DBGNode, strand bool) bool {
	var nd constructdbg.DBGNode
	if (e.StartNID == ne.StartNID) || (e.StartNID == ne.EndNID) {
		nd = nodesArr[e.StartNID]
	} else {
		nd = nodesArr[e.EndNID]
	}
	if e.StartNID == e.EndNID {
		if ne.StartNID != ne.EndNID {
			if constructdbg.CountEdgeIDComing(nd.EdgeIDOutcoming, e.ID) == 2 {
				if constructdbg.IsInComing(nd.EdgeIDIncoming, ne.ID) && (ne.EndNID == e.EndNID) {
					strand = !strand
				}
			} else if constructdbg.CountEdgeIDComing(nd.EdgeIDIncoming, e.ID) == 2 {
				if constructdbg.IsInComing(nd.EdgeIDOutcoming, ne.ID) && (ne.EndNID == e.EndNID) {
					strand = !strand
				}
			} else {
				if (constructdbg.IsInComing(nd.EdgeIDIncoming, ne.ID) && ne.StartNID == e.StartNID) || (constructdbg.IsInComing(nd.EdgeIDOutcoming, ne.ID) && ne.EndNID == e.EndNID) {
					strand = !strand
				}
			}
		}
	} else {
		if ne.StartNID == ne.EndNID {
			if constructdbg.CountEdgeIDComing(nd.EdgeIDOutcoming, ne.ID) == 2 {
				if constructdbg.IsInComing(nd.EdgeIDIncoming, e.ID) && ne.StartNID == e.StartNID {
					strand = !strand
				}
			} else if constructdbg.CountEdgeIDComing(nd.EdgeIDIncoming, ne.ID) == 2 {
				if constructdbg.IsInComing(nd.EdgeIDOutcoming, e.ID) && ne.StartNID == e.StartNID {
					strand = !strand
				}
			} else if constructdbg.IsInComing(nd.EdgeIDIncoming, e.ID) && (e.StartNID == ne.StartNID) || (constructdbg.IsInComing(nd.EdgeIDOutcoming, e.ID) && e.EndNID == ne.EndNID) {
				strand = !strand
			}
		} else {
			if (e.StartNID == ne.StartNID) || (e.EndNID == ne.EndNID) {
				strand = !strand
			}
		}
	}

	return strand
}

func GetSelfCycleEdgeNextEdgesInfo(direction uint8, e constructdbg.DBGEdge, strand bool, nodesArr []constructdbg.DBGNode, edgesArr []constructdbg.DBGEdge) (nID constructdbg.DBG_MAX_INT, ea []constructdbg.DBG_MAX_INT, strandArr []bool) {
	nID = e.EndNID
	nd := nodesArr[nID]
	if constructdbg.CountEdgeIDComing(nd.EdgeIDIncoming, e.ID) == 2 {
		for i := 0; i < bnt.BaseTypeNum; i++ {
			id := nd.EdgeIDOutcoming[i]
			if id > 1 {
				ea = append(ea, id)
				if edgesArr[id].StartNID != edgesArr[id].EndNID && edgesArr[id].EndNID == nID {
					strandArr = append(strandArr, !strand)
				} else {
					strandArr = append(strandArr, strand)
				}
			}
		}
	} else if constructdbg.CountEdgeIDComing(nd.EdgeIDOutcoming, e.ID) == 2 {
		for i := 0; i < bnt.BaseTypeNum; i++ {
			id := nd.EdgeIDIncoming[i]
			if id > 1 {
				ea = append(ea, id)
				if edgesArr[id].StartNID != edgesArr[id].EndNID && edgesArr[id].StartNID == nID {
					strandArr = append(strandArr, !strand)
				} else {
					strandArr = append(strandArr, strand)
				}
			}
		}
	} else { // eID contain two EdgeIDComing
		if (direction == constructdbg.FORWARD && strand == constructdbg.PLUS) || (direction == constructdbg.BACKWARD && strand == constructdbg.MINUS) {
			for i := 0; i < bnt.BaseTypeNum; i++ {
				id := nd.EdgeIDOutcoming[i]
				if id > 1 {
					ea = append(ea, id)
					if edgesArr[id].StartNID != edgesArr[id].EndNID && edgesArr[id].EndNID == nID {
						strandArr = append(strandArr, !strand)
					} else {
						strandArr = append(strandArr, strand)
					}
				}
			}
		} else {
			for i := 0; i < bnt.BaseTypeNum; i++ {
				id := nd.EdgeIDIncoming[i]
				if id > 1 {
					ea = append(ea, id)
					if edgesArr[id].StartNID != edgesArr[id].EndNID && edgesArr[id].StartNID == nID {
						strandArr = append(strandArr, !strand)
					} else {
						strandArr = append(strandArr, strand)
					}
				}
			}
		}
	}
	return
}

func GetNextEdgeInfo(direction uint8, e constructdbg.DBGEdge, strand bool, nodesArr []constructdbg.DBGNode, edgesArr []constructdbg.DBGEdge) (nID constructdbg.DBG_MAX_INT, ea []constructdbg.DBG_MAX_INT, strandArr []bool) {
	if e.StartNID == e.EndNID {
		return GetSelfCycleEdgeNextEdgesInfo(direction, e, strand, nodesArr, edgesArr)
	}

	if direction == constructdbg.FORWARD {
		if strand == constructdbg.PLUS {
			nID = e.EndNID
		} else {
			nID = e.StartNID
		}
	} else {
		if strand == constructdbg.PLUS {
			nID = e.StartNID
		} else {
			nID = e.EndNID
		}
	}

	nd := nodesArr[nID]
	if constructdbg.IsInComing(nd.EdgeIDIncoming, e.ID) {
		for i := 0; i < bnt.BaseTypeNum; i++ {
			id := nd.EdgeIDOutcoming[i]
			if id > 1 {
				ea = append(ea, id)
				ne := edgesArr[id]
				if e.EndNID == nID {
					if ne.StartNID != ne.EndNID && ne.EndNID == nID {
						strandArr = append(strandArr, !strand)
					} else {
						strandArr = append(strandArr, strand)
					}
				} else {
					if ne.StartNID == nID {
						strandArr = append(strandArr, !strand)
					} else {
						strandArr = append(strandArr, strand)
					}
				}
			}
		}
	} else {
		for i := 0; i < bnt.BaseTypeNum; i++ {
			id := nd.EdgeIDIncoming[i]
			if id > 1 {
				ea = append(ea, id)
				ne := edgesArr[id]
				if e.StartNID == nID {
					if ne.StartNID != ne.EndNID && ne.StartNID == nID {
						strandArr = append(strandArr, !strand)
					} else {
						strandArr = append(strandArr, strand)
					}
				} else {
					if ne.EndNID == nID {
						strandArr = append(strandArr, !strand)
					} else {
						strandArr = append(strandArr, strand)
					}
				}
			}
		}
	}
	return
}

func ResetDeduceReadPos(mi MapInfo, e constructdbg.DBGEdge, direction uint8) (deduceReadPos int) {
	deduceReadPos = mi.ReadPos
	if direction == constructdbg.FORWARD {
		if mi.Strand {
			deduceReadPos += len(e.Utg.Ks) - mi.EdgePos
		} else {
			deduceReadPos += mi.EdgePos
		}
	} else {
		if mi.Strand {
			deduceReadPos -= mi.EdgePos
		} else {
			deduceReadPos -= (len(e.Utg.Ks) - mi.EdgePos)
		}
	}
	return deduceReadPos
}

func IsBubble(e1, e2 constructdbg.DBGEdge, nodesArr []constructdbg.DBGNode) bool {
	if (e1.StartNID == e1.EndNID) || (e2.StartNID == e2.EndNID) || (e1.StartNID < 2) || (e1.EndNID < 2) || (e2.StartNID < 2) || (e2.EndNID < 2) {
		return false
	}
	// test startNID
	if constructdbg.IsInComing(nodesArr[e1.StartNID].EdgeIDIncoming, e1.ID) {
		if !constructdbg.IsInComing(nodesArr[e1.StartNID].EdgeIDIncoming, e2.ID) {
			return false
		}
	} else {
		if !constructdbg.IsInComing(nodesArr[e1.StartNID].EdgeIDOutcoming, e2.ID) {
			return false
		}
	}
	// test endNID
	if constructdbg.IsInComing(nodesArr[e1.EndNID].EdgeIDIncoming, e1.ID) {
		if !constructdbg.IsInComing(nodesArr[e1.EndNID].EdgeIDIncoming, e2.ID) {
			return false
		}
	} else {
		if !constructdbg.IsInComing(nodesArr[e1.EndNID].EdgeIDOutcoming, e2.ID) {
			return false
		}
	}

	return true
}

func IsBubbleEdge(e constructdbg.DBGEdge, nodesArr []constructdbg.DBGNode) bool {
	if (e.StartNID < 2) || (e.EndNID < 2) {
		return false
	}
	var ea1, ea2 []constructdbg.DBG_MAX_INT

	if constructdbg.IsInComing(nodesArr[e.StartNID].EdgeIDIncoming, e.ID) {
		ea1 = constructdbg.GetComingOtherEArr(nodesArr[e.StartNID].EdgeIDIncoming, e.ID)
	} else {
		ea1 = constructdbg.GetComingOtherEArr(nodesArr[e.StartNID].EdgeIDOutcoming, e.ID)
	}
	if constructdbg.IsInComing(nodesArr[e.EndNID].EdgeIDIncoming, e.ID) {
		ea2 = constructdbg.GetComingOtherEArr(nodesArr[e.EndNID].EdgeIDIncoming, e.ID)
	} else {
		ea2 = constructdbg.GetComingOtherEArr(nodesArr[e.EndNID].EdgeIDOutcoming, e.ID)
	}

	arr := constructdbg.InterSecDBG_MAX_INTArr(ea1, ea2)
	if len(arr) > 0 {
		return true
	} else {
		return false
	}

}

func EqualDBG_MAX_INTArr(arr1, arr2 []constructdbg.DBG_MAX_INT) bool {
	if len(arr1) != len(arr2) {
		return false
	}
	for i := 0; i < len(arr1); i++ {
		if arr1[i] != arr2[i] {
			return false
		}
	}
	return true
}

func MappingONTMostProbablePath(arr []LRRecord, edgesArr []constructdbg.DBGEdge, nodesArr []constructdbg.DBGNode, opt Options, flankAllow int) (extArr [2][]constructdbg.DBG_MAX_INT) {
	// check not select cycle edge
	pos := findSeedEID(arr, edgesArr, flankAllow)
	if pos < 0 || arr[pos].MapNum < opt.Kmer*2 {
		return
	}
	SeedRecord := arr[pos]
	//fmt.Printf("[MappingONTMostProbablePath] found seed Record: %v\n", SeedRecord)
	//fmt.Printf("[MappingONTMostProbablePath] len(ReadSeqBnt): %v\n", len(SeedRecord.ReadSeqBnt))

	extendMappingMinLen := 2000
	maxFlank := opt.Kmer - 2
	extArr[0] = append(extArr[0], SeedRecord.RefID)
	extArr[1] = append(extArr[1], 0)
	for i := 0; i < 2; i++ { // i == 0 note BACKWARD, i == 1 note FORWARD
		var flankLen int
		var mi MapInfo
		mi.EIDIdx = len(extArr[0]) - 1
		var direction uint8
		var deduceReadPos, extIdx int
		if i == 0 {
			if SeedRecord.Strand == constructdbg.PLUS {
				flankLen = SeedRecord.RefStart
				mi.EdgePos = SeedRecord.RefStart
			} else {
				flankLen = SeedRecord.RefLen - SeedRecord.RefEnd
				mi.EdgePos = SeedRecord.RefEnd
			}
			direction = constructdbg.BACKWARD
			mi.ReadPos = SeedRecord.Start
			deduceReadPos = mi.ReadPos - flankLen
		} else {
			if SeedRecord.Strand == constructdbg.PLUS {
				flankLen = SeedRecord.RefLen - SeedRecord.RefEnd
				mi.EdgePos = SeedRecord.RefEnd
			} else {
				flankLen = SeedRecord.RefStart
				mi.EdgePos = SeedRecord.RefStart
			}
			direction = constructdbg.FORWARD
			mi.ReadPos = SeedRecord.End
			deduceReadPos = mi.ReadPos + flankLen
		}

		mi.Strand = SeedRecord.Strand
		//mi.NID = nID
		extIdx = len(extArr[0])
		strand := SeedRecord.Strand
		for {
			//var alter constructdbg.DBG_MAX_INT // alternative edge, maybe cause by haplotype
			//fmt.Printf("[MappingMostProbablePath]mi: %v, flankLen: %v, deduceReadPos: %v\n\textArr: %v,  e.ID: %v, nID: %v\n", mi, flankLen, deduceReadPos, extArr, e.ID, nID)
			if flankLen >= maxFlank {
				break
			}
			e := edgesArr[extArr[0][len(extArr[0])-1]]
			//var nID constructdbg.DBG_MAX_INT
			//fmt.Printf("[MappingONTMostProbablePath]extArr: %v, e.ID: %v\n", extArr, e.ID)
			nID, ea, strArr := GetNextEdgeInfo(direction, e, strand, nodesArr, edgesArr)
			//fmt.Printf("[MappingONTMostProbablePath]next nID: %v, ea: %v, strandArr: %v\n", nID, ea, strArr)
			if nID == 0 || len(ea) == 0 {
				break
			}

			// test if readPos out of range
			for ; extIdx < len(extArr[0]); extIdx++ {
				leID := extArr[0][extIdx]
				if direction == constructdbg.FORWARD {
					deduceReadPos += (len(edgesArr[leID].Utg.Ks) - (opt.Kmer - 1))
				} else {
					deduceReadPos -= (len(edgesArr[leID].Utg.Ks) - (opt.Kmer - 1))
				}
			}
			if direction == constructdbg.FORWARD {
				if deduceReadPos >= SeedRecord.Len-opt.Kmer {
					break
				}
			} else {
				if deduceReadPos <= opt.Kmer {
					break
				}
			}

			//ea = constructdbg.GetNearEdgeIDArr(nodesArr[nID], e.ID, coming)
			//fmt.Printf("[MappingONTMostProbablePath]deduceReadPos: %v\n", deduceReadPos)
			if len(ea) == 1 {
				extArr[0] = append(extArr[0], ea[0])
				extArr[1] = append(extArr[1], 0)
				//ne := edgesArr[ea[0]]
				// reset strand
				strand = strArr[0]
				continue
			}

			{
				// case if len(ea) > 1, need long read extend mapping
				// set extMinLen
				extMinLen := extendMappingMinLen
				if direction == constructdbg.FORWARD {
					if deduceReadPos >= SeedRecord.Len-opt.Kmer*2 {
						break
					}
					if SeedRecord.Len-deduceReadPos < extendMappingMinLen {
						extMinLen = SeedRecord.Len - deduceReadPos
					}
				} else {
					if deduceReadPos <= opt.Kmer*2 {
						break
					}
					if deduceReadPos < extendMappingMinLen {
						extMinLen = deduceReadPos
					}
				}

				extendPathArr := GetExtendMappingPathArr(ea, strArr, direction, extMinLen, nodesArr, edgesArr, opt.Kmer)
				if len(extendPathArr) > 500 {
					fmt.Printf("[MappingONTMostProbablePath] need mapping path array is large >500, len(extendPathArr): %v\n", len(extendPathArr))
				} else if len(extendPathArr) == 0 {
					break
				}
				//fmt.Printf("[MappingONTMostProbablePath]mi: %v,opt.Kmer:%v, extMinLen: %v, extendPathArr: %v\n", mi, opt.Kmer, extMinLen, extendPathArr)
				pathSeqArr, edgeEndPosArr, strandArr, sharePathLen := GetPathSeqArr(mi, extArr[0][mi.EIDIdx:], direction, extendPathArr, edgesArr, nodesArr, opt.Kmer, extMinLen)
				extMinLen += sharePathLen
				if len(pathSeqArr) == 0 {
					break
				}
				//fmt.Printf("[MappingONTMostProbablePath]sharePathLen: %v, edgeEndPosArr: %v\n", sharePathLen, edgeEndPosArr)
				//fmt.Printf("[MappingONTMostProbablePath]strandArr: %v\n", strandArr)
				/*for m := 0; m < len(pathSeqArr); m++ {
					fmt.Printf("[MappingONTMostProbablePath]length of pathSeqArr[%v]: %v\n", m, len(pathSeqArr[m]))
				}*/
				var readSeqBnt []byte
				if direction == constructdbg.FORWARD {
					upper := mi.ReadPos + extMinLen
					if upper > SeedRecord.Len {
						upper = SeedRecord.Len
					}
					readSeqBnt = SeedRecord.ReadSeqBnt[mi.ReadPos:upper]
				} else {
					floor := mi.ReadPos - extMinLen
					if mi.ReadPos-extMinLen < floor {
						floor = 0
					}
					readSeqBnt = constructdbg.GetReverseByteArr(SeedRecord.ReadSeqBnt[floor:mi.ReadPos])
				}
				// mapping long read region
				cgArr := make([]CIGAR, len(pathSeqArr))
				chainScoreArr := make([]int, len(pathSeqArr))
				alignEdgePosArr := make([]int, len(pathSeqArr))
				alignReadPosArr := make([]int, len(pathSeqArr))
				for x, ps := range pathSeqArr {
					cgArr[x], chainScoreArr[x], alignEdgePosArr[x], alignReadPosArr[x] = DPAlignEdgePathExtend(ps, readSeqBnt, edgesArr, nodesArr)
				}

				// get align score, denote Match=1, Insert=-1, Delete=-1, MisMatch=-1
				sc := make([]int, len(pathSeqArr))
				for x, cg := range cgArr {
					sc[x] = int(cg.Mch)
				}
				fmt.Printf("[MappingONTMostProbablePath] cgArr: %v\n", cgArr)
				fmt.Printf("[MappingONTMostProbablePath]   chainScoreArr: %v\n", chainScoreArr)
				fmt.Printf("[MappingONTMostProbablePath]           score: %v\n", sc)
				fmt.Printf("[MappingONTMostProbablePath] alignEdgePosArr: %v\n", alignEdgePosArr)
				fmt.Printf("[MappingONTMostProbablePath] alignReadPosArr: %v\n", alignReadPosArr)
				sort.Ints(sc)
				max := sc[len(sc)-1]
				/*if max < extMinLen/2 {
					fmt.Printf("[MappingONTMostProbablePath] mapping score too low, ea: %v, cgArr: %v, <extMinLen/2: %v\n", ea, cgArr, extMinLen/2)
					break
				}*/
				// trace BACKWARD
				backLenArr := make([]int, len(extendPathArr))
				{
					notFoundNextPath := false
					for m := 0; m < len(extendPathArr); m++ {
						if int(cgArr[m].Mch) < max {
							continue
						}
						if alignEdgePosArr[m] < sharePathLen+opt.Kmer {
							notFoundNextPath = true
							break
						}
						notMappingEdgeLen := extMinLen - alignEdgePosArr[m]
						extendPathArr[m], edgeEndPosArr[m], strandArr[m], backLenArr[m] = TraceMappingEdgePos(extendPathArr[m], strandArr[m], edgeEndPosArr[m], direction, notMappingEdgeLen, edgesArr, nodesArr, opt.Kmer)
					}
					if notFoundNextPath {
						break
					}
				}
				if len(sc) == 1 || (len(sc) >= 2 && sc[len(sc)-1] > sc[len(sc)-2]) {
					x := IndexCigarArr(cgArr, max)
					pa := extendPathArr[x]
					if len(pa) <= 0 {
						break
					}
					if max < alignEdgePosArr[x]/2 {
						fmt.Printf("[MappingONTMostProbablePath] mapping score too low, cgArr: %v, < alignEdgePosArr[%v]: %v\n", cgArr[x], x, alignEdgePosArr[x]/2)
						break
					}
					// reset mi
					if direction == constructdbg.FORWARD {
						mi.ReadPos += (alignReadPosArr[x] - backLenArr[x])
					} else {
						mi.ReadPos -= (alignReadPosArr[x] - backLenArr[x])
					}
					//notMappingEdgeLen := extMinLen - alignEdgePosArr[x]
					//pa, mi.EdgePos, mi.Strand = TraceMappingEdgePos(pa, strandArr[x], edgeEndPosArr[x], direction, notMappingEdgeLen, edgesArr, nodesArr, opt.Kmer)

					for z := 0; z < len(pa); z++ {
						extArr[0] = append(extArr[0], pa[z])
						extArr[1] = append(extArr[1], 0)
					}
					mi.EIDIdx = len(extArr[0]) - 1
					mi.EdgePos = edgeEndPosArr[x]
					mi.Strand = strandArr[x]
					// set deduceReadPos
					deduceReadPos = ResetDeduceReadPos(mi, edgesArr[pa[len(pa)-1]], direction)
					extIdx = mi.EIDIdx
					strand = mi.Strand
				} else if len(sc) > 1 && sc[len(sc)-1] == sc[len(sc)-2] {
					samePath := true
					x1, x2 := -1, -1
					for m := 0; m < len(extendPathArr); m++ {
						if int(cgArr[m].Mch) != max {
							continue
						}
						if len(extendPathArr[m]) <= 0 {
							continue
						}
						if x1 >= 0 {
							if !EqualDBG_MAX_INTArr(extendPathArr[x1], extendPathArr[m]) {
								samePath = false
								x2 = m
							}
						} else {
							x1 = m
						}
					}
					if x1 >= 0 && max < alignEdgePosArr[x1]/2 {
						fmt.Printf("[MappingONTMostProbablePath] mapping score too low, cgArr: %v, < alignEdgePosArr[%v]: %v\n", cgArr[x1], x1, alignEdgePosArr[x1]/2)
						break
					}
					if x1 >= 0 && samePath {
						pa := extendPathArr[x1]
						// reset mi
						if direction == constructdbg.FORWARD {
							mi.ReadPos += (alignReadPosArr[x1] - backLenArr[x1])
						} else {
							mi.ReadPos -= (alignReadPosArr[x1] - backLenArr[x1])
						}
						for z := 0; z < len(pa); z++ {
							extArr[0] = append(extArr[0], pa[z])
							extArr[1] = append(extArr[1], 0)
						}
						mi.EIDIdx = len(extArr[0]) - 1
						mi.EdgePos = edgeEndPosArr[x1]
						mi.Strand = strandArr[x1]
						strand = mi.Strand
						// set deduceReadPos
						deduceReadPos = ResetDeduceReadPos(mi, edgesArr[pa[len(pa)-1]], direction)
						extIdx = mi.EIDIdx
					} else if x1 >= 0 && x2 >= 0 {
						if len(sc) > 2 {
							if sc[len(sc)-2] == sc[len(sc)-3] {
								break
							}
						}
						// maybe is bubble
						pa1, pa2 := extendPathArr[x1], extendPathArr[x2]
						if len(pa1) == 0 || len(pa2) == 0 || len(pa1) != len(pa2) {
							break
						}
						y := 0
						for ; y < len(pa1); y++ {
							if pa1[y] != pa2[y] {
								e1, e2 := edgesArr[pa1[y]], edgesArr[pa2[y]]
								if IsBubble(e1, e2, nodesArr) {

								} else {
									break
								}
							}
						}
						if y < len(pa1) {
							break
						}
						// reset mi
						if direction == constructdbg.FORWARD {
							mi.ReadPos += (alignReadPosArr[x1] - backLenArr[x1])
						} else {
							mi.ReadPos -= (alignReadPosArr[x1] - backLenArr[x1])
						}

						for z := 0; z < len(pa1); z++ {
							extArr[0] = append(extArr[0], pa1[z])
							if pa1[z] != pa2[z] {
								extArr[1] = append(extArr[1], pa2[z])
							} else {
								extArr[1] = append(extArr[1], 0)
							}
						}
						mi.EIDIdx = len(extArr[0]) - 1
						mi.EdgePos = edgeEndPosArr[x1]
						mi.Strand = strandArr[x1]
						strand = mi.Strand
						// set deduceReadPos
						deduceReadPos = ResetDeduceReadPos(mi, edgesArr[pa1[len(pa1)-1]], direction)
						extIdx = mi.EIDIdx
					} else {
						break
					}
				} else {
					break
				}
			}
		}

		// extArr process
		if direction == constructdbg.BACKWARD {
			constructdbg.ReverseDBG_MAX_INTArr(extArr[0])
			constructdbg.ReverseDBG_MAX_INTArr(extArr[1])
		}
	}

	return
}

/*{
			// find if can help by NGS mapping info
			eIDArr := GetNGSPath(extArr, edgesArr)
			if len(eIDArr) > 0 {
				for _, eID := range eIDArr {
					extArr[0] = append(extArr[0], eID)
					extArr[1] = append(extArr[1], alter)
					e = edgesArr[eID]
					if e.StartNID == nID {
						nID = e.EndNID
					} else {
						nID = e.StartNID
					}
				}
				fmt.Printf("[MappingMostProbablePath] extend path by NGS Path: %v\n", eIDArr)
				continue
			}
			fmt.Printf("[MappingMostProbablePath] extend path: %v\n", extArr)
			// mapping the ONT to the edges Path
			//minLen := GetMinLen(ea, edgesArr, opt.Kmer)
			var mappingPathArr [][]constructdbg.DBG_MAX_INT
			var lastEdgeArr []constructdbg.DBG_MAX_INT
			MinExtendLen := opt.Kmer / 2
			if minLen < MinExtendLen {
				var newml int
				mappingPathArr, lastEdgeArr, newml = GetMappingExtendPathArr(extArr[0], ea, nID, edgesArr, nodesArr, MinExtendLen, opt.Kmer)
				minLen = newml
			} else {
				for _, eID := range ea {
					var mp []constructdbg.DBG_MAX_INT
					mp = append(mp, extArr[0]...)
					mappingPathArr = append(mappingPathArr, mp)
					lastEdgeArr = append(lastEdgeArr, eID)
				}
			}
			extendLen := GetExtendLen(extArr[0], mi.EIDIdx+1, edgesArr, opt.Kmer)
			extendLen += minLen
			cgArr := make([]CIGAR, len(lastEdgeArr))
			chainScoreArr := make([]int, len(lastEdgeArr))
			notChainEdgeLenA := make([]int, len(lastEdgeArr))
			readPosA := make([]int, len(lastEdgeArr))
			edgesSeqLen := make([]int, len(lastEdgeArr))
			strandArr := make([]bool, len(lastEdgeArr))
			alignReadPosArr := make([]int, len(lastEdgeArr))
			fmt.Printf("[MappingMostProbablePath] mi: %v, direction: %v, minLen: %v, extendLen: %v\n", mi, direction, minLen, extendLen)
			for x, eID := range lastEdgeArr {
				fmt.Printf("[MappingMostProbablePath]**************** mappingPathArr[%v]:%v,lasteID: %v *****************\n", x, mappingPathArr[x], eID)
				cgArr[x], chainScoreArr[x], edgesSeqLen[x], notChainEdgeLenA[x], readPosA[x], strandArr[x], alignReadPosArr[x] = DPAlignEdgePath(mi, mappingPathArr[x], edgesArr, nodesArr, arr[pos].ReadSeqBnt, direction, eID, extendLen, opt.Kmer)
				fmt.Printf("[MappingMostProbablePath]################ cgArr[%v]: %v,edgesSeqLen[%v]: %v, notChainEdgeLenA[%v]: %v, readPosA[%v]: %v, alignReadPos[%v]: %v, strandArr[%v]: %v ###############\n", x, cgArr[x], x, edgesSeqLen[x], x, notChainEdgeLenA[x], x, readPosA[x], x, alignReadPosArr[x], x, strandArr[x])
				if notChainEdgeLenA[x] < 0 {
					log.Fatalf("[MappingMostProbablePath] notChainEdgeLenA[%v]: %v < 0\n", x, notChainEdgeLenA[x])
				}
			}

			sc := make([]int, len(lastEdgeArr))
			for x, cg := range cgArr {
				sc[x] = int(cg.Mch) - int(cg.Del) - int(cg.Ins) - int(cg.Mis)
			}
			sort.Ints(sc)
			x := IndexCigarArr(cgArr, sc[len(sc)-1])
			max := sc[len(sc)-1]
			if max < edgesSeqLen[x]*2/5 {
				fmt.Printf("[MappingMostProbablePath] mapping score too low, ea: %v, cgArr: %v, <edgeSeqLen[%v]/2\n", ea, cgArr, edgesSeqLen[x])
				break
			}
			if sc[len(sc)-1]-sc[len(sc)-2] < sc[len(sc)-1]/10 {
				// choose the readSeq mapping the most small one
				var diffArr []Diff
				for y, cg := range cgArr {
					if max-(int(cg.Mch)-int(cg.Del)-int(cg.Ins)-int(cg.Mis)) < max/10 {
						var diff Diff
						diff.Idx = y
						diff.Unidentity = (float64(cg.Del) + float64(cg.Ins) + float64(cg.Mis)) / float64(cg.Mch)
						diffArr = append(diffArr, diff)
					}
				}
				sort.Sort(DiffArr(diffArr))
				if diffArr[0].Unidentity < diffArr[1].Unidentity {
					sc[len(sc)-1] = math.MaxInt64
					x = diffArr[0].Idx
				}  else if (diffArr[0].Unidentity == diffArr[1].Unidentity) && (chainScoreArr[diffArr[0].Idx] != chainScoreArr[diffArr[1].Idx]) {
					if chainScoreArr[diffArr[0].Idx] > chainScoreArr[diffArr[1].Idx] {
						sc[len(sc)-1] = math.MaxInt64
						x = diffArr[0].Idx
					} else {
						sc[len(sc)-1] = math.MaxInt64
						x = diffArr[1].Idx
					}
				}
				readMappingMinLen := math.MaxInt64
				idx, count := -1, 0
				for y, s := range score {
					if s == max {
						if direction == constructdbg.BACKWARD {
							if mi.ReadPos-alignReadPosArr[y] < readMappingMinLen {
								readMappingMinLen = mi.ReadPos - alignReadPosArr[y]
								idx = y
								count = 1
							} else if mi.ReadPos-alignReadPosArr[y] == readMappingMinLen {
								count++
							}
						} else {
							if alignReadPosArr[y]-mi.ReadPos < readMappingMinLen {
								readMappingMinLen = alignReadPosArr[y] - mi.ReadPos
								idx = y
								count = 1
							} else if alignReadPosArr[y]-mi.ReadPos == readMappingMinLen {
								count++
							}
						}
					}
				}
				if count == 1 && idx >= 0 {
					score[idx] = math.MaxInt64
					sc[len(sc)-1] = math.MaxInt64
					x = idx
				}
				var idxArr []int
				for y, s := range score {
					if s == max {
						idxArr = append(idxArr, y)
					}
				}
				consisPath := GetLongConsisPath(mappingPathArr, idxArr)
				if len(consisPath) <= len(extArr[0]) {
					// test if a bubble edges
					if !(len(ea) == 2 && constructdbg.IsBubble(ea[0], ea[1], edgesArr)) {
						fmt.Printf("[MappingMostProbablePath] encounter not unsure case, ea: %v, score: %v", ea, score)
						break
					}
				} else {
					mappingPathArr[x] = consisPath[:len(consisPath)-1]
					lastEdgeArr[x] = consisPath[len(consisPath)-1]
					eID := lastEdgeArr[x]
					extendLen = GetExtendLen(mappingPathArr[x], mi.EIDIdx+1, edgesArr, opt.Kmer)
					extendLen += len(edgesArr[eID].Utg.Ks) - (opt.Kmer - 1)
					score[x], edgesSeqLen[x], notChainEdgeLenA[x], readPosA[x], strandArr[x] = DPAlignEdgePath(mi, mappingPathArr[x], edgesArr, nodesArr, arr[pos].ReadSeqBnt, direction, eID, extendLen, opt.Kmer)
					sc[len(sc)-1] = score[x]
					for y := 0; y < len(sc)-1; y++ {
						sc[y] = 0
					}
				}
			}
			eID := lastEdgeArr[x]
			e = edgesArr[eID]
			// process mi
			mi.Strand = strandArr[x]
			flankLen = notChainEdgeLenA[x]
			minLen = extendLen - GetExtendLen(mappingPathArr[x], mi.EIDIdx+1, edgesArr, opt.Kmer)
			mi.EdgePos = GetEdgePos(e, direction, mi.Strand, opt.Kmer, minLen, flankLen)

			mi.EIDIdx = len(mappingPathArr[x])
			mi.ReadPos = readPosA[x]
			if sc[len(sc)-1] > sc[len(sc)-2] {
				mappingPathArr[x] = append(mappingPathArr[x], lastEdgeArr[x])
				for y := len(extArr[0]); y < len(mappingPathArr[x]); y++ {
					e = edgesArr[mappingPathArr[x][y]]
					extArr[0] = append(extArr[0], e.ID)
					extArr[1] = append(extArr[1], alter)
					if e.StartNID == nID {
						nID = e.EndNID
					} else {
						nID = e.StartNID
					}
				}
				mi.NID = nID
				continue
			}
			if len(ea) == 2 && constructdbg.IsBubble(ea[0], ea[1], edgesArr) {
				alter = ea[1]
				eID := ea[0]
				extArr[0] = append(extArr[0], eID)
				extArr[1] = append(extArr[1], alter)
				e = edgesArr[eID]
				if e.StartNID == nID {
					nID = e.EndNID
				} else {
					nID = e.StartNID
				}
				mi.NID = nID
				continue
			}

			fmt.Printf("[MappingMostProbablePath] encounter not process case, ea: %v, cgArr: %v", ea, cgArr)
			break
		}

		// extArr process
		if direction == constructdbg.BACKWARD {
			constructdbg.ReverseDBG_MAX_INTArr(extArr[0])
			constructdbg.ReverseDBG_MAX_INTArr(extArr[1])
		}
	}
	return
} */

/*func findMostProbablePath(arr []LRRecord, edgesArr []constructdbg.DBGEdge, nodesArr []constructdbg.DBGNode, opt Options) (extArr [2][]constructdbg.DBG_MAX_INT) {
	// found seed edge ID
	MatchScore := 5
	pos := findSeedEID(arr, edgesArr)
	if pos < 0 || arr[pos].GapMapNum < opt.Kmer {
		return
	}
	fmt.Printf("[findMostProbablePath] found seed Record: %v\n", arr[pos])
	extArr = append(extArr, arr[pos].RefID)
	nID := edgesArr[arr[pos].RefID].StartNID
	if arr[pos].Strand {
		nID = edgesArr[arr[pos].RefID].EndNID
	}
	maxFlank := opt.Kmer * 4 / 5
	// extend start partition
	if nID > 0 && arr[pos].RefStart < maxFlank {
		e := edgesArr[arr[pos].RefID]
		nd := nodesArr[nID]
		si := pos
		seqLen := GetFlankLen(arr[si], constructdbg.BACKWARD)
		for {
			fmt.Printf("[findMostProbablePath]BACKWARD found e.ID: %v, nd.ID: %v, seqLen: %v, si: %v\n", e.ID, nd.ID, seqLen, si)
			var newArr []constructdbg.DBG_MAX_INT
			var puzzy bool // note if found a puzzy path, maybe Long Reads sequence error cause
			if len(e.PathMat) > 0 && len(e.PathMat[0].IDArr) > 0 {
				//var arrNGS []constructdbg.DBG_MAX_INT
				d := constructdbg.BACKWARD
				if e.EndNID == nd.ID {
					d = constructdbg.FORWARD
				}
				arrNGS := GetExtendPath(e.PathMat[0].IDArr, e.ID, d)
				if len(arrNGS) > 0 {
					//extArr = append(extArr, arrNGS...)
					var perr bool
					for _, eID := range arrNGS {
						if edgesArr[eID].ID == 0 || edgesArr[eID].GetDeleteFlag() > 0 {
							perr = true
						}
					}
					if !perr {
						newArr = append(newArr, arrNGS...)
					}
				}
			}

			if len(newArr) == 0 {
				ea := constructdbg.GetNearEdgeIDArr(nd, e.ID)
				if len(ea) == 1 {
					//extArr = append(extArr, ea[0])
					newArr = append(newArr, ea[0])
				} else if len(ea) > 1 {
					var a []LRRecord
					for _, eID := range ea {
						idx, ok := SearchLRRecordArr(arr, si, constructdbg.BACKWARD, eID, (seqLen+edgesArr[eID].GetSeqLen()-(opt.Kmer-1))*6/5)
						if ok && arr[idx].RefID == eID {
							a = append(a, arr[idx])
						}
					}

					minLen := GetMinLenEID(edgesArr, ea)
					if len(a) == 1 {
						if minLen > 2*opt.Kmer {
							if a[0].RefStart < a[0].Mch*1/5 && a[0].Mch > a[0].RefLen*4/5 { // if the long reads start position
								newArr = append(newArr, a[0].RefID)
							}
						}
					} else if len(a) == len(ea) {
						sort.Sort(LRRecordScoreArr(a))
						diffScore := a[0].Score - a[1].Score
						if diffScore > 0 {
							if a[0].RefLen-a[1].RefLen < diffScore/MatchScore {
								newArr = append(newArr, a[0].RefID)
							}
						} else { // diff score == 0
							if constructdbg.IsBubble(a[0].RefID, a[1].RefID, edgesArr) {
								puzzy = true
								newArr = append(newArr, a[0].RefID)
							}
						}
					}

					// extend path for  distinguish most likely path
					if len(newArr) == 0 {
						fmt.Printf("[findMostProbablePath]BACKWARD using Distinguish, ea: %v\n", ea)
						newArr = extendPathDistinguish(edgesArr, nd, nodesArr, ea, arr, si, seqLen, constructdbg.BACKWARD, opt)
						if len(newArr) == 0 {
							break
						}
					}

				} else {
					break
				}
			}

			// adjust si, e, nd, seqLen
			for _, item := range newArr {
				seqLen += (len(edgesArr[item].Utg.Ks) - opt.Kmer)
			}
			fmt.Printf("[findMostProbablePath]BACKWARD newArr: %v, extArr: %v\n", newArr, extArr)
			si, seqLen = RelocateLRRecordArr(arr, si, constructdbg.BACKWARD, newArr, seqLen, edgesArr, opt)
			e = edgesArr[newArr[len(newArr)-1]]

			for _, id := range newArr {
				if edgesArr[id].StartNID == nd.ID {
					nd = nodesArr[edgesArr[id].EndNID]
				} else {
					nd = nodesArr[edgesArr[id].StartNID]
				}
			}

			if puzzy {
				newArr[0] = 0
			}
			extArr = append(extArr, newArr...)
			if seqLen > arr[si].Start*6/5 || si == 0 || nd.ID == 0 {
				break
			}
		}
	}
	constructdbg.ReverseDBG_MAX_INTArr(extArr)

	nID = edgesArr[arr[pos].RefID].EndNID
	if arr[pos].Strand {
		nID = edgesArr[arr[pos].RefID].StartNID
	}
	if nID > 0 && pos < len(arr)-1 {
		e := edgesArr[arr[pos].RefID]
		nd := nodesArr[nID]
		si := pos
		seqLen := GetFlankLen(arr[si], constructdbg.FORWARD)
		for {
			fmt.Printf("[findMostProbablePath]FORWARD found e.ID: %v, nd.ID: %v, seqLen: %v, si: %v\n", e.ID, nd.ID, seqLen, si)
			var newArr []constructdbg.DBG_MAX_INT
			var puzzy bool
			if len(e.PathMat) > 0 && len(e.PathMat[0].IDArr) > 0 {
				d := constructdbg.FORWARD
				if e.StartNID == nd.ID {
					d = constructdbg.BACKWARD
				}
				arrNGS := GetExtendPath(e.PathMat[0].IDArr, e.ID, d)
				if len(arrNGS) > 0 {
					var perr bool
					for _, eID := range arrNGS {
						if edgesArr[eID].ID == 0 || edgesArr[eID].GetDeleteFlag() > 0 {
							perr = true
						}
					}
					if !perr {
						newArr = append(newArr, arrNGS...)
					}
				}
			}
			fmt.Printf("[findMostProbablePath]FORWARD after PathMat e.PathMat: %v\n\tnewArr: %v\n", e.PathMat, newArr)
			if len(newArr) == 0 {
				ea := constructdbg.GetNearEdgeIDArr(nd, e.ID)
				if len(ea) == 1 {
					newArr = append(newArr, ea[0])
				} else if len(ea) > 1 {
					var a []LRRecord
					for _, eID := range ea {
						idx, ok := SearchLRRecordArr(arr, si, constructdbg.FORWARD, eID, (seqLen+len(edgesArr[eID].Utg.Ks)-(opt.Kmer-1))*6/5)
						if ok && arr[idx].RefID == eID {
							a = append(a, arr[idx])
						}
					}
					minLen := GetMinLenEID(edgesArr, ea)
					if len(a) == 1 {
						if minLen > 2*opt.Kmer {
							if a[0].RefStart < a[0].Mch*1/5 && a[0].Mch > a[0].RefLen*4/5 {
								newArr = append(newArr, a[0].RefID)
							}
						}
					} else if len(a) == len(ea) {
						sort.Sort(LRRecordScoreArr(a))
						diffScore := a[0].Score - a[1].Score
						if diffScore > 0 {
							if a[0].RefLen-a[1].RefLen < diffScore/MatchScore {
								newArr = append(newArr, a[0].RefID)
							}
						} else { // diff score == 0
							if constructdbg.IsBubble(a[0].RefID, a[1].RefID, edgesArr) {
								puzzy = true
								newArr = append(newArr, a[0].RefID)
							}
						}
					}

					// extend path for  distinguish most likely path
					if len(newArr) == 0 {
						fmt.Printf("[findMostProbablePath]FORWARD using Distinguish, ea: %v\n", ea)
						newArr = extendPathDistinguish(edgesArr, nd, nodesArr, ea, arr, si, seqLen, constructdbg.FORWARD, opt)
						if len(newArr) == 0 {
							break
						}
					}

				} else {
					break
				}
			}

			// adjust si, e, nd, seqLen
			for _, item := range newArr {
				seqLen += (len(edgesArr[item].Utg.Ks) - opt.Kmer)
			}
			fmt.Printf("[findMostProbablePath] FORWARD newArr: %v, extArr: %v\n", newArr, extArr)
			si, seqLen = RelocateLRRecordArr(arr, si, constructdbg.FORWARD, newArr, seqLen, edgesArr, opt)
			e = edgesArr[newArr[len(newArr)-1]]
			for _, id := range newArr {
				if edgesArr[id].StartNID == nd.ID {
					nd = nodesArr[edgesArr[id].EndNID]
				} else {
					nd = nodesArr[edgesArr[id].StartNID]
				}
			}

			if puzzy { // need puzzy match
				newArr[0] = 0
			}
			extArr = append(extArr, newArr...)
			if seqLen > (arr[si].RefLen-arr[si].RefStart)*6/5 || si == len(arr)-1 || nd.ID == 0 {
				break
			}
		}
	}

	return extArr
} */

func Reverse(s string) string {
	r := []rune(s)
	for i, j := 0, len(r)-1; i < j; i, j = i+1, j-1 {
		r[i], r[j] = r[j], r[i]
	}
	return string(r)
}

func SplitBySpace(s string, num int) []string {
	arr := make([]string, num)
	count := 0
	for i := 0; i < len(s); {
		if s[i] == ' ' {
			i++
			continue
		}
		if count == num-1 {
			arr[count] = s[i:]
			count++
			break
		}
		j := i + 1
		for ; j < len(s); j++ {
			if s[j] == ' ' {
				break
			}
		}
		arr[count] = s[i:j]
		i = j
		count++
	}
	if count != num {
		log.Fatalf("[SplitBySpace] found count:%v != num: %v\n", count, num)
	}

	return arr
}

func ConvertLRRecord(pi PAFInfo) (lrd LRRecord) {
	//var err error
	R := pi.Sa
	tmp, err := strconv.Atoi(R[0])
	if err != nil {
		log.Fatalf("[ConvertLRRecord] ID: %v convert to int error: %v\n", R[0], err)
	}
	lrd.ID = constructdbg.DBG_MAX_INT(tmp)
	tmp, err = strconv.Atoi(R[5])
	if err != nil {
		log.Fatalf("[ConvertLRRecord] RefID: %v convert to int error: %v\n", R[5], err)
	}
	lrd.RefID = constructdbg.DBG_MAX_INT(tmp)

	lrd.RefLen, err = strconv.Atoi(R[6])
	if err != nil {
		log.Fatalf("[ConvertLRRecord] RefLen: %v convert to int error: %v\n", R[6], err)
	}

	lrd.RefStart, err = strconv.Atoi(R[7])
	if err != nil {
		log.Fatalf("[ConvertLRRecord] RefStart: %v convert to int error: %v\n", R[7], err)
	}

	lrd.RefEnd, err = strconv.Atoi(R[8])
	if err != nil {
		log.Fatalf("[ConvertLRRecord] RefEnd: %v convert to int error: %v\n", R[8], err)
	}

	if R[4] == "+" {
		lrd.Strand = constructdbg.PLUS
	} else {
		lrd.Strand = constructdbg.MINUS
	}

	lrd.Len, err = strconv.Atoi(R[1])
	if err != nil {
		log.Fatalf("[ConvertLRRecord] Len: %v convert to int error: %v\n", R[1], err)
	}

	lrd.Start, err = strconv.Atoi(R[2])
	if err != nil {
		log.Fatalf("[ConvertLRRecord] Start : %v convert to int error: %v\n", R[2], err)
	}

	lrd.End, err = strconv.Atoi(R[3])
	if err != nil {
		log.Fatalf("[ConvertLRRecord] End : %v convert to int error: %v\n", R[3], err)
	}

	lrd.MapNum, err = strconv.Atoi(R[9])
	if err != nil {
		log.Fatalf("[ConvertLRRecord] MapNum : %v convert to int error: %v\n", R[9], err)
	}

	lrd.GapMapNum, err = strconv.Atoi(R[10])
	if err != nil {
		log.Fatalf("[ConvertLRRecord] GapMapNum : %v convert to int error: %v\n", R[10], err)
	}
	lrd.ReadSeqBnt = pi.Seq
	return
}

func cleanLRRecordArr(arr []LRRecord, flankAllow int) []LRRecord {
	j := 0
	for _, r := range arr {
		if r.Start > flankAllow && r.End < r.Len-flankAllow {
			if r.RefEnd-r.RefStart > r.RefLen*3/4 && r.MapNum > r.GapMapNum*1/5 {
				arr[j] = r
				j++
			}
		} else if r.Start < flankAllow && r.End > r.Len-flankAllow {
			if r.End-r.Start > r.Len*3/4 {
				arr[j] = r
				j++
			}
		} else if r.Start <= flankAllow {
			if r.Strand == constructdbg.PLUS {
				if r.RefEnd > r.RefLen-flankAllow {
					arr[j] = r
					j++
				}
			} else {
				if r.RefStart < flankAllow {
					arr[j] = r
					j++
				}
			}
		} else { // r.End >= r.Len - flankAllow
			if r.Strand == constructdbg.PLUS {
				if r.RefStart < flankAllow {
					arr[j] = r
					j++
				}
			} else {
				if r.RefEnd > r.RefLen-flankAllow {
					arr[j] = r
					j++
				}
			}
		}
	}

	arr = arr[:j]
	return arr
}

func paraFindLongReadsMappingPath(rc <-chan []PAFInfo, wc chan LongReadMappingInfo, edgesArr []constructdbg.DBGEdge, nodesArr []constructdbg.DBGNode, opt Options) {
	//fragmentFreq := make([]int, 100)
	for {
		rArr, ok := <-rc
		if !ok {
			var guardPath LongReadMappingInfo
			wc <- guardPath
			break
		}

		// process PAFRecord array
		arr := make([]LRRecord, len(rArr))
		for i, R := range rArr {
			arr[i] = ConvertLRRecord(R)
		}

		// clean not whole length match
		/*flankdiff := 20
		j := 0
		for _, t := range arr {
			if t.RefConLen < t.RefLen*9/10 {
				if t.Start > flankdiff+t.RefConLen/10 &&
					(t.Len-(t.Start+t.ConLen)) > flankdiff+t.RefConLen/10 {
					continue
				}
			}
			arr[j] = t
			j++
		}
		arr = arr[:j] */

		//sort.Sort(LRRecordArr(arr))
		// debug code
		/*mapCount := make(map[constructdbg.DBG_MAX_INT]int)
		for _, R := range arr {
			if _, ok := mapCount[R.RefID]; ok {
				mapCount[R.RefID]++
			} else {
				mapCount[R.RefID] = 1
			}
		}
		for _, v := range mapCount {
			if v < 100 {
				fragmentFreq[v]++
			} else {
				fragmentFreq[99]++
			}
		}

		if len(arr) <= 1 {
			continue
		}

		fmt.Printf("[paraFindLongReadsMappingPath] read name: %v\n", strings.Split(rArr[0].Arr[2], " ")[1])
		for i, rd := range arr {
			arr[i].Mch, arr[i].Mis, arr[i].Ins, arr[i].Del = GetNM(rd)
			fmt.Printf("[paraFindLongReadsMappingPath] rd: %v\n", rd)
		}
		fmt.Printf("[paraFindLongReadsMappingPath] read name: %v\n", strings.Split(rArr[0].Arr[2], " ")[1])
		*/

		flankAllow := opt.Kmer * 2
		arr = cleanLRRecordArr(arr, flankAllow)
		/*if len(arr) > 0 {
			fmt.Printf("[paraFindLongReadsMappingPath]after clean arr[0].ID: %v\n", arr[0].ID)
		}*/
		// mapping ONT reads To the edge Path
		path := MappingONTMostProbablePath(arr, edgesArr, nodesArr, opt, flankAllow)
		// found most possible path
		//path := findMostProbablePath(arr, edgesArr, nodesArr, opt)
		//fmt.Printf("[paraFindLongReadsMappingPath]path: %v\n", path)
		if len(path[0]) >= 2 {
			var rmi LongReadMappingInfo
			rmi.ID = uint32(arr[0].ID)
			rmi.Path = path
			wc <- rmi
		}
	}
	/*for i, v := range fragmentFreq {
		fmt.Fprintf(os.Stderr, "%v\t%v\n", i, v)
	}*/
}

/*type MAFRecord struct {
	Arr [3]string
}*/

//type MAFRecordArr []MAFRecord

type LRInfo struct {
	ID  string
	Seq []byte
}

type PAFInfo struct {
	Sa  []string
	Seq []byte
}

func GetLRRead(ONTfp *bufio.Reader, format string) (LRInfo, error) {
	var li LRInfo
	s, err := constructcf.GetReadFileRecord(ONTfp, format, false)
	if err != nil {
		return li, err
	}
	//l := s.(*linear.Seq)
	li.ID = strconv.Itoa(int(s.ID))
	li.Seq = s.Seq
	//fmt.Printf("[GetLRRead]read ID: %v\n\tSeq: %v\n", l.ID, l.Seq)
	//fmt.Printf("[GetLRRead]len(li.Seq): %v\nli.Seq: %v\n", len(li.Seq), li.Seq)
	return li, err
}

func GetReadInfo(cs chan<- constructcf.ReadInfo, bufONTfp *bufio.Reader, format string) {
	for {
		s, err := constructcf.GetReadFileRecord(bufONTfp, format, false)
		if err == nil {
			cs <- s
		} else {
			if err == io.EOF {
				break
			} else {
				log.Fatalf("[GetReadInfo] read ONT file failed, err: %v\n", err)
			}
		}
	}

	close(cs)
}

func GetPAFRecord(paffn, ONTfn string, rc chan<- []PAFInfo, numCPU int) {
	ONTfp, err1 := os.Open(ONTfn)
	if err1 != nil {
		log.Fatalf("[GetPAFRecord] open ONT file: %s failed, err: %v\n", ONTfn, err1)
	}
	brONTfp := cbrotli.NewReaderSize(ONTfp, 1<<25)
	defer brONTfp.Close()
	bufONTfp := bufio.NewReader(brONTfp)
	paffp, err2 := os.Open(paffn)
	if err2 != nil {
		log.Fatalf("[GetPAFRecord] open PAF file: %s failed, err: %v\n", paffn, err2)
	}
	defer ONTfp.Close()
	defer paffp.Close()
	pafbuffp := bufio.NewReader(paffp)
	format := constructcf.GetReadsFileFormat(ONTfn)
	//ONTfafp := bufio.NewReader(ONTfp)
	//ONTfafp := fasta.NewReader(bufONTfp, linear.NewSeq("", nil, alphabet.DNA))
	cs := make(chan constructcf.ReadInfo, 10000)
	go GetReadInfo(cs, bufONTfp, format)

	var pa []PAFInfo
	var ri constructcf.ReadInfo
	for {
		line, err := pafbuffp.ReadString('\n')
		if err != nil {
			if err == io.EOF {
				break
			} else {
				log.Fatalf("[GetPAFRecord] Read line: %v, err: %v\n", line, err)
			}
		}

		//fmt.Printf("[GetMAFRecord] line: %v\n", line)
		sa := strings.Split(line[:len(line)-1], "\t")
		rID, err := strconv.Atoi(sa[0])
		if err != nil {
			log.Fatalf("[GetPAFReads] read file: %s, read ID not integer error: %v\n", paffn, err)
		}
		if rID != int(ri.ID) {
			for rID != int(ri.ID) {
				ri = <-cs
				if ri.ID == 0 {
					log.Fatalf("[GetPAFReads] read file: %s, error: %v\n", ONTfn, err)
				}
			}
		}
		var p PAFInfo
		p.Sa = sa
		p.Seq = ri.Seq
		//fmt.Printf("[GetPAFReads] p.Seq ID: %v, li.ID: %v\n", sa[5], li.ID)
		if len(pa) > 0 && sa[0] != pa[0].Sa[0] {

			rc <- pa
			//fmt.Printf("[GetPAFReads] pa[0]: %v\n", pa[0])
			var na []PAFInfo
			pa = na
		}
		pa = append(pa, p)
	}
	if len(pa) > 0 {
		rc <- pa
	}

	// notice para goroutinues the channel has not any more data
	close(rc)
}
