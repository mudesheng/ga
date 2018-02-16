package deconstructdbg

import (
	"bufio"
	"container/list"
	"fmt"
	"io"
	"log"
	"os"
	"sort"
	"strconv"
	"strings"

	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/sam"
	"github.com/mudesheng/GA/constructdbg"
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
	RefID              constructdbg.DBG_MAX_INT
	//Pair                        uint8 // 1 note Ref Name RefID/1, 2 note RefID/2
	RefStart, RefConLen, RefLen int
	Start, ConLen, Len          int
	RefSeq, Seq                 string
	Score                       int
	Strand                      bool // true note reverse strand
}

func (R LRRecord) GetNM() int {
	return R.Mis + R.Ins + R.Del
}

func (R LRRecord) GetRefEnd() int {
	return R.RefStart + R.RefConLen
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

type LRRecordScoreArr []LRRecord

func (a LRRecordScoreArr) Len() int {
	return len(a)
}

func (a LRRecordScoreArr) Less(i, j int) bool {
	return a[i].Score > a[j].Score
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

func findSeedEID(arr []LRRecord, edgesArr []constructdbg.DBGEdge) int {
	pos := -1
	var maxScore int
	for i, item := range arr {
		e := edgesArr[item.RefID]
		if e.GetUniqueFlag() > 0 && maxScore < item.Score {
			pos = i
			maxScore = item.Score
		}
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

func SearchLRRecordArr(arr []LRRecord, aidx int, direction uint8, eID constructdbg.DBG_MAX_INT, searchLen int) (idx int, ok bool) {
	if direction == constructdbg.FORWARD {
		max := arr[aidx].Start + arr[aidx].ConLen + searchLen
		for i := aidx + 1; i < len(arr); i++ {
			if arr[i].Start > max {
				break
			}
			if arr[i].RefID == eID {
				idx = i
				ok = true
				return idx, ok
			}
		}
	} else {
		min := arr[idx].Start - searchLen
		for i := aidx - 1; i >= 0; i-- {
			if arr[i].Start+arr[i].ConLen < min {
				break
			}
			if arr[i].RefID == eID {
				idx = i
				ok = true
				return idx, ok
			}
		}
	}

	return idx, ok
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
		if R.Strand {
			return R.RefStart
		} else {
			return R.RefLen - (R.RefStart + R.RefConLen)
		}
	} else { // BACKWARD
		if R.Strand {
			return R.RefLen - (R.RefStart + R.RefConLen)
		} else {
			return R.RefStart
		}
	}
}

func RelocateLRRecordArr(arr []LRRecord, si int, direction uint8, newArr []constructdbg.DBG_MAX_INT, searchLen int, edgesArr []constructdbg.DBGEdge, opt Options) (int, int) {
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
}

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

func GetNM(R LRRecord) (mch, mis, ins, del int) {
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
}

func Getscore(R LRRecord, begin, end int) (score int) {
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
}
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

// region special by Refence(Long Reads)
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
}

func findMostProbablePath(arr []LRRecord, edgesArr []constructdbg.DBGEdge, nodesArr []constructdbg.DBGNode, opt Options) (extArr []constructdbg.DBG_MAX_INT) {
	// found seed edge ID
	MatchScore := 5
	pos := findSeedEID(arr, edgesArr)
	if pos < 0 || arr[pos].RefConLen < opt.Kmer {
		return extArr
	}
	fmt.Printf("[findMostProbablePath] found seed Record: %v\n", arr[pos])
	extArr = append(extArr, arr[pos].RefID)
	nID := edgesArr[arr[pos].RefID].StartNID
	if arr[pos].Strand {
		nID = edgesArr[arr[pos].RefID].EndNID
	}
	if nID > 0 && pos > 0 {
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
			/*for _, item := range newArr {
				seqLen += (len(edgesArr[item].Utg.Ks) - opt.Kmer)
			}*/
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
			/*for _, item := range newArr {
				seqLen += (len(edgesArr[item].Utg.Ks) - opt.Kmer)
			}*/
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
}

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

func paraFindLongReadsMappingPath(rc chan []MAFRecord, wc chan []constructdbg.DBG_MAX_INT, edgesArr []constructdbg.DBGEdge, nodesArr []constructdbg.DBGNode, opt Options) {
	fragmentFreq := make([]int, 100)
	for {
		rArr := <-rc
		if len(rArr) == 0 {
			var guardPath []constructdbg.DBG_MAX_INT
			wc <- guardPath
			break
		}

		// process MAFRecord array
		var arr []LRRecord
		for _, R := range rArr {
			var lrd LRRecord
			sc := strings.Split(R.Arr[0], " ")[1]
			fmt.Printf("[paraFindLongReadsMappingPath]sc: %v\n", sc)
			var err error
			lrd.Score, err = strconv.Atoi(strings.Split(sc, "=")[1])
			if err != nil {
				log.Fatalf("[paraFindLongReadsMappingPath] score: %v convert to int error: %v\n", sc, err)
			}
			//sa := strings.Split(R.Arr[1], " ")
			sa := SplitBySpace(R.Arr[1], 7)
			fmt.Printf("[paraFindLongReadsMappingPath]sa: %v\n", sa)
			//ta := strings.Split(sa[1], "/")
			id, err := strconv.Atoi(sa[1])
			if err != nil {
				log.Fatalf("[paraFindLongReadsMappingPath] Ref Name: %v convert to int error: %v\n", sa, err)
			}
			lrd.RefID = constructdbg.DBG_MAX_INT(id)
			id, err = strconv.Atoi(sa[2])
			if err != nil {
				log.Fatalf("[paraFindLongReadsMappingPath] Ref Start position: %v convert to int error: %v\n", sa[2], err)
			}
			lrd.RefStart = id
			id, err = strconv.Atoi(sa[3])
			if err != nil {
				log.Fatalf("[paraFindLongReadsMappingPath] Consume length : %v convert to int error: %v\n", sa[3], err)
			}
			lrd.RefConLen = id
			id, err = strconv.Atoi(sa[5])
			if err != nil {
				log.Fatalf("[paraFindLongReadsMappingPath] Ref length : %v convert to int error: %v\n", sa[5], err)
			}
			lrd.RefLen = id
			lrd.RefSeq = sa[6]

			// Query parsing
			//sa = strings.Split(R.Arr[2], " ")
			sa = SplitBySpace(R.Arr[2], 7)
			fmt.Printf("[paraFindLongReadsMappingPath]sa: %v\n", sa)
			id, err = strconv.Atoi(sa[2])
			if err != nil {
				log.Fatalf("[paraFindLongReadsMappingPath] Query start : %v convert to int error: %v\n", sa[2], err)
			}
			lrd.Start = id
			id, err = strconv.Atoi(sa[3])
			if err != nil {
				log.Fatalf("[paraFindLongReadsMappingPath] Query consume length : %v convert to int error: %v\n", sa[3], err)
			}
			lrd.ConLen = id
			id, err = strconv.Atoi(sa[5])
			if err != nil {
				log.Fatalf("[paraFindLongReadsMappingPath] Query length : %v convert to int error: %v\n", sa[5], err)
			}
			lrd.Len = id
			lrd.Seq = sa[6]
			if sa[4] == "-" {
				lrd.Strand = true
				lrd.Start = lrd.Len - lrd.Start - lrd.ConLen
				lrd.RefSeq = Reverse(lrd.RefSeq)
				lrd.Seq = Reverse(lrd.Seq)
			}

			arr = append(arr, lrd)
		}

		// clean not whole length match
		flankdiff := 20
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
		arr = arr[:j]

		sort.Sort(LRRecordArr(arr))
		// debug code
		mapCount := make(map[constructdbg.DBG_MAX_INT]int)
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

		// found most possible path
		path := findMostProbablePath(arr, edgesArr, nodesArr, opt)
		if len(path) > 2 {
			wc <- path
		}
	}
	for i, v := range fragmentFreq {
		fmt.Fprintf(os.Stderr, "%v\t%v\n", i, v)
	}
}

type MAFRecord struct {
	Arr [3]string
}

//type MAFRecordArr []MAFRecord

func GetMAFRecord(maffn string, rc chan []MAFRecord, numCPU int) {
	fp, err := os.Open(maffn)
	if err != nil {
		log.Fatalf("[GetLastRecord] open file: %s failed, err: %v\n", maffn, err)
	}
	defer fp.Close()
	buffp := bufio.NewReader(fp)

	var mra []MAFRecord
	for {
		line, err := buffp.ReadString('\n')
		if err != nil {
			if err == io.EOF {
				rc <- mra
				break
			} else {
				log.Fatalf("[GetMAFRecord] Read line: %v, err: %v\n", line, err)
			}
		}
		if line[0] == '#' || line[0] == '\n' {
			continue
		}
		line2, err := buffp.ReadString('\n')
		if err != nil {
			log.Fatalf("[GetMAFRecord] Read line: %v, err: %v\n", line2, err)
		}
		line3, err := buffp.ReadString('\n')
		if err != nil {
			log.Fatalf("[GetMAFRecord] Read line: %v, err: %v\n", line3, err)
		}
		//fmt.Printf("[GetMAFRecord] line: %v\n", line)
		var mafR MAFRecord
		mafR.Arr[0] = line[:len(line)-1]
		mafR.Arr[1] = line2[:len(line2)-1]
		mafR.Arr[2] = line3[:len(line3)-1]
		fmt.Printf("[GetMAFRecord] mafR[0]: %v\n\tmafR[1]: %v\n\tmafR[2]: %v\n", mafR.Arr[0], mafR.Arr[1], mafR.Arr[2])
		if len(mra) > 0 {
			lrID1 := strings.Split(mra[0].Arr[2], " ")[1]
			lrID2 := strings.Split(mafR.Arr[2], " ")[1]
			if lrID1 != lrID2 {
				rc <- mra
				var n []MAFRecord
				mra = n
			}
		}
		mra = append(mra, mafR)
	}

	for i := 0; i < numCPU; i++ {
		var nilArr []MAFRecord
		rc <- nilArr
	}
}
