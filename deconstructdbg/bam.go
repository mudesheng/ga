package deconstructdbg

import (
	"bufio"
	"container/list"
	"fmt"
	"io"
	"log"
	"math"
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
	Ins, Del, Mis, Mch          int // Insertion, Deletion and Mismatch base number
	RefID                       constructdbg.DBG_MAX_INT
	Pair                        uint8 // 1 note Ref Name RefID/1, 2 note RefID/2
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

func findSeedEID(arr []LRRecord, edgesArr []constructdbg.DBGEdge) (int, int) {
	pos := -1
	var maxScore int
	for i, item := range arr {
		if edgesArr[item.RefID].GetUniqueFlag() > 0 {
			if item.Pair > 0 {
				pIdx := PairIndex(arr, i+1, item)
				if pIdx >= 0 {
					refMappingLen := edgesArr[item.RefID].GetSeqLen() - 2*item.RefLen + item.RefConLen + arr[pIdx].RefConLen
					readMappingLen := arr[pIdx].Start + arr[pIdx].ConLen - item.Start
					if math.Abs(float64(refMappingLen-readMappingLen)) < float64(readMappingLen/10) {
						if maxScore < item.Score+arr[pIdx].Score {
							pos = i
							maxScore = item.Score + arr[pIdx].Score
						}
					}
				} else {
					if pos >= 0 {
						if arr[pos].Pair > 0 {
							if maxScore < item.Score {
								pos = i
								maxScore = item.Score
							}
						} else {
							pos = i
							maxScore = item.Score
						}
					} else {
						pos = i
						maxScore = item.Score
					}
				}
			} else {
				if pos < 0 || (pos >= 0 && arr[pos].Pair == 0 && maxScore < item.Score) {
					pos = i
					maxScore = item.Score
				}
			}
		}
	}

	return pos, maxScore
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
		RefMax := arr[aidx].RefPos + searchLen
		for i := aidx + 1; i < len(arr); i++ {
			if arr[i].RefPos > RefMax {
				break
			}
			if arr[i].ID == eID {
				idx = i
				ok = true
				return idx, ok
			}
		}
	} else {
		RefMin := arr[idx].RefPos - searchLen
		for i := aidx - 1; i >= 0; i-- {
			if arr[i].RefPos < RefMin {
				break
			}
			if arr[i].ID == eID {
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
	minLen := len(edgesArr[eIDArr[0]].Utg.Ks)
	for _, eID := range eIDArr[1:] {
		if len(edgesArr[eID].Utg.Ks) < minLen {
			minLen = len(edgesArr[eID].Utg.Ks)
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

func GetFlankLen(R LRRecord, e constructdbg.DBGEdge, diretion uint8) int {
	if diretion == constructdbg.FORWARD {
		if R.Strand {
			return R.Pos
		} else {
			return e.GetSeqLen() - (R.Pos + R.Mch + R.Mis + R.Ins)
		}
	} else { // BACKWARD
		if R.Strand {
			return e.GetSeqLen() - (R.Pos + R.Mch + R.Mis + R.Ins)
		} else {
			return R.Pos
		}
	}
}

func RelocateLRRecordArr(arr []LRRecord, si int, direction uint8, newArr []constructdbg.DBG_MAX_INT, searchLen int, edgesArr []constructdbg.DBGEdge, opt Options) (int, int) {
	nl := searchLen
	for _, eID := range newArr {
		nl += edgesArr[eID].GetSeqLen() - (opt.Kmer - 1)
		idx, ok := SearchLRRecordArr(arr, si, direction, eID, nl*6/5)
		if ok && arr[idx].ID == eID {
			si = idx
			nl = GetFlankLen(arr[idx], edgesArr[eID], direction)
			if direction == constructdbg.FORWARD {
				if arr[idx].Strand {
					nl = arr[idx].Pos + arr[idx].GetQryConsumeLen()
				} else {
					nl = edgesArr[eID].GetSeqLen() - arr[idx].Pos
				}
			}
		}
		//searchLen = nl - (len(edgesArr[eID].Utg.Ks) - opt.Kmer)
		//nl -= (len(edgesArr[eID].Utg.Ks) - opt.Kmer)
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

func Getscore(R LRRecord, begin, end int) (score int) {
	if begin == R.RefPos && end == R.GetRefEnd() {
		score = R.GetScore()
		return score
	}
	indel := 0
	st, ed := begin-R.RefPos, end-R.RefPos
	pos := 0
	for _, op := range R.CIGAR {
		switch op.Type() {
		case sam.CigarMatch:
			pos += op.Len()
		case sam.CigarDeletion:
			if st <= pos && pos < ed {
				if pos+op.Len() <= ed {
					indel += op.Len()
				} else {
					indel += ed - pos
				}
			}
			pos += op.Len()
		case sam.CigarInsertion:
			if st <= pos && pos < ed {
				indel += op.Len()
			}
		}
		if pos >= ed {
			break
		}
	}
	pos, mis, mch := 0, 0, 0
	for i := 0; i < len(R.MD); {
		if R.MD[i] == '^' {
			for i++; i < len(R.MD); i++ {
				if 'A' <= R.MD[i] && R.MD[i] <= 'Z' {
					pos++
				} else {
					break
				}
			}
		} else if '0' <= R.MD[i] && R.MD[i] <= '9' {
			j := i + 1
			for ; j < len(R.MD); j++ {
				if R.MD[j] < '0' || R.MD[j] > '9' {
					break
				}
			}
			num, err := strconv.Atoi(R.MD[i:j])
			if err != nil {
				log.Fatalf("[GetScore] convert: %v to integar err: %v\n", R.MD[i:j], err)
			}
			if st <= pos && pos+num <= ed {
				mch += num
			} else if pos < st && pos+num > st && pos+num <= ed {
				mch += pos + num - st
			} else if pos >= st && pos < ed && pos+num > ed {
				mch += ed - pos
			} else if pos <= st && pos+num >= ed {
				mch += ed - st
			}
			//fmt.Printf("[GetScore] MD: %v, convert to num: %v\n\tst: %v, ed: %v, pos: %v, mch: %v\n", R.MD[i:j], num, st, ed, pos, mch)
			pos += num
			i = j
		} else if 'A' <= R.MD[i] && R.MD[i] <= 'Z' {
			if st <= pos && pos < ed {
				mis++
				pos++
			}
			i++
		} else {
			log.Fatalf("[Getscore] unknow char : %c\n", R.MD[i])
		}

		if pos >= ed {
			break
		}
	}
	score = mch - mis - indel
	fmt.Printf("[Getscore] RefPos: %v, RefEnd : %v\n\tbegin: %v, end: %v, mch: %v, indel: %v, mis: %v, score: %v\n\tMD: %v\n", R.RefPos, R.GetRefEnd(), begin, end, mch, indel, mis, score, R.MD)
	return score
}
func priorIndex(arr []LRRecord, idx int) int {
	pi := -1
	for i := idx - 1; i >= 0; i-- {
		if arr[i].ID != 0 {
			pi = i
			break
		}
	}

	return pi
}

func nextIndex(arr []LRRecord, idx int) int {
	ni := -1
	for i := idx + 1; i < len(arr); i++ {
		if arr[i].ID != 0 {
			ni = i
			break
		}
	}

	return ni
}

// region special by Refence(Long Reads)
func GetScoreFixedLen(edgesArr []constructdbg.DBGEdge, epi ExtPathInfo, flankRefPos int, direction uint8) int {
	var score int
	if direction == constructdbg.FORWARD {
		for i, _ := range epi.Path {
			if epi.LRRecordArr[i].ID != 0 {
				begin, end := epi.LRRecordArr[i].RefPos, epi.LRRecordArr[i].GetRefEnd()
				pi := priorIndex(epi.LRRecordArr, i)
				ni := nextIndex(epi.LRRecordArr, i)
				if pi >= 0 && epi.LRRecordArr[pi].GetRefEnd() > begin {
					begin += (epi.LRRecordArr[pi].GetRefEnd() - begin) / 2
				}
				if ni > i && epi.LRRecordArr[ni].RefPos < end {
					end -= (end - epi.LRRecordArr[ni].RefPos) / 2
				}
				if begin > flankRefPos {
					break
				}
				if end > flankRefPos {
					end = flankRefPos
				}
				score += Getscore(epi.LRRecordArr[i], begin, end)
			}
		}
	} else { // BACKWARD
		for i, _ := range epi.Path {
			if epi.LRRecordArr[i].ID != 0 {
				begin, end := epi.LRRecordArr[i].RefPos, epi.LRRecordArr[i].GetRefEnd()
				pi := priorIndex(epi.LRRecordArr, i)
				ni := nextIndex(epi.LRRecordArr, i)
				if pi >= 0 && epi.LRRecordArr[pi].RefPos < end {
					end -= (end - epi.LRRecordArr[pi].RefPos) / 2
				}
				if ni > i && epi.LRRecordArr[ni].GetRefEnd() > begin {
					begin += (epi.LRRecordArr[ni].GetRefEnd() - begin) / 2
				}
				if end < flankRefPos {
					break
				}
				if begin < flankRefPos {
					begin = flankRefPos
				}
				score += Getscore(epi.LRRecordArr[i], begin, end)
			}
		}
	}

	return score
}

func IsInRefRegionContain(a1, a2 LRRecord) (constructdbg.DBG_MAX_INT, bool) {
	a1E := a1.RefPos + a1.Mch + a1.Mis + a1.Del
	a2E := a2.RefPos + a2.Mch + a2.Mis + a2.Del
	if (a1.RefPos < a2.RefPos && a2E < a1E) || (a2.RefPos < a1.RefPos && a1E < a2E) || (a1.RefPos == a2.RefPos && a1E == a2E) {
		s1, s2 := a1.GetScore(), a2.GetScore()
		if s1 > s2 {
			return a1.ID, true
		} else {
			return a2.ID, true
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
			if ok && lrArr[idx].ID == eID {
				fmt.Printf("[extendPathDistinguish] lrArr[%v]: %v\n", idx, lrArr[idx])
				lrArrIdx = idx
				sl = GetFlankLen(lrArr[idx], edgesArr[eID], direction)
				if j > 0 && epi.LRRecordArr[j-1].ID != 0 { // if has a contain edge mapped, ignore...
					if id, ok := IsInRefRegionContain(epi.LRRecordArr[j-1], lrArr[idx]); ok {
						if id == lrArr[idx].ID {
							epi.LRRecordArr[j-1].ID = 0
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
		var flankRefPos int
		if direction == constructdbg.FORWARD {
			flankRefPos = lrArr[si].GetRefEnd() + opt.ExtLen
		} else {
			flankRefPos = lrArr[si].RefPos - opt.ExtLen
		}
		epi.Score = GetScoreFixedLen(edgesArr, epi, flankRefPos, direction) // region special by Refence(Long Reads)
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
	pos, maxScore := findSeedEID(arr, edgesArr)
	if pos < 0 || arr[pos].RefConLen < opt.Kmer {
		return extArr
	}

	fmt.Printf("[findMostProbablePath] found seed Record: %v\n", arr[pos])
	extArr = append(extArr, arr[pos].ID)
	nID := edgesArr[arr[pos].ID].StartNID
	if arr[pos].Strand {
		nID = edgesArr[arr[pos].ID].EndNID
	}
	if nID > 0 && pos > 0 {
		e := edgesArr[arr[pos].ID]
		nd := nodesArr[nID]
		si := pos
		seqLen := GetFlankLen(arr[si], e, constructdbg.BACKWARD)
		/*if arr[si].Strand {
			seqLen = len(e.Utg.Ks) - (arr[si].Pos + arr[si].Mch + arr[si].Mis + arr[si].Del)
		} else {
			seqLen = arr[si].Pos
		}*/
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
						if ok && arr[idx].ID == eID {
							a = append(a, arr[idx])
						}
					}

					minLen := GetMinLenEID(edgesArr, ea)
					if len(a) == 1 {
						if minLen > 2*opt.Kmer {
							if a[0].RefPos < a[0].Mch*1/5 { // if the long reads start position
								newArr = append(newArr, a[0].ID)
							} else if a[0].Mch > edgesArr[a[0].ID].GetSeqLen()*4/5 {
								newArr = append(newArr, a[0].ID)
							}
						}
					} else if len(a) == len(ea) {
						sort.Sort(LRRecordArr(a))
						diffScore := a[0].GetScore() - a[1].GetScore()
						if diffScore > 0 {
							if a[1].Pos+a[1].GetQryConsumeLen() < edgesArr[a[1].ID].GetSeqLen() {
								newArr = append(newArr, a[0].ID)
							} else if edgesArr[a[0].ID].GetSeqLen()-edgesArr[a[1].ID].GetSeqLen() < diffScore {
								newArr = append(newArr, a[0].ID)
							}
						} else { // diff score == 0
							if constructdbg.IsBubble(a[0].ID, a[1].ID, edgesArr) {
								puzzy = true
								newArr = append(newArr, a[0].ID)
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
			if len(newArr) > 1 {
				nID = FindShareNID(newArr[len(newArr)-2], newArr[len(newArr)-1], edgesArr)
				nd = nodesArr[nID]
			}
			if e.StartNID == nd.ID {
				nd = nodesArr[e.EndNID]
			} else {
				nd = nodesArr[e.StartNID]
			}

			if puzzy {
				newArr[0] = 0
			}
			extArr = append(extArr, newArr...)
			if seqLen > arr[si].RefPos*6/5 || si == 0 {
				break
			}
		}
	}
	constructdbg.ReverseDBG_MAX_INTArr(extArr)

	nID = edgesArr[arr[pos].ID].EndNID
	if arr[pos].Strand {
		nID = edgesArr[arr[pos].ID].StartNID
	}
	if nID > 0 && pos < len(arr)-1 {
		e := edgesArr[arr[pos].ID]
		nd := nodesArr[nID]
		si := pos
		seqLen := GetFlankLen(arr[si], e, constructdbg.FORWARD)
		/*if !arr[si].Strand {
			seqLen = len(e.Utg.Ks) - arr[si].Pos
		} else {
			seqLen = arr[si].Pos + arr[si].Mch + arr[si].Mis + arr[si].Del
		}*/
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
						if ok && arr[idx].ID == eID {
							a = append(a, arr[idx])
						}
					}
					minLen := GetMinLenEID(edgesArr, ea)
					if len(a) == 1 {
						if minLen > 2*opt.Kmer {
							if a[0].RefLen-a[0].RefPos < a[0].Mch*6/5 { // Long reads flank
								newArr = append(newArr, a[0].ID)
							} else {
								if a[0].Mch > edgesArr[a[0].ID].GetSeqLen()*4/5 {
									newArr = append(newArr, a[0].ID)
								}
							}
						}
					} else if len(a) == len(ea) {
						sort.Sort(LRRecordArr(a))
						diffScore := a[0].GetScore() - a[1].GetScore()
						if diffScore > 0 {
							if a[1].Pos+a[1].GetQryConsumeLen() < edgesArr[a[1].ID].GetSeqLen() {
								newArr = append(newArr, a[0].ID)
							} else if edgesArr[a[0].ID].GetSeqLen()-edgesArr[a[1].ID].GetSeqLen() < diffScore {
								newArr = append(newArr, a[0].ID)
							}
						} else { // diff score == 0
							if constructdbg.IsBubble(a[0].ID, a[1].ID, edgesArr) {
								puzzy = true
								newArr = append(newArr, a[0].ID)
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
			if len(newArr) > 1 {
				nID = FindShareNID(newArr[len(newArr)-2], newArr[len(newArr)-1], edgesArr)
				nd = nodesArr[nID]
			}
			if e.StartNID == nd.ID {
				nd = nodesArr[e.EndNID]
			} else {
				nd = nodesArr[e.StartNID]
			}
			if puzzy { // need puzzy match
				newArr[0] = 0
			}
			extArr = append(extArr, newArr...)
			if seqLen > (arr[si].RefLen-arr[si].RefPos)*6/5 || si == len(arr)-1 {
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
			var err error
			lrd.Score, err = strconv.Atoi(strings.Split(sc, "=")[1])
			if err != nil {
				log.Fatalf("[paraFindLongReadsMappingPath] score: %v convert to int error: %v\n", sc, err)
			}
			sa := strings.Split(R.Arr[1], " ")
			ta := strings.Split(sa[1], "/")
			id, err := strconv.Atoi(ta[0])
			if err != nil {
				log.Fatalf("[paraFindLongReadsMappingPath] Ref Name: %v convert to int error: %v\n", ta, err)
			}
			lrd.RefID = constructdbg.DBG_MAX_INT(id)
			if len(ta) > 1 {
				p, err := strconv.Atoi(ta[1])
				if err != nil {
					log.Fatalf("[paraFindLongReadsMappingPath] pair info: %v convert to int error: %v\n", ta, err)
				}
				if p != 1 && p != 2 {
					log.Fatalf("[paraFindLongReadsMappingPath] pair: %v must == 1 or ==2\n", p)
				}
				lrd.Pair = uint8(p)
			}
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
			sa = strings.Split(R.Arr[2], " ")
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
		for i, t := range arr {
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
			if _, ok := mapCount[R.ID]; ok {
				mapCount[R.ID]++
			} else {
				mapCount[R.ID] = 1
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

		// debug
		fmt.Printf("[paraFindLongReadsMappingPath] read name: %v\n", strings.Split(rArr[0].Arr[2])[1])
		for _, rd := range arr {
			fmt.Printf("[paraFindLongReadsMappingPath] rd: %v\n", rd)
		}
		fmt.Printf("[paraFindLongReadsMappingPath] read name: %v\n", strings.Split(rArr[0].Arr[2])[1])

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
		if line[0] == '#' {
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
		var mafR MAFRecord
		mafR.arr[0] = line[:len(line)-1]
		mafR.arr[1] = line2[:len(line2)-1]
		mafR.arr[2] = line3[:len(line3)-1]
		if len(mra) > 0 {
			lrID1 := strings.Split(mra[0].arr[2], " ")[1]
			lrID2 := strings.Split(mafR.arr[2], " ")[1]
			if lrID1 != lrID2 {
				rc <- mra
				var n []MAFRecord
				mra = n
			}
		}
		mra = append(mra, mafR)
	}

	for i := 0; i < numCPU; i++ {
		var nilArr MAFRecordArr
		rc <- nilArr
	}
}
