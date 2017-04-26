package deconstructdbg

import (
	"container/list"
	"fmt"
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
	RefPos, Pos        int
	Idty               float32
	RefLen             int
	MapQ               uint8
	Strand             bool
	ID                 constructdbg.DBG_MAX_INT
	CIGAR              sam.Cigar
	MD                 string
}

func (R LRRecord) GetScore() int {
	return R.Mch - R.Mis - R.Ins - R.Del
}

func (R LRRecord) GetRefEnd() int {
	return R.RefPos + R.Mch + R.Mis + R.Del
}

func (R LRRecord) GetQryConsumeLen() int {
	return R.Mch + R.Mis + R.Ins
}

func (R LRRecord) GetRefConsumeLen() int {
	return R.Mch + R.Mis + R.Del
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
	return a[i].GetScore() > a[j].GetScore()
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

func findSeedEID(arr []LRRecord, edgesArr []constructdbg.DBGEdge) int {
	pos := -1
	var max int
	for i, item := range arr {
		if edgesArr[item.ID].GetUniqueFlag() > 0 && item.Mch > max {
			max = item.Mch
			pos = i
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
	Path   []constructdbg.DBG_MAX_INT
	ExtLen int
	Score  int
	Nd     constructdbg.DBGNode
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
		return R.GetScore()
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
			} else if pos < st && pos+num <= ed {
				mch += pos + num - st
			} else if st >= pos && pos+num > ed {
				mch += ed - pos
			} else if pos <= st && pos+num >= ed {
				mch += ed - st
			}
			pos += num
			i = j
		} else if 'A' <= R.MD[i] && R.MD[i] <= 'Z' {
			if st <= pos && pos < ed {
				mis++
			}
			i++
		} else {
			log.Fatalf("[Getscore] unknow char : %c\n", R.MD[i])
		}

		if pos >= ed {
			break
		}
	}
	fmt.Printf("[Getscore] begin: %v, end: %v, indel: %v, mis: %v\n\tMD: %v\n", begin, end, indel, mis, R.MD)
	score = mch - mis - indel
	return score
}

// region special by Refence(Long Reads)
func GetScoreFixedLen(mapArr []LRRecord, edgesArr []constructdbg.DBGEdge, epi ExtPathInfo, opt Options, direction uint8) int {
	var score int
	var extLen int
	if direction == constructdbg.FORWARD {
		for i, R := range mapArr {
			begin, end := R.RefPos, R.GetRefEnd()
			if i > 0 && begin < mapArr[i-1].GetRefEnd() {
				begin += (mapArr[i-1].GetRefEnd() - R.RefPos) / 2
			}
			if i < len(mapArr)-1 && end > mapArr[i+1].RefPos {
				end -= (end - mapArr[i+1].RefPos) / 2
			}
			if extLen+(end-begin) > opt.ExtLen {
				end = opt.ExtLen - extLen + begin
			}
			//fmt.Printf("[GetScoreFixedLen] begin: %v, end: %v, indel: %v\n\tMD: %v\n", begin, end, indel, R.MD)
			score += Getscore(R, begin, end)
			extLen += end - begin
			if extLen >= opt.ExtLen {
				break
			}
		}
	} else { // BACKWARD
		for i, R := range mapArr {
			begin, end := R.RefPos, R.GetRefEnd()
			if i > 0 && end > mapArr[i-1].RefPos {
				end -= (end - mapArr[i-1].RefPos) / 2
			}
			if i < len(mapArr)-1 && begin < mapArr[i+1].GetRefEnd() {
				begin += (mapArr[i+1].GetRefEnd() - begin) / 2
			}
			if extLen+(end-begin) > opt.ExtLen {
				begin = end - (opt.ExtLen - extLen)
			}
			score += Getscore(R, begin, end)
			extLen += end - begin
			if extLen >= opt.ExtLen {
				break
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
		var mapArr []LRRecord
		fmt.Printf("[extendPathDistinguish] epi: %v\n", epi)
		for _, eID := range epi.Path {
			sl += edgesArr[eID].GetSeqLen() - (opt.Kmer - 1)
			idx, ok := SearchLRRecordArr(lrArr, lrArrIdx, direction, eID, sl*6/5)
			if ok && lrArr[idx].ID == eID {
				fmt.Printf("[extendPathDistinguish] lrArr[%v]: %v\n", idx, lrArr[idx])
				lrArrIdx = idx
				sl = GetFlankLen(lrArr[idx], edgesArr[eID], direction)
				if len(mapArr) > 0 { // if has a contain edge mapped, ignore...
					if id, ok := IsInRefRegionContain(mapArr[len(mapArr)-1], lrArr[idx]); ok {
						if id == lrArr[idx].ID {
							mapArr[len(mapArr)-1] = lrArr[idx]
						}
					} else {
						mapArr = append(mapArr, lrArr[idx])
					}

				} else {
					mapArr = append(mapArr, lrArr[idx])
				}
			}
		}
		epi.Score = GetScoreFixedLen(mapArr, edgesArr, epi, opt, direction) // region special by Refence(Long Reads)
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
	/*var refIDArr []constructdbg.DBG_MAX_INT
	var MchArr []int
	j := -1 // MchArr index
	for _, R := range arr {
		if j >= 0 && R.RefID == refIDArr[j] {
			MchArr[j] += R.Mch
		} else {
			refIDArr = append(refIDArr, R.RefID)
			MchArr = append(MchArr, R.Mch)
			j++
		}
	}*/
	pos := findSeedEID(arr, edgesArr)
	if pos < 0 || arr[pos].GetScore() < opt.Kmer {
		return extArr
		//log.Fatalf("[findMostProbablePath] findSeedEID pos: %v <= 0\n", pos)
	}
	//SeedEID := refIDArr[pos]
	// find SeedEID mapping strand, start and end index of arr
	/*strand, si, ei := findMapInfo(SeedEID, pos, arr)
	if strand {
		constructdbg.ReverseDBG_MAX_INTArr(refIDArr)
		ReverseLRRcordArr(arr)
		pos = len(refIDArr) - 1 - pos
	}*/

	// search possible path
	//seededge := edgesArr[SeedEID]
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

func paraFindLongReadsMappingPath(rc chan []sam.Record, wc chan []constructdbg.DBG_MAX_INT, edgesArr []constructdbg.DBGEdge, nodesArr []constructdbg.DBGNode, opt Options) {
	fragmentFreq := make([]int, 100)
	for {
		rArr := <-rc
		if len(rArr) <= 0 {
			var guardPath []constructdbg.DBG_MAX_INT
			wc <- guardPath
			break
		}

		// process sam.Record array
		NM := sam.Tag{'N', 'M'}
		MD := sam.Tag{'M', 'D'}
		//SA := sam.Tag{'S', 'A'}
		//SN := sam.Tag{'S', 'N'}
		var arr []LRRecord
		for _, R := range rArr {
			var lrd LRRecord
			//fmt.Printf("[paraFindLongReadsMappingPath] R.Ref: %v\tR.Pos: %v\tR.Start: %v\n", R.Ref, R.Pos, R.Start())
			nm := GetAuxUint(R.AuxFields.Get(NM).Value())
			lrd.Strand = (R.Strand() < 0)
			lrd.Mch, lrd.Ins, lrd.Del, lrd.Pos = AccumulateCigar(R.Cigar, lrd.Strand)
			lrd.Mis = nm - lrd.Ins - lrd.Del
			lrd.Mch -= lrd.Mis
			lrd.RefPos = R.Start()
			lrd.RefLen = R.Ref.Len()
			id, err := strconv.Atoi(strings.Split(R.Name, "/")[0])
			if err != nil {
				log.Fatalf("[paraFindLongReadsMappingPath] R.Name: %v convert to int error: %v\n", R.Name, err)
			}
			lrd.ID = constructdbg.DBG_MAX_INT(id)
			lrd.CIGAR = make(sam.Cigar, len(R.Cigar))
			copy(lrd.CIGAR, R.Cigar)
			//lrd.CIGAR = R.Cigar
			//fmt.Printf("lrd.CIGAR: %v\n\tR.Cigar: %v\n", lrd.CIGAR, R.Cigar)
			lrd.MD = R.AuxFields.Get(MD).Value().(string)
			lrd.MapQ = uint8(R.MapQ)
			lrd.Idty = float32(lrd.Mch) / float32(lrd.Mch+lrd.Mis+lrd.Ins+lrd.Del)

			// choose the most match length
			if len(arr) > 0 && arr[len(arr)-1].ID == lrd.ID {
				if arr[len(arr)-1].Mch < lrd.Mch {
					arr[len(arr)-1] = lrd
				}
			} else {
				if lrd.Mch+lrd.Mis+lrd.Ins >= opt.Kmer {
					arr = append(arr, lrd)
				}
			}

			/*sa := R.AuxFields.Get(SA)
			if sa == nil {
				continue
			}
			//sa := R.AuxFields.Get(SA).Value().(string)
			//fmt.Printf("[paraFindLongReadsMappingPath]sa: %v\n", sa.Value().(string))
			for _, rd := range strings.Split(sa.Value().(string), ";") {
				if len(rd) == 0 {
					continue
				}
				a := strings.Split(rd, ",")
				var t LRRecord
				id, err := strconv.Atoi(a[0])
				if err != nil {
					log.Fatalf("[paraFindLongReadsMappingPath] a: %v, Ref name: %v convert to int error: %v\n", a, a[0], err)
				}
				t.RefID = constructdbg.DBG_MAX_INT(id)
				if a[2] == "-" {
					t.Strand = true
				}

				t.CIGAR, err = sam.ParseCigar([]byte(a[3]))
				if err != nil {
					log.Fatalf("[paraFindLongReadsMappingPath] ParseCigar : %v error: %v\n", a[3], err)
				}

				nm, err := strconv.Atoi(a[5])
				if err != nil {
					log.Fatalf("[paraFindLongReadsMappingPath] NM : %v convert to int error: %v\n", a[5], err)
				}
				t.Mch, t.Ins, t.Del, t.Pos = AccumulateCigar(t.CIGAR, t.Strand)
				t.Mis = nm - t.Ins - t.Del
				t.Mch -= t.Mis

				RefPos, err := strconv.Atoi(a[1])
				if err != nil {
					log.Fatalf("[paraFindLongReadsMappingPath] Ref start position : %v convert to int error: %v\n", a[1], err)
				}
				if t.Strand {
					t.RefPos = RefPos - 1 + (t.Mch + t.Mis + t.Del) - 1
				} else {
					t.RefPos = RefPos - 1
				}
				fmt.Printf("[paraFindLongReadsMappingPath]a: %v\n", a)
				arr = append(arr, t)
			}*/
		}

		//sort.Sort(LRRecordArr(arr))
		//IDArr := GetRefIDArr(arr)

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

		// clean not whole length match
		i, j := 0, 0
		for ; i < len(arr) && j < len(arr); i++ {
			qml := arr[i].GetQryConsumeLen()
			maxQryLen := opt.MaxMapEdgeLen / 2
			if edgesArr[arr[i].ID].GetSeqLen() < opt.MaxMapEdgeLen {
				maxQryLen = edgesArr[arr[i].ID].GetSeqLen()
			}
			if qml >= opt.Kmer && maxQryLen-qml < 20 {
				arr[j] = arr[i]
				j++
			} else if (i == 0 || i == len(arr)-1) && qml > 2*opt.Kmer {
				if (i == 0 && arr[i].RefPos < 20) || (i == len(arr)-1 && arr[i].RefLen-(arr[i].RefPos+arr[i].GetRefConsumeLen()) < 20) {
					arr[j] = arr[i]
					j++
				}
			}
		}
		arr = arr[:j]

		if len(arr) <= 1 {
			continue
		}

		// debug
		fmt.Printf("[paraFindLongReadsMappingPath] read name: %v\n", rArr[0].Ref.Name())
		for _, rd := range arr {
			fmt.Printf("[paraFindLongReadsMappingPath] rd: %v\n", rd)
		}
		fmt.Printf("[paraFindLongReadsMappingPath] read name: %v\n", rArr[0].Ref.Name())

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
