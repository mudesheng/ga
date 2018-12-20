package deconstructdbg

import (
	"fmt"
	"log"
	"os"
	"reflect"
	"runtime/pprof"
	"sort"
	"strconv"

	"github.com/biogo/biogo/alphabet"
	"github.com/biogo/biogo/io/seqio/fasta"
	"github.com/biogo/biogo/seq/linear"
	"github.com/jwaldrip/odin/cli"
	"github.com/mudesheng/ga/constructdbg"
	"github.com/mudesheng/ga/utils"
)

func addPathToPathMat(edgesArr []constructdbg.DBGEdge, eID constructdbg.DBG_MAX_INT, path constructdbg.Path) bool {
	for i := 1; i < len(edgesArr[eID].PathMat); i++ {
		p := edgesArr[eID].PathMat[i]
		if reflect.DeepEqual(p.IDArr, path.IDArr) {
			edgesArr[eID].PathMat[i].Freq += path.Freq
			return true
		} else if reflect.DeepEqual(p.IDArr, constructdbg.GetReverseDBG_MAX_INTArr(path.IDArr)) {
			edgesArr[eID].PathMat[i].Freq += path.Freq
			return true
		}
	}
	edgesArr[eID].PathMat = append(edgesArr[eID].PathMat, path)
	return false
}

func WriteLongPathToDBG(wc chan [2][]constructdbg.DBG_MAX_INT, edgesArr []constructdbg.DBGEdge, numCPU int) [][2][]constructdbg.DBG_MAX_INT {
	var finishNum int
	var pathArr [][2][]constructdbg.DBG_MAX_INT
	for {
		extArr := <-wc
		fmt.Fprintf(os.Stderr, "[WriteLongPathToDBG] path: %v\n", extArr)
		p := extArr[0]
		if len(p) == 0 {
			finishNum++
			if finishNum == numCPU {
				break
			}
		}
		pathArr = append(pathArr, extArr)
		// write path to the DBG
		// all long reads path store this pathArr, and store index of pathArr to the edge PathMat[1]
		/*if len(p) > 2 {
			var path constructdbg.Path
			path.IDArr = p
			path.Freq = 1
			//pathArr = append(pathArr, p)
			//idx := constructdbg.DBG_MAX_INT(len(pathArr) - 1)
			for _, eID := range p {
				if edgesArr[eID].GetUniqueFlag() > 0 {
					if len(edgesArr[eID].PathMat) == 0 {
						var t constructdbg.Path
						edgesArr[eID].PathMat = append(edgesArr[eID].PathMat, t)
					}
					addPathToPathMat(edgesArr, eID, path)
					//edgesArr[eID].PathMat[1].IDArr = append(edgesArr[eID].PathMat[1].IDArr, idx)
				}
			}
		}*/
	}

	return pathArr
}

// set the unique edge of edgesArr
/*func setDBGEdgesUniqueFlag(edgesArr []constructdbg.DBGEdge, nodesArr []constructdbg.DBGNode) (uniqueNum int) {
	for i, e := range edgesArr {
		unique := true
		if e.ID == 0 || e.GetDeleteFlag() > 0 {
			continue
		}

		if e.StartNID > 0 {
			nd := nodesArr[e.StartNID]
			ea := constructdbg.GetNearEdgeIDArr(nd, e.ID)
			if len(ea) > 1 {
				unique = false
			}
		}

		if unique == true && e.EndNID > 0 {
			nd := nodesArr[e.EndNID]
			ea := constructdbg.GetNearEdgeIDArr(nd, e.ID)
			if len(ea) > 1 {
				unique = false
			}
		}

		if unique {
			edgesArr[i].SetUniqueFlag()
			uniqueNum++
		}
	}

	return uniqueNum
}*/

type IdxLen struct {
	Idx    constructdbg.DBG_MAX_INT
	Length int
}

type IdxLenArr []IdxLen

func (a IdxLenArr) Len() int {
	return len(a)
}

func (a IdxLenArr) Less(i, j int) bool {
	return a[i].Length > a[j].Length
}

func (a IdxLenArr) Swap(i, j int) {
	a[i], a[j] = a[j], a[i]
}

func sortIdxByUniqueEdgeLen(edgesArr []constructdbg.DBGEdge) (sortedEIDIdxArr []IdxLen) {
	for _, e := range edgesArr {
		if e.GetUniqueFlag() > 0 || e.GetTwoEdgesCycleFlag() == 0 {
			var il IdxLen
			il.Idx = e.ID
			il.Length = e.GetSeqLen()
			sortedEIDIdxArr = append(sortedEIDIdxArr, il)
		}
	}

	sort.Sort(IdxLenArr(sortedEIDIdxArr))
	return
}

/*func mergePathMat(edgesArr []constructdbg.DBGEdge, nodesArr []constructdbg.DBGNode) {
	for i, e := range edgesArr {
		if e.ID == 0 || e.GetDeleteFlag() > 0 || len(e.PathMat) == 0 {
			continue
		}

		if len(e.PathMat[0].IDArr) > 0 {
			var n constructdbg.Path
			edgesArr[i].PathMat[0] = n
		}
		if len(e.PathMat) <= 1 {
			continue
		}

		var consisPM [2]constructdbg.Path
		var SecondEID constructdbg.DBG_MAX_INT // if a two edges cycle, store e.EndNID second Edge
		if e.StartNID > 0 {
			ea := constructdbg.GetNearEdgeIDArr(nodesArr[e.StartNID], e.ID)
			if len(ea) == 1 {
				// check if a two edges cycle
				if e.EndNID > 0 {
					ea1 := constructdbg.GetNearEdgeIDArr(nodesArr[e.EndNID], e.ID)
					if len(ea1) == 1 && ea1[0] == ea[0] {
						a1 := constructdbg.GetOtherEArr(nodesArr[e.StartNID], e.ID)
						a2 := constructdbg.GetOtherEArr(nodesArr[e.EndNID], e.ID)
						if len(a1) == 1 && len(a2) == 1 {
							ea[0] = a2[0]
							SecondEID = a1[0]
							fmt.Printf("[mergePathMat]two edges cycle, a1: %v, a2: %v\n", a1, a2)
							//for j := 0; j < len(e.PathMat); j++ {
							//	fmt.Printf("[mergePathMat] edgesArr[%v].PathMat[%v]: %v\n", i, j, e.PathMat[j])
							//}
						}
					}
				}

				ca := constructdbg.FindConsisPath(ea[0], e)
				if len(ca.IDArr) > 1 {
					consisPM[0] = ca
				}
			}
		}

		if e.EndNID > 0 {
			ea := constructdbg.GetNearEdgeIDArr(nodesArr[e.EndNID], e.ID)
			if len(ea) == 1 {
				if SecondEID > 0 {
					ea[0] = SecondEID
				}
				ca := constructdbg.FindConsisPath(ea[0], e)
				if len(ca.IDArr) > 1 {
					consisPM[1] = ca
				}
			}
		}

		var mP constructdbg.Path // merge Path
		if len(consisPM[0].IDArr) > 1 {
			mP.IDArr = constructdbg.GetReverseDBG_MAX_INTArr(consisPM[0].IDArr)
			mP.Freq = consisPM[0].Freq
		}
		if len(consisPM[1].IDArr) > 1 {
			if len(mP.IDArr) > 0 {
				mP.IDArr = append(mP.IDArr, consisPM[1].IDArr[1:]...)
				if mP.Freq > consisPM[1].Freq {
					mP.Freq = consisPM[1].Freq
				}
			} else {
				mP = consisPM[1]
			}
		}
		if len(mP.IDArr) > 2 {
			edgesArr[i].PathMat = []constructdbg.Path{mP}
		} else {
			edgesArr[i].PathMat = nil
		}
		//fmt.Printf("[mergePathMat] mergePath[%v]: %v\nedge: %v\n", i, edgesArr[i].PathMat, e)

	}
} */

func CountNumEdgeIDInDBG_MAX_INTArr(arr []constructdbg.DBG_MAX_INT) (num int) {
	for _, eID := range arr {
		if eID > 0 {
			num++
		}
	}
	return
}

func reallocExtendPArr(extendPArr [][]constructdbg.DBG_MAX_INT, width int, direction uint8) {
	for i, arr := range extendPArr {
		na := make([]constructdbg.DBG_MAX_INT, width)
		if direction == constructdbg.FORWARD {
			copy(na, arr)
		} else {
			x := len(arr)
			copy(na[width-x:], arr)
		}
		extendPArr[i] = na
	}
}

func ChangeExtendPArr(extendPArr [][]constructdbg.DBG_MAX_INT, deleteCol int, width int, direction uint8) {
	if len(extendPArr) == 0 {
		return
	}
	if direction == constructdbg.FORWARD {
		for i, _ := range extendPArr {
			copy(extendPArr[i], extendPArr[i][deleteCol:])
		}
	} else {
		for i, _ := range extendPArr {
			copy(extendPArr[i][deleteCol:], extendPArr[i][:width-deleteCol])
		}
	}
}

func AppearEdgeInPathNum(arr []constructdbg.DBG_MAX_INT, eID constructdbg.DBG_MAX_INT) (count int) {
	if eID < 2 {
		fmt.Printf("[AppearEdgeInPathNum] eID: %v error, must >=2\n", eID)
	}
	for _, id := range arr {
		if id == eID {
			count++
		}
	}
	return
}

func ExtendPath(edgesArr []constructdbg.DBGEdge, nodesArr []constructdbg.DBGNode, e constructdbg.DBGEdge, pathArr [][2][]constructdbg.DBG_MAX_INT, edgesPathRelationArr [][]uint32, minMapingFreq int32, kmerlen int) constructdbg.Path {
	var p constructdbg.Path
	//minMapEdgesNum := 3
	p.IDArr = make([]constructdbg.DBG_MAX_INT, len(e.PathMat[0].IDArr))
	copy(p.IDArr, e.PathMat[0].IDArr)
	p.Freq = e.PathMat[0].Freq
	if int32(p.Freq) < minMapingFreq {
		return p
	}

	{
		idx := constructdbg.IndexEID(p.IDArr, e.ID)
		// found left partition path
		for i := idx - 1; i >= 0; i-- {
			e2 := edgesArr[p.IDArr[i]]
			//fmt.Printf("[ExtendPath] edgesArr[148] process Flag: %v, e2.ID: %d, i: %d, idx: %d\n", edgesArr[148].GetProcessFlag(), e2.ID, i, idx)
			if e2.GetUniqueFlag() == 0 || e2.GetDeleteFlag() > 0 || e2.GetTwoEdgesCycleFlag() > 0 || len(e2.PathMat) == 0 {
				continue
			}
			if e2.GetProcessFlag() > 0 {
				log.Fatalf("[ExtendPath] found e2.ID: %v have been processed in p: %v\n", e2.ID, p)
			}
			//if e2.GetTwoEdgesCycleFlag() > 0 {
			//	step = 2
			//}
			Idx := -1
			var matchPath constructdbg.Path
			var matchLen int
			matchPathNum := 0
			for j := 0; j < len(e2.PathMat); j++ {
				var ep constructdbg.Path
				ep.IDArr = make([]constructdbg.DBG_MAX_INT, len(e2.PathMat[j].IDArr))
				copy(ep.IDArr, e2.PathMat[j].IDArr)
				ep.Freq = e2.PathMat[j].Freq
				if AppearEdgeInPathNum(ep.IDArr, e2.ID) > 1 {
					fmt.Printf("[ExtendPath]e2.ID: %v appear more than one time in ep: %v\n", e2.ID, ep)
					continue
				}
				x := constructdbg.IndexEID(ep.IDArr, e2.ID)
				if x < len(ep.IDArr)-1 && ep.IDArr[x+1] != p.IDArr[i+1] {
					constructdbg.ReverseDBG_MAX_INTArr(ep.IDArr)
					x = len(ep.IDArr) - 1 - x
				}
				if (len(ep.IDArr)-x > len(p.IDArr)-i) || x <= i {
					continue
				}
				if reflect.DeepEqual(p.IDArr[:i+(len(ep.IDArr)-x)], ep.IDArr[x-i:]) {
					matchPathNum++
					matchPath = ep
					for y := 0; y < x-i; y++ {
						matchLen += (len(edgesArr[ep.IDArr[y]].Utg.Ks) - (kmerlen - 1))
					}
					matchLen += (kmerlen - 1)
					Idx = x - i
				}
			}
			if matchPathNum == 1 && matchLen > 2000 {
				p.IDArr = append(matchPath.IDArr[:Idx], p.IDArr...)
				if p.Freq > matchPath.Freq {
					p.Freq = matchPath.Freq
				}
				i += Idx
			}
			//ChangeExtendPArr(extendPArr, width-1-x, width, constructdbg.BACKWARD)
			fmt.Printf("[ExtendPath]extend left path p: %v\n", p)
		}
	}

	// add right partion
	{
		idx := constructdbg.IndexEID(p.IDArr, e.ID)
		for i := idx + 1; i < len(p.IDArr); i++ {
			e2 := edgesArr[p.IDArr[i]]
			if e2.GetUniqueFlag() == 0 || e2.GetDeleteFlag() > 0 || e2.GetTwoEdgesCycleFlag() > 0 || len(e2.PathMat) == 0 {
				continue
			}
			if e2.GetProcessFlag() > 0 {
				log.Fatalf("[ExtendPath] found e2.ID: %v have been processed in p: %v\n", e2.ID, p)
			}
			Idx := -1
			var matchPath constructdbg.Path
			var matchLen int
			matchPathNum := 0
			for j := 0; j < len(e2.PathMat); j++ {
				var ep constructdbg.Path
				ep.IDArr = make([]constructdbg.DBG_MAX_INT, len(e2.PathMat[j].IDArr))
				copy(ep.IDArr, e2.PathMat[j].IDArr)
				ep.Freq = e2.PathMat[j].Freq
				if AppearEdgeInPathNum(ep.IDArr, e2.ID) > 1 {
					fmt.Printf("[ExtendPath]e2.ID: %v appear more than one time in ep: %v\n", e2.ID, ep)
					continue
				}
				x := constructdbg.IndexEID(ep.IDArr, e2.ID)
				if x > 0 && ep.IDArr[x-1] != p.IDArr[i-1] {
					constructdbg.ReverseDBG_MAX_INTArr(ep.IDArr)
					x = len(ep.IDArr) - 1 - x
				}
				if (x > i) || (len(ep.IDArr)-x <= len(p.IDArr)-i) {
					continue
				}
				if reflect.DeepEqual(p.IDArr[i-x:], ep.IDArr[:x+(len(p.IDArr)-i)]) {
					matchPathNum++
					matchPath = ep
					for y := x + (len(p.IDArr) - i); y < len(ep.IDArr); y++ {
						matchLen += (len(edgesArr[ep.IDArr[y]].Utg.Ks) - (kmerlen - 1))
					}
					matchLen += (kmerlen - 1)
					Idx = x + (len(p.IDArr) - i)
				}
			}
			if matchPathNum == 1 && matchLen > 2000 {
				p.IDArr = append(p.IDArr, matchPath.IDArr[Idx:]...)
				if p.Freq > matchPath.Freq {
					p.Freq = matchPath.Freq
				}
				i += Idx
			}
			//ChangeExtendPArr(extendPArr, width-1-x, width, constructdbg.BACKWARD)
			fmt.Printf("[ExtendPath]after extend rigth path p: %v\n", p)
		}
	}

	return p
}

func findMaxPath(sortedEIDIdxArr []IdxLen, edgesArr []constructdbg.DBGEdge, nodesArr []constructdbg.DBGNode, pathArr [][2][]constructdbg.DBG_MAX_INT, edgesPathRelationArr [][]uint32, minMapingFreq int, kmerlen int) (pA []constructdbg.Path) {
	for _, item := range sortedEIDIdxArr {
		e := edgesArr[item.Idx]
		fmt.Printf("[findMaxPath]ProcFlag:%v, delFlag:%v, uniqFlag: %v,TwoEdgesCycleFlag:%v, len(e.PathMat): %v, eID: %v, length: %v\n", e.GetProcessFlag(), e.GetDeleteFlag(), e.GetUniqueFlag(), e.GetTwoEdgesCycleFlag(), len(e.PathMat), item.Idx, item.Length)
		fmt.Printf("[findMaxPath] e.PathMat: %v\n", e.PathMat)
		if e.GetDeleteFlag() > 0 || e.GetProcessFlag() > 0 || e.GetUniqueFlag() == 0 || e.GetTwoEdgesCycleFlag() > 0 || len(e.PathMat) != 1 || len(e.PathMat[0].IDArr) == 0 {
			continue
		}
		edgesArr[e.ID].SetProcessFlag()
		fmt.Printf("[findMaxPath] eID: %v, length: %v\n", item.Idx, item.Length)
		maxP := ExtendPath(edgesArr, nodesArr, e, pathArr, edgesPathRelationArr, int32(minMapingFreq), kmerlen)
		if len(maxP.IDArr) > 1 {
			pA = append(pA, maxP)
		}
	}
	return pA
}

func ConstructEdgesPathRelationship(edgesArr []constructdbg.DBGEdge, pathArr [][2][]constructdbg.DBG_MAX_INT, edgesPathRelationArr [][]uint32) {
	for i, extArr := range pathArr {
		ui := uint32(i)
		for j, eID := range extArr[0] {
			e := edgesArr[eID]
			if e.GetUniqueFlag() > 0 {
				edgesPathRelationArr[eID] = append(edgesPathRelationArr[eID], ui)
				if extArr[1][j] > 0 {
					aeID := extArr[1][j]
					edgesPathRelationArr[aeID] = append(edgesPathRelationArr[aeID], ui)
				}
			}
		}
	}
}

func cleanNGSPath(edgesArr []constructdbg.DBGEdge) {
	for i, e := range edgesArr {
		if i < 2 || e.GetDeleteFlag() > 0 {
			continue
		}
		edgesArr[i].PathMat = nil
	}
}

func copyPathArr(pathArr [][2][]constructdbg.DBG_MAX_INT, relationArr []uint32) [][2][]constructdbg.DBG_MAX_INT {
	subArr := make([][2][]constructdbg.DBG_MAX_INT, len(relationArr))
	for i, idx := range relationArr {
		if pathArr[idx][0] == nil {
			continue
		}
		srcArr := pathArr[idx]
		var extArr [2][]constructdbg.DBG_MAX_INT
		extArr[0] = make([]constructdbg.DBG_MAX_INT, len(srcArr[0]))
		extArr[1] = make([]constructdbg.DBG_MAX_INT, len(srcArr[1]))
		copy(extArr[0], srcArr[0])
		copy(extArr[1], srcArr[1])
		subArr[i] = extArr
	}
	return subArr
}

func IndexEID(extArr [2][]constructdbg.DBG_MAX_INT, eID constructdbg.DBG_MAX_INT) int {
	for i, id := range extArr[0] {
		if id == eID || extArr[1][i] == eID {
			return i
		}
	}
	return -1
}

func checkEdgePathArrDirection(edgesArr []constructdbg.DBGEdge, nodesArr []constructdbg.DBGNode, ePathArr [][2][]constructdbg.DBG_MAX_INT, eID constructdbg.DBG_MAX_INT) bool {
	e := edgesArr[eID]
	if (e.StartNID == e.EndNID) || (e.GetTwoEdgesCycleFlag() > 0) {
		log.Fatalf("[checkEdgePathArrDirection] encounter cycle edge ID: %v\n", e.ID)
	}
	var ea1, ea2 []constructdbg.DBG_MAX_INT
	step := 1
	if e.StartNID > 0 {
		if constructdbg.IsInComing(nodesArr[e.StartNID].EdgeIDIncoming, e.ID) {
			ea1 = constructdbg.GetNearEdgeIDArr(nodesArr[e.StartNID], e.ID, true)
		} else {
			ea1 = constructdbg.GetNearEdgeIDArr(nodesArr[e.StartNID], e.ID, false)
		}
	}
	if e.EndNID > 0 {
		if constructdbg.IsInComing(nodesArr[e.EndNID].EdgeIDIncoming, e.ID) {
			ea2 = constructdbg.GetNearEdgeIDArr(nodesArr[e.EndNID], e.ID, true)
		} else {
			ea2 = constructdbg.GetNearEdgeIDArr(nodesArr[e.EndNID], e.ID, false)
		}
	}

	if constructdbg.IsIntersection(ea1, ea2) {
		edgesArr[eID].ResetUniqueFlag()
		return false
	}

	for i, extArr := range ePathArr {
		idx := IndexEID(extArr, eID)
		if idx < len(extArr[0])-step {
			if constructdbg.IsInDBG_MAX_INTArr(ea1, extArr[0][idx+step]) || constructdbg.IsInDBG_MAX_INTArr(ea1, extArr[1][idx+step]) {
				constructdbg.ReverseDBG_MAX_INTArr(ePathArr[i][0])
				constructdbg.ReverseDBG_MAX_INTArr(ePathArr[i][1])
			}
		} else {
			//fmt.Printf("[checkEdgePathArrDirection] idx: %v, step: %v, extArr: %v, ea1: %v, ea2: %v\n", idx, step, extArr, ea1, ea2)
			if idx-step >= 0 && (constructdbg.IsInDBG_MAX_INTArr(ea2, extArr[0][idx-step]) || constructdbg.IsInDBG_MAX_INTArr(ea2, extArr[1][idx-step])) {
				constructdbg.ReverseDBG_MAX_INTArr(ePathArr[i][0])
				constructdbg.ReverseDBG_MAX_INTArr(ePathArr[i][1])
			}
		}
	}
	return true
}

type IDFreq struct {
	ID   constructdbg.DBG_MAX_INT
	Freq int32
}

type IDFreqArr []IDFreq

func (arr IDFreqArr) Len() int {
	return len(arr)
}

func (arr IDFreqArr) Less(i, j int) bool {
	return arr[i].Freq > arr[j].Freq
}

func (arr IDFreqArr) Swap(i, j int) {
	arr[i], arr[j] = arr[j], arr[i]
}

func AddIDFreqArr(idFreqArr []IDFreq, eID constructdbg.DBG_MAX_INT, freq int) []IDFreq {
	if eID == 0 {
		return idFreqArr
	}
	add := false
	for j, idFreq := range idFreqArr {
		if idFreq.ID == eID {
			idFreqArr[j].Freq += int32(freq)
			add = true
			break
		}
	}
	if add {
		return idFreqArr
	}
	var tmp IDFreq
	tmp.ID, tmp.Freq = eID, int32(freq)
	idFreqArr = append(idFreqArr, tmp)

	return idFreqArr
}

/*func AddIDFreqArr(idFreqArr []IDFreq, ea []constructdbg.DBG_MAX_INT) []IDFreq {
	for _, eID := range ea {
		if eID == 0 {
			continue
		}
		add := false
		for j, idFreq := range idFreqArr {
			if idFreq.ID == eID {
				idFreqArr[j].Freq++
				add = true
				break
			}
		}
		if add {
			continue
		}
		var tmp IDFreq
		tmp.ID, tmp.Freq = eID, 1
		idFreqArr = append(idFreqArr, tmp)
	}

	return idFreqArr
}*/

func GetTotalFreq(idFreqArr []IDFreq) (totalFreq int32) {
	for _, idFreq := range idFreqArr {
		totalFreq += idFreq.Freq
	}
	return
}

func IsInIDFreqArr(idFreqArr []IDFreq, eID constructdbg.DBG_MAX_INT) bool {
	for _, idFeq := range idFreqArr {
		if idFeq.ID == eID {
			return true
		}
	}
	return false
}

func MergeLongEdgePathMat(e constructdbg.DBGEdge, pm []constructdbg.Path, edgesArr []constructdbg.DBGEdge, nodesArr []constructdbg.DBGNode) (mergePath constructdbg.Path) {
	var ea1, ea2 []constructdbg.DBG_MAX_INT
	if constructdbg.IsInComing(nodesArr[e.StartNID].EdgeIDIncoming, e.ID) {
		ea1 = constructdbg.GetNearEdgeIDArr(nodesArr[e.StartNID], e.ID, true)
	} else {
		ea1 = constructdbg.GetNearEdgeIDArr(nodesArr[e.StartNID], e.ID, false)
	}
	if constructdbg.IsInComing(nodesArr[e.EndNID].EdgeIDIncoming, e.ID) {
		ea2 = constructdbg.GetNearEdgeIDArr(nodesArr[e.EndNID], e.ID, true)
	} else {
		ea2 = constructdbg.GetNearEdgeIDArr(nodesArr[e.EndNID], e.ID, false)
	}
	if len(ea1) != 1 || len(ea2) != 1 {
		return
	}
	idx1 := constructdbg.IndexEID(pm[0].IDArr, e.ID)
	idx2 := constructdbg.IndexEID(pm[1].IDArr, e.ID)
	ok1, ok2 := false, false
	if idx1 == 0 {
		if pm[0].IDArr[idx1+1] == ea1[0] {
			pm[0].IDArr = constructdbg.ReverseDBG_MAX_INTArr(pm[0].IDArr)
			idx1 = len(pm[0].IDArr) - 1
			ok1 = true
		} else if pm[0].IDArr[idx1+1] == ea2[0] {
			pm[0], pm[1] = pm[1], pm[0]
			idx1, idx2 = idx2, idx1
			if idx1 == 0 {
				if pm[0].IDArr[idx1+1] == ea1[0] {
					pm[0].IDArr = constructdbg.ReverseDBG_MAX_INTArr(pm[0].IDArr)
					idx1 = len(pm[0].IDArr) - 1
					ok1 = true
				}
			} else if idx1 == len(pm[0].IDArr)-1 {
				ok1 = true
			}
		}
	} else if idx1 == len(pm[0].IDArr)-1 {
		if pm[0].IDArr[idx1-1] == ea1[0] {
			ok1 = true
		} else if pm[0].IDArr[idx1-1] == ea2[0] {
			pm[0], pm[1] = pm[1], pm[0]
			idx1, idx2 = idx2, idx1
			if idx1 == 0 {
				if pm[0].IDArr[idx1+1] == ea1[0] {
					pm[0].IDArr = constructdbg.ReverseDBG_MAX_INTArr(pm[0].IDArr)
					idx1 = len(pm[0].IDArr) - 1
					ok1 = true
				}
			} else if idx1 == len(pm[0].IDArr)-1 {
				ok1 = true
			}
		}
	} else {
		return
	}
	if ok1 {
		if idx2 == 0 {
			if pm[1].IDArr[idx2+1] == ea2[0] {
				pm[1].IDArr = constructdbg.ReverseDBG_MAX_INTArr(pm[1].IDArr)
				idx2 = len(pm[1].IDArr) - 1
				ok2 = true
			}
		} else if idx2 == len(pm[1].IDArr)-1 {
			if pm[1].IDArr[idx2-1] == ea2[0] {
				ok2 = true
			}
		}
	}

	if ok1 == true && ok2 == true {
		mergePath.IDArr = append(pm[0].IDArr, pm[1].IDArr[1:]...)
		if pm[0].Freq > pm[1].Freq {
			mergePath.Freq = pm[1].Freq
		} else {
			mergePath.Freq = pm[0].Freq
		}
	}
	return
}

type PathArrInfo struct {
	PathArr [][2][]constructdbg.DBG_MAX_INT
	Idx     int
}

func MergePathArr(edgesArr []constructdbg.DBGEdge, nodesArr []constructdbg.DBGNode, pathArr [][2][]constructdbg.DBG_MAX_INT, edgesPathRelationArr [][]uint32, minMapFreq int) {
	for i, e := range edgesArr {
		if i < 2 || e.GetDeleteFlag() > 0 || len(edgesPathRelationArr[i]) < minMapFreq || e.GetUniqueFlag() == 0 || e.GetTwoEdgesCycleFlag() > 0 {
			continue
		}

		ePathArr := copyPathArr(pathArr, edgesPathRelationArr[i])
		if !checkEdgePathArrDirection(edgesArr, nodesArr, ePathArr, e.ID) {
			continue
		}

		// merge process
		var leftMax, rightMax int
		for _, extArr := range ePathArr {
			idx := IndexEID(extArr, e.ID)
			if idx > leftMax {
				leftMax = idx
			}
			if len(extArr[0])-idx > rightMax {
				rightMax = len(extArr[0]) - idx
			}
		}
		leftMax++
		rightMax++
		al := leftMax + rightMax
		// adjust length of ePathArr
		for j, extArr := range ePathArr {
			a1 := make([]constructdbg.DBG_MAX_INT, al)
			a2 := make([]constructdbg.DBG_MAX_INT, al)
			idx := IndexEID(extArr, e.ID)
			copy(a1[leftMax-idx:], extArr[0])
			copy(a2[leftMax-idx:], extArr[1])
			ePathArr[j][0] = a1
			ePathArr[j][1] = a2
		}
		fmt.Printf("[MergePathArr]: merge eID: %v\n", e.ID)
		for x, p := range ePathArr {
			fmt.Printf("[MergePathArr]: ePathArr[%v]: %v\n", x, p)
		}

		// find consis path, allow many consis path, if Freq > minMapFreq
		//var cleanPathNum int
		//var path constructdbg.Path
		//path.Freq = math.MaxInt64

		// process left partition
		var leftStack, rightStack []PathArrInfo
		var pai PathArrInfo
		pai.PathArr = ePathArr
		pai.Idx = leftMax - 1
		leftStack = append(leftStack, pai)
		for len(leftStack) > 0 {
			p := leftStack[len(leftStack)-1]
			leftStack = leftStack[:len(leftStack)-1]

			// process just path in the right
			if p.Idx == leftMax-1 {
				var pai PathArrInfo
				pai.Idx = leftMax
				freq := 0
				for w := 0; w < len(p.PathArr); w++ {
					if p.PathArr[w][0][p.Idx] == 0 {
						pai.PathArr = append(pai.PathArr, p.PathArr[w])
						freq++
					}
				}
				if freq >= minMapFreq {
					fmt.Printf("[MergePathArr]: add rightStack pai: %v\n", pai)
					rightStack = append(rightStack, pai)
				}
			}

			for j := p.Idx; j >= 0; j-- {

				var idFreqArr []IDFreq
				for k := 0; k < len(p.PathArr); k++ {
					if p.PathArr[k][0][j] > 0 {
						if p.PathArr[k][1][j] > 0 {
							idFreqArr = AddIDFreqArr(idFreqArr, p.PathArr[k][0][j], 1)
							idFreqArr = AddIDFreqArr(idFreqArr, p.PathArr[k][1][j], 1)
						} else {
							idFreqArr = AddIDFreqArr(idFreqArr, p.PathArr[k][0][j], 2)
						}
					}
				}
				fmt.Printf("[MergePathArr]j: %v, idFreqArr: %v\n", j, idFreqArr)
				{
					var y int
					if len(idFreqArr) >= 1 {
						sort.Sort(IDFreqArr(idFreqArr))
						y = len(idFreqArr) - 1
						for ; y >= 0; y-- {
							if idFreqArr[y].Freq/2 >= int32(minMapFreq) {
								break
							}
						}
					} else {
						y = -1
					}
					// clean low freq edgeID
					for z := y + 1; z < len(idFreqArr); z++ {
						id := idFreqArr[z].ID
						for w := 0; w < len(p.PathArr); w++ {
							if p.PathArr[w][0][j] == id && p.PathArr[w][1][j] == 0 {
								for m := j; m >= 0; m-- {
									if p.PathArr[w][0][m] > 0 {
										p.PathArr[w][0][m] = 0
										p.PathArr[w][1][m] = 0
									} else {
										break
									}
								}
							} else if p.PathArr[w][0][j] == id && p.PathArr[w][1][j] > 0 {
								p.PathArr[w][0][j] = p.PathArr[w][1][j]
								p.PathArr[w][1][j] = 0
							} else if p.PathArr[w][0][j] > 0 && p.PathArr[w][1][j] == id {
								p.PathArr[w][1][j] = 0
							}
						}
					}
					if y == -1 || y >= 1 { // need stack
						if y == -1 {
							var pai PathArrInfo
							pai.Idx = leftMax
							for w := 0; w < len(p.PathArr); w++ {
								if p.PathArr[w][0][j+1] > 0 {
									pai.PathArr = append(pai.PathArr, p.PathArr[w])
								}
							}
							fmt.Printf("[MergePathArr]: add rightStack pai: %v\n", pai)
							rightStack = append(rightStack, pai)
						} else {
							for w := 0; w < y+1; w++ {
								var pai PathArrInfo
								pai.Idx = j - 1
								id := idFreqArr[w].ID
								for m := 0; m < len(p.PathArr); m++ {
									if p.PathArr[m][0][j] == id && p.PathArr[m][1][j] == 0 {
										pai.PathArr = append(pai.PathArr, p.PathArr[m])
									} else if p.PathArr[m][0][j] == id && p.PathArr[m][1][j] > 0 {
										var pa [2][]constructdbg.DBG_MAX_INT
										pa[0] = make([]constructdbg.DBG_MAX_INT, al)
										pa[1] = make([]constructdbg.DBG_MAX_INT, al)
										copy(pa[0], p.PathArr[m][0])
										copy(pa[1], p.PathArr[m][1])
										pa[1][j] = 0
										pai.PathArr = append(pai.PathArr, pa)
										p.PathArr[m][0][j] = p.PathArr[m][1][j]
										p.PathArr[m][1][j] = 0
									} else if p.PathArr[m][0][j] > 0 && p.PathArr[m][1][j] == id {
										var pa [2][]constructdbg.DBG_MAX_INT
										pa[0] = make([]constructdbg.DBG_MAX_INT, al)
										pa[1] = make([]constructdbg.DBG_MAX_INT, al)
										copy(pa[0], p.PathArr[m][0])
										copy(pa[1], p.PathArr[m][1])
										pa[0][j] = pa[1][j]
										pa[1][j] = 0
										pai.PathArr = append(pai.PathArr, pa)
										p.PathArr[m][1][j] = 0
									}
								}
								fmt.Printf("[MergePathArr]: add leftStack pai: %v\n", pai)
								leftStack = append(leftStack, pai)
							}
						}
						break
					}
				}
			}
		}

		// process right partition
		for len(rightStack) > 0 {
			p := rightStack[len(rightStack)-1]
			rightStack = rightStack[:len(rightStack)-1]

			for j := p.Idx; j < al; j++ {
				var idFreqArr []IDFreq
				for k := 0; k < len(p.PathArr); k++ {
					if p.PathArr[k][0][j] > 0 {
						if p.PathArr[k][1][j] > 0 {
							idFreqArr = AddIDFreqArr(idFreqArr, p.PathArr[k][0][j], 1)
							idFreqArr = AddIDFreqArr(idFreqArr, p.PathArr[k][1][j], 1)
						} else {
							idFreqArr = AddIDFreqArr(idFreqArr, p.PathArr[k][0][j], 2)
						}
					}
				}
				fmt.Printf("[MergePathArr]j: %v, idFreqArr: %v\n", j, idFreqArr)
				{
					var y int
					if len(idFreqArr) >= 1 {
						sort.Sort(IDFreqArr(idFreqArr))
						y = len(idFreqArr) - 1
						for ; y >= 0; y-- {
							if idFreqArr[y].Freq/2 >= int32(minMapFreq) {
								break
							}
						}
					} else {
						y = -1
					}
					// clean low freq edgeID
					for z := y + 1; z < len(idFreqArr); z++ {
						id := idFreqArr[z].ID
						for w := 0; w < len(p.PathArr); w++ {
							if p.PathArr[w][0][j] == id && p.PathArr[w][1][j] == 0 {
								for m := j; m < al; m++ {
									if p.PathArr[w][0][m] > 0 {
										p.PathArr[w][0][m] = 0
										p.PathArr[w][1][m] = 0
									} else {
										break
									}
								}
							} else if p.PathArr[w][0][j] == id && p.PathArr[w][1][j] > 0 {
								p.PathArr[w][0][j] = p.PathArr[w][1][j]
								p.PathArr[w][1][j] = 0
							} else if p.PathArr[w][0][j] > 0 && p.PathArr[w][1][j] == id {
								p.PathArr[w][1][j] = 0
							}
						}
					}

					if y == -1 || y >= 1 { // need stack
						if y == -1 {
							var path constructdbg.Path
							ok := false
							for w := 0; w < len(p.PathArr); w++ {
								if p.PathArr[w][0][j-1] > 0 {
									if !ok {
										m := j - 2
										for ; m >= 0; m-- {
											if p.PathArr[w][0][m] == 0 {
												break
											}
										}
										path.IDArr = p.PathArr[w][0][m+1 : j]
										ok = true
									}
									path.Freq++
								}
							}
							if len(path.IDArr) > 1 {
								edgesArr[i].PathMat = append(edgesArr[i].PathMat, path)
							}
						} else {
							for w := 0; w < y+1; w++ {
								var pai PathArrInfo
								pai.Idx = j + 1
								id := idFreqArr[w].ID
								for m := 0; m < len(p.PathArr); m++ {
									if p.PathArr[m][0][j] == id && p.PathArr[m][1][j] == 0 {
										pai.PathArr = append(pai.PathArr, p.PathArr[m])
									} else if p.PathArr[m][0][j] == id && p.PathArr[m][1][j] > 0 {
										var pa [2][]constructdbg.DBG_MAX_INT
										pa[0] = make([]constructdbg.DBG_MAX_INT, al)
										pa[1] = make([]constructdbg.DBG_MAX_INT, al)
										copy(pa[0], p.PathArr[m][0])
										copy(pa[1], p.PathArr[m][1])
										pa[1][j] = 0
										pai.PathArr = append(pai.PathArr, pa)
										p.PathArr[m][0][j] = p.PathArr[m][1][j]
										p.PathArr[m][1][j] = 0
									} else if p.PathArr[m][0][j] > 0 && p.PathArr[m][1][j] == id {
										var pa [2][]constructdbg.DBG_MAX_INT
										pa[0] = make([]constructdbg.DBG_MAX_INT, al)
										pa[1] = make([]constructdbg.DBG_MAX_INT, al)
										copy(pa[0], p.PathArr[m][0])
										copy(pa[1], p.PathArr[m][1])
										pa[0][j] = pa[1][j]
										pa[1][j] = 0
										pai.PathArr = append(pai.PathArr, pa)
										p.PathArr[m][1][j] = 0
									}
								}
								fmt.Printf("[MergePathArr]: add rightStack pai: %v\n", pai)
								rightStack = append(rightStack, pai)
							}
						}
						break
					}
				}

			}
		}
		for x, p := range edgesArr[i].PathMat {
			fmt.Printf("[MergePathArr]edgesArr[%v]PathMat[%v] : %v\n", i, x, p)
		}
		// if long edge, can merge
		if len(edgesArr[i].PathMat) == 2 {
			mergePath := MergeLongEdgePathMat(edgesArr[i], edgesArr[i].PathMat, edgesArr, nodesArr)
			if len(mergePath.IDArr) > 0 {
				edgesArr[i].PathMat[0] = mergePath
				edgesArr[i].PathMat = edgesArr[i].PathMat[:1]
			}
		}
	}
	/*// add left==0 and right==1 partition
		for z := 0; z < 2; z++ {
			var j, step int
			if z == 0 {
				j = leftMax
				step = -1
			} else {
				j = leftMax + 1
				step = 1
			}
			for ; ; j += step {
				if z == 0 {
					if j < 0 {
						break
					}
				} else {
					if j >= al {
						break
					}
				}
				//var id1, id2 constructdbg.DBG_MAX_INT
				//var freq1, freq2 int
				var idFreqArr []IDFreq
				var uncertainFreq int
				for k := 0; k < len(ePathArr); k++ {
					if len(ePathArr[k][0]) > 0 && ePathArr[k][0][j] > 0 {
						if ePathArr[k][1][j] > 0 {
							uncertainFreq++
						} else {
							idFreqArr = AddIDFreqArr(idFreqArr, ePathArr[k][0][j])
						}
						//fmt.Printf("[MergePathArr]k: %v, eID1:%v, eID2: %v, j: %v id1: %v, id2: %v, freq1: %v, freq2: %v\n", k, eID1, eID2, j, id1, id2, freq1, freq2)
					}
				}
				fmt.Printf("[MergePathArr] idFreqArr: %v\n", idFreqArr)
				if len(idFreqArr) > 1 {
					//totalFreq := GetTotalFreq(idFreqArr)
					sort.Sort(IDFreqArr(idFreqArr))
					if int(idFreqArr[0].Freq)+uncertainFreq < minMapFreq {
						break
					}
					if float32(idFreqArr[0].Freq)/float32(2) < float32(idFreqArr[1].Freq) {
						if IsBubble(edgesArr[idFreqArr[0].ID], edgesArr[idFreqArr[1].ID], nodesArr) && (float32(idFreqArr[0].Freq) >= float32(idFreqArr[1].Freq*2)) {
							// this is a bubble cause mapping uncertain
						} else {
							notConsis = true
							break
						}
					}

					{ // this maybe wrong mapping path, need clean up
						for x, p := range ePathArr {
							if len(p[0]) > 0 && IsInIDFreqArr(idFreqArr[1:], p[0][j]) {
								idx := edgesPathRelationArr[e.ID][x]
								if len(pathArr[idx][0]) > 0 {
									fmt.Fprintf(os.Stderr, "[MergePathArr]clean path idx: %v, path: %v\n", idx, pathArr[idx])
									pathArr[idx][0] = nil
									pathArr[idx][1] = nil
									cleanPathNum++
								}
							}
						}
					}
				} else {
					if len(idFreqArr) == 0 || idFreqArr[0].Freq+int32(uncertainFreq) < int32(minMapFreq) {
						break
					}
				}
				path.IDArr = append(path.IDArr, idFreqArr[0].ID)
				if path.Freq > int(idFreqArr[0].Freq)+uncertainFreq {
					path.Freq = int(idFreqArr[0].Freq) + uncertainFreq
				}
				fmt.Printf("[MergePathArr]j: %v,  path: %v\n", j, path)
			}
			if notConsis {
				break
			}
			if len(path.IDArr) < 1 {
				break
			}
			if z == 0 {
				path.IDArr = constructdbg.GetReverseDBG_MAX_INTArr(path.IDArr)
			}
		}
		fmt.Printf("[MergePathArr]notConsis: %v, path: %v\n", notConsis, path)

		if !notConsis && len(path.IDArr) >= 2 {
			var pm []constructdbg.Path
			pm = append(pm, path)
			edgesArr[i].PathMat = pm
			fmt.Printf("[MergePathArr]eID: %v,  pm: %v\n", i, pm)
		}
	}*/
}

func GetContigSeq(p constructdbg.Path, edgesArr []constructdbg.DBGEdge, nodesArr []constructdbg.DBGNode, kmerlen int) (seq []byte) {
	var strand bool
	e := edgesArr[p.IDArr[0]]
	e1 := edgesArr[p.IDArr[1]]
	if e.StartNID == e.EndNID {
		log.Fatalf("[GetContigSeq] in path: %v, start edge is self cycle: %v\n", p, e)
	}
	if (e.EndNID == e1.StartNID) || (e.EndNID == e1.EndNID) {
		strand = true
	} else {
		strand = false
	}
	for i, id := range p.IDArr {
		e := edgesArr[id]
		pos := 0
		if i > 0 {
			pos = kmerlen - 1
		}
		if strand {
			seq = append(seq, e.Utg.Ks[pos:]...)
		} else {
			seq = append(seq, constructdbg.GetReverseCompByteArr(e.Utg.Ks[:len(e.Utg.Ks)-pos])...)
		}
		// reset strand
		if i < len(p.IDArr)-1 {
			ne := edgesArr[p.IDArr[i+1]]
			strand = GetNextEdgeStrand(e, ne, nodesArr, strand)
		}
	}

	return
}

func ExtractSeq(edgesArr []constructdbg.DBGEdge, nodesArr []constructdbg.DBGNode, joinPathArr []constructdbg.Path, DcDBGEdgesfn string, kmerlen int) {
	fp, err := os.Create(DcDBGEdgesfn)
	if err != nil {
		log.Fatalf("[ExtractSeq] create file: %s failed, err: %v\n", DcDBGEdgesfn, err)
	}
	defer fp.Close()

	fafp := fasta.NewWriter(fp, 80)
	seqID := int(1)

	for _, p := range joinPathArr {
		if len(p.IDArr) <= 1 {
			continue
		}
		fmt.Printf("[ReconstructDBG] p: %v\n", p)
		//if IsTwoEdgeCyclePath(path) { joinPathArr[i].IDArr = }

		//e := constructdbg.CascadePath(p, edgesArr, nodesArr, kmerlen, false)
		seq := GetContigSeq(p, edgesArr, nodesArr, kmerlen)
		contig := linear.NewSeq("", alphabet.BytesToLetters(seq), alphabet.DNA)
		contig.ID = strconv.Itoa(int(seqID))
		seqID++
		// Add start and end adapter seq
		//qs := constructdbg.Transform2QSeq(e.Utg)
		//seq.AppendQLetters(qs...)
		var ps string
		for _, eID := range p.IDArr {
			ps += strconv.Itoa(int(eID)) + "-"
		}
		ps = ps[:len(ps)-1]
		ans := "path: " + ps + "\tFreq: " + strconv.Itoa(p.Freq)
		contig.Annotation.SetDescription(ans)
		_, err := fafp.Write(contig)
		if err != nil {
			log.Fatalf("[StoreEdgesToFn] write seq: %v; err: %v\n", contig, err)
		}
	}
}

func SimplifyByLongReadsPath(edgesArr []constructdbg.DBGEdge, nodesArr []constructdbg.DBGNode, pathArr [][2][]constructdbg.DBG_MAX_INT, opt Options) []constructdbg.Path {

	sortedEIDIdxArr := sortIdxByUniqueEdgeLen(edgesArr)
	edgesPathRelationArr := make([][]uint32, len(edgesArr))
	ConstructEdgesPathRelationship(edgesArr, pathArr, edgesPathRelationArr)

	cleanNGSPath(edgesArr)
	MergePathArr(edgesArr, nodesArr, pathArr, edgesPathRelationArr, opt.MinMapFreq)
	//MergePathArr(edgesArr, nodesArr, pathArr, edgesPathRelationArr, opt.MinMapFreq)
	//constructdbg.MergePathMat(edgesArr, nodesArr, 2)
	//fmt.Printf("[SimplifyByLongReadsPath] sortedEIDIdxArr: %v\n", sortedEIDIdxArr)
	mergePathArr := findMaxPath(sortedEIDIdxArr, edgesArr, nodesArr, pathArr, edgesPathRelationArr, opt.MinMapFreq, opt.Kmer)
	// constuct map edge ID to the path
	//IDMapPath := constructdbg.ConstructIDMapPath(pathArr)
	//constructdbg.DeleteJoinPathArrEnd(edgesArr, pathArr)
	return mergePathArr
}

type Options struct {
	utils.ArgsOpt
	//MaxMapEdgeLen int
	MinCov        int
	WinSize       int
	MaxNGSReadLen int
	MinMapFreq    int
	ExtLen        int
	ONTFn         string
	Correct       bool
}

func checkArgs(c cli.Command) (opt Options, succ bool) {

	/*var ok bool
	opt.TipMaxLen, ok = c.Flag("tipMaxLen").Get().(int)
	if !ok {
		log.Fatalf("[checkArgs] argument 'tipMaxLen': %v set error\n ", c.Flag("tipMaxlen").String())
	}
	//opt.TipMaxLen = tmp

	//tmp, err = strconv.Atoi(c.Flag("WinSize").String())
	//tmp = c.Flag("WinSize")
	opt.WinSize, ok = c.Flag("WinSize").Get().(int)
	if !ok {
		log.Fatalf("[checkArgs] argument 'WinSize': %v set error\n ", c.Flag("WinSize").String())
	}
	if opt.WinSize < 1 || opt.WinSize > 100 {
		log.Fatalf("[checkArgs] argument 'WinSize': %v must between 1~100\n", c.Flag("WinSize"))
	} */
	var ok bool
	/*opt.MaxMapEdgeLen, ok = c.Flag("MaxMapEdgeLen").Get().(int)
	if !ok {
		log.Fatalf("[checkArgs] argument 'MaxMapEdgeLen': %v set error\n ", c.Flag("MaxMapEdgeLen").String())
	}
	if opt.MaxMapEdgeLen < 2000 {
		log.Fatalf("[checkArgs] argument 'MaxMapEdgeLen': %v must bigger than 2000\n", c.Flag("MaxMapEdgeLen").String())
	}*/
	opt.MinCov, ok = c.Flag("MinCov").Get().(int)
	if !ok {
		log.Fatalf("[checkArgs] argument 'MinCov': %v set error\n ", c.Flag("MinCov").String())
	}
	if opt.MinCov < 2 || opt.MinCov > 10 {
		log.Fatalf("[checkArgs] argument 'MinCov': %v must between 2~10\n", c.Flag("MinCov"))
	}
	opt.WinSize, ok = c.Flag("WinSize").Get().(int)
	if !ok {
		log.Fatalf("[checkArgs] argument 'WinSize': %v set error\n ", c.Flag("WinSize").String())
	}
	if opt.WinSize < 1 || opt.WinSize > 100 {
		log.Fatalf("[checkArgs] argument 'WinSize': %v must between 1~100\n", c.Flag("WinSize"))
	}
	opt.MaxNGSReadLen, ok = c.Flag("MaxNGSReadLen").Get().(int)
	if !ok {
		log.Fatalf("[checkArgs] argument 'MaxNGSReadLen': %v set error\n ", c.Flag("MaxNGSReadLen").String())
	}
	if opt.MaxNGSReadLen < opt.Kmer+50 {
		log.Fatalf("[checkArgs] argument 'MaxNGSReadLen': %v must bigger than K+50\n", c.Flag("MaxNGSReadLen").String())
	}

	opt.MinMapFreq, ok = c.Flag("MinMapFreq").Get().(int)
	if !ok {
		log.Fatalf("[checkArgs] argument 'MinMapFreq': %v set error\n ", c.Flag("MinMapFreq").String())
	}
	if opt.MinMapFreq < 5 && opt.MinMapFreq >= 20 {
		log.Fatalf("[checkArgs] argument 'MinMapFreq': %v must 5 <= MinMapFreq < 20\n", c.Flag("MinMapFreq").String())
	}

	opt.ExtLen, ok = c.Flag("ExtLen").Get().(int)
	if !ok {
		log.Fatalf("[checkArgs] argument 'ExtLen': %v set error\n ", c.Flag("ExtLen").String())
	}
	if opt.ExtLen < 1000 {
		log.Fatalf("[checkArgs] argument 'ExtLen': %v must bigger than 1000\n", c.Flag("ExtLen").String())
	}
	opt.Correct, ok = c.Flag("Correct").Get().(bool)
	if !ok {
		log.Fatalf("[checkArgs] argument 'Correct': %v set error\n ", c.Flag("Correct").String())
	}

	opt.ONTFn = c.Flag("LongReadFile").String()

	succ = true
	return opt, succ
}

func DeconstructDBG(c cli.Command) {
	// check arguments
	gOpt, suc := utils.CheckGlobalArgs(c.Parent())
	if suc == false {
		log.Fatalf("[Smfy] check global Arguments error, opt: %v\n", gOpt)
	}

	opt := Options{gOpt, 0, 0, 0, 0, 0, "", false}
	tmp, suc := checkArgs(c)
	if suc == false {
		log.Fatalf("[Smfy] check Arguments error, opt: %v\n", tmp)
	}
	//opt.TipMaxLen = tmp.TipMaxLen
	opt.MinCov = tmp.MinCov
	opt.WinSize = tmp.WinSize
	opt.MaxNGSReadLen = tmp.MaxNGSReadLen
	opt.MinMapFreq = tmp.MinMapFreq
	//opt.MaxMapEdgeLen = tmp.MaxMapEdgeLen
	opt.ExtLen = tmp.ExtLen
	opt.ONTFn = tmp.ONTFn
	opt.Correct = tmp.Correct
	//constructdbg.Kmerlen = opt.Kmer
	fmt.Printf("Arguments: %v\n", opt)

	profileFn := opt.Prefix + ".decdbg.prof"
	cpuprofilefp, err := os.Create(profileFn)
	if err != nil {
		log.Fatalf("[CCF] open cpuprofile file: %v failed\n", profileFn)
	}
	pprof.StartCPUProfile(cpuprofilefp)
	defer pprof.StopCPUProfile()

	// read files and construt DBG
	DBGInfofn := opt.Prefix + ".smfy.DBGInfo"
	eSize, nSize := constructdbg.DBGInfoReader(DBGInfofn)
	nodesfn := opt.Prefix + ".nodes.smfy.Arr"
	nodesArr := constructdbg.NodesArrReader(nodesfn)
	if len(nodesArr) != nSize {
		log.Fatalf("[DeconstructDBG] len(nodesArr): %v != nodesArr Size: %v in file: %v\n", len(nodesArr), nSize, DBGInfofn)
	}
	edgesArr := make([]constructdbg.DBGEdge, eSize)
	edgesfn := opt.Prefix + ".edges.smfy.fq"
	constructdbg.LoadEdgesfqFromFn(edgesfn, edgesArr, true)

	constructdbg.CheckInterConnectivity(edgesArr, nodesArr)

	uniqueNum, semiUniqueNum, twoCycleNum, selfCycle := constructdbg.SetDBGEdgesUniqueFlag(edgesArr, nodesArr)
	fmt.Printf("[DeconstructDBG] unique edge number is : %v, semiUniqueNum: %v, twoCycleNum:%v, selfCycle:%v\n", uniqueNum, semiUniqueNum, twoCycleNum, selfCycle)

	// remap NGS reads to the new samplify DBG
	/*var copt constructdbg.Options
	copt.Kmer = opt.Kmer
	copt.NumCPU = opt.NumCPU
	copt.WinSize = opt.WinSize
	copt.MaxNGSReadLen = opt.MaxNGSReadLen
	copt.CfgFn = opt.CfgFn
	wrFn := opt.Prefix + ".decdbg.NGSAlignment"
	constructdbg.MapNGS2DBG(copt, nodesArr, edgesArr, wrFn)
	constructdbg.AddPathToDBGEdge(edgesArr, wrFn)
	constructdbg.MergePathMat(edgesArr, nodesArr, opt.MinMapFreq)*/

	// get ont reads mapping info by minimap2
	//paffn := opt.Prefix + ".paf"

	// get ont Long reads Mapping info by minimap2, must use ont or other Long reads as reference, and smfy edges as query
	paffn := opt.Prefix + ".paf"
	rc := make(chan []PAFInfo, opt.NumCPU)
	wc := make(chan [2][]constructdbg.DBG_MAX_INT, opt.NumCPU)

	go GetPAFRecord(paffn, opt.ONTFn, rc, opt.NumCPU)

	for i := 0; i < opt.NumCPU; i++ {
		go paraFindLongReadsMappingPath(rc, wc, edgesArr, nodesArr, opt)
	}

	pathArr := WriteLongPathToDBG(wc, edgesArr, opt.NumCPU)

	// Simplify using Long Reads Mapping info
	joinPathArr := SimplifyByLongReadsPath(edgesArr, nodesArr, pathArr, opt)

	graphfn := opt.Prefix + ".afterLR.dot"
	constructdbg.GraphvizDBGArr(nodesArr, edgesArr, graphfn)
	DcDBGEdgesfn := opt.Prefix + ".edges.DcDBG.fq"
	ExtractSeq(edgesArr, nodesArr, joinPathArr, DcDBGEdgesfn, opt.Kmer)
	//constructdbg.StoreEdgesToFn(DcDBGEdgesfn, edgesArr)
}
