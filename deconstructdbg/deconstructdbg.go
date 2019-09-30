package deconstructdbg

import (
	"bufio"
	"fmt"
	"io"
	"io/ioutil"
	"log"
	"math"
	"os"
	"reflect"
	"sort"
	"strconv"
	"time"

	"github.com/biogo/biogo/alphabet"
	"github.com/biogo/biogo/io/seqio/fasta"
	"github.com/biogo/biogo/seq/linear"
	"github.com/gogo/protobuf/proto"
	"github.com/jwaldrip/odin/cli"
	"github.com/mudesheng/ga/bnt"
	"github.com/mudesheng/ga/cbrotli"
	"github.com/mudesheng/ga/constructcf"
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

func WriteLongPathToFile(wc chan LongReadMappingInfo, pathfn string, numCPU int) {
	pathfp, err := os.Create(pathfn)
	if err != nil {
		log.Fatalf("[WriteLongPathToFile] file %s create error, err: %v\n", pathfn, err)
	}
	defer pathfp.Close()
	cbrofp := cbrotli.NewWriter(pathfp, cbrotli.WriterOptions{Quality: 1})
	defer cbrofp.Close()
	buffp := bufio.NewWriterSize(cbrofp, 1<<25) // 1<<25 == 2**25
	var finishNum int
	readNum := 0
	var pathArr constructdbg.PathArr
	//var rpArr []constructdbg.ReadPath
	for {
		rmi := <-wc
		extArr := rmi.Path
		fmt.Fprintf(os.Stderr, "[WriteLongPathToFile] path: %v\n", extArr)
		//p := extArr[0]
		if rmi.ID == 0 {
			finishNum++
			if finishNum == numCPU {
				break
			}
		}
		p0 := make([]uint32, len(extArr[0]))
		p1 := make([]uint32, len(extArr[1]))
		for j := 0; j < len(extArr[0]); j++ {
			p0[j] = uint32(extArr[0][j])
			p1[j] = uint32(extArr[1][j])
		}
		var rp constructdbg.ReadPath
		rp = constructdbg.ReadPath{
			ReadID: rmi.ID,
			Qual:   rmi.Qual,
			Path0:  p0,
			Path1:  p1,
		}
		//rpArr = append(rpArr, rp)

		pathArr.Arr = append(pathArr.Arr, &rp)
		readNum++
		//pathArr = append(pathArr, extArr)
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
	data, err := proto.Marshal(&pathArr)
	if err != nil {
		log.Fatalf("[WriteLongPathToFile] Mashal data error: %v\n", err)
	}

	buffp.Write(data)
	/*data, err := proto.Marshal(&pathArr)
	if err != nil {
		log.Fatalf("[WriteLongPathToFile] Mashal data error: %v\n", err)
	}
	buffp.Write(data)*/

	if err := buffp.Flush(); err != nil {
		log.Fatalf("[WriteLongPathToFile] failed to flush file: %s, err: %v\n", pathfn, err)
	}

	if err := cbrofp.Flush(); err != nil {
		log.Fatalf("[WriteLongPathToFile] failed to flush file: %s, err: %v\n", pathfn, err)
	}
}

func LoadLongPathFromFile(pathfn string) (readPathArr constructdbg.PathArr) {
	pathfp, err := os.Open(pathfn)
	if err != nil {
		log.Fatalf("[LoadLongPathFromFile] file %s create error, err: %v\n", pathfn, err)
	}
	defer pathfp.Close()
	cbrofp := cbrotli.NewReaderSize(pathfp, 1<<25)
	defer cbrofp.Close()
	buffp := bufio.NewReader(cbrofp) // 1<<25 == 2**25
	buf, err := ioutil.ReadAll(buffp)
	if err != nil {
		log.Fatalf("[LoadLongPathFromFile] read file %s failed, err:%v\n", pathfn, err)
	}
	var a constructdbg.PathArr
	err = proto.Unmarshal(buf, &a)
	if err != nil {
		log.Fatalf("[LoadLongPathFromFile] proto.Unmarshal() err:%v\n", err)
	}
	readPathArr = a
	/*readPathArr = make([]constructdbg.ReadPath, len(arr.Arr))
	for i, path := range arr.Arr {
		readPathArr[i] = path
		var p [2][]constructdbg.DBG_MAX_INT
		p[0] = make([]constructdbg.DBG_MAX_INT, len(path.Path0))
		p[1] = make([]constructdbg.DBG_MAX_INT, len(path.Path1))
		for j := 0; j < len(path.Path0); j++ {
			p[0][j] = constructdbg.DBG_MAX_INT(path.Path0[j])
			p[1][j] = constructdbg.DBG_MAX_INT(path.Path1[j])
		}
		pathArr[i] = p
	}*/

	return
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

func SortIdxByUniqueEdgeLen(edgesArr []constructdbg.DBGEdge) (sortedEIDIdxArr []IdxLen) {
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

func IsInPathMat(ePathArr [][2][]constructdbg.DBG_MAX_INT, matchPath [2][]constructdbg.DBG_MAX_INT) bool {
	ml := len(matchPath[0])
	for _, p := range ePathArr {
		if len(p[0]) < ml {
			continue
		}
		if EqualDBG_MAX_INTArr(p[0][:ml], matchPath[0]) {
			return true
		}
	}
	return false
}

func GetOverlapPath(matchPath [2][]constructdbg.DBG_MAX_INT, eID constructdbg.DBG_MAX_INT, readPathArr []*constructdbg.ReadPath, edgesPathRelation []uint32, direction uint8, minMapFreq int) (extendP [2][]constructdbg.DBG_MAX_INT, freq int) {
	Idx := IndexEID(matchPath, eID)
	ePathArr := CopyReadPathArr(readPathArr, edgesPathRelation)
	for i := 0; i < len(ePathArr); i++ {
		p := ePathArr[i]
		pl, ml := len(p[0]), len(matchPath[0])
		if pl < ml {
			ePathArr[i][0] = nil
			ePathArr[i][1] = nil
			continue
		}
		idx := IndexEID(p, eID)
		if Idx > 0 && idx > 0 {
			if matchPath[0][Idx-1] != p[0][idx-1] {
				constructdbg.ReverseDBG_MAX_INTArr(p[0])
				constructdbg.ReverseDBG_MAX_INTArr(p[1])
				idx = pl - 1 - idx
			}
		} else if Idx < ml-1 && idx < pl-1 {
			if matchPath[0][Idx+1] != p[0][idx+1] {
				constructdbg.ReverseDBG_MAX_INTArr(p[0])
				constructdbg.ReverseDBG_MAX_INTArr(p[1])
				idx = pl - 1 - idx
			}
		} else if Idx > 0 && idx < pl-1 {
			if matchPath[0][Idx-1] == p[0][idx+1] {
				constructdbg.ReverseDBG_MAX_INTArr(p[0])
				constructdbg.ReverseDBG_MAX_INTArr(p[1])
				idx = pl - 1 - idx
			}
		} else {
			if matchPath[0][Idx+1] == p[0][idx-1] {
				constructdbg.ReverseDBG_MAX_INTArr(p[0])
				constructdbg.ReverseDBG_MAX_INTArr(p[1])
				idx = pl - 1 - idx
			}
		}
		if idx > Idx && ml-Idx < pl-idx {
			if EqualDBG_MAX_INTArr(matchPath[0], p[0][idx-Idx:idx+(ml-Idx)]) {
				if direction == constructdbg.FORWARD {
					ePathArr[i][0] = p[0][idx-Idx:]
					ePathArr[i][1] = p[1][idx-Idx:]
					//fmt.Printf("[GetOverlapPath]Idx: %v, idx:%v,i:%v, p:%v\n", Idx, idx, i, p)
				} else {
					ePathArr[i][0] = constructdbg.ReverseDBG_MAX_INTArr(p[0][:idx+(ml-Idx)])
					ePathArr[i][1] = constructdbg.ReverseDBG_MAX_INTArr(p[1][:idx+(ml-Idx)])
					//fmt.Printf("[GetOverlapPath]Idx: %v, idx:%v,i:%v, p:%v\n", Idx, idx, i, p)
				}
			} else {
				ePathArr[i][0] = nil
				ePathArr[i][1] = nil
			}
		} else {
			ePathArr[i][0] = nil
			ePathArr[i][1] = nil
		}
	}

	/*for m := 0; m < len(ePathArr); m++ {
		if ePathArr[m][0] != nil {
			fmt.Printf("[GetOverlapPath] ePathArr[%v]:%v\n", m, ePathArr[m])
		}
	}*/

	// found consis path
	if direction == constructdbg.BACKWARD {
		matchPath[0] = constructdbg.ReverseDBG_MAX_INTArr(matchPath[0])
		matchPath[1] = constructdbg.ReverseDBG_MAX_INTArr(matchPath[1])
	}

	// find Max Freq edge path
	for j := len(matchPath[0]); ; j++ {
		var idFreqArr []IDFreq
		for k := 0; k < len(ePathArr); k++ {
			p := ePathArr[k]
			if p[0] == nil || len(p[0]) <= j {
				continue
			}
			if p[0][j] > 0 && p[1][j] == 0 {
				idFreqArr = AddIDFreqArr(idFreqArr, p[0][j], 1)
			}
		}

		if len(idFreqArr) == 0 {
			break
		}

		sort.Sort(IDFreqArr(idFreqArr))
		y := len(idFreqArr) - 1
		for ; y >= 0; y-- {
			if idFreqArr[y].Freq >= int32(minMapFreq) {
				break
			}
		}
		//fmt.Printf("[GetOverlapPath] y:%v,j: %v, idFreqArr: %v\n", y, j, idFreqArr)
		if y >= 1 {
			if idFreqArr[0].Freq > idFreqArr[1].Freq {
				y = 0
			}
			//ePathArr = CleanLowFreq(idFreqArr[y+1:], j, ePathArr, constructdbg.FORWARD)
		} else if y < 0 {
			break
		}
		if y == 0 {
			matchPath[0] = append(matchPath[0], idFreqArr[0].ID)
			matchPath[1] = append(matchPath[1], 0)
			freq = int(idFreqArr[0].Freq)
		} else {
			break
		}
	}
	//fmt.Printf("[GetOverlapPath] matchPath: %v\n", matchPath)

	for j := len(matchPath[0]); j > 0; j-- {
		matchPath[0] = matchPath[0][:j]
		matchPath[1] = matchPath[1][:j]
		if IsInPathMat(ePathArr, matchPath) {
			break
		}
	}

	/*for j := len(matchPath[0]); ; j++ {
		var idFreqArr []IDFreq
		for k := 0; k < len(ePathArr); k++ {
			p := ePathArr[k]
			if p[0] == nil || len(p[0]) <= j {
				continue
			}
			if p[0][j] > 0 && p[1][j] == 0 {
				idFreqArr = AddIDFreqArr(idFreqArr, p[0][j], 1)
			}
		}

		if len(idFreqArr) == 0 {
			break
		}

		sort.Sort(IDFreqArr(idFreqArr))
		y := len(idFreqArr) - 1
		for ; y >= 0; y-- {
			if idFreqArr[y].Freq >= int32(minMapFreq) {
				break
			}
		}
		//fmt.Printf("[GetOverlapPath] y:%v,j: %v, idFreqArr: %v\n", y, j, idFreqArr)
		if y == 1 {
			if idFreqArr[0].Freq > idFreqArr[1].Freq*2 {
				y = 0
			}
			ePathArr = CleanLowFreq(idFreqArr[y+1:], j, ePathArr, constructdbg.FORWARD)
		} else if y < 0 || y > 1 {
			break
		}
		if y == 0 {
			matchPath[0] = append(matchPath[0], idFreqArr[0].ID)
			matchPath[1] = append(matchPath[1], 0)
			freq = int(idFreqArr[0].Freq)
		} else {
			break
		}
	}*/

	if direction == constructdbg.BACKWARD {
		matchPath[0] = constructdbg.ReverseDBG_MAX_INTArr(matchPath[0])
		matchPath[1] = constructdbg.ReverseDBG_MAX_INTArr(matchPath[1])
	}
	extendP = matchPath
	return
}

func ExtendPath(edgesArr []constructdbg.DBGEdge, nodesArr []constructdbg.DBGNode, e constructdbg.DBGEdge, readPathArr []*constructdbg.ReadPath, edgesPathRelationArr [][]uint32, minMapingFreq, minMappingLen, kmerlen int) constructdbg.Path {
	var p, p1 constructdbg.Path
	//minMapEdgesNum := 3
	p.IDArr = make([]constructdbg.DBG_MAX_INT, len(e.PathMat[0].IDArr))
	p1.IDArr = make([]constructdbg.DBG_MAX_INT, len(e.PathMat1[0].IDArr))
	copy(p.IDArr, e.PathMat[0].IDArr)
	copy(p1.IDArr, e.PathMat1[0].IDArr)
	p.Freq = e.PathMat[0].Freq
	p1.Freq = e.PathMat1[0].Freq
	if p.Freq < minMapingFreq {
		return p
	}

	{
		idx := constructdbg.IndexEID(p.IDArr, e.ID)
		// found left partition path
		lastEID := constructdbg.DBG_MAX_INT(0)
		for {
			i := idx - 1
			for ; i >= 0; i-- {
				e2 := edgesArr[p.IDArr[i]]
				//fmt.Printf("[ExtendPath]e2.ID: %v\n", e2.ID)
				//fmt.Printf("[ExtendPath] edgesArr[148] process Flag: %v, e2.ID: %d, i: %d, idx: %d\n", edgesArr[148].GetProcessFlag(), e2.ID, i, idx)
				if e2.GetMergedFlag() == 0 || len(e2.PathMat) == 0 {
					continue
				}
				/*if e2.GetProcessFlag() > 0 {
					log.Fatalf("[ExtendPath] found e2.ID: %v have been processed in p: %v\n", e2.ID, p)
				}*/
				//if e2.GetTwoEdgesCycleFlag() > 0 {
				//	step = 2
				//}
				Idx := -1
				var matchPath, matchPath1 constructdbg.Path
				var matchLen int
				matchPathNum := 0
				fmt.Printf("[ExtendPath]e2.ID: %v, e2.PathMat:%v\n", e2.ID, e2.PathMat)
				for j := 0; j < len(e2.PathMat); j++ {
					var ep, ep1 constructdbg.Path
					ep.IDArr = make([]constructdbg.DBG_MAX_INT, len(e2.PathMat[j].IDArr))
					ep1.IDArr = make([]constructdbg.DBG_MAX_INT, len(e2.PathMat1[j].IDArr))
					copy(ep.IDArr, e2.PathMat[j].IDArr)
					copy(ep1.IDArr, e2.PathMat1[j].IDArr)
					ep.Freq = e2.PathMat[j].Freq
					ep1.Freq = e2.PathMat1[j].Freq
					if AppearEdgeInPathNum(ep.IDArr, e2.ID) > 1 {
						fmt.Printf("[ExtendPath]e2.ID: %v appear more than one time in ep: %v\n", e2.ID, ep)
						continue
					}
					x := constructdbg.IndexEID(ep.IDArr, e2.ID)
					if x < len(ep.IDArr)-1 && ep.IDArr[x+1] != p.IDArr[i+1] {
						constructdbg.ReverseDBG_MAX_INTArr(ep.IDArr)
						constructdbg.ReverseDBG_MAX_INTArr(ep1.IDArr)
						x = len(ep.IDArr) - 1 - x
					}
					if (len(ep.IDArr)-x > len(p.IDArr)-i) || x <= i {
						continue
					}
					if EqualDBG_MAX_INTArr(p.IDArr[:i+(len(ep.IDArr)-x)], ep.IDArr[x-i:]) {
						matchPathNum++
						matchPath = ep
						matchPath1 = ep1
						for y := x - i; y < len(ep.IDArr); y++ {
							matchLen += (len(edgesArr[ep.IDArr[y]].Utg.Ks) - (kmerlen - 1))
						}
						matchLen += (kmerlen - 1)
						Idx = x - i
					}
				}
				//fmt.Printf("[ExtendPath]Idx: %v, matchPathNum: %v, matchLen: %v\n", Idx, matchPathNum, matchLen)
				if matchPathNum == 1 && matchLen > minMappingLen {
					p.IDArr = append(matchPath.IDArr[:Idx], p.IDArr...)
					p1.IDArr = append(matchPath1.IDArr[:Idx], p1.IDArr...)
					p.Freq = matchPath.Freq
					p1.Freq = matchPath1.Freq
					i += Idx
					edgesArr[e2.ID].SetProcessFlag()
				}
				//ChangeExtendPArr(extendPArr, width-1-x, width, constructdbg.BACKWARD)
				//fmt.Printf("[ExtendPath]i: %v, extend left path p: %v\n", i, p)
			}
			if edgesArr[p.IDArr[0]].StartNID == 0 || edgesArr[p.IDArr[0]].EndNID == 0 {
				break
			}
			// if not found merged edge, maybe can found unique edge path
			{
				var matchPath [2][]constructdbg.DBG_MAX_INT
				var freq int
				var eID constructdbg.DBG_MAX_INT // anchor edge
				pathLen := kmerlen - 1
				for m := 0; m < len(p.IDArr); m++ {
					te := edgesArr[p.IDArr[m]]
					if eID < 2 && te.GetMergedFlag() == 0 && len(edgesPathRelationArr[te.ID]) > 0 {
						eID = te.ID
					}
					matchPath[0] = append(matchPath[0], p.IDArr[m])
					matchPath[1] = append(matchPath[1], p1.IDArr[m])
					pathLen += (len(te.Utg.Ks) - (kmerlen - 1))
					if pathLen > minMappingLen && eID >= 2 {
						break
					}
				}
				if pathLen < minMappingLen || pathLen > 10000 || eID < 2 || eID == lastEID {
					break
				}
				lastEID = eID
				//freq = p.Freq
				ml := len(matchPath[0])
				//fmt.Printf("[ExtendPath]BACKWARD anchor edge ID: %v, len(edgesPathRelationArr[eID]): %v, pathLen: %v, matchPath: %v\n", eID, len(edgesPathRelationArr[eID]), pathLen, matchPath)
				extendP, freq := GetOverlapPath(matchPath, eID, readPathArr, edgesPathRelationArr[eID], constructdbg.BACKWARD, minMapingFreq)
				//fmt.Printf("[ExtendPath]BACKWARD extendP: %v, freq:%v\n", extendP, freq)
				if len(extendP[0]) > ml && freq >= minMapingFreq {
					p.IDArr = append(extendP[0][:ml], p.IDArr...)
					p1.IDArr = append(extendP[1][:ml], p1.IDArr...)
					p.Freq = freq
					p1.Freq = freq
					idx = ml
				} else {
					break
				}
			}
			//fmt.Printf("[ExtendPath]BACKWARD after extend unique edge path :%v\n", p)
		}
	}

	// add right partion
	{
		idx := constructdbg.IndexEID(p.IDArr, e.ID)
		lastEID := constructdbg.DBG_MAX_INT(0)
		for {
			i := idx + 1
			for ; i < len(p.IDArr); i++ {
				e2 := edgesArr[p.IDArr[i]]
				//fmt.Printf("[ExtendPath]e2.ID: %v\n", e2.ID)
				if e2.GetMergedFlag() == 0 || len(e2.PathMat) == 0 {
					continue
				}
				/*if e2.GetProcessFlag() > 0 {
					log.Fatalf("[ExtendPath] found e2.ID: %v have been processed in p: %v\n", e2.ID, p)
				}*/
				Idx := -1
				var matchPath, matchPath1 constructdbg.Path
				var matchLen int
				matchPathNum := 0
				fmt.Printf("[ExtendPath]i: %v, e2.ID: %v, e2.PathMat:%v\n", i, e2.ID, e2.PathMat)
				for j := 0; j < len(e2.PathMat); j++ {
					var ep, ep1 constructdbg.Path
					ep.IDArr = make([]constructdbg.DBG_MAX_INT, len(e2.PathMat[j].IDArr))
					ep1.IDArr = make([]constructdbg.DBG_MAX_INT, len(e2.PathMat1[j].IDArr))
					copy(ep.IDArr, e2.PathMat[j].IDArr)
					copy(ep1.IDArr, e2.PathMat1[j].IDArr)
					ep.Freq = e2.PathMat[j].Freq
					ep1.Freq = e2.PathMat1[j].Freq
					if AppearEdgeInPathNum(ep.IDArr, e2.ID) > 1 {
						fmt.Printf("[ExtendPath]e2.ID: %v appear more than one time in ep: %v\n", e2.ID, ep)
						continue
					}
					x := constructdbg.IndexEID(ep.IDArr, e2.ID)
					if x > 0 && ep.IDArr[x-1] != p.IDArr[i-1] {
						constructdbg.ReverseDBG_MAX_INTArr(ep.IDArr)
						constructdbg.ReverseDBG_MAX_INTArr(ep1.IDArr)
						x = len(ep.IDArr) - 1 - x
					}
					if (x > i) || (len(ep.IDArr)-x <= len(p.IDArr)-i) {
						continue
					}
					if EqualDBG_MAX_INTArr(p.IDArr[i-x:], ep.IDArr[:x+(len(p.IDArr)-i)]) {
						matchPathNum++
						matchPath = ep
						matchPath1 = ep1
						for y := 0; y < x+(len(p.IDArr)-i); y++ {
							matchLen += (len(edgesArr[ep.IDArr[y]].Utg.Ks) - (kmerlen - 1))
						}
						matchLen += (kmerlen - 1)
						Idx = x + (len(p.IDArr) - i)
					}
				}
				//fmt.Printf("[ExtendPath]Idx: %v, matchPathNum: %v, matchLen: %v\n", Idx, matchPathNum, matchLen)
				if matchPathNum == 1 && matchLen > minMappingLen {
					p.IDArr = append(p.IDArr, matchPath.IDArr[Idx:]...)
					p1.IDArr = append(p1.IDArr, matchPath1.IDArr[Idx:]...)
					p.Freq = matchPath.Freq
					p1.Freq = matchPath1.Freq
					edgesArr[e2.ID].SetProcessFlag()
					//i += Idx
				}
				//ChangeExtendPArr(extendPArr, width-1-x, width, constructdbg.BACKWARD)
				//fmt.Printf("[ExtendPath]i: %v,len(p): %v, after extend rigth path p: %v\n", i, len(p.IDArr), p)
			}
			// if not found merged edge, maybe can found unique edge path
			{
				var matchPath [2][]constructdbg.DBG_MAX_INT
				var freq int
				var eID constructdbg.DBG_MAX_INT // anchor edge
				pathLen := kmerlen - 1
				for m := len(p.IDArr) - 1; m > 0; m-- {
					te := edgesArr[p.IDArr[m]]
					if eID < 2 && te.GetMergedFlag() == 0 && len(edgesPathRelationArr[te.ID]) > 0 {
						eID = te.ID
					}
					matchPath[0] = append(matchPath[0], p.IDArr[m])
					matchPath[1] = append(matchPath[1], p1.IDArr[m])
					pathLen += (len(te.Utg.Ks) - (kmerlen - 1))
					if pathLen > minMappingLen && eID >= 2 {
						break
					}
				}
				if pathLen < minMappingLen || pathLen > 10000 || eID < 2 || lastEID == eID {
					break
				}
				lastEID = eID
				matchPath[0] = constructdbg.ReverseDBG_MAX_INTArr(matchPath[0])
				matchPath[1] = constructdbg.ReverseDBG_MAX_INTArr(matchPath[1])
				ml := len(matchPath[0])
				//fmt.Printf("[ExtendPath]FORWARD anchor edge ID: %v,len(edgesPathRelationArr[eID]): %v, pathLen: %v, matchPath: %v\n", eID, len(edgesPathRelationArr[eID]), pathLen, matchPath)
				//freq = p.Freq
				extendP, freq := GetOverlapPath(matchPath, eID, readPathArr, edgesPathRelationArr[eID], constructdbg.FORWARD, minMapingFreq)
				//fmt.Printf("[ExtendPath]FORWARD extendP: %v, freq:%v\n", extendP, freq)
				if len(extendP[0]) > ml && freq >= minMapingFreq {
					p.IDArr = append(p.IDArr, extendP[0][ml:]...)
					p1.IDArr = append(p1.IDArr, extendP[1][ml:]...)
					p.Freq = freq
					p1.Freq = freq
					idx = i - 1
				} else {
					break
				}
				//fmt.Printf("[ExtendPath]FORWARD after extend unique edge path :%v\n", p)
			}
		}
	}

	return p
}

func FindMaxPath(sortedEIDIdxArr []IdxLen, edgesArr []constructdbg.DBGEdge, nodesArr []constructdbg.DBGNode, readPathArr []*constructdbg.ReadPath, edgesPathRelationArr [][]uint32, minMapingFreq, minMappingLen int, kmerlen int) (pA []constructdbg.Path) {
	for _, item := range sortedEIDIdxArr {
		e := edgesArr[item.Idx]
		if e.GetMergedFlag() == 0 || e.GetProcessFlag() > 0 || len(e.PathMat) != 1 || len(e.PathMat[0].IDArr) == 0 {
			continue
		}
		//fmt.Printf("[findMaxPath]ProcFlag:%v, delFlag:%v, uniqFlag: %v,TwoEdgesCycleFlag:%v, len(e.PathMat): %v, eID: %v, length: %v\n", e.GetProcessFlag(), e.GetDeleteFlag(), e.GetUniqueFlag(), e.GetTwoEdgesCycleFlag(), len(e.PathMat), item.Idx, item.Length)
		//fmt.Printf("[findMaxPath] e.PathMat: %v\n", e.PathMat)
		edgesArr[e.ID].SetProcessFlag()
		fmt.Printf("[findMaxPath] eID: %v, length: %v\n", item.Idx, item.Length)
		maxP := ExtendPath(edgesArr, nodesArr, e, readPathArr, edgesPathRelationArr, minMapingFreq, minMappingLen, kmerlen)
		if len(maxP.IDArr) > 1 {
			fmt.Printf("[findMaxPath] maxP: %v\n", maxP)
			pA = append(pA, maxP)
		}
	}
	return pA
}

func ConstructEdgesPathRelationship(edgesArr []constructdbg.DBGEdge, nodesArr []constructdbg.DBGNode, readPathArr []*constructdbg.ReadPath, edgesPathRelationArr [][]uint32, minLen, kmerlen int) [][]uint32 {
	for i, rp := range readPathArr {
		ui := uint32(i)
		var extArr [2][]uint32
		extArr[0], extArr[1] = rp.GetPath0(), rp.GetPath1()
		for j, eID := range extArr[0] {
			e := edgesArr[eID]
			if e.GetUniqueFlag() > 0 && extArr[1][j] == 0 {
				if IsBubbleEdge(e, nodesArr) {
					//fmt.Printf("[ConstructEdgesPathRelationship] e.ID: %v is a buble edge\n", e.ID)
					if len(e.Utg.Ks) > kmerlen*2+10 {
						edgesPathRelationArr[eID] = append(edgesPathRelationArr[eID], ui)
					}
				} else {
					if len(e.Utg.Ks) > minLen {
						edgesPathRelationArr[eID] = append(edgesPathRelationArr[eID], ui)
					}
				}
			}
		}
	}

	return edgesPathRelationArr
}

func CleanNGSPath(edgesArr []constructdbg.DBGEdge) {
	for i, e := range edgesArr {
		if i < 2 || e.GetDeleteFlag() > 0 {
			continue
		}
		edgesArr[i].PathMat = nil
	}
}

func CopyReadPathArr(readPathArr []*constructdbg.ReadPath, relationArr []uint32) [][2][]constructdbg.DBG_MAX_INT {
	subArr := make([][2][]constructdbg.DBG_MAX_INT, len(relationArr))
	for i, idx := range relationArr {
		path0, path1 := readPathArr[idx].GetPath0(), readPathArr[idx].GetPath1()
		if len(path0) == 0 {
			continue
		}
		//srcArr := pathArr[idx]
		var extArr [2][]constructdbg.DBG_MAX_INT
		extArr[0] = make([]constructdbg.DBG_MAX_INT, len(path0))
		extArr[1] = make([]constructdbg.DBG_MAX_INT, len(path1))
		for j, id := range path0 {
			extArr[0][j] = constructdbg.DBG_MAX_INT(id)
		}
		for j, id := range path1 {
			extArr[1][j] = constructdbg.DBG_MAX_INT(id)
		}
		subArr[i] = extArr
	}
	return subArr
}

func CopyPathArr(pathArr [][2][]constructdbg.DBG_MAX_INT, relationArr []uint32) [][2][]constructdbg.DBG_MAX_INT {
	subArr := make([][2][]constructdbg.DBG_MAX_INT, len(relationArr))
	for i, idx := range relationArr {
		//path0, path1 := readPathArr[idx].GetPath0(), readPathArr[idx].GetPath1()
		if len(pathArr[idx][0]) == 0 {
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

func MergeLongEdgePathMat(e constructdbg.DBGEdge, edgesArr []constructdbg.DBGEdge, nodesArr []constructdbg.DBGNode) (mergePath, mergePath1 constructdbg.Path) {
	var ea1, ea2 []constructdbg.DBG_MAX_INT
	var pm, pm1 [2]constructdbg.Path
	pm[0].IDArr = make([]constructdbg.DBG_MAX_INT, len(e.PathMat[0].IDArr))
	pm[1].IDArr = make([]constructdbg.DBG_MAX_INT, len(e.PathMat[1].IDArr))
	pm1[0].IDArr = make([]constructdbg.DBG_MAX_INT, len(e.PathMat1[0].IDArr))
	pm1[1].IDArr = make([]constructdbg.DBG_MAX_INT, len(e.PathMat1[1].IDArr))
	copy(pm[0].IDArr, e.PathMat[0].IDArr)
	copy(pm[1].IDArr, e.PathMat[1].IDArr)
	copy(pm1[0].IDArr, e.PathMat1[0].IDArr)
	copy(pm1[1].IDArr, e.PathMat1[1].IDArr)
	pm[0].Freq, pm[1].Freq = e.PathMat[0].Freq, e.PathMat[1].Freq
	pm1[0].Freq, pm1[1].Freq = e.PathMat1[0].Freq, e.PathMat1[1].Freq
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
	if idx1 == 0 && len(pm[1].IDArr)-idx2 == len(pm[0].IDArr) {
		if EqualDBG_MAX_INTArr(pm[0].IDArr, pm[1].IDArr[:len(pm[1].IDArr)-idx2]) {
			mergePath.IDArr = pm[1].IDArr
			mergePath1.IDArr = pm1[1].IDArr
			mergePath.Freq += pm[0].Freq
			mergePath1.Freq += pm1[0].Freq
		}
	} else if idx2 == 0 && len(pm[0].IDArr)-idx1 == len(pm[1].IDArr) {
		if EqualDBG_MAX_INTArr(pm[1].IDArr, pm[0].IDArr[:len(pm[0].IDArr)-idx1]) {
			mergePath.IDArr = pm[0].IDArr
			mergePath1.IDArr = pm1[0].IDArr
			mergePath.Freq += pm[1].Freq
			mergePath1.Freq += pm1[1].Freq
		}
	}
	ok1, ok2 := false, false
	if idx1 == 0 {
		if pm[0].IDArr[idx1+1] == ea1[0] {
			pm[0].IDArr = constructdbg.ReverseDBG_MAX_INTArr(pm[0].IDArr)
			pm1[0].IDArr = constructdbg.ReverseDBG_MAX_INTArr(pm1[0].IDArr)
			idx1 = len(pm[0].IDArr) - 1
			ok1 = true
		} else if pm[0].IDArr[idx1+1] == ea2[0] {
			pm[0], pm[1] = pm[1], pm[0]
			pm1[0], pm1[1] = pm1[1], pm1[0]
			idx1, idx2 = idx2, idx1
			if idx1 == 0 {
				if pm[0].IDArr[idx1+1] == ea1[0] {
					pm[0].IDArr = constructdbg.ReverseDBG_MAX_INTArr(pm[0].IDArr)
					pm1[0].IDArr = constructdbg.ReverseDBG_MAX_INTArr(pm1[0].IDArr)
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
			pm1[0], pm1[1] = pm1[1], pm1[0]
			idx1, idx2 = idx2, idx1
			if idx1 == 0 {
				if pm[0].IDArr[idx1+1] == ea1[0] {
					pm[0].IDArr = constructdbg.ReverseDBG_MAX_INTArr(pm[0].IDArr)
					pm1[0].IDArr = constructdbg.ReverseDBG_MAX_INTArr(pm1[0].IDArr)
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
				ok2 = true
			}
		} else if idx2 == len(pm[1].IDArr)-1 {
			if pm[1].IDArr[idx2-1] == ea2[0] {
				pm[1].IDArr = constructdbg.ReverseDBG_MAX_INTArr(pm[1].IDArr)
				pm1[1].IDArr = constructdbg.ReverseDBG_MAX_INTArr(pm1[1].IDArr)
				idx2 = len(pm[1].IDArr) - 1
				ok2 = true
			}
		}
	}

	if ok1 == true && ok2 == true {
		mergePath.IDArr = append(pm[0].IDArr, pm[1].IDArr[1:]...)
		mergePath1.IDArr = append(pm1[0].IDArr, pm1[1].IDArr[1:]...)
		if pm[0].Freq > pm[1].Freq {
			mergePath.Freq = pm[1].Freq
			mergePath1.Freq = pm1[1].Freq
		} else {
			mergePath.Freq = pm[0].Freq
			mergePath1.Freq = pm1[0].Freq
		}
	}
	return
}

type PathArrInfo struct {
	PathArr [][2][]constructdbg.DBG_MAX_INT
	Idx     int
}

func CleanLowFreq(idFreqArr []IDFreq, j int, pathArr [][2][]constructdbg.DBG_MAX_INT, direction uint8) [][2][]constructdbg.DBG_MAX_INT {
	for z := 0; z < len(idFreqArr); z++ {
		id := idFreqArr[z].ID
		for w := 0; w < len(pathArr); w++ {
			p := pathArr[w]
			if direction == constructdbg.BACKWARD {
				if len(p[0]) < j {
					continue
				}
				if pathArr[w][0][j] == id && pathArr[w][1][j] == 0 {
					for m := j; m >= 0; m-- {
						if pathArr[w][0][m] > 0 {
							pathArr[w][0][m] = 0
							pathArr[w][1][m] = 0
						} else {
							break
						}
					}
				} else if pathArr[w][0][j] == id && pathArr[w][1][j] > 0 {
					pathArr[w][0][j] = pathArr[w][1][j]
					pathArr[w][1][j] = 0
				} else if pathArr[w][0][j] > 0 && pathArr[w][1][j] == id {
					pathArr[w][1][j] = 0
				}
			} else { // FORWARD
				if len(p[0]) <= j {
					continue
				}
				if pathArr[w][0][j] == id && pathArr[w][1][j] == 0 {
					for m := j; m < len(pathArr[w][0]); m++ {
						if pathArr[w][0][m] > 0 {
							pathArr[w][0][m] = 0
							pathArr[w][1][m] = 0
						} else {
							break
						}
					}
				} else if pathArr[w][0][j] == id && pathArr[w][1][j] > 0 {
					pathArr[w][0][j] = pathArr[w][1][j]
					pathArr[w][1][j] = 0
				} else if pathArr[w][0][j] > 0 && pathArr[w][1][j] == id {
					pathArr[w][1][j] = 0
				}
			}
		}
	}
	return pathArr
}

func GetIDFreqArr(ePathArr [][2][]constructdbg.DBG_MAX_INT, j int) (idFreqArr []IDFreq) {
	for k := 0; k < len(ePathArr); k++ {
		p := ePathArr[k]
		if p[0][j] > 0 && p[1][j] == 0 {
			idFreqArr = AddIDFreqArr(idFreqArr, p[0][j], 1)
		}
	}
	return
}

func ReplaceEPathMat(ePathMatArr [][2][]constructdbg.DBG_MAX_INT, idx int, direction uint8, mp, oldP, newP []constructdbg.DBG_MAX_INT) {
	fmt.Printf("[ReplaceEPathMat]mp: %v\n\t\toldP: %v\n\t\tnewP: %v\n", mp, oldP, newP)

	if direction == constructdbg.BACKWARD {
		oldP = constructdbg.GetReverseDBG_MAX_INTArr(oldP)
		newP = constructdbg.GetReverseDBG_MAX_INTArr(newP)
		for i, pm := range ePathMatArr {
			if reflect.DeepEqual(pm[0][idx+1:idx+1+len(mp)], mp) && reflect.DeepEqual(pm[0][idx+1-len(oldP):idx+1], oldP) {
				bl := idx + 1 - len(newP)
				np := make([]constructdbg.DBG_MAX_INT, len(pm[0]))
				x := idx + 1 - len(oldP) - 1
				for j := bl - 1; j >= 0 && x >= 0; j-- {
					np[j] = pm[0][x]
					x--
				}
				copy(np[idx+1-len(newP):idx+1], newP)
				copy(np[idx+1:], pm[0][idx+1:])
				ePathMatArr[i][0] = np
			}
		}
	} else {
		for i, pm := range ePathMatArr {
			if reflect.DeepEqual(pm[0][idx-len(mp):idx], mp) && reflect.DeepEqual(pm[0][idx:idx+len(oldP)], oldP) {
				np := make([]constructdbg.DBG_MAX_INT, len(pm[0]))
				copy(np[:idx], pm[0][:idx])
				copy(np[idx:idx+len(newP)], newP)
				bl := idx + len(newP)
				x := idx + len(oldP)
				for j := bl; j < len(np) && x < len(pm[0]); j++ {
					np[j] = pm[0][x]
					x++
				}
				ePathMatArr[i][0] = np
			}
		}
	}
	return
}

func ReplaceEPath(pm []constructdbg.DBG_MAX_INT, idx int, direction uint8, mp, oldP, newP []constructdbg.DBG_MAX_INT) (path []constructdbg.DBG_MAX_INT) {
	fmt.Printf("[ReplaceEPath]mp: %v\n\t\toldP: %v\n\t\tnewP: %v\n", mp, oldP, newP)
	if direction == constructdbg.BACKWARD {
		oldP = constructdbg.GetReverseDBG_MAX_INTArr(oldP)
		newP = constructdbg.GetReverseDBG_MAX_INTArr(newP)
		if reflect.DeepEqual(pm[idx+1:idx+1+len(mp)], mp) && reflect.DeepEqual(pm[idx+1-len(oldP):idx+1], oldP) {
			bl := idx + 1 - len(newP)
			path = make([]constructdbg.DBG_MAX_INT, len(pm))
			x := idx + 1 - len(oldP) - 1
			for j := bl - 1; j >= 0 && x >= 0; j-- {
				path[j] = pm[x]
				x--
			}
			copy(path[idx+1-len(newP):idx+1], newP)
			copy(path[idx+1:], pm[idx+1:])
		}
	} else {
		if reflect.DeepEqual(pm[idx-len(mp):idx], mp) && reflect.DeepEqual(pm[idx:idx+len(oldP)], oldP) {
			path = make([]constructdbg.DBG_MAX_INT, len(pm))
			copy(path[:idx], pm[:idx])
			copy(path[idx:idx+len(newP)], newP)
			bl := idx + len(newP)
			x := idx + len(oldP)
			for j := bl; j < len(path) && x < len(pm); j++ {
				path[j] = pm[x]
				x++
			}
		}
	}
	fmt.Printf("[ReplaceEPath]changed path: %v\n", path)

	return
}

func GetEPathSectionPath(pm []constructdbg.DBG_MAX_INT, idx int, direction uint8, mp []constructdbg.DBG_MAX_INT) (path []constructdbg.DBG_MAX_INT) {
	if direction == constructdbg.FORWARD {
		if reflect.DeepEqual(pm[idx-len(mp):idx], mp) {
			path = make([]constructdbg.DBG_MAX_INT, len(pm[idx:]))
			copy(path, pm[idx:])
		}
	} else {
		if reflect.DeepEqual(pm[idx+1:idx+1+len(mp)], mp) {
			path = constructdbg.GetReverseDBG_MAX_INTArr(pm[:idx+1])
		}
	}
	return
}

func GetEPathMatPath(ePathMatArr [][2][]constructdbg.DBG_MAX_INT, idx int, direction uint8, mp []constructdbg.DBG_MAX_INT, eID constructdbg.DBG_MAX_INT) (path []constructdbg.DBG_MAX_INT) {
	slen := len(ePathMatArr[0][0])
	path = append(path, eID)
	if direction == constructdbg.FORWARD {
		for j := idx + 1; j < slen; j++ {
			//idFreqArr := GetIDFreqArr(ePathMatArr, j)
			var idFreqArr []IDFreq
			for _, pm := range ePathMatArr {
				if pm[0][idx] == eID {
					idFreqArr = AddIDFreqArr(idFreqArr, pm[0][j], 1)
				}
			}
			if len(idFreqArr) > 0 {
				sort.Sort(IDFreqArr(idFreqArr))
				path = append(path, idFreqArr[0].ID)
			} else {
				break
			}
		}
	} else {
		for j := idx - 1; j >= 0; j-- {
			//idFreqArr := GetIDFreqArr(ePathMatArr, j)
			var idFreqArr []IDFreq
			for _, pm := range ePathMatArr {
				if pm[0][idx] == eID {
					idFreqArr = AddIDFreqArr(idFreqArr, pm[0][j], 1)
				}
			}
			if len(idFreqArr) > 0 {
				sort.Sort(IDFreqArr(idFreqArr))
				path = append(path, idFreqArr[0].ID)
			} else {
				break
			}
		}
	}
	return
}

func AdjustEPathMat(ePathMatArr [][2][]constructdbg.DBG_MAX_INT, idx int, idFreqArr []IDFreq, direction uint8, edgesArr []constructdbg.DBGEdge, mp []constructdbg.DBG_MAX_INT, kmerlen int) {
	mp = mp[:len(mp)-1]
	if direction == constructdbg.BACKWARD {
		mp = constructdbg.GetReverseDBG_MAX_INTArr(mp)
	}
	path0 := GetEPathMatPath(ePathMatArr, idx, direction, mp, idFreqArr[0].ID)
	path1 := GetEPathMatPath(ePathMatArr, idx, direction, mp, idFreqArr[1].ID)
	fmt.Printf("[AdjustEPathMat]mp: %v\n\t\tpath0: %v\n\t\tpath1: %v\n", mp, path0, path1)
	if len(path0) > 1 && len(path1) > 1 {
		l0 := len(edgesArr[path0[0]].Utg.Ks) - (kmerlen - 1)
		//l0 += len(edgesArr[path0[1]].Utg.Ks) - (kmerlen-1)
		l1 := len(edgesArr[path1[0]].Utg.Ks) - (kmerlen - 1)
		//l1 += len(edgesArr[path1[1]].Utg.Ks) - (kmerlen-1)
		j := 1
		for i := 1; i < len(path0) && j < len(path1); i++ {
			if path0[i] == path1[j] {
				if utils.AbsInt(l0-l1) < constructdbg.Min(l0, l1)/10 {
					ReplaceEPathMat(ePathMatArr, idx, direction, mp, path1[:j], path0[:i])
					break
				}
			}
			if l1 < l0 {
				ok := false
				for ; j < len(path1); j++ {
					if path0[i] == path1[j] {
						if utils.AbsInt(l0-l1) < constructdbg.Min(l0, l1)/10 {
							ReplaceEPathMat(ePathMatArr, idx, direction, mp, path1[:j], path0[:i])
							ok = true
							break
						}
					} else if l1 > l0 {
						break
					}
					l1 += len(edgesArr[path1[j]].Utg.Ks) - (kmerlen - 1)
				}
				if ok {
					break
				}
			}
			l0 += len(edgesArr[path0[i]].Utg.Ks) - (kmerlen - 1)
		}
		// print
		for i, pm := range ePathMatArr {
			fmt.Printf("[AdjustEPathMat]ePathMatArr[%v]: %v\n", i, pm)
		}
	}
	return
}

func AdjustEPathMat2(ePathMatArr [][2][]constructdbg.DBG_MAX_INT, idx int, idFreqArr []IDFreq, direction uint8, edgesArr []constructdbg.DBGEdge, mp []constructdbg.DBG_MAX_INT, kmerlen int) {
	mp = mp[:len(mp)-1]
	if direction == constructdbg.BACKWARD {
		mp = constructdbg.GetReverseDBG_MAX_INTArr(mp)
	}
	path0 := GetEPathMatPath(ePathMatArr, idx, direction, mp, idFreqArr[0].ID)
	for i, pm := range ePathMatArr {
		path1 := GetEPathSectionPath(pm[0], idx, direction, mp)
		//fmt.Printf("[AdjustEPathMat2]\t\tpath[%v]: %v\n\t\tpath0: %v\n\t\tpath1: %v\n\t\tmp: %v\n", i, pm[0], path0, path1, mp)
		if len(path0) > 1 && len(path1) > 1 && path0[0] != path1[0] && path1[0] > 0 {
			fmt.Printf("[AdjustEPathMat2]need change path[%v]: %v\n\t\tpath0: %v\n\t\tpath1: %v\n\t\tmp: %v\n", i, pm[0], path0, path1, mp)
			l0 := len(edgesArr[path0[0]].Utg.Ks) - (kmerlen - 1)
			l1 := len(edgesArr[path1[0]].Utg.Ks) - (kmerlen - 1)
			y := 1
			for x := 1; x < len(path0) && y < len(path1); x++ {
				if path0[x] == path1[y] {
					if utils.AbsInt(l0-l1) < constructdbg.Min(l0, l1)/10 {
						ePathMatArr[i][0] = ReplaceEPath(pm[0], idx, direction, mp, path1[:y], path0[:x])
						break
					}
				}
				if l1 < l0 {
					ok := false
					for ; y < len(path1); y++ {
						if path0[x] == path1[y] {
							if utils.AbsInt(l0-l1) < constructdbg.Min(l0, l1)/10 {
								ePathMatArr[i][0] = ReplaceEPath(pm[0], idx, direction, mp, path1[:y], path0[:x])
								ok = true
								break
							}
						} else if l1 > l0 {
							break
						}
						l1 += len(edgesArr[path1[y]].Utg.Ks) - (kmerlen - 1)
					}
					if ok {
						break
					}
				}
				l0 += len(edgesArr[path0[x]].Utg.Ks) - (kmerlen - 1)
			}
		}
	}
	// print
	for i, pm := range ePathMatArr {
		fmt.Printf("[AdjustEPathMat2]ePathMatArr[%v]: %v\n", i, pm)
	}
	return
}

func MergePathArr(edgesArr []constructdbg.DBGEdge, nodesArr []constructdbg.DBGNode, readPathArr []*constructdbg.ReadPath, edgesPathRelationArr [][]uint32, minMapFreq, kmerlen int, minUniqueEdgeLen, depth, avgReadLen int) {
	for i, e := range edgesArr {
		if i < 2 || e.GetDeleteFlag() > 0 || e.GetUniqueFlag() == 0 || e.GetTwoEdgesCycleFlag() > 0 {
			continue
		}
		if IsBubbleEdge(e, nodesArr) && len(e.Utg.Ks) < kmerlen*2+5 {
			continue
		}

		// check is a unique edge in the genome
		Unique := false
		el := len(e.Utg.Ks)
		if el > minUniqueEdgeLen {
			Unique = true
		} else {
			freq := len(edgesPathRelationArr[i])
			if el > kmerlen*3 && freq < depth*7/5 { // (el / avgReadLen)*Depth = one long read cov this edge probility
				Unique = true
			}
		}
		// if Unique == true, start merge , other just add path to the edge PathMat
		if Unique == false {
			continue
		}
		// merge process
		ePathArr := CopyReadPathArr(readPathArr, edgesPathRelationArr[i])
		if !checkEdgePathArrDirection(edgesArr, nodesArr, ePathArr, e.ID) {
			continue
		}

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

		// every position just choose max freq edge
		{
			var mp []constructdbg.DBG_MAX_INT
			freq := math.MaxInt32
			// BACKWARD
			{
				for j := leftMax; j >= 0; j-- {
					idFreqArr := GetIDFreqArr(ePathArr, j)
					if len(idFreqArr) == 0 {
						break
					}
					sort.Sort(IDFreqArr(idFreqArr))
					fmt.Printf("[MergePathArr][%v] idFreqArr: %v\n", j, idFreqArr)
					//fmt.Printf("[MergePathArr]j: %v, idFreqArr: %v\n", j, idFreqArr)
					y := len(idFreqArr) - 1
					for ; y >= 0; y-- {
						if idFreqArr[y].Freq >= int32(minMapFreq) {
							break
						}
					}
					if y < 0 {
						break
					}

					if y >= 1 {
						if el > minUniqueEdgeLen && idFreqArr[0].Freq > idFreqArr[1].Freq {
							y = 0
						}
					}

					if y == 0 {
						mp = append(mp, idFreqArr[0].ID)
						if idFreqArr[0].Freq < int32(freq) {
							freq = int(idFreqArr[0].Freq)
						}
						te := edgesArr[idFreqArr[0].ID]
						// found if encounter complex bubble
						if te.GetTwoEdgesCycleFlag() == 0 && j < leftMax && j > 0 && len(idFreqArr) > 1 {
							AdjustEPathMat2(ePathArr, j, idFreqArr, constructdbg.BACKWARD, edgesArr, mp, kmerlen)
						}
					} else {
						break
					}
				}
			}
			if len(mp) >= 2 { // may found boundary uncertain
				if mp[len(mp)-1] == mp[len(mp)-2] {
					mp = mp[:len(mp)-1]
				}
			}
			mp = constructdbg.ReverseDBG_MAX_INTArr(mp)
			leftLen := len(mp) - 1
			// FORWARD
			{

				for j := leftMax + 1; j < leftMax+rightMax; j++ {
					idFreqArr := GetIDFreqArr(ePathArr, j)
					if len(idFreqArr) == 0 {
						break
					}
					sort.Sort(IDFreqArr(idFreqArr))
					fmt.Printf("[MergePathArr][%v] idFreqArr: %v\n", j, idFreqArr)
					//fmt.Printf("[MergePathArr]j: %v, idFreqArr: %v\n", j, idFreqArr)
					y := len(idFreqArr) - 1
					for ; y >= 0; y-- {
						if idFreqArr[y].Freq >= int32(minMapFreq) {
							break
						}
					}
					if y < 0 {
						break
					}

					if y >= 1 {
						if el > minUniqueEdgeLen && idFreqArr[0].Freq > idFreqArr[1].Freq {
							y = 0
						}
					}

					if y == 0 {
						mp = append(mp, idFreqArr[0].ID)
						if idFreqArr[0].Freq < int32(freq) {
							freq = int(idFreqArr[0].Freq)
						}
						te := edgesArr[idFreqArr[0].ID]
						// found if encounter complex bubble
						if te.GetTwoEdgesCycleFlag() == 0 && j < leftMax+rightMax-1 && len(idFreqArr) > 1 {
							AdjustEPathMat2(ePathArr, j, idFreqArr, constructdbg.FORWARD, edgesArr, mp[leftLen:], kmerlen)
						}
					} else {
						break
					}
				}
			}
			if len(mp) >= 2 {
				if mp[len(mp)-1] == mp[len(mp)-2] {
					mp = mp[:len(mp)-1]
				}
			}
			if len(mp) >= 3 {
				var path1, path2 constructdbg.Path
				path1.IDArr = mp
				path2.IDArr = make([]constructdbg.DBG_MAX_INT, len(mp))
				path1.Freq, path2.Freq = freq, freq
				edgesArr[i].PathMat = append(edgesArr[i].PathMat, path1)
				edgesArr[i].PathMat1 = append(edgesArr[i].PathMat1, path2)
				// clean relationship array
				edgesPathRelationArr[i] = nil
				edgesArr[i].SetMergedFlag()
				fmt.Printf("[MergePathArr]mergeed successed, edgesArr[%v]PathMat : %v, edgesArr[%v]PathMat1: %v\n", i, edgesArr[i].PathMat, i, edgesArr[i].PathMat1)
			}
		}
	}

	/*// process left partition
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
				//uncertainNum := 0
				var idFreqArr []IDFreq
				for k := 0; k < len(p.PathArr); k++ {
					if p.PathArr[k][0][j] > 0 {
						if p.PathArr[k][1][j] > 0 {
							//uncertainNum++
							//idFreqArr = AddIDFreqArr(idFreqArr, p.PathArr[k][0][j], 1)
							//idFreqArr = AddIDFreqArr(idFreqArr, p.PathArr[k][1][j], 1)
						} else {
							idFreqArr = AddIDFreqArr(idFreqArr, p.PathArr[k][0][j], 1)
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
							if idFreqArr[y].Freq >= int32(minMapFreq) {
								break
							}
						}
						// process cycle edges path, boundary uncertain problem
						if y == 1 {
							te0 := edgesArr[idFreqArr[0].ID]
							te1 := edgesArr[idFreqArr[1].ID]
							if (te0.GetTwoEdgesCycleFlag() > 0 || te0.StartNID == te0.EndNID) && te1.GetUniqueFlag() > 0 && idFreqArr[0].Freq > idFreqArr[1].Freq {
								if idFreqArr[0].Freq > idFreqArr[1].Freq*2 {

								} else {
									idFreqArr[0], idFreqArr[1] = idFreqArr[1], idFreqArr[0]
								}
								y = 0
							} else if Unique {
								if idFreqArr[0].Freq*4/5 > idFreqArr[1].Freq {
									y = 0
								}
							}
							// if is a bubble
							if IsBubble(te0, te1, nodesArr) {
								if Unique {
									if idFreqArr[0].Freq > idFreqArr[1].Freq*2 {
										y = 0
									} else { // edge allow bubble
										for k := 0; k < len(p.PathArr); k++ {
											if p.PathArr[k][0][j] > 0 {
												p.PathArr[k][0][j] = idFreqArr[0].ID
												if p.PathArr[k][1][j] > 0 {
													p.PathArr[k][1][j] = idFreqArr[1].ID
												}
											}
										}
										y = 0
										idFreqArr = idFreqArr[:1]
									}
								}
							}
						}
						// clean low freq edgeID
						p.PathArr = CleanLowFreq(idFreqArr[y+1:], j, p.PathArr, constructdbg.BACKWARD)
					} else {
						y = -1
					}

					if y == -1 || y >= 1 { // need stack
						if y == -1 {
							var pai PathArrInfo
							pai.Idx = leftMax
							for w := 0; w < len(p.PathArr); w++ {
								if p.PathArr[w][0][j+1] > 0 {
									for z := j; z >= 0; z-- {
										if p.PathArr[w][0][z] > 0 {
											p.PathArr[w][0][z] = 0
											p.PathArr[w][1][z] = 0
										} else {
											break
										}
									}
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
				//uncertainNum := 0
				for k := 0; k < len(p.PathArr); k++ {
					if p.PathArr[k][0][j] > 0 {
						if p.PathArr[k][1][j] > 0 {
							//uncertainNum++
						} else {
							idFreqArr = AddIDFreqArr(idFreqArr, p.PathArr[k][0][j], 1)
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
							if idFreqArr[y].Freq >= int32(minMapFreq) {
								break
							}
						}
						// process cycle edges path, boundary uncertain problem
						if y == 1 {
							te0 := edgesArr[idFreqArr[0].ID]
							te1 := edgesArr[idFreqArr[1].ID]
							if (te0.GetTwoEdgesCycleFlag() > 0 || te0.StartNID == te0.EndNID) && te1.GetUniqueFlag() > 0 && idFreqArr[0].Freq > idFreqArr[1].Freq {
								if idFreqArr[0].Freq > idFreqArr[1].Freq*2 {

								} else {
									idFreqArr[0], idFreqArr[1] = idFreqArr[1], idFreqArr[0]
								}
								y = 0
							} else if Unique {
								if idFreqArr[0].Freq*4/5 > idFreqArr[1].Freq {
									y = 0
								}
							}
							// if is a bubble
							if IsBubble(te0, te1, nodesArr) {
								if Unique { //avgReadLen/(avgReadLen+el*2) one long read cov this edge probility
									if idFreqArr[0].Freq > idFreqArr[1].Freq*2 {
										y = 0
									} else { // edge allow bubble
										for k := 0; k < len(p.PathArr); k++ {
											if p.PathArr[k][0][j] > 0 {
												p.PathArr[k][0][j] = idFreqArr[0].ID
												if p.PathArr[k][1][j] > 0 {
													p.PathArr[k][1][j] = idFreqArr[1].ID
												}
											}
										}
										y = 0
										idFreqArr = idFreqArr[:1]
									}
								}
							}
						}
					} else {
						y = -1
					}
					// clean low freq edgeID
					p.PathArr = CleanLowFreq(idFreqArr[y+1:], j, p.PathArr, constructdbg.FORWARD)

					if y == -1 || y >= 1 { // need stack
						if y == -1 {
							var path1, path2 constructdbg.Path
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
										path1.IDArr = p.PathArr[w][0][m+1 : j]
										path2.IDArr = p.PathArr[w][1][m+1 : j]
										ok = true
									}
									path1.Freq++
									path2.Freq++
								}
							}
							if len(path1.IDArr) > 1 {
								edgesArr[i].PathMat = append(edgesArr[i].PathMat, path1)
								edgesArr[i].PathMat1 = append(edgesArr[i].PathMat1, path2)
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
			mergePath0, mergePath1 := MergeLongEdgePathMat(edgesArr[i], edgesArr, nodesArr)
			if len(mergePath0.IDArr) > 0 {
				edgesArr[i].PathMat[0] = mergePath0
				edgesArr[i].PathMat1[0] = mergePath1
				edgesArr[i].PathMat = edgesArr[i].PathMat[:1]
				edgesArr[i].PathMat1 = edgesArr[i].PathMat1[:1]
				fmt.Printf("[MergePathArr]mergeed successed, edgesArr[%v]PathMat : %v, edgesArr[%v]PathMat1: %v\n", i, edgesArr[i].PathMat, i, edgesArr[i].PathMat1)
			}
		}
	}*/
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
	buffp := bufio.NewWriterSize(fp, 1<<25)

	fafp := fasta.NewWriter(buffp, 100)
	seqID := int(1)

	for _, p := range joinPathArr {
		if len(p.IDArr) <= 1 {
			continue
		}
		fmt.Printf("[ReconstructDBG] p: %v\n", p)
		//if IsTwoEdgeCyclePath(path) { joinPathArr[i].IDArr = }

		//e := constructdbg.CascadePath(p, edgesArr, nodesArr, kmerlen, false)
		seq := GetContigSeq(p, edgesArr, nodesArr, kmerlen)

		contig := linear.NewSeq("", constructdbg.Transform2Letters(seq), alphabet.DNA)
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
	if err := buffp.Flush(); err != nil {
		log.Fatalf("[StoreEdgesToFn] failed to flush file: %s, err: %v\n", DcDBGEdgesfn, err)
	}
}

func SimplifyByLongReadsPath(edgesArr []constructdbg.DBGEdge, nodesArr []constructdbg.DBGNode, pathArr constructdbg.PathArr, opt Options, minUniqueEdgeLen, minMappingLen, depth, avgReadLen int) (mergePathArr []constructdbg.Path) {
	readPathArr := pathArr.Arr
	sortedEIDIdxArr := SortIdxByUniqueEdgeLen(edgesArr)
	edgesPathRelationArr := make([][]uint32, len(edgesArr))
	edgesPathRelationArr = ConstructEdgesPathRelationship(edgesArr, nodesArr, readPathArr, edgesPathRelationArr, opt.Kmer*3/2, opt.Kmer)

	CleanNGSPath(edgesArr)
	MergePathArr(edgesArr, nodesArr, readPathArr, edgesPathRelationArr, opt.MinMapFreq, opt.Kmer, minUniqueEdgeLen, depth, avgReadLen)
	//MergePathArr(edgesArr, nodesArr, pathArr, edgesPathRelationArr, opt.MinMapFreq)
	//constructdbg.MergePathMat(edgesArr, nodesArr, 2)
	fmt.Printf("[SimplifyByLongReadsPath] sortedEIDIdxArr: %v\n", sortedEIDIdxArr)
	mergePathArr = FindMaxPath(sortedEIDIdxArr, edgesArr, nodesArr, readPathArr, edgesPathRelationArr, opt.MinMapFreq, minMappingLen, opt.Kmer)
	// constuct map edge ID to the path
	//IDMapPath := constructdbg.ConstructIDMapPath(pathArr)
	//constructdbg.DeleteJoinPathArrEnd(edgesArr, pathArr)
	return
	//return mergePathArr
}

func CompressToUint64(ks []byte) (cb uint64) {
	for _, b := range ks {
		cb <<= bnt.NumBitsInBase
		cb |= uint64(b)
	}
	return
}

func GetKmerInfo(ks []byte, ID, pos uint32, strand bool) (ki KmerInfo) {
	ki.Kmer = CompressToUint64(ks)
	ki.ID = ID
	ki.SetPos(pos)
	ki.SetStrand(strand)
	return
}

func FindMinKmerInfo(ka []KmerInfo) (min KmerInfo, idx int) {
	min = ka[0]
	idx = 0
	for i := 1; i < len(ka); i++ {
		ki := ka[i]
		if min.Kmer > ki.Kmer {
			min = ki
			idx = i
		}
	}
	return
}

func GetSeqMiniKmerInfoArr(ks []byte, ID uint32, SeedLen, WinSize int) []KmerInfo {
	var kmerInfoArr []KmerInfo
	rs := constructdbg.GetReverseCompByteArr(ks)
	sl := len(ks)
	bufSize := 100
	buf := make([]KmerInfo, bufSize)
	idx := 0
	last := -1
	var ki, rki KmerInfo
	for i := 0; i < WinSize-1; i++ {
		ki = GetKmerInfo(ks[i:i+SeedLen], ID, uint32(i), constructdbg.PLUS)
		rki = GetKmerInfo(rs[sl-i-SeedLen:sl-i], ID, uint32(i), constructdbg.MINUS)
		if ki.Kmer < rki.Kmer {
			buf[idx] = ki
		} else {
			buf[idx] = rki
		}
		idx++
	}
	//ki := GetKmerInfo(ks[WinSize-1:WinSize-1+SeedLen], ID, uint32(WinSize-1), true)
	//rki := GetKmerInfo(rs[sl-(WinSize-1)-SeedLen:sl-(WinSize-1)], ID, uint32(WinSize-1), false)
	BITNUM := (uint64(SeedLen) - 1) * bnt.NumBitsInBase
	MUSK := (uint64(1) << (uint64(SeedLen) * bnt.NumBitsInBase)) - 1
	for i := WinSize - 1; i < sl-(SeedLen-1); i++ {
		ki.Kmer <<= bnt.NumBitsInBase
		ki.Kmer |= uint64(ks[i+SeedLen-1])
		ki.Kmer &= MUSK //rki.Kmer &= MUSK
		ki.SetPos(uint32(i))
		ki.SetStrand(constructdbg.PLUS)

		rki.Kmer >>= bnt.NumBitsInBase
		var bt uint64
		bt = uint64(rs[sl-i-SeedLen])
		rki.Kmer |= (bt << BITNUM)
		rki.SetPos(uint32(i))
		rki.SetStrand(constructdbg.MINUS)
		if ki.Kmer < rki.Kmer {
			buf[idx] = ki
		} else {
			buf[idx] = rki
		}
		idx++

		min, x := FindMinKmerInfo(buf[idx-WinSize : idx])
		j := idx - WinSize + x
		//fmt.Printf("[GetSeqMiniKmerInfos]last: %v, j: %v, min: %v\n", last, j, min)
		if last < 0 {
			kmerInfoArr = append(kmerInfoArr, min)
			last = j
		} else {
			if j > last {
				kmerInfoArr = append(kmerInfoArr, min)
				last = j
			}
		}

		// adjust buf
		if idx == bufSize {
			y := 0
			for x := idx - WinSize; x < idx; x++ {
				buf[y] = buf[x]
				y++
			}
			last = last - (idx - WinSize)
			idx = y
		}
	}

	return kmerInfoArr
}

type KIArr []KmerInfo

func (arr KIArr) Len() int {
	return len(arr)
}

func (arr KIArr) Less(i, j int) bool {
	return arr[i].Kmer < arr[j].Kmer
}
func (arr KIArr) Swap(i, j int) {
	arr[i], arr[j] = arr[j], arr[i]
}

func SortDBGEdgesKmer(edgesArr []constructdbg.DBGEdge, SeedLen, WinSize int) (sortKmerArr []KmerInfo) {
	var kmerArr []KmerInfo
	for i := 0; i < len(edgesArr); i++ {
		e := edgesArr[i]
		if e.ID < 2 || e.GetDeleteFlag() > 0 {
			continue
		}
		kmerArr = append(kmerArr, GetSeqMiniKmerInfoArr(e.Utg.Ks, uint32(e.ID), SeedLen, WinSize)...)
	}
	sort.Sort(KIArr(kmerArr))

	sortKmerArr = kmerArr
	return
}

func SortLongReadsKmer(rc <-chan []constructcf.ReadInfo, SeedLen, WinSize int, bufSize int64, cs chan<- ReadsBucketInfo, finishT chan<- int) {
	//var addReadKmerLen int64
	for {
		ra := <-rc
		if len(ra) == 0 {
			break
		}
		var kmerArr []KmerInfo
		//kmerArr := make([]KmerInfo, 0, bufSize/int64(WinSize)+(1<<20))
		for _, ri := range ra {
			kmerArr = append(kmerArr, GetSeqMiniKmerInfoArr(ri.Seq, uint32(ri.ID), SeedLen, WinSize)...)
		}
		sort.Sort(KIArr(kmerArr))
		var rbi ReadsBucketInfo
		rbi.ReadsArr = ra
		rbi.KmerSortArr = kmerArr
		cs <- rbi
	}
	finishT <- 1
}

type Options struct {
	utils.ArgsOpt
	//MaxMapEdgeLen int
	MinDepth         int
	AvgDepth         int
	MinUniqueEdgeLen int
	AvgReadLen       int
	WinSize          int
	MaxNGSReadLen    int
	MinMapFreq       int
	ExtLen           int
	ONTFn            string
	Correct          bool
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

	opt.MinUniqueEdgeLen, ok = c.Flag("MinUniqueEdgeLen").Get().(int)
	if !ok {
		log.Fatalf("[checkArgs] argument 'MinUniqueEdgeLen': %v set error\n ", c.Flag("MinUniqueEdgeLen").String())
	}
	if opt.MinUniqueEdgeLen < 1000 || opt.MinUniqueEdgeLen > 1000000 {
		log.Fatalf("[checkArgs] argument 'MinUniqueEdgeLen': %v must between [1k~1m]\n", c.Flag("MinUniqueEdgeLen"))
	}
	opt.AvgReadLen, ok = c.Flag("AvgReadLen").Get().(int)
	if !ok {
		log.Fatalf("[checkArgs] argument 'AvgReadLen': %v set error\n ", c.Flag("AvgReadLen").String())
	}
	if opt.AvgReadLen < 1000 || opt.AvgReadLen > 1000000 {
		log.Fatalf("[checkArgs] argument 'AvgReadLen': %v must between [1k~1m]\n", c.Flag("AvgReadLen"))
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

type MapingKmerInfo struct {
	EID          uint32
	RInfo, EInfo uint32
	Len          uint16 // mapping base length
}

func (m MapingKmerInfo) GetRPos() uint32 {
	return m.RInfo >> 1
}
func (m *MapingKmerInfo) SetRPos(p uint32) {
	m.RInfo = (p << 1) | (m.RInfo & 0x1)
}
func (m MapingKmerInfo) GetEPos() uint32 {
	return m.EInfo >> 1
}
func (m *MapingKmerInfo) SetEPos(p uint32) {
	m.EInfo = (p << 1) | (m.EInfo & 0x1)
}
func (m MapingKmerInfo) GetRStrand() bool {
	if (m.RInfo & 0x1) > 0 {
		return constructdbg.PLUS
	} else {
		return constructdbg.MINUS
	}
}
func (m MapingKmerInfo) GetEStrand() bool {
	if (m.EInfo & 0x1) > 0 {
		return constructdbg.PLUS
	} else {
		return constructdbg.MINUS
	}
}

func (m MapingKmerInfo) GetStrand() bool {
	if m.GetEStrand() == m.GetRStrand() {
		return constructdbg.PLUS
	} else {
		return constructdbg.MINUS
	}
}

type MapingKmerInfoArr []MapingKmerInfo

func (arr MapingKmerInfoArr) Len() int {
	return len(arr)
}

func (arr MapingKmerInfoArr) Less(i, j int) bool {
	return (arr[i].RInfo >> 1) < (arr[j].RInfo >> 1)
}
func (arr MapingKmerInfoArr) Swap(i, j int) {
	arr[i], arr[j] = arr[j], arr[i]
}

func GetLongReadsFile(cfgInfo constructcf.CfgInfo) (fnArr []string) {
	for _, lib := range cfgInfo.Libs {
		if lib.SeqProfile == 3 || lib.SeqProfile == 4 {
			fnArr = append(fnArr, lib.FnName...)
		}
	}
	return
}

type ReadsBucketInfo struct {
	ReadsArr    []constructcf.ReadInfo
	KmerSortArr []KmerInfo
}

func LoadLongReads(ONTFnbr string, rc chan<- []constructcf.ReadInfo, bufSize int64, processT chan<- int) {
	var processNumReads, totalBases int
	format := constructcf.GetReadsFileFormat(ONTFnbr)
	fp, err := os.Open(ONTFnbr)
	if err != nil {
		log.Fatalf("[LoadLongReads] open file: %v failed..., err: %v\n", ONTFnbr, err)
	}
	defer fp.Close()
	fmt.Printf("[LoadLongReads] processing reads in file: %v\n", ONTFnbr)
	brfp := cbrotli.NewReaderSize(fp, 1<<25)
	defer brfp.Close()
	buffp := bufio.NewReader(brfp)
	var err1 error
	var readInfoArr []constructcf.ReadInfo
	var processBases int64
	for err1 != io.EOF {
		ri, err1 := constructcf.GetReadFileRecord(buffp, format, false)
		if ri.ID == 0 {
			if err1 == io.EOF {
				break
			} else {
				log.Fatalf("[LoadLongReads] file: %s encounter err: %v\n", ONTFnbr, err1)
			}
		}
		if ri.ID > math.MaxUint32 && ri.ID <= 0 {
			log.Fatalf("[LoadLongReads] read ID must [1~%v], ri.ID: %v\n", math.MaxUint32, ri.ID)
		}
		readInfoArr = append(readInfoArr, ri)
		processNumReads++
		totalBases += len(ri.Seq)
		processBases += int64(len(ri.Seq))
		if processBases >= bufSize {
			rc <- readInfoArr
			var na []constructcf.ReadInfo
			readInfoArr = na
			processBases = 0
		}
	}
	if len(readInfoArr) > 0 {
		rc <- readInfoArr
	}

	fmt.Printf("[LoadLongReads] processed reads number is: %d, totalBases: %v, finished processed file: %v\n", processNumReads, totalBases, ONTFnbr)
	processT <- 1
}

func LoadLongReadsAndSort(cfgInfo constructcf.CfgInfo, NumCPU int, cs chan<- ReadsBucketInfo, SeedLen, WinSize int, bufSize int64) {
	fnArr := GetLongReadsFile(cfgInfo)
	rc := make(chan []constructcf.ReadInfo, NumCPU)
	processCPUNum := NumCPU * 2 / 3
	LoadCPUNum := NumCPU - processCPUNum
	if processCPUNum <= 2 {
		processCPUNum = 1
		LoadCPUNum = 1
	}
	loadT := make(chan int, LoadCPUNum)
	finishedT := make(chan int, processCPUNum)
	for i := 0; i < processCPUNum; i++ {
		go SortLongReadsKmer(rc, SeedLen, WinSize, bufSize, cs, finishedT)
	}
	for i := 0; i < LoadCPUNum; i++ {
		loadT <- 1
	}

	for i := 0; i < len(fnArr); i++ {
		<-loadT
		go LoadLongReads(fnArr[i], rc, bufSize, loadT)
	}

	for i := 0; i < LoadCPUNum; i++ {
		<-loadT
	}
	close(rc)
	for i := 0; i < processCPUNum; i++ {
		<-finishedT
	}
	close(cs)
}

func GetShareKmer(edgesKmerSortArr, kiArr []KmerInfo, firstReadID uint32, readsNum int) [][]MapingKmerInfo {
	mkiArr := make([][]MapingKmerInfo, readsNum)
	i, j := 0, 0
	for ; i < len(kiArr) && j < len(edgesKmerSortArr); i++ {
		kiR, kiE := kiArr[i], edgesKmerSortArr[j]
		if kiR.Kmer < kiE.Kmer {
			continue
		} else if kiR.Kmer > kiE.Kmer {
			for ; j < len(edgesKmerSortArr); j++ {
				kiE = edgesKmerSortArr[j]
				if kiR.Kmer <= kiE.Kmer {
					break
				}
			}
		}
		if kiR.Kmer < kiE.Kmer {
			continue
		} else if kiR.Kmer > kiE.Kmer {
			break
		}
		// kiR.Kmer == kiE.Kmer
		for m := j; m < len(edgesKmerSortArr); m++ {
			tmp := edgesKmerSortArr[m]
			if tmp.Kmer > kiR.Kmer {
				break
			}
			var mki MapingKmerInfo
			mki.EID, mki.RInfo, mki.EInfo = kiE.ID, kiR.Info, kiE.Info
			mkiArr[kiE.ID-firstReadID] = append(mkiArr[kiE.ID-firstReadID], mki)
		}
	}

	return mkiArr
}

func TransformToChainArr(mkByEdgeMap map[uint32][]MapingKmerInfo, edgeNum, MinSeedNum int, SeedLen uint32) (edgesChainArr [][]Chain, edgeMapInfoArr []EdgeMapInfo) {
	edgesChainArr = make([][]Chain, edgeNum)
	edgeMapInfoArr = make([]EdgeMapInfo, edgeNum)
	idx := 0
	for k, v := range mkByEdgeMap {
		if len(v) >= MinSeedNum {
			edgeMapInfoArr[idx].EdgeID = k
			// choose best mapping strand
			var PLUSNum, MINUSNum int
			for _, mki := range v {
				if mki.GetEStrand() == mki.GetRStrand() {
					PLUSNum++
				} else {
					MINUSNum++
				}
			}
			var strand bool
			if PLUSNum > MINUSNum {
				strand = constructdbg.PLUS
			} else {
				strand = constructdbg.MINUS
			}
			edgeMapInfoArr[idx].Strand = strand
			// choose chain by strand
			for _, mki := range v {
				if strand {
					if mki.GetEStrand() != mki.GetRStrand() {
						continue
					}
				} else {
					if mki.GetEStrand() == mki.GetRStrand() {
						continue
					}
				}
				var c Chain
				c.X, c.Y = mki.GetEPos(), mki.GetRPos()
				c.Len = SeedLen
				edgesChainArr[idx] = append(edgesChainArr[idx], c)
			}
			idx++
		}
	}
	return
}

func GetFlankReadPos(maxChainArr []Chain) (startP, endP uint32) {
	fc, lc := maxChainArr[0], maxChainArr[len(maxChainArr)-1]
	if fc.Y < lc.Y {
		startP = fc.Y
		endP = lc.Y + lc.Len
	} else {
		startP = lc.Y
		endP = fc.Y + fc.Len
	}
	return
}

type EdgeMapInfo struct {
	EdgeID uint32
	Strand bool
}

func GetRegionChainArrScore(ca []Chain, sp, ep uint32) (sc int) {
	if ca[0].Y < ca[1].Y {
		for _, c := range ca {
			if c.Y+c.Len < sp {
				continue
			} else if c.Y+c.Len > ep {
				break
			}
			sc += int(c.Len)
		}
	} else {
		for i := len(ca) - 1; i >= 0; i-- {
			c := ca[i]
			if c.Y+c.Len < sp {
				continue
			} else if c.Y+c.Len > ep {
				break
			}
			sc += int(c.Len)
		}
	}
	return
}

func FindReadRegionHighQualityEdges(karr []MapingKmerInfo, startP, endP, edgeID uint32, maxChainArr []Chain, SeedLen, MinSeedNum int) (highMapingQualityEdgesArr []uint32) {
	// classify by edge
	mkByEdgeMap, edgeNum := ClassifyByEdge(karr, startP, endP, uint32(SeedLen))
	delete(mkByEdgeMap, edgeID)
	// transform to Chain array
	edgesChainArr, edgeMapInfoArr := TransformToChainArr(mkByEdgeMap, edgeNum, MinSeedNum, uint32(SeedLen))

	// find max score by every edge
	edgesMaxChainArr := make([][]Chain, edgeNum)
	edgesMaxScore := make([]int, edgeNum)
	for i, carr := range edgesChainArr {
		edgesMaxChainArr[i], edgesMaxScore[i] = GetMaxChainArr(carr, edgeMapInfoArr[i].Strand, uint32(SeedLen))
	}
	for i, mc := range edgesMaxChainArr {
		sp, ep := GetFlankReadPos(mc)
		sc := GetRegionChainArrScore(maxChainArr, sp, ep)
		if sc <= edgesMaxScore[i] {
			highMapingQualityEdgesArr = append(highMapingQualityEdgesArr, edgeMapInfoArr[i].EdgeID)
		}
	}
	return
}

func ClassifyByEdge(karr []MapingKmerInfo, startP, endP, SeedLen uint32) (mkByEdgeMap map[uint32][]MapingKmerInfo, edgeNum int) {
	mkByEdgeMap = make(map[uint32][]MapingKmerInfo)
	for _, mki := range karr {
		if mki.GetEPos() < startP {
			continue
		} else if mki.GetEPos()+SeedLen > endP {
			break
		}
		if _, ok := mkByEdgeMap[mki.EID]; ok {
			mkByEdgeMap[mki.EID] = append(mkByEdgeMap[mki.EID], mki)
		} else {
			mkByEdgeMap[mki.EID] = []MapingKmerInfo{mki}
			edgeNum++
		}
	}
	return
}

type UniqueRegion struct {
	StartP, EndP int
	EdgeID       uint32
	MaxScore     int
	ChainArr     []Chain
}
type UniqueRegionArr []UniqueRegion

func (arr UniqueRegionArr) Len() int {
	return len(arr)
}

func (arr UniqueRegionArr) Less(i, j int) bool {
	return arr[i].MaxScore < arr[j].MaxScore
}
func (arr UniqueRegionArr) Swap(i, j int) {
	arr[i], arr[j] = arr[j], arr[i]
}

func FindPrimaryMappingEdges(karr []MapingKmerInfo, ReadLen, SeedLen int) {
	MinSeedNum := 4 // minimum seed mapping number for one hit
	MinUniqueRegionLen := 500
	// classify by edge
	mkByEdgeMap, edgeNum := ClassifyByEdge(karr, 0, uint32(ReadLen), uint32(SeedLen))

	// transform to Chain array
	edgesChainArr, edgeMapInfoArr := TransformToChainArr(mkByEdgeMap, edgeNum, MinSeedNum, uint32(SeedLen))

	// find max score by every edge
	edgesMaxChainArr := make([][]Chain, edgeNum)
	edgesMaxScore := make([]int, edgeNum)
	for i, carr := range edgesChainArr {
		edgesMaxChainArr[i], edgesMaxScore[i] = GetMaxChainArr(carr, edgeMapInfoArr[i].Strand, uint32(SeedLen))
	}

	// confirm long read unique mapping region
	scArr := make([]Score, edgeNum)
	for i, ms := range edgesMaxScore {
		scArr[i].Sc = ms
		scArr[i].Idx = i
	}
	sort.Sort(ScoreArr(scArr))
	var uniqArr []UniqueRegion
	for i := edgeNum - 1; i >= 0; i-- {
		sc := scArr[i]
		maxChainArr := edgesMaxChainArr[sc.Idx]
		startP, endP := GetFlankReadPos(maxChainArr)
		if int(endP)-int(startP) < MinUniqueRegionLen {
			break
		}
		highMapingQualityEdgesArr := FindReadRegionHighQualityEdges(karr, startP, endP, edgeMapInfoArr[sc.Idx].EdgeID, maxChainArr, SeedLen, MinSeedNum)
		if len(highMapingQualityEdgesArr) > 0 {
			continue
		}
		// add to unique region
		var ur UniqueRegion
		ur.StartP, ur.EndP = int(startP), int(endP)
		ur.EdgeID = edgeMapInfoArr[sc.Idx].EdgeID
		ur.MaxScore = edgesMaxScore[sc.Idx]
		ur.ChainArr = maxChainArr
		uniqArr = append(uniqArr, ur)
	}

	// Extend path by unique region
	sort.Sort(UniqueRegionArr(uniqArr))
	//extendPathArr := MappingMostProbablePath()
}

/*func FindONTReadsPath(edgesKmerSortArr []KmerInfo, cs <-chan ReadsBucketInfo, SeedLen, WinSize int, wc chan<- [2][]constructdbg.DBG_MAX_INT) {
	for {
		rbi := <-cs
		if len(rbi.ReadsArr) == 0 {
			var tmp [2][]constructdbg.DBG_MAX_INT
			wc <- tmp
			break
		}
		rArr := rbi.ReadsArr
		kIArr := rbi.KmerSortArr
		// first find ONT reads mapping seed chains
		// read ID must successive uint32
		kmerShareArrByReadID := GetShareKmer(edgesKmerSortArr, kIArr, uint32(rArr[0].ID), len(rArr))
		for i, karr := range kmerShareArrByReadID {
			// sort by long read Postion
			sort.Sort(MapingKmerInfoArr(karr))
			// find primary mapping Edges and clean secondary edges
			edgesMapingInfo := FindPrimaryMappingEdges(karr, len(rArr[i].Seq), SeedLen)
			// find seed edges and local align uncertain edges

			// output edge mapping path

		}

	}
}*/

type LongReadMappingInfo struct {
	ID           uint32 // Read ID
	Qual         uint32 // mapping Quality [0~100]
	RStart, REnd uint32 // start and end position of long Read sequence
	EStart, EEnd uint32 // start and end position of edges
	AnchorEdgeID uint32 // anchor(alignment start edge) edge ID
	AnchorStrand bool   //
	Score        int    // total seed chain maping score
	Path         [2][]constructdbg.DBG_MAX_INT
}

func GetInterSectionKmerInfo(kmerInfoArrA, kmerInfoArrB []KmerInfo, BLen int, StartID uint32, SeedLen uint16) (bucketInterSecKmerInfoArr [][]MapingKmerInfo) {
	//startReadID := uint32(rb.ReadsArr[0].ID)
	bucketInterSecKmerInfoArr = make([][]MapingKmerInfo, BLen)
	i, j := 0, 0
	ek := kmerInfoArrA[j]
	for ; i < len(kmerInfoArrB) && j < len(kmerInfoArrA); i++ {
		//fmt.Printf("[GetInterSectionKmerInfo]kmerInfoArrB[%v]: %v, kmerInfoArrA[%v]: %v\n", i, kmerInfoArrB[i], j, kmerInfoArrA[j])
		rk := kmerInfoArrB[i]
		if rk.Kmer < ek.Kmer {
			continue
		} else if rk.Kmer > ek.Kmer {
			for j = j + 1; j < len(kmerInfoArrA); j++ {
				ek = kmerInfoArrA[j]
				if rk.Kmer <= ek.Kmer {
					break
				}
			}
		}

		if rk.Kmer < ek.Kmer {
			continue
		} else if rk.Kmer > ek.Kmer {
			break
		}
		// rk.Kmer == ek.Kmer
		for m := j; m < len(kmerInfoArrA); m++ {
			tk := kmerInfoArrA[m]
			if tk.Kmer > rk.Kmer {
				break
			}
			var mk MapingKmerInfo
			mk.EID, mk.EInfo, mk.RInfo = tk.ID, tk.Info, rk.Info
			mk.Len = SeedLen
			//fmt.Printf("[GetInterSectionKmerInfo]rk: %v, tk: %v\n", rk, tk)

			//fmt.Printf("[GetInterSectionKmerInfo]edgeID: %v: ReadPos: %v, EdgePos: %v,Kmer: %v\n", rk.ID, tk.GetPos(), rk.GetPos(), tk.Kmer)
			bucketInterSecKmerInfoArr[rk.ID-StartID] = append(bucketInterSecKmerInfoArr[rk.ID-StartID], mk)
		}
	}

	return
}

func GetChainBlocks(ka []MapingKmerInfo) (kb []MapingKmerInfo) {
	flagArr := make([]bool, len(ka))
	loopNum := 0
	for { // until to the no flag false
		loopNum++
		kbl := len(kb)
		for i := 0; i < len(flagArr); i++ {
			if flagArr[i] == true {
				//i++
				continue
			}
			mk := ka[i]
			flagArr[i] = true
			j := i + 1
			for ; j < len(flagArr); j++ {
				nk := ka[j]
				if flagArr[j] == true {
					continue
				}
				if nk.GetRPos() > mk.GetRPos()+uint32(mk.Len) {
					break
				}
				if mk.EID != nk.EID || mk.GetStrand() != nk.GetStrand() {
					continue
				}
				//fmt.Printf("[GetChainBlocks] EID: %v, mk RPos: %v, EPos: %v,Len: %v, nk RPos: %v, EPos: %v, Len:%v\n", mk.EID, mk.GetRPos(), mk.GetEPos(), mk.Len, nk.GetRPos(), nk.GetEPos(), nk.Len)
				strand := mk.GetStrand()
				if strand == constructdbg.PLUS {
					if nk.GetRPos()-mk.GetRPos() == nk.GetEPos()-mk.GetEPos() {
						mk.Len = uint16(nk.GetRPos() + uint32(nk.Len) - mk.GetRPos())
						//mk.Len += uint16(((nk.GetRPos() + uint32(nk.Len)) - (mk.GetRPos() + uint32(mk.Len))))
						flagArr[j] = true
					}
				} else {
					if nk.GetRPos()-mk.GetRPos() == (mk.GetEPos()+uint32(mk.Len))-(nk.GetEPos()+uint32(nk.Len)) {
						//fmt.Printf("[GetChainBlocks] EID: %v, mk RPos: %v, EPos: %v,Len: %v, nk RPos: %v, EPos: %v, Len:%v\n", mk.EID, mk.GetRPos(), mk.GetEPos(), mk.Len, nk.GetRPos(), nk.GetEPos(), nk.Len)
						mk.Len = uint16(nk.GetRPos() + uint32(nk.Len) - mk.GetRPos())
						mk.SetEPos(nk.GetEPos())
						//mk.Len += uint16(((nk.GetRPos() + uint32(nk.Len)) - (mk.GetRPos() + uint32(mk.Len))))
						flagArr[j] = true
						//fmt.Printf("[GetChainBlocks] mk RPos: %v, EPos: %v,Len: %v\n", mk.GetRPos(), mk.GetEPos(), mk.Len)
					}
				}
			}

			kb = append(kb, mk)
			//i = j
		}
		// no new added block
		if kbl == len(kb) {
			break
		}
	}
	fmt.Printf("[GetChainBlocks] loop num: %v\n", loopNum)

	return
}

type EdgeKmerInfo struct {
	PlusNum, MinusNum int
	ScoreP, ScoreM    int
	IdxP, IdxM        int
	MkInfoArr         []MapingKmerInfo
}

type ChainBlockMapLongReadInfo struct {
	//EID        constructdbg.DBG_MAX_INT
	Start, End int
	Score      int
}

type EdgeKI struct {
	Strand    bool
	Score     int
	Idx       int
	MkInfoArr []MapingKmerInfo
}

type EdgeKIArr []EdgeKI

func (arr EdgeKIArr) Len() int {
	return len(arr)
}

func (arr EdgeKIArr) Less(i, j int) bool {
	return arr[i].Score < arr[j].Score
}

func (arr EdgeKIArr) Swap(i, j int) {
	arr[i], arr[j] = arr[j], arr[i]
}

func GetEdgesKIArr(kb []MapingKmerInfo, edgesArr []constructdbg.DBGEdge, MinKmerScore int) (edgesKIArr []EdgeKI) {
	MinScore := MinKmerScore
	edgeKIMap := make(map[uint32]EdgeKmerInfo)
	for i, mk := range kb {
		var eki EdgeKmerInfo
		eki, _ = edgeKIMap[mk.EID]
		eki.MkInfoArr = append(eki.MkInfoArr, mk)
		if mk.GetStrand() == constructdbg.PLUS {
			eki.PlusNum++
			eki.ScoreP += int(mk.Len)
			if mk.Len > kb[eki.IdxP].Len {
				eki.IdxP = i
			}
		} else {
			eki.MinusNum++
			eki.ScoreM += int(mk.Len)
			if mk.Len > kb[eki.IdxM].Len {
				eki.IdxM = i
			}
		}
		edgeKIMap[mk.EID] = eki
	}

	edgesKIArr = make([]EdgeKI, 0, len(edgeKIMap))
	for _, v := range edgeKIMap {
		if v.ScoreM < MinScore && v.ScoreP < MinScore {
			continue
		}
		edgeLen := len(edgesArr[kb[v.IdxM].EID].Utg.Ks)
		if edgeLen < 500 {
			if v.ScoreM < edgeLen/15 && v.ScoreP < edgeLen/15 {
				continue
			}
		}

		var eKI EdgeKI
		if v.ScoreP > v.ScoreM {
			eKI.Strand = constructdbg.PLUS
			eKI.Score = v.ScoreP
			eKI.Idx = v.IdxP
		} else {
			eKI.Strand = constructdbg.MINUS
			eKI.Score = v.ScoreM
			eKI.Idx = v.IdxM
		}
		eKI.MkInfoArr = v.MkInfoArr
		edgesKIArr = append(edgesKIArr, eKI)
	}
	return
}

func GetTwoBlockScore(cen, mk MapingKmerInfo, GapCost int, edgesArr []constructdbg.DBGEdge, flankLen, kmerlen int) (sc int) {
	var gapLen, dX, dY int
	if mk.GetStrand() == constructdbg.PLUS {
		if cen.GetRPos() < mk.GetRPos() {
			dX = int(mk.GetRPos()) - (int(cen.GetRPos()) + int(cen.Len))
			if cen.EID != mk.EID {
				dY = int(mk.GetEPos()) - (kmerlen - 1) + flankLen
			} else {
				dY = int(mk.GetEPos()) - (int(cen.GetEPos()) + int(cen.Len))
			}
		} else {
			dX = int(cen.GetRPos()) - (int(mk.GetRPos()) + int(mk.Len))
			if cen.EID != mk.EID {
				dY = len(edgesArr[mk.EID].Utg.Ks) - (kmerlen - 1) + flankLen - (int(mk.GetEPos()) + int(mk.Len))
			} else {
				dY = int(cen.GetEPos()) - (int(mk.GetEPos()) + int(mk.Len))
			}
		}
	} else {
		if cen.GetRPos() < mk.GetRPos() {
			dX = int(mk.GetRPos()) - (int(cen.GetRPos()) + int(cen.Len))
			if cen.EID != mk.EID {
				dY = len(edgesArr[mk.EID].Utg.Ks) - (kmerlen - 1) + flankLen - (int(mk.GetEPos()) + int(mk.Len))
			} else {
				dY = int(cen.GetEPos()) - (int(mk.GetEPos()) + int(mk.Len))
			}
		} else {
			dX = int(cen.GetRPos()) - (int(mk.GetRPos()) + int(mk.Len))
			if cen.EID != mk.EID {
				dY = int(mk.GetEPos()) - ((kmerlen - 1) - flankLen)
			} else {
				dY = int(mk.GetEPos()) - (int(cen.GetEPos()) + int(cen.Len))
			}
		}
	}
	/*if dX > 2000 || dY > 2000 {
		sc = 0
		return
	}*/
	gapLen = utils.AbsInt(dX - dY)
	min := constructdbg.Min(dX, dY)
	if gapLen > min/5+20 {
		sc = 0
		return
	}
	a := int(mk.Len)
	// a gap length of gapLen costs
	var b float64
	if gapLen > 1 {
		b = float64(0.01)*float64(mk.Len)*float64(gapLen) + float64(0.5)*math.Log2(float64(gapLen))
	} else {
		if dX == 0 && dY == 0 {
			b = 0
		} else {
			b = 1
		}
	}
	sc = a - int(b)
	//sc = int(mk.Len) - gapLen*GapCost
	//fmt.Printf("[GetTwoBlockScore]dX: %v, dY: %v, sc: %v, a: %v, b: %v\n", dX, dY, sc, a, b)
	return
}

func CutMapingKmerBoundary(cen, mk MapingKmerInfo, direction uint8, strand bool) MapingKmerInfo {
	if cen.EID != mk.EID {
		log.Fatalf("[CutMapingKmerBoundary] cen.EID: %v != mk.EID: %v\n", cen.EID, mk.EID)
	}
	if direction == constructdbg.FORWARD {
		dR := int(cen.GetRPos()) + int(cen.Len) - int(mk.GetRPos())
		if strand == constructdbg.PLUS {
			dE := int(cen.GetEPos()) + int(cen.Len) - int(mk.GetEPos())
			max := utils.MaxInt(dR, dE)
			if max > 0 {
				if max > int(mk.Len) {
					max = int(mk.Len)
				}
				mk.SetRPos(mk.GetRPos() + uint32(max))
				mk.SetEPos(mk.GetEPos() + uint32(max))
				mk.Len -= uint16(max)
			}
		} else {
			dE := int(mk.GetEPos()) + int(mk.Len) - int(cen.GetEPos())
			max := utils.MaxInt(dR, dE)
			if max > 0 {
				if max > int(mk.Len) {
					max = int(mk.Len)
				}
				mk.SetRPos(mk.GetRPos() + uint32(max))
				//mk.SetEPos(mk.GetEPos() + max)
				mk.Len -= uint16(max)
			}
		}
	} else {
		dR := int(mk.GetRPos()) + int(mk.Len) - int(cen.GetRPos())
		if strand == constructdbg.PLUS {
			dE := int(mk.GetEPos()) + int(mk.Len) - int(cen.GetEPos())
			max := utils.MaxInt(dR, dE)
			if max > 0 {
				if max > int(mk.Len) {
					max = int(mk.Len)
				}
				mk.Len -= uint16(max)
			}
		} else {
			dE := int(cen.GetEPos()) + int(cen.Len) - int(mk.GetEPos())
			max := utils.MaxInt(dR, dE)
			if max > 0 {
				if max > int(mk.Len) {
					max = int(mk.Len)
				}
				mk.SetEPos(mk.GetEPos() + uint32(max))
				mk.Len -= uint16(max)
			}
		}
	}
	return mk
}

func BlockToNearBlocksScoreArr(kb []MapingKmerInfo, Idx int, strand bool, visitArr []bool, GapCost, MaxGapLen int, SeqType uint8, bsA []int, edgesArr []constructdbg.DBGEdge, flankLen, kmerlen int, SeedLen uint16) {
	cen := kb[Idx]
	// BACKWARD
	for i := Idx - 1; i >= 0; i-- {
		mk := kb[i]
		if visitArr[i] == true {
			break
		}
		if bsA[i] == math.MinInt32 {
			continue
		}
		if mk.GetStrand() != strand {
			continue
		}
		/*if int(cen.GetRPos())-int(mk.GetRPos()) > MaxGapLen {
			break
		}*/
		if strand == constructdbg.PLUS {
			if mk.GetRPos()+uint32(mk.Len) > cen.GetRPos() || mk.GetEPos()+uint32(mk.Len) > cen.GetEPos() {
				mk = CutMapingKmerBoundary(cen, mk, constructdbg.BACKWARD, strand)
				if mk.Len >= SeedLen {
					kb[i] = mk
				} else {
					bsA[i] = math.MinInt32
					continue
				}
			}
			sc := GetTwoBlockScore(cen, mk, GapCost, edgesArr, flankLen, kmerlen)
			if sc > bsA[i] {
				bsA[i] = sc
			}
		} else { // strand == constructdbg.MINUS
			if mk.GetRPos()+uint32(mk.Len) > cen.GetRPos() || mk.GetEPos() < cen.GetEPos()+uint32(cen.Len) {
				mk = CutMapingKmerBoundary(cen, mk, constructdbg.BACKWARD, strand)
				if mk.Len >= SeedLen {
					kb[i] = mk
				} else {
					bsA[i] = math.MinInt32
					continue
				}
			}
			sc := GetTwoBlockScore(cen, mk, GapCost, edgesArr, flankLen, kmerlen)
			if sc > bsA[i] {
				bsA[i] = sc
			}
		}
	}

	// FORWARD
	for i := Idx + 1; i < len(kb); i++ {
		mk := kb[i]
		//fmt.Printf("[BlockToNearBlocksScoreArr] i: %v, mk ReadPos: %v, cen ReadPos: %v\n", i, mk.GetRPos(), cen.GetRPos())
		if visitArr[i] == true {
			break
		}
		if bsA[i] == math.MinInt32 {
			continue
		}
		if mk.GetStrand() != strand {
			continue
		}
		/*if int(mk.GetRPos())-(int(cen.GetRPos())+int(mk.Len)) > MaxGapLen {
			break
		}*/
		var sc int
		if strand == constructdbg.PLUS {
			if cen.GetRPos()+uint32(cen.Len) > mk.GetRPos() || cen.GetEPos()+uint32(cen.Len) > mk.GetEPos() {
				mk = CutMapingKmerBoundary(cen, mk, constructdbg.FORWARD, strand)
				if mk.Len >= SeedLen {
					kb[i] = mk
				} else {
					bsA[i] = math.MinInt32
					continue
				}
			}
			sc = GetTwoBlockScore(cen, mk, GapCost, edgesArr, flankLen, kmerlen)

			if sc > bsA[i] {
				bsA[i] = sc
			}
		} else { //strand == constructdbg.MINUS
			if cen.GetRPos()+uint32(cen.Len) > mk.GetRPos() || mk.GetEPos()+uint32(mk.Len) > cen.GetEPos() {
				mk = CutMapingKmerBoundary(cen, mk, constructdbg.FORWARD, strand)
				if mk.Len >= SeedLen {
					kb[i] = mk
				} else {
					bsA[i] = math.MinInt32
					continue
				}
			}
			sc = GetTwoBlockScore(cen, mk, GapCost, edgesArr, flankLen, kmerlen)

			if sc > bsA[i] {
				bsA[i] = sc
			}
		}
		//fmt.Printf("[BlockToNearBlocksScoreArr] sc: %v\n", sc)
	}

	//sort.Sort(ScoreArr(bsA))
	return
}

func BlockToNearBlocksScoreArrCrossEdges(kb []MapingKmerInfo, Idx int, strand bool, visitArr []bool, GapCost, MaxGapLen int, SeqType uint8, bsA []int, direction uint8, flankLen, kmerlen int, edgesArr []constructdbg.DBGEdge, SeedLen uint16) {
	cen := kb[Idx]
	// BACKWARD
	for i := Idx - 1; i >= 0; i-- {
		mk := kb[i]
		if visitArr[i] == true {
			break
		}
		if bsA[i] == math.MinInt32 {
			continue
		}
		if mk.GetStrand() != strand {
			continue
		}
		/*if int(cen.GetRPos())-int(mk.GetRPos()) > MaxGapLen {
			break
		}*/
		if strand == constructdbg.PLUS {
			if mk.EID != cen.EID {
				if mk.GetRPos()+uint32(mk.Len) > cen.GetRPos() || int(mk.GetEPos())+int(mk.Len) > len(edgesArr[mk.EID].Utg.Ks)-((kmerlen-1)-flankLen) {
					bsA[i] = math.MinInt32
					continue
				}
			} else {
				if mk.GetRPos()+uint32(mk.Len) > cen.GetRPos() || mk.GetEPos()+uint32(mk.Len) > cen.GetEPos() {
					mk = CutMapingKmerBoundary(cen, mk, constructdbg.BACKWARD, strand)
					if mk.Len >= SeedLen {
						kb[Idx] = mk
					} else {
						bsA[i] = math.MinInt32
						continue
					}
				}
			}

			sc := GetTwoBlockScore(cen, mk, GapCost, edgesArr, flankLen, kmerlen)
			if sc > bsA[i] {
				bsA[i] = sc
			}
		} else { // strand == constructdbg.MINUS
			if mk.EID != cen.EID {
				if mk.GetRPos()+uint32(mk.Len) > cen.GetRPos() || int(mk.GetEPos()) < (kmerlen-1)-flankLen {
					bsA[i] = math.MinInt32
					continue
				}
			} else {
				if mk.GetRPos()+uint32(mk.Len) > cen.GetRPos() || mk.GetEPos() < cen.GetEPos()+uint32(cen.Len) {
					mk = CutMapingKmerBoundary(cen, mk, constructdbg.BACKWARD, strand)
					if mk.Len >= SeedLen {
						kb[Idx] = mk
					} else {
						bsA[i] = math.MinInt32
						continue
					}
				}
			}
			sc := GetTwoBlockScore(cen, mk, GapCost, edgesArr, flankLen, kmerlen)
			if sc > bsA[i] {
				bsA[i] = sc
			}
		}
	}

	// FORWARD
	for i := Idx + 1; i < len(kb); i++ {
		mk := kb[i]
		//fmt.Printf("[BlockToNearBlocksScoreArr] i: %v, mk ReadPos: %v, cen ReadPos: %v\n", i, mk.GetRPos(), cen.GetRPos())
		if visitArr[i] == true {
			break
		}
		if bsA[i] == math.MinInt32 {
			continue
		}
		if mk.GetStrand() != strand {
			continue
		}
		/*if int(mk.GetRPos())-(int(cen.GetRPos())+int(mk.Len)) > MaxGapLen {
			break
		}*/
		var sc int
		if strand == constructdbg.PLUS {
			if cen.EID != mk.EID {
				if cen.GetRPos()+uint32(cen.Len) > mk.GetRPos() || (kmerlen-1)-flankLen > int(mk.GetEPos()) {
					bsA[i] = math.MinInt32
					continue
				}
			} else {
				if cen.GetRPos()+uint32(cen.Len) > mk.GetRPos() || cen.GetEPos()+uint32(cen.Len) > mk.GetEPos() {
					mk = CutMapingKmerBoundary(cen, mk, constructdbg.FORWARD, strand)
					if mk.Len >= SeedLen {
						kb[Idx] = mk
					} else {
						bsA[i] = math.MinInt32
						continue
					}
				}
			}
			sc = GetTwoBlockScore(cen, mk, GapCost, edgesArr, flankLen, kmerlen)

			if sc > bsA[i] {
				bsA[i] = sc
			}
		} else { //strand == constructdbg.MINUS
			if cen.EID != mk.EID {
				if cen.GetRPos()+uint32(cen.Len) > mk.GetRPos() || int(mk.GetEPos())+int(mk.Len) > len(edgesArr[mk.EID].Utg.Ks)-((kmerlen-1)-flankLen) {
					bsA[i] = math.MinInt32
					continue
				}
			} else {
				if cen.GetRPos()+uint32(cen.Len) > mk.GetRPos() || mk.GetEPos()+uint32(mk.Len) > cen.GetEPos() {
					mk = CutMapingKmerBoundary(cen, mk, constructdbg.FORWARD, strand)
					if mk.Len >= SeedLen {
						kb[Idx] = mk
					} else {
						bsA[i] = math.MinInt32
						continue
					}
				}
			}
			sc = GetTwoBlockScore(cen, mk, GapCost, edgesArr, flankLen, kmerlen)

			if sc > bsA[i] {
				bsA[i] = sc
			}
		}
		//fmt.Printf("[BlockToNearBlocksScoreArr] sc: %v\n", sc)
	}

	//sort.Sort(ScoreArr(bsA))
	return
}

func GetMaxScoreFromBaSA(blockToAnchorsScoreArr []int, visitArr []bool) (max Score) {
	for i, sc := range blockToAnchorsScoreArr {
		if visitArr[i] == true {
			continue
		}
		if max.Sc < sc {
			max.Sc = sc
			max.Idx = i
		}
	}
	return
}

func SearchNearMK(mkArr []MapingKmerInfo, mk MapingKmerInfo, direction uint8) (nk MapingKmerInfo) {
	if direction == constructdbg.BACKWARD {
		for _, tk := range mkArr {
			if tk.GetRPos() < mk.GetRPos() {
				if tk.GetRPos() > nk.GetRPos() {
					nk = tk
				}
			}
		}
	} else {
		nk.SetRPos(math.MaxUint32 >> 2)
		for _, tk := range mkArr {
			if tk.GetRPos() > mk.GetRPos() {
				if tk.GetRPos() < nk.GetRPos() {
					nk = tk
				}
			}
		}
	}
	return
}

func CheckChainBoundary(mkArr []MapingKmerInfo, mk MapingKmerInfo, maxSc Score) (MapingKmerInfo, Score) {
	sc := maxSc
	lmk := SearchNearMK(mkArr, mk, constructdbg.BACKWARD)
	rmk := SearchNearMK(mkArr, mk, constructdbg.FORWARD)
	if lmk.Len > 0 {
		if lmk.GetRPos()+uint32(lmk.Len) > mk.GetRPos() {
			if mk.GetRPos()+uint32(mk.Len) <= lmk.GetRPos()+uint32(lmk.Len) {
				sc.Sc = -1
				return mk, sc
			}
			mk.Len -= uint16(lmk.GetRPos() + uint32(lmk.Len) - mk.GetRPos())
			sc.Sc -= int(lmk.GetRPos() + uint32(lmk.Len) - mk.GetRPos())
			if mk.GetStrand() == constructdbg.PLUS {
				mk.SetEPos(mk.GetEPos() + (lmk.GetRPos() + uint32(lmk.Len) - mk.GetRPos()))
			} else {
				//
			}
			mk.SetRPos(lmk.GetRPos() + uint32(lmk.Len))
		}
		if mk.GetStrand() == constructdbg.PLUS {
			if mk.GetEPos() < lmk.GetEPos() || mk.GetEPos()+uint32(mk.Len) <= lmk.GetEPos()+uint32(lmk.Len) {
				sc.Sc = -1
				return mk, sc
			} else if mk.GetEPos() < lmk.GetEPos()+uint32(lmk.Len) {
				mk.Len -= uint16(lmk.GetEPos() + uint32(lmk.Len) - mk.GetEPos())
				sc.Sc -= int(lmk.GetEPos() + uint32(lmk.Len) - mk.GetEPos())
				mk.SetRPos(mk.GetRPos() + (lmk.GetEPos() + uint32(lmk.Len) - mk.GetEPos()))
				mk.SetEPos(lmk.GetEPos() + uint32(lmk.Len))
			}
		} else {
			if mk.GetEPos() >= lmk.GetEPos() {
				sc.Sc = -1
				return mk, sc
			} else if mk.GetEPos()+uint32(mk.Len) > lmk.GetEPos() {
				mk.Len -= uint16(mk.GetEPos() + uint32(mk.Len) - lmk.GetRPos())
				sc.Sc -= int(mk.GetEPos() + uint32(mk.Len) - lmk.GetRPos())
				mk.SetRPos(mk.GetRPos() + (mk.GetEPos() + uint32(mk.Len) - lmk.GetRPos()))
				//mk.SetEPos()
			}
		}
	}

	if rmk.Len > 0 {
		if mk.GetRPos()+uint32(mk.Len) > rmk.GetRPos() {
			if mk.GetRPos() >= rmk.GetRPos() {
				sc.Sc = -1
				return mk, sc
			}
			mk.Len -= uint16(mk.GetRPos() + uint32(mk.Len) - rmk.GetRPos())
			sc.Sc -= int(mk.GetRPos() + uint32(mk.Len) - rmk.GetRPos())
			if mk.GetStrand() == constructdbg.PLUS {

			} else {
				mk.SetEPos(mk.GetEPos() + (mk.GetRPos() + uint32(mk.Len) - rmk.GetRPos()))
			}
		}
		if mk.GetStrand() == constructdbg.PLUS {
			if mk.GetEPos() >= rmk.GetEPos() || mk.GetEPos()+uint32(mk.Len) >= rmk.GetEPos()+uint32(rmk.Len) {
				sc.Sc = -1
				return mk, sc
			} else if mk.GetEPos()+uint32(mk.Len) > rmk.GetEPos() {
				mk.Len -= uint16(mk.GetEPos() + uint32(mk.Len) - rmk.GetEPos())
				sc.Sc -= int(mk.GetEPos() + uint32(mk.Len) - rmk.GetEPos())
			}
		} else {
			if mk.GetEPos() < rmk.GetEPos() || mk.GetEPos()+uint32(mk.Len) < rmk.GetEPos()+uint32(rmk.Len) {
				sc.Sc = -1
				return mk, sc
			} else if mk.GetEPos() < rmk.GetEPos()+uint32(rmk.Len) {
				mk.Len -= uint16(rmk.GetEPos() + uint32(rmk.Len) - mk.GetEPos())
				sc.Sc -= int(rmk.GetEPos() + uint32(rmk.Len) - mk.GetEPos())
				mk.SetEPos(mk.GetEPos() + (rmk.GetEPos() + uint32(rmk.Len) - mk.GetEPos()))
			}
		}
	}

	return mk, sc
}

func GetMappingKIArrScore(mkArr []MapingKmerInfo, GapCost int, edgesArr []constructdbg.DBGEdge, flankLen, kmerlen int) (score int) {
	score = int(mkArr[0].Len)
	for i := 1; i < len(mkArr); i++ {
		lk := mkArr[i-1]
		nk := mkArr[i]
		sc := GetTwoBlockScore(lk, nk, GapCost, edgesArr, flankLen, kmerlen)
		//fmt.Printf("[GetMappingKIArrScore] nk EdgePos: %v, ReadPos: %v, Len: %v, score: %v\n", nk.GetEPos(), nk.GetRPos(), nk.Len, sc)
		score += sc
	}

	return
}

type MaxPathInfo struct {
	Arr   []MapingKmerInfo
	Path  []uint32
	Score int
	//DistE int
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

type MaxPathInfoTmpArr []MaxPathInfo

func (arr MaxPathInfoTmpArr) Len() int {
	return len(arr)
}

func (arr MaxPathInfoTmpArr) Less(i, j int) bool {
	if arr[i].Arr[0].GetRPos() < arr[j].Arr[0].GetRPos() {
		return true
	} else if arr[i].Arr[0].GetRPos() > arr[j].Arr[0].GetRPos() {
		return false
	} else {
		return arr[i].Score < arr[j].Score
	}
}

func (arr MaxPathInfoTmpArr) Swap(i, j int) {
	arr[i], arr[j] = arr[j], arr[i]
}

type BoundScore struct {
	RPos  uint32
	Score int
}

func GetAnchorMKI(ki EdgeKI) (anchorMKI MapingKmerInfo, idx int) {
	for i, mki := range ki.MkInfoArr {
		if ki.Strand == mki.GetStrand() {
			if anchorMKI.Len < mki.Len {
				anchorMKI = mki
				idx = i
			}
		}
	}
	return
}

func GetMaxScoreArr(ki EdgeKI, SeedLen, GapCost, MaxGapLen int, SeqType uint8, edgesArr []constructdbg.DBGEdge, kmerlen int) (mp MaxPathInfo) {
	var maxA []MapingKmerInfo
	//var blockToBlocksScorePat [][]Score
	sort.Sort(MapingKmerInfoArr(ki.MkInfoArr))
	/*for i, mki := range ki.MkInfoArr {
		fmt.Printf("[GetMaxScoreArr] ki.MkInfoArr[%v] EdgePos: %v, ReadPos : %v, Len: %v\n", i, mki.GetEPos(), mki.GetRPos(), mki.Len)
	}*/
	anchorMKI, idx := GetAnchorMKI(ki)
	visitArr := make([]bool, len(ki.MkInfoArr))
	fmt.Printf("[GetMaxScoreArr] anchor[%v] edgeID: %v, mki ReadPos : %v,EdgePos: %v, Len: %v\n", idx, anchorMKI.EID, anchorMKI.GetRPos(), anchorMKI.GetEPos(), anchorMKI.Len)
	visitArr[idx] = true
	maxA = append(maxA, anchorMKI)
	bsA := make([]int, len(ki.MkInfoArr))
	bsA[idx] = int(anchorMKI.Len)
	BlockToNearBlocksScoreArr(ki.MkInfoArr, idx, ki.Strand, visitArr, GapCost, MaxGapLen, SeqType, bsA, edgesArr, 0, kmerlen, uint16(SeedLen))
	//blockToBlocksScorePat = append(blockToBlocksScorePat, baSA)
	for {
		//fmt.Printf("[GetMaxScoreArr]bsA: %v\n", bsA)
		maxSc := GetMaxScoreFromBaSA(bsA, visitArr)
		if maxSc.Sc <= 0 {
			break
		}
		/*tk, sc := CheckChainBoundary(maxA, ki.MkInfoArr[maxSc.Idx], maxSc)
		if sc.Sc < maxSc.Sc {
			continue
		}*/
		// add to Max Chain
		visitArr[maxSc.Idx] = true
		bsA[maxSc.Idx] = int(ki.MkInfoArr[maxSc.Idx].Len)
		maxA = append(maxA, ki.MkInfoArr[maxSc.Idx])
		//fmt.Printf("[GetMaxScoreArr]added mki ReadPos: %v, EdgePos: %v, Len: %v, Score: %v\n", ki.MkInfoArr[maxSc.Idx].GetRPos(), ki.MkInfoArr[maxSc.Idx].GetEPos(), ki.MkInfoArr[maxSc.Idx].Len, maxSc.Sc)
		BlockToNearBlocksScoreArr(ki.MkInfoArr, maxSc.Idx, ki.Strand, visitArr, GapCost, MaxGapLen, SeqType, bsA, edgesArr, 0, kmerlen, uint16(SeedLen))
		//blockToBlocksScorePat = append(blockToBlocksScorePat, baSA)
	}
	sort.Sort(MapingKmerInfoArr(maxA))
	arrScore := GetMappingKIArrScore(maxA, GapCost, edgesArr, 0, kmerlen)
	mp.Arr = maxA
	mp.Score = arrScore
	//fmt.Printf("[GetMaxScoreArr]mp.Score: %v, mp.Arr: %v\n", mp.Score, mp.Arr)
	/*if arrScore > maxScore {
			maxScore = arrScore
			maxArr, maxScArr = CleanLowScore(maxArr, maxScArr, maxScore*98/100)
			maxArr = append(maxArr, maxA)
			maxScArr = append(maxScArr, arrScore)
			loopNum = 0
	} else if arrScore > maxScore *98/100{
		maxArr = append(maxArr, maxA)
		maxScArr = append(maxScArr, arrScore)
		loopNum++
	}
	if loopNum >= MinLoopNum {
		break
	}*/
	return
}

func IndexMapingKmerInfoArr(kb []MapingKmerInfo, mk MapingKmerInfo) (idx int) {
	idx = -1
	for i, tk := range kb {
		if tk.GetRPos() == mk.GetRPos() && tk.Len == mk.Len && tk.EID == mk.EID {
			idx = i
			break
		}
	}
	if idx < 0 {
		log.Fatalf("[IndexMapingKmerInfoArr] mk: %v not include kb\n", mk)
	}
	return
}

func GetReverseMapingKmerInfoArr(arr []MapingKmerInfo) []MapingKmerInfo {
	rev := make([]MapingKmerInfo, len(arr))
	for i, mk := range arr {
		rev[len(rev)-1-i] = mk
	}
	return rev
}

/*func GetNextEdgeMaxScore(kb, shareArr []MapingKmerInfo, direction uint8, SeedLen, GapCost, MaxGapLen, MinLoopNum int, SeqType uint8) (mpi MaxPathInfo) {
	var maxA []MapingKmerInfo
	var blockToBlocksScorePat [][]Score
	var arr []int
	mk := shareArr[len(shareArr)-1]
	Idx := IndexMapingKmerInfoArr(kb, mk)
	arr = append(arr, Idx)
	maxA = append(maxA, mk)
	baSA := BlockToNearBlocksScoreArr(kb, Idx, arr, GapCost, MaxGapLen, SeqType)
	blockToBlocksScorePat = append(blockToBlocksScorePat, baSA)
	for {
		maxSc := GetMaxScoreFromBaSAPat(blockToBlocksScorePat, arr)
		if maxSc.Sc <= 0 {
			break
		}
		tk, sc := CheckChainBoundary(maxA, kb[maxSc.Idx], maxSc)
		if sc.Sc < maxSc.Sc {
			continue
		}
		// add to Max Chain
		arr = append(arr, sc.Idx)
		maxA = append(maxA, tk)
		baSA := BlockToNearBlocksScoreArr(kb, sc.Idx, arr, GapCost, MaxGapLen, SeqType)
		blockToBlocksScorePat = append(blockToBlocksScorePat, baSA)
	}
	sort.Sort(MapingKmerInfoArr(maxA))
	arrScore := GetMappingKIArrScore(maxA, MaxGapLen, GapCost)
	if direction == constructdbg.BACKWARD {
		maxA = GetReverseMapingKmerInfoArr(maxA)
	}

	if !reflect.DeepEqual(maxA[len(maxA)-len(shareArr):], shareArr) {
		fmt.Printf("[GetNextEdgeMaxScore] shareArr: %v, not consis with maxA: %v\n", shareArr, maxA)
	} else {
		mpi.Arr = maxA
		mpi.Score = arrScore
	}

	return
}*/

func CheckBound(mkArr []MapingKmerInfo, ReadLen, EdgeLen, kmerlen int, direction uint8) (distR, distE int, shareArr []MapingKmerInfo) {
	mk := mkArr[len(mkArr)-1]
	if direction == constructdbg.FORWARD {
		distR = ReadLen - int(mk.GetRPos()+uint32(mk.Len))
		if mk.GetStrand() == constructdbg.PLUS {
			distE = EdgeLen - int(mk.GetEPos()+uint32(mk.Len))
			for i := len(mkArr) - 1; i >= 0; i-- {
				tk := mkArr[i]
				if int(tk.GetEPos()) > EdgeLen-(kmerlen-1) {
					shareArr = append(shareArr, tk)
				} else {
					break
				}
			}
		} else {
			distE = int(mk.GetEPos())
			for i := len(mkArr) - 1; i >= 0; i-- {
				tk := mkArr[i]
				if tk.GetEPos()+uint32(tk.Len) < uint32(kmerlen-1) {
					shareArr = append(shareArr, tk)
				} else {
					break
				}
			}
		}
	} else {
		distR = int(mk.GetRPos())
		if mk.GetStrand() == constructdbg.PLUS {
			distE = int(mk.GetEPos())
			for i := len(mkArr) - 1; i >= 0; i-- {
				tk := mkArr[i]
				if tk.GetEPos()+uint32(tk.Len) < uint32(kmerlen-1) {
					shareArr = append(shareArr, tk)
				} else {
					break
				}
			}
		} else {
			distE = EdgeLen - int(mk.GetEPos()+uint32(mk.Len))
			for i := len(mkArr) - 1; i >= 0; i-- {
				tk := mkArr[i]
				if int(tk.GetEPos()) > EdgeLen-(kmerlen-1) {
					shareArr = append(shareArr, tk)
				} else {
					break
				}
			}
		}
	}
	return
}

func GetReadBound(mpArr []MaxPathInfo, direction uint8) (ok bool, bd uint32) {
	if direction == constructdbg.FORWARD {
		for _, mp := range mpArr {
			if bd == 0 {
				bd = mp.Arr[len(mp.Arr)-1].GetRPos()
				ok = true
			} else {
				if bd != mp.Arr[len(mp.Arr)-1].GetRPos() {
					ok = false
					break
				}
			}
		}
	} else {
		for _, mp := range mpArr {
			if bd == 0 {
				bd = mp.Arr[0].GetRPos()
				ok = true
			} else {
				if bd != mp.Arr[0].GetRPos() {
					ok = false
					break
				}
			}
		}
	}
	return
}

func CheckMaxScoreBound(RPos uint32, Score int, direction uint8, bsArr []BoundScore) ([]BoundScore, bool) {
	for i := len(bsArr); i >= 0; i-- {
		bs := bsArr[i]
		if direction == constructdbg.FORWARD {
			if RPos < bs.RPos {
				continue
			}
		} else {
			if RPos > bs.RPos {
				continue
			}
		}
		if RPos == bs.RPos {
			if bs.Score <= Score {
				bsArr[i].Score = Score
				return bsArr, true
			} else {
				return bsArr, false
			}
		} else {
			if Score <= bs.Score {
				return bsArr, false
			} else {
				na := make([]BoundScore, len(bsArr)+1)
				copy(na[:i+1], bsArr[:i+1])
				na[i+1].RPos, na[i+1].Score = RPos, Score
				copy(na[i+2:], bsArr[i+1:])
				bsArr = na
				return bsArr, true
			}
		}
	}
	return bsArr, false
}

func InsertSharePat(sharePat [][]MapingKmerInfo, mk MapingKmerInfo) [][]MapingKmerInfo {
	ok := false
	for i, arr := range sharePat {
		if arr[0].EID == mk.EID {
			sharePat[i] = append(sharePat[i], mk)
			ok = true
			break
		}
	}
	if !ok {
		var arr []MapingKmerInfo
		arr = append(arr, mk)
		sharePat = append(sharePat, arr)
	}
	return sharePat
}

func GetShareKBArr(kb, shareArr []MapingKmerInfo, direction uint8) (sharePat [][]MapingKmerInfo) {
	if direction == constructdbg.FORWARD {
		low, high := shareArr[0].GetRPos(), shareArr[len(shareArr)-1].GetRPos()
		idx := 0
		for _, mk := range kb {
			if mk.GetRPos() < low {
				continue
			} else if mk.GetRPos() > high {
				break
			}
			if mk.GetRPos() > shareArr[idx].GetRPos() {
				idx++
				if idx >= len(shareArr) {
					break
				}
			}
			if mk.GetRPos() == shareArr[idx].GetRPos() && mk.Len == shareArr[idx].Len {
				sharePat = InsertSharePat(sharePat, mk)
			}
		}
	} else {
		high, low := shareArr[0].GetRPos(), shareArr[len(shareArr)-1].GetRPos()
		idx := len(shareArr) - 1
		for _, mk := range kb {
			if mk.GetRPos() < low {
				continue
			} else if mk.GetRPos() > high {
				break
			}
			if mk.GetRPos() > shareArr[idx].GetRPos() {
				idx--
				if idx < 0 {
					break
				}
			}
			if mk.GetRPos() == shareArr[idx].GetRPos() && mk.Len == shareArr[idx].Len {
				sharePat = InsertSharePat(sharePat, mk)
			}
		}
	}
	// check sharePat
	idx := 0
	for _, arr := range sharePat {
		if len(arr) != len(shareArr) {
			continue
		}
		ok := true
		strand := arr[0].GetStrand()
		for j := 1; j < len(arr); j++ {
			if strand != arr[j].GetStrand() {
				ok = false
				break
			}
		}
		if ok {
			sharePat[idx] = arr
			idx++
		}
	}
	return
}

func IsConnective(e, ne constructdbg.DBGEdge, nd constructdbg.DBGNode) (ok bool) {
	if constructdbg.IsInComing(nd.EdgeIDIncoming, e.ID) {
		if constructdbg.IsInComing(nd.EdgeIDOutcoming, ne.ID) {
			ok = true
		} else {
			ok = false
		}
	} else {
		if constructdbg.IsInComing(nd.EdgeIDIncoming, ne.ID) {
			ok = true
		} else {
			ok = false
		}
	}
	return ok
}

func CheckDBG(sharePat [][]MapingKmerInfo, edgesArr []constructdbg.DBGEdge, nodesArr []constructdbg.DBGNode, mk MapingKmerInfo, direction uint8) [][]MapingKmerInfo {
	e := edgesArr[mk.EID]
	strand := mk.GetStrand()
	idx := 0
	if direction == constructdbg.FORWARD {
		if strand == constructdbg.PLUS {
			for _, arr := range sharePat {
				ne := edgesArr[arr[0].EID]
				nstrd := arr[0].GetStrand()
				if !IsConnective(e, ne, nodesArr[e.EndNID]) {
					fmt.Printf("[CheckDBG] eID: %v, and neID: %v can't connective\n", e.ID, ne.ID)
					continue
				}
				if GetNextEdgeStrand(e, ne, nodesArr, strand) == nstrd {
					sharePat[idx] = arr
					idx++
				}
			}
		} else {
			for _, arr := range sharePat {
				ne := edgesArr[arr[0].EID]
				nstrd := arr[0].GetStrand()
				if !IsConnective(e, ne, nodesArr[e.StartNID]) {
					fmt.Printf("[CheckDBG] eID: %v, and neID: %v can't connective\n", e.ID, ne.ID)
					continue
				}
				if GetNextEdgeStrand(e, ne, nodesArr, strand) == nstrd {
					sharePat[idx] = arr
					idx++
				}
			}
		}
	} else {
		if strand == constructdbg.PLUS {
			for _, arr := range sharePat {
				ne := edgesArr[arr[0].EID]
				nstrd := arr[0].GetStrand()
				if !IsConnective(e, ne, nodesArr[e.StartNID]) {
					fmt.Printf("[CheckDBG] eID: %v, and neID: %v can't connective\n", e.ID, ne.ID)
					continue
				}
				if GetNextEdgeStrand(e, ne, nodesArr, strand) == nstrd {
					sharePat[idx] = arr
					idx++
				}
			}
		} else {
			for _, arr := range sharePat {
				ne := edgesArr[arr[0].EID]
				nstrd := arr[0].GetStrand()
				if !IsConnective(e, ne, nodesArr[e.EndNID]) {
					fmt.Printf("[CheckDBG] eID: %v, and neID: %v can't connective\n", e.ID, ne.ID)
					continue
				}
				if GetNextEdgeStrand(e, ne, nodesArr, strand) == nstrd {
					sharePat[idx] = arr
					idx++
				}
			}
		}
	}

	return sharePat
}

func AnchorFilter(cBMRIArr []ChainBlockMapLongReadInfo, cb ChainBlockMapLongReadInfo) bool {
	filter := false
	for _, ele := range cBMRIArr {
		if ele.Start <= cb.Start && cb.End <= cb.End {
			alen, blen := int(ele.End)-int(ele.Start), int(cb.End)-int(cb.Start)
			if cb.Score*100/blen < ele.Score*100/alen {
				filter = true
				break
			}
		}
	}
	return filter
}

func GetMaxChainBlockArr(kb []MapingKmerInfo, edgesArr []constructdbg.DBGEdge, nodesArr []constructdbg.DBGNode, ReadLen, SeedLen, GapCost, MaxGapLen, MinKmerScore, kmerlen int, SeqType uint8) (maxPIArr []MaxPathInfo) {
	edgesKIArr := GetEdgesKIArr(kb, edgesArr, MinKmerScore)
	sort.Sort(EdgeKIArr(edgesKIArr))
	//var maxScore int
	//var maxScArr []int
	//var loopNum int
	//var cBMRIArr []ChainBlockMapLongReadInfo
	var eIDArr []constructdbg.DBG_MAX_INT
	var maxPathInfoArr []MaxPathInfo
	//var tmpArr []MaxPathInfo
	// search long read local best mapping to the edge
	for i := len(edgesKIArr) - 1; i >= 0; i-- {
		ki := edgesKIArr[i]
		e := edgesArr[ki.MkInfoArr[0].EID]
		// filter
		/*if e.StartNID == e.EndNID || e.GetTwoEdgesCycleFlag() > 0 {
			continue
		}
		if len(e.Utg.Ks) < kmerlen*2-3 {
			continue
		}*/
		maxPI := GetMaxScoreArr(ki, SeedLen, GapCost, MaxGapLen, SeqType, edgesArr, kmerlen)
		var cb ChainBlockMapLongReadInfo
		lastMK := maxPI.Arr[len(maxPI.Arr)-1]
		cb.Start, cb.End, cb.Score = int(maxPI.Arr[0].GetRPos()), int(lastMK.GetRPos())+int(lastMK.Len), maxPI.Score
		//tmpArr = append(tmpArr, maxPI)
		//if AnchorFilter(cBMRIArr, cb) == false {
		eIDArr = append(eIDArr, e.ID)
		maxPathInfoArr = append(maxPathInfoArr, maxPI)
		//}
	}

	sort.Sort(MaxPathInfoTmpArr(maxPathInfoArr))
	fmt.Printf("[FindONTReadsPath] ReadLen: %v\n", ReadLen)
	for i, mpi := range maxPathInfoArr {
		start, end := mpi.Arr[0].GetRPos(), mpi.Arr[len(mpi.Arr)-1].GetRPos()+uint32(mpi.Arr[len(mpi.Arr)-1].Len)
		fmt.Printf("[FindONTReadsPath] maxPathInfoArr[%v]: Edge ID: %v, EdgeLen: %v, Read region[%v---%v] Len:%v, Edge region[%v---%v], score: %v, Percent: %v\n", i, mpi.Arr[0].EID, len(edgesArr[mpi.Arr[0].EID].Utg.Ks), start, end, end-start, mpi.Arr[0].GetEPos(), mpi.Arr[len(mpi.Arr)-1].GetEPos(), mpi.Score, mpi.Score*100/(int(end)-int(start)))
	}
	/*for y, mki := range ele.Arr {
		var cb ChainBlockMapLongReadInfo
		lastMK := mpi.Arr[len(mpi.Arr)-1]
		cb.Start, cb.End, cb.Score = int(mpi.Arr[0].GetRPos()), int(lastMK.GetRPos())+int(lastMK.Len), mpi.Score
		fmt.Printf("[GetMaxChainBlockArr][%v] EID: %v, start: %v, end: %v, score: %v\n", i, lastMK.EID, cb.Start, cb.End, cb.Score)
	}*/

	// filter low mapping quality edge and filter noly maping middle region of edge
	var highQMaxPathArr []MaxPathInfo
	for i, mpi := range maxPathInfoArr {
		filter := false
		mka, mkb := mpi.Arr[0], mpi.Arr[len(mpi.Arr)-1]
		Start, End := int(mka.GetRPos()), int(mkb.GetRPos())+int(mkb.Len)
		if mpi.Score < SeedLen*2 || End-Start < kmerlen/2 {
			filter = true
			continue
		}
		boundLen := kmerlen
		e := edgesArr[mpi.Arr[0].EID]
		if len(e.Utg.Ks) > 1000 {
			boundLen = 3 * kmerlen
		}
		strand := mpi.Arr[0].GetStrand()
		if strand == constructdbg.PLUS {
			eS, eE := int(mka.GetEPos()), int(mkb.GetEPos())+int(mkb.Len)
			if Start > boundLen && eS > boundLen {
				filter = true
				continue
			}
			if End < ReadLen-boundLen && eE < len(e.Utg.Ks)-boundLen {
				filter = true
				continue
			}
		} else {
			eS, eE := int(mkb.GetEPos()), int(mka.GetEPos())+int(mka.Len)
			if Start > boundLen && len(e.Utg.Ks)-eE > boundLen {
				filter = true
				continue
			}
			if ReadLen-End > boundLen && eS > boundLen {
				filter = true
				continue
			}
		}
		for j, hmpi := range maxPathInfoArr {
			if i == j {
				continue
			}
			hStart, hEnd := int(hmpi.Arr[0].GetRPos()), int(hmpi.Arr[len(hmpi.Arr)-1].GetRPos())+int(hmpi.Arr[len(hmpi.Arr)-1].Len)
			if hStart <= Start && End <= hEnd && mpi.Score*100/(End-Start) < hmpi.Score*100/(hEnd-hStart) {
				filter = true
				fmt.Printf("[GetMaxChainBlockArr] filter EID: %v, start: %v, end: %v\n", mpi.Arr[0].EID, Start, End)
				break
			} else if Start == hStart && hEnd == End && hmpi.Score*100/(hEnd-hStart) == mpi.Score*100/(End-Start) {
				filter = true
				break
			}
		}
		if filter == false {
			highQMaxPathArr = append(highQMaxPathArr, mpi)
		}
	}
	maxPIArr = highQMaxPathArr
	return
}

/*{
	for i := len(edgesKIArr) - 1; i >= 0; i-- {
		ki := edgesKIArr[i]
		mk := kb[ki.Idx]
		var lstack, rstack []MaxPathInfo
		maxPI := GetMaxScoreArr(ki, SeedLen, GapCost, MaxGapLen, MinLoopNum, SeqType)
		maxPI.Arr = GetReverseMapingKmerInfoArr(maxPI.Arr)
		maxPI.Path = append(maxPI.Path, maxPI.Arr[0].EID)
		lstack = append(lstack, maxPI)

		// BACKWARD
		{
			var bsArr []BoundScore
			for len(lstack) > 0 {
				pi := lstack[len(lstack)-1]
				lstack = lstack[:len(lstack)-1]
				for {
					distR, distE, shareArr := CheckBound(pi.Arr, ReadLen, len(edgesArr[pi.Arr[0].EID].Utg.Ks), kmerlen, constructdbg.BACKWARD)
					if distE > kmerlen-1 || distR < kmerlen+distE || len(shareArr) == 0 {
						// add rstack
						pi.Arr = GetReverseMapingKmerInfoArr(pi.Arr)
						rstack = append(rstack, pi)
						break
					}
					sharePat := GetShareKBArr(kb, shareArr, constructdbg.BACKWARD)
					sharePat = CheckDBG(sharePat, edgesArr, nodesArr, shareArr[len(shareArr)-1], constructdbg.BACKWARD)
					if len(sharePat) == 0 {
						// add rstack
						pi.Arr = GetReverseMapingKmerInfoArr(pi.Arr)
						rstack = append(rstack, pi)
						break
					}
					pIArr := make([]MaxPathInfo, len(sharePat))
					//eIDArr := make([]uint32, len(sharePat))
					for j, a := range sharePat {
						pIArr[j] = GetNextEdgeMaxScore(kb, a, constructdbg.BACKWARD, SeedLen, GapCost, MaxGapLen, MinLoopNum, SeqType)
					}
					if ok, bd := GetReadBound(pIArr, constructdbg.BACKWARD); ok {
						maxSc := pIArr[0].Score
						count := 1
						idx := 0
						for m := 1; m < len(pIArr); m++ {
							if pIArr[m].Score > maxSc {
								maxSc = pIArr[m].Score
								count = 1
								idx = m
							} else if pIArr[m].Score == maxSc {
								count++
							}
						}
						if count == 1 {
							pi.Arr = append(pi.Arr, pIArr[idx].Arr...)
							pi.Score = pi.Score + pIArr[idx].Score
							//pi.DistE = pIArr[idx].DistE
							var ok bool
							bsArr, ok = CheckMaxScoreBound(pi.Arr[len(pi.Arr)-1].GetRPos(), pi.Score, constructdbg.BACKWARD, bsArr)
							if ok {
								continue
							} else {
								break
							}
						} else {
							for n, p := range pIArr {
								if p.Score == maxSc {
									var mpi MaxPathInfo
									mpi.Arr = make([]MapingKmerInfo, len(pi.Arr)+len(p.Arr))
									copy(mpi.Arr[:len(pi.Arr)], pi.Arr)
									copy(mpi.Arr[len(pi.Arr):], p.Arr)
									//mpi.DistE = p.DistE
									mpi.Score = pi.Score + p.Score
									//distR, distE := CheckBound(mpi.Arr, ReadLen, constructdbg.FORWARD)
									var ok bool
									bsArr, ok = CheckMaxScoreBound(mpi.Arr[len(mpi.Arr)-1].GetRPos(), mpi.Score, constructdbg.BACKWARD, bsArr)
									if ok {
										lstack = append(lstack, mpi)
									} else {
										continue
									}
								}
							}
							break
						}
					} else {
						for n, p := range pIArr {
							var mpi MaxPathInfo
							mpi.Arr = make([]MapingKmerInfo, len(pi.Arr)+len(p.Arr))
							copy(mpi.Arr[:len(pi.Arr)], pi.Arr)
							copy(mpi.Arr[len(pi.Arr):], p.Arr)
							//mpi.DistE = p.DistE
							mpi.Score = pi.Score + p.Score
							//distR, distE := CheckBound(mpi.Arr, ReadLen, constructdbg.BACKWARD)
							var ok bool
							bsArr, ok = CheckMaxScoreBound(mpi.Arr[len(mpi.Arr)-1].GetRPos(), mpi.Score, constructdbg.BACKWARD, bsArr)
							if ok {
								lstack = append(lstack, mpi)
							} else {
								continue
							}
						}
						break
					}
				}
			}
		}
		//distR, distE := CheckBound(maxPI.Arr,ReadLen, constructdbg.BACKWARD)

		//FORWARD
		{
			var bsArr []BoundScore
			//maxPI.DistE = distE
			//stack = append(stack, maxPI)
			//var bs BoundScore
			//bs.RPos = maxPI.Arr[len(maxPI.Arr)-1].GetRPos() + uint32(maxPI.Arr[len(maxPI.Arr)-1].Len)
			//bs.Score = maxPI.Score
			//bsArr = append(bsArr, bs)
			for len(rstack) > 0 {
				pi := rstack[len(rstack)-1]
				rstack = rstack[:len(rstack)-1]
				for {
					distR, distE, shareArr := CheckBound(pi.Arr, ReadLen, len(edgesArr[pi.Arr[len(pi.Arr)-1].EID].Utg.Ks), kmerlen, constructdbg.FORWARD)
					if distE > kmerlen-1 || distR < kmerlen+distE || len(shareArr) == 0 {
						// add to the maxPIArr
						maxPIArr = append(maxPIArr, pi)
						break
					}
					sharePat := GetShareKBArr(kb, shareArr, constructdbg.FORWARD)
					sharePat = CheckDBG(sharePat, edgesArr, nodesArr, shareArr[len(shareArr)-1], constructdbg.FORWARD)
					if len(sharePat) == 0 {
						// add to maxPIArr
						maxPIArr = append(maxPIArr, pi)
						break
					}
					pIArr := make([]MaxPathInfo, len(sharePat))
					for j, a := range sharePat {
						pIArr[j] = GetNextEdgeMaxScore(kb, a, constructdbg.FORWARD, SeedLen, GapCost, MaxGapLen, MinLoopNum, SeqType)
					}
					if ok, bd := GetReadBound(pIArr, constructdbg.FORWARD); ok {
						maxSc := pIArr[0].Score
						count := 1
						idx := 0
						for m := 1; m < len(pIArr); m++ {
							if pIArr[m].Score > maxSc {
								maxSc = pIArr[m].Score
								count = 1
								idx = m
							} else if pIArr[m].Score == maxSc {
								count++
							}
						}
						if count == 1 {
							pi.Arr = append(pi.Arr, pIArr[idx].Arr...)
							pi.Score = pi.Score + pIArr[idx].Score
							//pi.DistE = pIArr[idx].DistE
							var ok bool
							bsArr, ok = CheckMaxScoreBound(pi.Arr[len(pi.Arr)-1].GetRPos()+uint32(pi.Arr[len(pi.Arr)-1].Len), pi.Score, constructdbg.FORWARD, bsArr)
							if ok {
								continue
							} else {
								break
							}
						} else {
							for n, p := range pIArr {
								if p.Score == maxSc {
									var mpi MaxPathInfo
									mpi.Arr = make([]MapingKmerInfo, len(pi.Arr)+len(p.Arr))
									copy(mpi.Arr[:len(pi.Arr)], pi.Arr)
									copy(mpi.Arr[len(pi.Arr):], p.Arr)
									//mpi.DistE = p.DistE
									mpi.Score = pi.Score + p.Score
									//distR, distE := CheckBound(mpi.Arr, ReadLen, constructdbg.FORWARD)
									var ok bool
									bsArr, ok = CheckMaxScoreBound(mpi.Arr[len(mpi.Arr)-1].GetRPos()+uint32(mpi.Arr[len(mpi.Arr)-1].Len), mpi.Score, constructdbg.FORWARD, bsArr)
									if ok {
										rstack = append(rstack, mpi)
									} else {
										continue
									}
								}
							}
							break
						}
					} else {
						for n, p := range pIArr {
							var mpi MaxPathInfo
							mpi.Arr = make([]MapingKmerInfo, len(pi.Arr)+len(p.Arr))
							copy(mpi.Arr[:len(pi.Arr)], pi.Arr)
							copy(mpi.Arr[len(pi.Arr):], p.Arr)
							//mpi.DistE = p.DistE
							mpi.Score = pi.Score + p.Score
							//distR, distE := CheckBound(mpi.Arr, ReadLen, constructdbg.FORWARD)
							var ok bool
							bsArr, ok = CheckMaxScoreBound(mpi.Arr[len(mpi.Arr)-1].GetRPos()+uint32(mpi.Arr[len(mpi.Arr)-1].Len), mpi.Score, constructdbg.FORWARD, bsArr)
							if ok {
								rstack = append(rstack, mpi)
							} else {
								continue
							}
						}
						break
					}
				}
			}
		}
	}
}*/

func GetEdgeRegionSeq(e constructdbg.DBGEdge, mpi MaxPathInfo) (edgeRegionSeq []byte) {
	// check mpi.Arr
	if mpi.Arr[0].EID != uint32(e.ID) || mpi.Arr[len(mpi.Arr)-1].EID != uint32(e.ID) {
		log.Fatalf("[GetEdgeRegionSeq] e.ID: %v != mpi.Arr: %v\n", e.ID, mpi.Arr)
	}
	strand := mpi.Arr[0].GetStrand()
	last := mpi.Arr[len(mpi.Arr)-1]
	fmt.Printf("[GetEdgeRegionSeq]eID: %v, arr[0] ReadPos: %v, EdgePos: %v, strand: %v, Len: %v, arr[last] ReadPos: %v, EdgePos: %v, strand: %v, Len: %v, len(edge): %v\n", e.ID, mpi.Arr[0].GetRPos(), mpi.Arr[0].GetEPos(), mpi.Arr[0].GetStrand(), mpi.Arr[0].Len, last.GetRPos(), last.GetEPos(), last.GetStrand(), last.Len, len(e.Utg.Ks))
	if strand == constructdbg.PLUS {
		start, end := mpi.Arr[0].GetEPos(), last.GetEPos()+uint32(last.Len)
		edgeRegionSeq = append(edgeRegionSeq, e.Utg.Ks[start:end]...)
	} else {
		start, end := last.GetEPos(), mpi.Arr[0].GetEPos()+uint32(mpi.Arr[0].Len)
		edgeRegionSeq = append(edgeRegionSeq, constructdbg.GetReverseCompByteArr(e.Utg.Ks[start:end])...)
	}
	return
}
func GetDirectionEdgesSeq(edgesArr []constructdbg.DBGEdge, nodesArr []constructdbg.DBGNode, direction uint8, kb []MapingKmerInfo, kmerlen, extendLen int) (edgesSeq []byte) {
	mkS := kb[0]
	mkE := kb[len(kb)-1]
	fmt.Printf("[GetDirectionEdgesSeq] mkE EID: %v, EdgePos: %v,Len: %v, strand: %v\n", mkE.EID, mkE.GetEPos(), mkE.Len, mkE.GetStrand())
	if direction == constructdbg.FORWARD {
		// start edge
		e := edgesArr[mkS.EID]
		var finished bool
		if mkS.GetStrand() == constructdbg.PLUS {
			edgesSeq = append(edgesSeq, e.Utg.Ks[mkS.GetEPos():]...)
		} else {
			edgesSeq = append(edgesSeq, constructdbg.GetReverseCompByteArr(e.Utg.Ks[:mkS.GetEPos()+uint32(mkS.Len)])...)
		}

		// middle edges
		for i := 1; i < len(kb)-1; i++ {
			mk := kb[i]
			if mk.EID == uint32(e.ID) {
				continue
			}
			if mk.EID == mkE.EID {
				break
			}
			e = edgesArr[mk.EID]
			extLen := len(e.Utg.Ks) - (kmerlen - 1)
			finished = false
			if len(edgesSeq)+(len(e.Utg.Ks)-(kmerlen-1)) >= extendLen {
				finished = true
				extLen = extendLen - len(edgesSeq)
			}
			if extLen <= 0 {
				finished = true
				break
			}
			if mk.GetStrand() == constructdbg.PLUS {
				edgesSeq = append(edgesSeq, e.Utg.Ks[kmerlen-1:kmerlen-1+extLen]...)
			} else {
				edgesSeq = append(edgesSeq, constructdbg.GetReverseCompByteArr(e.Utg.Ks[len(e.Utg.Ks)-(kmerlen-1)-extLen:len(e.Utg.Ks)-(kmerlen-1)])...)
			}
			if finished {
				break
			}
		}
		if finished {
			return
		}
		// add last edge
		e = edgesArr[mkE.EID]
		extLen := extendLen - len(edgesSeq)
		if extLen > len(e.Utg.Ks)-(kmerlen-1) {
			extLen = len(e.Utg.Ks) - (kmerlen - 1)
		}
		if mkE.GetStrand() == constructdbg.PLUS {
			/*if int(mkE.GetEPos())+int(mkE.Len) < kmerlen-1 {
				dl := kmerlen - 1 - (int(mkE.GetEPos()) + int(mkE.Len))
				edgesSeq = edgesSeq[:len(edgesSeq)-dl]
			} else {*/

			edgesSeq = append(edgesSeq, e.Utg.Ks[kmerlen-1:kmerlen-1+extLen]...)
			//}
		} else {
			/*if int(mkE.GetEPos()) > len(e.Utg.Ks)-(kmerlen-1) {
				dl := int(mkE.GetEPos()) - (len(e.Utg.Ks) - (kmerlen - 1))
				edgesSeq = edgesSeq[:len(edgesSeq)-dl]
			} else {*/

			edgesSeq = append(edgesSeq, constructdbg.GetReverseCompByteArr(e.Utg.Ks[len(e.Utg.Ks)-(kmerlen-1)-extLen:len(e.Utg.Ks)-(kmerlen-1)])...)
			//}
		}
	} else {
		// start edge
		e := edgesArr[mkS.EID]
		var finished bool
		if mkS.GetStrand() == constructdbg.PLUS {
			edgesSeq = append(edgesSeq, constructdbg.GetReverseByteArr(e.Utg.Ks[:mkS.GetEPos()+uint32(mkS.Len)])...)
		} else {
			edgesSeq = append(edgesSeq, constructdbg.GetCompByteArr(e.Utg.Ks[mkS.GetEPos():])...)
		}

		// middle edges
		for i := 1; i < len(kb)-1; i++ {
			mk := kb[i]
			if mk.EID == uint32(e.ID) {
				continue
			}
			if mk.EID == mkE.EID {
				break
			}
			e = edgesArr[mk.EID]
			extLen := len(e.Utg.Ks) - (kmerlen - 1)
			finished = false
			if len(edgesSeq)+(len(e.Utg.Ks)-(kmerlen-1)) >= extendLen {
				finished = true
				extLen = extendLen - len(edgesSeq)
			}
			if extLen <= 0 {
				finished = true
				break
			}
			if mkS.GetStrand() == constructdbg.PLUS {
				edgesSeq = append(edgesSeq, constructdbg.GetReverseByteArr(e.Utg.Ks[len(e.Utg.Ks)-(kmerlen-1)-extLen:len(e.Utg.Ks)-(kmerlen-1)])...)
			} else {
				edgesSeq = append(edgesSeq, constructdbg.GetCompByteArr(e.Utg.Ks[(kmerlen-1):(kmerlen-1)+extLen])...)
			}
			if finished {
				break
			}
		}
		if finished {
			return
		}

		// add last edge
		e = edgesArr[mkE.EID]
		extLen := extendLen - len(edgesSeq)
		if extLen > len(e.Utg.Ks)-(kmerlen-1) {
			extLen = len(e.Utg.Ks) - (kmerlen - 1)
		}
		if mkE.GetStrand() == constructdbg.PLUS {
			/*if int(mkE.GetEPos()) > len(e.Utg.Ks)-(kmerlen-1) {
				dl := int(mkE.GetEPos()) - (len(e.Utg.Ks) - (kmerlen - 1))
				edgesSeq = edgesSeq[:len(edgesSeq)-dl]
			} else {*/
			edgesSeq = append(edgesSeq, constructdbg.GetReverseByteArr(e.Utg.Ks[len(e.Utg.Ks)-(kmerlen-1)-extLen:len(e.Utg.Ks)-(kmerlen-1)])...)
			//}
		} else {
			/*if (kmerlen - 1) > int(mkE.GetEPos())+int(mkE.Len) {
				dl := (kmerlen - 1) - (int(mkE.GetEPos()) + int(mkE.Len))
				edgesSeq = edgesSeq[:len(edgesSeq)-dl]
			} else {*/
			edgesSeq = append(edgesSeq, constructdbg.GetCompByteArr(e.Utg.Ks[(kmerlen-1):(kmerlen-1)+extLen])...)
			//}
		}
	}
	return
}

func GetEdgesSeq(edgesArr []constructdbg.DBGEdge, nodesArr []constructdbg.DBGNode, kb []MapingKmerInfo, kmerlen int) (edgesSeq []byte, path []constructdbg.DBG_MAX_INT, chainArr []Chain) {
	chainArr = make([]Chain, len(kb))
	mkS := kb[0]
	mkE := kb[len(kb)-1]
	var startY, startX uint32
	// start Anchor
	{
		e := edgesArr[mkS.EID]
		if mkS.GetStrand() == constructdbg.PLUS {
			edgesSeq = append(edgesSeq, e.Utg.Ks[mkS.GetEPos():]...)
			startX = mkS.GetEPos()
		} else {
			edgesSeq = append(edgesSeq, constructdbg.GetReverseCompByteArr(e.Utg.Ks[:mkS.GetEPos()+uint32(mkS.Len)])...)
			startX = mkS.GetEPos() + uint32(mkS.Len)
		}
		path = append(path, e.ID)
		chainArr[0].Len = uint32(mkS.Len)
		startY = mkS.GetRPos()
	}
	// middle
	{
		// determine boundary
		i := 1
		for ; i < len(kb)-1; i++ {
			mk := kb[i]
			if mk.EID == uint32(path[len(path)-1]) {
				if mk.GetStrand() == constructdbg.PLUS {
					chainArr[i].X = mk.GetEPos() - startX
				} else {
					chainArr[i].X = startX - (mk.GetEPos() + uint32(mk.Len))
				}
				chainArr[i].Len = uint32(mk.Len)
				chainArr[i].Y = mk.GetRPos() - startY
				continue
			} else {
				break
			}
		}
		j := len(kb) - 1
		for ; j > 0; j-- {
			mk := kb[j]
			if mk.EID == mkE.EID {
				continue
			} else {
				break
			}
		}
		for m := i; m <= j; m++ {
			mk := kb[m]
			{
				if mk.GetStrand() == constructdbg.PLUS {
					chainArr[i].X = uint32(len(edgesSeq)) + (mk.GetEPos() - (uint32(kmerlen) - 1))
				} else {
					chainArr[i].X = startX - (mk.GetEPos() + uint32(mk.Len))
				}
				chainArr[i].Len = uint32(mk.Len)
				chainArr[i].Y = mk.GetRPos() - startY
			}
			if mk.EID == uint32(path[len(path)-1]) {
				continue
			}
			e := edgesArr[mk.EID]
			if mk.GetStrand() == constructdbg.PLUS {
				edgesSeq = append(edgesSeq, e.Utg.Ks[kmerlen-1:]...)
			} else {
				edgesSeq = append(edgesSeq, constructdbg.GetReverseCompByteArr(e.Utg.Ks[:len(e.Utg.Ks)-(kmerlen-1)])...)
			}
			path = append(path, e.ID)
		}
	}
	// end Anthor
	{
		e := edgesArr[mkE.EID]
		if mkE.GetStrand() == constructdbg.PLUS {
			edgesSeq = append(edgesSeq, e.Utg.Ks[kmerlen-1:mkE.GetEPos()+uint32(mkE.Len)]...)
		} else {
			edgesSeq = append(edgesSeq, constructdbg.GetReverseCompByteArr(e.Utg.Ks[mkE.GetEPos():len(e.Utg.Ks)-(kmerlen-1)])...)
		}
		path = append(path, e.ID)
	}

	return
}
func GetReadSeq(seq []byte, start, end uint32) (readSeq []byte) {
	readSeq = seq[start:end]
	return
}

type ExtendPathInfo struct {
	MpArr               []MapingKmerInfo
	Path                []constructdbg.DBG_MAX_INT
	Path1               []constructdbg.DBG_MAX_INT
	Strand              bool
	Coming              bool
	ExtendLen, FlankLen int // FlankLen note this edge unseed flank boundary length
	Score               int
}

func GetEdgeComing(e constructdbg.DBGEdge, direction uint8, strand bool, nodesArr []constructdbg.DBGNode, coming bool) bool {
	if e.StartNID > 0 && e.StartNID == e.EndNID {
		nd := nodesArr[e.StartNID]
		if constructdbg.IsInComing(nd.EdgeIDIncoming, e.ID) && constructdbg.IsInComing(nd.EdgeIDOutcoming, e.ID) {
			return coming
		} else {
			return !coming
		}
	}
	var nd constructdbg.DBGNode
	if strand == constructdbg.PLUS {
		if direction == constructdbg.FORWARD {
			if e.EndNID > 0 {
				nd = nodesArr[e.EndNID]
			}
		} else {
			if e.StartNID > 0 {
				nd = nodesArr[e.StartNID]
			}
		}
	} else {
		if direction == constructdbg.FORWARD {
			if e.StartNID > 0 {
				nd = nodesArr[e.StartNID]
			}
		} else {
			if e.EndNID > 0 {
				nd = nodesArr[e.EndNID]
			}
		}
	}
	if constructdbg.IsInComing(nd.EdgeIDIncoming, e.ID) {
		coming = true
	} else {
		coming = false
	}
	return coming
}

func GetNearDirectionEdgeIDArr(e constructdbg.DBGEdge, direction uint8, strand, coming bool, nodesArr []constructdbg.DBGNode) (eArr []constructdbg.DBG_MAX_INT) {
	var nd constructdbg.DBGNode
	if strand == constructdbg.PLUS {
		if direction == constructdbg.FORWARD {
			if e.EndNID <= 0 {
				return
			}
			nd = nodesArr[e.EndNID]
		} else {
			if e.StartNID <= 0 {
				return
			}
			nd = nodesArr[e.StartNID]
		}
	} else {
		if direction == constructdbg.FORWARD {
			if e.StartNID <= 0 {
				return
			}
			nd = nodesArr[e.StartNID]
		} else {
			if e.EndNID <= 0 {
				return
			}
			nd = nodesArr[e.EndNID]
		}
	}
	eArr = constructdbg.GetNearEdgeIDArr(nd, e.ID, coming)
	return
}

func GetMappingDBGFlankLen(epi ExtendPathInfo, direction uint8, eLen int) (flankLen int) {
	anchor := epi.MpArr[len(epi.MpArr)-1]
	if epi.Strand == constructdbg.PLUS {
		if direction == constructdbg.FORWARD {
			flankLen = eLen - (int(anchor.GetEPos()) + int(anchor.Len))
		} else {
			flankLen = int(anchor.GetEPos())
		}
	} else {
		if direction == constructdbg.FORWARD {
			flankLen = int(anchor.GetEPos())
		} else {
			flankLen = eLen - (int(anchor.GetEPos()) + int(anchor.Len))
		}
	}
	return
}

func IsInChainBlockMapLongReadInfoArr(cBMRIArr []ChainBlockMapLongReadInfo, start, end int) bool {
	for _, cb := range cBMRIArr {
		if cb.Start < start && start < cb.End {
			return true
		} else if cb.Start < end && end < cb.End {
			return false
		}
	}

	return false
}

func LocalMapingKmerInfoArr(kb []MapingKmerInfo, mki MapingKmerInfo) int {
	idx := -1
	for i, c := range kb {
		if c.EID == mki.EID && c.GetRPos() <= mki.GetRPos() && mki.GetRPos()+uint32(mki.Len) <= c.GetRPos()+uint32(c.Len) {
			if c.GetRStrand() == mki.GetRStrand() && c.GetEStrand() == mki.GetEStrand() {
				if c.GetEPos() <= mki.GetEPos() && mki.GetEPos()+uint32(mki.Len) <= c.GetEPos()+uint32(c.Len) {
					idx = i
					break
				}
			}
		}
	}
	if idx < 0 {
		log.Fatalf("[LocalMapingKmerInfoArr] not found mki: %v in kb\n", mki)
	}
	return idx
}

func GetEdgeMKIArr(kb []MapingKmerInfo, idx int, eID uint32, direction uint8, strand bool, distance int) (mpArr []MapingKmerInfo) {
	cen := kb[idx]
	mpArr = append(mpArr, cen)
	//fmt.Printf("[GetEdgeMKIArr]eID: %v, idx: %v,cen RPos: %v, direction: %v, strand: %v, distance: %v\n", eID, idx, cen.GetRPos(), direction, strand, distance)
	if direction == constructdbg.FORWARD {
		for i := idx + 1; i < len(kb); i++ {
			mk := kb[i]
			if int(mk.GetRPos()) > int(cen.GetRPos())+int(cen.Len)+distance {
				break
			}
			//fmt.Printf("[GetEdgeMKIArr] kb[%v]: %v, RPos: %v, strand: %v, EID: %v\n", i, mk, mk.GetRPos(), mk.GetStrand(), mk.EID)
			if mk.EID != eID || mk.GetStrand() != strand {
				continue
			}
			mpArr = append(mpArr, mk)
		}
	} else {
		for i := idx - 1; i > 0; i-- {
			mk := kb[i]
			//fmt.Printf("[GetEdgeMKIArr] kb[%v]: %v, RPos: %v, EID: %v\n", i, mk, mk.GetRPos(), mk.EID)
			if int(mk.GetRPos()) < int(cen.GetRPos())-distance {
				break
			}
			if mk.EID != eID || mk.GetStrand() != strand {
				continue
			}
			mpArr = append(mpArr, mk)
		}
	}
	sort.Sort(MapingKmerInfoArr(mpArr))
	return
}

func RevMappingKmerInfoArr(mpArr []MapingKmerInfo) {
	for i := 0; i < len(mpArr)/2; i++ {
		tmp := mpArr[i]
		mpArr[i] = mpArr[len(mpArr)-1-i]
		mpArr[len(mpArr)-1-i] = tmp
	}
	return
}

func ExtendPathChain(mpArr []MapingKmerInfo, flankLen, kmerlen int, direction uint8, ne constructdbg.DBGEdge, kb []MapingKmerInfo, strand bool, score int, edgesArr []constructdbg.DBGEdge, SeedLen uint16) (newArr []MapingKmerInfo, newScore int) {
	//extendMpArr = append(extendMpArr, mpArr...)
	newScore = score
	var extendMpArr []MapingKmerInfo
	lastMki := mpArr[len(mpArr)-1]
	idx := LocalMapingKmerInfoArr(kb, lastMki)
	ekb := GetEdgeMKIArr(kb, idx, uint32(ne.ID), direction, strand, flankLen+(len(ne.Utg.Ks)-(kmerlen-1))*7/5+50)
	if len(ekb) <= 1 {
		newArr = append(newArr, mpArr...)
		return
	}
	//fmt.Printf("[ExtendPathChain] ekb: %v\n", ekb)
	sort.Sort(MapingKmerInfoArr(ekb))
	var ekbIdx int
	if direction == constructdbg.FORWARD {
		ekbIdx = 0
	} else {
		ekbIdx = len(ekb) - 1
	}
	//PrintAddedMpArr(ekb)
	visitArr := make([]bool, len(ekb))
	bsA := make([]int, len(ekb))
	//fmt.Printf("[ExtendPathChain] visitArr: %v\n\tekbIdx: %v\n", visitArr, ekbIdx)
	visitArr[ekbIdx] = true
	extendMpArr = append(extendMpArr, lastMki)
	bsA[ekbIdx] = int(lastMki.Len)
	BlockToNearBlocksScoreArrCrossEdges(ekb, ekbIdx, strand, visitArr, 1, 500, uint8(3), bsA, direction, flankLen, kmerlen, edgesArr, SeedLen)
	for {
		maxSc := GetMaxScoreFromBaSA(bsA, visitArr)
		if maxSc.Sc <= 0 {
			break
		}
		visitArr[maxSc.Idx] = true
		bsA[maxSc.Idx] = int(ekb[maxSc.Idx].Len)
		extendMpArr = append(extendMpArr, ekb[maxSc.Idx])
		//fmt.Printf("[ExtendPathChain] bsA: %v\n\t\tadded mk: %v\n", bsA, ekb[maxSc.Idx])
		BlockToNearBlocksScoreArrCrossEdges(ekb, maxSc.Idx, strand, visitArr, 1, 500, uint8(3), bsA, direction, flankLen, kmerlen, edgesArr, SeedLen)
	}
	sort.Sort(MapingKmerInfoArr(extendMpArr))
	//PrintAddedMpArr(extendMpArr)
	if direction == constructdbg.BACKWARD {
		RevMappingKmerInfoArr(extendMpArr)
	}
	//PrintAddedMpArr(extendMpArr)
	//fmt.Printf("[ExtendPathChain] mpArr: %v\n", mpArr)
	//fmt.Printf("[ExtendPathChain] exendMpArr: %v\n", extendMpArr)

	newScore += GetMappingKIArrScore(extendMpArr, 1, edgesArr, flankLen, kmerlen) - int(lastMki.Len)
	newArr = append(newArr, mpArr...)
	newArr = append(newArr, extendMpArr[1:]...)
	//fmt.Printf("[ExtendPathChain] exendMpArr: %v\n", extendMpArr)
	//fmt.Printf("[ExtendPathChain] newArr: %v\n", newArr)

	/*//startRPos := int(lastMki.GetRPos()) + int(lastMki.Len)
		for i := idx + 1; i < len(kb); i++ {
			mki := kb[i]
			lRPos := int(lastMki.GetRPos())+int(lastMki.Len)
			if mki.EID != uint32(ne.ID) || strand != mki.GetStrand() {
				continue
			}
			if int(mki.GetRPos()) < lRPos {
				continue
			}
			if int(mki.GetRPos())-lRPos > 500 {
				break
			}
			dR := int(mki.GetRPos()) - lRPos
			var dE int
			if strand == constructdbg.PLUS {
				dE = int(mki.GetEPos()) - (int(lastMki.GetEPos()) + int(lastMki.Len))
			} else {
				dE = int(lastMki.GetEPos()) - (int(mki.GetEPos()) + int(mki.Len))
			}
			if dE <= 0 {
				continue
			}
			min := constructdbg.Min(dR, dE)
			gapLen := AbsInt(dR-dE)
			if gapLen > min/5+10 {
				continue
			}
			a := constructdbg.Min(min, int(mki.Len))
			var b float64
			if gapLen > 1 {
				b = float64(0.005)*float64(mki.Len)*float64(gapLen) + float64(0.5)*math.Log2(float64(gapLen))
			} else {
				b = 1
			}
			sc = constructdbg.Min(a-int(b), int(nc.Len))
				sc := int(mki.Len) - AbsInt(dR-dE)
				if sc > 0 {
					extendMpArr = append(extendMpArr, mki)
					lastMki = mki
					newScore += sc
				}
			}
		}
	} else { // direction == constructdbg.BACKWARD
		lastMki := extendMpArr[len(extendMpArr)-1]
		idx := LocalMapingKmerInfoArr(kb, lastMki)
		startRPos := int(lastMki.GetRPos())
		for i := idx - 1; i >= 0; i-- {
			mki := kb[i]
			if int(mki.GetRPos()) < startRPos-(flankLen+len(ne.Utg.Ks)-(kmerlen-1)+len(ne.Utg.Ks)/10) {
				break
			}
			if mki.EID != uint32(ne.ID) || strand != mki.GetStrand() {
				continue
			}
			dR := int(lastMki.GetRPos()) - (int(mki.GetRPos()) + int(mki.Len))
			var dE int
			if strand == constructdbg.PLUS {
				dE = int(lastMki.GetEPos()) - (int(mki.GetEPos()) + int(mki.Len))
			} else {
				dE = int(mki.GetEPos()) - (int(lastMki.GetEPos()) + int(lastMki.Len))
			}
			if dE <= 0 {
				continue
			}
			min := constructdbg.Min(dR, dE)
			if AbsInt(dR-dE) < min/10+8 {
				sc := int(mki.Len) - AbsInt(dR-dE)
				if sc > 0 {
					extendMpArr = append(extendMpArr, mki)
					lastMki = mki
					newScore += sc
				}
			}
		}
	}*/
	return
}

func GetBestAlignmentPath(bestArr []ExtendPathInfo, direction uint8, edgesArr []constructdbg.DBGEdge, nodesArr []constructdbg.DBGNode, read constructcf.ReadInfo, SeedLen, kmerlen, minMapLen int) (bestEPI ExtendPathInfo) {
	var bestNum int
	for i, epi := range bestArr {
		fmt.Printf("[GetBestAlignmentPath] bestArr[%v]: score: %v, mapLen: %v,RPos: %v, path: %v\n", i, epi.Score, epi.ExtendLen-epi.FlankLen, epi.MpArr[len(epi.MpArr)-1].GetRPos(), epi.Path)
		mkia, mkib := epi.MpArr[0], epi.MpArr[len(epi.MpArr)-1]
		var readRegionSeq []byte
		if direction == constructdbg.FORWARD {
			start, _ := int(mkia.GetRPos()), int(mkib.GetRPos())+int(mkib.Len)
			if start+minMapLen > len(read.Seq) {
				minMapLen = len(read.Seq) - start
			}
			readRegionSeq = read.Seq[start : start+minMapLen]
		} else {
			start, _ := int(mkia.GetRPos())+int(mkia.Len), int(mkib.GetRPos())
			if minMapLen > start {
				minMapLen = start
			}
			readRegionSeq = constructdbg.GetReverseByteArr(read.Seq[start-minMapLen : start])
		}

		edgesSeq := GetDirectionEdgesSeq(edgesArr, nodesArr, direction, epi.MpArr, kmerlen, minMapLen*6/5)
		chainArr, _ := ComputingChainArr(edgesSeq, readRegionSeq, uint32(SeedLen))
		fmt.Printf("[GetBestAlignmentPath] base alignment found chainArr length: %v\n", len(chainArr))
		if len(chainArr) == 0 {
			continue
		}
		cg := DPLocalAlign(edgesSeq, readRegionSeq, chainArr)
		fmt.Printf("[GetBestAlignmentPath][%v] cg Mch: %v, Mis: %v, Del: %v, Ins: %v\n", i, cg.Mch, cg.Mis, cg.Del, cg.Ins)
		sc := int(cg.Mch) - int(cg.Mis) - int(cg.Del) - int(cg.Ins)
		if sc > bestEPI.Score {
			bestEPI = epi
			bestEPI.Score = sc
			bestNum = 1
		} else if sc == bestEPI.Score {
			a, b := epi.ExtendLen-epi.FlankLen, bestEPI.ExtendLen-bestEPI.FlankLen
			if a != b {
				if a > b {
					bestEPI = epi
					bestEPI.Score = sc
				}
			} else {
				bestNum++
			}
		}
	}

	if bestNum > 1 {
		var t ExtendPathInfo
		bestEPI = t
	}
	return
}

func GetBestEPathInfo(highQualEPathInfoArr []ExtendPathInfo, edgesArr []constructdbg.DBGEdge, nodesArr []constructdbg.DBGNode, read constructcf.ReadInfo, SeedLen, kmerlen int, direction uint8) (bestEPI ExtendPathInfo) {
	if len(highQualEPathInfoArr) == 0 {
		return
	}
	if len(highQualEPathInfoArr) == 1 {
		bestEPI = highQualEPathInfoArr[0]
		return
	}
	for i, epi := range highQualEPathInfoArr {
		l := epi.ExtendLen - epi.FlankLen
		fmt.Printf("[GetBestEPathInfo]high[%v]len(epi.MpArr): %v, Score: %v,Len: %v, Percent: %v, ExtendLen: %v, FlankLen: %v, Strand: %v, Coming: %v, Path: %v\n", i, len(epi.MpArr), epi.Score, l, epi.Score*100/l, epi.ExtendLen, epi.FlankLen, epi.Strand, epi.Coming, epi.Path)
	}
	var bestArr []ExtendPathInfo
	bestScoreNum := 0
	var maxScoreArr []ScorePos
	for _, epi := range highQualEPathInfoArr {
		var sp ScorePos
		sp.Pos = epi.ExtendLen - epi.FlankLen
		sp.Score = epi.Score
		maxScoreArr = AddedToMaxScoreArr(maxScoreArr, sp)
	}
	//fmt.Printf("[GetBestEPathInfo]maxScoreArr: %v\n", maxScoreArr)

	minMaplen := math.MaxInt32
	for _, epi := range highQualEPathInfoArr {
		var sp ScorePos
		sp.Pos = epi.ExtendLen - epi.FlankLen
		sp.Score = epi.Score
		if BiggerThanMaxScoreArr(maxScoreArr, sp) && len(epi.MpArr) > 1 {
			epi = DeleteUnmapPath(epi, direction, edgesArr, kmerlen)
			bestEPI = epi
			bestScoreNum++
			bestArr = append(bestArr, epi)
			if sp.Pos < minMaplen {
				minMaplen = sp.Pos
			}
		}
	}
	if minMaplen < 2*kmerlen {
		var tmp ExtendPathInfo
		bestEPI = tmp
		return
	}

	if bestScoreNum > 1 {
		fmt.Printf("[GetBestEPathInfo] minMapLen: %v\n", minMaplen)
		var maxArr []ExtendPathInfo
		for j, epi := range bestArr {
			//fmt.Printf("[GetBestEPathInfo] bestArr[%v]: score: %v, mapLen: %v,RPos: %v, path: %v\n", j, epi.Score, epi.ExtendLen-epi.FlankLen, epi.MpArr[len(epi.MpArr)-1].GetRPos(), epi.Path)
			ok := true
			epiMapLen := epi.ExtendLen - epi.FlankLen
			for x, ele := range bestArr {
				eleMapLen := ele.ExtendLen - ele.FlankLen
				if x != j && epiMapLen < eleMapLen && len(epi.Path) < len(ele.Path) {
					if reflect.DeepEqual(epi.Path, ele.Path[:len(epi.Path)]) {
						ok = false
						break
					}
				}
				if x != j && epiMapLen <= eleMapLen {
					if epi.Score*100/epiMapLen < ele.Score*100/eleMapLen {
						ok = false
						break
					}
				}
			}
			if ok {
				maxArr = append(maxArr, epi)
			}
		}
		if len(maxArr) == 1 {
			bestEPI = maxArr[0]
		} else if len(maxArr) > 1 && len(maxArr) < 20 {
			bestEPI = GetBestAlignmentPath(maxArr, direction, edgesArr, nodesArr, read, SeedLen, kmerlen, minMaplen)
		} else {
			var tmp ExtendPathInfo
			bestEPI = tmp
		}
	}
	return
}

func MuskBubble(eArr []constructdbg.DBG_MAX_INT, stk []ExtendPathInfo) {
	if len(stk) < 2 {
		return
	}
	s1, s2 := stk[len(stk)-2], stk[len(stk)-1]
	if s1.Path[len(s1.Path)-1] == eArr[0] && s2.Path[len(s2.Path)-1] == eArr[1] {
		if s1.Score < s2.Score {
			stk[len(stk)-2].Score = 0
		} else {
			stk[len(stk)-1].Score = 0
		}
	}
	return
}

func DeleteUnmapPath(bestEPI ExtendPathInfo, direction uint8, edgesArr []constructdbg.DBGEdge, kmerlen int) ExtendPathInfo {
	fl := bestEPI.FlankLen
	el := bestEPI.ExtendLen
	if len(bestEPI.MpArr) == 0 {
		return bestEPI
	}
	lastMki := bestEPI.MpArr[len(bestEPI.MpArr)-1]
	eID := lastMki.EID
	i := len(bestEPI.Path) - 1
	for ; i > 0; i-- {
		e := edgesArr[bestEPI.Path[i]]
		if len(e.Utg.Ks)-(kmerlen-1) < fl {
			fl -= (len(e.Utg.Ks) - (kmerlen - 1))
			el -= (len(e.Utg.Ks) - (kmerlen - 1))
			continue
		}

		if uint32(e.ID) == eID {
			if len(e.Utg.Ks)-fl < 2*kmerlen {
				i--
			}
			break
		}
	}
	bestEPI.Path = bestEPI.Path[:i+1]
	if len(bestEPI.Path1) > len(bestEPI.Path) {
		bestEPI.Path1 = bestEPI.Path1[:len(bestEPI.Path)]
	}
	bestEPI.FlankLen = fl
	bestEPI.ExtendLen = el

	return bestEPI
}

func MergeTwoFlankMapingPath(bestEPIB, bestEPIF ExtendPathInfo, mpi MaxPathInfo, eID uint32) (longRMInfo LongReadMappingInfo) {
	var path, path1 []constructdbg.DBG_MAX_INT
	cb1 := mpi.Arr[0]
	cb2 := mpi.Arr[len(mpi.Arr)-1]
	if len(bestEPIB.Path) <= 1 && len(bestEPIF.Path) <= 1 {
		path = append(path, constructdbg.DBG_MAX_INT(eID))
		path1 = append(path1, constructdbg.DBG_MAX_INT(0))
		longRMInfo.RStart, longRMInfo.REnd = cb1.GetRPos(), cb2.GetRPos()+uint32(cb2.Len)
		if cb1.GetStrand() == constructdbg.PLUS {
			longRMInfo.EStart, longRMInfo.EEnd = cb1.GetEPos(), cb2.GetEPos()+uint32(cb2.Len)
		} else {
			longRMInfo.EStart, longRMInfo.EEnd = cb1.GetEPos()+uint32(cb1.Len), cb2.GetEPos()
		}
	} else {
		if len(bestEPIB.Path) > 1 {
			lastB := bestEPIB.MpArr[len(bestEPIB.MpArr)-1]
			longRMInfo.RStart = lastB.GetRPos()
			if lastB.GetStrand() == constructdbg.PLUS {
				longRMInfo.EStart = lastB.GetEPos()
			} else {
				longRMInfo.EStart = lastB.GetEPos() + uint32(lastB.Len)
			}
			path = append(path, constructdbg.GetReverseDBG_MAX_INTArr(bestEPIB.Path)...)
			if len(bestEPIB.Path1) > 0 {
				path1 = append(path1, constructdbg.GetReverseDBG_MAX_INTArr(bestEPIB.Path1)...)
			} else {
				p := make([]constructdbg.DBG_MAX_INT, len(bestEPIB.Path))
				path1 = append(path1, p...)
			}
		} else {
			path = append(path, constructdbg.DBG_MAX_INT(eID))
			path1 = append(path1, constructdbg.DBG_MAX_INT(0))
			longRMInfo.RStart = cb1.GetRPos()
			if cb1.GetStrand() == constructdbg.PLUS {
				longRMInfo.EStart = cb1.GetEPos()
			} else {
				longRMInfo.EStart = cb1.GetEPos() + uint32(cb1.Len)
			}
		}

		if len(bestEPIF.Path) > 1 {
			lastF := bestEPIF.MpArr[len(bestEPIF.MpArr)-1]
			longRMInfo.REnd = lastF.GetRPos() + uint32(lastF.Len)
			if lastF.GetStrand() == constructdbg.PLUS {
				longRMInfo.EEnd = lastF.GetEPos() + uint32(lastF.Len)
			} else {
				longRMInfo.EEnd = lastF.GetEPos()
			}
			path = append(path, bestEPIF.Path[1:]...)
			if len(bestEPIF.Path1) > 0 {
				path1 = append(path1, bestEPIF.Path1[1:]...)
			} else {
				p := make([]constructdbg.DBG_MAX_INT, len(bestEPIF.Path))
				path1 = append(path1, p[1:]...)
			}
		} else {
			longRMInfo.REnd = cb2.GetRPos() + uint32(cb2.Len)
			if cb1.GetStrand() == constructdbg.PLUS {
				longRMInfo.EEnd = cb2.GetEPos() + uint32(cb2.Len)
			} else {
				longRMInfo.EEnd = cb2.GetEPos()
			}
		}
	}
	longRMInfo.Path[0] = path
	longRMInfo.Path[1] = path1
	return
}

//"denote max score before Pos"
type ScorePos struct {
	Pos   int
	Score int
}

func BiggerThanMaxScoreArr(maxScoreArr []ScorePos, sa ScorePos) bool {
	ok := true
	for i := len(maxScoreArr) - 1; i >= 0; i-- {
		t := maxScoreArr[i]
		if sa.Pos >= t.Pos {
			if sa.Pos > t.Pos {
				if sa.Score <= t.Score {
					ok = false
				} /*else if i < len(maxScoreArr) - 1 && sa.Score > 100 {
					next := maxScoreArr[i+1]
					fraction := float64(next.Score)/float64(next.Pos)
					if float64(sa.Score)/float64(sa.Pos) < fraction * 0.95 {
						ok = false
					}
				}*/
			} else if sa.Pos == t.Pos {
				if sa.Score < t.Score {
					ok = false
				} /* else if i < len(maxScoreArr) - 1 && sa.Score > 100 {
					next := maxScoreArr[i+1]
					fraction := float64(next.Score)/float64(next.Pos)
					if float64(sa.Score)/float64(sa.Pos) < fraction * 0.95 {
						ok = false
					}
				} */
			}
			break
		}
	}
	return ok

}

func AddedToMaxScoreArr(maxScoreArr []ScorePos, sp ScorePos) []ScorePos {
	added := false
	max := math.MinInt32
	for i, t := range maxScoreArr {
		if t.Score > max {
			max = t.Score
		}
		if t.Pos >= sp.Pos {
			if t.Score < sp.Score {
				added = true
				maxScoreArr[i] = sp
				j := i + 1
				for j < len(maxScoreArr) {
					nt := maxScoreArr[j]
					if nt.Score <= sp.Score {
						j++
						continue
					} else {
						break
					}
				}
				if j > i+1 {
					var a []ScorePos
					a = append(a, maxScoreArr[:i+1]...)
					a = append(a, maxScoreArr[j:]...)
					maxScoreArr = a
				}
			} else if t.Score > sp.Score {
				if i > 0 {
					if maxScoreArr[i-1].Score < sp.Score {
						added = true
					}
				} else {
					added = true
				}
				if added {
					var a []ScorePos
					a = append(a, maxScoreArr[:i]...)
					a = append(a, sp)
					a = append(a, maxScoreArr[i:]...)
					maxScoreArr = a
				}
			}
			break
		}
	}
	if added == false && sp.Score > max {
		maxScoreArr = append(maxScoreArr, sp)
		//added = true
	}
	return maxScoreArr
}

func GetUniqueMaxPathInfo(maxPathInfoArr []MaxPathInfo, kmerlen int) (uniqueMaxPathInfo MaxPathInfo) {
	for i, mpi := range maxPathInfoArr {
		a := mpi.Arr[0]
		b := mpi.Arr[len(mpi.Arr)-1]
		if i > 0 {
			mki := maxPathInfoArr[i-1].Arr[len(maxPathInfoArr[i-1].Arr)-1]
			last := int(mki.GetRPos()) + int(mki.Len)
			if last-int(a.GetRPos()) > kmerlen {
				continue
			}
		}
		if i < len(maxPathInfoArr)-1 {
			mki := maxPathInfoArr[i+1].Arr[0]
			first := int(mki.GetRPos())
			if int(b.GetRPos())+int(b.Len)-first > kmerlen {
				continue
			}
		}
		fmt.Printf("[GetUniqueMaxPathInfo] unique [%v] EID: %v, score: %v\n", i, mpi.Arr[0].EID, mpi.Score)
		if mpi.Score > uniqueMaxPathInfo.Score {
			uniqueMaxPathInfo = mpi
		}
	}
	return
}

func PrintAddedMpArr(mpArr []MapingKmerInfo) {
	for i, mki := range mpArr {
		fmt.Printf("[PrintAddedMpArr] [%v] EID: %v, strand: %v, EPos: %v, RPos: %v, Len: %v\n", i, mki.EID, mki.GetStrand(), mki.GetEPos(), mki.GetRPos(), mki.Len)
	}
}

func GetBubbleEdgeIDArr(path []constructdbg.DBG_MAX_INT, edgesArr []constructdbg.DBGEdge, nodesArr []constructdbg.DBGNode) (bubArr []constructdbg.DBG_MAX_INT, idxArr []int) {
	for i, id := range path {
		e := edgesArr[id]
		if e.GetBubbleFlag() > 0 {
			bubArr = append(bubArr, id)
			eArr := GetBubbleEIDArr(e, edgesArr, nodesArr)
			for _, eID := range eArr {
				if eID != id {
					bubArr = append(bubArr, eID)
				}
				idxArr = append(idxArr, i)
			}
		}
	}
	return
}

func IsBubbleRepeatInEPIPath(bubArr []constructdbg.DBG_MAX_INT) bool {
	for i, eID := range bubArr {
		for j := i + 1; j < len(bubArr); j++ {
			if eID == bubArr[j] {
				return true
			}
		}
	}
	return false
}

func GetEdgeMappingStrand(eID constructdbg.DBG_MAX_INT, mkArr []MapingKmerInfo, path []constructdbg.DBG_MAX_INT, direction uint8, edgesArr []constructdbg.DBGEdge, nodesArr []constructdbg.DBGNode) (strand bool) {
	strand = mkArr[0].GetStrand()
	e := edgesArr[path[0]]
	ok := false
	if eID == path[0] {
		return strand
	}
	for i := 1; i < len(path); i++ {
		ne := edgesArr[path[i]]
		strand = GetNextEdgeStrand(e, ne, nodesArr, strand)
		if ne.ID == eID {
			ok = true
			break
		}
		e = ne
	}
	if ok == false {
		log.Fatalf("[GetEdgeMappingStrand] not found edgeID: %v in the path: %v\n", eID, path)
	}
	return
}

func BubbleAlignment(epi ExtendPathInfo, edgesArr []constructdbg.DBGEdge, nodesArr []constructdbg.DBGNode, rd constructcf.ReadInfo, SeedLen, kmerlen int, direction uint8) ExtendPathInfo {
	bubArr, idxArr := GetBubbleEdgeIDArr(epi.Path, edgesArr, nodesArr)
	if IsBubbleRepeatInEPIPath(bubArr) {
		fmt.Printf("[BubbleAlignment] this Path: %v encounter one bubble display >2 times\n", epi.Path)
		return epi
	}
	if len(bubArr) == 0 {
		return epi
	}
	fmt.Printf("[BubbleAlignment] this Path: %v has bubble: %v, idxArr: %v\n", epi.Path, bubArr, idxArr)
	path1 := make([]constructdbg.DBG_MAX_INT, len(epi.Path))
	mkia, mkib := epi.MpArr[0], epi.MpArr[len(epi.MpArr)-1]
	var readRegionSeq []byte
	if direction == constructdbg.FORWARD {
		start, end := int(mkia.GetRPos()), int(mkib.GetRPos())+int(mkib.Len)
		readRegionSeq = rd.Seq[start:end]
	} else {
		start, end := int(mkia.GetRPos())+int(mkia.Len), int(mkib.GetRPos())
		readRegionSeq = rd.Seq[end:start]
	}
	edgesSeqArr := make([][]byte, len(bubArr))
	var strand bool
	for i, eID := range bubArr {
		if eID == epi.Path[idxArr[i]] {
			strand = GetEdgeMappingStrand(eID, epi.MpArr, epi.Path, direction, edgesArr, nodesArr)
		} else {
			e := edgesArr[epi.Path[idxArr[i]]]
			ne := edgesArr[eID]
			if ne.StartNID == e.EndNID && ne.EndNID == e.StartNID {
				strand = !strand
			}
		}
		if strand == constructdbg.PLUS {
			edgesSeqArr[i] = append(edgesSeqArr[i], edgesArr[eID].Utg.Ks...)
		} else {
			edgesSeqArr[i] = append(edgesSeqArr[i], constructdbg.GetReverseCompByteArr(edgesArr[eID].Utg.Ks)...)
		}
	}
	//fmt.Printf("[BubbleAlignment] readRegionSeq: %v\n", readRegionSeq)
	//fmt.Printf("[BubbleAlignment] edgesSeqArr: %v\n", edgesSeqArr)
	cgArr := AlignmentBlocks(readRegionSeq, edgesSeqArr, uint32(SeedLen))
	for i := 0; i < len(bubArr); {
		//e1 := edgesArr[bubArr[i]]
		j := i + 1
		for ; j < len(bubArr); j++ {
			//e2 := edgesArr[bubArr[j]]
			if idxArr[j] == idxArr[i] { // test if is a same bubble edge
				continue
			} else {
				break
			}
		}
		maxScore := math.MinInt32
		var eIDArr []constructdbg.DBG_MAX_INT
		for x := i; x < j; x++ {
			cg := cgArr[x]
			score := int(cg.Mch) - int(cg.Mis) - int(cg.Ins) - int(cg.Del)
			if score > maxScore {
				maxScore = score
				eIDArr = eIDArr[:0]
				eIDArr = append(eIDArr, bubArr[x])
			} else if score == maxScore {
				eIDArr = append(eIDArr, bubArr[x])
			}
		}
		if len(eIDArr) == 1 {
			epi.Path[idxArr[i]] = eIDArr[0]
		} else {
			epi.Path[idxArr[i]] = eIDArr[0]
			path1[idxArr[i]] = eIDArr[1]
		}
		i = j
	}
	epi.Path1 = path1
	return epi
}

func FindONTReadsPath(edgesKmerSortArr []KmerInfo, cs <-chan ReadsBucketInfo, SeedLen, winSize, GapCost, MaxGapLen, kmerlen int, wc chan<- LongReadMappingInfo, edgesArr []constructdbg.DBGEdge, nodesArr []constructdbg.DBGNode, SeqType uint8) {
	extendLen := 2 * kmerlen
	var totalReadNum, notFoundSeedEdgeNum, mapOneEdgeNum, ExtendNum int
	for {
		rb := <-cs
		if len(rb.ReadsArr) == 0 {
			var tmp LongReadMappingInfo
			wc <- tmp
			break
		}
		totalReadNum += len(rb.ReadsArr)
		startReadID := rb.ReadsArr[0].ID
		readsBucketInterSecKmerInfoArr := GetInterSectionKmerInfo(edgesKmerSortArr, rb.KmerSortArr, len(rb.ReadsArr), uint32(startReadID), uint16(SeedLen))
		// loop process every read
		for j, ka := range readsBucketInterSecKmerInfoArr {
			rID := uint32(startReadID) + uint32(j)
			rl := len(rb.ReadsArr[j].Seq)
			sort.Sort(MapingKmerInfoArr(ka))
			kb := GetChainBlocks(ka)
			sort.Sort(MapingKmerInfoArr(kb))
			fmt.Printf("[FindONTReadsPath] read[%v] length: %v\n", rID, rl)
			//kbEdgeArr := GetEdgeChainBlocks(kb, MaxGapLen)
			/*for y, mki := range kb {
				fmt.Printf("[FindONTReadsPath] kb[%v]EID: %v, ReadPos : %v, EdgePos: %v,  Len: %v\n", y, mki.EID, mki.GetRPos(), mki.GetEPos(), mki.Len)
			}*/
			MinKmerScore := SeedLen * 2
			maxPathInfoArr := GetMaxChainBlockArr(kb, edgesArr, nodesArr, rl, SeedLen, GapCost, MaxGapLen, MinKmerScore, kmerlen, SeqType)

			sort.Sort(MaxPathInfoTmpArr(maxPathInfoArr))
			//var cBMRIArr []ChainBlockMapLongReadInfo
			for x, ele := range maxPathInfoArr {
				//if ele.Score > 5*SeedLen {
				start, end := ele.Arr[0].GetRPos(), ele.Arr[len(ele.Arr)-1].GetRPos()+uint32(ele.Arr[len(ele.Arr)-1].Len)
				fmt.Printf("[FindONTReadsPath] maxPathInfoArr[%v]: Edge ID: %v, Read region[%v---%v] Len:%v, score: %v, Percent: %v\n", x, ele.Arr[0].EID, start, end, end-start, ele.Score, ele.Score*100/(int(end)-int(start)))
				/*for y, mki := range ele.Arr {
					fmt.Printf("[FindONTReadsPath]  [%v] EdgePos: %v, ReadPos : %v, Len: %v\n", y, mki.GetEPos(), mki.GetRPos(), mki.Len)
				}*/
				//}
			}
			uniqueMaxPathInfo := GetUniqueMaxPathInfo(maxPathInfoArr, kmerlen)
			if uniqueMaxPathInfo.Score == 0 {
				notFoundSeedEdgeNum++
				fmt.Printf("[FindONTReadsPath] read ID[%v] not found unique boundary mapping region\n", rID)
				continue
			}
			start, end := int(uniqueMaxPathInfo.Arr[0].GetRPos()), int(uniqueMaxPathInfo.Arr[len(uniqueMaxPathInfo.Arr)-1].GetRPos())+int(uniqueMaxPathInfo.Arr[len(uniqueMaxPathInfo.Arr)-1].Len)
			fmt.Printf("[FindONTReadsPath] uniqueMaxPathInfo: Edge ID: %v, Read region[%v---%v] Len:%v, score: %v, Percent: %v\n", uniqueMaxPathInfo.Arr[0].EID, start, end, end-start, uniqueMaxPathInfo.Score, uniqueMaxPathInfo.Score*100/(int(end)-int(start)))
			if int(start) < 2*kmerlen && int(end) > rl-2*kmerlen {
				mapOneEdgeNum++
				fmt.Printf("[FindONTReadsPath] read ID[%v] only map one long edge[%v]\n", rID, uniqueMaxPathInfo.Arr[0].EID)
				continue
			}
			//for x := len(maxPathInfoArr) - 1; x >= 0; x-- {
			mpi := uniqueMaxPathInfo
			eID := mpi.Arr[0].EID
			edge := edgesArr[eID]
			if edge.StartNID == edge.EndNID || edge.GetTwoEdgesCycleFlag() > 0 || len(edge.Utg.Ks) < 3*kmerlen || end-start < 2*kmerlen {
				fmt.Printf("[FindONTReadsPath] read ID[%v] unique map edge[%v] not allow extend path\n", rID, uniqueMaxPathInfo.Arr[0].EID)
				continue
			}
			// check this region has been processed
			last := mpi.Arr[len(mpi.Arr)-1]
			//start, end := int(mpi.Arr[0].GetRPos()), int(last.GetRPos())+int(last.Len)
			fmt.Printf("[FindONTReadsPath] Read region[%v---%v], readLen: %v, Edge region[%v---%v], EdgeLen: %v, score: %v\n", mpi.Arr[0].GetRPos(), last.GetRPos(), rl, mpi.Arr[0].GetEPos(), last.GetEPos(), len(edgesArr[eID].Utg.Ks), mpi.Score)
			/*if mpi.Score < 5*SeedLen || IsInChainBlockMapLongReadInfoArr(cBMRIArr, start, end) {
				continue
			}*/
			// get base alignment for seed chain
			readRegionSeq := rb.ReadsArr[j].Seq[start:end]
			edgeRegionSeq := GetEdgeRegionSeq(edgesArr[eID], mpi)
			chainArr, _ := ComputingChainArr(edgeRegionSeq, readRegionSeq, uint32(SeedLen))
			fmt.Printf("[FindONTReadsPath] base alignment found chainArr length: %v\n", len(chainArr))
			//fmt.Printf(">edge\n%v\n", constructdbg.Transform2Letters(edgeRegionSeq))
			//fmt.Printf(">read\n%v\n", constructdbg.Transform2Letters(readRegionSeq))
			if len(chainArr) == 0 {
				continue
			}
			cg := DPLocalAlign(edgeRegionSeq, readRegionSeq, chainArr)
			if int(cg.Mch) < len(edgeRegionSeq)*82/100 {
				continue
			}
			//avgSeqIdentityPercent := mpi.Score * 100 / len(edgeRegionSeq)
			var bestEPIF, bestEPIB ExtendPathInfo
			fmt.Printf("[FindONTReadsPath] Seed EdgeID: %v, mpi.Score: %v, cg.Mch: %v, len(edgeRegionSeq): %v\n", eID, mpi.Score, cg.Mch, len(edgeRegionSeq))
			// read FORWARD mapping
			{
				// check if is boundary of long read
				pos := end
				if pos+extendLen < rl {
					anchor := last
					var mpArr []MapingKmerInfo
					mpArr = append(mpArr, anchor)
					var epi ExtendPathInfo
					id := eID
					epi.Path = append(epi.Path, constructdbg.DBG_MAX_INT(id))
					epi.Strand = anchor.GetStrand()
					var coming bool
					epi.Coming = GetEdgeComing(edgesArr[id], constructdbg.FORWARD, epi.Strand, nodesArr, coming)
					if anchor.GetStrand() == constructdbg.PLUS {
						epi.ExtendLen = len(edgesArr[id].Utg.Ks) - int(anchor.GetEPos())
						epi.FlankLen = len(edgesArr[id].Utg.Ks) - (int(anchor.GetEPos()) + int(anchor.Len))
					} else {
						epi.ExtendLen = int(anchor.GetEPos()) + int(anchor.Len)
						epi.FlankLen = int(anchor.GetEPos())
					}
					epi.MpArr = mpArr
					epi.Score = int(anchor.Len)
					//maxScore := epi.Score
					//scoreLen := int(anchor.Len)
					//maxMapLen := int(anchor.Len)
					var stk []ExtendPathInfo
					stk = append(stk, epi)
					var highQualEPathInfoArr []ExtendPathInfo
					var maxScoreArr []ScorePos
					for {
						// found min ExtendLen path
						var t ExtendPathInfo
						max := float64(0) // max score in the stack
						idx := -1
						//var ele ExtendPathInfo
						for x, ele := range stk {
							if ele.Score <= 0 {
								continue
							}
							mapLen := ele.ExtendLen - ele.FlankLen
							fraction := float64(ele.Score) / float64(mapLen)
							if fraction > max {
								max = fraction
								t = ele
								idx = x
							}
						}
						if idx >= 0 {
							stk[idx].Score = 0
						}
						if t.Score <= 0 {
							break
						}
						var sa ScorePos
						sa.Score, sa.Pos = t.Score, t.ExtendLen-t.FlankLen
						if BiggerThanMaxScoreArr(maxScoreArr, sa) == false {
							continue
						}
						// add next near edge info
						id := t.Path[len(t.Path)-1]
						e := edgesArr[id]
						eArr := GetNearDirectionEdgeIDArr(e, constructdbg.FORWARD, t.Strand, t.Coming, nodesArr)
						fmt.Printf("[FindONTReadsPath]idx: %v, next eArr: %v, len(t.MpArr): %v, Score: %v, ExtendLen: %v, FlankLen: %v, Strand: %v, Coming: %v, Path: %v\n", idx, eArr, len(t.MpArr), t.Score, t.ExtendLen, t.FlankLen, t.Strand, t.Coming, t.Path)
						for _, eID := range eArr {
							ne := edgesArr[eID]
							var nt ExtendPathInfo
							strand := GetNextEdgeStrand(e, ne, nodesArr, t.Strand)
							nt.Path = append(nt.Path, t.Path...)
							nt.Path = append(nt.Path, eID)
							nt.Strand = strand
							// extend chains
							nt.MpArr, nt.Score = ExtendPathChain(t.MpArr, t.FlankLen, kmerlen, constructdbg.FORWARD, ne, kb, strand, t.Score, edgesArr, uint16(SeedLen))
							//fmt.Printf("[FindONTReadsPath]eID: %v, nt.Score: %v\n", eID, nt.Score)
							nt.Coming = GetEdgeComing(edgesArr[eID], constructdbg.FORWARD, strand, nodesArr, t.Coming)
							nt.ExtendLen += t.ExtendLen + (len(ne.Utg.Ks) - (kmerlen - 1))
							if len(nt.MpArr) > len(t.MpArr) {
								nt.FlankLen = GetMappingDBGFlankLen(nt, constructdbg.FORWARD, len(ne.Utg.Ks))
							} else {
								nt.FlankLen = t.FlankLen + len(ne.Utg.Ks) - (kmerlen - 1)
							}
							fmt.Printf("[FindONTReadsPath]eID: %v, len(nt.MpArr): %v, Score: %v, ExtendLen: %v, FlankLen: %v, Strand: %v, Coming: %v, Path: %v\n", eID, len(nt.MpArr), nt.Score, nt.ExtendLen, nt.FlankLen, nt.Strand, nt.Coming, nt.Path)
							// delete low score extend paths
							var sp ScorePos
							sp.Pos = nt.ExtendLen - nt.FlankLen
							sp.Score = nt.Score
							if BiggerThanMaxScoreArr(maxScoreArr, sp) {
								if nt.FlankLen > 500 {
									highQualEPathInfoArr = append(highQualEPathInfoArr, nt)
								} else {
									p := int(nt.MpArr[len(nt.MpArr)-1].GetRPos()) + int(nt.MpArr[len(nt.MpArr)-1].Len)
									if p > rl-kmerlen {
										highQualEPathInfoArr = append(highQualEPathInfoArr, nt)
									} else {
										stk = append(stk, nt)
									}
								}
								maxScoreArr = AddedToMaxScoreArr(maxScoreArr, sp)
								fmt.Printf("[FindONTReadsPath]sp: %v, maxScoreArr: %v\n", sp, maxScoreArr)
							}
						}
						fmt.Printf("[FindONTReadsPath] len(stk): %v\n", len(stk))
						// process bubble edges, just choose best score edge
						if len(eArr) == 2 && len(stk) > 2 && edgesArr[eArr[0]].GetBubbleFlag() > 0 {
							MuskBubble(eArr, stk)
						}
						//stk[idx].Score = 0
						//constructdbg.GetNearEdgeIDArr(nd, eID, coming)
						//constructdbg.GetNeighbourEID(eID, nID, edgesArr, nodesArr, max_insert_size)
						//constructdbg.GetNextDirection(eID, edge, nodesArr)
						//constructdbg.GetNextEID(eID, node)
						//constructdbg.GetNextMappingEID(nd, e, base)
					}

					// process highQualEPathInfoArr
					bestEPIF = GetBestEPathInfo(highQualEPathInfoArr, edgesArr, nodesArr, rb.ReadsArr[j], SeedLen, kmerlen, constructdbg.FORWARD)
					fmt.Printf("[FindONTReadsPath] bestEPIF path: %v\n", bestEPIF.Path)
					if len(bestEPIF.Path) > 1 {
						bestEPIF = BubbleAlignment(bestEPIF, edgesArr, nodesArr, rb.ReadsArr[j], SeedLen, kmerlen, constructdbg.FORWARD)
					}
				}
			}
			// read BACKWARD mapping
			{
				// check if is boundary of long read
				pos := start
				if pos-extendLen > 0 {
					anchor := mpi.Arr[0]
					var mpArr []MapingKmerInfo
					mpArr = append(mpArr, anchor)
					var epi ExtendPathInfo
					id := eID
					epi.Path = append(epi.Path, constructdbg.DBG_MAX_INT(id))
					epi.Strand = anchor.GetStrand()
					var coming bool
					epi.Coming = GetEdgeComing(edgesArr[id], constructdbg.BACKWARD, epi.Strand, nodesArr, coming)
					if anchor.GetStrand() == constructdbg.PLUS {
						epi.ExtendLen = int(anchor.GetEPos()) + int(anchor.Len)
						epi.FlankLen = int(anchor.GetEPos())
					} else {
						epi.ExtendLen = len(edgesArr[id].Utg.Ks) - int(anchor.GetEPos())
						epi.FlankLen = epi.ExtendLen - int(anchor.Len)
					}
					epi.MpArr = mpArr
					epi.Score = int(anchor.Len)
					//maxMapLen := int(anchor.Len)
					var stk []ExtendPathInfo
					stk = append(stk, epi)
					var highQualEPathInfoArr []ExtendPathInfo
					var maxScoreArr []ScorePos
					for {
						// found min ExtendLen path
						var t ExtendPathInfo
						max := float64(0) // max score in the stack
						idx := -1
						for x, ele := range stk {
							if ele.Score <= 0 {
								continue
							}
							mapLen := ele.ExtendLen - ele.FlankLen
							fraction := float64(ele.Score) / float64(mapLen)
							if fraction > max {
								max = fraction
								t = ele
								idx = x
							}
						}
						if idx >= 0 {
							stk[idx].Score = 0
						}
						if t.Score <= 0 {
							break
						}
						var sa ScorePos
						sa.Score, sa.Pos = t.Score, t.ExtendLen-t.FlankLen
						if BiggerThanMaxScoreArr(maxScoreArr, sa) == false {
							continue
						}
						// add next near edge info
						id := t.Path[len(t.Path)-1]
						e := edgesArr[id]
						eArr := GetNearDirectionEdgeIDArr(e, constructdbg.BACKWARD, t.Strand, t.Coming, nodesArr)
						fmt.Printf("[FindONTReadsPath]idx: %v, next eArr: %v, len(t.MpArr): %v, Score: %v, ExtendLen: %v, FlankLen: %v, Strand: %v, Coming: %v, Path: %v\n", idx, eArr, len(t.MpArr), t.Score, t.ExtendLen, t.FlankLen, t.Strand, t.Coming, t.Path)
						for _, eID := range eArr {
							ne := edgesArr[eID]
							var nt ExtendPathInfo
							strand := GetNextEdgeStrand(e, ne, nodesArr, t.Strand)
							nt.Path = append(nt.Path, t.Path...)
							nt.Path = append(nt.Path, eID)
							nt.Strand = strand
							nt.MpArr, nt.Score = ExtendPathChain(t.MpArr, t.FlankLen, kmerlen, constructdbg.BACKWARD, ne, kb, strand, t.Score, edgesArr, uint16(SeedLen))
							//PrintAddedMpArr(nt.MpArr[len(t.MpArr):])
							nt.Coming = GetEdgeComing(edgesArr[eID], constructdbg.BACKWARD, strand, nodesArr, t.Coming)
							nt.ExtendLen += t.ExtendLen + (len(ne.Utg.Ks) - (kmerlen - 1))
							if len(nt.MpArr) > len(t.MpArr) {
								nt.FlankLen = GetMappingDBGFlankLen(nt, constructdbg.BACKWARD, len(ne.Utg.Ks))
							} else {
								nt.FlankLen = t.FlankLen + len(ne.Utg.Ks) - (kmerlen - 1)
							}
							fmt.Printf("[FindONTReadsPath]eID: %v, len(nt.MpArr): %v, Score: %v, ExtendLen: %v, FlankLen: %v, Strand: %v, Coming: %v, Path: %v\n", eID, len(nt.MpArr), nt.Score, nt.ExtendLen, nt.FlankLen, nt.Strand, nt.Coming, nt.Path)
							// delete low score extend paths
							//if nt.Score*100/(nt.ExtendLen-nt.FlankLen) < avgSeqIdentityPercent*9/10 {
							var sp ScorePos
							sp.Pos = nt.ExtendLen - nt.FlankLen
							sp.Score = nt.Score
							//fmt.Printf("[FindONTReadsPath]sp: %v, maxScoreArr: %v\n", sp, maxScoreArr)
							if BiggerThanMaxScoreArr(maxScoreArr, sp) {
								if nt.FlankLen > 500 {
									highQualEPathInfoArr = append(highQualEPathInfoArr, nt)
								} else {
									p := int(nt.MpArr[len(nt.MpArr)-1].GetRPos())
									if p < kmerlen {
										highQualEPathInfoArr = append(highQualEPathInfoArr, nt)
									} else {
										stk = append(stk, nt)
									}
								}
								maxScoreArr = AddedToMaxScoreArr(maxScoreArr, sp)
								fmt.Printf("[FindONTReadsPath]sp: %v, maxScoreArr: %v\n", sp, maxScoreArr)
							}

						}
						fmt.Printf("[FindONTReadsPath] len(stk): %v\n", len(stk))
						// process bubble edges, just choose best score edge
						if len(eArr) == 2 && len(stk) > 2 && edgesArr[eArr[0]].GetBubbleFlag() > 0 {
							MuskBubble(eArr, stk)
						}
						//stk[idx].Score = 0
					}

					// process highQualEPathInfoArr
					bestEPIB = GetBestEPathInfo(highQualEPathInfoArr, edgesArr, nodesArr, rb.ReadsArr[j], SeedLen, kmerlen, constructdbg.BACKWARD)
					fmt.Printf("[FindONTReadsPath] bestEPIB path: %v\n", bestEPIB.Path)
					if len(bestEPIB.Path) > 1 {
						bestEPIB = BubbleAlignment(bestEPIB, edgesArr, nodesArr, rb.ReadsArr[j], SeedLen, kmerlen, constructdbg.BACKWARD)
					}
				}
			}

			// extend path by DBG
			var longRMInfo LongReadMappingInfo
			bl, fl := len(bestEPIB.Path), len(bestEPIF.Path)
			bestEPIB = DeleteUnmapPath(bestEPIB, constructdbg.BACKWARD, edgesArr, kmerlen)
			bestEPIF = DeleteUnmapPath(bestEPIF, constructdbg.FORWARD, edgesArr, kmerlen)
			if len(bestEPIB.Path) < bl || len(bestEPIF.Path) < fl {
				fmt.Printf("[FindONTReadsPath]len(bestEPIB.Path): %v < bl: %v, len(bestEPIF.Path):%v < fl: %v\n", len(bestEPIB.Path), bl, len(bestEPIF.Path), fl)
			}
			fmt.Printf("[FindONTReadsPath] bestEPIB: %v\n", bestEPIB)
			fmt.Printf("[FindONTReadsPath] bestEPIF: %v\n", bestEPIF)
			longRMInfo = MergeTwoFlankMapingPath(bestEPIB, bestEPIF, mpi, eID)
			longRMInfo.ID = rID
			longRMInfo.Score = bestEPIB.Score + mpi.Score + bestEPIF.Score
			longRMInfo.AnchorEdgeID = eID
			longRMInfo.AnchorStrand = last.GetStrand()
			longRMInfo.Qual = uint32(longRMInfo.Score) * 100 / (longRMInfo.REnd - longRMInfo.RStart)
			/*// add read region to the cBMRIArr
			var cbmr ChainBlockMapLongReadInfo
			cbmr.Start, cbmr.End, cbmr.Score = int(longRMInfo.RStart), int(longRMInfo.REnd), int(longRMInfo.Score)
			cBMRIArr = append(cBMRIArr, cbmr)*/
			ok := true
			if bestEPIB.Score == 0 {
				if start > extendLen {
					ok = false
				}
			}
			if bestEPIF.Score == 0 {
				if end < rl-extendLen {
					ok = false
				}
			}
			if ok {
				ExtendNum++
			}
			if len(longRMInfo.Path) > 0 {
				wc <- longRMInfo
				fmt.Printf("[FindONTReadsPath] longRMInfo: %v\n", longRMInfo)
			}

			/*//kbArr := CheckByDBG(maxKbArr, nodesArr, edgesArr)
			kbArr := maxKbArr
			scArr := make([]int, len(kbArr))
			idx := 0
			for m, kb := range kbArr {
				scArr[m] = GetMappingKIArrScore(kb.Arr, MaxGapLen, GapCost)
				if scArr[m] < rl/15 {
					fmt.Printf("[FindONTReadsPath] chainBlocksScore: %v < readLen/15: %v\n", scArr[m], rl/15)
				} else {
					kbArr[idx] = kbArr[m]
					idx++
				}
			}
			kbArr = kbArr[:idx]
			cgArr := make([]CIGAR, len(kbArr))
			pathArr := make([][]constructdbg.DBG_MAX_INT, len(kbArr))
			for m, mp := range kbArr {
				kb := mp.Arr
				var edgesSeq []byte
				var chainArr []Chain
				edgesSeq, pathArr[m], chainArr = GetEdgesSeq(edgesArr, nodesArr, kb, kmerlen)
				readSeq := GetReadSeq(rb.ReadsArr[j].Seq, kb[0].GetRPos(), kb[len(kb)-1].GetRPos()+uint32(kb[len(kb)-1].Len))
				cgArr[m] = DPLocalAlign(edgesSeq, readSeq, chainArr)
			}
			maxScore := math.MinInt32
			var count int
			maxArr := make([][]constructdbg.DBG_MAX_INT, 0, len(cgArr))
			for m, cg := range cgArr {
				sc := int(cg.Mch) - int(cg.Mis) - int(cg.Del) - int(cg.Ins)
				if sc > maxScore {
					maxScore = sc
					maxArr[0] = pathArr[m]
					maxArr = maxArr[:1]
				} else if sc == maxScore {
					maxArr = append(maxArr, pathArr[m])
				}
			}

			// found trustful path
			var longRMInfo LongReadMappingInfo
			longRMInfo.ID = rID
			for n := 0; n < len(maxArr[0]); n++ {
				var eID1, eID2 constructdbg.DBG_MAX_INT
				stop := false
				for m, path := range maxArr {
					if len(path) <= n {
						stop = true
						break
					}
					if eID1 < 2 {
						eID1 = path[n]
					} else if path[n] != eID1 {
						if eID2 < 2 {
							eID2 = path[n]
						} else if path[n] != eID2 {
							stop = true
							break
						}
					}
				}
				if stop {
					break
				}
				if eID1 > 1 {
					if eID2 > 1 {
						if IsBubble(edgesArr[eID1], edgesArr[eID2], nodesArr) {
							longRMInfo.Path[0] = append(longRMInfo.Path[0], eID1)
							longRMInfo.Path[1] = append(longRMInfo.Path[1], eID2)
						}
					} else {
						longRMInfo.Path[0] = append(longRMInfo.Path[0], eID1)
					}
				} else {
					break
				}
			}
			if len(longRMInfo.Path[0]) > 2 {
				wc <- longRMInfo
			}*/

		}
	}
	fmt.Printf("[FindONTReadsPath] total reads num: %v, not found seed edge num: %v, map one edge num: %v, extend read num: %v\n", totalReadNum, notFoundSeedEdgeNum, mapOneEdgeNum, ExtendNum)
	return
}

func DeconstructDBG(c cli.Command) {
	// check arguments
	gOpt, suc := utils.CheckGlobalArgs(c.Parent())
	if suc == false {
		log.Fatalf("[DeconstructDBG] check global Arguments error, opt: %v\n", gOpt)
	}

	opt := Options{gOpt, 0, 0, 0, 0, 0, 0, 0, 0, "", false}
	tmp, suc := checkArgs(c)
	if suc == false {
		log.Fatalf("[DeconstructDBG] check Arguments error, opt: %v\n", tmp)
	}
	//opt.TipMaxLen = tmp.TipMaxLen
	opt.MinDepth = tmp.MinDepth
	opt.AvgDepth = tmp.AvgDepth
	opt.MinUniqueEdgeLen = tmp.MinUniqueEdgeLen
	opt.AvgReadLen = tmp.AvgReadLen
	opt.WinSize = tmp.WinSize
	opt.MaxNGSReadLen = tmp.MaxNGSReadLen
	opt.MinMapFreq = tmp.MinMapFreq
	//opt.MaxMapEdgeLen = tmp.MaxMapEdgeLen
	opt.ExtLen = tmp.ExtLen
	opt.ONTFn = tmp.ONTFn
	opt.Correct = tmp.Correct
	GapCost := 1
	MaxGapLen := 2000
	SeqType := uint8(1)
	//constructdbg.Kmerlen = opt.Kmer
	fmt.Printf("Arguments: %v\n", c.Flags())

	/*profileFn := opt.Prefix + ".decdbg.prof"
	cpuprofilefp, err := os.Create(profileFn)
	if err != nil {
		log.Fatalf("[DeconstructDBG] open cpuprofile file: %v failed\n", profileFn)
	}
	pprof.StartCPUProfile(cpuprofilefp)
	defer pprof.StopCPUProfile()*/

	// read files and construt DBG
	DBGInfofn := opt.Prefix + ".smfy.DBGInfo"
	eSize, nSize := constructdbg.DBGInfoReader(DBGInfofn)
	nodesfn := opt.Prefix + ".nodes.smfy.Arr.br"
	nodesArr := make([]constructdbg.DBGNode, nSize)
	constructdbg.NodesArrReader(nodesfn, nodesArr, opt.Kmer)
	if len(nodesArr) != nSize {
		log.Fatalf("[DeconstructDBG] len(nodesArr): %v != nodesArr Size: %v in file: %v\n", len(nodesArr), nSize, DBGInfofn)
	}
	edgesArr := make([]constructdbg.DBGEdge, eSize)
	edgesfn := opt.Prefix + ".edges.smfy.fq.br"
	constructdbg.LoadEdgesfqFromFn(edgesfn, edgesArr, true)

	constructdbg.CheckInterConnectivity(edgesArr, nodesArr)

	uniqueNum, semiUniqueNum, twoCycleNum, selfCycle, selfCycleSameComingNum, bubbleEdgeNum := constructdbg.SetDBGEdgesUniqueFlag(edgesArr, nodesArr, opt.MinUniqueEdgeLen)
	fmt.Printf("[DeconstructDBG] unique edge number is : %v, semiUniqueNum: %v, twoCycleNum:%v, selfCycle:%v, selfCycleSameComingNum: %v, bubbleEdgeNum: %v\n", uniqueNum, semiUniqueNum, twoCycleNum, selfCycle, selfCycleSameComingNum, bubbleEdgeNum)
	t0 := time.Now()
	pathFn := opt.Prefix + ".Path.protobuf"
	// no any other align tools for seed edge by ONT reads, use same as daligner,
	// first sort DBG edges by seed kmer, and split ONT reads by bucket, every bucket sort by seed kmer, and found share kmers from edges sort kmersArr
	if _, err := os.Stat(pathFn); err != nil {
		SeedLen := 15
		//BucketSize := (1 << 28) // ==2**28
		cfgInfo, err := constructcf.ParseCfg(opt.CfgFn, false)
		if err != nil {
			log.Fatalf("[DeconstructDBG] ParseCfg 'C': %v err :%v\n", opt.CfgFn, err)
		}
		fmt.Printf("[DeconstructDBG] opt: %v\n\tcfgInfo: %v\n", opt, cfgInfo)

		edgesKmerSortArr := SortDBGEdgesKmer(edgesArr, SeedLen, opt.WinSize)
		fmt.Printf("[DeconstructDBG]total found edges minimizers number is: %v\n", len(edgesKmerSortArr))
		bufSize := int64(len(edgesKmerSortArr)) * int64(opt.WinSize)
		var processedCPUNum int
		if opt.NumCPU <= 2 {
			processedCPUNum = 1
		} else {
			processedCPUNum = opt.NumCPU * 3 / 4
		}
		cs := make(chan ReadsBucketInfo, opt.NumCPU)
		wc := make(chan LongReadMappingInfo, processedCPUNum*100)
		go LoadLongReadsAndSort(cfgInfo, opt.NumCPU-processedCPUNum, cs, SeedLen, opt.WinSize, bufSize)
		for i := 0; i < processedCPUNum; i++ {
			go FindONTReadsPath(edgesKmerSortArr, cs, SeedLen, opt.WinSize, GapCost, MaxGapLen, opt.Kmer, wc, edgesArr, nodesArr, SeqType)
		}
		//func FindONTReadsPath(edgesKmerSortArr []KmerInfo, cs <-chan ReadsBucketInfo, SeedLen, winSize, GapCost, MaxGapLen, kmerlen int, wc chan<- LongReadMappingInfo, edgesArr []constructdbg.DBGEdge, nodesArr []constructdbg.DBGNode, SeqType uint8) {

		//pathFn := opt.Prefix + ".Path.protobuf"
		WriteLongPathToFile(wc, pathFn, processedCPUNum)
		fmt.Printf("[DeconstructDBG] Maping long reads to DBG used: %v\n", time.Now().Sub(t0))
	}
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
	t0 = time.Now()
	readPathArr := LoadLongPathFromFile(pathFn)
	// Simplify using Long Reads Mapping info
	joinPathArr := SimplifyByLongReadsPath(edgesArr, nodesArr, readPathArr, opt, opt.MinUniqueEdgeLen, 2000, opt.AvgDepth, opt.AvgReadLen)

	//graphfn := opt.Prefix + ".afterLR.dot"
	//constructdbg.GraphvizDBGArr(nodesArr, edgesArr, graphfn)
	DcDBGEdgesfn := opt.Prefix + ".edges.DcDBG.fq"
	ExtractSeq(edgesArr, nodesArr, joinPathArr, DcDBGEdgesfn, opt.Kmer)
	fmt.Printf("[DeconstructDBG] find max edge path and produce sequence used: %v\n", time.Now().Sub(t0))
	//constructdbg.StoreEdgesToFn(DcDBGEdgesfn, edgesArr)
}
