package main

import (
	"bufio"
	"fmt"
	"io"
	"log"
	"math"
	"os"
	"reflect"
	"runtime"
	"runtime/pprof"
	"sort"
	"strconv"
	"time"

	"github.com/google/brotli/go/cbrotli"
	"github.com/jwaldrip/odin/cli"
)

//"Debug const set Debug model if true"
var DebugModel bool = false
var BuckReadsNum int = 100
var EdgeMapKmerMin int = 60 // minimum kmer number of edge used map Long Reads
var MinEdgeLen int = 240    // (kmerlen-1)*2
var MinMinimerNum int = 5

const SeedLen = 21 // SeedLen allow <=23
const MaxAllowSeqLen = (1 << seqPositionBitNum) - 1

type ReadPath struct {
	ReadID       uint32
	Path0, Path1 []uint32
}

func addPathToPathMat(edgesArr []DBGEdge, eID uint32, path Path) bool {
	for i := 1; i < len(edgesArr[eID].PathMat); i++ {
		p := edgesArr[eID].PathMat[i]
		if reflect.DeepEqual(p.IDArr, path.IDArr) {
			edgesArr[eID].PathMat[i].Freq += path.Freq
			return true
		} else if reflect.DeepEqual(p.IDArr, GetReverseUint32Arr(path.IDArr)) {
			edgesArr[eID].PathMat[i].Freq += path.Freq
			return true
		}
	}
	edgesArr[eID].PathMat = append(edgesArr[eID].PathMat, path)
	return false
}

func WriteLongPathToFile(wc chan string, pathfn, ONTMapInfofn string, numCPU int) {
	pathfp, err := os.Create(pathfn)
	if err != nil {
		log.Fatalf("[WriteLongPathToFile] file %s create error, err: %v\n", pathfn, err)
	}
	defer pathfp.Close()
	cbrofp := cbrotli.NewWriter(pathfp, cbrotli.WriteroptionsDDBG{Quality: 1})
	defer cbrofp.Close()
	buffp := bufio.NewWriterSize(cbrofp, 1<<14) // 1<<25 == 2**25
	var finishNum int
	mapRecordNum := 0
	//var pathArr PathArr
	//pathArr.Arr = make([]*ReadPath, 0, 10000)
	//var rpArr []ReadPath
	for {
		ms := <-wc
		//extArr := rmi.Path
		//fmt.Fprintf(os.Stderr, "[WriteLongPathToFile] ms: %v\n", ms)
		//p := extArr[0]
		if len(ms) == 0 {
			finishNum++
			if finishNum == numCPU {
				break
			}
			continue
		}

		fmt.Fprintf(buffp, "%s\n", ms)
		//buffp.WriteString(ms)
		mapRecordNum++
		//pathArr = append(pathArr, extArr)
		// write path to the DBG
		// all long reads path store this pathArr, and store index of pathArr to the edge PathMat[1]
		/*if len(p) > 2 {
			var path Path
			path.IDArr = p
			path.Freq = 1
			//pathArr = append(pathArr, p)
			//idx := uint32(len(pathArr) - 1)
			for _, eID := range p {
				if edgesArr[eID].GetUniqueFlag() > 0 {
					if len(edgesArr[eID].PathMat) == 0 {
						var t Path
						edgesArr[eID].PathMat = append(edgesArr[eID].PathMat, t)
					}
					addPathToPathMat(edgesArr, eID, path)
					//edgesArr[eID].PathMat[1].IDArr = append(edgesArr[eID].PathMat[1].IDArr, idx)
				}
			}
		}*/
	}
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

	mapfp, err := os.Create(ONTMapInfofn)
	if err != nil {
		log.Fatalf("[WriteLongPathToFile] file %s create error, err: %v\n", ONTMapInfofn, err)
	}
	defer mapfp.Close()
	fmt.Fprintf(mapfp, "MapRecordSize:%d\n", mapRecordNum)
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

// set the unique edge of edgesArr
/*func setDBGEdgesUniqueFlag(edgesArr []DBGEdge, nodesArr []DBGNode) (uniqueNum int) {
	for i, e := range edgesArr {
		unique := true
		if e.ID == 0 || e.GetDeleteFlag() > 0 {
			continue
		}

		if e.StartNID > 0 {
			nd := nodesArr[e.StartNID]
			ea := GetNearEdgeIDArr(nd, e.ID)
			if len(ea) > 1 {
				unique = false
			}
		}

		if unique == true && e.EndNID > 0 {
			nd := nodesArr[e.EndNID]
			ea := GetNearEdgeIDArr(nd, e.ID)
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

type EIDLen struct {
	EID     uint32
	NextEID uint32
	Length  int
}

type EIDLenArr []EIDLen

func (a EIDLenArr) Len() int {
	return len(a)
}

func (a EIDLenArr) Less(i, j int) bool {
	return a[i].Length > a[j].Length
}

func (a EIDLenArr) Swap(i, j int) {
	a[i], a[j] = a[j], a[i]
}

func GetDirectionBubbleRepeatArr(e *DBGEdge, edgesArr []DBGEdge, direction uint8) (bubbleRepeatArr []uint32) {
	strand := PLUS
	le := e
	var ea []uint32
	for {
		if direction == FORWARD {
			if strand {
				ea = GetComingEdgeArr(le.EdgeIDOutcoming)
			} else {
				ea = GetComingEdgeArr(le.EdgeIDIncoming)
			}
		} else {
			if strand {
				ea = GetComingEdgeArr(le.EdgeIDIncoming)
			} else {
				ea = GetComingEdgeArr(le.EdgeIDOutcoming)
			}
		}
		if len(ea) == 0 {
			break
		}
		ne := &edgesArr[ea[0]]
		strand = GetNextEdgeStrand2(le, ne, strand, direction)
		le = ne
		if ne.GetBubbleFlag() > 0 {
			continue
		} else if ne.GetBubbleRepeatFlag() > 0 {
			bubbleRepeatArr = append(bubbleRepeatArr, uint32(ne.ID))
		} else {
			break
		}
	}
	return
}

func SortByUniqueEdgeLen(edgesArr []DBGEdge, nodesArr []DBGNode, Haplotype bool, kmerlen int) (sortedEIDLenArr []EIDLen) {
	bubbleRepeatEdgePathMinLen := 3 * kmerlen
	for i := range edgesArr {
		e := &edgesArr[i]
		if e.ID <= 0 || e.GetProcessFlag() > 0 {
			continue
		}

		if e.GetBubbleRepeatFlag() > 0 {
			var eIDLen EIDLen
			eIDLen.Length = e.GetSeqLen()
			eIDLen.EID = e.ID
			eIDLen.NextEID = e.ID
			var bubbleRepeatArr []uint32
			bubbleRepeatArr = append(bubbleRepeatArr, e.ID)
			// BACKWARD
			{
				ne := e
				for ne.ID > 0 {
					eIDArr := GetNextEIDArr(ne, PLUS, BACKWARD)
					if len(eIDArr) == 0 {
						break
					}
					ne = &edgesArr[eIDArr[0]]
					if ne.GetBubbleFlag() > 0 {
						eIDLen.Length += len(ne.Ks) - (kmerlen - 1)
					} else if ne.GetBubbleRepeatFlag() > 0 {
						eIDLen.Length += len(ne.Ks) - (kmerlen - 1)
						eIDLen.EID = ne.ID
						bubbleRepeatArr = append(bubbleRepeatArr, ne.ID)
					} else {
						break
					}
				}
			}

			// FORWARD
			{
				ne := e
				for ne.ID > 0 {
					eIDArr := GetNextEIDArr(ne, PLUS, FORWARD)
					if len(eIDArr) == 0 {
						break
					}
					ne = &edgesArr[eIDArr[0]]
					if ne.GetBubbleFlag() > 0 {
						eIDLen.Length += len(ne.Ks) - (kmerlen - 1)
					} else if ne.GetBubbleRepeatFlag() > 0 {
						eIDLen.Length += len(ne.Ks) - (kmerlen - 1)
						eIDLen.NextEID = ne.ID
						bubbleRepeatArr = append(bubbleRepeatArr, ne.ID)
					} else {
						break
					}
				}
			}

			if eIDLen.Length > bubbleRepeatEdgePathMinLen {
				for _, id := range bubbleRepeatArr {
					edgesArr[id].SetProcessFlag()
				}
				sortedEIDLenArr = append(sortedEIDLenArr, eIDLen)
			}
		}

		if edgesArr[e.ID].GetProcessFlag() == 0 && e.GetUniqueFlag() > 0 && e.StartNID != e.EndNID && e.GetTwoEdgesCycleFlag() == 0 {
			var il EIDLen
			il.EID = e.ID
			il.Length = e.GetSeqLen()
			edgesArr[e.ID].SetProcessFlag()
			sortedEIDLenArr = append(sortedEIDLenArr, il)
		}
	}

	sort.Sort(EIDLenArr(sortedEIDLenArr))
	return
}

/*func mergePathMat(edgesArr []DBGEdge, nodesArr []DBGNode) {
	for i, e := range edgesArr {
		if e.ID == 0 || e.GetDeleteFlag() > 0 || len(e.PathMat) == 0 {
			continue
		}

		if len(e.PathMat[0].IDArr) > 0 {
			var n Path
			edgesArr[i].PathMat[0] = n
		}
		if len(e.PathMat) <= 1 {
			continue
		}

		var consisPM [2]Path
		var SecondEID uint32 // if a two edges cycle, store e.EndNID second Edge
		if e.StartNID > 0 {
			ea := GetNearEdgeIDArr(nodesArr[e.StartNID], e.ID)
			if len(ea) == 1 {
				// check if a two edges cycle
				if e.EndNID > 0 {
					ea1 := GetNearEdgeIDArr(nodesArr[e.EndNID], e.ID)
					if len(ea1) == 1 && ea1[0] == ea[0] {
						a1 := GetOtherEArr(nodesArr[e.StartNID], e.ID)
						a2 := GetOtherEArr(nodesArr[e.EndNID], e.ID)
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

				ca := FindConsisPath(ea[0], e)
				if len(ca.IDArr) > 1 {
					consisPM[0] = ca
				}
			}
		}

		if e.EndNID > 0 {
			ea := GetNearEdgeIDArr(nodesArr[e.EndNID], e.ID)
			if len(ea) == 1 {
				if SecondEID > 0 {
					ea[0] = SecondEID
				}
				ca := FindConsisPath(ea[0], e)
				if len(ca.IDArr) > 1 {
					consisPM[1] = ca
				}
			}
		}

		var mP Path // merge Path
		if len(consisPM[0].IDArr) > 1 {
			mP.IDArr = GetReverseUint32Arr(consisPM[0].IDArr)
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
			edgesArr[i].PathMat = []Path{mP}
		} else {
			edgesArr[i].PathMat = nil
		}
		//fmt.Printf("[mergePathMat] mergePath[%v]: %v\nedge: %v\n", i, edgesArr[i].PathMat, e)

	}
} */

func CountNumEdgeIDInDbgMaxIntArr(arr []uint32) (num int) {
	for _, eID := range arr {
		if eID > 0 {
			num++
		}
	}
	return
}

func reallocExtendPArr(extendPArr [][]uint32, width int, direction uint8) {
	for i, arr := range extendPArr {
		na := make([]uint32, width)
		if direction == FORWARD {
			copy(na, arr)
		} else {
			x := len(arr)
			copy(na[width-x:], arr)
		}
		extendPArr[i] = na
	}
}

func ChangeExtendPArr(extendPArr [][]uint32, deleteCol int, width int, direction uint8) {
	if len(extendPArr) == 0 {
		return
	}
	if direction == FORWARD {
		for i, _ := range extendPArr {
			copy(extendPArr[i], extendPArr[i][deleteCol:])
		}
	} else {
		for i, _ := range extendPArr {
			copy(extendPArr[i][deleteCol:], extendPArr[i][:width-deleteCol])
		}
	}
}

func AppearEdgeInPathNum(arr []uint32, eID uint32) (count int) {
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

func IsInPathMat(ePathArr [][2][]uint32, matchPath [2][]uint32) bool {
	ml := len(matchPath[0])
	for _, p := range ePathArr {
		if len(p[0]) < ml {
			continue
		}
		if EqualUint32Arr(p[0][:ml], matchPath[0]) {
			return true
		}
	}
	return false
}

/*func GetOverlapPath(matchPath [2][]uint32, eID uint32, readPathArr []*ReadPath, edgesPathRelation []int, direction uint8, minMapFreq int) (extendP [2][]uint32, freq int) {
	Idx := IndexUint32(matchPath, eID)
	ePathArr := CopyReadPathArr(readPathArr, edgesPathRelation)
	for i := 0; i < len(ePathArr); i++ {
		p := ePathArr[i]
		pl, ml := len(p[0]), len(matchPath[0])
		if pl < ml {
			ePathArr[i][0] = nil
			ePathArr[i][1] = nil
			continue
		}
		idx := IndexUint32(p, eID)
		if Idx > 0 && idx > 0 {
			if Idx > 0 && idx > 0 && matchPath[0][Idx-1] != p[0][idx-1] {
				ReverseUint32Arr(p[0])
				ReverseUint32Arr(p[1])
				idx = pl - 1 - idx
			}
		} else if Idx < ml-1 && idx < pl-1 {
			if Idx+1 < len(matchPath[0]) && idx+1 < len(p[0]) && matchPath[0][Idx+1] != p[0][idx+1] {
				ReverseUint32Arr(p[0])
				ReverseUint32Arr(p[1])
				idx = pl - 1 - idx
			}
		} else if Idx > 0 && idx < pl-1 {
			if Idx > 0 && idx+1 < len(p[0]) && matchPath[0][Idx-1] == p[0][idx+1] {
				ReverseUint32Arr(p[0])
				ReverseUint32Arr(p[1])
				idx = pl - 1 - idx
			}
		} else {
			if Idx+1 < len(matchPath[0]) && idx-1 > 0 && matchPath[0][Idx+1] == p[0][idx-1] {
				ReverseUint32Arr(p[0])
				ReverseUint32Arr(p[1])
				idx = pl - 1 - idx
			}
		}
		if idx > Idx && ml-Idx < pl-idx {
			if EqualUint32Arr(matchPath[0], p[0][idx-Idx:idx+(ml-Idx)]) {
				if direction == FORWARD {
					ePathArr[i][0] = p[0][idx-Idx:]
					ePathArr[i][1] = p[1][idx-Idx:]
					//fmt.Printf("[GetOverlapPath]Idx: %v, idx:%v,i:%v, p:%v\n", Idx, idx, i, p)
				} else {
					ePathArr[i][0] = ReverseUint32Arr(p[0][:idx+(ml-Idx)])
					ePathArr[i][1] = ReverseUint32Arr(p[1][:idx+(ml-Idx)])
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
	}

	// found consis path
	if direction == BACKWARD {
		matchPath[0] = ReverseUint32Arr(matchPath[0])
		matchPath[1] = ReverseUint32Arr(matchPath[1])
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
			//ePathArr = CleanLowFreq(idFreqArr[y+1:], j, ePathArr, FORWARD)
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
			ePathArr = CleanLowFreq(idFreqArr[y+1:], j, ePathArr, FORWARD)
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
	}

	if direction == BACKWARD {
		matchPath[0] = ReverseUint32Arr(matchPath[0])
		matchPath[1] = ReverseUint32Arr(matchPath[1])
	}
	extendP = matchPath
	return
} */

type MatchInfo struct {
	Path        Path
	EID         uint32 // path cen EID
	PathPos     int
	MatchSeqLen int
	DiffSeqLen  int
	ExtendLen   int
}

type MatchArr []MatchInfo

func (arr MatchArr) Len() int {
	return len(arr)
}

func (arr MatchArr) Less(i, j int) bool {
	return arr[i].MatchSeqLen-arr[i].DiffSeqLen+arr[i].ExtendLen > arr[j].MatchSeqLen-arr[j].DiffSeqLen+arr[j].ExtendLen
}

func (arr MatchArr) Swap(i, j int) {
	arr[i], arr[j] = arr[j], arr[i]
}

func GetBestMatchPath(matchArr []MatchInfo, minMappingLen, MinExtendLen int) (bestmi MatchInfo) {
	sort.Sort(MatchArr(matchArr))
	for i, mi := range matchArr {
		fmt.Printf("[GetBestMatchPath]matchArr[%d]: EID: %v, MatchSeqLen: %d, DiffSeqLen: %d, ExtendSeqLen: %d, PathPos: %d\n", i, mi.EID, mi.MatchSeqLen, mi.DiffSeqLen, mi.ExtendLen, mi.PathPos)
	}
	for _, mi := range matchArr {
		//fmt.Printf("[GetBestMatchPath]matchArr[%d]: EID: %v, MatchSeqLen: %d, DiffSeqLen: %d, ExtendSeqLen: %d, PathPos: %d\n", i, mi.EID, mi.MatchSeqLen, mi.DiffSeqLen, mi.ExtendLen, mi.PathPos)
		if mi.MatchSeqLen-mi.DiffSeqLen > minMappingLen && mi.ExtendLen > MinExtendLen {
			bestmi = mi
			break
		}
	}
	return
}

func ExtendPathDDBG(edgesArr []DBGEdge, nodesArr []DBGNode, readPathArr []ReadPath, edgesPathRelationArr [][]uint32, minMapingFreq, minMappingLen, kmerlen int, eIDLen EIDLen, MaxPathLen, MaxBubbleSeqLen int) Path {
	e := edgesArr[eIDLen.EID]
	MinExtendLen := 2000
	var p Path
	//minMapEdgesNum := 3
	p.IDArr = make([]uint32, len(e.PathMat[0].IDArr))
	p.AltArr = make([]uint32, len(e.PathMat[0].IDArr))
	copy(p.IDArr, e.PathMat[0].IDArr)
	copy(p.AltArr, e.PathMat[0].AltArr)
	p.Freq = e.PathMat[0].Freq
	if p.Freq < minMapingFreq {
		return p
	}

	fmt.Printf("[ExtendPath]e.ID: %v, p.IDArr:%v\n\t\t\tp.AltArr: %v\n", e.ID, p.IDArr, p.AltArr)

	{
		idx := IndexUint32(p.IDArr, eIDLen.EID)
		// found left partition path
		//lastEID := uint32(0)
		//fmt.Printf("[ExtendPath]idx: %v\n", idx)
		for idx >= 0 && idx < len(p.IDArr) {
			var matchArr []MatchInfo
			i := idx - 1
			for ; i >= 0; i-- {
				e2 := edgesArr[p.IDArr[i]]
				/*if e2.GetUniqueFlag() > 0 {
				}*/
				//fmt.Printf("[ExtendPath] edgesArr[148] process Flag: %v, e2.ID: %d, i: %d, idx: %d\n", edgesArr[148].GetProcessFlag(), e2.ID, i, idx)
				if len(e2.PathMat) == 0 {
					continue
				}
				fmt.Printf("[ExtendPath]e2.ID:%v,  e2.BubbleRepeatFlag: %v, UniqueFlag: %v, e2 Len: %d, PathLen: %d\n", e2.ID, e2.GetBubbleRepeatFlag(), e2.GetUniqueFlag(), e2.GetSeqLen(), GetPathSeqLen(e2.PathMat[0].IDArr, edgesArr, kmerlen))

				/*if e2.GetProcessFlag() > 0 {
					log.Fatalf("[ExtendPath] found e2.ID: %v have been processed in p: %v\n", e2.ID, p)
				}*/
				//if e2.GetTwoEdgesCycleFlag() > 0 {
				//	step = 2
				//}
				for j := 0; j < 2; j++ {
					var tp Path
					if j == 0 {
						tp = e2.PathMat[0]
					} else {
						if len(e2.PathMat1) > 0 {
							tp = e2.PathMat1[0]
						}
					}
					if len(tp.IDArr) == 0 {
						continue
					}
					var ep Path
					ep.IDArr = make([]uint32, len(tp.IDArr))
					ep.AltArr = make([]uint32, len(tp.IDArr))
					copy(ep.IDArr, tp.IDArr)
					copy(ep.AltArr, tp.AltArr)
					ep.Freq = tp.Freq
					if AppearEdgeInPathNum(ep.IDArr, e2.ID) > 1 {
						fmt.Printf("[ExtendPath]e2.ID: %v appear more than one time in ep: %v\n", e2.ID, ep)
						continue
					}
					x := IndexUint32(ep.IDArr, e2.ID)
					if x < len(ep.IDArr)-1 && i < len(p.IDArr)-1 && ep.IDArr[x+1] != p.IDArr[i+1] {
						ReverseUint32Arr(ep.IDArr)
						ReverseUint32Arr(ep.AltArr)
						x = len(ep.IDArr) - 1 - x
					} else if x > 0 && i > 0 && ep.IDArr[x-1] != p.IDArr[i-1] {
						ReverseUint32Arr(ep.IDArr)
						ReverseUint32Arr(ep.AltArr)
						x = len(ep.IDArr) - 1 - x
					}
					if (len(ep.IDArr)-x > len(p.IDArr)-i) || x <= i {
						continue
					}
					fmt.Printf("[ExtendPath]i: %v, e2.ID: %v, ep.IDArr:%v\n\t\t\tep.AltArr: %v\n", i, e2.ID, ep.IDArr, ep.AltArr)
					//if EqualUint32Arr(p.IDArr[:i+(len(ep.IDArr)-x)], ep.IDArr[x-i:]) {
					if epPos, matchLen, diffLen, ok := EqualPath2(p, ep, e2.ID, BACKWARD, edgesArr, nodesArr, kmerlen, MaxPathLen, MaxBubbleSeqLen); ok {
						var mi MatchInfo
						mi.EID = e2.ID
						mi.Path = ep
						mi.DiffSeqLen = diffLen
						mi.MatchSeqLen = matchLen
						mi.PathPos = epPos
						mi.ExtendLen = GetPathSeqLen(ep.IDArr[:epPos], edgesArr, kmerlen)
						matchArr = append(matchArr, mi)
					}
				}
			}

			//fmt.Printf("[ExtendPath]Idx: %v, matchPathNum: %v, matchLen: %v\n", Idx, matchPathNum, matchLen)
			mok := false
			if len(matchArr) > 0 {
				bestmi := GetBestMatchPath(matchArr, minMappingLen, MinExtendLen)
				if bestmi.Path.Freq >= minMapingFreq {
					mok = true
					fmt.Printf("[ExtendPath]i: %v, extend left path: %v\n", i, bestmi.Path.IDArr[:bestmi.PathPos])
					p.IDArr = append(bestmi.Path.IDArr[:bestmi.PathPos], p.IDArr...)
					p.AltArr = append(bestmi.Path.AltArr[:bestmi.PathPos], p.AltArr...)
					p.Freq = bestmi.Path.Freq
					idx = IndexUint32(p.IDArr, bestmi.EID)
					edgesArr[bestmi.EID].SetProcessFlag()
				}
			}

			if !mok {
				idx = -1
			}
		}
		//ChangeExtendPArr(extendPArr, width-1-x, width, BACKWARD)
		// if not found merged edge, maybe can found unique edge path
		/*{
			var matchPath [2][]uint32
			var freq int
			var eID uint32 // anchor edge
			pathLen := kmerlen - 1
			for m := 0; m < len(p.IDArr); m++ {
				te := edgesArr[p.IDArr[m]]
				if eID < 2 && te.GetMergedFlag() == 0 && len(edgesPathRelationArr[te.ID]) > 0 {
					eID = te.ID
				}
				matchPath[0] = append(matchPath[0], p.IDArr[m])
				matchPath[1] = append(matchPath[1], p1.IDArr[m])
				pathLen += (len(te.Ks) - (kmerlen - 1))
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
			extendP, freq := GetOverlapPath(matchPath, eID, readPathArr, edgesPathRelationArr[eID], BACKWARD, minMapingFreq)
			//fmt.Printf("[ExtendPath]BACKWARD extendP: %v, freq:%v\n", extendP, freq)
			if len(extendP[0]) > ml && freq >= minMapingFreq {
				p.IDArr = append(extendP[0][:ml], p.IDArr...)
				p1.IDArr = append(extendP[1][:ml], p1.IDArr...)
				p.Freq = freq
				p1.Freq = freq
				idx = ml
			} else {
				//idx = i
				break
			}
		}*/
		//fmt.Printf("[ExtendPath]BACKWARD after extend unique edge path :%v\n", p)
	}

	// add right partion
	{
		idx := IndexUint32(p.IDArr, eIDLen.NextEID)
		//ne := edgesArr[eIDLen.NextEID]
		//lastEID := uint32(0)
		for idx >= 0 && idx < len(p.IDArr) {
			var matchArr []MatchInfo
			i := idx + 1
			for ; i < len(p.IDArr); i++ {
				e2 := edgesArr[p.IDArr[i]]
				//fmt.Printf("[ExtendPath]e2.ID: %v\n", e2.ID)
				if len(e2.PathMat) == 0 {
					continue
				}
				fmt.Printf("[ExtendPath]e2.ID:%v,  e2.BubbleRepeatFlag: %v, UniqueFlag: %v, e2 Len: %d, PathLen: %d\n", e2.ID, e2.GetBubbleRepeatFlag(), e2.GetUniqueFlag(), e2.GetSeqLen(), GetPathSeqLen(e2.PathMat[0].IDArr, edgesArr, kmerlen))
				/*if e2.GetProcessFlag() > 0 {
					log.Fatalf("[ExtendPath] found e2.ID: %v have been processed in p: %v\n", e2.ID, p)
				}*/
				for j := 0; j < 2; j++ {
					var tp Path
					if j == 0 {
						tp = e2.PathMat[0]
					} else {
						if len(e2.PathMat1) > 0 {
							tp = e2.PathMat1[0]
						}
					}
					if len(tp.IDArr) == 0 {
						continue
					}
					var ep Path
					ep.IDArr = make([]uint32, len(tp.IDArr))
					ep.AltArr = make([]uint32, len(tp.IDArr))
					copy(ep.IDArr, tp.IDArr)
					copy(ep.AltArr, tp.AltArr)
					ep.Freq = tp.Freq
					if AppearEdgeInPathNum(ep.IDArr, e2.ID) > 1 {
						fmt.Printf("[ExtendPath]i: %v, e2.ID: %v appear more than one time in ep: %v\n", i, e2.ID, ep)
						continue
					}
					x := IndexUint32(ep.IDArr, e2.ID)
					if x > 0 && i > 0 && ep.IDArr[x-1] != p.IDArr[i-1] {
						ReverseUint32Arr(ep.IDArr)
						ReverseUint32Arr(ep.AltArr)
						x = len(ep.IDArr) - 1 - x
					} else if x < len(ep.IDArr)-1 && i < len(p.IDArr)-1 && ep.IDArr[x+1] != p.IDArr[i+1] {
						ReverseUint32Arr(ep.IDArr)
						ReverseUint32Arr(ep.AltArr)
						x = len(ep.IDArr) - 1 - x
					}
					if (x > i) || (len(ep.IDArr)-x <= len(p.IDArr)-i) {
						continue
					}
					fmt.Printf("[ExtendPath]i: %v, e2.ID: %v, ep.IDArr: %v\n\t\t\tep.AltArr: %v\n", i, e2.ID, ep.IDArr, ep.AltArr)
					//if EqualUint32Arr(p.IDArr[i-x:], ep.IDArr[:x+(len(p.IDArr)-i)]) {
					if epPos, matchLen, diffLen, ok := EqualPath2(p, ep, e2.ID, FORWARD, edgesArr, nodesArr, kmerlen, MaxPathLen, MaxBubbleSeqLen); ok {
						var mi MatchInfo
						mi.EID = e2.ID
						mi.Path = ep
						mi.DiffSeqLen = diffLen
						mi.MatchSeqLen = matchLen
						mi.PathPos = epPos
						mi.ExtendLen = GetPathSeqLen(ep.IDArr[epPos:], edgesArr, kmerlen)
						matchArr = append(matchArr, mi)
					}
				}
			}
			//fmt.Printf("[ExtendPath]Idx: %v, matchPathNum: %v, matchLen: %v\n", Idx, matchPathNum, matchLen)
			mok := false
			if len(matchArr) > 0 {
				bestmi := GetBestMatchPath(matchArr, minMappingLen, MinExtendLen)
				if bestmi.Path.Freq >= minMapingFreq {
					mok = true
					p.IDArr = append(p.IDArr, bestmi.Path.IDArr[bestmi.PathPos:]...)
					p.AltArr = append(p.AltArr, bestmi.Path.AltArr[bestmi.PathPos:]...)
					p.Freq = bestmi.Path.Freq
					idx = IndexUint32(p.IDArr, bestmi.EID)
					edgesArr[bestmi.EID].SetProcessFlag()
					fmt.Printf("[ExtendPath]i: %v, extend right path: %v\n", i, bestmi.Path.IDArr[bestmi.PathPos:])
				}
			}

			if !mok {
				idx = -1
			}
			//ChangeExtendPArr(extendPArr, width-1-x, width, BACKWARD)
			//fmt.Printf("[ExtendPath]i: %v,len(p): %v, after extend rigth path p: %v\n", i, len(p.IDArr), p)
		}
		// if not found merged edge, maybe can found unique edge path
		/*{
			var matchPath [2][]uint32
			var freq int
			var eID uint32 // anchor edge
			pathLen := kmerlen - 1
			for m := len(p.IDArr) - 1; m > 0; m-- {
				te := edgesArr[p.IDArr[m]]
				if eID < 2 && te.GetMergedFlag() == 0 && len(edgesPathRelationArr[te.ID]) > 0 {
					eID = te.ID
				}
				matchPath[0] = append(matchPath[0], p.IDArr[m])
				matchPath[1] = append(matchPath[1], p1.IDArr[m])
				pathLen += (len(te.Ks) - (kmerlen - 1))
				if pathLen > minMappingLen && eID >= 2 {
					break
				}
			}
			if pathLen < minMappingLen || pathLen > 10000 || eID < 2 || lastEID == eID {
				break
			}
			lastEID = eID
			matchPath[0] = ReverseUint32Arr(matchPath[0])
			matchPath[1] = ReverseUint32Arr(matchPath[1])
			ml := len(matchPath[0])
			//fmt.Printf("[ExtendPath]FORWARD anchor edge ID: %v,len(edgesPathRelationArr[eID]): %v, pathLen: %v, matchPath: %v\n", eID, len(edgesPathRelationArr[eID]), pathLen, matchPath)
			//freq = p.Freq
			extendP, freq := GetOverlapPath(matchPath, eID, readPathArr, edgesPathRelationArr[eID], FORWARD, minMapingFreq)
			//fmt.Printf("[ExtendPath]FORWARD extendP: %v, freq:%v\n", extendP, freq)
			if len(extendP[0]) > ml && freq >= minMapingFreq {
				p.IDArr = append(p.IDArr, extendP[0][ml:]...)
				p1.IDArr = append(p1.IDArr, extendP[1][ml:]...)
				p.Freq = freq
				p1.Freq = freq
				//idx = i - 1
			} else {
				//idx = i
				break
			}
			//fmt.Printf("[ExtendPath]FORWARD after extend unique edge path :%v\n", p)
		}*/
	}

	return p
}

func FindMaxPath(sortedEIDLenArr []EIDLen, edgesArr []DBGEdge, nodesArr []DBGNode, readPathArr []ReadPath, edgesPathRelationArr [][]uint32, minMapingFreq, minMappingLen int, kmerlen int, Haplotype bool, MaxPathLen, MaxBubbleSeqLen int) (pA []Path) {
	//bubbleNearEdgeMinLen := 2000
	for _, item := range sortedEIDLenArr {
		e := edgesArr[item.EID]
		ne := edgesArr[item.NextEID]
		if e.GetProcessFlag() > 0 || ne.GetProcessFlag() > 0 || len(e.PathMat) == 0 || len(e.PathMat[0].IDArr) < minMapingFreq*2 || item.Length < 500 {
			continue
		}

		edgesArr[e.ID].SetProcessFlag()
		if item.NextEID > 0 {
			edgesArr[item.NextEID].SetProcessFlag()
		}

		/*if Haplotype {
			if e.GetBubbleFlag() > 0 {
				ok := true
				if e.StartNID > 0 {
					seID := GetNextEID(e.ID, nodesArr[e.StartNID])
					if seID > 0 && len(edgesArr[seID].Ks) < bubbleNearEdgeMinLen {
						ok = false
					}
				}
				if e.EndNID > 0 {
					eeID := GetNextEID(e.ID, nodesArr[e.EndNID])
					if eeID > 0 && len(edgesArr[eeID].Ks) < bubbleNearEdgeMinLen {
						ok = false
					}
				}
				if ok == false {
					continue
				}
			}
		}*/
		//fmt.Printf("[findMaxPath]ProcFlag:%v, delFlag:%v, uniqFlag: %v,TwoEdgesCycleFlag:%v, len(e.PathMat): %v, eID: %v, length: %v\n", e.GetProcessFlag(), e.GetDeleteFlag(), e.GetUniqueFlag(), e.GetTwoEdgesCycleFlag(), len(e.PathMat), item.Idx, item.Length)
		//fmt.Printf("[FindMaxPath] e.PathMat: %v\n", e.PathMat)
		//fmt.Printf("[findMaxPath] e.PathMat: %v\n", e.PathMat)
		fmt.Printf("[FindMaxPath] eID: %v, length: %v, NextEID: %v, path: %v\n\tpath1: %v\n", item.EID, item.Length, item.NextEID, e.PathMat, e.PathMat1)
		if Haplotype {
			maxP := ExtendPath(edgesArr, nodesArr, readPathArr, edgesPathRelationArr, minMapingFreq, minMappingLen, kmerlen, item, MaxPathLen, MaxBubbleSeqLen)
			if len(maxP.IDArr) > 2 {
				fmt.Printf("[FindMaxPath] maxP: %v\n", maxP)
				pA = append(pA, maxP)
			}
			if len(e.PathMat1) > 0 {
				e.PathMat, e.PathMat1 = e.PathMat1, e.PathMat
				maxP := ExtendPath(edgesArr, nodesArr, readPathArr, edgesPathRelationArr, minMapingFreq, minMappingLen, kmerlen, item, MaxPathLen, MaxBubbleSeqLen)
				if len(maxP.IDArr) > 2 {
					fmt.Printf("[FindMaxPath] maxP: %v\n", maxP)
					pA = append(pA, maxP)
				}
			}
		} else {
			if len(e.PathMat1) > 0 {
				continue
			}
			maxP := ExtendPath(edgesArr, nodesArr, readPathArr, edgesPathRelationArr, minMapingFreq, minMappingLen, kmerlen, item, MaxPathLen, MaxBubbleSeqLen)
			if len(maxP.IDArr) > 2 {
				fmt.Printf("[FindMaxPath] maxP: %v\n", maxP)
				pA = append(pA, maxP)
			}
		}
	}
	return pA
}

func Transform2ReadPath(mi MapingInfo) (rp ReadPath) {
	rp.ReadID = mi.ID
	for _, pa := range mi.PathMat {
		if len(pa) == 1 {
			arr := pa[0]
			for _, id := range arr {
				rp.Path0 = append(rp.Path0, id)
				rp.Path1 = append(rp.Path1, 0)
			}
		} else if len(pa) == 2 {
			if len(pa[1]) > len(pa[0]) {
				pa[0], pa[1] = pa[1], pa[0]
			}
			for _, id := range pa[0] {
				rp.Path0 = append(rp.Path0, id)
			}
			for _, id := range pa[1] {
				rp.Path1 = append(rp.Path1, id)
			}
			if len(pa[0]) > len(pa[1]) {
				diff := len(pa[0]) - len(pa[1])
				for j := 0; j < diff; j++ {
					rp.Path1 = append(rp.Path1, 0)
				}
				len1 := len(rp.Path1)
				rp.Path1[len1-1] = rp.Path1[len1-1-diff]
				for x := len1 - 1 - diff; x < len1-1; x++ {
					rp.Path1[x] = 1
				}
			}
		} else {
			log.Fatalf("[Transform2ReadPath]len(pa):%d pa:%v case unprocessed\n", len(pa), pa)
		}
	}
	return
}

func ConstructEdgesPathRelationship(edgesArr []DBGEdge, nodesArr []DBGNode, cs <-chan MapingInfo, mapRecordSize int) (readPathArr []ReadPath, edgesPathRelationArr [][]uint32) {
	readPathArr = make([]ReadPath, 0, mapRecordSize)
	edgesPathRelationArr = make([][]uint32, len(edgesArr))
	var mi MapingInfo
	var ok bool
	for {
		mi, ok = <-cs
		if !ok {
			break
		}
		rp := Transform2ReadPath(mi)
		readPathArr = append(readPathArr, rp)
		i := uint32(len(readPathArr) - 1)
		var extArr [2][]uint32
		extArr[0], extArr[1] = rp.Path0, rp.Path1
		if len(extArr[0]) < 3 {
			continue
		}
		for j, eID := range extArr[0] {
			e := &edgesArr[eID]
			if e.GetBubbleRepeatFlag() > 0 {
				edgesPathRelationArr[eID] = append(edgesPathRelationArr[eID], i)
			} else if e.GetUniqueFlag() > 0 && extArr[1][j] == 0 {
				edgesPathRelationArr[eID] = append(edgesPathRelationArr[eID], i)
			}
		}
	}

	if mapRecordSize != len(readPathArr) {
		fmt.Fprintf(os.Stderr, "[ConstructEdgesPathRelationship]mapRecordSize:%d != len(readPathArr):%d\n", mapRecordSize, len(readPathArr))
	}

	return
}

func CleanNGSPath(edgesArr []DBGEdge) {
	for i, e := range edgesArr {
		if i < 2 || e.GetDeleteFlag() > 0 {
			continue
		}
		edgesArr[i].PathMat = nil
	}
}

func ResetProcessFlagAll(edgesArr []DBGEdge) {
	for i, e := range edgesArr {
		if i < 2 || e.GetDeleteFlag() > 0 {
			continue
		}
		edgesArr[i].ResetProcessFlag()
	}
}

func CopyReadPathArr(readPathArr []*ReadPath, relationArr []int) [][2][]uint32 {
	subArr := make([][2][]uint32, len(relationArr))
	for i, idx := range relationArr {
		path0, path1 := readPathArr[idx].Path0, readPathArr[idx].Path1
		if len(path0) == 0 {
			continue
		}
		//srcArr := pathArr[idx]
		var extArr [2][]uint32
		extArr[0] = make([]uint32, len(path0))
		extArr[1] = make([]uint32, len(path1))
		for j, id := range path0 {
			extArr[0][j] = uint32(id)
		}
		for j, id := range path1 {
			extArr[1][j] = uint32(id)
		}
		subArr[i] = extArr
	}
	return subArr
}

func CopyPathArr(pathArr [][2][]uint32, relationArr []int) [][2][]uint32 {
	subArr := make([][2][]uint32, len(relationArr))
	for i, idx := range relationArr {
		//path0, path1 := readPathArr[idx].GetPath0(), readPathArr[idx].GetPath1()
		if len(pathArr[idx][0]) == 0 {
			continue
		}
		srcArr := pathArr[idx]
		var extArr [2][]uint32
		extArr[0] = make([]uint32, len(srcArr[0]))
		extArr[1] = make([]uint32, len(srcArr[1]))
		copy(extArr[0], srcArr[0])
		copy(extArr[1], srcArr[1])
		subArr[i] = extArr
	}
	return subArr
}

func IndexPathEID(pathArr [2][]uint32, eID uint32) int {
	for i, id := range pathArr[0] {
		if id == eID {
			return i
		}
	}
	for i, id := range pathArr[1] {
		if id == eID {
			return i
		}
	}
	return -1
}

func IndexBubbleRepeatArr(pathArr [2][]uint32, bubbleRepeatArr []uint32) (int, int) {
	idxP, idxB := -1, -1
	eID1, eID2 := bubbleRepeatArr[0], bubbleRepeatArr[len(bubbleRepeatArr)-1]
	for i, id := range pathArr[0] {
		if id == eID1 {
			idxP = i
			idxB = 0
			break
		} else if id == eID2 {
			idxP = i
			idxB = len(bubbleRepeatArr) - 1
			break
		} else {
			idx := IndexUint32(bubbleRepeatArr, id)
			if idx >= 0 {
				idxP = i
				idxB = idx
				break
			}
		}
	}

	return idxP, idxB
}

func IsInDirectionBubbleRepeatPath(eID1, eID2 uint32, edgesArr []DBGEdge, direction uint8) bool {
	ok := false
	if direction == FORWARD {
		ea := GetComingEdgeArr(edgesArr[eID1].EdgeIDOutcoming)
		if len(ea) == 0 {
			return false
		}
		ne := &edgesArr[ea[0]]
		if ne.GetBubbleFlag() == 0 {
			return false
		}
		if ne.StartNID == edgesArr[eID1].EndNID {
			ea = GetComingEdgeArr(ne.EdgeIDOutcoming)
		} else {
			ea = GetComingEdgeArr(ne.EdgeIDIncoming)
		}
		if len(ea) != 1 {
			return false
		}
		if uint32(ea[0]) == eID2 {
			return true
		}
	} else {
		ea := GetComingEdgeArr(edgesArr[eID1].EdgeIDIncoming)
		if len(ea) == 0 {
			return false
		}
		ne := &edgesArr[ea[0]]
		if ne.GetBubbleFlag() == 0 {
			return false
		}
		if ne.EndNID == edgesArr[eID1].StartNID {
			ea = GetComingEdgeArr(ne.EdgeIDIncoming)
		} else {
			ea = GetComingEdgeArr(ne.EdgeIDOutcoming)
		}
		if len(ea) != 1 {
			return false
		}
		if uint32(ea[0]) == eID2 {
			return true
		}
	}
	return ok
}

func GetNextEIDArr(e *DBGEdge, strand bool, direction uint8) (arr []uint32) {
	if direction == FORWARD {
		if strand {
			arr = GetComingEdgeArr(e.EdgeIDOutcoming)
		} else {
			arr = GetComingEdgeArr(e.EdgeIDIncoming)
		}
	} else {
		if strand {
			arr = GetComingEdgeArr(e.EdgeIDIncoming)
		} else {
			arr = GetComingEdgeArr(e.EdgeIDOutcoming)
		}
	}
	return
}

func IsInUint32Arr(arr []uint32, eID uint32) bool {
	for _, id := range arr {
		if id == eID {
			return true
		}
	}
	return false
}

func ReverseUint32Arr(arr []uint32) []uint32 {
	al := len(arr)
	dl := al / 2
	for i := 0; i < dl; i++ {
		arr[i], arr[al-1-i] = arr[al-1-i], arr[i]
	}

	return arr
}

func checkEdgePathArrDirection(edgesArr []DBGEdge, ePathArr [][2][]uint32, bubbleRepeatArr []uint32) [][2][]uint32 {
	var eID, nextEID uint32
	var strand1, strand2 bool
	var eIDLeftArr, eIDRightArr, nextEIDLeftArr, nextEIDRightArr []uint32
	eID = bubbleRepeatArr[0]
	e := &edgesArr[eID]
	if len(bubbleRepeatArr) > 1 {
		if IsInDirectionBubbleRepeatPath(eID, bubbleRepeatArr[1], edgesArr, FORWARD) {
			strand1 = PLUS
		} else {
			strand1 = MINUS
		}
		nextEID = bubbleRepeatArr[len(bubbleRepeatArr)-1]
		if IsInDirectionBubbleRepeatPath(nextEID, bubbleRepeatArr[len(bubbleRepeatArr)-2], edgesArr, BACKWARD) {
			strand2 = PLUS
		} else {
			strand2 = MINUS
		}
	} else {
		strand1 = PLUS
	}
	eIDLeftArr = GetNextEIDArr(e, strand1, BACKWARD)
	eIDRightArr = GetNextEIDArr(e, strand1, FORWARD)
	if len(bubbleRepeatArr) > 1 {
		ne := &edgesArr[nextEID]
		nextEIDLeftArr = GetNextEIDArr(ne, strand2, BACKWARD)
		nextEIDRightArr = GetNextEIDArr(ne, strand2, FORWARD)
	}

	if (e.StartNID == e.EndNID) || (e.GetTwoEdgesCycleFlag() > 0) {
		log.Fatalf("[checkEdgePathArrDirection] encounter cycle edge ID: %v\n", e.ID)
	}

	for i, extArr := range ePathArr {
		idx := IndexPathEID(extArr, eID)
		if idx >= 0 {
			if idx > 0 {
				if IsInUint32Arr(eIDRightArr, extArr[0][idx-1]) {
					ReverseUint32Arr(ePathArr[i][0])
					ReverseUint32Arr(ePathArr[i][1])
				}
			} else {
				if IsInUint32Arr(eIDLeftArr, extArr[0][idx+1]) {
					ReverseUint32Arr(ePathArr[i][0])
					ReverseUint32Arr(ePathArr[i][1])
				}
			}
			continue
		}

		idx = IndexPathEID(extArr, nextEID)
		if idx >= 0 {
			if idx > 0 {
				if IsInUint32Arr(nextEIDRightArr, extArr[0][idx-1]) {
					ReverseUint32Arr(ePathArr[i][0])
					ReverseUint32Arr(ePathArr[i][1])
				}
			} else {
				if IsInUint32Arr(nextEIDLeftArr, extArr[0][idx+1]) {
					ReverseUint32Arr(ePathArr[i][0])
					ReverseUint32Arr(ePathArr[i][1])
				}
			}
			continue
		}

		// path between [eID ~~~~ nextEID]
		for j, id := range extArr[0] {
			te := &edgesArr[id]
			if te.GetBubbleRepeatFlag() == 0 {
				continue
			}
			idx1 := IndexUint32(bubbleRepeatArr, uint32(te.ID))
			if idx1 < 0 {
				continue
				//log.Fatalf("[checkEdgePathArrDirection]bubbleRepeatArr:%v id:%d extArr[0]:%v\n", bubbleRepeatArr, id, extArr[0])
			}
			var strand bool
			if idx1 < len(bubbleRepeatArr)-1 {
				if IsInDirectionBubbleRepeatPath(bubbleRepeatArr[idx1], bubbleRepeatArr[idx1+1], edgesArr, FORWARD) {
					strand = PLUS
				} else {
					strand = MINUS
				}
			} else {
				if IsInDirectionBubbleRepeatPath(bubbleRepeatArr[idx1], bubbleRepeatArr[idx1-1], edgesArr, BACKWARD) {
					strand = PLUS
				} else {
					strand = MINUS
				}
			}
			eIDRightArr = GetNextEIDArr(te, strand, FORWARD)
			if j > 0 {
				if IsInUint32Arr(eIDRightArr, extArr[0][j-1]) {
					ReverseUint32Arr(ePathArr[i][0])
					ReverseUint32Arr(ePathArr[i][1])
				}
			} else {
				if !IsInUint32Arr(eIDRightArr, extArr[0][j+1]) {
					ReverseUint32Arr(ePathArr[i][0])
					ReverseUint32Arr(ePathArr[i][1])
				}
			}
			break
		}
	}

	return ePathArr
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

/*func AddIDFreqArr(idFreqArr []IDFreq, ea []uint32) []IDFreq {
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

func GetTotalFreq(idFreqArr []IDFreq) (totalFreq uint32) {
	for _, idFreq := range idFreqArr {
		totalFreq += idFreq.Freq
	}
	return
}

func IsInIDFreqArr(idFreqArr []IDFreq, eID uint32) bool {
	for _, idFeq := range idFreqArr {
		if idFeq.ID == eID {
			return true
		}
	}
	return false
}

func MergeLongEdgePathMat(e DBGEdge, edgesArr []DBGEdge, nodesArr []DBGNode) (mergePath, mergePath1 Path) {
	var ea1, ea2 []uint32
	var pm, pm1 [2]Path
	pm[0].IDArr = make([]uint32, len(e.PathMat[0].IDArr))
	pm[1].IDArr = make([]uint32, len(e.PathMat[1].IDArr))
	pm1[0].IDArr = make([]uint32, len(e.PathMat1[0].IDArr))
	pm1[1].IDArr = make([]uint32, len(e.PathMat1[1].IDArr))
	copy(pm[0].IDArr, e.PathMat[0].IDArr)
	copy(pm[1].IDArr, e.PathMat[1].IDArr)
	copy(pm1[0].IDArr, e.PathMat1[0].IDArr)
	copy(pm1[1].IDArr, e.PathMat1[1].IDArr)
	pm[0].Freq, pm[1].Freq = e.PathMat[0].Freq, e.PathMat[1].Freq
	pm1[0].Freq, pm1[1].Freq = e.PathMat1[0].Freq, e.PathMat1[1].Freq
	if IsInComing(nodesArr[e.StartNID].EdgeIDIncoming, e.ID) {
		ea1 = GetNearEdgeIDArr(nodesArr[e.StartNID], e.ID, true)
	} else {
		ea1 = GetNearEdgeIDArr(nodesArr[e.StartNID], e.ID, false)
	}
	if IsInComing(nodesArr[e.EndNID].EdgeIDIncoming, e.ID) {
		ea2 = GetNearEdgeIDArr(nodesArr[e.EndNID], e.ID, true)
	} else {
		ea2 = GetNearEdgeIDArr(nodesArr[e.EndNID], e.ID, false)
	}
	if len(ea1) != 1 || len(ea2) != 1 {
		return
	}
	idx1 := IndexUint32(pm[0].IDArr, e.ID)
	idx2 := IndexUint32(pm[1].IDArr, e.ID)
	if idx1 == 0 && len(pm[1].IDArr)-idx2 == len(pm[0].IDArr) {
		if EqualUint32Arr(pm[0].IDArr, pm[1].IDArr[:len(pm[1].IDArr)-idx2]) {
			mergePath.IDArr = pm[1].IDArr
			mergePath1.IDArr = pm1[1].IDArr
			mergePath.Freq += pm[0].Freq
			mergePath1.Freq += pm1[0].Freq
		}
	} else if idx2 == 0 && len(pm[0].IDArr)-idx1 == len(pm[1].IDArr) {
		if EqualUint32Arr(pm[1].IDArr, pm[0].IDArr[:len(pm[0].IDArr)-idx1]) {
			mergePath.IDArr = pm[0].IDArr
			mergePath1.IDArr = pm1[0].IDArr
			mergePath.Freq += pm[1].Freq
			mergePath1.Freq += pm1[1].Freq
		}
	}
	ok1, ok2 := false, false
	if idx1 == 0 {
		if pm[0].IDArr[idx1+1] == ea1[0] {
			pm[0].IDArr = ReverseUint32Arr(pm[0].IDArr)
			pm1[0].IDArr = ReverseUint32Arr(pm1[0].IDArr)
			idx1 = len(pm[0].IDArr) - 1
			ok1 = true
		} else if pm[0].IDArr[idx1+1] == ea2[0] {
			pm[0], pm[1] = pm[1], pm[0]
			pm1[0], pm1[1] = pm1[1], pm1[0]
			idx1, idx2 = idx2, idx1
			if idx1 == 0 {
				if pm[0].IDArr[idx1+1] == ea1[0] {
					pm[0].IDArr = ReverseUint32Arr(pm[0].IDArr)
					pm1[0].IDArr = ReverseUint32Arr(pm1[0].IDArr)
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
					pm[0].IDArr = ReverseUint32Arr(pm[0].IDArr)
					pm1[0].IDArr = ReverseUint32Arr(pm1[0].IDArr)
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
				pm[1].IDArr = ReverseUint32Arr(pm[1].IDArr)
				pm1[1].IDArr = ReverseUint32Arr(pm1[1].IDArr)
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
	PathArr [][2][]uint32
	Idx     int
}

func CleanLowFreq(idFreqArr []IDFreq, j int, pathArr [][2][]uint32, direction uint8) [][2][]uint32 {
	for z := 0; z < len(idFreqArr); z++ {
		id := idFreqArr[z].ID
		for w := 0; w < len(pathArr); w++ {
			p := pathArr[w]
			if direction == BACKWARD {
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

func GetIDFreqArr(ePathArr [][2][]uint32, j int) (idFreqArr []IDFreq) {
	for k := 0; k < len(ePathArr); k++ {
		p := ePathArr[k]
		if p[0][j] > 0 && p[1][j] == 0 {
			idFreqArr = AddIDFreqArr(idFreqArr, p[0][j], 1)
		}
	}
	return
}

func GetIDFreqArrHaplotype(ePathArr [][2][]uint32, j int, direction uint8, diffLen int, path, altPath []uint32, UniqueScore int) (idFreqArr []IDFreq) {
	if direction == FORWARD {
		for k := 0; k < len(ePathArr); k++ {
			p := ePathArr[k]
			ok := true
			y := len(path) - 1
			for x := j - 1; x >= j-diffLen && y >= 0; x-- {
				if (x < len(p[0]) && p[0][x] == path[y]) || (x < len(p[1]) && p[1][x] == path[y]) {

				} else if y < len(altPath) {
					if p[0][x] == altPath[y] || p[1][x] == altPath[y] {

					} else {
						ok = false
						break
					}
				} else {
					ok = false
					break
				}
				y--
			}
			if ok && len(p[0]) > j && len(p[1]) > j {
				if p[0][j] > 0 && p[1][j] == 0 {
					idFreqArr = AddIDFreqArr(idFreqArr, p[0][j], uint32(UniqueScore))
				} else if p[0][j] > 0 && p[1][j] > 0 {
					idFreqArr = AddIDFreqArr(idFreqArr, p[0][j], 1)
					idFreqArr = AddIDFreqArr(idFreqArr, p[1][j], 1)
				}
			}
		}
	} else {
		for k := 0; k < len(ePathArr); k++ {
			p := ePathArr[k]
			ok := true
			y := 0
			for x := j + 1; x < j+1+diffLen && y < len(path); x++ {
				if p[0][x] == path[y] || p[1][x] == path[y] {

				} else if y >= len(path)-len(altPath) {
					sz := len(path) - len(altPath)
					if p[0][x] == altPath[y-sz] || p[1][x] == altPath[y-sz] {

					} else {
						ok = false
						break
					}
				} else {
					ok = false
					break
				}
				y++
			}
			if ok && len(p[0]) > j && len(p[1]) > j {
				if p[0][j] > 0 && p[1][j] == 0 {
					idFreqArr = AddIDFreqArr(idFreqArr, p[0][j], uint32(UniqueScore))
				} else if p[0][j] > 0 && p[1][j] > 0 {
					idFreqArr = AddIDFreqArr(idFreqArr, p[0][j], 1)
					idFreqArr = AddIDFreqArr(idFreqArr, p[1][j], 1)
				}
			}
		}
	}

	return
}

func ReplaceEPathMat(ePathMatArr [][2][]uint32, idx int, direction uint8, mp, oldP, newP []uint32) {
	fmt.Printf("[ReplaceEPathMat]mp: %v\n\t\toldP: %v\n\t\tnewP: %v\n", mp, oldP, newP)

	if direction == BACKWARD {
		oldP = GetReverseUint32Arr(oldP)
		newP = GetReverseUint32Arr(newP)
		for i, pm := range ePathMatArr {
			if reflect.DeepEqual(pm[0][idx+1:idx+1+len(mp)], mp) && reflect.DeepEqual(pm[0][idx+1-len(oldP):idx+1], oldP) {
				bl := idx + 1 - len(newP)
				np := make([]uint32, len(pm[0]))
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
				np := make([]uint32, len(pm[0]))
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

func AdjustAltPath(pm []uint32, idx int, direction uint8, oldP, newP []uint32) (path []uint32) {
	if direction == BACKWARD {
		dl := len(newP) - len(oldP)
		if idx <= MaxInt(len(newP), len(oldP)) {
			path = pm
			return
		}
		path = make([]uint32, 0, len(pm))
		if dl > 0 {
			path = append(path, pm[dl:idx+1-len(oldP)]...)
			path = path[:len(path)+dl]
			path = append(path, pm[idx+1-len(oldP):idx+1-len(oldP)+(len(pm)-len(path))]...)
		} else {
			path = append(path, pm[:idx+1-len(oldP)]...)
			path = append(path, pm[idx+1-len(newP):]...)
		}
	} else {
		if idx >= len(pm)-MaxInt(len(newP), len(oldP)) {
			path = pm
			return
		}
		dl := len(newP) - len(oldP)
		path = make([]uint32, 0, len(pm))
		if dl > 0 {
			path = append(path, pm[:idx+len(oldP)]...)
			path = path[:len(path)+dl]
			path = append(path, pm[idx+len(oldP):idx+len(oldP)+(len(pm)-len(path))]...)
		} else {
			path = append(path, pm[:idx+len(newP)]...)
			path = append(path, pm[idx+len(oldP):]...)
		}
	}
	path = path[:len(pm)]

	return
}

func ReplaceEPath(pm []uint32, idx int, direction uint8, oldP, newP []uint32) (path []uint32) {
	//fmt.Printf("[ReplaceEPath]idx: %v, oldP: %v\n\t\tnewP: %v\n", idx, oldP, newP)
	if direction == BACKWARD {
		oldP = GetReverseUint32Arr(oldP)
		newP = GetReverseUint32Arr(newP)
		bl := idx + 1 - len(newP)
		if bl < 0 {
			return pm
		}
		path = make([]uint32, len(pm))
		x := idx + 1 - len(oldP) - 1
		for j := bl - 1; j >= 0 && x >= 0; j-- {
			path[j] = pm[x]
			x--
		}
		copy(path[idx+1-len(newP):idx+1], newP)
		copy(path[idx+1:], pm[idx+1:])
	} else {
		path = make([]uint32, len(pm))
		copy(path[:idx], pm[:idx])
		copy(path[idx:idx+len(newP)], newP)
		bl := idx + len(newP)
		if bl >= len(pm) {
			return pm
		}
		x := idx + len(oldP)
		for j := bl; j < len(path) && x < len(pm); j++ {
			path[j] = pm[x]
			x++
		}
	}
	//fmt.Printf("[ReplaceEPath]changed path: %v\n", path)
	return
}

func FoundFirstEID(arr []uint32) (idx int) {
	idx = -1
	for i, id := range arr {
		if id > 0 {
			idx = i
			break
		}
	}
	return
}

func FoundLastEID(arr []uint32) (idx int) {
	idx = -1
	for i := len(arr) - 1; i >= 0; i-- {
		if arr[i] > 0 {
			idx = i
			break
		}
	}
	return
}

func GetEPathSectionPath(pm []uint32, idx int, direction uint8, mp []uint32, diffLen, MaxPathLen int) (path []uint32) {
	if direction == FORWARD {
		pl := MaxPathLen
		if pl > len(pm)-idx {
			pl = len(pm) - idx
		}
		if diffLen > 0 {
			if idx-diffLen < 0 {
				return
			}
			if reflect.DeepEqual(pm[idx-diffLen:idx], mp[len(mp)-diffLen:]) {
				path = make([]uint32, pl)
				copy(path, pm[idx:idx+pl])
			}
		} else {
			path = make([]uint32, pl)
			copy(path, pm[idx:idx+pl])
		}
	} else {
		pl := MaxPathLen
		if pl > idx+1 {
			pl = idx + 1
		}
		if diffLen > 0 {
			if idx+1+diffLen > len(pm) {
				return
			}
			if reflect.DeepEqual(pm[idx+1:idx+1+diffLen], mp[:diffLen]) {
				path = GetReverseUint32Arr(pm[idx+1-pl : idx+1])
			}
		} else {
			path = GetReverseUint32Arr(pm[idx+1-pl : idx+1])
		}
	}
	return
}

func GetEPathMatPath(ePathMatArr [][2][]uint32, idx int, direction uint8, mp []uint32, eID uint32, diffLen int) (path []uint32) {
	slen := len(ePathMatArr[0][0])
	path = append(path, eID)
	if direction == FORWARD {
		for j := idx + 1; j < slen; j++ {
			//idFreqArr := GetIDFreqArr(ePathMatArr, j)
			var idFreqArr []IDFreq
			for _, pm := range ePathMatArr {
				pmIdx := FoundFirstEID(pm[0])
				if idx-pmIdx < diffLen {
					continue
				}
				if pm[0][idx] == eID {
					if diffLen > 0 {
						if reflect.DeepEqual(pm[0][idx-diffLen:idx], mp[len(mp)-diffLen:]) {
							idFreqArr = AddIDFreqArr(idFreqArr, pm[0][j], 1)
						}
					} else {
						idFreqArr = AddIDFreqArr(idFreqArr, pm[0][j], 1)
					}
				}
			}
			if len(idFreqArr) > 0 {
				sort.Sort(IDFreqArr(idFreqArr))
				if len(idFreqArr) > 1 && idFreqArr[0].Freq == idFreqArr[1].Freq {
					break
				}
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
				pmIdx := FoundLastEID(pm[0])
				if pmIdx+1-(idx+1) < diffLen {
					continue
				}
				if pm[0][idx] == eID {
					if diffLen > 0 {
						if reflect.DeepEqual(pm[0][idx+1:idx+1+diffLen], mp[:diffLen]) {
							idFreqArr = AddIDFreqArr(idFreqArr, pm[0][j], 1)
						}
					} else {
						idFreqArr = AddIDFreqArr(idFreqArr, pm[0][j], 1)
					}
				}
			}
			if len(idFreqArr) > 0 {
				sort.Sort(IDFreqArr(idFreqArr))
				if len(idFreqArr) > 1 && idFreqArr[0].Freq == idFreqArr[1].Freq {
					break
				}
				path = append(path, idFreqArr[0].ID)
			} else {
				break
			}
		}
	}
	return
}

func AdjustEPathMat(ePathMatArr [][2][]uint32, idx int, idFreqArr []IDFreq, direction uint8, edgesArr []DBGEdge, mp []uint32, kmerlen int) {
	mp = mp[:len(mp)-1]
	if direction == BACKWARD {
		mp = GetReverseUint32Arr(mp)
	}
	path0 := GetEPathMatPath(ePathMatArr, idx, direction, mp, idFreqArr[0].ID, len(mp))
	path1 := GetEPathMatPath(ePathMatArr, idx, direction, mp, idFreqArr[1].ID, len(mp))
	fmt.Printf("[AdjustEPathMat]mp: %v\n\t\tpath0: %v\n\t\tpath1: %v\n", mp, path0, path1)
	if len(path0) > 1 && len(path1) > 1 {
		l0 := len(edgesArr[path0[0]].Ks) - (kmerlen - 1)
		//l0 += len(edgesArr[path0[1]].Ks) - (kmerlen-1)
		l1 := len(edgesArr[path1[0]].Ks) - (kmerlen - 1)
		//l1 += len(edgesArr[path1[1]].Ks) - (kmerlen-1)
		j := 1
		for i := 1; i < len(path0) && j < len(path1); i++ {
			if path0[i] == path1[j] {
				if AbsInt(l0-l1) < Min(l0, l1)/10 {
					ReplaceEPathMat(ePathMatArr, idx, direction, mp, path1[:j], path0[:i])
					break
				}
			}
			if l1 < l0 {
				ok := false
				for ; j < len(path1); j++ {
					if path0[i] == path1[j] {
						if AbsInt(l0-l1) < Min(l0, l1)/10 {
							ReplaceEPathMat(ePathMatArr, idx, direction, mp, path1[:j], path0[:i])
							ok = true
							break
						}
					} else if l1 > l0 {
						break
					}
					l1 += len(edgesArr[path1[j]].Ks) - (kmerlen - 1)
				}
				if ok {
					break
				}
			}
			l0 += len(edgesArr[path0[i]].Ks) - (kmerlen - 1)
		}
		// print
		for i, pm := range ePathMatArr {
			fmt.Printf("[AdjustEPathMat]ePathMatArr[%v]: %v\n", i, pm)
		}
	}
	return
}

func GetShareDBGNode(nodesArr []DBGNode, le, ne *DBGEdge) (nd DBGNode) {
	if le.StartNID == ne.StartNID || le.StartNID == ne.EndNID {
		nd = nodesArr[le.StartNID]
	} else if le.EndNID == ne.StartNID || le.EndNID == ne.EndNID {
		nd = nodesArr[le.EndNID]
	}
	return
}
func GetShareDBGNID(le, ne *DBGEdge) (nID uint32) {
	if le.StartNID == ne.StartNID || le.StartNID == ne.EndNID {
		nID = le.StartNID
	} else if le.EndNID == ne.StartNID || le.EndNID == ne.EndNID {
		nID = le.EndNID
	}
	return
}

func SortByFreq(pathArr []BubblePath, idFreqArr []IDFreq) []BubblePath {
	count := 0
	maxLen := 0
	for i, bp := range pathArr {
		if bp.Path[0] == idFreqArr[0].ID {
			if len(bp.Path) > maxLen {
				pathArr[0], pathArr[i] = pathArr[i], pathArr[0]
				maxLen = len(bp.Path)
			}
			count++
		}
	}
	if count != 1 {
		fmt.Fprintf(os.Stderr, "[SortByFreq] sort by freq error.... idFreqArr:%v, pathArr:%v\n", idFreqArr, pathArr)
		//return nil
	}
	return pathArr
}

func IsInPathArr(path []uint32, pathArr []BubblePath, edgesArr []DBGEdge, kmerlen int) (np []uint32, ok bool) {
	for _, bp := range pathArr {
		bl := len(bp.Path)
		if len(path) >= bl && reflect.DeepEqual(path[:bl], bp.Path) {
			np = path[:bl]
			ok = true
			break
		}
	}

	return
}

func CleanPathFlank(pm []uint32, idx int, direction uint8) []uint32 {
	if direction == FORWARD {
		for i := idx; i < len(pm); i++ {
			if pm[i] == 0 {
				break
			}
			pm[i] = 0
		}
	} else {
		for i := idx; i >= 0; i-- {
			if pm[i] == 0 {
				break
			}
			pm[i] = 0
		}
	}
	return pm
}

func AdjustEPathMat2(ePathMatArr [][2][]uint32, idx int, idFreqArr []IDFreq, direction uint8, edgesArr []DBGEdge, nodesArr []DBGNode, mp []uint32, kmerlen int, diffLen int, MaxPathLen, MaxBubbleSeqLen int) {
	/*if direction == FORWARD {
		mp = mp[:len(mp)-1]
		path0 := GetEPathMatPath(ePathMatArr, idx, direction, mp, eID, diffLen)
		for i, pm := range ePathMatArr {
			path1 := GetEPathSectionPath(pm[0], idx, direction, mp, diffLen)

			pmIdx := FoundFirstEID(pm[0])
			if idx-pmIdx >= diffNum {
				if reflect.DeepEqual(pm[0][pmIdx:pmIdx+(idx-pmIdx)], mp[len(mp)-(idx-pmIdx):]) {
					if pm[0][idx] != eID {
						ePathMatArr[i][0][idx] = eID
					}
					ePathMatArr[i][1][idx] = 0
				}
			}
		}
	} else { // direction == BACKWARD

	} */
	if direction == FORWARD {
		mp = mp[:len(mp)-1]
	} else {
		mp = mp[1:]
	}
	/*if direction == BACKWARD {
		mp = GetReverseUint32Arr(mp)
	}*/
	var pathArr []BubblePath
	e1, e2 := edgesArr[idFreqArr[0].ID], edgesArr[idFreqArr[1].ID]
	var le *DBGEdge
	if direction == FORWARD {
		le = &edgesArr[mp[len(mp)-1]]
		//nd := GetShareDBGNode(nodesArr, le, edgesArr[idFreqArr[0].ID])
	} else {
		le = &edgesArr[mp[0]]
	}
	var path0, path1 []uint32
	if e1.GetTwoEdgesCycleFlag() > 0 || e2.GetTwoEdgesCycleFlag() > 0 || e1.StartNID == e1.EndNID || e2.StartNID == e2.EndNID {
		if e1.GetTwoEdgesCycleFlag() > 0 {
			path0 = append(path0, e1.ID, le.ID, e2.ID)
			path1 = append(path1, e2.ID)
		} else if e2.GetTwoEdgesCycleFlag() > 0 {
			path0 = append(path0, e1.ID)
			path1 = append(path1, e2.ID, le.ID, e1.ID)
		} else if e1.StartNID == e1.EndNID {
			path0 = append(path0, e1.ID, e2.ID)
			path1 = append(path1, e2.ID)
		} else if e2.StartNID == e2.EndNID {
			path0 = append(path0, e1.ID)
			path1 = append(path1, e2.ID, e1.ID)
		}
		fmt.Printf("[AdjustEPathMat2] path0: %v, path1:%v\n", path0, path1)
	} else if len(mp) >= 1 {
		nd := GetShareDBGNode(nodesArr, le, &edgesArr[idFreqArr[0].ID])
		pa := GetNextPathArr(le, nd, MaxPathLen+1, MaxBubbleSeqLen, edgesArr, nodesArr, kmerlen)

		pathArr, _ = GetBubblePathArr(pa, edgesArr, kmerlen, MaxPathLen, MaxBubbleSeqLen)
		for x, p := range pathArr {
			fmt.Printf("[GetDiffPathLenBubbleArr]pathArr[%v]: %v\n", x, p)
		}
		pathArr = SortByFreq(pathArr, idFreqArr)
		if len(pathArr) < 2 {
			fmt.Printf("[AdjustEPathMat2] not found pathArr: %v\n", pathArr)
			//return
		}
		fmt.Printf("[AdjustEPathMat2] pathArr: %v\n", pathArr)
	}
	//path0 := GetEPathMatPath(ePathMatArr, idx, direction, mp, eID, diffLen)
	for i, pm := range ePathMatArr {
		if pm[0][idx] == 0 || pm[0][idx] == idFreqArr[0].ID {
			continue
		}
		if e1.GetTwoEdgesCycleFlag() > 0 || e2.GetTwoEdgesCycleFlag() > 0 || e1.StartNID == e1.EndNID || e2.StartNID == e2.EndNID {
			ePathMatArr[i][0] = ReplaceEPath(pm[0], idx, direction, path1, path0)
			ePathMatArr[i][1] = AdjustAltPath(pm[1], idx, direction, path1, path0)
		} else if len(mp) >= 1 {
			if len(pathArr) > 1 {
				path1 := GetEPathSectionPath(pm[0], idx, direction, mp, diffLen, MaxPathLen)
				if path, ok := IsInPathArr(path1, pathArr[1:], edgesArr, kmerlen); ok {
					ePathMatArr[i][0] = ReplaceEPath(pm[0], idx, direction, path, pathArr[0].Path)
					ePathMatArr[i][1] = AdjustAltPath(pm[1], idx, direction, path, pathArr[0].Path)
				}
			} else {
				ePathMatArr[i][0] = CleanPathFlank(pm[0], idx, direction)
				//ePathMatArr[i][1] = CleanPathFlank(pm[1], idx, direction)
			}
		}

	}
	// print
	/*for i, pm := range ePathMatArr {
		fmt.Printf("[AdjustEPathMat2]ePathMatArr[%v]: %v\n", i, pm)
	}*/
	return
}

func GetBubbleRepeatArr(e *DBGEdge, edgesArr []DBGEdge) (bubbleRepeatArr []uint32) {
	ea1 := GetDirectionBubbleRepeatArr(e, edgesArr, BACKWARD)
	ea2 := GetDirectionBubbleRepeatArr(e, edgesArr, FORWARD)
	bubbleRepeatArr = make([]uint32, 0, len(ea1)+len(ea2)+1)
	bubbleRepeatArr = append(bubbleRepeatArr, uint32(e.ID))
	if len(ea1) > 0 && len(ea2) > 0 {
		log.Fatalf("[GetBubbleRepeatArr]ea1:%v ea2:%v ea1 and ea2 must one empty\n", ea1, ea2)
	} else if len(ea1) > 0 {
		//bubbleRepeatArr = append(bubbleRepeatArr, GetReverseUint32Arr(ea1)...)
		bubbleRepeatArr = append(bubbleRepeatArr, ea1...)
	} else if len(ea2) > 0 {
		bubbleRepeatArr = append(bubbleRepeatArr, ea2...)
	}
	return
}

func GetPathRelationArr(rp ReadPath) [2][]uint32 {
	var extArr [2][]uint32
	extArr[0] = make([]uint32, len(rp.Path0))
	extArr[1] = make([]uint32, len(rp.Path1))
	copy(extArr[0], rp.Path0)
	copy(extArr[1], rp.Path1)
	return extArr
}

func GetBubbleRepeatReadPathArr(bubbleRepeatArr []uint32, edgesPathRelationArr [][]uint32, readPathArr []ReadPath) (ePathArr [][2][]uint32) {
	readPathMap := make(map[uint32]bool)
	for _, eID := range bubbleRepeatArr {
		for _, pr := range edgesPathRelationArr[eID] {
			if _, ok := readPathMap[pr]; !ok {
				ePathArr = append(ePathArr, GetPathRelationArr(readPathArr[pr]))
				readPathMap[pr] = true
			}
		}
	}
	return
}

func GetCorrelationEID(ePathArr [][2][]uint32, eID1, eID2 uint32, idx1, idx2 int) (correlationNum int) {
	for _, path := range ePathArr {
		if path[0][idx1] == eID1 && path[1][idx1] == 0 && path[0][idx2] == eID2 && path[1][idx2] == 0 {
			correlationNum++
		}
	}
	return
}

func GetHaplotypePath(ePathArr [][2][]uint32, idx int, mp [2][]uint32, direction uint8, edgesArr []DBGEdge, eID1, eID2 uint32) (h1, h2 uint32, ok bool) {
	ok = true
	if direction == FORWARD {
		px := -1
		for i := len(mp[0]) - 1; i >= 0; i-- {
			if mp[0][i] != mp[1][i] {
				px = i
				break
			}
		}
		if px < 0 {
			ok = false
			return
		}
		mp0TeID1 := GetCorrelationEID(ePathArr, mp[0][px], eID1, idx-(len(mp[0])-px), idx)
		mp0TeID2 := GetCorrelationEID(ePathArr, mp[0][px], eID2, idx-(len(mp[0])-px), idx)
		mp1TeID1 := GetCorrelationEID(ePathArr, mp[1][px], eID1, idx-(len(mp[1])-px), idx)
		mp1TeID2 := GetCorrelationEID(ePathArr, mp[1][px], eID2, idx-(len(mp[1])-px), idx)
		if mp0TeID1*2 < mp0TeID2 {
			h1 = eID2
		} else if mp0TeID1 > mp0TeID2*2 {
			h1 = eID1
		} else {
			ok = false
		}
		if mp1TeID1*2 < mp1TeID2 {
			h2 = eID2
		} else if mp1TeID1 > mp1TeID2*2 {
			h2 = eID1
		} else {
			ok = false
		}
		if h1 == h2 {
			ok = false
		}
	} else { // direction == BACKWARD
		px := -1
		for i := 0; i < len(mp[0]); i++ {
			if mp[0][i] != mp[1][i] {
				px = i
				break
			}
		}
		if px < 0 {
			ok = false
			return
		}
		mp0TeID1 := GetCorrelationEID(ePathArr, eID1, mp[0][px], idx, idx+1+px)
		mp0TeID2 := GetCorrelationEID(ePathArr, eID2, mp[0][px], idx, idx+1+px)
		mp1TeID1 := GetCorrelationEID(ePathArr, eID1, mp[1][px], idx, idx+1+px)
		mp1TeID2 := GetCorrelationEID(ePathArr, eID2, mp[1][px], idx, idx+1+px)
		if mp0TeID1*2 < mp0TeID2 {
			h1 = eID2
		} else if mp0TeID1 > mp0TeID2*2 {
			h1 = eID1
		} else {
			ok = false
		}
		if mp1TeID1*2 < mp1TeID2 {
			h2 = eID2
		} else if mp1TeID1 > mp1TeID2*2 {
			h2 = eID1
		} else {
			ok = false
		}
		if h1 == h2 {
			ok = false
		}
	}

	return
}

func DiffLenFirst(a1, a2 []uint32, diffLen int) int {
	if len(a1) != len(a2) {
		fmt.Printf("[DiffLenFirst] len(a1): %v != len(a2): %v\n", len(a1), len(a2))
	}
	min := MinInt(len(a1), len(a2))
	//max := MaxInt(len(a1), len(a2))
	i := 0
	for ; i < min; i++ {
		if a1[i] != a2[i] {
			break
		}
	}
	if i == min {
		return 0
	} else {
		return i + 1
	}
}

func DiffLenLast(a1, a2 []uint32) int {
	if len(a1) != len(a2) {
		fmt.Printf("[DiffLenLast] len(a1): %v != len(a2): %v\n", len(a1), len(a2))
	}
	min := MinInt(len(a1), len(a2))
	//max := MaxInt(len(a1), len(a2))
	i := min - 1
	for ; i >= 0; i-- {
		if a1[i] != a2[i] {
			break
		}
	}
	if i < 0 {
		return 0
	} else {
		return min - i
	}
}

func CheckPathConnectity(path Path, edgesArr []DBGEdge, nodesArr []DBGNode) Path {
	for i := 1; i < len(path.IDArr)-1; i++ {
		id := path.IDArr[i]
		e := edgesArr[id]
		le := edgesArr[path.IDArr[i-1]]
		ne := edgesArr[path.IDArr[i+1]]
		if ((e.StartNID == le.StartNID) || (e.StartNID == le.EndNID)) && ((e.EndNID == ne.StartNID) || (e.EndNID == ne.EndNID)) {
			continue
		} else if (e.EndNID == le.StartNID || e.EndNID == le.EndNID) && (e.StartNID == ne.StartNID || e.StartNID == ne.EndNID) {
			continue
		}

		if path.AltArr[i] > 0 {
			path.IDArr[i], path.AltArr[i] = path.AltArr[i], path.IDArr[i]
			fmt.Printf("[CheckPathConnectity] changed IDArrr[%d]: %d to AltArr[%d]: %d in path: %v\n", i, path.IDArr[i], i, path.AltArr[i], path.IDArr)
			continue
		}

		var nd DBGNode
		if (e.StartNID == le.StartNID) || (e.StartNID == le.EndNID) {
			nd = nodesArr[e.StartNID]
		} else if (e.EndNID == le.StartNID) || (e.EndNID == le.EndNID) {
			nd = nodesArr[e.EndNID]
		}
		if nd.ID > 0 {
			eIDArr := GetOtherEArr(nd, e.ID)
			for _, id := range eIDArr {
				te := edgesArr[id]
				if te.StartNID == nd.ID {
					if te.EndNID == ne.StartNID || te.EndNID == ne.EndNID {
						path.AltArr[i] = path.IDArr[i]
						path.IDArr[i] = te.ID
						fmt.Printf("[CheckPathConnectity] changed IDArrr[%d]: %d to EID: %d in path: %v\n", i, path.IDArr[i], te.ID, path.IDArr)
						break
					}
				} else {
					if te.StartNID == ne.StartNID || te.StartNID == ne.EndNID {
						path.AltArr[i] = path.IDArr[i]
						path.IDArr[i] = te.ID
						fmt.Printf("[CheckPathConnectity] changed IDArrr[%d]: %d to EID: %d in path: %v\n", i, path.IDArr[i], te.ID, path.IDArr)
						break
					}
				}
			}
		}

		nd.ID = 0
		if (e.StartNID == ne.StartNID) || (e.StartNID == ne.EndNID) {
			nd = nodesArr[e.StartNID]
		} else if (e.EndNID == ne.StartNID) || (e.EndNID == ne.EndNID) {
			nd = nodesArr[e.EndNID]
		}
		if nd.ID > 0 {
			eIDArr := GetOtherEArr(nd, e.ID)
			for _, id := range eIDArr {
				te := edgesArr[id]
				if te.StartNID == nd.ID {
					if te.EndNID == le.StartNID || te.EndNID == le.EndNID {
						path.AltArr[i] = path.IDArr[i]
						path.IDArr[i] = te.ID
						fmt.Printf("[CheckPathConnectity] changed IDArrr[%d]: %d to EID: %d in path: %v\n", i, path.IDArr[i], te.ID, path.IDArr)
						break
					}
				} else {
					if te.StartNID == le.StartNID || te.StartNID == le.EndNID {
						path.AltArr[i] = path.IDArr[i]
						path.IDArr[i] = te.ID
						fmt.Printf("[CheckPathConnectity] changed IDArrr[%d]: %d to EID: %d in path: %v\n", i, path.IDArr[i], te.ID, path.IDArr)
						break
					}
				}
			}
		}
	}

	return path
}

func IsBubblePath(leID uint32, idFreqArr []IDFreq, edgesArr []DBGEdge, nodesArr []DBGNode, kmerlen int) (bubOk bool) {
	MaxPathLen := 10
	MaxBubbleSeqLen := kmerlen * 3
	var nd DBGNode
	le := &edgesArr[leID]
	e := &edgesArr[idFreqArr[0].ID]
	if le.StartNID == e.StartNID || le.StartNID == e.EndNID {
		nd = nodesArr[le.StartNID]
	} else {
		nd = nodesArr[le.EndNID]
	}
	pathArr := GetNextPathArr(le, nd, MaxPathLen+1, MaxBubbleSeqLen, edgesArr, nodesArr, kmerlen)
	/*for x, p := range pathArr {
		fmt.Printf("[IsBubblePath]pathArr[%v]: %v\n", x, p)
	}*/
	bubPathArr, _ := GetBubblePathArr(pathArr, edgesArr, kmerlen, MaxPathLen, MaxBubbleSeqLen)
	//fmt.Printf("[IsBubblePath]bubPathArr: %v\n", bubPathArr)
	var bubIDArr []uint32
	for _, bp := range bubPathArr {
		bubIDArr = append(bubIDArr, bp.Path[0])
	}
	bubOk = true
	for _, idf := range idFreqArr {
		if IsInUint32Arr(bubIDArr, idf.ID) {

		} else {
			bubOk = false
			break
		}
	}
	fmt.Printf("[IsBubblePath]leID: %v, idFreqArr: %v, bubOk: %v, bubPathArr: %v\n", leID, idFreqArr, bubOk, bubPathArr)
	return
}

func GetMaxEPathArr(ePathArr [][2][]uint32, eIDLen EIDLen, bubbleRepeatArr []uint32) (leftMax, rightMax int) {
	for _, extArr := range ePathArr {
		idx := IndexPathEID(extArr, uint32(eIDLen.EID))
		if idx > leftMax {
			leftMax = idx
		}
		if len(bubbleRepeatArr) > 1 {
			idx2 := IndexPathEID(extArr, uint32(eIDLen.NextEID))
			if len(extArr[0])-idx2 > rightMax {
				rightMax = len(extArr[0]) - idx2
			}
		} else {
			if len(extArr[0])-idx > rightMax {
				rightMax = len(extArr[0]) - idx
			}
		}
	}
	return
}

func ReAdjustEPathArr(ePathArr [][2][]uint32, eIDLen EIDLen, bubbleRepeatArr []uint32, edgesArr []DBGEdge, leftMax, al int) {
	for j, extArr := range ePathArr {
		a1 := make([]uint32, al)
		a2 := make([]uint32, al)
		idx := IndexPathEID(extArr, uint32(eIDLen.EID))
		if idx >= 0 {
			copy(a1[leftMax-idx:], extArr[0])
			copy(a2[leftMax-idx:], extArr[1])
			ePathArr[j][0] = a1
			ePathArr[j][1] = a2
			continue
		}

		idx = IndexPathEID(extArr, uint32(eIDLen.NextEID))
		if idx >= 0 {
			pos := leftMax + 2*(len(bubbleRepeatArr)-1)
			copy(a1[pos-idx:], extArr[0])
			copy(a2[pos-idx:], extArr[1])
			ePathArr[j][0] = a1
			ePathArr[j][1] = a2
			continue
		}

		for x, id := range extArr[0] {
			if edgesArr[id].GetBubbleRepeatFlag() > 0 {
				idx = IndexUint32(bubbleRepeatArr, id)
				pos := leftMax + 2*idx
				if pos-x < 0 {
					fmt.Printf("[ReAdjustEPathArr]extArr[0]:%v\n\t\tx:%d id:%d bubblereapeatArr:%v idx:%d pos:%d\n", extArr[0], x, id, bubbleRepeatArr, idx, pos)
				}
				copy(a1[pos-x:], extArr[0])
				copy(a2[pos-x:], extArr[1])
				ePathArr[j][0] = a1
				ePathArr[j][1] = a2
				break
			}
		}
	}
}

type NextInfo struct {
	ID     uint16
	Weight uint16
}

type DAGEdge struct {
	EIDArr     []uint32
	Weight     uint16
	NextArr    []NextInfo
	InEdgesID  []uint16
	OutEdgesID []uint16
}

func GetPathEndIdx(pathArr [2][]uint32, bubbleRepeatArr []uint32, direction uint8, j int) (int, int) {
	var end int
	var bub bool
	if direction == FORWARD {
		end = len(pathArr[0])
		if pathArr[1][j] > 0 {
			bub = true
		} else {
			idx := IndexUint32(bubbleRepeatArr, pathArr[0][j])
			if idx >= 0 {
				return j + 1, idx
			}
		}
		for x := j + 1; x < end; x++ {
			eID := pathArr[0][x]
			idx := IndexUint32(bubbleRepeatArr, eID)
			if idx >= 0 {
				return x, -1
			}
			if bub {
				if pathArr[1][x] == 0 {
					return x, -1
				}
			} else {
				if pathArr[1][x] > 0 {
					return x, -1
				}
			}
		}
	} else {
		end = 0
		if pathArr[1][j] > 0 {
			bub = true
		} else {
			idx := IndexUint32(bubbleRepeatArr, pathArr[0][j])
			if idx >= 0 {
				return j, idx
			}
		}
		for x := j; x >= end; x-- {
			eID := pathArr[0][x]
			idx := IndexUint32(bubbleRepeatArr, eID)
			if idx >= 0 {
				return x, idx
			}
			if bub {
				if pathArr[1][x] == 0 {
					return x + 1, -1
				}
			} else {
				if pathArr[1][x] > 0 {
					return x + 1, -1
				}
			}
		}
	}

	return end, -1
}

func AddOutID(ge *DAGEdge, geID uint16) bool {
	var ok bool
	for _, id := range ge.OutEdgesID {
		if id == geID {
			ok = true
			break
		}
	}
	if !ok {
		ge.OutEdgesID = append(ge.OutEdgesID, geID)
	}
	return ok
}

func AddInID(ge *DAGEdge, geID uint16) bool {
	var ok bool
	for _, id := range ge.InEdgesID {
		if id == geID {
			ok = true
			break
		}
	}
	if !ok {
		ge.InEdgesID = append(ge.InEdgesID, geID)
	}
	return ok
}

func AddNextInfo(ge *DAGEdge, geID uint16) bool {
	var ok bool
	for i, ni := range ge.NextArr {
		if ni.ID == geID {
			ok = true
			ge.NextArr[i].Weight++
			break
		}
	}
	if !ok {
		var ni NextInfo
		ni.ID = geID
		ni.Weight = 1
		ge.NextArr = append(ge.NextArr, ni)
	}
	return ok
}

func SplitEIDArrAddGraph(graph []DAGEdge, id uint16, splitLen int, direction uint8) []DAGEdge {
	if len(graph[id].EIDArr) == splitLen {
		return graph
	}
	var te DAGEdge
	if direction == FORWARD {
		te.EIDArr = make([]uint32, len(graph[id].EIDArr)-splitLen)
		te.Weight = graph[id].Weight
		te.OutEdgesID = graph[id].OutEdgesID
		copy(te.EIDArr, graph[id].EIDArr[splitLen:])
		AddInID(&te, id)
		graph[id].EIDArr = graph[id].EIDArr[:splitLen]
		graph = append(graph, te)
		var n []uint16
		graph[id].OutEdgesID = n
		AddOutID(&graph[id], uint16(len(graph)-1))
	} else {
		te.EIDArr = make([]uint32, splitLen)
		te.Weight = graph[id].Weight
		te.InEdgesID = graph[id].InEdgesID
		copy(te.EIDArr, graph[id].EIDArr[:splitLen])
		AddOutID(&te, id)
		graph[id].EIDArr = graph[id].EIDArr[splitLen:]
		graph = append(graph, te)
		var n []uint16
		graph[id].InEdgesID = n
		AddInID(&graph[id], uint16(len(graph)-1))
	}
	return graph
}

func GetExtPath(eIDArr []uint32) []uint32 {
	var arr []uint32
	for _, id := range eIDArr {
		if id > 1 {
			arr = append(arr, id)
		} else if id == 0 {
			log.Fatalf("[GetExtPath]eIDArr:%v contain 0\n", eIDArr)
		}
	}

	return arr
}

func GetNextNextGEID(graph []DAGEdge, geID uint16, extArr [2][]uint32, j int, direction uint8) (int, int) {
	nID, nnID := -1, -1
	ge := &graph[geID]
	if direction == FORWARD {
		//j += len(ge.EIDArr)
		for i := 0; i < 2 && j < len(extArr[0]); i++ {
			for _, id := range ge.OutEdgesID {
				if graph[id].EIDArr[0] == extArr[0][j] {
					ge = &graph[id]
					if i == 0 {
						nID = int(id)
					} else {
						nnID = int(id)
					}
					break
				}
			}
			if j+len(ge.EIDArr) <= len(extArr[0]) && !EqualUint32Arr(ge.EIDArr, extArr[0][j:j+len(ge.EIDArr)]) {
				log.Fatalf("[GetNextNextGEID]ge.EIDArr:%v extArr[0][j:j+len(ge.EIDArr)]:%v\n", ge.EIDArr, extArr[0][j:j+len(ge.EIDArr)])
			}
			j += len(ge.EIDArr)
		}
	} else {
		//j--
		for i := 0; i < 2 && j >= 0; i++ {
			for _, id := range ge.InEdgesID {
				if graph[id].EIDArr[len(graph[id].EIDArr)-1] == extArr[0][j] {
					ge = &graph[id]
					if i == 0 {
						nID = int(id)
					} else {
						nnID = int(id)
					}
					break
				}
			}
			if j+1-len(ge.EIDArr) >= 0 && !EqualUint32Arr(ge.EIDArr, extArr[0][j+1-len(ge.EIDArr):j+1]) {
				log.Fatalf("[GetNextNextGEID]ge.EIDArr:%v extArr[0][j+1-len(ge.EIDArr):j+1]:%v\n", ge.EIDArr, extArr[0][j+1-len(ge.EIDArr):j+1])
			}
			j -= len(ge.EIDArr)
		}
	}

	return nID, nnID
}

func FoundFirstGeID(graph []DAGEdge, gID int, extArr [2][]uint32, j int) (geID int) {
	geID = gID
	if j < 0 {
		return geID
	}
	for x := j; x >= 0; {
		nID, nnID := GetNextNextGEID(graph, uint16(gID), extArr, x, BACKWARD)
		if nID >= 0 {
			geID = nID
			gID = nID
			x -= len(graph[nID].EIDArr)
		}
		if nnID >= 0 {
			geID = nnID
			gID = nnID
			x -= len(graph[nnID].EIDArr)
		}
	}
	return geID
}

func CheckINOUT(graph []DAGEdge) {
	for i := range graph {
		ge := &graph[i]
		for _, id := range ge.InEdgesID {
			AddOutID(&graph[id], uint16(i))
		}
		for _, id := range ge.OutEdgesID {
			AddInID(&graph[id], uint16(i))
		}
	}
}

func GetPathArr(graph []DAGEdge, geID uint16, direction uint8) (path []uint32) {
	if direction == FORWARD {
		path = append(path, graph[geID].EIDArr...)
		if len(graph[geID].OutEdgesID) == 0 {
			return
		}
		ngeID := graph[geID].OutEdgesID[0]
		path = append(path, graph[ngeID].EIDArr...)
		for {
			fmt.Printf("geID:%d\n", geID)
			ge := &graph[geID]
			geID = ngeID
			if len(ge.NextArr) == 0 {
				break
			} else if len(ge.NextArr) == 1 {
				ngeID = ge.NextArr[0].ID
				path = append(path, graph[ngeID].EIDArr...)
			} else {
				maxW := uint16(0)
				maxID := uint16(0)
				for _, ni := range ge.NextArr {
					if ni.Weight > maxW {
						maxW = ni.Weight
						maxID = ni.ID
					}
				}
				ngeID = maxID
				path = append(path, graph[maxID].EIDArr...)
			}
		}
	} else {
		if len(graph[geID].InEdgesID) == 0 {
			return
		}
		ngeID := graph[geID].InEdgesID[0]
		path = append(path, graph[ngeID].EIDArr...)
		for {
			ge := &graph[geID]
			fmt.Printf("geID:%d\n", geID)
			geID = ngeID
			if len(ge.NextArr) == 0 {
				break
			} else if len(ge.NextArr) == 1 {
				ngeID = ge.NextArr[0].ID
				path = append(path, graph[ngeID].EIDArr...)
			} else {
				maxW := uint16(0)
				maxID := uint16(0)
				for _, ni := range ge.NextArr {
					if ni.Weight > maxW {
						maxW = ni.Weight
						maxID = ni.ID
					}
				}
				ngeID = maxID
				path = append(path, graph[maxID].EIDArr...)
			}
		}
	}

	return
}

func GetBubbleArr(edgesArr []DBGEdge, eID1, eID2 uint32) (arr []uint32) {
	e1, e2 := &edgesArr[eID1], &edgesArr[eID2]
	for _, id := range e1.EdgeIDOutcoming {
		if id < 2 {
			continue
		}
		if IsInComing(e2.EdgeIDIncoming, id) || IsInComing(e2.EdgeIDOutcoming, id) {
			arr = append(arr, uint32(id))
		}
	}
	if len(arr) > 0 {
		return
	}

	for _, id := range e1.EdgeIDIncoming {
		if id < 2 {
			continue
		}
		if IsInComing(e2.EdgeIDIncoming, id) || IsInComing(e2.EdgeIDOutcoming, id) {
			arr = append(arr, uint32(id))
		}
	}

	return
}

func MergePathArr(edgesArr []DBGEdge, nodesArr []DBGNode, readPathArr []ReadPath, sortEIDLenArr []EIDLen, edgesPathRelationArr [][]uint32, opt optionsDDBG, MaxPathLen, MaxBubbleSeqLen int) {
	for _, eIDLen := range sortEIDLenArr {
		e := &edgesArr[eIDLen.EID]
		//el := len(e.Ks)
		//if IsBubbleEdge(e, nodesArr) && len(e.Ks) < kmerlen*2+5 {
		//	continue
		//}

		// check is a unique edge in the genome
		//Unique := true
		freq := len(edgesPathRelationArr[eIDLen.EID])
		if freq > opt.AvgDepth*2 { // (el / avgReadLen)*Depth = one long read cov this edge probility
			continue
		}

		// merge process
		var bubbleRepeatArr []uint32
		if e.GetBubbleRepeatFlag() > 0 && eIDLen.NextEID != eIDLen.EID {
			bubbleRepeatArr = GetBubbleRepeatArr(e, edgesArr)
		} else {
			bubbleRepeatArr = append(bubbleRepeatArr, uint32(e.ID))
		}
		//ePathArr := CopyReadPathArr(readPathArr, edgesPathRelationArr[e.ID])
		ePathArr := GetBubbleRepeatReadPathArr(bubbleRepeatArr, edgesPathRelationArr, readPathArr)
		ePathArr = checkEdgePathArrDirection(edgesArr, ePathArr, bubbleRepeatArr)

		if DebugModel {
			fmt.Printf("[MergePathArr]eIDLen:%v,bubbleRepeatArr:%v\n", eIDLen, bubbleRepeatArr)
			for x, p := range ePathArr {
				fmt.Printf("ePathArr[%d]:%v\n", x, p)
			}
		}

		// Directed Acyclic Graph and Topological order
		{
			//eID := uint32(eIDLen.EID)
			graph := make([]DAGEdge, len(bubbleRepeatArr), len(bubbleRepeatArr)*3+8)
			for i := range graph {
				graph[i].EIDArr = []uint32{bubbleRepeatArr[i]}
			}
			for _, extArr := range ePathArr {
				idxP, idxB := IndexBubbleRepeatArr(extArr, bubbleRepeatArr)
				if idxB < 0 {
					log.Fatalf("[MergePathArr]pathArr:%v not share bubbleRepeat in arr:%v\n", extArr, bubbleRepeatArr)
				}
				//fmt.Printf("[MergePathArr]idxP:%d idxB:%d ePathArr[%d]:%v\n", idxP, idxB, i, extArr)
				// FORWARD
				{
					geID := uint16(idxB)
					for j := idxP; j < len(extArr[0]); {
						ge := &graph[geID]
						if extArr[0][j] != ge.EIDArr[0] {
							log.Fatalf("[MergePathArr]ge:%v extArr[0]:%v extArr[1]:%v\n", *ge, extArr[0][j:], extArr[1][j:])
						}
						weight := uint16(1)
						if extArr[1][j] > 0 {
							weight = 0
						}
						//fmt.Printf("j:%d\n", j)
						//idx := IndexUint32(bubbleRepeatArr, extArr[0][j])
						//if idx >= 0 {
						//	graph[idx].Weight++
						//	j++
						//}
						end, idx := GetPathEndIdx(extArr, bubbleRepeatArr, FORWARD, j)
						if end-j >= len(ge.EIDArr) {
							if EqualUint32Arr(extArr[0][j:j+len(ge.EIDArr)], ge.EIDArr) {
								j += len(ge.EIDArr)
							} else {
								splitLen := len(ge.EIDArr)
								for x, id := range ge.EIDArr {
									if id != extArr[0][j+x] {
										splitLen = x
										break
									}
								}
								graph = SplitEIDArrAddGraph(graph, geID, splitLen, FORWARD)
								ge = &graph[geID]
								j += splitLen
							}
						} else {
							splitLen := end - j
							for x := j; x < end; x++ {
								if extArr[0][x] != ge.EIDArr[x-j] {
									splitLen = x - j
									break
								}
							}
							graph = SplitEIDArrAddGraph(graph, geID, splitLen, FORWARD)
							ge = &graph[geID]
							j += splitLen
							//log.Fatalf("[MergePathArr]j:%d extArr[0]:%v ge.EIDArr:%v\n", j, extArr[0][j:end], ge.EIDArr)
						}
						ge.Weight += weight
						if j >= len(extArr[0]) {
							break
						}
						if j < idxP {
							break
						}
						end, idx = GetPathEndIdx(extArr, bubbleRepeatArr, FORWARD, j)
						//fmt.Printf("j:%d end:%d\n", j, end)
						if idx >= 0 {
							geID = uint16(idx)
							AddOutID(ge, geID)
						} else {
							eIDArr := make([]uint32, end-j)
							copy(eIDArr, extArr[0][j:end])
							var added bool
							for _, id := range ge.OutEdgesID {
								if graph[id].EIDArr[0] != eIDArr[0] {
									continue
								}
								added = true
								geID = id
								break
							}
							if !added {
								var te DAGEdge
								te.EIDArr = eIDArr
								AddInID(&te, geID)
								geID = uint16(len(graph))
								AddOutID(ge, uint16(len(graph)))
								graph = append(graph, te)
								//fmt.Printf("graph[%d]:%v\n", len(graph)-1, te)
							}
						}
					}
				}

				// BACKWARD
				{
					geID := uint16(idxB)
					for j := idxP + len(graph[geID].EIDArr) - 1; j >= 0; {
						ge := &graph[geID]
						if extArr[0][j] != ge.EIDArr[len(ge.EIDArr)-1] {
							log.Fatalf("[MergePathArr]ge:%v extArr[0]:%v extArr[1]:%v\n", *ge, extArr[0][:j+1], extArr[1][:j+1])
						}
						weight := uint16(1)
						if extArr[1][j] > 0 {
							weight = 0
						}
						//idx := IndexUint32(bubbleRepeatArr, extArr[0][j])
						//if idx >= 0 {
						//	graph[idx].Weight++
						//	j++
						//}
						end, idx := GetPathEndIdx(extArr, bubbleRepeatArr, BACKWARD, j)
						//fmt.Printf("j:%d ge:%v\n", j, *ge)
						if j+1-end >= len(ge.EIDArr) {
							if EqualUint32Arr(extArr[0][j+1-len(ge.EIDArr):j+1], ge.EIDArr) {
								j -= len(ge.EIDArr)
							} else {
								splitLen := len(ge.EIDArr)
								for x := len(ge.EIDArr) - 1; x >= 0; x-- {
									y := j - (len(ge.EIDArr) - 1 - x)
									if ge.EIDArr[x] != extArr[0][y] {
										splitLen = x + 1
										break
									}
								}
								graph = SplitEIDArrAddGraph(graph, geID, splitLen, BACKWARD)
								ge = &graph[geID]
								j -= len(ge.EIDArr)
							}
						} else {
							splitLen := len(ge.EIDArr) - (j + 1 - end)
							for x := j; x >= end; x-- {
								y := len(ge.EIDArr) - 1 - (j - x)
								if extArr[0][x] != ge.EIDArr[y] {
									splitLen = y + 1
									break
								}
							}
							graph = SplitEIDArrAddGraph(graph, geID, splitLen, BACKWARD)
							ge = &graph[geID]
							j -= len(ge.EIDArr)
							//log.Fatalf("[MergePathArr]j:%d extArr[0]:%v ge.EIDArr:%v\n", j, extArr[0][j:end], ge.EIDArr)
						}
						ge.Weight += weight
						if j < 0 {
							break
						}
						if j > idxP {
							break
						}
						end, idx = GetPathEndIdx(extArr, bubbleRepeatArr, BACKWARD, j)
						//fmt.Printf("end:%d j:%d ge:%v\n", end, j, *ge)
						if idx >= 0 {
							geID = uint16(idx)
							AddInID(ge, geID)
						} else {
							eIDArr := make([]uint32, j+1-end)
							copy(eIDArr, extArr[0][end:j+1])
							var added bool
							for _, id := range ge.InEdgesID {
								tge := &graph[id]
								if tge.EIDArr[len(tge.EIDArr)-1] != eIDArr[len(eIDArr)-1] {
									continue
								}
								added = true
								geID = id
								break
							}
							if !added {
								var te DAGEdge
								te.EIDArr = eIDArr
								AddOutID(&te, geID)
								geID = uint16(len(graph))
								AddInID(ge, uint16(len(graph)))
								graph = append(graph, te)
								//fmt.Printf("added graph[%d]:%v\n", len(graph)-1, te)
							}
						}
					}
				}
				CheckINOUT(graph)
			}

			// set NextInfo
			for _, extArr := range ePathArr {
				//fmt.Printf("extArr:%v\n", extArr)
				idxP, idxB := IndexBubbleRepeatArr(extArr, bubbleRepeatArr)
				geID := -1
				start := -1
				if idxB < 0 {
					log.Fatalf("[MergePathArr]pathArr:%v not share bubbleRepeat in arr:%v\n", extArr, bubbleRepeatArr)
				} else if idxB > 0 {
					geID = FoundFirstGeID(graph, idxB, extArr, idxP-1)
					start = 0
				} else {
					geID = idxB
					start = idxP
					if e.GetBubbleRepeatFlag() > 0 {
						start += len(graph[geID].EIDArr)
						geID, _ = GetNextNextGEID(graph, uint16(geID), extArr, start, FORWARD)
					}
				}
				// FORWARD
				gID := geID
				for j := start; j < len(extArr[0]); {
					ge := &graph[gID]
					if !EqualUint32Arr(ge.EIDArr, extArr[0][j:j+len(ge.EIDArr)]) {
						log.Fatalf("[MergePathArr]ge.EIDArr:%v, extArr[0][j:j+len(ge.EIDArr)]:%v\n", ge.EIDArr, extArr[0][j:j+len(ge.EIDArr)])
					}
					j += len(ge.EIDArr)
					nID, nnID := GetNextNextGEID(graph, uint16(gID), extArr, j, FORWARD)
					//fmt.Printf("j:%d nID:%d nnID:%d\n", j, nID, nnID)
					if nnID >= 0 {
						AddNextInfo(ge, uint16(nnID))
					} else {
						break
					}
					gID = nID
				}

				// BACKWARD
				if idxB == 0 {
					gID := geID
					for j := start + len(graph[gID].EIDArr) - 1; j >= 0; {
						ge := &graph[gID]
						if !EqualUint32Arr(ge.EIDArr, extArr[0][j+1-len(ge.EIDArr):j+1]) {
							log.Fatalf("[MergePathArr]ge.EIDArr:%v, extArr[0][j+1-len(ge.EIDArr):j+1]:%v\n", ge.EIDArr, extArr[0][j+1-len(ge.EIDArr):j+1])
						}
						j -= len(ge.EIDArr)
						nID, nnID := GetNextNextGEID(graph, uint16(gID), extArr, j, BACKWARD)
						//fmt.Printf("j:%d nID:%d nnID:%d\n", j, nID, nnID)
						if nnID >= 0 {
							AddNextInfo(ge, uint16(nnID))
						} else {
							break
						}
						gID = nID
					}
				}
			}
			// check code
			if DebugModel {
				for i, ge := range graph {
					fmt.Printf("graph[%d]:%v\n", i, ge)
				}
			}

			// get path
			geID := -1
			for i, ge := range graph {
				if ge.Weight > 8 && len(ge.NextArr) > 0 {
					geID = i
					break
				}
			}
			//if geID < 0 {
			//	continue
			//}
			fmt.Printf("seed geID:%d\n", geID)
			var p [2][]uint32
			if geID == 0 {
				if e.GetUniqueFlag() > 0 {
					path := GetPathArr(graph, uint16(geID), BACKWARD)
					ReverseUint32Arr(path)
					path = append(path, GetPathArr(graph, uint16(geID), FORWARD)...)
					p[0] = Transform2Path(path)
				} else if e.GetBubbleRepeatFlag() > 0 {
					for i, gID := range graph[geID].OutEdgesID {
						if i > 1 {
							break
						}
						path := GetPathArr(graph, gID, BACKWARD)
						ReverseUint32Arr(path)
						path = append(path, GetPathArr(graph, gID, FORWARD)...)
						p[i] = Transform2Path(path)
					}
				}
			} else if geID > 0 {
				for i, gID := range graph[geID].OutEdgesID {
					if i > 1 {
						break
					}
					path := GetPathArr(graph, gID, FORWARD)
					p[i] = Transform2Path(path)
				}
			}

			for i := 0; i < 2; i++ {
				var tp []uint32
				for j, id := range bubbleRepeatArr {
					tp = append(tp, id)
					if j < len(bubbleRepeatArr)-1 {
						arr := GetBubbleArr(edgesArr, id, bubbleRepeatArr[j+1])
						if len(arr) == 0 {
							break
						}
						tp = append(tp, arr[time.Now().Nanosecond()%len(arr)])
					}
				}
				fmt.Printf("tp:%v\n", tp)
				var pm Path
				pm.IDArr = Transform2Path(tp)
				for _, eID := range tp {
					e := &edgesArr[eID]
					if e.GetBubbleRepeatFlag() > 0 || e.GetUniqueFlag() > 0 {
						e.PathMat = append(e.PathMat, pm)
					}
				}
			}

			/*for i := range p {
				fmt.Printf("path[%d]:%v\n", i, p[i])
				var pm Path
				pm.IDArr = p[i]
				for _, eID := range p[i] {
					e := &edgesArr[eID]
					if e.GetBubbleRepeatFlag() > 0 || e.GetUniqueFlag() > 0 {
						e.PathMat = append(e.PathMat, pm)
					}
				}
			}*/
		}

		/*leftMax, rightMax := GetMaxEPathArr(ePathArr, eIDLen, bubbleRepeatArr)
		var ne DBGEdge
		if len(bubbleRepeatArr) > 1 {
			ne = edgesArr[eIDLen.NextEID]
		}

		var al int
		al = leftMax + 2*(len(bubbleRepeatArr)-1) + rightMax
		fmt.Printf("[MergePathArr]leftMax:%d rightMax:%d al:%d\n", leftMax, rightMax, al)
		// adjust length of ePathArr
		ReAdjustEPathArr(ePathArr, eIDLen, bubbleRepeatArr, edgesArr, leftMax, al)

		fmt.Printf("[MergePathArr]merge eID:%d\n", e.ID)
		if DebugModel {
			for x, p := range ePathArr {
				fmt.Printf("ePathArr[%d]:%v\n", x, p)
			}
		}

		// find consis path, allow many consis path, if Freq > minMapFreq
		//var cleanPathNum int
		//var path Path
		//path.Freq = math.MaxInt64

		var mp [2][]uint32
		var altmp [2][]uint32
		MinDeleteFreqDiff := 5 // delete low freq eID if is MinDeleteFreqDiff*freq[i] < freq[i-1]
		UniqueScore := 5       // if path just one edge, Path.AltArr is nil, path is unique, not ambiguous
		freq0, freq1 := int32(math.MaxInt32), int32(math.MaxInt32)
		// FORWARD
		{
			mp[0] = append(mp[0], e.ID)
			if opt.Haplotype {
				mp[1] = append(mp[1], e.ID)
			}
			//mp0Idx, mp1Idx := len(mp[0]), len(mp[1])
			//freq := math.MaxInt32
			ok0, ok1 := true, true
			for j := leftMax + 1; j < al; j++ {
				diffLen := 0
				if opt.Haplotype {
					diffLen = DiffLenLast(mp[0], mp[1])
				}
				var idFreqArr0, idFreqArr1 []IDFreq
				if ok0 {
					idFreqArr0 = GetIDFreqArrHaplotype(ePathArr, j, FORWARD, diffLen, mp[0], altmp[0], UniqueScore)
				}
				if diffLen > 0 && ok1 {
					idFreqArr1 = GetIDFreqArrHaplotype(ePathArr, j, FORWARD, diffLen, mp[1], altmp[1], UniqueScore)
				}
				sort.Sort(IDFreqArr(idFreqArr0))
				sort.Sort(IDFreqArr(idFreqArr1))
				num0, num1 := len(idFreqArr0), len(idFreqArr1)
				fmt.Printf("[MergePathArr][%v] idFreqArr0: %v, idFreqArr1: %v\n", j, idFreqArr0, idFreqArr1)
				idFreqArr0 = DeleteLowFreq(idFreqArr0, int32(opt.MinMapFreq*UniqueScore), int32(MinDeleteFreqDiff))
				idFreqArr1 = DeleteLowFreq(idFreqArr1, int32(opt.MinMapFreq*UniqueScore), int32(MinDeleteFreqDiff))

				if len(idFreqArr0) == 0 && len(idFreqArr1) == 0 {
					break
				}

				if len(idFreqArr0) <= 1 && len(idFreqArr1) <= 1 {
					if len(idFreqArr0) == 1 {
						mp[0] = append(mp[0], idFreqArr0[0].ID)
						if freq0 > idFreqArr0[0].Freq {
							freq0 = idFreqArr0[0].Freq
						}
						//AdjustEPathMat2(ePathArr, j, idFreqArr0[0].ID, FORWARD, edgesArr, mp[0], kmerlen, diffLen)
					}
					if len(idFreqArr1) == 1 {
						mp[1] = append(mp[1], idFreqArr1[0].ID)
						if freq1 > idFreqArr1[0].Freq {
							freq1 = idFreqArr1[0].Freq
						}
						//AdjustEPathMat2(ePathArr, j, idFreqArr1[0].ID, FORWARD, edgesArr, mp[1], kmerlen, diffLen)
					} else if opt.Haplotype {
						mp[1] = append(mp[1], idFreqArr0[0].ID)
						if freq1 > idFreqArr0[0].Freq {
							freq1 = idFreqArr0[0].Freq
						}
						//AdjustEPathMat2(ePathArr, j, idFreqArr0[0].ID, FORWARD, edgesArr, mp[1], kmerlen, diffLen)
					}
				} else {
					if diffLen == 0 {
						if idFreqArr0[0].Freq > idFreqArr0[1].Freq*5 {
							mp[0] = append(mp[0], idFreqArr0[0].ID)
							if freq0 > idFreqArr0[0].Freq {
								freq0 = idFreqArr0[0].Freq
							}
							//AdjustEPathMat2(ePathArr, j, idFreqArr0[0].ID, FORWARD, edgesArr, mp[0], kmerlen, diffLen)
							if opt.Haplotype {
								mp[1] = append(mp[1], idFreqArr0[0].ID)
								if freq1 > idFreqArr0[0].Freq {
									freq1 = idFreqArr0[0].Freq
								}
								//AdjustEPathMat2(ePathArr, j, idFreqArr0[0].ID, FORWARD, edgesArr, mp[1], kmerlen, diffLen)
							}
						} else if !opt.Haplotype && idFreqArr0[0].Freq > idFreqArr0[1].Freq*2 {
							mp[0] = append(mp[0], idFreqArr0[0].ID)
							if freq0 > idFreqArr0[0].Freq {
								freq0 = idFreqArr0[0].Freq
							}
							//AdjustEPathMat2(ePathArr, j, idFreqArr0[0].ID, FORWARD, edgesArr, mp[0], kmerlen, diffLen)
						} else { // len(idFreqArr0) > 1, low quality edge
							fmt.Fprintf(os.Stdout, "[MergePathArr]eIDLen: %v, idx: %d has low quality edge: {eID1[%d] : Freq[%d]}, {eID2[%d] : Freq[%d]}\n", eIDLen, j, idFreqArr0[0].ID, idFreqArr0[0].Freq, idFreqArr0[1].ID, idFreqArr0[1].Freq)
							e1, e2 := edgesArr[idFreqArr0[0].ID], edgesArr[idFreqArr0[1].ID]
							if e1.GetTwoEdgesCycleFlag() > 0 || e2.GetTwoEdgesCycleFlag() > 0 || e1.StartNID == e1.EndNID || e2.StartNID == e2.EndNID {
								mp[0] = append(mp[0], idFreqArr0[0].ID)
								if freq0 > idFreqArr0[0].Freq {
									freq0 = idFreqArr0[0].Freq
								}
							} else if opt.Haplotype {
								mp[0] = append(mp[0], idFreqArr0[0].ID)
								if freq0 > idFreqArr0[0].Freq {
									freq0 = idFreqArr0[0].Freq
								}
								mp[1] = append(mp[1], idFreqArr0[1].ID)
								if freq1 > idFreqArr0[1].Freq {
									freq1 = idFreqArr0[1].Freq
								}
								//AdjustEPathMat2(ePathArr, j, idFreqArr1[1].ID, FORWARD, edgesArr, mp[1], kmerlen, diffLen)
							} else if IsBubblePath(mp[0][len(mp[0])-1], idFreqArr0, edgesArr, nodesArr, opt.Kmer) {
								dl := len(mp[0]) - len(altmp[0])
								ta := make([]uint32, dl)
								altmp[0] = append(altmp[0], ta...)
								mp[0] = append(mp[0], idFreqArr0[0].ID)
								if freq0 > idFreqArr0[0].Freq {
									freq0 = idFreqArr0[0].Freq
								}
								altmp[0] = append(altmp[0], idFreqArr0[1].ID)
								//AdjustEPathMat2(ePathArr, j, idFreqArr0[0].ID, FORWARD, edgesArr, mp[0], kmerlen, diffLen)
							} else {
								//ok0 = false
								break
							}
						}
					} else { // diffLen > 0
						if len(idFreqArr0) == 1 {
							mp[0] = append(mp[0], idFreqArr0[0].ID)
							if freq0 > idFreqArr0[0].Freq {
								freq0 = idFreqArr0[0].Freq
							}
							//AdjustEPathMat2(ePathArr, j, idFreqArr0[0].ID, FORWARD, edgesArr, mp[0], kmerlen, diffLen)
						} else if len(idFreqArr0) > 1 {
							fmt.Fprintf(os.Stderr, "[MergePathArr]diffLen: %v, eIDLen: %v, idx: %d has low quality edge: {eID1[%d] : Freq[%d]}, {eID2[%d] : Freq[%d]}\n", diffLen, eIDLen, j, idFreqArr0[0].ID, idFreqArr0[0].Freq, idFreqArr0[1].ID, idFreqArr0[1].Freq)
							if idFreqArr0[0].Freq > 2*idFreqArr0[1].Freq || IsBubblePath(mp[0][len(mp[0])-1], idFreqArr0, edgesArr, nodesArr, opt.Kmer) {
								dl := len(mp[0]) - len(altmp[0])
								ta := make([]uint32, dl)
								altmp[0] = append(altmp[0], ta...)
								altmp[0] = append(altmp[0], idFreqArr0[1].ID)
								mp[0] = append(mp[0], idFreqArr0[0].ID)
								if freq0 > idFreqArr0[0].Freq {
									freq0 = idFreqArr0[0].Freq
								}
							} else {
								e1, e2 := edgesArr[idFreqArr0[0].ID], edgesArr[idFreqArr0[1].ID]
								if e1.GetTwoEdgesCycleFlag() > 0 || e2.GetTwoEdgesCycleFlag() > 0 || e1.StartNID == e1.EndNID || e2.StartNID == e2.EndNID {
									mp[0] = append(mp[0], idFreqArr0[0].ID)
									if freq0 > idFreqArr0[0].Freq {
										freq0 = idFreqArr0[0].Freq
									}
								} else {
									ok0 = false
								}
							}
						}

						if len(idFreqArr1) == 1 {
							mp[1] = append(mp[1], idFreqArr1[0].ID)
							if freq1 > idFreqArr1[0].Freq {
								freq1 = idFreqArr1[0].Freq
							}
							//AdjustEPathMat2(ePathArr, j, idFreqArr1[0].ID, FORWARD, edgesArr, mp[1], kmerlen, diffLen)
						} else if len(idFreqArr1) > 1 {
							fmt.Fprintf(os.Stderr, "[MergePathArr]diffLen: %v, eIDLen: %v, idx: %d has low quality edge: {eID1[%d] : Freq[%d]}, {eID2[%d] : Freq[%d]}\n", diffLen, eIDLen, j, idFreqArr1[0].ID, idFreqArr1[0].Freq, idFreqArr1[1].ID, idFreqArr1[1].Freq)
							if idFreqArr1[0].Freq > 2*idFreqArr1[1].Freq || IsBubblePath(mp[1][len(mp[1])-1], idFreqArr1, edgesArr, nodesArr, opt.Kmer) {
								dl := len(mp[1]) - len(altmp[1])
								ta := make([]uint32, dl)
								altmp[1] = append(altmp[1], ta...)
								altmp[1] = append(altmp[1], idFreqArr1[1].ID)
								mp[1] = append(mp[1], idFreqArr1[0].ID)
								if freq1 > idFreqArr1[0].Freq {
									freq1 = idFreqArr1[0].Freq
								}
							} else {
								e1, e2 := edgesArr[idFreqArr1[0].ID], edgesArr[idFreqArr1[1].ID]
								if e1.GetTwoEdgesCycleFlag() > 0 || e2.GetTwoEdgesCycleFlag() > 0 || e1.StartNID == e1.EndNID || e2.StartNID == e2.EndNID {
									mp[1] = append(mp[1], idFreqArr0[0].ID)
									if freq1 > idFreqArr1[0].Freq {
										freq1 = idFreqArr1[0].Freq
									}
								} else {
									ok1 = false
								}
							}
						}
					}
				}

				//adjust EPathMat
				{
					if opt.Haplotype {
						dl := DiffLenLast(mp[0], mp[1])
						if dl == 1 {
							continue
						}
					}

					if num0 > 1 {
						AdjustEPathMat2(ePathArr, j, idFreqArr0[:num0], FORWARD, edgesArr, nodesArr, mp[0], opt.Kmer, diffLen, MaxPathLen, MaxBubbleSeqLen)
					}

					if num1 > 1 {
						AdjustEPathMat2(ePathArr, j, idFreqArr1[:num1], FORWARD, edgesArr, nodesArr, mp[1], opt.Kmer, diffLen, MaxPathLen, MaxBubbleSeqLen)
					}
				}
			}
		}

		if len(altmp[0]) < len(mp[0]) {
			dl := len(mp[0]) - len(altmp[0])
			ta := make([]uint32, dl)
			altmp[0] = append(altmp[0], ta...)
		}
		if opt.Haplotype && len(altmp[1]) < len(mp[1]) {
			dl := len(mp[1]) - len(altmp[1])
			ta := make([]uint32, dl)
			altmp[1] = append(altmp[1], ta...)
		}
		// BACKWARD
		{
			forwardDiffLen := AbsInt(len(mp[0]) - len(mp[1]))
			//freq := math.MaxInt32
			ok0, ok1 := true, true
			for j := leftMax - 1; j >= 0; j-- {
				diffLen := 0
				if opt.Haplotype {
					diffLen = DiffLenFirst(mp[0], mp[1], forwardDiffLen)
				}
				GetIDFreqArr(ePathArr, j)
				var idFreqArr0, idFreqArr1 []IDFreq
				if ok0 {
					idFreqArr0 = GetIDFreqArrHaplotype(ePathArr, j, BACKWARD, diffLen, mp[0], altmp[0], UniqueScore)
				}
				if diffLen > 0 && ok1 {
					idFreqArr1 = GetIDFreqArrHaplotype(ePathArr, j, BACKWARD, diffLen, mp[1], altmp[1], UniqueScore)
				}

				sort.Sort(IDFreqArr(idFreqArr0))
				sort.Sort(IDFreqArr(idFreqArr1))
				num0, num1 := len(idFreqArr0), len(idFreqArr1)
				fmt.Printf("[MergePathArr][%v] idFreqArr0: %v, idFreqArr1: %v\n", j, idFreqArr0, idFreqArr1)
				idFreqArr0 = DeleteLowFreq(idFreqArr0, int32(opt.MinMapFreq*UniqueScore), int32(MinDeleteFreqDiff))
				idFreqArr1 = DeleteLowFreq(idFreqArr1, int32(opt.MinMapFreq*UniqueScore), int32(MinDeleteFreqDiff))
				if len(idFreqArr0) == 0 && len(idFreqArr1) == 0 {
					break
				}

				if len(idFreqArr0) <= 1 && len(idFreqArr1) <= 1 {
					var n0, n1 []uint32
					if len(idFreqArr0) == 1 {
						n0 = append(n0, idFreqArr0[0].ID)
						n0 = append(n0, mp[0]...)
						mp[0] = n0
						if freq0 > idFreqArr0[0].Freq {
							freq0 = idFreqArr0[0].Freq
						}
						//AdjustEPathMat2(ePathArr, j, idFreqArr0[0].ID, BACKWARD, edgesArr, mp[0], kmerlen, diffLen)
					}
					if len(idFreqArr1) == 1 {
						n1 = append(n1, idFreqArr1[0].ID)
						n1 = append(n1, mp[1]...)
						mp[1] = n1
						if freq1 > idFreqArr1[0].Freq {
							freq1 = idFreqArr1[0].Freq
						}
						//AdjustEPathMat2(ePathArr, j, idFreqArr1[0].ID, BACKWARD, edgesArr, mp[1], kmerlen, diffLen)
					} else if opt.Haplotype {
						n1 = append(n1, idFreqArr0[0].ID)
						n1 = append(n1, mp[1]...)
						mp[1] = n1
						if freq1 > idFreqArr0[0].Freq {
							freq1 = idFreqArr0[0].Freq
						}
						//AdjustEPathMat2(ePathArr, j, idFreqArr0[0].ID, BACKWARD, edgesArr, mp[1], kmerlen, diffLen)
					}
				} else {
					if diffLen == 0 {
						if idFreqArr0[0].Freq > idFreqArr0[1].Freq*5 {
							var n0, n1 []uint32
							n0 = append(n0, idFreqArr0[0].ID)
							n0 = append(n0, mp[0]...)
							mp[0] = n0
							if freq0 > idFreqArr0[0].Freq {
								freq0 = idFreqArr0[0].Freq
							}
							//AdjustEPathMat2(ePathArr, j, idFreqArr0[0].ID, BACKWARD, edgesArr, mp[0], kmerlen, diffLen)
							if opt.Haplotype {
								n1 = append(n1, idFreqArr0[0].ID)
								n1 = append(n1, mp[1]...)
								mp[1] = n1
								if freq1 > idFreqArr0[0].Freq {
									freq1 = idFreqArr0[0].Freq
								}
								//AdjustEPathMat2(ePathArr, j, idFreqArr0[0].ID, BACKWARD, edgesArr, mp[1], kmerlen, diffLen)
							}
						} else if !opt.Haplotype && idFreqArr0[0].Freq > idFreqArr0[1].Freq*2 {
							var n0 []uint32
							n0 = append(n0, idFreqArr0[0].ID)
							n0 = append(n0, mp[0]...)
							mp[0] = n0
							if freq0 > idFreqArr0[0].Freq {
								freq0 = idFreqArr0[0].Freq
							}
							//AdjustEPathMat2(ePathArr, j, idFreqArr0[0].ID, BACKWARD, edgesArr, mp[0], kmerlen, diffLen)
						} else { // len(idFreqArr0) > 1, low quality edge
							fmt.Fprintf(os.Stderr, "[MergePathArr]eIDLen: %v, idx: %d has low quality edge: {eID1[%d] : Freq[%d]}, {eID2[%d] : Freq[%d]}\n", eIDLen, j, idFreqArr0[0].ID, idFreqArr0[0].Freq, idFreqArr0[1].ID, idFreqArr0[1].Freq)
							e1, e2 := edgesArr[idFreqArr0[0].ID], edgesArr[idFreqArr0[1].ID]
							if e1.GetTwoEdgesCycleFlag() > 0 || e2.GetTwoEdgesCycleFlag() > 0 || e1.StartNID == e1.EndNID || e2.StartNID == e2.EndNID {
								var n0, alt []uint32
								dl := len(mp[0]) - len(altmp[0])
								n0 = append(n0, idFreqArr0[0].ID)
								n0 = append(n0, mp[0]...)
								mp[0] = n0
								alt = append(alt, idFreqArr0[1].ID)
								//AdjustEPathMat2(ePathArr, j, idFreqArr0[0].ID, BACKWARD, edgesArr, mp[0], kmerlen, diffLen)
								ta := make([]uint32, dl)
								alt = append(alt, ta...)
								alt = append(alt, altmp[0]...)
								altmp[0] = alt
							} else if opt.Haplotype {
								var n0, n1 []uint32
								n0 = append(n0, idFreqArr0[0].ID)
								n0 = append(n0, mp[0]...)
								mp[0] = n0
								if freq0 > idFreqArr0[0].Freq {
									freq0 = idFreqArr0[0].Freq
								}
								n1 = append(n1, idFreqArr0[1].ID)
								n1 = append(n1, mp[1]...)
								mp[1] = n1
								if freq1 > idFreqArr0[1].Freq {
									freq1 = idFreqArr0[1].Freq
								}
							} else if IsBubblePath(mp[0][0], idFreqArr0, edgesArr, nodesArr, opt.Kmer) {
								var n0, alt []uint32
								dl := len(mp[0]) - len(altmp[0])
								n0 = append(n0, idFreqArr0[0].ID)
								n0 = append(n0, mp[0]...)
								mp[0] = n0
								alt = append(alt, idFreqArr0[1].ID)
								//AdjustEPathMat2(ePathArr, j, idFreqArr0[0].ID, BACKWARD, edgesArr, mp[0], kmerlen, diffLen)
								ta := make([]uint32, dl)
								alt = append(alt, ta...)
								alt = append(alt, altmp[0]...)
								altmp[0] = alt
							} else {
								break
							}
						}
					} else { // diffLen > 0
						if len(idFreqArr0) == 1 {
							var n0 []uint32
							n0 = append(n0, idFreqArr0[0].ID)
							n0 = append(n0, mp[0]...)
							mp[0] = n0
							if freq0 > idFreqArr0[0].Freq {
								freq0 = idFreqArr0[0].Freq
							}
							//AdjustEPathMat2(ePathArr, j, idFreqArr0[0].ID, BACKWARD, edgesArr, mp[0], kmerlen, diffLen)
						} else if len(idFreqArr0) > 1 {
							fmt.Fprintf(os.Stderr, "[MergePathArr]diffLen: %v, eIDLen: %v, idx: %d has low quality edge: {eID1[%d] : Freq[%d]}, {eID2[%d] : Freq[%d]}\n", diffLen, eIDLen, j, idFreqArr0[0].ID, idFreqArr0[0].Freq, idFreqArr0[1].ID, idFreqArr0[1].Freq)
							var n0, alt []uint32
							if idFreqArr0[0].Freq > 2*idFreqArr0[1].Freq || IsBubblePath(mp[0][0], idFreqArr0, edgesArr, nodesArr, opt.Kmer) {
								dl := len(mp[0]) - len(altmp[0])
								n0 = append(n0, idFreqArr0[0].ID)
								n0 = append(n0, mp[0]...)
								mp[0] = n0
								if freq0 > idFreqArr0[0].Freq {
									freq0 = idFreqArr0[0].Freq
								}
								alt = append(alt, idFreqArr0[1].ID)
								ta := make([]uint32, dl)
								alt = append(alt, ta...)
								altmp[0] = append(alt, altmp[0]...)
								//AdjustEPathMat2(ePathArr, j, idFreqArr0[0].ID, BACKWARD, edgesArr, mp[0], kmerlen, diffLen)
							} else {
								e1, e2 := edgesArr[idFreqArr0[0].ID], edgesArr[idFreqArr0[1].ID]
								if e1.GetTwoEdgesCycleFlag() > 0 || e2.GetTwoEdgesCycleFlag() > 0 || e1.StartNID == e1.EndNID || e2.StartNID == e2.EndNID {
									dl := len(mp[0]) - len(altmp[0])
									n0 = append(n0, idFreqArr0[0].ID)
									n0 = append(n0, mp[0]...)
									mp[0] = n0
									if freq0 > idFreqArr0[0].Freq {
										freq0 = idFreqArr0[0].Freq
									}
									alt = append(alt, idFreqArr0[1].ID)
									ta := make([]uint32, dl)
									alt = append(alt, ta...)
									altmp[0] = append(alt, altmp[0]...)
								} else {
									ok0 = false
								}
							}
						}

						if len(idFreqArr1) == 1 {
							var n1 []uint32
							n1 = append(n1, idFreqArr1[0].ID)
							n1 = append(n1, mp[1]...)
							mp[1] = n1
							if freq1 > idFreqArr1[0].Freq {
								freq1 = idFreqArr1[0].Freq
							}
							//AdjustEPathMat2(ePathArr, j, idFreqArr1[0].ID, BACKWARD, edgesArr, mp[1], kmerlen, diffLen)
						} else if len(idFreqArr1) > 1 {
							fmt.Fprintf(os.Stderr, "[MergePathArr]diffLen: %v, eIDLen: %v, idx: %d has low quality edge: {eID1[%d] : Freq[%d]}, {eID2[%d] : Freq[%d]}\n", diffLen, eIDLen, j, idFreqArr1[0].ID, idFreqArr1[0].Freq, idFreqArr1[1].ID, idFreqArr1[1].Freq)
							var n1, alt []uint32
							if idFreqArr1[0].Freq > 2*idFreqArr1[1].Freq || IsBubblePath(mp[1][0], idFreqArr1, edgesArr, nodesArr, opt.Kmer) {
								dl := len(mp[1]) - len(altmp[1])
								ta := make([]uint32, dl)
								alt = append(alt, idFreqArr1[1].ID)
								alt = append(alt, ta...)
								altmp[1] = append(alt, altmp[1]...)
								n1 = append(n1, idFreqArr1[0].ID)
								n1 = append(n1, mp[1]...)
								mp[1] = n1
								if freq1 > idFreqArr1[0].Freq {
									freq1 = idFreqArr1[0].Freq
								}
								//AdjustEPathMat2(ePathArr, j, idFreqArr1[0].ID, BACKWARD, edgesArr, mp[1], kmerlen, diffLen)
							} else {
								e1, e2 := edgesArr[idFreqArr1[0].ID], edgesArr[idFreqArr1[1].ID]
								if e1.GetTwoEdgesCycleFlag() > 0 || e2.GetTwoEdgesCycleFlag() > 0 || e1.StartNID == e1.EndNID || e2.StartNID == e2.EndNID {
									dl := len(mp[1]) - len(altmp[1])
									n1 = append(n1, idFreqArr1[0].ID)
									n1 = append(n1, mp[1]...)
									mp[1] = n1
									if freq1 > idFreqArr1[0].Freq {
										freq1 = idFreqArr1[0].Freq
									}
									alt = append(alt, idFreqArr1[1].ID)
									ta := make([]uint32, dl)
									alt = append(alt, ta...)
									altmp[1] = append(alt, altmp[1]...)
								} else {
									ok1 = false
								}
							}
						}
					}
				}

				//adjust EPathMat
				{
					if opt.Haplotype {
						dl := DiffLenFirst(mp[0], mp[1], diffLen)
						if dl == 1 {
							continue
						}
					}
					if num0 > 1 {
						AdjustEPathMat2(ePathArr, j, idFreqArr0[:num0], BACKWARD, edgesArr, nodesArr, mp[0], opt.Kmer, diffLen, MaxPathLen, MaxBubbleSeqLen)
					}

					if num1 > 1 {
						AdjustEPathMat2(ePathArr, j, idFreqArr1[:num1], BACKWARD, edgesArr, nodesArr, mp[1], opt.Kmer, diffLen, MaxPathLen, MaxBubbleSeqLen)
					}
				}
			}
		}

		if len(altmp[0]) < len(mp[0]) {
			dl := len(mp[0]) - len(altmp[0])
			ta := make([]uint32, dl)
			altmp[0] = append(ta, altmp[0]...)
		}
		if opt.Haplotype && len(altmp[1]) < len(mp[1]) {
			dl := len(mp[1]) - len(altmp[1])
			ta := make([]uint32, dl)
			altmp[1] = append(ta, altmp[1]...)
		}

		if len(mp[0]) >= 3 {
			var path1, path2 Path
			path1.IDArr = mp[0]
			path1.AltArr = altmp[0]
			path1.Freq = int(freq0)
			path1 = CheckPathConnectity(path1, edgesArr, nodesArr)
			edgesArr[e.ID].PathMat = append(edgesArr[e.ID].PathMat, path1)
			if len(bubbleRepeatArr) > 1 {
				ne := edgesArr[eIDLen.NextEID]
				idx2 := IndexUint32(mp[0], eIDLen.NextEID)
				if idx2 >= 0 {
					edgesArr[ne.ID].PathMat = append(edgesArr[ne.ID].PathMat, path1)
				} else {
					fmt.Printf("[MergePathArr] NextEID: %v not include mp[0]: %v\n", eIDLen.NextEID, mp[0])
				}
			}
			if len(mp[1]) > 0 && reflect.DeepEqual(mp[0], mp[1]) == false {
				path2.IDArr = mp[1]
				path2.AltArr = altmp[1]
				path2.Freq = int(freq1)
				path2 = CheckPathConnectity(path2, edgesArr, nodesArr)
				edgesArr[e.ID].PathMat1 = append(edgesArr[e.ID].PathMat1, path2)
				if len(bubbleRepeatArr) > 1 {
					ne := edgesArr[eIDLen.NextEID]
					idx2 := IndexUint32(mp[1], eIDLen.NextEID)
					if idx2 >= 0 {
						edgesArr[ne.ID].PathMat1 = append(edgesArr[ne.ID].PathMat1, path2)
					} else {
						fmt.Printf("[MergePathArr] NextEID: %v not include mp[1]: %v\n", eIDLen.NextEID, mp[1])
					}
				}
			}
			// clean relationship array
			//edgesPathRelationArr[e.ID] = nil
			edgesArr[e.ID].SetMergedFlag()
			fmt.Printf("[MergePathArr]mergeed successed, edgesArr[%v]PathMat  : %v\n", e.ID, edgesArr[e.ID].PathMat)
			fmt.Printf("[MergePathArr]mergeed successed, edgesArr[%v]PathMat1 : %v\n", e.ID, edgesArr[e.ID].PathMat1)
		} */
	}
	// every position just choose max freq edge
	/*{
		var mp []uint32
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
					if idFreqArr[y].Freq >= int32(opt.MinMapFreq) {
						break
					}
				}
				if y < 0 {
					break
				}

				if y >= 1 {
					if idFreqArr[0].Freq > idFreqArr[1].Freq {
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
						AdjustEPathMat2(ePathArr, j, idFreqArr, BACKWARD, edgesArr,nodesArr, mp, opt.Kmer,0, MaxPathLen, MaxBubbleSeqLen)
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
		mp = ReverseUint32Arr(mp)
		//mp = append(mp, e.ID)
		var leftLen int
		if len(mp) > 1 {
			leftLen = len(mp) - 1
		}
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
					if idFreqArr[y].Freq >= int32(opt.MinMapFreq) {
						break
					}
				}
				if y < 0 {
					break
				}

				if y >= 1 {
					if idFreqArr[0].Freq > idFreqArr[1].Freq {
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
						AdjustEPathMat2(ePathArr, j, idFreqArr, FORWARD, edgesArr, mp[leftLen:], opt.Kmer)
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
			var path1, path2 Path
			path1.IDArr = mp
			path2.IDArr = make([]uint32, len(mp))
			path1.Freq, path2.Freq = freq, freq
			edgesArr[e.ID].PathMat = append(edgesArr[e.ID].PathMat, path1)
			edgesArr[e.ID].PathMat1 = append(edgesArr[e.ID].PathMat1, path2)
			// clean relationship array
			//edgesPathRelationArr[e.ID] = nil
			edgesArr[e.ID].SetMergedFlag()
			fmt.Printf("[MergePathArr]mergeed successed, edgesArr[%v]PathMat : %v, edgesArr[%v]PathMat1: %v\n", e.ID, edgesArr[e.ID].PathMat, e.ID, edgesArr[e.ID].PathMat1)
		}
	}*/

	/*	// process left partition
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
					if freq >= opt.MinMapFreq {
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
								if idFreqArr[y].Freq >= int32(opt.MinMapFreq) {
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
							p.PathArr = CleanLowFreq(idFreqArr[y+1:], j, p.PathArr, BACKWARD)
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
											var pa [2][]uint32
											pa[0] = make([]uint32, al)
											pa[1] = make([]uint32, al)
											copy(pa[0], p.PathArr[m][0])
											copy(pa[1], p.PathArr[m][1])
											pa[1][j] = 0
											pai.PathArr = append(pai.PathArr, pa)
											p.PathArr[m][0][j] = p.PathArr[m][1][j]
											p.PathArr[m][1][j] = 0
										} else if p.PathArr[m][0][j] > 0 && p.PathArr[m][1][j] == id {
											var pa [2][]uint32
											pa[0] = make([]uint32, al)
											pa[1] = make([]uint32, al)
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
								if idFreqArr[y].Freq >= int32(opt.MinMapFreq) {
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
						p.PathArr = CleanLowFreq(idFreqArr[y+1:], j, p.PathArr, FORWARD)

						if y == -1 || y >= 1 { // need stack
							if y == -1 {
								var path1, path2 Path
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
									e.PathMat = append(e.PathMat, path1)
									e.PathMat1 = append(e.PathMat1, path2)
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
											var pa [2][]uint32
											pa[0] = make([]uint32, al)
											pa[1] = make([]uint32, al)
											copy(pa[0], p.PathArr[m][0])
											copy(pa[1], p.PathArr[m][1])
											pa[1][j] = 0
											pai.PathArr = append(pai.PathArr, pa)
											p.PathArr[m][0][j] = p.PathArr[m][1][j]
											p.PathArr[m][1][j] = 0
										} else if p.PathArr[m][0][j] > 0 && p.PathArr[m][1][j] == id {
											var pa [2][]uint32
											pa[0] = make([]uint32, al)
											pa[1] = make([]uint32, al)
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
			for x, p := range e.PathMat {
				fmt.Printf("[MergePathArr]edgesArr[%d]PathMat[%d]:%v\n", e.ID, x, p)
			}
			// if long edge, can merge
			if len(e.PathMat) == 2 {
				mergePath0, mergePath1 := MergeLongEdgePathMat(edgesArr[i], edgesArr, nodesArr)
				if len(mergePath0.IDArr) > 0 {
					e.PathMat[0] = mergePath0
					e.PathMat1[0] = mergePath1
					e.PathMat = e.PathMat[:1]
					e.PathMat1 = e.PathMat1[:1]
					fmt.Printf("[MergePathArr]mergeed successed, edgesArr[%v]PathMat : %v, edgesArr[%v]PathMat1: %v\n", i, edgesArr[i].PathMat, i, edgesArr[i].PathMat1)
				}
			}
		}
		// add left==0 and right==1 partition
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
				//var id1, id2 uint32
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
				path.IDArr = GetReverseUint32Arr(path.IDArr)
			}
		}
		fmt.Printf("[MergePathArr]notConsis: %v, path: %v\n", notConsis, path)

		if !notConsis && len(path.IDArr) >= 2 {
			var pm []Path
			pm = append(pm, path)
			edgesArr[i].PathMat = pm
			fmt.Printf("[MergePathArr]eID: %v,  pm: %v\n", i, pm)
		}*/
}

func GetContigSeq(p Path, edgesArr []DBGEdge, nodesArr []DBGNode, kmerlen int) (seq []byte) {
	var strand bool
	e := &edgesArr[p.IDArr[0]]
	e1 := &edgesArr[p.IDArr[1]]
	if e.StartNID == e.EndNID {
		log.Fatalf("[GetContigSeq] in path: %v, start edge is self cycle: %v\n", p, e)
	}
	if (e.EndNID == e1.StartNID) || (e.EndNID == e1.EndNID) {
		strand = true
	} else {
		strand = false
	}
	for i, id := range p.IDArr {
		e := &edgesArr[id]
		pos := 0
		if i > 0 {
			pos = kmerlen - 1
		}
		if strand {
			seq = append(seq, e.Ks[pos:]...)
		} else {
			seq = append(seq, GetReverseCompByteArr(e.Ks[:len(e.Ks)-pos])...)
		}
		// reset strand
		if i < len(p.IDArr)-1 {
			ne := &edgesArr[p.IDArr[i+1]]
			strand = GetNextEdgeStrand2(e, ne, strand, FORWARD)
		}
	}

	return
}

func ExtractSeq(edgesArr []DBGEdge, nodesArr []DBGNode, joinPathArr []Path, DcDBGEdgesfn string, kmerlen int) {
	fp, err := os.Create(DcDBGEdgesfn)
	if err != nil {
		log.Fatalf("[ExtractSeq] create file: %s failed, err: %v\n", DcDBGEdgesfn, err)
	}
	defer fp.Close()
	buffp := bufio.NewWriterSize(fp, 1<<20)

	//fafp := fasta.NewWriter(buffp, 100)
	seqID := int(1)

	for _, p := range joinPathArr {
		if len(p.IDArr) <= 1 {
			continue
		}
		fmt.Printf("[ReconstructDBG] p: %v\n", p)
		//if IsTwoEdgeCyclePath(path) { joinPathArr[i].IDArr = }

		//e := CascadePath(p, edgesArr, nodesArr, kmerlen, false)
		seq := GetContigSeq(p, edgesArr, nodesArr, kmerlen)

		//contig := linear.NewSeq("", Transform2Letters(seq), alphabet.DNA)
		//contig.ID = strconv.Itoa(int(seqID))
		seqID++
		// Add start and end adapter seq
		//qs := Transform2QSeq(e.
		//seq.AppendQLetters(qs...)
		var ps string
		for _, eID := range p.IDArr {
			ps += strconv.Itoa(int(eID)) + "-"
		}
		ps = ps[:len(ps)-1]
		ans := ">" + strconv.Itoa(int(seqID)) + "\tpath: " + ps + "\tFreq: " + strconv.Itoa(p.Freq) + "\n"
		//contig.Annotation.SetDescription(ans)
		//_, err := fafp.Write(contig)
		buffp.WriteString(ans)
		buffp.Write(Transform2Char(seq))
		buffp.WriteString("\n")

	}
	if err := buffp.Flush(); err != nil {
		log.Fatalf("[StoreEdgesToFn] failed to flush file: %s, err: %v\n", DcDBGEdgesfn, err)
	}
}

func SimplifyByLongReadsPath(edgesArr []DBGEdge, nodesArr []DBGNode, cs <-chan MapingInfo, opt optionsDDBG, MaxPathLen, MaxBubbleSeqLen, mapRecordSize int) (mergePathArr []Path) {
	//readPathArr := pathArr.Arr
	sortedEIDLenArr := SortByUniqueEdgeLen(edgesArr, nodesArr, opt.Haplotype, opt.Kmer)
	fmt.Printf("[SimplifyByLongReadsPath]len(sortedEIDLenArr):%d\n", len(sortedEIDLenArr))
	for i, el := range sortedEIDLenArr {
		fmt.Printf("[SimplifyByLongReadsPath] EIDLenArr[%d]:%d\n", i, el)
	}
	/*var readPathArr [][][][]uint32
	for {
		pm, ok := <-cs
		if !ok {
			break
		}
		readPathArr = append(readPathArr, pm)
	}*/
	//edgesPathRelationArr := make([][]uint32, len(sortedEIDLenArr))
	readPathArr, edgesPathRelationArr := ConstructEdgesPathRelationship(edgesArr, nodesArr, cs, mapRecordSize)

	//CleanNGSPath(edgesArr)
	t0 := time.Now()
	MergePathArr(edgesArr, nodesArr, readPathArr, sortedEIDLenArr, edgesPathRelationArr, opt, MaxPathLen, MaxBubbleSeqLen)
	fmt.Printf("[SimplifyByLongReadsPath] MergePathArr used: %v\n", time.Now().Sub(t0))

	//MergePathArr(edgesArr, nodesArr, pathArr, edgesPathRelationArr, opt.MinMapFreq)
	//MergePathMat(edgesArr, nodesArr, 2)
	//fmt.Printf("[SimplifyByLongReadsPath] sortedEIDIdxArr: %v\n", sortedEIDLenArr)
	ResetProcessFlagAll(edgesArr)
	mergePathArr = FindMaxPath(sortedEIDLenArr, edgesArr, nodesArr, readPathArr, edgesPathRelationArr, opt.MinMapFreq, 500, opt.Kmer, opt.Haplotype, MaxPathLen, MaxBubbleSeqLen)
	// constuct map edge ID to the path
	//IDMapPath := ConstructIDMapPath(pathArr)
	//DeleteJoinPathArrEnd(edgesArr, pathArr)
	return
	//return mergePathArr
}

func CompressToUint64(ks []byte) (cb uint64) {
	for _, b := range ks {
		cb <<= NumBitsInBase
		cb |= uint64(b)
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

func FindMinKmerID(ka []KmerID) (min KmerID, idx int) {
	min = ka[0]
	idx = 0
	for i, ki := range ka[1:] {
		if min.GetKmer() > ki.GetKmer() {
			min = ki
			idx = i + 1
		}
	}
	return
}

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

func GetSeqMiniKmerInfoArr(ks []byte, buf, kmerInfoArr []KI, ID uint32, SeedLen, WinSize int) []KI {
	if len(ks) < SeedLen+WinSize {
		var t []KI
		return t
	}
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
	fmt.Printf("[GetSeqMiniKmerInfoArr]bs:%X\n", bs)
	fmt.Printf("[GetSeqMiniKmerInfoArr]rs:%X\n", rs)
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
		fmt.Printf("[GetSeqMiniKmerInfoArr]bs:%X\n", bs)
		fmt.Printf("[GetSeqMiniKmerInfoArr]rs:%X\n", rs)
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

func RadixSortKIArr(kmerArr []KI, bulkArr [][]KI, SeedLen, BaseNum, NumBulk uint64) ([]KI, [][]KI) {
	if len(kmerArr) < 2 {
		return kmerArr, bulkArr
	}
	if len(kmerArr) < 10000 {
		sort.Sort(KIArr(kmerArr))
		return kmerArr, bulkArr
	}
	MUSK := NumBulk - 1
	overflowNum := 0
	countArr := make([]int, NumBulk)
	for i := uint64(0); i < (SeedLen+BaseNum-1)/BaseNum; i++ {
		// clean countArr
		for j := range countArr {
			countArr[j] = 0
		}
		for _, ki := range kmerArr {
			idx := (ki.GetKmer() >> (i * BaseNum * NumBitsInBase)) & MUSK
			if countArr[idx] >= len(bulkArr[idx]) {
				overflowNum++
				extSize := countArr[idx] * 3 / 2
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
		ct := 0
		for j := uint64(0); j < NumBulk; j++ {
			for x := 0; x < countArr[j]; x++ {
				kmerArr[ct] = bulkArr[j][x]
				ct++
			}
		}
		if ct != len(kmerArr) {
			log.Fatalf("[RadixSortKIArr] ct:%d != len(kmerArr):%d\n", ct, len(kmerArr))
		}
	}
	fmt.Printf("[RadixSortKIArr] overflowNum:%d\n", overflowNum)
	return kmerArr, bulkArr
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
	for i := uint64(0); i < sz; i++ {
		// clean countArr
		for j := range countArr {
			countArr[j] = 0
		}
		shiftNum := i * BitNum
		for _, ki := range kmerArr {
			idx := (ki.GetKmer() >> shiftNum) & MUSK
			if countArr[idx] >= len(bulkArr[idx]) {
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
		ct := 0
		for j := uint64(0); j < NumBulk; j++ {
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

type Uint32Slice []uint32

func (p Uint32Slice) Len() int           { return len(p) }
func (p Uint32Slice) Less(i, j int) bool { return p[i] < p[j] }
func (p Uint32Slice) Swap(i, j int)      { p[i], p[j] = p[j], p[i] }

type Uint16Slice []uint16

func (p Uint16Slice) Len() int           { return len(p) }
func (p Uint16Slice) Less(i, j int) bool { return p[i] < p[j] }
func (p Uint16Slice) Swap(i, j int)      { p[i], p[j] = p[j], p[i] }

func IsRepeatEdge(e *DBGEdge) bool {
	if e.GetUniqueFlag() == 0 && e.GetBubbleFlag() == 0 && e.GetBubbleRepeatFlag() == 0 && e.GetTwoEdgesCycleFlag() == 0 {
		return true
	}
	return false
}

func IsRepeatEdgeNoDB(e *DBGEdge, kmerlen int) bool {
	if e.GetUniqueFlag() == 0 && e.GetBubbleFlag() == 0 && e.GetBubbleRepeatFlag() == 0 && e.GetSeqLen() < kmerlen+EdgeMapKmerMin {
		return true
	}
	return false
}

func IsEdgeNoDB(e *DBGEdge, kmerlen int) bool {
	if IsRepeatEdge(e) {
		if e.GetSeqLen() < (kmerlen-1)*3/2 {
			return true
		}
	} else {
		if e.GetSeqLen() < kmerlen*5/4 {
			return true
		}
	}

	return false
}

type FreqInfo struct {
	Pos  uint32
	Freq uint16
}

func IndexFreqInfo(freqArr []FreqInfo, WinSize int, pos uint32) int {
	idx := -1
	step := WinSize * 4 / 7
	x := int(pos) / step
	if x >= len(freqArr) {
		x = len(freqArr) - 1
	}
	if freqArr[x].Pos == pos {
		idx = x
		return idx
	} else if freqArr[x].Pos < pos {
		x += int(pos-freqArr[x].Pos) / step
		if x >= len(freqArr) {
			x = len(freqArr) - 1
		}
	} else {
		x -= int(freqArr[x].Pos-pos) / step
		if x < 0 {
			x = 0
		}
	}

	if freqArr[x].Pos == pos {
		idx = x
	} else if freqArr[x].Pos < pos {
		for i := x + 1; i < len(freqArr); i++ {
			if freqArr[i].Pos == pos {
				idx = i
				break
			}
		}
	} else {
		for i := x - 1; i >= 0; i-- {
			if freqArr[i].Pos == pos {
				idx = i
				break
			}
		}
	}
	if idx < 0 {
		fmt.Printf("[IndexEdgeFreqArr]step:%d pos:%d freqArr:%v", step, pos, freqArr)
	}
	return idx
}

func CountEdgeMiniNum(freqArr []FreqInfo) (count int) {
	for _, fi := range freqArr {
		if fi.Freq == 0 || fi.Freq == math.MaxUint16 {
			count++
		}
	}

	return
}

func SortDBGEdgesKmer(dbgCS chan<- DBGKmerInfoPac, edgesArr []DBGEdge, SeedLen, WinSize, HighFreqKmerMin, HighFreqKmerMax, MinimerKmerNum, kmerlen int) {
	var dbgKmerInfoPac DBGKmerInfoPac
	var deleteKmerNum int
	var edgeNum, mapEdgeNum int
	MinEdgeLen = (kmerlen - 1) * 3 / 2
	kmerArr := make([]KI, 0, MinimerKmerNum)
	bufSize := WinSize * 3
	buf := make([]KI, bufSize)
	ka := make([]KI, 5000)
	//rs := make([]byte, 10000)
	for i := 0; i < len(edgesArr); i++ {
		e := &edgesArr[i]
		if e.ID < 2 || e.GetDeleteFlag() > 0 {
			continue
		}
		edgeNum++
		if e.GetSeqLen() > MaxAllowSeqLen {
			log.Fatalf("[SortDBGEdgesKmer]eID:%d seqLen:%d > MaxAllowSeqLen:%d\n", e.ID, e.GetSeqLen(), MaxAllowSeqLen)
		}
		ka = GetSeqMiniKmerInfoArr(e.Ks, buf, ka, uint32(e.ID), SeedLen, WinSize)
		if IsEdgeNoDB(e, kmerlen) {
			sort.Sort(KIArr(ka))
			e.SeedInfoArr = make([]SeedInfo, len(ka))
			for x, ki := range ka {
				e.SeedInfoArr[x] = ki.SeedInfo
			}
			edgesArr[e.ID].SeedInfoArr = e.SeedInfoArr
			continue
		}

		mapEdgeNum++
		//var ka []KI

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

	BitNum := uint64(3)
	NumBulk := uint64(1) << BitNum
	bulkArr := make([][]KI, NumBulk)
	for i := range bulkArr {
		bulkArr[i] = make([]KI, len(kmerArr)/int(NumBulk)*3/2)
	}
	kmerArr, _ = RadixSortKIArr3Bit(kmerArr, bulkArr, uint64(SeedLen), BitNum, NumBulk)
	totalKALen := len(kmerArr)
	//sort.Sort(KIArr(kmerArr))
	count := 0
	deleteCount := 0
	UniqueNum := 0
	highFreqKmerMap := make(map[uint64][]IDInfo, 2000)

	//sortKmerArr := make([]KI, count)
	//lowPercentMinimersEdgesMap := make(map[uint64][]uint32, 2000) // key is kmer, value eID array
	//lowFreqNumArr := make([]uint8, len(edgesArr))
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
		if j-i == 1 {
			UniqueNum++
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
	sortKmerArr := kmerArr[:skaIdx]

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
	for i := 0; i < len(edgesArr); i++ {
		e := &edgesArr[i]
		if e.ID < 2 || e.GetDeleteFlag() > 0 {
			continue
		}
		if IsEdgeNoDB(e, kmerlen) {
			continue
		}
		edgesArr[i].CovD = miniNumArr[i]
		if int(miniNumArr[i]) >= MinMinimerNum {
			edgesArr[i].SeedInfoArr = nil
		}

		//if int(miniNumArr[i]) < MinMinimerNum {
		//fmt.Printf("[SortDBGEdgesKmer]eID:%d eLen:%d miniNum:%d repeat:%t\n", e.ID, e.GetSeqLen(), miniNumArr[i], IsRepeatEdge(e))
		//}

	}

	totalHighFreq := 0
	for _, v := range highFreqKmerMap {
		totalHighFreq += len(v)
	}

	//fmt.Printf("[SortDBGEdgesKmer]highFreqInDBNum:%d percent:%.3f\n", highFreqInDBNum, float64(highFreqInDBNum)/float64(totalKALen))
	fmt.Printf("[SortDBGEdgesKmer]edgeNum:%d mapEdgeNum:%d\n", edgeNum, mapEdgeNum)
	fmt.Printf("[SortDBGEdgesKmer]deleteCount:%d deleteKmerNum:%d percent:%.4f\n", deleteCount, deleteKmerNum, float64(deleteKmerNum)/float64(totalKALen))
	fmt.Printf("[SortDBGEdgesKmer]HighFreqKmerMin:%d MinimerKmerNum:%d kmerArr len:%d count:%d \n", HighFreqKmerMin, MinimerKmerNum, len(kmerArr), count)
	fmt.Printf("[SortDBGEdgesKmer]UniqueNum:%d Unique percent:%v high Freq kmer num:%d avg high freq:%d high freq percent:%v\n", UniqueNum, float32(UniqueNum)/float32(totalKALen), len(highFreqKmerMap), totalHighFreq/(len(highFreqKmerMap)+1), float32(totalHighFreq)/float32(totalKALen))
	dbgKmerInfoPac.KmerinfoArr = sortKmerArr
	dbgKmerInfoPac.HighFreqMap = highFreqKmerMap
	//dbgKmerInfoPac.LowFreqPercentEdgeMap = lowPercentMinimersEdgesMap
	dbgCS <- dbgKmerInfoPac
	close(dbgCS)
	//time.Sleep(time.Second)
}

func GetReadArrSeqLen(ra []ReadInfo) (sl int) {
	for _, ri := range ra {
		sl += len(ri.Seq)
	}
	return
}

func SortLongReadsKmer(rc <-chan []ReadInfo, SeedLen, WinSize int, buckReadsNum int, cs chan<- ReadsBucketInfo, finishT chan<- int) {
	//var addReadKmerLen int64
	//rs := make([]byte, 10000)
	bufSize := WinSize * 10
	buf := make([]KI, bufSize)
	BitNum := uint64(3)
	NumBulk := uint64(1) << BitNum
	bulkArr := make([][]KI, NumBulk)
	for i := range bulkArr {
		bulkArr[i] = make([]KI, buckReadsNum*10000/int(NumBulk)*3/2)
	}
	ka := make([]KI, 5000)
	for {
		ra := <-rc
		if len(ra) == 0 {
			break
		}
		tl := GetReadArrSeqLen(ra)
		//var kmerArr []KI
		kmerArr := make([]KI, 0, tl*2/WinSize)
		seedInfoMat := make([][]SeedInfo, len(ra))
		for i, ri := range ra {
			if len(ri.Seq) > MaxAllowSeqLen {
				log.Fatalf("[SortLongReadsKmer]read ID:%d seqLen:%d > MaxAllowSeqLen:%d\n", ri.ID, len(ri.Seq), MaxAllowSeqLen)
			}
			ka = GetSeqMiniKmerInfoArr(ri.Seq, buf, ka, uint32(ri.ID), SeedLen, WinSize)
			kmerArr = append(kmerArr, ka...)
			seedInfoArr := make([]SeedInfo, len(ka))
			for x, ki := range ka {
				seedInfoArr[x] = ki.SeedInfo
			}
			seedInfoMat[i] = seedInfoArr
		}
		kmerArr, bulkArr = RadixSortKIArr3Bit(kmerArr, bulkArr, uint64(SeedLen), BitNum, NumBulk)
		//sort.Sort(KIArr(kmerArr))
		var rbi ReadsBucketInfo
		rbi.ReadsArr = ra
		rbi.SeedInfoMat = seedInfoMat
		rbi.KmerSortArr = kmerArr
		cs <- rbi
		fmt.Printf("[SortLongReadsKmer]finished sort reads num : %d, miniKmer num: %d\n", len(ra), len(kmerArr))
	}
	finishT <- 1
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
		log.Fatalf("[DeconstructDBG] check global Arguments error, opt: %v\n", gOpt)
	}
	opt.ArgsOpt = gOpt
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
	if opt.MinUniqueEdgeLen < 200 || opt.MinUniqueEdgeLen > 10000 {
		log.Fatalf("[checkArgs] argument 'MinUniqueEdgeLen': %v must between [200~10k]\n", c.Flag("MinUniqueEdgeLen"))
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
	var minPerc float64
	minPerc, ok = c.Flag("MinChainScoreIdentityPercent").Get().(float64)
	if !ok {
		log.Fatalf("[checkArgs] argument 'MinChainScoreIdentityPercent': %v set error\n ", c.Flag("MinChainScoreIdentityPercent").String())
	}
	if minPerc <= 0 || minPerc >= 1 {
		log.Fatalf("[checkArgs] argument 'MinChainScoreIdentityPercent': %v must between [1~100]\n", c.Flag("MinChainScoreIdentityPercent").String())
	}
	opt.MinChainScoreIdentityPercent = float32(minPerc)

	opt.Correct, ok = c.Flag("Correct").Get().(bool)
	if !ok {
		log.Fatalf("[checkArgs] argument 'Correct': %v set error\n ", c.Flag("Correct").String())
	}
	opt.Haplotype, ok = c.Flag("Haplotype").Get().(bool)
	if !ok {
		log.Fatalf("[checkArgs] argument 'Haplotype': %v set error\n ", c.Flag("Haplotype").String())
	}
	opt.Debug, ok = c.Flag("Debug").Get().(bool)
	if !ok {
		log.Fatalf("[checkArgs] argument 'Debug': %v set error\n ", c.Flag("Debug").String())
	}

	opt.ONTFn = c.Flag("LongReadFile").String()

	succ = true
	return opt, succ
}

//const seqPositionBitNum = 64 - 1 - SeedLen*NumBitsInBase
//const MaxAllowSeqLen = (1 << seqPositionBitNum) - 1
//const KmerMask = (1 << (seqPositionBitNum + 1)) - 1
//const PositionMask = (1 << seqPositionBitNum) - 1
//const KmerReSetPosition = (((1 << SeedLen * NumBitsInBase) - 1) << (seqPositionBitNum + 1)) | 0x1
//const KmerResetStrand = (1 << (SeedLen*NumBitsInBase + seqPositionBitNum)) - 1
const mapLenBitNum = 32 - 1 - seqPositionBitNum
const mapingKmerInfoMaxMapLen = (1 << mapLenBitNum) - 1
const infoReSetMapLen = (1 << (seqPositionBitNum + 1)) - 1
const ePositionMask = (1 << seqPositionBitNum) - 1
const infoReSetPosition = ((1<<(mapLenBitNum) - 1) << (seqPositionBitNum + 1)) | 0x1

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
		m.Info = m.Info & ePositionMask
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
		m.Info = m.Info & StrandMask
	}
}

/*const UsedFlag = 1 << 63
const UsedFlagMusk = 1<<63 - 1

func (m MapingKmerInfo) GetUsedFlag() bool {
	return m.SeedInfo&UsedFlag > 0
}

func (m *MapingKmerInfo) SetUsedFlag() {
	m.SeedInfo |= UsedFlag
}

func (m *MapingKmerInfo) ResetUsedFlag() {
	m.SeedInfo &= UsedFlagMusk
} */

/*func (m MapingKmerInfo) GetKmer() uint64 {
	var km uint64
	km = (uint64(m.RInfo) << 32) | uint64(m.EInfo)
	return km
}
func (m *MapingKmerInfo) SetKmer(km uint64) {
	m.RInfo = uint32(km >> 32)
	m.EInfo = uint32(km & (1<<32 - 1))
}  */

/*func (m MapingKmerInfo) GetRPos() uint32 {
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
	return (m.RInfo & 0x1) > 0
}
func (m MapingKmerInfo) GetEStrand() bool {
	return (m.EInfo & 0x1) > 0
}

func (m MapingKmerInfo) GetStrand() bool {
	return (m.RInfo & 0x1) == (m.EInfo & 0x1)
}*/

type MapingKmerInfoEPosArr []MapingKmerInfo

func (arr MapingKmerInfoEPosArr) Len() int {
	return len(arr)
}

func (arr MapingKmerInfoEPosArr) Less(i, j int) bool {
	return arr[i].GetEPosition() < arr[j].GetEPosition()
}

func (arr MapingKmerInfoEPosArr) Swap(i, j int) {
	arr[i], arr[j] = arr[j], arr[i]
}

type Uint32Arr []uint32 // uint32 array used sort

func (arr Uint32Arr) Len() int {
	return len(arr)
}

func (arr Uint32Arr) Less(i, j int) bool {
	return arr[i] < arr[j]
}

func (arr Uint32Arr) Swap(i, j int) {
	arr[i], arr[j] = arr[j], arr[i]
}

func GetLongReadsFile(cfgInfo CfgInfo) (fnArr []string) {
	for _, lib := range cfgInfo.Libs {
		if lib.SeqProfile == 3 || lib.SeqProfile == 4 {
			fnArr = append(fnArr, lib.FnName...)
		}
	}
	return
}

type ReadsBucketInfo struct {
	ReadsArr    []ReadInfo
	SeedInfoMat [][]SeedInfo
	KmerSortArr []KI
}

func LoadLongReads(ONTFnbr string, rc chan<- []ReadInfo, buckReadsNum int, processT chan<- int) {
	var processNumReads, totalBases int
	format := GetReadsFileFormat(ONTFnbr)
	bufSize := (1 << 20)
	size := 1000
	fncs1 := make(chan []byte, 3)
	recordChan1 := make(chan ReadInfo, size)

	fmt.Printf("[LoadLongReads] begin processed file: %v\n", ONTFnbr)
	go ReadBrFile(ONTFnbr, fncs1, bufSize)
	go GetReadFileRecord(fncs1, recordChan1, format, bufSize, true)
	readInfoArr := make([]ReadInfo, 0, 1000)
	var processCount int
	for {
		ri, ok := <-recordChan1
		if !ok {
			break
		}
		if ri.ID > math.MaxUint32 && ri.ID <= 0 {
			log.Fatalf("[LoadLongReads] read ID must [1~%v], ri.ID: %v\n", math.MaxUint32, ri.ID)
		}
		readInfoArr = append(readInfoArr, ri)
		processNumReads++
		totalBases += len(ri.Seq)
		processCount++
		if processCount >= buckReadsNum {
			rc <- readInfoArr
			var na []ReadInfo
			readInfoArr = na
			processCount = 0
		}
	}
	if len(readInfoArr) > 0 {
		rc <- readInfoArr
	}

	fmt.Printf("[LoadLongReads] processed reads number is: %d, totalBases: %v, finished processed file: %v\n", processNumReads, totalBases, ONTFnbr)
	processT <- 1
}

func LoadLongReadsAndSort(cfgInfo CfgInfo, NumCPU int, cs chan<- ReadsBucketInfo, SeedLen, WinSize int, buckReadsNum int) {
	fnArr := GetLongReadsFile(cfgInfo)
	rc := make(chan []ReadInfo, NumCPU)
	processCPUNum := NumCPU * 4 / 5
	LoadCPUNum := NumCPU - processCPUNum
	if processCPUNum <= 2 {
		processCPUNum = 1
		LoadCPUNum = 1
	}
	loadT := make(chan int, LoadCPUNum)
	finishedT := make(chan int, processCPUNum)
	for i := 0; i < processCPUNum; i++ {
		go SortLongReadsKmer(rc, SeedLen, WinSize, buckReadsNum, cs, finishedT)
	}
	for i := 0; i < LoadCPUNum; i++ {
		loadT <- 1
	}

	for i := 0; i < len(fnArr); i++ {
		<-loadT
		go LoadLongReads(fnArr[i], rc, buckReadsNum, loadT)
	}

	for i := 0; i < LoadCPUNum; i++ {
		<-loadT
	}
	close(rc)
	//time.Sleep(time.Second * 60)
	for i := 0; i < processCPUNum; i++ {
		<-finishedT
	}
	close(cs)
	//time.Sleep(time.Second * 60)
}

func GetShareKmer(edgesKmerSortArr, kiArr []KI, firstReadID uint32, readsNum int) [][]MapingKmerInfo {
	mkiArr := make([][]MapingKmerInfo, readsNum)
	i, j := 0, 0
	for ; i < len(kiArr) && j < len(edgesKmerSortArr); i++ {
		kiR, kiE := kiArr[i], edgesKmerSortArr[j]
		if kiR.GetKmer() < kiE.GetKmer() {
			continue
		} else if kiR.GetKmer() > kiE.GetKmer() {
			for ; j < len(edgesKmerSortArr); j++ {
				kiE = edgesKmerSortArr[j]
				if kiR.GetKmer() <= kiE.GetKmer() {
					break
				}
			}
		}
		if kiR.GetKmer() < kiE.GetKmer() {
			continue
		} else if kiR.GetKmer() > kiE.GetKmer() {
			break
		}
		// kiR.GetKmer() == kiE.GetKmer()
		for m := j; m < len(edgesKmerSortArr); m++ {
			tmp := edgesKmerSortArr[m]
			if tmp.GetKmer() > kiR.GetKmer() {
				break
			}
			var mki MapingKmerInfo
			mki.SetKmer(kiR.GetKmer())
			mki.EID = kiE.ID
			mki.SetPosition(kiR.GetPosition())
			mki.SetStrand(kiR.GetStrand())
			mki.SetEPosition(uint32(kiE.GetPosition()))
			mki.SetEStrand(kiE.GetStrand())
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
				if mki.GetEStrand() == mki.GetStrand() {
					PLUSNum++
				} else {
					MINUSNum++
				}
			}
			var strand bool
			if PLUSNum > MINUSNum {
				strand = PLUS
			} else {
				strand = MINUS
			}
			edgeMapInfoArr[idx].Strand = strand
			// choose chain by strand
			for _, mki := range v {
				if strand {
					if mki.GetEStrand() != mki.GetStrand() {
						continue
					}
				} else {
					if mki.GetEStrand() == mki.GetStrand() {
						continue
					}
				}
				var c Chain
				c.X, c.Y = uint32(mki.GetEPosition()), uint32(mki.GetPosition())
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
	EdgeID   uint32
	EdgeLen  uint32
	MPIIdx   uint16
	Distance uint16
	Strand   bool
	MPI      *MaxPathInfo
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
		if mki.GetEPosition() < startP {
			continue
		} else if mki.GetEPosition()+SeedLen > endP {
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
		scArr[i].Sc = int(ms)
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

/*func FindONTReadsPath(edgesKmerSortArr []KI, cs <-chan ReadsBucketInfo, SeedLen, WinSize int, wc chan<- [2][]uint32) {
	for {
		rbi := <-cs
		if len(rbi.ReadsArr) == 0 {
			var tmp [2][]uint32
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
	Path         [2][]uint32
}

func GetInterSectionKmerInfo(dbgKmerInfoPac DBGKmerInfoPac, kmerInfoArrB []KI, BLen int, StartID uint32, SeedLen uint16) (bucketInterSecKmerInfoArr [][]MapingKmerInfo) {
	//startReadID := uint32(rb.ReadsArr[0].ID)
	kmerInfoArrA := dbgKmerInfoPac.KmerinfoArr
	bucketInterSecKmerInfoArr = make([][]MapingKmerInfo, BLen)
	for i := 0; i < BLen; i++ {
		bucketInterSecKmerInfoArr[i] = make([]MapingKmerInfo, 0, 400000)
	}
	i, j := 0, 0
	ek := kmerInfoArrA[j]
	for ; i < len(kmerInfoArrB) && j < len(kmerInfoArrA); i++ {
		//fmt.Printf("[GetInterSectionKmerInfo]kmerInfoArrB[%v]: %v, kmerInfoArrA[%v]: %v\n", i, kmerInfoArrB[i], j, kmerInfoArrA[j])
		rk := kmerInfoArrB[i]
		if rk.GetKmer() < ek.GetKmer() {
			continue
		} else if rk.GetKmer() > ek.GetKmer() {
			for j = j + 1; j < len(kmerInfoArrA); j++ {
				ek = kmerInfoArrA[j]
				if rk.GetKmer() <= ek.GetKmer() {
					break
				}
			}
		}

		if rk.GetKmer() < ek.GetKmer() {
			continue
		} else if rk.GetKmer() > ek.GetKmer() {
			break
		}
		// rk.GetKmer() == ek.GetKmer()
		for m := j; m < len(kmerInfoArrA); m++ {
			tk := kmerInfoArrA[m]
			if tk.GetKmer() > rk.GetKmer() {
				break
			}
			var mk MapingKmerInfo
			//mk.SetKmer(tk.GetKmer())
			mk.SeedInfo = rk.SeedInfo
			mk.EID = tk.ID
			mk.SetEPosition(uint32(tk.GetPosition()))
			mk.SetEStrand(tk.GetStrand())
			mk.SetMapLen(uint32(SeedLen))
			//fmt.Printf("[GetInterSectionKmerInfo]rk: %v, tk: %v\n", rk, tk)
			//fmt.Printf("[GetInterSectionKmerInfo]edgeID: %v: ReadPos: %v, EdgePos: %v.SeedInfo: %v\n", rk.ID, tk.GetPos(), rk.GetPos(), tk.GetKmer())
			bucketInterSecKmerInfoArr[rk.ID-StartID] = append(bucketInterSecKmerInfoArr[rk.ID-StartID], mk)
		}
	}

	return
}

func GetInterSectionKmerInfoB(dbgKmerInfoPac DBGKmerInfoPac, kmerInfoArrB []KI, StartID uint32, SeedLen uint16, bucketInterSecKmerInfoArr [][]MapingKmerInfo) [][]MapingKmerInfo {
	//startReadID := uint32(rb.ReadsArr[0].ID)
	// clean bucketInterSecKmerInfoArr
	for i := range bucketInterSecKmerInfoArr {
		bucketInterSecKmerInfoArr[i] = bucketInterSecKmerInfoArr[i][:0]
	}
	kmerInfoArrA := dbgKmerInfoPac.KmerinfoArr
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
		mk.SetEPosition(uint32(ek.GetPosition()))
		mk.SetEStrand(ek.GetStrand())
		mk.SetMapLen(uint32(SeedLen))
		for m := j; m < len(kmerInfoArrB); m++ {
			tk := kmerInfoArrB[m]
			if tk.GetKmer() > ek.GetKmer() {
				break
			}
			mk.SeedInfo = tk.SeedInfo
			//fmt.Printf("[GetInterSectionKmerInfo]rk: %v, tk: %v\n", rk, tk)
			//fmt.Printf("[GetInterSectionKmerInfo]edgeID: %v: ReadPos: %v, EdgePos: %v.SeedInfo: %v\n", rk.ID, tk.GetPos(), rk.GetPos(), tk.GetKmer())
			bucketInterSecKmerInfoArr[tk.ID-StartID] = append(bucketInterSecKmerInfoArr[tk.ID-StartID], mk)
		}
	}

	return bucketInterSecKmerInfoArr
}

func GetInterSectionKmerInfo2(kmerInfoArrA, kmerInfoArrB []KI, BLen int, StartID uint32, SeedLen uint16) (bucketInterSecKmerInfoArr [][]MapingKmerInfo) {
	//startReadID := uint32(rb.ReadsArr[0].ID)
	bucketInterSecKmerInfoArr = make([][]MapingKmerInfo, BLen)
	for i := 0; i < BLen; i++ {
		bucketInterSecKmerInfoArr[i] = make([]MapingKmerInfo, 0, 20)
	}
	i, j := 0, 0
	ek := kmerInfoArrA[j]
	for ; i < len(kmerInfoArrB) && j < len(kmerInfoArrA); i++ {
		//fmt.Printf("[GetInterSectionKmerInfo]kmerInfoArrB[%v]: %v, kmerInfoArrA[%v]: %v\n", i, kmerInfoArrB[i], j, kmerInfoArrA[j])
		rk := kmerInfoArrB[i]
		if rk.GetKmer() < ek.GetKmer() {
			continue
		} else if rk.GetKmer() > ek.GetKmer() {
			for j = j + 1; j < len(kmerInfoArrA); j++ {
				ek = kmerInfoArrA[j]
				if rk.GetKmer() <= ek.GetKmer() {
					break
				}
			}
		}

		if rk.GetKmer() < ek.GetKmer() {
			continue
		} else if rk.GetKmer() > ek.GetKmer() {
			break
		}
		// rk.GetKmer() == ek.GetKmer()
		for m := j; m < len(kmerInfoArrA); m++ {
			tk := kmerInfoArrA[m]
			if tk.GetKmer() > rk.GetKmer() {
				break
			}
			var mk MapingKmerInfo
			mk.SeedInfo = tk.SeedInfo
			mk.EID = rk.ID
			mk.SetEPosition(uint32(rk.GetPosition()))
			mk.SetEStrand(rk.GetStrand())
			mk.SetMapLen(uint32(SeedLen))
			//fmt.Printf("[GetInterSectionKmerInfo]rk: %v, tk: %v\n", rk, tk)

			//fmt.Printf("[GetInterSectionKmerInfo]edgeID: %v: ReadPos: %v, EdgePos: %v.SeedInfo: %v\n", rk.ID, tk.GetPos(), rk.GetPos(), tk.GetKmer())
			bucketInterSecKmerInfoArr[rk.ID-StartID] = append(bucketInterSecKmerInfoArr[rk.ID-StartID], mk)
		}
	}

	return
}

func GetChainBlocks(ka, buf []MapingKmerInfo, esArr []EdgeStrand) ([]MapingKmerInfo, []MapingKmerInfo, []EdgeMKINum) {
	if cap(buf) < len(ka) {
		buf = make([]MapingKmerInfo, len(ka))
	} else {
		buf = buf[:len(ka)]
	}
	enArr := Transform2EdgeMKINumArr(esArr)
	for _, mki := range ka {
		enArr[IndexEdgeMKINumArr(enArr, mki.EID)].Num++
	}
	sum := uint32(0)
	for i, en := range enArr {
		enArr[i].Num = sum
		sum += en.Num
	}

	if int(sum) != len(ka) {
		log.Fatalf("[GetChainBlocks]sum:%d != len(ka):%d\n", sum, len(ka))
	}

	for _, mki := range ka {
		buf[enArr[IndexEdgeMKINumArr(enArr, mki.EID)].Num] = mki
		enArr[IndexEdgeMKINumArr(enArr, mki.EID)].Num++
	}
	// sort by mapping edge position
	{
		lastIdx := 0
		eID := buf[lastIdx].EID
		i := 1
		for ; i < len(buf); i++ {
			if buf[i].EID != eID {
				sort.Sort(MapingKmerInfoEPosArr(buf[lastIdx:i]))
				lastIdx = i
				eID = buf[lastIdx].EID
			}
		}
		if i-lastIdx > 0 {
			sort.Sort(MapingKmerInfoEPosArr(buf[lastIdx:i]))
		}
	}
	//PrintMKIArr(buf)
	// clean edgeCountMap

	//kb = make([]MapingKmerInfo, 0, len(ka)/5+2)
	ka = ka[:0]
	for i := 0; i < len(buf); i++ {
		mk := buf[i]
		if mk.GetMapLen() == 0 {
			continue
		}
		//fmt.Printf("[GetChainBlocks] %v\n", PrintMKI(mk))
		mE := mk.GetEPosition()
		last := mk
		lR := uint32(last.GetPosition())
		lE := mE
		strand := last.GetMapStrand()
		j := i + 1
		for ; j < len(buf); j++ {
			nk := buf[j]
			nR := uint32(nk.GetPosition())
			nE := nk.GetEPosition()
			if nk.EID != mk.EID || nE > lE+last.GetMapLen() {
				break
			}

			if nk.GetMapStrand() != strand {
				continue
			}
			//fmt.Printf("[GetChainBlocks] %v\n", PrintMKI(nk))
			if strand == PLUS {
				if nR-lR == nE-lE {
					mk.SetMapLen(nE + nk.GetMapLen() - mE)
					buf[j].SetMapLen(0)
					last = nk
					lR, lE = nR, nE
				}
			} else {
				if nE-lE == (lR+last.GetMapLen())-(nR+nk.GetMapLen()) {
					//fmt.Printf("[GetChainBlocks] EID: %v, mk RPos: %v, EPos: %v,Len: %v, nk RPos: %v, EPos: %v, Len:%v\n", mk.EID, mk.GetPosition(), mk.GetEPosition(), mk.Len, nk.GetPosition(), nk.GetEPosition(), nk.GetMapLen())
					mk.SetPosition(uint64(nR))
					mk.SetMapLen(nE + nk.GetMapLen() - mE)
					//mk.SetEPosition(nE)
					//mk.GetMapLen() += uint16(((nk.GetPosition() + uint32(nk.GetMapLen())) - (mk.GetPosition() + uint32(mk.GetMapLen()))))
					buf[j].SetMapLen(0)
					last = nk
					lR, lE = nR, nE
					//fmt.Printf("[GetChainBlocks] mk RPos: %v, EPos: %v,Len: %v\n", mk.GetPosition(), mk.GetEPosition(), mk.GetMapLen())
				}
			}
		}
		//buf[i].GetMapLen() = 0
		//fmt.Printf("[GetChainBlocks]added mk: %v\n", PrintMKI(mk))
		ka = append(ka, mk)
		if len(ka) == 1 {
			enArr[IndexEdgeMKINumArr(enArr, mk.EID)].Num = 0
		} else if mk.EID != ka[len(ka)-2].EID {
			enArr[IndexEdgeMKINumArr(enArr, mk.EID)].Num = uint32(len(ka) - 1)
		}
	}
	return ka, buf, enArr
}

func GetChainBlocksOneEdge(ka, kb []MapingKmerInfo) []MapingKmerInfo {
	kb = kb[:0]
	for i := 0; i < len(ka); i++ {
		mk := ka[i]
		if mk.GetMapLen() == 0 {
			continue
		}
		mR := uint32(mk.GetPosition())
		mE := mk.GetEPosition()
		j := i + 1
		last := mk
		lR := mR
		lE := mE
		for ; j < len(ka); j++ {
			nk := ka[j]
			nR := uint32(nk.GetPosition())
			nE := nk.GetEPosition()
			if nR > lR+last.GetMapLen() {
				break
			}
			if nk.GetStrand() != last.GetStrand() {
				continue
			}
			strand := last.GetStrand() == last.GetEStrand()
			if strand == PLUS {
				if nR-lR == nE-lE {
					mk.SetMapLen(nR + nk.GetMapLen() - mR)
					ka[j].SetMapLen(0)
					last = nk
					lR, lE = nR, nE
				}
			} else {
				if nR-lR == (lE+last.GetMapLen())-(nE+nk.GetMapLen()) {
					//fmt.Printf("[GetChainBlocks] EID: %v, mk RPos: %v, EPos: %v,Len: %v, nk RPos: %v, EPos: %v, Len:%v\n", mk.EID, mk.GetPosition(), mk.GetEPosition(), mk.GetMapLen(), nk.GetPosition(), nk.GetEPosition(), nk.GetMapLen())
					mk.SetMapLen(nR + nk.GetMapLen() - mR)
					mk.SetEPosition(nE)
					//mk.GetMapLen() += uint16(((nk.GetPosition() + uint32(nk.GetMapLen())) - (mk.GetPosition() + uint32(mk.GetMapLen()))))
					ka[j].SetMapLen(0)
					last = nk
					lR, lE = nR, nE
					//fmt.Printf("[GetChainBlocks] mk RPos: %v, EPos: %v,Len: %v\n", mk.GetPosition(), mk.GetEPosition(), mk.GetMapLen())
				}
			}
		}
		ka[i].SetMapLen(0)
		kb = append(kb, mk)
	}
	return kb
}

/*func GetChainBlocks(ka []MapingKmerInfo) (kb []MapingKmerInfo) {
	flagArr := make([]bool, len(ka))
	kb = make([]MapingKmerInfo, 0, len(ka))
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
				if nk.GetPosition() > mk.GetPosition()+uint32(mk.GetMapLen()) {
					break
				}
				if mk.EID != nk.EID || mk.GetStrand() != nk.GetStrand() {
					continue
				}
				//fmt.Printf("[GetChainBlocks] EID: %v, mk RPos: %v, EPos: %v,Len: %v, nk RPos: %v, EPos: %v, Len:%v\n", mk.EID, mk.GetPosition(), mk.GetEPosition(), mk.GetMapLen(), nk.GetPosition(), nk.GetEPosition(), nk.GetMapLen())
				strand := mk.GetStrand()
				if strand == PLUS {
					if nk.GetPosition()-mk.GetPosition() == nk.GetEPosition()-mk.GetEPosition() {
						mk.GetMapLen() = uint16(nk.GetPosition() + uint32(nk.GetMapLen()) - mk.GetPosition())
						//mk.GetMapLen() += uint16(((nk.GetPosition() + uint32(nk.GetMapLen())) - (mk.GetPosition() + uint32(mk.GetMapLen()))))
						flagArr[j] = true
					}
				} else {
					if nk.GetPosition()-mk.GetPosition() == (mk.GetEPosition()+uint32(mk.GetMapLen()))-(nk.GetEPosition()+uint32(nk.GetMapLen())) {
						//fmt.Printf("[GetChainBlocks] EID: %v, mk RPos: %v, EPos: %v,Len: %v, nk RPos: %v, EPos: %v, Len:%v\n", mk.EID, mk.GetPosition(), mk.GetEPosition(), mk.GetMapLen(), nk.GetPosition(), nk.GetEPosition(), nk.GetMapLen())
						mk.GetMapLen() = uint16(nk.GetPosition() + uint32(nk.GetMapLen()) - mk.GetPosition())
						mk.SetEPos(nk.GetEPosition())
						//mk.GetMapLen() += uint16(((nk.GetPosition() + uint32(nk.GetMapLen())) - (mk.GetPosition() + uint32(mk.GetMapLen()))))
						flagArr[j] = true
						//fmt.Printf("[GetChainBlocks] mk RPos: %v, EPos: %v,Len: %v\n", mk.GetPosition(), mk.GetEPosition(), mk.GetMapLen())
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
	if loopNum > 2 {
		fmt.Printf("[GetChainBlocks] loop num: %v\n", loopNum)
		sort.Sort(MapingKmerInfoArr(kb))
	}

	return
}*/

type EdgeKmerInfo struct {
	PlusNum, MinusNum int
	ScoreP, ScoreM    int
	IdxP, IdxM        int
	MkInfoArr         []MapingKmerInfo
}

type ChainBlockMapLongReadInfo struct {
	//EID        uint32
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

func GetEdgesKIArr(kb []MapingKmerInfo, edgesArr []DBGEdge, MinKmerScore, LongEdgeMinLen int) (edgesKIArr []EdgeKI) {
	MinScore := MinKmerScore
	edgeKIMap := make(map[uint32]EdgeKmerInfo)
	for i, mk := range kb {
		//var eki EdgeKmerInfo
		eki, ok := edgeKIMap[mk.EID]
		if !ok {
			eki.MkInfoArr = make([]MapingKmerInfo, 0, 20)
		}
		eki.MkInfoArr = append(eki.MkInfoArr, mk)
		if mk.GetStrand() == PLUS {
			eki.PlusNum++
			eki.ScoreP += int(mk.GetMapLen())
			if mk.GetMapLen() > kb[eki.IdxP].GetMapLen() {
				eki.IdxP = i
			}
		} else {
			eki.MinusNum++
			eki.ScoreM += int(mk.GetMapLen())
			if mk.GetMapLen() > kb[eki.IdxM].GetMapLen() {
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
		edgeLen := len(edgesArr[kb[v.IdxM].EID].Ks)
		if edgeLen < LongEdgeMinLen {
			if v.ScoreM < edgeLen/10 && v.ScoreP < edgeLen/10 {
				continue
			}
		} else {
			if v.ScoreM < LongEdgeMinLen/10 && v.ScoreP < LongEdgeMinLen/10 {
				continue
			}
		}

		var eKI EdgeKI
		if v.ScoreP > v.ScoreM {
			eKI.Strand = PLUS
			eKI.Score = v.ScoreP
			eKI.Idx = v.IdxP
		} else {
			eKI.Strand = MINUS
			eKI.Score = v.ScoreM
			eKI.Idx = v.IdxM
		}
		eKI.MkInfoArr = v.MkInfoArr
		edgesKIArr = append(edgesKIArr, eKI)
	}
	return
}

func GetMKIDistance(cen, mk MapingKmerInfo, edgesArr []DBGEdge, distanceE int) (dX, dY int) {
	if mk.GetStrand() == PLUS {
		if cen.GetPosition() <= mk.GetPosition() {
			dX = int(mk.GetPosition()) - (int(cen.GetPosition()) + int(cen.GetMapLen()))
			if cen.EID != mk.EID {
				dY = int(mk.GetEPosition()) + distanceE
			} else {
				dY = int(mk.GetEPosition()) - (int(cen.GetEPosition()) + int(cen.GetMapLen()))
			}
		} else {
			dX = int(cen.GetPosition()) - (int(mk.GetPosition()) + int(mk.GetMapLen()))
			if cen.EID != mk.EID {
				dY = len(edgesArr[mk.EID].Ks) - (int(mk.GetEPosition()) + int(mk.GetMapLen())) + distanceE
			} else {
				dY = int(cen.GetEPosition()) - (int(mk.GetEPosition()) + int(mk.GetMapLen()))
			}
		}
	} else {
		if cen.GetPosition() <= mk.GetPosition() {
			dX = int(mk.GetPosition()) - (int(cen.GetPosition()) + int(cen.GetMapLen()))
			if cen.EID != mk.EID {
				dY = len(edgesArr[mk.EID].Ks) - (int(mk.GetEPosition()) + int(mk.GetMapLen())) + distanceE
			} else {
				dY = int(cen.GetEPosition()) - (int(mk.GetEPosition()) + int(mk.GetMapLen()))
			}
		} else {
			dX = int(cen.GetPosition()) - (int(mk.GetPosition()) + int(mk.GetMapLen()))
			if cen.EID != mk.EID {
				dY = int(mk.GetEPosition()) + distanceE
			} else {
				dY = int(mk.GetEPosition()) - (int(cen.GetEPosition()) + int(cen.GetMapLen()))
			}
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

func GetTwoBlockScore(cen, mk MapingKmerInfo, GapCost int, edgesArr []DBGEdge, distanceE int) (sc int) {
	var gapLen int
	dX, dY := GetMKIDistance(cen, mk, edgesArr, distanceE)

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
	if gapLen > min/6+20 {
		sc = 0
		return
	}
	a := int(mk.GetMapLen())
	// a gap length of gapLen costs
	var b float64
	if gapLen > 1 {
		//b = float64(0.01)*float64(mk.GetMapLen())*float64(gapLen) + float64(0.5)*math.Log2(float64(gapLen))
		b = float64(0.01)*float64(mk.GetMapLen())*float64(gapLen) + float64(0.5)*float64(gapLen/10)
	} else {
		if dX == 0 && dY == 0 {
			b = 0
		} else {
			b = 1
		}
	}
	sc = a - int(b)
	//sc = int(mk.GetMapLen()) - gapLen*GapCost
	//fmt.Printf("[GetTwoBlockScore]dX: %v, dY: %v, sc: %v, a: %v, b: %v\n", dX, dY, sc, a, b)
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

func CutMapingKmerBoundary(cen, mk MapingKmerInfo, direction uint8, strand bool) MapingKmerInfo {
	if cen.EID != mk.EID {
		log.Fatalf("[CutMapingKmerBoundary] cen.EID: %v != mk.EID: %v\n", cen.EID, mk.EID)
	}
	if direction == FORWARD {
		dR := int(cen.GetPosition()) + int(cen.GetMapLen()) - int(mk.GetPosition())
		if strand == PLUS {
			dE := int(cen.GetEPosition()) + int(cen.GetMapLen()) - int(mk.GetEPosition())
			max := MaxInt(dR, dE)
			if max > 0 {
				if max > int(mk.GetMapLen()) {
					max = int(mk.GetMapLen())
				}
				mk.SetPosition(mk.GetPosition() + uint64(max))
				mk.SetEPosition(mk.GetEPosition() + uint32(max))
				mk.SetMapLen(mk.GetMapLen() - uint32(max))
			}
		} else {
			dE := int(mk.GetEPosition()) + int(mk.GetMapLen()) - int(cen.GetEPosition())
			max := MaxInt(dR, dE)
			if max > 0 {
				if max > int(mk.GetMapLen()) {
					max = int(mk.GetMapLen())
				}
				mk.SetPosition(mk.GetPosition() + uint64(max))
				//mk.SetEPos(mk.GetEPosition() + max)
				mk.SetMapLen(mk.GetMapLen() - uint32(max))
			}
		}
	} else {
		dR := int(mk.GetPosition()) + int(mk.GetMapLen()) - int(cen.GetPosition())
		if strand == PLUS {
			dE := int(mk.GetEPosition()) + int(mk.GetMapLen()) - int(cen.GetEPosition())
			max := MaxInt(dR, dE)
			if max > 0 {
				if max > int(mk.GetMapLen()) {
					max = int(mk.GetMapLen())
				}
				mk.SetMapLen(mk.GetMapLen() - uint32(max))
			}
		} else {
			dE := int(cen.GetEPosition()) + int(cen.GetMapLen()) - int(mk.GetEPosition())
			max := MaxInt(dR, dE)
			if max > 0 {
				if max > int(mk.GetMapLen()) {
					max = int(mk.GetMapLen())
				}
				mk.SetEPosition(mk.GetEPosition() + uint32(max))
				mk.SetMapLen(mk.GetMapLen() - uint32(max))
			}
		}
	}
	return mk
}

func AddToIdxScArr(idxScArr IdxScArr, idxSc IdxScore) IdxScArr {
	if idxScArr.ArrCount < len(idxScArr.Arr) {
		idxScArr.Arr[idxScArr.ArrCount] = idxSc
		idxScArr.ArrCount++
		iS, idx := GetMinIdxSc(idxScArr)
		idxScArr.MinIdx = idx
		idxScArr.MinScore = iS.Sc
	} else {
		if idxSc.Sc > idxScArr.MinScore {
			idxScArr.Arr[idxScArr.MinIdx] = idxSc
			iS, idx := GetMinIdxSc(idxScArr)
			idxScArr.MinIdx = idx
			idxScArr.MinScore = iS.Sc
		}
	}

	return idxScArr
}

//const MaxUint14 uint16 = (math.MaxUint16 >> 2)

func BlockToNearBlocksScoreArr(kb []MapingKmerInfo, Idx int, strand bool, GapCost int, bsA []uint8, SeedLen uint16) {
	cen := kb[Idx]
	//maxExtendNum := 6
	//count := 0
	extendSize := 200
	//rightExtendSize := 100
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

func BlockToNearBlocksScoreArrCrossEdges(kb []MapingKmerInfo, Idx int, strand bool, visitArr []bool, GapCost, MaxGapLen int, SeqType uint8, bsA []int, direction uint8, flankLen, kmerlen int, edgesArr []DBGEdge, SeedLen uint16, distanceE int) {
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
		/*if int(cen.GetPosition())-int(mk.GetPosition()) > MaxGapLen {
			break
		}*/
		if strand == PLUS {
			if mk.EID != cen.EID {
				dY := distanceE + (len(edgesArr[mk.EID].Ks) - (int(mk.GetEPosition()) + int(mk.GetMapLen())))
				if dY < 0 || mk.GetPosition()+uint64(mk.GetMapLen()) > cen.GetPosition() {
					bsA[i] = math.MinInt32
					continue
				}
			} else {
				if mk.GetPosition()+uint64(mk.GetMapLen()) > cen.GetPosition() || mk.GetEPosition()+uint32(mk.GetMapLen()) > cen.GetEPosition() {
					mk = CutMapingKmerBoundary(cen, mk, BACKWARD, strand)
					if mk.GetMapLen() >= uint32(SeedLen) {
						kb[i] = mk
						fmt.Printf("[BlockToNearBlocksScoreArr] i: %v, mk ReadPos: %v, EPos: %v\n", i, mk.GetPosition(), mk.GetEPosition())
					} else {
						bsA[i] = math.MinInt32
						continue
					}
				}
			}

			sc := GetTwoBlockScore(cen, mk, GapCost, edgesArr, distanceE)
			if sc > bsA[i] {
				bsA[i] = sc
			}
		} else { // strand == MINUS
			if mk.EID != cen.EID {
				dY := distanceE + int(mk.GetEPosition())
				if dY < 0 || mk.GetPosition()+uint64(mk.GetMapLen()) > cen.GetPosition() {
					bsA[i] = math.MinInt32
					continue
				}
			} else {
				if mk.GetPosition()+uint64(mk.GetMapLen()) > cen.GetPosition() || mk.GetEPosition() < cen.GetEPosition()+uint32(cen.GetMapLen()) {
					mk = CutMapingKmerBoundary(cen, mk, BACKWARD, strand)
					if mk.GetMapLen() >= uint32(SeedLen) {
						kb[i] = mk
						fmt.Printf("[BlockToNearBlocksScoreArr] i: %v, mk ReadPos: %v, EPos: %v\n", i, mk.GetPosition(), mk.GetEPosition())
					} else {
						bsA[i] = math.MinInt32
						continue
					}
				}
			}
			sc := GetTwoBlockScore(cen, mk, GapCost, edgesArr, distanceE)
			if sc > bsA[i] {
				bsA[i] = sc
			}
		}
	}

	// FORWARD
	for i := Idx + 1; i < len(kb); i++ {
		mk := kb[i]
		//fmt.Printf("[BlockToNearBlocksScoreArr] i: %v, mk ReadPos: %v, cen ReadPos: %v\n", i, mk.GetPosition(), cen.GetPosition())
		if visitArr[i] == true {
			break
		}
		if bsA[i] == math.MinInt32 {
			continue
		}
		if mk.GetStrand() != strand {
			continue
		}
		/*if int(mk.GetPosition())-(int(cen.GetPosition())+int(mk.GetMapLen())) > MaxGapLen {
			break
		}*/
		var sc int
		if strand == PLUS {
			if cen.EID != mk.EID {
				dY := distanceE + int(mk.GetEPosition())
				if dY < 0 || cen.GetPosition()+uint64(cen.GetMapLen()) > mk.GetPosition() {
					bsA[i] = math.MinInt32
					continue
				}
			} else {
				if cen.GetPosition()+uint64(cen.GetMapLen()) > mk.GetPosition() || cen.GetEPosition()+uint32(cen.GetMapLen()) > mk.GetEPosition() {
					mk = CutMapingKmerBoundary(cen, mk, FORWARD, strand)
					if mk.GetMapLen() >= uint32(SeedLen) {
						kb[i] = mk
						fmt.Printf("[BlockToNearBlocksScoreArr] i: %v, mk ReadPos: %v, EPos: %v\n", i, mk.GetPosition(), mk.GetEPosition())
					} else {
						bsA[i] = math.MinInt32
						continue
					}
				}
			}
			sc = GetTwoBlockScore(cen, mk, GapCost, edgesArr, distanceE)

			if sc > bsA[i] {
				bsA[i] = sc
			}
		} else { //strand == MINUS
			if cen.EID != mk.EID {
				dY := distanceE + len(edgesArr[mk.EID].Ks) - (int(mk.GetEPosition()) + int(mk.GetMapLen()))
				if dY < 0 || cen.GetPosition()+uint64(cen.GetMapLen()) > mk.GetPosition() {
					bsA[i] = math.MinInt32
					continue
				}
			} else {
				if cen.GetPosition()+uint64(cen.GetMapLen()) > mk.GetPosition() || mk.GetEPosition()+uint32(mk.GetMapLen()) > cen.GetEPosition() {
					mk = CutMapingKmerBoundary(cen, mk, FORWARD, strand)
					if mk.GetMapLen() >= uint32(SeedLen) {
						kb[i] = mk
						fmt.Printf("[BlockToNearBlocksScoreArr] i: %v, mk ReadPos: %v, EPos: %v\n", i, mk.GetPosition(), mk.GetEPosition())
					} else {
						bsA[i] = math.MinInt32
						continue
					}
				}
			}
			sc = GetTwoBlockScore(cen, mk, GapCost, edgesArr, distanceE)
			if sc > bsA[i] {
				bsA[i] = sc
			}
		}
	}

	if distanceE != 0 {
		//fmt.Printf("[BlockToNearBlocksScoreArrCrossEdges] bsA: %v\n", bsA)
	}
	//sort.Sort(ScoreArr(bsA))
	return
}

func GetMinIdxSc(idxScArr IdxScArr) (idxSc IdxScore, idx int) {
	idx = -1
	for i, is := range idxScArr.Arr[:idxScArr.ArrCount] {
		if is.Sc > 0 && (idxSc.Sc == 0 || is.Sc < idxSc.Sc) {
			idxSc = is
			idx = i
		}
	}
	return
}

func GetMaxIdxSc(idxScArr IdxScArr) (idxSc IdxScore, idx int) {
	idx = -1
	for i, is := range idxScArr.Arr[:idxScArr.ArrCount] {
		if is.Sc > 0 && (idxSc.Sc == 0 || is.Sc > idxSc.Sc) {
			idxSc = is
			idx = i
		}
	}
	return
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

func SearchNearMK(mkArr []MapingKmerInfo, mk MapingKmerInfo, direction uint8) (nk MapingKmerInfo) {
	if direction == BACKWARD {
		for _, tk := range mkArr {
			if tk.GetPosition() < mk.GetPosition() {
				if tk.GetPosition() > nk.GetPosition() {
					nk = tk
				}
			}
		}
	} else {
		nk.SetPosition(math.MaxUint32 >> 2)
		for _, tk := range mkArr {
			if tk.GetPosition() > mk.GetPosition() {
				if tk.GetPosition() < nk.GetPosition() {
					nk = tk
				}
			}
		}
	}
	return
}

func CheckChainBoundary(mkArr []MapingKmerInfo, mk MapingKmerInfo, maxSc Score) (MapingKmerInfo, Score) {
	sc := maxSc
	lmk := SearchNearMK(mkArr, mk, BACKWARD)
	rmk := SearchNearMK(mkArr, mk, FORWARD)
	if lmk.GetMapLen() > 0 {
		if lmk.GetPosition()+uint64(lmk.GetMapLen()) > mk.GetPosition() {
			if mk.GetPosition()+uint64(mk.GetMapLen()) <= lmk.GetPosition()+uint64(lmk.GetMapLen()) {
				sc.Sc = -1
				return mk, sc
			}
			mk.SetMapLen(mk.GetMapLen() - (uint32(lmk.GetPosition()) + lmk.GetMapLen() - uint32(mk.GetPosition())))
			sc.Sc -= int(lmk.GetPosition() + uint64(lmk.GetMapLen()) - mk.GetPosition())
			if mk.GetStrand() == PLUS {
				mk.SetEPosition(mk.GetEPosition() + uint32(lmk.GetPosition()+uint64(lmk.GetMapLen())-mk.GetPosition()))
			} else {
				//
			}
			mk.SetPosition(lmk.GetPosition() + uint64(lmk.GetMapLen()))
		}
		if mk.GetStrand() == PLUS {
			if mk.GetEPosition() < lmk.GetEPosition() || mk.GetEPosition()+uint32(mk.GetMapLen()) <= lmk.GetEPosition()+uint32(lmk.GetMapLen()) {
				sc.Sc = -1
				return mk, sc
			} else if mk.GetEPosition() < lmk.GetEPosition()+uint32(lmk.GetMapLen()) {
				mk.SetMapLen(mk.GetMapLen() - (lmk.GetEPosition() + lmk.GetMapLen() - mk.GetEPosition()))
				sc.Sc -= int(lmk.GetEPosition() + uint32(lmk.GetMapLen()) - mk.GetEPosition())
				mk.SetPosition(mk.GetPosition() + uint64(lmk.GetEPosition()+lmk.GetMapLen()-mk.GetEPosition()))
				mk.SetEPosition(lmk.GetEPosition() + lmk.GetMapLen())
			}
		} else {
			if mk.GetEPosition() >= lmk.GetEPosition() {
				sc.Sc = -1
				return mk, sc
			} else if mk.GetEPosition()+uint32(mk.GetMapLen()) > lmk.GetEPosition() {
				mk.SetMapLen(mk.GetMapLen() - (mk.GetEPosition() + mk.GetMapLen() - uint32(lmk.GetPosition())))
				sc.Sc -= int(mk.GetEPosition() + mk.GetMapLen() - uint32(lmk.GetPosition()))
				mk.SetPosition(mk.GetPosition() + uint64(mk.GetEPosition()+mk.GetMapLen()-uint32(lmk.GetPosition())))
				//mk.SetEPos()
			}
		}
	}

	if rmk.GetMapLen() > 0 {
		if mk.GetPosition()+uint64(mk.GetMapLen()) > rmk.GetPosition() {
			if mk.GetPosition() >= rmk.GetPosition() {
				sc.Sc = -1
				return mk, sc
			}
			mk.SetMapLen(mk.GetMapLen() - uint32(mk.GetPosition()+uint64(mk.GetMapLen())-rmk.GetPosition()))
			sc.Sc -= int(mk.GetPosition() + uint64(mk.GetMapLen()) - rmk.GetPosition())
			if mk.GetStrand() == PLUS {

			} else {
				mk.SetEPosition(mk.GetEPosition() + uint32(mk.GetPosition()+uint64(mk.GetMapLen())-rmk.GetPosition()))
			}
		}
		if mk.GetStrand() == PLUS {
			if mk.GetEPosition() >= rmk.GetEPosition() || mk.GetEPosition()+uint32(mk.GetMapLen()) >= rmk.GetEPosition()+uint32(rmk.GetMapLen()) {
				sc.Sc = -1
				return mk, sc
			} else if mk.GetEPosition()+uint32(mk.GetMapLen()) > rmk.GetEPosition() {
				mk.SetMapLen(mk.GetMapLen() - (mk.GetEPosition() + mk.GetMapLen() - rmk.GetEPosition()))
				sc.Sc -= int(mk.GetEPosition() + mk.GetMapLen() - rmk.GetEPosition())
			}
		} else {
			if mk.GetEPosition() < rmk.GetEPosition() || mk.GetEPosition()+uint32(mk.GetMapLen()) < rmk.GetEPosition()+uint32(rmk.GetMapLen()) {
				sc.Sc = -1
				return mk, sc
			} else if mk.GetEPosition() < rmk.GetEPosition()+uint32(rmk.GetMapLen()) {
				mk.SetMapLen(mk.GetMapLen() - (rmk.GetEPosition() + rmk.GetMapLen() - mk.GetEPosition()))
				sc.Sc -= int(rmk.GetEPosition() + uint32(rmk.GetMapLen()) - mk.GetEPosition())
				mk.SetEPosition(mk.GetEPosition() + (rmk.GetEPosition() + rmk.GetMapLen() - mk.GetEPosition()))
			}
		}
	}

	return mk, sc
}

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

const ArrSize = 16

type MaxPathInfo struct {
	EID   uint32
	Base  uint32
	Score float32
	//Flag  uint8
	//Score2             uint16
	//FirstHalfArrScore  uint16 // if LenArr == ArrSize
	//SecondHalfArrScore uint16 // if LenArr == ArrSize
	//Path               []uint32
	Arr    [ArrSize]uint8 // Arr just store base index of Kb, if need mkiArr[i] = Kb[Base+Arr[i]]
	LenArr uint8
	Kb     []MapingKmerInfo // source mkiArr
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

type BoundScore struct {
	RPos  uint32
	Score int
}

func GetAnchorMKI(ki EdgeKI) (anchorMKI MapingKmerInfo, idx int) {
	for i, mki := range ki.MkInfoArr {
		if ki.Strand == mki.GetStrand() {
			if anchorMKI.GetMapLen() < mki.GetMapLen() {
				anchorMKI = mki
				idx = i
			}
		}
	}
	return
}

func CleanBoolArr(arr []bool) []bool {
	for i := range arr {
		arr[i] = false
	}
	return arr
}

func GetMKIArrByIndex(kb []MapingKmerInfo, indexArr []uint16) []MapingKmerInfo {
	ka := make([]MapingKmerInfo, len(indexArr))
	for i, idx := range indexArr {
		ka[i] = kb[idx]
	}
	return ka
}

//const Uint16MaxFlag = 1 << 15

func CleanUint8Arr(arr []uint8) {
	for i, c := range arr {
		if c != math.MaxUint8 {
			arr[i] = 0
		}
	}
	return
}

func CleanInt16Arr(arr []int16) {
	for i := range arr {
		arr[i] = 0
	}
	return
}

type IdxScore struct {
	Idx uint16
	Sc  uint16
}

//const IdxScArrLen = 10

type IdxScArr struct {
	Arr      []IdxScore
	ArrCount int
	MinScore uint16
	MinIdx   int
}

func GetVisitArr(visitArr []uint8, i int) bool {
	if i < 0 || i >= len(visitArr)*8 {
		log.Fatalf("[GetVisitArr] i:%d not allow in len(visitArr):%d array!\n", i, len(visitArr))
	}
	b := visitArr[i/8]
	idx := uint64(i) % 8
	if (b >> idx & 0x1) > 0 {
		return true
	} else {
		return false
	}
}

func SetVisitArr(visitArr []uint8, i int) {
	if i < 0 || i >= len(visitArr)*8 {
		log.Fatalf("[SetVisitArr] i:%d not allow in len(visitArr):%d array!\n", i, len(visitArr))
	}
	visitArr[i/8] |= 0x1 << (uint64(i) % 8)
}

func GetMaxScoreArr(ka []MapingKmerInfo, strand bool, SeedLen, GapCost, MaxGapLen int, SeqType uint8, el int, bsA []uint8, mkiAllFlag bool) (mp MaxPathInfo) {
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
		for j := 0; j < kaLen; j++ {
			cb := ka[j]
			if maxLen < uint32(cb.GetMapLen()) && (cb.GetStrand() == cb.GetEStrand()) == strand && !IsInIntArr(usedArr, j) {
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
		//visitArr[idx] = true
		//anchorMKI := ka[idx]
		bCArr = bCArr[:0]
		bCArr = append(bCArr, uint16(idx))
		bsA[idx] = math.MaxUint8 - 1
		BlockToNearBlocksScoreArr(ka, idx, strand, GapCost, bsA, uint16(SeedLen))
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
			BlockToNearBlocksScoreArr(ka, int(maxSc.Idx), strand, GapCost, bsA, uint16(SeedLen))
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

	// reset UsedFlag
	//for _, idx := range usedArr {
	//	ka[idx].ResetUsedFlag()
	//}
	//fmt.Printf("[GetMaxScoreArr]mp.Score: %v, mp.Arr: %v\n", mp.Score, mp.Arr[:mp.GetMapLen()Arr])
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
	return mp
}

func IndexMapingKmerInfoArr(kb []MapingKmerInfo, mk MapingKmerInfo) (idx int) {
	idx = -1
	for i, tk := range kb {
		if tk.GetPosition() == mk.GetPosition() && tk.GetMapLen() == mk.GetMapLen() && tk.EID == mk.EID {
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
	if direction == BACKWARD {
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
	if direction == FORWARD {
		distR = ReadLen - int(mk.GetPosition()+uint64(mk.GetMapLen()))
		if mk.GetStrand() == PLUS {
			distE = EdgeLen - int(mk.GetEPosition()+uint32(mk.GetMapLen()))
			for i := len(mkArr) - 1; i >= 0; i-- {
				tk := mkArr[i]
				if int(tk.GetEPosition()) > EdgeLen-(kmerlen-1) {
					shareArr = append(shareArr, tk)
				} else {
					break
				}
			}
		} else {
			distE = int(mk.GetEPosition())
			for i := len(mkArr) - 1; i >= 0; i-- {
				tk := mkArr[i]
				if tk.GetEPosition()+uint32(tk.GetMapLen()) < uint32(kmerlen-1) {
					shareArr = append(shareArr, tk)
				} else {
					break
				}
			}
		}
	} else {
		distR = int(mk.GetPosition())
		if mk.GetStrand() == PLUS {
			distE = int(mk.GetEPosition())
			for i := len(mkArr) - 1; i >= 0; i-- {
				tk := mkArr[i]
				if tk.GetEPosition()+uint32(tk.GetMapLen()) < uint32(kmerlen-1) {
					shareArr = append(shareArr, tk)
				} else {
					break
				}
			}
		} else {
			distE = EdgeLen - int(mk.GetEPosition()+uint32(mk.GetMapLen()))
			for i := len(mkArr) - 1; i >= 0; i-- {
				tk := mkArr[i]
				if int(tk.GetEPosition()) > EdgeLen-(kmerlen-1) {
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
	if direction == FORWARD {
		for _, mp := range mpArr {
			mkb := mp.Kb[mp.Base+uint32(mp.Arr[len(mp.Arr)-1])]
			if bd == 0 {
				bd = uint32(mkb.GetPosition())
				ok = true
			} else {
				if bd != uint32(mkb.GetPosition()) {
					ok = false
					break
				}
			}
		}
	} else {
		for _, mp := range mpArr {
			mka := mp.Kb[mp.Base+uint32(mp.Arr[0])]
			if bd == 0 {
				bd = uint32(mka.GetPosition())
				ok = true
			} else {
				if bd != uint32(mka.GetPosition()) {
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
		if direction == FORWARD {
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
	if direction == FORWARD {
		low, high := shareArr[0].GetPosition(), shareArr[len(shareArr)-1].GetPosition()
		idx := 0
		for _, mk := range kb {
			if mk.GetPosition() < low {
				continue
			} else if mk.GetPosition() > high {
				break
			}
			if mk.GetPosition() > shareArr[idx].GetPosition() {
				idx++
				if idx >= len(shareArr) {
					break
				}
			}
			if mk.GetPosition() == shareArr[idx].GetPosition() && mk.GetMapLen() == shareArr[idx].GetMapLen() {
				sharePat = InsertSharePat(sharePat, mk)
			}
		}
	} else {
		high, low := shareArr[0].GetPosition(), shareArr[len(shareArr)-1].GetPosition()
		idx := len(shareArr) - 1
		for _, mk := range kb {
			if mk.GetPosition() < low {
				continue
			} else if mk.GetPosition() > high {
				break
			}
			if mk.GetPosition() > shareArr[idx].GetPosition() {
				idx--
				if idx < 0 {
					break
				}
			}
			if mk.GetPosition() == shareArr[idx].GetPosition() && mk.GetMapLen() == shareArr[idx].GetMapLen() {
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
		strand := arr[0].GetStrand() == arr[0].GetEStrand()
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

func IsConnective(e, ne *DBGEdge, nd DBGNode) (ok bool) {
	if IsInComing(nd.EdgeIDIncoming, e.ID) {
		if IsInComing(nd.EdgeIDOutcoming, ne.ID) {
			ok = true
		} else {
			ok = false
		}
	} else {
		if IsInComing(nd.EdgeIDIncoming, ne.ID) {
			ok = true
		} else {
			ok = false
		}
	}
	return ok
}

func CheckDBG(sharePat [][]MapingKmerInfo, edgesArr []DBGEdge, nodesArr []DBGNode, mk MapingKmerInfo, direction uint8) [][]MapingKmerInfo {
	e := &edgesArr[mk.EID]
	strand := mk.GetStrand() == mk.GetEStrand()
	idx := 0
	if direction == FORWARD {
		if strand == PLUS {
			for _, arr := range sharePat {
				ne := &edgesArr[arr[0].EID]
				nstrd := arr[0].GetStrand()
				if !IsConnective(e, ne, nodesArr[e.EndNID]) {
					fmt.Printf("[CheckDBG] eID: %v, and neID: %v can't connective\n", e.ID, ne.ID)
					continue
				}
				if GetNextEdgeStrand2(e, ne, strand, direction) == nstrd {
					sharePat[idx] = arr
					idx++
				}
			}
		} else {
			for _, arr := range sharePat {
				ne := &edgesArr[arr[0].EID]
				nstrd := arr[0].GetStrand()
				if !IsConnective(e, ne, nodesArr[e.StartNID]) {
					fmt.Printf("[CheckDBG] eID: %v, and neID: %v can't connective\n", e.ID, ne.ID)
					continue
				}
				if GetNextEdgeStrand2(e, ne, strand, direction) == nstrd {
					sharePat[idx] = arr
					idx++
				}
			}
		}
	} else {
		if strand == PLUS {
			for _, arr := range sharePat {
				ne := &edgesArr[arr[0].EID]
				nstrd := arr[0].GetStrand()
				if !IsConnective(e, ne, nodesArr[e.StartNID]) {
					fmt.Printf("[CheckDBG] eID: %v, and neID: %v can't connective\n", e.ID, ne.ID)
					continue
				}
				if GetNextEdgeStrand2(e, ne, strand, direction) == nstrd {
					sharePat[idx] = arr
					idx++
				}
			}
		} else {
			for _, arr := range sharePat {
				ne := &edgesArr[arr[0].EID]
				nstrd := arr[0].GetStrand()
				if !IsConnective(e, ne, nodesArr[e.EndNID]) {
					fmt.Printf("[CheckDBG] eID: %v, and neID: %v can't connective\n", e.ID, ne.ID)
					continue
				}
				if GetNextEdgeStrand2(e, ne, strand, direction) == nstrd {
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

const BaseNumInByte = 4 // one bit note BaseNumInBit base

func CountBit(covArr []uint8) (sum int) {
	for _, c := range covArr {
		for j := 0; j < 8; j++ {
			sum += int(c & 0x1)
			c >>= 1
		}
	}
	return
}

func GetSumMappingKIArrLen(mkInfoArr []MapingKmerInfo, strand bool, el int) (int, int, int) {
	num := (el + 8 - 1) / 8
	covArr := make([]uint8, num)
	a, b := uint32(math.MaxInt32), uint32(0)
	for _, mki := range mkInfoArr {
		if (mki.GetStrand() == mki.GetEStrand()) != strand {
			continue
		}
		start, end := mki.GetEPosition(), mki.GetEPosition()+uint32(mki.GetMapLen())
		for j := start; j < end; j++ {
			covArr[j/8] |= (0x1 << (j % 8))
		}
		if start < a {
			a = start
		}
		if b < end {
			b = end
		}
	}
	sum := CountBit(covArr)
	//sum *= BaseNumInBit
	return sum, int(a), int(b)
}

func InitStrandInfo(kb []MapingKmerInfo, bsA []uint8, strand bool) []uint8 {
	for i, mki := range kb {
		if (mki.GetStrand() == mki.GetEStrand()) != strand {
			bsA[i] = math.MaxUint8
		} else {
			bsA[i] = 0
		}
	}
	return bsA
}

func GetMaxChainBlockArr(kb []MapingKmerInfo, maxPathInfoArr []MaxPathInfo, edgesArr []DBGEdge, nodesArr []DBGNode, SeedLen, GapCost, MaxGapLen, kmerlen int, SeqType uint8, LongEdgeMinLen int, MinChainScoreIdentityPercent, seedNumPercent float32, logfp io.Writer) []MaxPathInfo {
	//edgesKIArr := GetEdgesKIArr(kb, edgesArr, MinKmerScore, LongEdgeMinLen)
	//sort.Sort(EdgeKIArr(edgesKIArr))
	//var maxScore int
	//var maxScArr []int
	//var loopNum int
	//var cBMRIArr []ChainBlockMapLongReadInfo
	//var eIDArr []uint32
	//MinChainScoreIdentityPercent = seedNumPercent * 3 / 2
	MinChainScoreIdentityPercent = seedNumPercent
	if MinChainScoreIdentityPercent > 0.9 {
		MinChainScoreIdentityPercent = 0.9
	}
	maxPathInfoArr = maxPathInfoArr[:0]
	//EdgeMinLen := (kmerlen-1)*2 - 10
	bsShare := make([]uint8, 30) // visit just used lowest-bit of bsA element &0x1
	//var tmpArr []MaxPathInfo
	// search long read local best mapping to the edge
	totalNum, GetMaxScoreNum := 0, 0
	for i := 0; i < len(kb); {
		scP, scM := 0, 0
		eID := kb[i].EID
		totalNum++
		j := i
		//fmt.Printf("[GetMaxChainBlockArr]i:%d eID:%d\t", i, eID)
		for ; j < len(kb); j++ {
			mki := kb[j]
			if eID != mki.EID {
				break
			}
			if mki.GetStrand() == mki.GetEStrand() {
				scP += int(mki.GetMapLen())
			} else {
				scM += int(mki.GetMapLen())
			}
		}
		//fmt.Printf("[GetMaxChainBlockArr]i:%d j:%d\n", i, j)
		var strand bool
		if scP > scM {
			strand = PLUS
		} else {
			strand = MINUS
		}
		//num := MaxInt(num1, num2)
		//e := edgesArr[ki.MkInfoArr[0].EID]
		// filter
		/*if e.StartNID == e.EndNID || e.GetTwoEdgesCycleFlag() > 0 {
			continue
		}
		if len(e.Ks) < kmerlen*2-3 {
			continue
		}*/
		//PrintAddedMpArr(ki.MkInfoArr)
		el := edgesArr[eID].GetSeqLen()
		if j-i > el*5 {
			fmt.Fprintf(os.Stderr, "[GetMaxChainBlockArr]repeat edge? eID:%d el:%d EdgeType:%s kb num:%d\n", eID, el, GetEdgeType(&edgesArr[eID]), j-i)
			//continue
		}
		sc, a, b := GetSumMappingKIArrLen(kb[i:j], strand, el)
		if DebugModel {
			fmt.Fprintf(logfp, "eID:%d el:%d edgeMinNum:%d len(ka):%d mapping region:%d sc:%d\n", eID, el, edgesArr[eID].CovD, j-i, b-a, sc)
		}
		/*if sc <= SeedLen {q
			fmt.Fprintf(os.Stderr, "[GetMaxChainBlockArr]eID:%d el:%d edgeMinNum:%d len(ka):%d mapping region:%d sc:%d\n", eID, el, edgesArr[eID].CovD, j-i, b-a, sc)
			for _, mki := range kb[i:j] {
				fmt.Fprintf(os.Stderr, "%s\n", PrintMKI(mki))
			}
		}*/
		if el >= LongEdgeMinLen {
			if float32(sc) < float32(LongEdgeMinLen)*MinChainScoreIdentityPercent || b-a < LongEdgeMinLen*3/4 {
				i = j
				continue
			}
		} else {
			/*var minL int
			if el > EdgeMinLen {
				minL = el * (seedNumPercent - 5 + 67) / 100
			} else {
				minL = el * (seedNumPercent - 5 + 75) / 100
			}
			if minL > el*9/10 {
				minL = el * 9 / 10
			}*/
			if float32(sc) < float32(el)*MinChainScoreIdentityPercent || b-a < el/2 {
				i = j
				continue
			}
		}
		//PrintMKIArr(kb[i:j])
		if j-i > cap(bsShare) {
			bsShare = make([]uint8, (j - i))
		}
		bsA := bsShare[:j-i]
		bsA = InitStrandInfo(kb[i:j], bsA, strand)
		//PrintMKIArr(kb[i:j])
		//PrintMKI(mki)
		maxPI := GetMaxScoreArr(kb[i:j], strand, SeedLen, GapCost, MaxGapLen, SeqType, el, bsA, false)
		//fmt.Printf("[GetMaxChainBlockArr]maxPI.Arr:%v score:%.1f\n", maxPI.Arr[:maxPI.LenArr], maxPI.Score)
		if maxPI.LenArr < 1 {
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
		if el >= LongEdgeMinLen {
			if float32(maxPI.Score) > float32(LongEdgeMinLen)*MinChainScoreIdentityPercent && AbsInt(readRegLen-edgeRegLen) < min/10 && min > LongEdgeMinLen*3/4 {
				maxPI.Base = uint32(i)
				maxPI.EID = eID
				maxPI.Kb = kb
				maxPathInfoArr = append(maxPathInfoArr, maxPI)
				if DebugModel {
					fmt.Fprintf(logfp, "add mpi LenArr:%d score:%.1f\n", maxPI.LenArr, maxPI.Score)
				}
			}
		} else {
			/*var minL int
			if el > EdgeMinLen {
				minL = el * (seedNumPercent - 5 + 67) / 100
			} else {
				minL = el * (seedNumPercent - 5 + 75) / 100
			}
			if minL > el*9/10 {
				minL = el * 9 / 10
			}*/
			if float32(maxPI.Score) > float32(el)*MinChainScoreIdentityPercent && AbsInt(readRegLen-edgeRegLen) < min/10 && min > el/2 {
				maxPI.Base = uint32(i)
				maxPI.EID = eID
				maxPI.Kb = kb
				maxPathInfoArr = append(maxPathInfoArr, maxPI)
				if DebugModel {
					fmt.Fprintf(logfp, "add mpi LenArr:%d score:%.1f\n", maxPI.LenArr, maxPI.Score)
				}
			}
		}
		//if maxPI.EID > 0 {
		//	PrintMPI(maxPI)
		//}
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
					distR, distE, shareArr := CheckBound(pi.Arr, ReadLen, len(edgesArr[pi.Arr[0].EID].Ks), kmerlen, BACKWARD)
					if distE > kmerlen-1 || distR < kmerlen+distE || len(shareArr) == 0 {
						// add rstack
						pi.Arr = GetReverseMapingKmerInfoArr(pi.Arr)
						rstack = append(rstack, pi)
						break
					}
					sharePat := GetShareKBArr(kb, shareArr, BACKWARD)
					sharePat = CheckDBG(sharePat, edgesArr, nodesArr, shareArr[len(shareArr)-1], BACKWARD)
					if len(sharePat) == 0 {
						// add rstack
						pi.Arr = GetReverseMapingKmerInfoArr(pi.Arr)
						rstack = append(rstack, pi)
						break
					}
					pIArr := make([]MaxPathInfo, len(sharePat))
					//eIDArr := make([]uint32, len(sharePat))
					for j, a := range sharePat {
						pIArr[j] = GetNextEdgeMaxScore(kb, a, BACKWARD, SeedLen, GapCost, MaxGapLen, MinLoopNum, SeqType)
					}
					if ok, bd := GetReadBound(pIArr, BACKWARD); ok {
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
							bsArr, ok = CheckMaxScoreBound(pi.Arr[len(pi.Arr)-1].GetPosition(), pi.Score, BACKWARD, bsArr)
							if ok {
								continue
							} else {
								break˝
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
									//distR, distE := CheckBound(mpi.Arr, ReadLen, FORWARD)
									var ok bool
									bsArr, ok = CheckMaxScoreBound(mpi.Arr[len(mpi.Arr)-1].GetPosition(), mpi.Score, BACKWARD, bsArr)
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
							//distR, distE := CheckBound(mpi.Arr, ReadLen, BACKWARD)
							var ok bool
							bsArr, ok = CheckMaxScoreBound(mpi.Arr[len(mpi.Arr)-1].GetPosition(), mpi.Score, BACKWARD, bsArr)
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
		//distR, distE := CheckBound(maxPI.Arr,ReadLen, BACKWARD)

		//FORWARD
		{
			var bsArr []BoundScore
			//maxPI.DistE = distE
			//stack = append(stack, maxPI)
			//var bs BoundScore
			//bs.RPos = maxPI.Arr[len(maxPI.Arr)-1].GetPosition() + uint32(maxPI.Arr[len(maxPI.Arr)-1].GetMapLen())
			//bs.Score = maxPI.Score
			//bsArr = append(bsArr, bs)
			for len(rstack) > 0 {
				pi := rstack[len(rstack)-1]
				rstack = rstack[:len(rstack)-1]
				for {
					distR, distE, shareArr := CheckBound(pi.Arr, ReadLen, len(edgesArr[pi.Arr[len(pi.Arr)-1].EID].Ks), kmerlen, FORWARD)
					if distE > kmerlen-1 || distR < kmerlen+distE || len(shareArr) == 0 {
						// add to the maxPIArr
						maxPIArr = append(maxPIArr, pi)
						break
					}
					sharePat := GetShareKBArr(kb, shareArr, FORWARD)
					sharePat = CheckDBG(sharePat, edgesArr, nodesArr, shareArr[len(shareArr)-1], FORWARD)
					if len(sharePat) == 0 {
						// add to maxPIArr
						maxPIArr = append(maxPIArr, pi)
						break
					}
					pIArr := make([]MaxPathInfo, len(sharePat))
					for j, a := range sharePat {
						pIArr[j] = GetNextEdgeMaxScore(kb, a, FORWARD, SeedLen, GapCost, MaxGapLen, MinLoopNum, SeqType)
					}
					if ok, bd := GetReadBound(pIArr, FORWARD); ok {
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
							bsArr, ok = CheckMaxScoreBound(pi.Arr[len(pi.Arr)-1].GetPosition()+uint32(pi.Arr[len(pi.Arr)-1].GetMapLen()), pi.Score, FORWARD, bsArr)
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
									//distR, distE := CheckBound(mpi.Arr, ReadLen, FORWARD)
									var ok bool
									bsArr, ok = CheckMaxScoreBound(mpi.Arr[len(mpi.Arr)-1].GetPosition()+uint32(mpi.Arr[len(mpi.Arr)-1].GetMapLen()), mpi.Score, FORWARD, bsArr)
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
							//distR, distE := CheckBound(mpi.Arr, ReadLen, FORWARD)
							var ok bool
							bsArr, ok = CheckMaxScoreBound(mpi.Arr[len(mpi.Arr)-1].GetPosition()+uint32(mpi.Arr[len(mpi.Arr)-1].GetMapLen()), mpi.Score, FORWARD, bsArr)
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

func GetEdgeRegionSeq(e DBGEdge, mpi MaxPathInfo) (edgeRegionSeq []byte) {
	// check mpi.Arr
	mka, mkb := mpi.Kb[mpi.Base+uint32(mpi.Arr[0])], mpi.Kb[mpi.Base+uint32(mpi.Arr[len(mpi.Arr)-1])]
	if mka.EID != uint32(e.ID) || mkb.EID != uint32(e.ID) {
		log.Fatalf("[GetEdgeRegionSeq] e.ID: %v != mpi.Arr: %v\n", e.ID, mpi.Arr)
	}
	strand := mka.GetStrand() == mka.GetEStrand()
	//last := mpi.Arr[len(mpi.Arr)-1]
	//fmt.Printf("[GetEdgeRegionSeq]eID: %v, arr[0] ReadPos: %v, EdgePos: %v, strand: %v, Len: %v, arr[last] ReadPos: %v, EdgePos: %v, strand: %v, Len: %v, len(edge): %v\n", e.ID, mpi.Arr[0].GetPosition(), mpi.Arr[0].GetEPosition(), mpi.Arr[0].GetStrand(), mpi.Arr[0].GetMapLen(), last.GetPosition(), last.GetEPosition(), last.GetStrand(), last.GetMapLen(), len(e.Ks))
	if strand == PLUS {
		start, end := mka.GetEPosition(), mkb.GetEPosition()+uint32(mkb.GetMapLen())
		edgeRegionSeq = append(edgeRegionSeq, e.Ks[start:end]...)
	} else {
		start, end := mkb.GetEPosition(), mka.GetEPosition()+uint32(mka.GetMapLen())
		edgeRegionSeq = append(edgeRegionSeq, GetReverseCompByteArr(e.Ks[start:end])...)
	}
	return
}

type EdgeRegInfo struct {
	EID        uint32
	Strand     bool
	Start, End uint32
}

func RevEdgeRegInfoArr(edgeRegInfoArr []EdgeRegInfo) {
	al := len(edgeRegInfoArr)
	dl := al / 2
	for i := 0; i < dl; i++ {
		edgeRegInfoArr[i], edgeRegInfoArr[al-1-i] = edgeRegInfoArr[al-1-i], edgeRegInfoArr[i]
	}
}

func GetDirectionEdgesSeq(edgesArr []DBGEdge, nodesArr []DBGNode, kmerlen int, strand bool, path []uint32) (edgesSeq []byte) {
	extendLen := kmerlen - 1
	for _, eID := range path {
		extendLen += edgesArr[eID].GetSeqLen() - (kmerlen - 1)
	}
	edgesSeq = make([]byte, 0, extendLen)
	seqLen := 0
	var e *DBGEdge
	var startPos int
	for i, eID := range path {
		ne := &edgesArr[eID]
		if i == 0 {
			startPos = 0
		} else {
			strand = GetNextEdgeStrand2(e, ne, strand, FORWARD)
			startPos = kmerlen - 1
		}
		if strand == PLUS {
			edgesSeq = append(edgesSeq, ne.Ks[startPos:]...)
		} else {
			edgesSeq = append(edgesSeq, ne.Ks[:ne.GetSeqLen()-startPos]...)
			ReverseCompByteArr(edgesSeq[seqLen:])
		}
		e = ne
		seqLen += (ne.GetSeqLen() - startPos)
	}
	return
}

func GetDirectionEdgesSeqLen(edgesArr []DBGEdge, nodesArr []DBGNode, kmerlen int, strand bool, startPos, extendLen int, path []uint32, eID uint32) (edgesSeq []byte) {
	if len(path) == 0 || path[0] != eID {
		log.Fatalf("[GetDirectionEdgesSeqLen]path:%v, eID:%d\n", path, eID)
	}
	e := &edgesArr[eID]
	var finished bool
	edgesSeq = make([]byte, 0, extendLen)
	seqLen := 0
	if strand == PLUS {
		edgesSeq = append(edgesSeq, e.Ks[startPos:]...)
	} else {
		edgesSeq = append(edgesSeq, e.Ks[startPos:]...)
		ReverseCompByteArr(edgesSeq[seqLen:])
	}
	seqLen += len(edgesSeq)
	for i := 1; i < len(path); i++ {
		eID := path[i]
		ne := &edgesArr[eID]
		strand = GetNextEdgeStrand2(e, ne, strand, FORWARD)
		extLen := len(ne.Ks) - (int(kmerlen) - 1)
		if seqLen+extLen >= extendLen {
			finished = true
			extLen = extendLen - seqLen
		}
		if extLen <= 0 {
			//finished = true
			break
		}
		if strand == PLUS {
			edgesSeq = append(edgesSeq, ne.Ks[kmerlen-1:kmerlen-1+extLen]...)
		} else {
			edgesSeq = append(edgesSeq, ne.Ks[ne.GetSeqLen()-(kmerlen-1)-extLen:ne.GetSeqLen()-(kmerlen-1)]...)
			ReverseCompByteArr(edgesSeq[seqLen:])
		}
		if finished {
			break
		}
		e = ne
		seqLen += extLen
	}

	if len(edgesSeq) != extendLen {
		//fmt.Fprintf(os.Stderr, "[GetDirectionEdgesSeqLen] not found extendLen:%d seq, len(edgesSeq):%d \n", extendLen, len(edgesSeq))
		log.Fatalf("[GetDirectionEdgesSeqLen] not found extendLen:%d seq, len(edgesSeq):%d \n", extendLen, len(edgesSeq))
	}
	return
}

func GetDirectionEdgesSeqRegion(edgesArr []DBGEdge, nodesArr []DBGNode, direction uint8, kmerlen uint32, extendLen int, path []uint32, anchor MapingKmerInfo) (edgesSeq []byte) {
	if len(path) == 0 || uint32(path[0]) != anchor.EID {
		log.Fatalf("[GetDirectionEdgesSeq]path:%v, anchor.EID:%d\n", path, anchor.EID)
	}
	eID := anchor.EID
	e := &edgesArr[eID]
	var finished bool
	edgeRegInfoArr := make([]EdgeRegInfo, 0, len(path))
	//seqLen := 0
	var eri EdgeRegInfo
	eri.EID = path[0]
	eri.Strand = anchor.GetStrand()
	if direction == FORWARD {
		if eri.Strand == PLUS {
			eri.Start, eri.End = anchor.GetEPosition(), uint32(edgesArr[eri.EID].GetSeqLen())
		} else {
			eri.Start, eri.End = 0, anchor.GetEPosition()+uint32(anchor.GetMapLen())
		}
	} else {
		if eri.Strand == PLUS {
			eri.Start, eri.End = 0, anchor.GetEPosition()+uint32(anchor.GetMapLen())
		} else {
			eri.Start, eri.End = anchor.GetEPosition(), uint32(edgesArr[eri.EID].GetSeqLen())
		}
	}
	edgeRegInfoArr = append(edgeRegInfoArr, eri)
	seqLen := int(eri.End) - int(eri.Start)
	for i := 1; i < len(path); i++ {
		eID := path[i]
		ne := &edgesArr[eID]
		eri.Strand = GetNextEdgeStrand2(e, ne, eri.Strand, FORWARD)
		eri.EID = eID
		extLen := len(ne.Ks) - (int(kmerlen) - 1)
		if seqLen+extLen >= extendLen {
			finished = true
			extLen = extendLen - seqLen
		}
		if extLen <= 0 {
			//finished = true
			break
		}
		seqLen += extLen
		if direction == FORWARD {
			if eri.Strand == PLUS {
				eri.Start, eri.End = kmerlen-1, kmerlen-1+uint32(extLen)
			} else {
				eri.Start, eri.End = uint32(edgesArr[eri.EID].GetSeqLen())-(kmerlen-1)-uint32(extLen), uint32(edgesArr[eri.EID].GetSeqLen())-(kmerlen-1)
			}
		} else {
			if eri.Strand == PLUS {
				eri.Start, eri.End = uint32(edgesArr[eri.EID].GetSeqLen())-(kmerlen-1)-uint32(extLen), uint32(edgesArr[eri.EID].GetSeqLen())-(kmerlen-1)
			} else {
				eri.Start, eri.End = kmerlen-1, kmerlen-1+uint32(extLen)
			}
		}
		edgeRegInfoArr = append(edgeRegInfoArr, eri)
		if finished {
			break
		}
		e = ne
	}

	// get seq
	edgesSeq = make([]byte, 0, extendLen)
	if direction == BACKWARD {
		RevEdgeRegInfoArr(edgeRegInfoArr)
	}

	seqLen = 0
	for _, eri := range edgeRegInfoArr {
		if eri.Strand == PLUS {
			edgesSeq = append(edgesSeq, edgesArr[eri.EID].Ks[eri.Start:eri.End]...)
		} else {
			edgesSeq = append(edgesSeq, edgesArr[eri.EID].Ks[eri.Start:eri.End]...)
			ReverseCompByteArr(edgesSeq[seqLen:])
		}
		seqLen += int(eri.End - eri.Start)
	}

	if len(edgesSeq) != extendLen {
		//fmt.Fprintf(os.Stderr, "[GetDirectionEdgesSeq] not found extendLen:%d seq, len(edgesSeq):%d \n", extendLen, len(edgesSeq))
		log.Fatalf("[GetDirectionEdgesSeq] not found extendLen:%d seq, len(edgesSeq):%d \n", extendLen, len(edgesSeq))
	}
	return
}

/*func GetDirectionEdgesSeq(edgesArr []DBGEdge, nodesArr []DBGNode, direction uint8, kb []MapingKmerInfo, kmerlen, extendLen int) (edgesSeq []byte) {
	mkS := kb[0]
	mkE := kb[len(kb)-1]
	fmt.Printf("[GetDirectionEdgesSeq] mkE EID: %v, EdgePos: %v,Len: %v, strand: %v\n", mkE.EID, mkE.GetEPosition(), mkE.GetMapLen(), mkE.GetStrand())
	if direction == FORWARD {
		// start edge
		e := edgesArr[mkS.EID]
		var finished bool
		if mkS.GetStrand() == PLUS {
			edgesSeq = append(edgesSeq, e.Ks[mkS.GetEPosition():]...)
		} else {
			edgesSeq = append(edgesSeq, GetReverseCompByteArr(e.Ks[:mkS.GetEPosition()+uint32(mkS.GetMapLen())])...)
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
			extLen := len(e.Ks) - (kmerlen - 1)
			finished = false
			if len(edgesSeq)+(len(e.Ks)-(kmerlen-1)) >= extendLen {
				finished = true
				extLen = extendLen - len(edgesSeq)
			}
			if extLen <= 0 {
				finished = true
				break
			}
			if mk.GetStrand() == PLUS {
				edgesSeq = append(edgesSeq, e.Ks[kmerlen-1:kmerlen-1+extLen]...)
			} else {
				edgesSeq = append(edgesSeq, GetReverseCompByteArr(e.Ks[len(e.Ks)-(kmerlen-1)-extLen:len(e.Ks)-(kmerlen-1)])...)
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
		if extLen > len(e.Ks)-(kmerlen-1) {
			extLen = len(e.Ks) - (kmerlen - 1)
		}
		if mkE.GetStrand() == PLUS {
			if int(mkE.GetEPosition())+int(mkE.GetMapLen()) < kmerlen-1 {
				dl := kmerlen - 1 - (int(mkE.GetEPosition()) + int(mkE.GetMapLen()))
				edgesSeq = edgesSeq[:len(edgesSeq)-dl]
			} else {

			edgesSeq = append(edgesSeq, e.Ks[kmerlen-1:kmerlen-1+extLen]...)
			//}
		} else {
			if int(mkE.GetEPosition()) > len(e.Ks)-(kmerlen-1) {
				dl := int(mkE.GetEPosition()) - (len(e.Ks) - (kmerlen - 1))
				edgesSeq = edgesSeq[:len(edgesSeq)-dl]
			} else {

			edgesSeq = append(edgesSeq, GetReverseCompByteArr(e.Ks[len(e.Ks)-(kmerlen-1)-extLen:len(e.Ks)-(kmerlen-1)])...)
			//}
		}
	} else {
		// start edge
		e := edgesArr[mkS.EID]
		var finished bool
		if mkS.GetStrand() == PLUS {
			edgesSeq = append(edgesSeq, GetReverseByteArr(e.Ks[:mkS.GetEPosition()+uint32(mkS.GetMapLen())])...)
		} else {
			edgesSeq = append(edgesSeq, GetCompByteArr(e.Ks[mkS.GetEPosition():])...)
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
			extLen := len(e.Ks) - (kmerlen - 1)
			finished = false
			if len(edgesSeq)+(len(e.Ks)-(kmerlen-1)) >= extendLen {
				finished = true
				extLen = extendLen - len(edgesSeq)
			}
			if extLen <= 0 {
				finished = true
				break
			}
			if mkS.GetStrand() == PLUS {
				edgesSeq = append(edgesSeq, GetReverseByteArr(e.Ks[len(e.Ks)-(kmerlen-1)-extLen:len(e.Ks)-(kmerlen-1)])...)
			} else {
				edgesSeq = append(edgesSeq, GetCompByteArr(e.Ks[(kmerlen-1):(kmerlen-1)+extLen])...)
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
		if extLen > len(e.Ks)-(kmerlen-1) {
			extLen = len(e.Ks) - (kmerlen - 1)
		}
		if mkE.GetStrand() == PLUS {
			if int(mkE.GetEPosition()) > len(e.Ks)-(kmerlen-1) {
				dl := int(mkE.GetEPosition()) - (len(e.Ks) - (kmerlen - 1))
				edgesSeq = edgesSeq[:len(edgesSeq)-dl]
			} else {
			edgesSeq = append(edgesSeq, GetReverseByteArr(e.Ks[len(e.Ks)-(kmerlen-1)-extLen:len(e.Ks)-(kmerlen-1)])...)
			//}
		} else {
			if (kmerlen - 1) > int(mkE.GetEPosition())+int(mkE.GetMapLen()) {
				dl := (kmerlen - 1) - (int(mkE.GetEPosition()) + int(mkE.GetMapLen()))
				edgesSeq = edgesSeq[:len(edgesSeq)-dl]
			} else {
			edgesSeq = append(edgesSeq, GetCompByteArr(e.Ks[(kmerlen-1):(kmerlen-1)+extLen])...)
			//}
		}
	}
	return
}*/

func GetEdgesSeq(edgesArr []DBGEdge, nodesArr []DBGNode, kb []MapingKmerInfo, kmerlen int) (edgesSeq []byte, path []uint32, chainArr []Chain) {
	chainArr = make([]Chain, len(kb))
	mkS := kb[0]
	mkE := kb[len(kb)-1]
	var startY, startX uint32
	// start Anchor
	{
		e := edgesArr[mkS.EID]
		if mkS.GetStrand() == PLUS {
			edgesSeq = append(edgesSeq, e.Ks[mkS.GetEPosition():]...)
			startX = mkS.GetEPosition()
		} else {
			edgesSeq = append(edgesSeq, GetReverseCompByteArr(e.Ks[:mkS.GetEPosition()+uint32(mkS.GetMapLen())])...)
			startX = mkS.GetEPosition() + uint32(mkS.GetMapLen())
		}
		path = append(path, e.ID)
		chainArr[0].Len = mkS.GetMapLen()
		startY = uint32(mkS.GetPosition())
	}
	// middle
	{
		// determine boundary
		i := 1
		for ; i < len(kb)-1; i++ {
			mk := kb[i]
			if mk.EID == uint32(path[len(path)-1]) {
				if mk.GetStrand() == PLUS {
					chainArr[i].X = mk.GetEPosition() - startX
				} else {
					chainArr[i].X = startX - (mk.GetEPosition() + uint32(mk.GetMapLen()))
				}
				chainArr[i].Len = mk.GetMapLen()
				chainArr[i].Y = uint32(mk.GetPosition()) - startY
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
				if mk.GetStrand() == PLUS {
					chainArr[i].X = uint32(len(edgesSeq)) + (mk.GetEPosition() - (uint32(kmerlen) - 1))
				} else {
					chainArr[i].X = startX - (mk.GetEPosition() + uint32(mk.GetMapLen()))
				}
				chainArr[i].Len = mk.GetMapLen()
				chainArr[i].Y = uint32(mk.GetPosition()) - startY
			}
			if mk.EID == uint32(path[len(path)-1]) {
				continue
			}
			e := edgesArr[mk.EID]
			if mk.GetStrand() == PLUS {
				edgesSeq = append(edgesSeq, e.Ks[kmerlen-1:]...)
			} else {
				edgesSeq = append(edgesSeq, GetReverseCompByteArr(e.Ks[:len(e.Ks)-(kmerlen-1)])...)
			}
			path = append(path, e.ID)
		}
	}
	// end Anthor
	{
		e := edgesArr[mkE.EID]
		if mkE.GetStrand() == PLUS {
			edgesSeq = append(edgesSeq, e.Ks[kmerlen-1:mkE.GetEPosition()+uint32(mkE.GetMapLen())]...)
		} else {
			edgesSeq = append(edgesSeq, GetReverseCompByteArr(e.Ks[mkE.GetEPosition():len(e.Ks)-(kmerlen-1)])...)
		}
		path = append(path, e.ID)
	}

	return
}
func GetReadSeq(seq []byte, start, end uint32) (readSeq []byte) {
	readSeq = seq[start:end]
	return
}

type ExtendEdgePathInfo struct {
	Last *ExtendEdgePathInfo
	EID  uint32
}

type ExtendPathInfo struct {
	MPI       MaxPathInfo
	ExtendEPI ExtendEdgePathInfo
	Path      []uint32
	Path1     []uint32
	Flag      uint8
	//Strand              bool
	//Coming              bool
	ExtendLen, FlankLen uint32 // FlankLen note this edge unseed flank boundary length
	Score               uint32
}

func (epi *ExtendPathInfo) GetStrandFlag() bool {
	return (epi.Flag & 0x1) > 0
}

func (epi *ExtendPathInfo) SetStrandFlag(strand bool) {
	if strand == PLUS {
		epi.Flag = epi.Flag | 0x1
	} else {
		epi.Flag = epi.Flag & (0xFF - 0x1)
	}
}

func (epi *ExtendPathInfo) GetComingFlag() uint8 {
	if (epi.Flag & 0x2) > 0 {
		return FORWARD
	} else {
		return BACKWARD
	}
}

func (epi *ExtendPathInfo) SetComingFlag(direction uint8) {
	if direction == FORWARD {
		epi.Flag = epi.Flag | 0x2
	} else {
		epi.Flag = epi.Flag & (0xFF - 0x2)
	}
}

func GetEdgeComing(e DBGEdge, direction uint8, strand bool, nodesArr []DBGNode, coming uint8) uint8 {
	if e.StartNID > 0 && e.StartNID == e.EndNID {
		nd := nodesArr[e.StartNID]
		if IsInComing(nd.EdgeIDIncoming, e.ID) && IsInComing(nd.EdgeIDOutcoming, e.ID) {
			return coming
		} else {
			if coming > 0 {
				return 0
			} else {
				return 1
			}
		}
	}
	var nd DBGNode
	if strand == PLUS {
		if direction == FORWARD {
			if e.EndNID > 0 {
				nd = nodesArr[e.EndNID]
			}
		} else {
			if e.StartNID > 0 {
				nd = nodesArr[e.StartNID]
			}
		}
	} else {
		if direction == FORWARD {
			if e.StartNID > 0 {
				nd = nodesArr[e.StartNID]
			}
		} else {
			if e.EndNID > 0 {
				nd = nodesArr[e.EndNID]
			}
		}
	}
	if IsInComing(nd.EdgeIDIncoming, e.ID) {
		coming = FORWARD
	} else {
		coming = BACKWARD
	}
	return coming
}

func GetNearDirectionEdgeIDArr(e DBGEdge, direction uint8, strand bool, coming uint8, nodesArr []DBGNode) (eArr []uint32) {
	var nd DBGNode
	if strand == PLUS {
		if direction == FORWARD {
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
		if direction == FORWARD {
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
	c := (coming > 0)
	eArr = GetNearEdgeIDArr(nd, e.ID, c)
	return
}

func GetMappingDBGExtendInfo(mpi, extendMPI MaxPathInfo, score, flankLen uint32, strand bool, direction uint8, eLen, distanceE, kmerlen int, edgesArr []DBGEdge) (uint32, uint32) {
	same := false
	var fl int
	//var last2MKI MapingKmerInfo
	//fl = flankLen
	newScore := float32(score)
	lastMKI := mpi.Kb[mpi.Base+uint32(mpi.Arr[len(mpi.Arr)-1])]
	//fmt.Printf("[GetMappingDBGExtendInfo] : %v dY: %v, lastMKI: %v, empi: %v\n", dX, dY, PrintMKI(lastMKI), PrintMKI(empi))
	if direction == FORWARD {
		for i, idx := range extendMPI.Arr {
			empi := extendMPI.Kb[extendMPI.Base+uint32(idx)]
			if empi.GetPosition() > lastMKI.GetPosition() {
				break
			}
			/*if mpi.EID != lastEID {
				break
			}*/

			if empi.GetPosition() == lastMKI.GetPosition() && empi.GetMapLen() > lastMKI.GetMapLen() {
				_, dY := GetMKIDistance(lastMKI, empi, edgesArr, distanceE)
				if dY != 0-int(lastMKI.GetMapLen()) {
					//fmt.Printf("[GetMappingDBGExtendInfo] dX: %v dY: %v, lastMKI: %v, empi: %v\n", dX, dY, PrintMKI(lastMKI), PrintMKI(empi))
					continue
				}
				//score += int(empi.GetMapLen()) - int(lastMKI.GetMapLen())
				sc := float32(lastMKI.GetMapLen())
				last := empi
				for x := i - 1; x >= 0; x-- {
					mk := extendMPI.Kb[extendMPI.Base+uint32(extendMPI.Arr[x])]
					sc += GetTwoBlockScore2(last, mk, 1)
					last = mk
				}
				newScore = float32(score) + extendMPI.Score - sc
				same = true
				break
			}
		}

		if !same {
			var i int
			var empi MapingKmerInfo
			for x, idx := range extendMPI.Arr {
				empi = extendMPI.Kb[extendMPI.Base+uint32(idx)]
				i = x
				if empi.GetPosition() >= lastMKI.GetPosition()+uint64(lastMKI.GetMapLen()) {
					break
				}
			}
			if i == 0 {
				dX, dY := GetMKIDistance(lastMKI, empi, edgesArr, distanceE)
				min := MinInt(dX, dY)
				if min < 0-int(lastMKI.GetMapLen()) {
					//fmt.Printf("[GetMappingDBGExtendInfo] dX: %v dY: %v, cen: %v, mk: %v\n", dX, dY, PrintMKI(lastMKI), PrintMKI(extendMpArr[i]))
					//sc := GetMappingKIArrScore(extendMpArr, 1, edgesArr, 0)
					//score += sc * 9 / 10
				} else {
					/*if min < 0 {
						extendMpArr[i].SetRPos(extendMpArr[i].GetPosition() + uint32(AbsInt(min)))
						if extendMpArr[i].GetStrand() == PLUS {
							extendMpArr[i].SetEPos(extendMpArr[i].GetEPosition() + uint32(AbsInt(min)))
						}
						extendMpArr[i].GetMapLen() -= uint16(AbsInt(min))
					}*/

					gapLen := AbsInt(dX - dY)
					if gapLen < min/4+20 {
						//newArr = append(newArr, extendMpArr[i:]...)
						sc := GetTwoBlockScore(lastMKI, empi, 1, edgesArr, distanceE)
						sc -= int(empi.GetMapLen())
						newScore = float32(score) + extendMPI.Score + float32(sc)
						/*for x := len(mpArr) - 1; x < len(newArr)-1; x++ {
							score += GetTwoBlockScore(newArr[x], newArr[x+1], 1, edgesArr, distanceE)
						}*/
					}
				}
			}
		}

		if newScore > float32(score) {
			newLast := extendMPI.Kb[extendMPI.Base+uint32(extendMPI.Arr[len(extendMPI.Arr)-1])]
			if strand == PLUS {
				fl = eLen - (int(newLast.GetEPosition()) + int(newLast.GetMapLen()))
			} else {
				fl = int(newLast.GetEPosition())
			}
		} else {
			fl = int(flankLen) + eLen - (kmerlen - 1)
		}

	} else { // BACKWARD
		for i, idx := range extendMPI.Arr {
			empi := extendMPI.Kb[extendMPI.Base+uint32(idx)]
			if empi.GetPosition()+uint64(empi.GetMapLen()) < lastMKI.GetPosition()+uint64(lastMKI.GetMapLen()) {
				break
			}

			if empi.GetPosition() < lastMKI.GetPosition() && empi.GetPosition()+uint64(empi.GetMapLen()) == lastMKI.GetPosition()+uint64(lastMKI.GetMapLen()) {
				_, dY := GetMKIDistance(lastMKI, empi, edgesArr, distanceE)
				if dY != 0-int(lastMKI.GetMapLen()) {
					//fmt.Printf("[GetMappingDBGExtendInfo] dX: %v dY: %v, lastMKI: %v, empi: %v\n", dX, dY, PrintMKI(lastMKI), PrintMKI(empi))
					continue
				}
				sc := float32(lastMKI.GetMapLen())
				last := empi
				for x := i - 1; x >= 0; x-- {
					mk := extendMPI.Kb[extendMPI.Base+uint32(extendMPI.Arr[x])]
					sc += GetTwoBlockScore2(last, mk, 1)
					last = mk
				}
				newScore = float32(score) + extendMPI.Score - sc
				same = true
				break
			}
		}
		if !same {
			var i int
			var empi MapingKmerInfo
			for x, idx := range extendMPI.Arr {
				empi = extendMPI.Kb[extendMPI.Base+uint32(idx)]
				i = x
				if empi.GetPosition()+uint64(empi.GetMapLen()) <= lastMKI.GetPosition() {
					break
				}
			}
			if i == 0 {
				dX, dY := GetMKIDistance(lastMKI, empi, edgesArr, distanceE)
				min := MinInt(dX, dY)
				if min < 0-int(lastMKI.GetMapLen()) {
					//fmt.Printf("[GetMappingDBGExtendInfo] dX: %v dY: %v, cen: %v, mk: %v\n", dX, dY, PrintMKI(lastMKI), PrintMKI(extendMpArr[i]))
					//sc := GetMappingKIArrScore(extendMpArr, 1, edgesArr, 0)
					//score += sc * 9 / 10
				} else {
					/*if min < 0 {
						if extendMpArr[i].GetStrand() == MINUS {
							extendMpArr[i].SetEPos(extendMpArr[i].GetEPosition() + uint32(AbsInt(min)))
						}
						extendMpArr[i].GetMapLen() -= uint16(AbsInt(min))
					}*/

					gapLen := AbsInt(dX - dY)
					if gapLen < min/4+20 {
						sc := GetTwoBlockScore(lastMKI, empi, 1, edgesArr, distanceE)
						sc -= int(empi.GetMapLen())
						newScore = float32(score) + extendMPI.Score + float32(sc)
						//newArr = append(newArr, extendMpArr[i:]...)
						/*for x := len(mpArr) - 1; x < len(newArr)-1; x++ {
							score += GetTwoBlockScore(newArr[x], newArr[x+1], 1, edgesArr, distanceE)
						}*/
					}
				}
			}
		}

		if newScore-float32(score) > 0.5 {
			newLast := extendMPI.Kb[extendMPI.Base+uint32(extendMPI.Arr[len(extendMPI.Arr)-1])]
			if strand == PLUS {
				fl = int(newLast.GetEPosition())
			} else {
				fl = eLen - (int(newLast.GetEPosition()) + int(newLast.GetMapLen()))
			}
		} else {
			fl = int(flankLen) + eLen - (kmerlen - 1)
		}
	}
	if newScore < 0 {
		newScore = 0
	}
	return uint32(newScore), uint32(fl)
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
	if len(kb) < 1 {
		return idx
	}
	a, b := int(kb[0].GetPosition()), int(kb[len(kb)-1].GetPosition())
	mR := mki.GetPosition()
	mE := mki.GetEPosition()
	pos := len(kb) * int(mR) / (b - a)
	if pos < 0 {
		pos = 0
	} else if pos >= len(kb) {
		pos = len(kb) - 1
	}
	for j := pos; j >= 0; j-- {
		if kb[j].GetPosition()+uint64(kb[j].GetMapLen()) >= mR-uint64(mki.GetMapLen()) {
			pos = j
		} else {
			break
		}
	}
	if kb[pos].GetPosition() <= mR {
		for i := pos; i < len(kb); i++ {
			c := kb[i]
			if c.EID == mki.EID && c.GetPosition() <= mR && mR+uint64(mki.GetMapLen()) <= c.GetPosition()+uint64(c.GetMapLen()) {
				if c.GetStrand() == mki.GetStrand() && c.GetEStrand() == mki.GetEStrand() {
					if c.GetEPosition() <= mE && mE+uint32(mki.GetMapLen()) <= c.GetEPosition()+uint32(c.GetMapLen()) {
						idx = i
						break
					}
				}
			}
		}
		if idx < 0 {
			for i := pos - 1; i >= 0; i-- {
				c := kb[i]
				if c.EID == mki.EID && c.GetPosition() <= mR && mR+uint64(mki.GetMapLen()) <= c.GetPosition()+uint64(c.GetMapLen()) {
					if c.GetStrand() == mki.GetStrand() && c.GetEStrand() == mki.GetEStrand() {
						if c.GetEPosition() <= mE && mE+mki.GetMapLen() <= c.GetEPosition()+c.GetMapLen() {
							idx = i
							break
						}
					}
				}
			}
		}
	} else {
		for i := pos; i >= 0; i-- {
			c := kb[i]
			cR := c.GetPosition()
			cE := c.GetEPosition()
			if c.EID == mki.EID && cR <= mR && mR+uint64(mki.GetMapLen()) <= cR+uint64(c.GetMapLen()) {
				if c.GetStrand() == mki.GetStrand() && c.GetEStrand() == mki.GetEStrand() {
					if cE <= mE && mE+uint32(mki.GetMapLen()) <= cE+uint32(c.GetMapLen()) {
						idx = i
						break
					}
				}
			}
		}

		if idx < 0 {
			for i := pos + 1; i < len(kb); i++ {
				c := kb[i]
				if c.EID == mki.EID && c.GetPosition() <= mR && mR+uint64(mki.GetMapLen()) <= c.GetPosition()+uint64(c.GetMapLen()) {
					if c.GetStrand() == mki.GetStrand() && c.GetEStrand() == mki.GetEStrand() {
						if c.GetEPosition() <= mE && mE+uint32(mki.GetMapLen()) <= c.GetEPosition()+uint32(c.GetMapLen()) {
							idx = i
							break
						}
					}
				}
			}
		}
	}
	if idx < 0 {
		log.Fatalf("[LocalMapingKmerInfoArr] not found mki: %v in kb\n", mki)
	}
	return idx
}

func GetEdgeMKIArr(kb []MapingKmerInfo, idx int, eID uint32, direction uint8, strand bool, start, end int) (mpArr []MapingKmerInfo) {
	//mpArr = append(mpArr, cen)
	kl := AbsInt(end - start)
	mpArr = make([]MapingKmerInfo, 0, kl/40+2)
	//fmt.Printf("[GetEdgeMKIArr]eID: %v, idx: %v,cen RPos: %v, direction: %v, strand: %v, distance: %v\n", eID, idx, cen.GetPosition(), direction, strand, end)
	if direction == FORWARD {
		for i := idx - 1; i >= 0; i-- {
			mk := kb[i]
			if int(mk.GetPosition()) < start {
				break
			}
			if mk.EID != eID || mk.GetStrand() != strand {
				continue
			}
			mpArr = append(mpArr, mk)
		}
		RevMappingKmerInfoArr(mpArr)

		for i := idx + 1; i < len(kb); i++ {
			mk := kb[i]
			if int(mk.GetPosition())+int(mk.GetMapLen()) > end {
				break
			}
			//fmt.Printf("[GetEdgeMKIArr] kb[%v]: %v, RPos: %v, strand: %v, EID: %v\n", i, mk, mk.GetPosition(), mk.GetStrand(), mk.EID)
			if mk.EID != eID || mk.GetStrand() != strand {
				continue
			}
			mpArr = append(mpArr, mk)
		}
	} else {
		for i := idx - 1; i > 0; i-- {
			mk := kb[i]
			//fmt.Printf("[GetEdgeMKIArr] kb[%v]: %v, RPos: %v, EID: %v\n", i, mk, mk.GetPosition(), mk.EID)
			if int(mk.GetPosition()) < end {
				break
			}
			if mk.EID != eID || mk.GetStrand() != strand {
				continue
			}
			mpArr = append(mpArr, mk)
		}
		RevMappingKmerInfoArr(mpArr)

		for i := idx + 1; i < len(kb); i++ {
			mk := kb[i]
			if int(mk.GetPosition())+int(mk.GetMapLen()) > start {
				break
			}
			if mk.EID != eID || mk.GetStrand() != strand {
				continue
			}
			mpArr = append(mpArr, mk)
		}

	}
	//sort.Sort(MapingKmerInfoArr(mpArr))
	return
}

func GetEdgeMKIArr2(kb []MapingKmerInfo, kbIdx int, eID uint32, strand bool, start, end int) (mpArr []MapingKmerInfo, indexArr []uint16) {
	//mpArr = append(mpArr, cen)
	kl := AbsInt(end - start)
	mpArr = make([]MapingKmerInfo, 0, kl/40+2)
	indexArr = make([]uint16, 0, kl/40+2)
	//fmt.Printf("[GetEdgeMKIArr]eID: %v, idx: %v,cen RPos: %v, direction: %v, strand: %v, distance: %v\n", eID, idx, cen.GetPosition(), direction, strand, end)
	ka := kb[kbIdx:]
	for i, mk := range ka {
		if int(mk.GetPosition()) < start {
			continue
		}
		if mk.EID != eID || int(mk.GetPosition())+int(mk.GetMapLen()) > end {
			break
		}
		//fmt.Printf("[GetEdgeMKIArr] kb[%v]: %v, RPos: %v, strand: %v, EID: %v\n", i, mk, mk.GetPosition(), mk.GetStrand(), mk.EID)
		if mk.GetStrand() != strand {
			continue
		}
		mpArr = append(mpArr, mk)
		indexArr = append(indexArr, uint16(i))
	}
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

func RevUint8Arr(arr [ArrSize]uint8, al int) {
	dl := al / 2
	for i := 0; i < dl; i++ {
		arr[i], arr[al-1-i] = arr[al-1-i], arr[i]
	}
	return
}

func GetRevUint8Arr(arr []uint8, al int) []uint8 {
	ra := make([]uint8, al)
	for i, c := range arr {
		ra[al-1-i] = c
	}
	return ra
}

func PrintMKI(mki MapingKmerInfo) string {
	return fmt.Sprintf("mki EID:%d RPos:%d EPos:%d Strand:%v Len:%d", mki.EID, mki.GetPosition(), mki.GetEPosition(), mki.GetStrand() == mki.GetEStrand(), mki.GetMapLen())
}

func GetMKI2EdgeDistance(extendEPI ExtendEdgePathInfo, neID uint32, lastMki MapingKmerInfo, direction uint8, edgesArr []DBGEdge, kmerlen int) (distance int) {
	if direction == FORWARD {
		if lastMki.EID == uint32(neID) {
			if lastMki.GetStrand() == PLUS {
				distance = 0 - (int(lastMki.GetEPosition()) + int(lastMki.GetMapLen()))
			} else {
				distance = int(lastMki.GetEPosition()) + int(lastMki.GetMapLen()) - len(edgesArr[lastMki.EID].Ks)
			}
		} else {
			//idx := GetIndexDbgMaxIntArrDirection(path, uint32(neID), BACKWARD)
			epi := extendEPI
			for epi.EID > 0 {
				eID := epi.EID
				if uint32(eID) != lastMki.EID {
					distance += len(edgesArr[eID].Ks) - (kmerlen - 1)
				} else {
					if lastMki.GetStrand() == PLUS {
						distance += len(edgesArr[eID].Ks) - (kmerlen - 1) - (int(lastMki.GetEPosition()) + int(lastMki.GetMapLen()))
					} else {
						distance += int(lastMki.GetEPosition()) - (kmerlen - 1)
					}
					break
				}
				if epi.Last != nil {
					epi = *epi.Last
				} else {
					epi.EID = 0
					epi.Last = nil
				}
			}
		}
	} else { // BACKWARD
		if lastMki.EID == uint32(neID) {
			if lastMki.GetStrand() == PLUS {
				distance = (int(lastMki.GetEPosition()) + int(lastMki.GetMapLen())) - len(edgesArr[lastMki.EID].Ks)
			} else {
				distance = 0 - (int(lastMki.GetEPosition()) + int(lastMki.GetMapLen()))
			}
		} else {
			epi := extendEPI
			for epi.EID > 0 {
				eID := epi.EID
				if uint32(eID) != lastMki.EID {
					distance += len(edgesArr[eID].Ks) - (kmerlen - 1)
				} else {
					if lastMki.GetStrand() == PLUS {
						distance += int(lastMki.GetEPosition()) - (kmerlen - 1)
					} else {
						distance += len(edgesArr[eID].Ks) - (kmerlen - 1) - (int(lastMki.GetEPosition()) + int(lastMki.GetMapLen()))
					}
					break
				}
				if epi.Last != nil {
					epi = *epi.Last
				} else {
					epi.EID = 0
					epi.Last = nil
				}
			}
		}
	}
	//fmt.Printf("[GetMKI2EdgeDistance] lastMki: %v\t\t: neID: %v, distance: %v, path: %v\n", PrintMKI(lastMki), neID, distance, path)

	return
}

func GetEdgeAnchorMKI(ekb []MapingKmerInfo) (max, count int) {
	count = 0
	max = math.MinInt32
	for _, mki := range ekb {
		if int(mki.GetMapLen()) > max {
			max = int(mki.GetMapLen())
			count++
		} else if int(mki.GetMapLen()) == max {
			count++
		}
	}
	return
}

/*func GetExtendMapingKmerInfoArr(ekb []MapingKmerInfo, ekbIdx int, strand bool, direction uint8, edgesArr []DBGEdge, SeedLen uint16) (extendMpArr []MapingKmerInfo) {
	anchor := ekb[ekbIdx]
	visitArr := make([]bool, len(ekb))
	bsA := make([]int, len(ekb))
	//fmt.Printf("[ExtendPathChain] visitArr: %v\n\tekbIdx: %v\n", visitArr, ekbIdx)
	visitArr[ekbIdx] = true
	extendMpArr = append(extendMpArr, anchor)
	bsA[ekbIdx] = int(anchor.GetMapLen())
	BlockToNearBlocksScoreArr(ekb, ekbIdx, strand, visitArr, 1, bsA, edgesArr, SeedLen)
	for {
		maxSc := GetMaxScoreFromBaSA(bsA, visitArr)
		if maxSc.Sc <= 0 {
			break
		}
		//fmt.Printf("[ExtendPathChain] added mk[%v]: %v maxSc.Sc: %v\n", maxSc.Idx, PrintMKI(ekb[maxSc.Idx]), maxSc.Sc)
		visitArr[maxSc.Idx] = true
		bsA[maxSc.Idx] = int(ekb[maxSc.Idx].GetMapLen())
		extendMpArr = append(extendMpArr, ekb[maxSc.Idx])
		BlockToNearBlocksScoreArr(ekb, maxSc.Idx, strand, visitArr, 1, bsA, edgesArr, SeedLen)
	}
	sort.Sort(MapingKmerInfoArr(extendMpArr))
	return
}*/

func ExtendPathChain(mpi MaxPathInfo, extendEPI ExtendEdgePathInfo, flankLen uint32, kmerlen int, direction uint8, ne DBGEdge, kbIdx int, strand bool, score uint32, edgesArr []DBGEdge, nodesArr []DBGNode, SeedLen uint16, maxPathInfoArr []MaxPathInfo, LongEdgeMinLen, MinChainScoreIdentityPercent int) (nmpi MaxPathInfo, newScore, fl uint32) {
	//extendMpArr = append(extendMpArr, mpArr...)
	//fmt.Printf("[ExtendPathChain] mpArr: \n")
	//PrintAddedMpArr(mpArr)
	kb := mpi.Kb
	newScore = score
	var extendMPI MaxPathInfo
	lastMki := kb[mpi.Base+uint32(mpi.Arr[len(mpi.Arr)-1])]
	/*var distance int
	leID := path[len(path)-2]
	if len(path) > 1 && edgesArr[leID].GetTwoEdgesCycleFlag() > 0 {
		//lastPathPos := IndexLastEID(path, uint32(lastMki.EID))
		tl := SumTwoEdgeCycleLen(edgesArr[leID], edgesArr, nodesArr)
		distance = tl - 2*(kmerlen-1)
	}*/
	var start, end int
	SD := 50
	if direction == FORWARD {
		start = int(lastMki.GetPosition()) + int(lastMki.GetMapLen()) + int(flankLen)*9/10 - (int(kmerlen) - 1) - SD
		end = start + ne.GetSeqLen()*7/5 + 2*SD
	} else {
		start = int(lastMki.GetPosition()) - int(flankLen)*9/10 + (int(kmerlen) - 1) + SD
		end = start - ne.GetSeqLen()*7/5 - 2*SD
	}

	if start > end {
		start, end = end, start
	}

	mpiArrOk := false
	mpiArrIndex := GetIndexMaxPathInfoArr(maxPathInfoArr, uint32(ne.ID))
	if ne.GetTwoEdgesCycleFlag() == 0 && ne.StartNID != ne.EndNID && mpiArrIndex >= 0 {
		//extendMpArr = make([]MapingKmerInfo, len(maxPathInfoArr[mpiArrIndex].Arr))
		//copy(extendMpArr, maxPathInfoArr[mpiArrIndex].Arr)
		t := maxPathInfoArr[mpiArrIndex]
		mka, mkb := t.Kb[t.Base+uint32(t.Arr[0])], t.Kb[t.Base+uint32(t.Arr[len(t.Arr)-1])]
		if int(mka.GetPosition()) > start && int(mkb.GetPosition())+int(mkb.GetMapLen()) < end {
			extendMPI = t
			mpiArrOk = true
			if direction == BACKWARD {
				RevUint8Arr(extendMPI.Arr, int(extendMPI.LenArr))
			}
		}
	}

	if !mpiArrOk {
		//idx := LocalMapingKmerInfoArr(kb, lastMki)
		//ekb := GetEdgeMKIArr(kb, idx, uint32(ne.ID), direction, strand, flankLen+(len(ne.Ks)-(kmerlen-1))*7/5+50)

		//fmt.Printf("[ExtendPathChain]start:%d end:%d distance:%d lastMKI:%v\n", start, end, distance, PrintMKI(lastMki))
		ekb, indexArr := GetEdgeMKIArr2(kb, kbIdx, uint32(ne.ID), strand, start, end)
		if len(ekb) == 0 {
			nmpi = mpi
			//copy(nmpi.Arr[:mpi.GetMapLen()Arr], mpi.Arr[:mpi.GetMapLen()Arr])
			//nmpi = append(newArr, mpArr...)
			fl = flankLen + uint32(len(ne.Ks)-(kmerlen-1))
			return
		}
		//fmt.Printf("[ExtendPathChain]indexArr:%v\n", indexArr)
		bsA := make([]uint8, len(ekb)) // visit just used lowest-bit of bsA element &0x1
		maxPI := GetMaxScoreArr(ekb, strand, int(SeedLen), 1, 500, 1, ne.GetSeqLen(), bsA, false)
		el := ne.GetSeqLen()
		if el > LongEdgeMinLen {
			el = LongEdgeMinLen
		}
		if int(maxPI.Score) < el*MinChainScoreIdentityPercent/100 {
			nmpi = mpi
			fl = flankLen + uint32(len(ne.Ks)-(kmerlen-1))
			return
		}
		//fmt.Printf("[ExtendPathChain]eID:%d maxPI.Arr:%v\n", ne.ID, maxPI.Arr[:maxPI.GetMapLen()Arr])
		maxPI.Base = uint32(kbIdx)
		maxPI.EID = uint32(ne.ID)
		maxPI.Kb = kb
		for i, idx := range maxPI.Arr {
			maxPI.Arr[i] = uint8(indexArr[idx])
		}
		extendMPI = maxPI
		if direction == BACKWARD {
			RevUint8Arr(extendMPI.Arr, int(extendMPI.LenArr))
		}
		//ekb = append(ekb, lastMki)
		//sort.Sort(MapingKmerInfoArr(ekb))
		//fmt.Printf("[ExtendPathChain] ekb: %v\n", ekb)
		//ekbIdx := LocalMapingKmerInfoArr(ekb, lastMki)
		/*var ekbIdx int
		if direction == FORWARD {
			ekbIdx = 0
		} else {
			ekbIdx = len(ekb) - 1
		}*/
		//fmt.Printf("[ExtendPathChain] ekb: \n")
		//PrintAddedMpArr(ekb)
	}

	//fmt.Printf("[ExtendPathChain]extendMpArr: \n")
	//PrintAddedMpArr(extendMpArr)

	//fmt.Printf("[ExtendPathChain]extendMPI:%v\n", extendMPI.Arr[:extendMPI.GetMapLen()Arr])
	//fmt.Printf("[ExtendPathChain] mpArr: %v\n", mpArr)
	//fmt.Printf("[ExtendPathChain] exendMpArr: %v\n", extendMpArr)

	if len(extendMPI.Arr) > 0 {
		distanceE := GetMKI2EdgeDistance(extendEPI, ne.ID, lastMki, direction, edgesArr, int(kmerlen))
		//fmt.Printf("[ExtendPathChain] distanceE: %v\n", distanceE)
		newScore, fl = GetMappingDBGExtendInfo(mpi, extendMPI, score, flankLen, strand, direction, len(ne.Ks), distanceE, int(kmerlen), edgesArr)
		if newScore > score {
			nmpi = extendMPI
		} else {
			nmpi = mpi
		}
	} else {
		nmpi = mpi
		//copy(nmpi.Arr[:mpi.GetMapLen()Arr], mpi.Arr[:mpi.GetMapLen()Arr])
		fl = flankLen + uint32(len(ne.Ks)-(kmerlen-1))
	}

	/*fmt.Printf("[ExtendPathChain] lastMKI: %v extendMpArr: \n", PrintMKI(lastMki))
	extL := len(newArr) - len(extendMpArr) - 3
	if extL < 0 {
		extL = 0
	}
	PrintAddedMpArr(newArr[extL:])*/
	//fmt.Printf("[ExtendPathChain] flankLen: %v, added mki num: %d, score: %d\n", fl, len(extendMpArr)-1, newScore)
	//fmt.Printf("[ExtendPathChain] newArr: \n")
	//PrintAddedMpArr(newArr)
	//fmt.Printf("[ExtendPathChain] ekb: \n")
	//PrintAddedMpArr(ekb)
	//fmt.Printf("[ExtendPathChain] exendMpArr: %v\n", extendMpArr)
	//fmt.Printf("[ExtendPathChain] newArr: %v\n", newArr)

	/*//startRPos := int(lastMki.GetPosition()) + int(lastMki.GetMapLen())
		for i := idx + 1; i < len(kb); i++ {
			mki := kb[i]
			lRPos := int(lastMki.GetPosition())+int(lastMki.GetMapLen())
			if mki.EID != uint32(ne.ID) || strand != mki.GetStrand() {
				continue
			}
			if int(mki.GetPosition()) < lRPos {
				continue
			}
			if int(mki.GetPosition())-lRPos > 500 {
				break
			}
			dR := int(mki.GetPosition()) - lRPos
			var dE int
			if strand == PLUS {
				dE = int(mki.GetEPosition()) - (int(lastMki.GetEPosition()) + int(lastMki.GetMapLen()))
			} else {
				dE = int(lastMki.GetEPosition()) - (int(mki.GetEPosition()) + int(mki.GetMapLen()))
			}
			if dE <= 0 {
				continue
			}
			min := Min(dR, dE)
			gapLen := AbsInt(dR-dE)
			if gapLen > min/5+10 {
				continue
			}
			a := Min(min, int(mki.GetMapLen()))
			var b float64
			if gapLen > 1 {
				b = float64(0.005)*float64(mki.GetMapLen())*float64(gapLen) + float64(0.5)*math.Log2(float64(gapLen))
			} else {
				b = 1
			}
			sc = Min(a-int(b), int(nc.GetMapLen()))
				sc := int(mki.GetMapLen()) - AbsInt(dR-dE)
				if sc > 0 {
					extendMpArr = append(extendMpArr, mki)
					lastMki = mki
					newScore += sc
				}
			}
		}
	} else { // direction == BACKWARD
		lastMki := extendMpArr[len(extendMpArr)-1]
		idx := LocalMapingKmerInfoArr(kb, lastMki)
		startRPos := int(lastMki.GetPosition())
		for i := idx - 1; i >= 0; i-- {
			mki := kb[i]
			if int(mki.GetPosition()) < startRPos-(flankLen+len(ne.Ks)-(kmerlen-1)+len(ne.Ks)/10) {
				break
			}
			if mki.EID != uint32(ne.ID) || strand != mki.GetStrand() {
				continue
			}
			dR := int(lastMki.GetPosition()) - (int(mki.GetPosition()) + int(mki.GetMapLen()))
			var dE int
			if strand == PLUS {
				dE = int(lastMki.GetEPosition()) - (int(mki.GetEPosition()) + int(mki.GetMapLen()))
			} else {
				dE = int(mki.GetEPosition()) - (int(lastMki.GetEPosition()) + int(lastMki.GetMapLen()))
			}
			if dE <= 0 {
				continue
			}
			min := Min(dR, dE)
			if AbsInt(dR-dE) < min/10+8 {
				sc := int(mki.GetMapLen()) - AbsInt(dR-dE)
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

func SharePathLen(path, path1 []uint32, edgesArr []DBGEdge, kmerlen, MaxPathLen, MaxBubbleSeqLen int) (sl int) {
	i, j := 0, 0
	//min := MinInt(len(path), len(path1))
	for i < len(path) && j < len(path1) {
		if path[i] != path1[i] {
			if i > 0 && j > 0 {
				var pathArr [][]uint32
				pathArr = append(pathArr, path[i-1:], path1[j-1:])
				bubPathArr, num := GetBubblePathArr(pathArr, edgesArr, kmerlen, MaxPathLen, MaxBubbleSeqLen)
				if num == 2 {
					if reflect.DeepEqual(path[i:i+len(bubPathArr[0].Path)], bubPathArr[0].Path) {
						i += len(bubPathArr[0].Path)
						j += len(bubPathArr[1].Path)
					} else {
						i += len(bubPathArr[1].Path)
						j += len(bubPathArr[0].Path)
					}
					continue
				}
			}
			break
		}
		i++
		j++
	}
	sl = i
	return
}

func GetExtendPathInfoPath(extendEPI ExtendEdgePathInfo) (path []uint32) {
	path = make([]uint32, 0, 10)
	epi := extendEPI
	path = append(path, epi.EID)
	for epi.Last != nil {
		epi = *epi.Last
		path = append(path, epi.EID)
	}
	ReverseUint32Arr(path)
	return path
}

func GetBestAlignmentPath(bestArr []ExtendPathInfo, direction uint8, edgesArr []DBGEdge, nodesArr []DBGNode, read ReadInfo, SeedLen, kmerlen, minMapLen, MaxPathLen, MaxBubbleSeqLen, WinSize int, maxPathInfoArr []MaxPathInfo, anchor MapingKmerInfo) (bestEPI ExtendPathInfo) {
	var readRegionSeq []byte
	{
		readRegionLen := minMapLen * 6 / 5
		//mpi := bestArr[0].MPI
		//mkia := anchor
		if direction == FORWARD {
			start := int(anchor.GetPosition())
			if start+readRegionLen > len(read.Seq) {
				readRegionLen = len(read.Seq) - start
			}
			readRegionSeq = read.Seq[start : start+readRegionLen]
		} else {
			start := int(anchor.GetPosition()) + int(anchor.GetMapLen())
			if readRegionLen > start {
				readRegionLen = start
			}
			readRegionSeq = read.Seq[start-readRegionLen : start]
		}
	}

	edgesSeqArr := make([][]byte, len(bestArr))
	for i, epi := range bestArr {
		if DebugModel {
			mpi := &epi.MPI
			mkb := mpi.Kb[mpi.Base+uint32(mpi.Arr[len(mpi.Arr)-1])]
			fmt.Printf("[GetBestAlignmentPath] bestArr[%v]: score: %v, mapLen: %v,RPos: %v, path: %v\n", i, epi.Score, epi.ExtendLen-epi.FlankLen, mkb.GetPosition(), epi.Path)
		}
		//epi.Path = GetExtendPathInfoPath(epi.ExtendEPI)
		//mkia := mpi.Kb[mpi.Base+uint32(mpi.Arr[0])]
		edgesSeq := GetDirectionEdgesSeqRegion(edgesArr, nodesArr, direction, uint32(kmerlen), minMapLen, epi.Path, anchor)
		edgesSeqArr[i] = edgesSeq
	}

	// alignment
	cgArr := AlignmentBlocks(readRegionSeq, edgesSeqArr, SeedLen, WinSize)
	maxScore := math.MinInt32
	idx := -1
	maxCount := 0
	for i, cg := range cgArr {
		score := int(cg.Mch) - int(cg.Mis) - int(cg.Ins) - int(cg.Del)
		if score > maxScore {
			maxScore = score
			idx = i
			maxCount = 1
		} else if score == maxScore {
			maxCount++
		}
	}

	if maxCount == 1 && maxScore > minMapLen*60/100 {
		bestEPI = bestArr[idx]
	} else if maxCount > 1 {
		var shareEPI ExtendPathInfo
		shareEPI = bestArr[idx]
		for i, cg := range cgArr {
			score := int(cg.Mch) - int(cg.Mis) - int(cg.Ins) - int(cg.Del)
			if i == idx || score != maxScore {
				continue
			}
			sl := SharePathLen(shareEPI.Path, bestArr[i].Path, edgesArr, kmerlen, MaxPathLen, MaxBubbleSeqLen)
			shareEPI.Path = shareEPI.Path[:sl]
			if sl == 0 {
				break
			}
		}
		bestEPI = shareEPI
	} else {
		var t ExtendPathInfo
		bestEPI = t
	}
	fmt.Printf("[GetBestAlignmentPath] len(bestArr):%d best score:%d \n", len(bestArr), bestEPI.Score)
	return
}

func IsInExtendpathInfoArr(epiArr []ExtendPathInfo, epi ExtendPathInfo) (cok bool) {
	for _, t := range epiArr {
		if len(t.Path) == len(epi.Path) && reflect.DeepEqual(t.Path, epi.Path) {
			cok = true
			break
		}
	}
	return
}

func GetBestEPathInfo(highQualEPathInfoArr []ExtendPathInfo, edgesArr []DBGEdge, nodesArr []DBGNode, read ReadInfo, SeedLen, kmerlen int, direction uint8, MaxPathLen, MaxBubbleSeqLen, WinSize int, maxPathInfoArr []MaxPathInfo, maxScoreArr []ScorePos, anchor MapingKmerInfo, ScoreWinSize int) (bestEPI ExtendPathInfo) {
	if len(highQualEPathInfoArr) == 0 {
		return
	}
	if len(highQualEPathInfoArr) == 1 {
		bestEPI = highQualEPathInfoArr[0]
		return
	}
	var maxLen uint32
	for _, epi := range highQualEPathInfoArr {
		l := epi.ExtendLen - epi.FlankLen
		if maxLen < l {
			maxLen = l
		}
	}
	//var bestArr []ExtendPathInfo
	bestScoreNum := 0
	//fmt.Printf("[GetBestEPathInfo]maxScoreArr: %v\n", maxScoreArr)
	if int(maxLen) < kmerlen*2 {
		return
	}

	var startPos int
	if direction == FORWARD {
		startPos = int(anchor.GetPosition())
	} else {
		startPos = int(anchor.GetPosition()) + int(anchor.GetMapLen())
	}

	for i, epi := range highQualEPathInfoArr {
		mpi := epi.MPI
		if epi.ExtendLen-epi.FlankLen < maxLen*8/10 || len(mpi.Arr) == 0 {
			highQualEPathInfoArr[i].Score = 0
			continue
		}
		mki := mpi.Kb[mpi.Base+uint32(mpi.Arr[len(mpi.Arr)-1])]
		var sp ScorePos
		if direction == FORWARD {
			sp.Pos = uint32(mki.GetPosition()) + mki.GetMapLen()
		} else {
			sp.Pos = uint32(mki.GetPosition())
		}
		sp.Score = uint32(epi.Score)
		if BiggerThanMaxScoreArr(maxScoreArr, startPos, sp, 0, direction, ScoreWinSize) {
			epi.Path = GetExtendPathInfoPath(epi.ExtendEPI)
			epi = DeleteUnmapPath(epi, direction, edgesArr, kmerlen)
			bestEPI = epi
			bestScoreNum++
			highQualEPathInfoArr[i] = epi
			if DebugModel {
				mpi := epi.MPI
				l := epi.ExtendLen - epi.FlankLen
				//mkia := mpi.Kb[mpi.Base+int(mpi.Arr[0])]
				fmt.Printf("[GetBestEPathInfo]high[%v]len(mpi.Arr): %v, Score: %v,Len: %v, Percent: %v, ExtendLen: %v, FlankLen: %v, Strand: %v, Coming: %v, Path: %v\n", i, len(mpi.Arr), epi.Score, l, epi.Score*100/l, epi.ExtendLen, epi.FlankLen, epi.GetStrandFlag(), epi.GetComingFlag(), epi.Path)
			}
		} else {
			highQualEPathInfoArr[i].Score = 0
		}
	}

	if bestScoreNum > 1 {
		/*if minMaplen < kmerlen {
			return
		}*/
		//fmt.Printf("[GetBestEPathInfo] minMapLen: %v\n", minMaplen)
		maxArr := make([]ExtendPathInfo, 0, 10)
		minMaplen := uint32(math.MaxInt32)
		for j, epi := range highQualEPathInfoArr {
			//fmt.Printf("[GetBestEPathInfo] bestArr[%v]: score: %v, mapLen: %v,RPos: %v, path: %v\n", j, epi.Score, epi.ExtendLen-epi.FlankLen, epi.MpArr[len(epi.MpArr)-1].GetPosition(), epi.Path)
			if epi.Score == 0 {
				continue
			}
			ok := true
			epiMapLen := epi.ExtendLen - epi.FlankLen
			for x := j + 1; x < len(highQualEPathInfoArr); x++ {
				ele := highQualEPathInfoArr[x]
				if epi.Score == 0 {
					continue
				}
				eleMapLen := ele.ExtendLen - ele.FlankLen
				if epiMapLen <= eleMapLen && len(epi.Path) <= len(ele.Path) {
					if reflect.DeepEqual(epi.Path, ele.Path[:len(epi.Path)]) {
						ok = false
						break
					}
				} else if eleMapLen < epiMapLen && len(epi.Path) >= len(ele.Path) {
					if reflect.DeepEqual(epi.Path[:len(ele.Path)], ele.Path) {
						highQualEPathInfoArr[x].Score = 0
					}
				} else {
					if epiMapLen <= eleMapLen {
						if epi.Score > ele.Score {
							highQualEPathInfoArr[x].Score = 0
						}
					} else { // epiMapLen > eleMapLen
						if epi.Score <= ele.Score {
							ok = false
							break
						}
					}
				}
			}
			if ok {
				maxArr = append(maxArr, epi)
				if epi.ExtendLen-epi.FlankLen < minMaplen {
					minMaplen = epi.ExtendLen - epi.FlankLen
				}
			}
		}
		if len(maxArr) == 1 {
			bestEPI = maxArr[0]
		} else if minMaplen > uint32(kmerlen)*2 && len(maxArr) > 1 && len(maxArr) < 20 {
			bestEPI = GetBestAlignmentPath(maxArr, direction, edgesArr, nodesArr, read, SeedLen, kmerlen, int(minMaplen), MaxPathLen, MaxBubbleSeqLen, WinSize, maxPathInfoArr, anchor)
		} else {
			var tmp ExtendPathInfo
			bestEPI = tmp
		}
	}
	return
}

func MuskBubble(eArr []uint32, stk []ExtendPathInfo) ([]ExtendPathInfo, int) {
	var deleteNum int
	if len(stk) < 2 {
		return stk, deleteNum
	}
	s1, s2 := stk[len(stk)-2], stk[len(stk)-1]
	if s1.ExtendEPI.EID == eArr[0] && s2.ExtendEPI.EID == eArr[1] {
		if s1.Score < s2.Score {
			stk[len(stk)-2] = s2
			stk = stk[:len(stk)-1]
			deleteNum++
		} else {
			stk = stk[:len(stk)-1]
			deleteNum++
		}
	}
	return stk, deleteNum
}

func MuskDiffLenBubble(le *DBGEdge, eArr []uint32, stk []ExtendPathInfo, edgesArr []DBGEdge, nodesArr []DBGNode, kmerlen, MaxPathLen, MaxBubbleSeqLen int) ([]ExtendPathInfo, int) {
	var deleteNum int
	if len(stk) < len(eArr) || len(eArr) < 2 {
		return stk, deleteNum
	}
	if le.GetTwoEdgesCycleFlag() > 0 || le.StartNID == le.EndNID {
		return stk, deleteNum
	}

	maxIdx := -1
	max, sMax := math.MinInt32, math.MinInt32
	for i, id := range eArr {
		e := &edgesArr[id]
		if e.GetSeqLen() > max {
			max = e.GetSeqLen()
			maxIdx = i
		} else if e.GetSeqLen() > sMax {
			sMax = e.GetSeqLen()
		}
	}
	e := &edgesArr[eArr[maxIdx]]
	if 2*(kmerlen-1) <= max && max < MaxBubbleSeqLen && sMax < 2*(kmerlen-1) && e.GetTwoEdgesCycleFlag() == 0 && e.StartNID != e.EndNID {

	} else {
		return stk, deleteNum
	}
	var nd DBGNode
	if le.StartNID == e.StartNID || le.StartNID == e.EndNID {
		nd = nodesArr[le.StartNID]
	} else {
		nd = nodesArr[le.EndNID]
	}
	pathArr := GetNextPathArr(le, nd, MaxPathLen+1, max*11/10, edgesArr, nodesArr, kmerlen)
	_, num := GetBubblePathArr(pathArr, edgesArr, kmerlen, MaxPathLen, MaxBubbleSeqLen)
	if num != len(pathArr) {
		return stk, deleteNum
	}

	for i := 0; i < len(eArr); i++ {
		epi := stk[len(stk)-len(eArr)+i]
		if i == maxIdx {
			if epi.ExtendEPI.EID != eArr[i] {
				log.Fatalf("[MuskDiffLenBubble] eArr[%v]: %v not concensis with stk\n", i, eArr[i])
			}
			continue
		} else {
			stk[len(stk)-len(eArr)+i].Score = 0
			deleteNum++
		}
	}
	stk[len(stk)-len(eArr)] = stk[len(stk)-len(eArr)+maxIdx]
	stk = stk[:len(stk)-len(eArr)+1]
	return stk, deleteNum
}

func DeleteUnmapPath(bestEPI ExtendPathInfo, direction uint8, edgesArr []DBGEdge, kmerlen int) ExtendPathInfo {
	fl := int(bestEPI.FlankLen)
	el := int(bestEPI.ExtendLen)
	if len(bestEPI.MPI.Arr) == 0 {
		return bestEPI
	}
	mpi := bestEPI.MPI
	//lastMki := mpi.Kb[mpi.Base+int(mpi.Arr[mpi.GetMapLen()Arr-1])]
	eID := mpi.EID
	i := len(bestEPI.Path) - 1
	for ; i > 0; i-- {
		e := edgesArr[bestEPI.Path[i]]
		if len(e.Ks)-(kmerlen-1) < fl {
			fl -= (len(e.Ks) - (kmerlen - 1))
			el -= (len(e.Ks) - (kmerlen - 1))
			continue
		}

		if uint32(e.ID) == eID {
			if len(e.Ks)-fl < kmerlen*2 || e.GetBubbleFlag() > 0 {
				i--
			}
			break
		}
	}
	bestEPI.Path = bestEPI.Path[:i+1]
	if len(bestEPI.Path1) > len(bestEPI.Path) {
		bestEPI.Path1 = bestEPI.Path1[:len(bestEPI.Path)]
	}
	bestEPI.FlankLen = uint32(fl)
	bestEPI.ExtendLen = uint32(el)

	return bestEPI
}

func MergeTwoFlankMapingPath(bestEPIB, bestEPIF ExtendPathInfo, mpi MaxPathInfo, eID uint32) (longRMInfo LongReadMappingInfo) {
	path := make([]uint32, 0, len(bestEPIB.Path)+len(bestEPIF.Path)+1)
	path1 := make([]uint32, 0, len(bestEPIB.Path)+len(bestEPIF.Path)+1)
	cb1 := mpi.Kb[mpi.Base+uint32(mpi.Arr[0])]
	cb2 := mpi.Kb[mpi.Base+uint32(mpi.Arr[len(mpi.Arr)-1])]
	if len(bestEPIB.Path) <= 1 && len(bestEPIF.Path) <= 1 {
		path = append(path, uint32(eID))
		path1 = append(path1, uint32(0))
		longRMInfo.RStart, longRMInfo.REnd = uint32(cb1.GetPosition()), uint32(cb2.GetPosition())+cb2.GetMapLen()
		if cb1.GetStrand() == PLUS {
			longRMInfo.EStart, longRMInfo.EEnd = cb1.GetEPosition(), cb2.GetEPosition()+uint32(cb2.GetMapLen())
		} else {
			longRMInfo.EStart, longRMInfo.EEnd = cb1.GetEPosition()+uint32(cb1.GetMapLen()), cb2.GetEPosition()
		}
	} else {
		if len(bestEPIB.Path) > 1 {
			mpi := bestEPIB.MPI
			lastB := mpi.Kb[mpi.Base+uint32(mpi.Arr[len(mpi.Arr)-1])]
			//lastB := bestEPIB.MpArr[len(bestEPIB.MpArr)-1]
			longRMInfo.RStart = uint32(lastB.GetPosition())
			if lastB.GetStrand() == PLUS {
				longRMInfo.EStart = lastB.GetEPosition()
			} else {
				longRMInfo.EStart = lastB.GetEPosition() + uint32(lastB.GetMapLen())
			}
			path = append(path, GetReverseUint32Arr(bestEPIB.Path)...)
			if len(bestEPIB.Path1) > 0 {
				path1 = append(path1, GetReverseUint32Arr(bestEPIB.Path1)...)
			} else {
				p := make([]uint32, len(bestEPIB.Path))
				path1 = append(path1, p...)
			}
		} else {
			path = append(path, uint32(eID))
			path1 = append(path1, uint32(0))
			longRMInfo.RStart = uint32(cb1.GetPosition())
			if cb1.GetStrand() == PLUS {
				longRMInfo.EStart = cb1.GetEPosition()
			} else {
				longRMInfo.EStart = cb1.GetEPosition() + uint32(cb1.GetMapLen())
			}
		}

		if len(bestEPIF.Path) > 1 {
			mpi := bestEPIF.MPI
			lastF := mpi.Kb[mpi.Base+uint32(mpi.Arr[len(mpi.Arr)-1])]
			//lastF := bestEPIF.MpArr[len(bestEPIF.MpArr)-1]
			longRMInfo.REnd = uint32(lastF.GetPosition()) + lastF.GetMapLen()
			if lastF.GetStrand() == PLUS {
				longRMInfo.EEnd = lastF.GetEPosition() + uint32(lastF.GetMapLen())
			} else {
				longRMInfo.EEnd = lastF.GetEPosition()
			}
			path = append(path, bestEPIF.Path[1:]...)
			if len(bestEPIF.Path1) > 0 {
				path1 = append(path1, bestEPIF.Path1[1:]...)
			} else {
				p := make([]uint32, len(bestEPIF.Path))
				path1 = append(path1, p[1:]...)
			}
		} else {
			longRMInfo.REnd = uint32(cb2.GetPosition()) + cb2.GetMapLen()
			if cb1.GetStrand() == PLUS {
				longRMInfo.EEnd = cb2.GetEPosition() + uint32(cb2.GetMapLen())
			} else {
				longRMInfo.EEnd = cb2.GetEPosition()
			}
		}
	}
	longRMInfo.Path[0] = path
	longRMInfo.Path[1] = path1
	return
}

//"denote max score before Pos"
type ScorePos struct {
	Pos   uint32 // ReadPos
	Score uint32
}

func FindLastSP(maxScoreArr []ScorePos, pos int) (sp ScorePos) {
	for i := pos; i >= 0; i-- {
		sp = maxScoreArr[i]
		if sp.Score > 0 {
			break
		}
	}
	return
}

func BiggerThanMaxScoreArr(maxScoreArr []ScorePos, startPos int, sa ScorePos, extendLen int, direciton uint8, ScoreWinSize int) bool {
	ok := true
	AllowDiffLen := int(ScoreWinSize) + 10
	if direciton == FORWARD {
		idx := (int(sa.Pos) - startPos) / ScoreWinSize
		if int(idx) >= len(maxScoreArr) {
			idx = len(maxScoreArr) - 1
		}
		if idx < 0 {
			if len(maxScoreArr) == 0 {
				return ok
			}
			idx = 0
		}
		sp := maxScoreArr[idx]
		if sp.Score == 0 {
			sp = FindLastSP(maxScoreArr, int(idx)-1)
		}
		if sp.Pos < sa.Pos {
			if sp.Score >= sa.Score {
				ok = false
			}
		} else if sp.Pos == sa.Pos {
			if sp.Score > sa.Score {
				ok = false
			}
		} else { // sp.Pos > sa.Pos
			if sp.Score > sa.Score && idx > 0 {
				sp = maxScoreArr[idx-1]
				if sp.Score == 0 {
					sp = FindLastSP(maxScoreArr, int(idx)-2)
				}
				if sp.Score >= sa.Score {
					ok = false
				}
			}
		}

		// test extendLen
		if ok && extendLen > 0 {
			idx := (extendLen - AllowDiffLen - ScoreWinSize + 1) / ScoreWinSize
			if idx >= len(maxScoreArr) {
				idx = len(maxScoreArr) - 1
			}
			if idx < 0 {
				if len(maxScoreArr) == 0 {
					return ok
				}
				idx = 0
			}
			sp := maxScoreArr[idx]
			if sp.Score == 0 {
				sp = FindLastSP(maxScoreArr, int(idx)-1)
			}
			if sp.Score > sa.Score {
				ok = false
			}
		}
	} else {
		idx := (startPos - int(sa.Pos)) / ScoreWinSize
		if int(idx) >= len(maxScoreArr) {
			idx = len(maxScoreArr) - 1
		}
		if idx < 0 {
			if len(maxScoreArr) == 0 {
				return ok
			}
			idx = 0
		}
		sp := maxScoreArr[idx]
		if sp.Score == 0 {
			sp = FindLastSP(maxScoreArr, int(idx)-1)
		}
		if sp.Pos > sa.Pos {
			if sp.Score >= sa.Score {
				ok = false
			}
		} else if sp.Pos == sa.Pos {
			if sp.Score > sa.Score {
				ok = false
			}
		} else { // sp.Pos > sa.Pos
			if sp.Score > sa.Score && idx > 0 {
				sp = maxScoreArr[idx-1]
				if sp.Score == 0 {
					sp = FindLastSP(maxScoreArr, int(idx)-2)
				}
				if sp.Score >= sa.Score {
					ok = false
				}
			}
		}

		// test extendLen
		if ok && extendLen > 0 {
			idx := (extendLen - AllowDiffLen - ScoreWinSize + 1) / ScoreWinSize
			if idx >= len(maxScoreArr) {
				idx = len(maxScoreArr) - 1
			}
			if idx < 0 {
				if len(maxScoreArr) == 0 {
					return ok
				}
				idx = 0
			}
			sp := maxScoreArr[idx]
			if sp.Score == 0 {
				sp = FindLastSP(maxScoreArr, idx-1)
			}
			if sp.Score > sa.Score {
				ok = false
			}
		}
	}
	return ok
}

func AddedToMaxScoreArr(maxScoreArr []ScorePos, startPos int, sa ScorePos, direction uint8, ScoreWinSize int) ([]ScorePos, bool) {
	added := false
	if direction == FORWARD {
		regLen := int(sa.Pos) - startPos
		idx := regLen / ScoreWinSize
		if idx >= len(maxScoreArr) {
			if cap(maxScoreArr) > idx {
				maxScoreArr = maxScoreArr[:idx+1]
				maxScoreArr[idx] = sa
			} else {
				na := make([]ScorePos, idx*6/5+10)
				copy(na, maxScoreArr)
				na[idx] = sa
				maxScoreArr = na[:idx+1]
			}
			added = true
		} else {
			sp := maxScoreArr[idx]
			if sp.Score == 0 {
				maxScoreArr[idx] = sa
				added = true
			} else {
				if sp.Pos >= sa.Pos && sp.Score <= sa.Score {
					maxScoreArr[idx] = sa
					added = true
				} else {
					if float32(sa.Score)/float32(regLen) > float32(sp.Score)/(float32(sp.Pos)-float32(startPos)) {
						maxScoreArr[idx] = sa
						added = true
					}
				}
			}
		}
	} else {
		regLen := startPos - int(sa.Pos)
		idx := regLen / ScoreWinSize
		if idx >= len(maxScoreArr) {
			if cap(maxScoreArr) > idx {
				maxScoreArr = maxScoreArr[:idx+1]
				maxScoreArr[idx] = sa
			} else {
				na := make([]ScorePos, idx*6/5+10)
				copy(na, maxScoreArr)
				na[idx] = sa
				maxScoreArr = na[:idx+1]
			}
			added = true
		} else {
			sp := maxScoreArr[idx]
			if sp.Score == 0 {
				maxScoreArr[idx] = sa
				added = true
			} else {
				if sp.Pos <= sa.Pos && sp.Score <= sa.Score {
					maxScoreArr[idx] = sa
					added = true
				} else {
					if float32(sa.Score)/float32(regLen) > float32(sp.Score)/(float32(startPos)-float32(sp.Pos)) {
						maxScoreArr[idx] = sa
						added = true
					}
				}
			}
		}
	}
	//fmt.Printf("[AddedToMaxScoreArr]sa:%v added:%v\n", sa, added)
	return maxScoreArr, added
}

func AddedToMaxScoreArrMPI(maxScoreArr []ScorePos, startPos int, mpi MaxPathInfo, score int, direction uint8, ScoreWinSize int) ([]ScorePos, int) {
	addNum := 0
	sc := float32(score)
	sz := len(mpi.Arr)
	if sz >= ArrSize {
		sz = ArrSize / 2
		sc = float32(score) - mpi.Score
	}
	var ok bool
	for i := sz - 1; i >= 0; i-- {
		var sa ScorePos
		mki := mpi.Kb[mpi.Base+uint32(mpi.Arr[i])]
		if direction == FORWARD {
			sa.Pos = uint32(mki.GetPosition()) + mki.GetMapLen()
			if int(sa.Pos) <= startPos {
				continue
			}
		} else {
			sa.Pos = uint32(mki.GetPosition())
			if int(sa.Pos) >= startPos {
				continue
			}
		}
		sa.Score = uint32(sc)
		maxScoreArr, ok = AddedToMaxScoreArr(maxScoreArr, startPos, sa, direction, ScoreWinSize)
		if ok {
			addNum++
		}
		if i > 0 {
			sc -= GetTwoBlockScore2(mpi.Kb[mpi.Base+uint32(mpi.Arr[i-1])], mki, 1)
		}
	}
	//fmt.Printf("[AddedToMaxScoreArrMPI]addNum:%d\n", addNum)
	return maxScoreArr, addNum
}

func GetUniqueMaxPathInfo(maxPathInfoArr []MaxPathInfo, kmerlen, SeedLen, ReadLen int, edgesArr []DBGEdge, MinChainScoreIdentityPercent int) (uniqueMaxPathInfo MaxPathInfo) {
	// filter low mapping quality edge and filter only maping middle region of edge
	filterFlag := make([]bool, len(maxPathInfoArr))
	boundLen := kmerlen / 2
	for i, mpi := range maxPathInfoArr {
		if filterFlag[i] {
			continue
		}
		mka, mkb := mpi.Kb[mpi.Base+uint32(mpi.Arr[0])], mpi.Kb[mpi.Base+uint32(mpi.Arr[len(mpi.Arr)-1])]
		Start, End := int(mka.GetPosition()), int(mkb.GetPosition())+int(mkb.GetMapLen())
		strand := mka.GetStrand() == mka.GetEStrand()
		e := edgesArr[mpi.EID]
		if strand == PLUS {
			eS, eE := int(mka.GetEPosition()), int(mkb.GetEPosition())+int(mkb.GetMapLen())
			if Start > boundLen && eS > boundLen {
				filterFlag[i] = true
				continue
			}
			if End < ReadLen-boundLen && eE < len(e.Ks)-boundLen {
				filterFlag[i] = true
				continue
			}
		} else {
			eS, eE := int(mkb.GetEPosition()), int(mka.GetEPosition())+int(mka.GetMapLen())
			if Start > boundLen && len(e.Ks)-eE > boundLen {
				filterFlag[i] = true
				continue
			}
			if ReadLen-End > boundLen && eS > boundLen {
				filterFlag[i] = true
				continue
			}
		}

		for j := i + 1; j < len(maxPathInfoArr); j++ {
			if filterFlag[j] {
				continue
			}
			hmpi := maxPathInfoArr[j]
			mka, mkb := hmpi.Kb[hmpi.Base+uint32(hmpi.Arr[0])], hmpi.Kb[hmpi.Base+uint32(hmpi.Arr[len(hmpi.Arr)-1])]
			hStart, hEnd := int(mka.GetPosition()), int(mkb.GetPosition())+int(mkb.GetMapLen())
			if End <= hStart {
				break
			}
			if hStart <= Start && End <= hEnd && int(mpi.Score)*100/(End-Start) < int(hmpi.Score)*100/(hEnd-hStart) {
				filterFlag[i] = true
				break
			} else if hStart <= Start && End <= hEnd && mpi.Score > hmpi.Score {
				filterFlag[j] = true
			} else if hStart == Start && End == hEnd && mpi.Score < hmpi.Score {
				filterFlag[i] = true
				break
			} else if hStart <= Start && End <= hEnd && End-Start <= hEnd-hStart && int(mpi.Score)*100/(End-Start) > (int(hmpi.Score)*100/(hEnd-hStart)) {
				filterFlag[j] = true
			}
		}
	}

	if DebugModel {
		for x, mpi := range maxPathInfoArr {
			if filterFlag[x] {
				continue
			}
			mka, mkb := mpi.Kb[mpi.Base+uint32(mpi.Arr[0])], mpi.Kb[mpi.Base+uint32(mpi.Arr[len(mpi.Arr)-1])]
			start, end := uint32(mka.GetPosition()), uint32(mkb.GetPosition())+mkb.GetMapLen()
			if start == end {
				end++
			}
			fmt.Printf("[GetUMPI]MPIArr[%v]: EdgeID:%v ReadRegion[%v--%v] eLen:%d Len:%v score:%v Percent:%v\n", x, mka.EID, start, end, edgesArr[mka.EID].GetSeqLen(), end-start, mpi.Score, int(mpi.Score)*100/(int(end)-int(start)))
		}
	}

	// prior choose unique edge or bubble and bubbleRepeat edge
	for i, mpi := range maxPathInfoArr {
		if filterFlag[i] {
			continue
		}
		mka := mpi.Kb[mpi.Base+uint32(mpi.Arr[0])]
		//mka, mkb := mpi.Kb[mpi.Base+int(mpi.Arr[0])], mpi.Kb[mpi.Base+int(mpi.Arr[mpi.GetMapLen()Arr-1])]
		eID := mka.EID
		e := edgesArr[eID]
		if int(mpi.Score) < e.GetSeqLen()*MinChainScoreIdentityPercent/100 {
			filterFlag[i] = true
			continue
		}
		if e.GetUniqueFlag() == 0 && e.GetBubbleFlag() == 0 && e.GetBubbleRepeatFlag() == 0 {
			continue
		}
		if mpi.Score > uniqueMaxPathInfo.Score {
			//fmt.Printf("[GetUniqueMaxPathInfo]MPIArr[%v]: EdgeID:%v  score:%v\n", i, eID, mpi.Score)
			uniqueMaxPathInfo = mpi
		}
	}

	if uniqueMaxPathInfo.Score > 0 {
		return
	}

	// second choose largest score
	for i, mpi := range maxPathInfoArr {
		if filterFlag[i] {
			continue
		}
		if mpi.Score > uniqueMaxPathInfo.Score {
			uniqueMaxPathInfo = mpi
		}
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

func IsSameDBGNode(e, e2 DBGEdge, strand, strand2 bool) bool {
	if strand == strand2 {
		return e.StartNID == e2.StartNID || e.EndNID == e2.EndNID
	} else {
		return e.StartNID == e2.EndNID || e.EndNID == e2.StartNID
	}
}

func DiffFloat32(f1, f2 float32) (diff float32) {
	if f1 > f2 {
		diff = f1 - f2
	} else {
		diff = f2 - f1
	}
	return diff
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

	/*// filter highQMPIArr[i], highQMPIArr[i+1] share region > kmerlen
	filterFlag = filterFlag[:len(highQMPIArr)]
	for i := range filterFlag {
		filterFlag[i] = false
	}

	for i, mpi := range highQMPIArr {
		if filterFlag[i] {
			continue
		}
		mka, mkb := mpi.Kb[mpi.Base+uint32(mpi.Arr[0])], mpi.Kb[mpi.Base+uint32(mpi.Arr[len(mpi.Arr)-1])]
		start, end := uint32(mka.GetPosition()), uint32(mkb.GetPosition())+mkb.GetMapLen()
		for j := i + 1; j < len(highQMPIArr); j++ {
			if filterFlag[j] {
				continue
			}
			nm := highQMPIArr[j]
			a, b := nm.Kb[nm.Base+uint32(nm.Arr[0])], nm.Kb[nm.Base+uint32(nm.Arr[len(nm.Arr)-1])]
			nstart, nend := uint32(a.GetPosition()), uint32(b.GetPosition())+b.GetMapLen()
			if int(nstart) >= int(end)-kmerlen {
				break
			}
			if start == nstart || end >= nend || (int(end)-int(start)-(int(nend)-int(nstart)) == int(mpi.Score)-int(nm.Score)) {
				continue
			}

			sp, ep := nstart, end
			sz := len(mpi.Arr)
			dr := BACKWARD
			if sz == ArrSize {
				sz = sz / 2
			}
			sc, _ := GetMPIRegionScore(mpi, sz, dr, sp, ep)

			sz = len(nm.Arr)
			dr = FORWARD
			if sz == ArrSize {
				sz = sz / 2
			}
			sc2, _ := GetMPIRegionScore(nm, sz, dr, sp, ep)
			if sc < sc2 {
				filterFlag[i] = true
				break
			} else if sc > sc2 {
				filterFlag[j] = true
			}
		}
	}

	count = 0
	for i := range highQMPIArr {
		if filterFlag[i] {
			continue
		}
		highQMPIArr[count] = highQMPIArr[i]
		count++
	}
	fmt.Printf("[GetMaxScoreRegionMPI]len(highQMPIArr):%d after filter share Region count:%d\n", len(highQMPIArr), count)
	highQMPIArr = highQMPIArr[:count]

	if DebugModel {
		for x, mpi := range highQMPIArr {
			mka, mkb := mpi.Kb[mpi.Base+uint32(mpi.Arr[0])], mpi.Kb[mpi.Base+uint32(mpi.Arr[len(mpi.Arr)-1])]
			start, end := uint32(mka.GetPosition()), uint32(mkb.GetPosition())+mkb.GetMapLen()
			if start == end {
				end++
			}
			fmt.Printf("[GetUMPI]highQMPIArr[%v]: EdgeID:%v ReadRegion[%v--%v] eLen:%d Len:%v score:%v Percent:%v\n", x, mka.EID, start, end, edgesArr[mka.EID].GetSeqLen(), end-start, mpi.Score, int(mpi.Score)*100/(int(end)-int(start)))
		}
	}

	if len(highQMPIArr) > math.MaxUint16 {
		fmt.Fprintf(os.Stderr, "[Warning][GetMaxScoreRegionMPI]len(highQMPIArr):%d > %d\n", len(highQMPIArr), math.MaxUint16)
	} */
	return maxPathInfoArr
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
		strand := mka.GetStrand() == mka.GetEStrand()
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
			nstrand := mka.GetStrand() == mka.GetEStrand()
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

			/*else if nstart == start || nend == end {
				if nstart == start {
					if j+1 < len(highQMPIArr) {
						ta := highQMPIArr[j+1].Kb[highQMPIArr[j+1].Base+uint32(highQMPIArr[j+1].Arr[0])]
						if uint32(ta.GetPosition()) == start {
							continue
						}
					}
					strand, strand2 := mka.GetStrand(), a.GetStrand()
					if IsSameDBGNode(e, e2, strand, strand2) {
						priorIdxArr[i] = uint16(i) + 1
						priorIdxArr[j] = uint16(j) + 1
						continue
					}
				} else {
					if j+1 < len(highQMPIArr) {
						tn := highQMPIArr[j+1]
						tb := tn.Kb[tn.Base+uint32(tn.Arr[len(tn.Arr)-1])]
						if uint32(tb.GetPosition())+tb.GetMapLen() == end {
							continue
						}
					}
					strand, strand2 := mka.GetStrand(), a.GetStrand()
					if IsSameDBGNode(e, e2, strand, strand2) {
						priorIdxArr[i] = uint16(i) + 1
						priorIdxArr[j] = uint16(j) + 1
						continue
					}
				}
			} else if nstart < end-uint32(kmerlen) {
				if mpi.Score > nm.Score*2 {
					priorIdxArr[i] = uint16(i) + 1
				} else {
					priorIdxArr[i] = math.MaxUint16
				}
				priorIdxArr[j] = math.MaxUint16
			}*/
		}

		/*if priorIdxArr[i] == math.MaxUint16 {
			continue
		}
		if bubbleCount == 1 && bubbleIdx >= 0 && priorIdxArr[bubbleIdx] != math.MaxUint16 {
			e2 := edgesArr[highQMPIArr[bubbleIdx].EID]
			if IsBubble(e, e2, nodesArr) {
				priorIdxArr[i] = uint16(i) + 1
				priorIdxArr[bubbleIdx] = uint16(bubbleIdx) + 1
			} else {
				priorIdxArr[i] = math.MaxUint16
				priorIdxArr[bubbleIdx] = math.MaxUint16
			}
			continue
		}

		priorIdxArr[i] = uint16(i) + 1
		if e.GetUniqueFlag() == 0 && e.GetBubbleFlag() == 0 && e.GetBubbleRepeatFlag() == 0 {
			if regl > MinRegLen && e.StartNID != e.EndNID {
				if mpi.Score > score {
					if repeatOK {
						if mpi.Score > score {
							anchorIdx = i
							score = mpi.Score
						}
					} else {
						if mpi.Score > score*2 {
							anchorIdx = i
							score = mpi.Score
							repeatOK = true
						}
					}
				}
			}
		} else {
			if repeatOK {
				if mpi.Score > score/2 {
					anchorIdx = i
					score = mpi.Score
					repeatOK = false
				}
			} else {
				if mpi.Score > score {
					anchorIdx = i
					score = mpi.Score
				}
			}
		}*/
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

/*func PrintMPIArr(mpi MaxPathInfo) {
	for i, idx := range mpi.Arr {
		mki := mpi.Kb[mpi.Base+uint32(idx)]
		fmt.Printf("[PrintMPI]arr[%d] EID:%v strand:%v RPos:%d EPos:%d Len:%d\n", i, mki.EID, mki.GetStrand(), mki.GetPosition(), mki.GetEPosition(), mki.GetMapLen())
	}
}*/
func PrintMPI(mpi *MaxPathInfo, edgesArr []DBGEdge) string {
	var s string
	if mpi.LenArr == 0 {
		return s
	}
	start, end, strand, startE, endE := GetMPIInfo(mpi)
	if startE == endE {
		endE++
	}
	s = fmt.Sprintf("EdgeID:%d R[%d--%d] E[%d--%d] std:%t eLen:%d Len:%d score:%0.1f Percent:%0.3f", mpi.EID, start, end, startE, endE, strand, edgesArr[mpi.EID].GetSeqLen(), endE-startE, mpi.Score, mpi.Score/float32(endE-startE))
	return s
}

func PrintMKIArr(mkiArr []MapingKmerInfo) {
	for i, mki := range mkiArr {
		fmt.Printf("[PrintMKIArr]arr[%d] EID:%v strand:%v RPos:%d EPos:%d Len:%d\n", i, mki.EID, mki.GetStrand() == mki.GetEStrand(), mki.GetPosition(), mki.GetEPosition(), mki.GetMapLen())
	}
}

func GetBubbleEdgeIDArr(path []uint32, edgesArr []DBGEdge, nodesArr []DBGNode) (bubArr []uint32, idxArr []int) {
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

func IsRepeatInEPIPath(bubArr []uint32) bool {
	for i, eID := range bubArr {
		for j := i + 1; j < len(bubArr); j++ {
			if eID == bubArr[j] {
				return true
			}
		}
	}
	return false
}

func GetEdgeMappingStrand(eID uint32, strd bool, path []uint32, direction uint8, edgesArr []DBGEdge, nodesArr []DBGNode) (strand bool) {
	strand = strd
	e := &edgesArr[path[0]]
	ok := false
	if eID == path[0] {
		return strand
	}
	for i := 1; i < len(path); i++ {
		ne := &edgesArr[path[i]]
		strand = GetNextEdgeStrand2(e, ne, strand, direction)
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

type CyclePath struct {
	EID  uint32
	Idx  int
	Path []uint32
}

func GetTwoEdgeCycleArr(path []uint32, edgesArr []DBGEdge, nodesArr []DBGNode, cycleSize, TwoEdgeCyleMaxLen int) (cyclePathArr []CyclePath) {
	for i, eID := range path {
		e := edgesArr[eID]
		if e.GetTwoEdgesCycleFlag() > 0 && SumTwoEdgeCycleLen(e, edgesArr, nodesArr) < TwoEdgeCyleMaxLen {
			if i > 0 && edgesArr[path[i-1]].GetTwoEdgesCycleFlag() > 0 {
				continue
			}
			if i < len(path)-1 && edgesArr[path[i+1]].GetTwoEdgesCycleFlag() > 0 {
				continue
			}

			neID := GetOtherEdgeInTwoEdgesCycle(edgesArr, nodesArr, eID)
			var cp CyclePath
			cp.EID = eID
			cp.Idx = i
			cp.Path = append(cp.Path, eID)
			cyclePathArr = append(cyclePathArr, cp)
			for j := 0; j < cycleSize; j++ {
				var tcp CyclePath
				tcp.EID = eID
				tcp.Idx = i
				tcp.Path = append(tcp.Path, cp.Path...)
				tcp.Path = append(tcp.Path, neID)
				tcp.Path = append(tcp.Path, eID)
				cyclePathArr = append(cyclePathArr, tcp)
				cp = tcp
			}
		}
	}
	return
}

func GetEqualPath(path []uint32, bubPathArr []BubblePath) (idx int) {
	idx = -1
	for i, bp := range bubPathArr {
		if len(bp.Path) <= len(path) && reflect.DeepEqual(path[:len(bp.Path)], bp.Path) {
			idx = i
			break
		}
	}
	return
}

func GetDiffPathLenBubbleArr(path []uint32, edgesArr []DBGEdge, nodesArr []DBGNode, MaxPathLen, MaxBubbleSeqLen int, kmerlen int) (diffLenBubArr []BubblePath) {
	//MaxBubbleSeqLen := 2*kmerlen + kmerlen/2
	//fmt.Printf("[GetDiffPathLenBubbleArr]path: %v\n", path)
	for i := 1; i < len(path); {
		eID := path[i]
		e := &edgesArr[eID]
		if e.StartNID == e.EndNID || e.GetTwoEdgesCycleFlag() > 0 || e.GetBubbleFlag() > 0 || e.GetBubbleRepeatFlag() > 0 || e.GetSeqLen() > MaxBubbleSeqLen {
			i++
			continue
		}
		le := &edgesArr[path[i-1]]
		var nd DBGNode
		if le.EndNID == e.StartNID || le.EndNID == e.EndNID {
			nd = nodesArr[le.EndNID]
		} else {
			nd = nodesArr[le.StartNID]
		}
		pathArr := GetNextPathArr(le, nd, MaxPathLen+1, MaxBubbleSeqLen, edgesArr, nodesArr, kmerlen)
		/*for x, p := range pathArr {
			fmt.Printf("[GetDiffPathLenBubbleArr]pathArr[%v]: %v\n", x, p)
		}*/
		bubPathArr, _ := GetBubblePathArr(pathArr, edgesArr, kmerlen, MaxPathLen, MaxBubbleSeqLen)
		/*for x, bp := range bubPathArr {
			fmt.Printf("[GetDiffPathLenBubbleArr]bubPathArr[%v]: %v\n", x, bp)
		}*/
		if len(bubPathArr) > 1 {
			//fmt.Printf("[GetDiffPathLenBubbleArr]path[%v:]: %v\n", i, path[i:])
			idx := GetEqualPath(path[i:], bubPathArr)
			if idx >= 0 {
				var bp BubblePath
				bp.Idx = i
				bp.Path = bubPathArr[idx].Path
				diffLenBubArr = append(diffLenBubArr, bp)
				for x, ele := range bubPathArr {
					if x == idx {
						continue
					}
					var tbp BubblePath
					tbp.Idx = i
					tbp.Path = ele.Path
					diffLenBubArr = append(diffLenBubArr, tbp)
				}
				i += len(bp.Path)
				continue
			}
		}
		i++
	}
	return
}

func GetExtendFlank(le DBGEdge, direction uint8, strand bool, flankExtendLen int) (startPos, extLen int) {
	if le.GetSeqLen() > flankExtendLen {
		extLen = flankExtendLen
	} else {
		extLen = le.GetSeqLen()
	}
	if direction == FORWARD {
		if strand == PLUS {
			startPos = le.GetSeqLen() - extLen
		} else {
			startPos = extLen
		}
	} else {
		if strand == PLUS {
			startPos = extLen
		} else {
			startPos = le.GetSeqLen() - extLen
		}
	}
	return
}

func BubbleAlignment(epi ExtendPathInfo, edgesArr []DBGEdge, nodesArr []DBGNode, rd ReadInfo, SeedLen, kmerlen int, direction uint8, TwoEdgeCyleMaxLen, MaxPathLen, MaxBubbleSeqLen, WinSize int, maxPathInfoArr MaxPathInfoArr) ExtendPathInfo {
	// Get bubble edges
	bubArr, idxArr := GetBubbleEdgeIDArr(epi.Path, edgesArr, nodesArr)
	if IsRepeatInEPIPath(bubArr) {
		fmt.Printf("[BubbleAlignment] this Path: %v encounter one bubble display >2 times\n", epi.Path)
		return epi
	}

	// Get two Edge cycle
	CycleSize := 1
	flankExtendLen := 800
	var cyclePathArr []CyclePath
	cyclePathArr = GetTwoEdgeCycleArr(epi.Path, edgesArr, nodesArr, CycleSize, TwoEdgeCyleMaxLen)

	// Get Diff Edge length bubble, for example , one edge and two edge produce complex bubble
	//MaxPathLen := 3
	diffLenBubArr := GetDiffPathLenBubbleArr(epi.Path, edgesArr, nodesArr, MaxPathLen, MaxBubbleSeqLen, kmerlen)

	bubNum := len(bubArr)
	if len(bubArr) == 0 && len(cyclePathArr) == 0 && len(diffLenBubArr) == 0 {
		return epi
	}
	if DebugModel {
		fmt.Printf("[BubbleAlignment] this Path: %v has bubble: %v, idxArr: %v\n", epi.Path, bubArr, idxArr)
		if len(cyclePathArr) > 0 {
			fmt.Printf("[BubbleAlignment] cyclePathArr: %v\n", cyclePathArr)
		}
		if len(diffLenBubArr) > 0 {
			fmt.Printf("[BubbleAlignment] diffLenBubArr: %v\n", diffLenBubArr)
		}
	}
	path1 := make([]uint32, len(epi.Path))
	mpi := epi.MPI
	mkia, mkib := mpi.Kb[mpi.Base+uint32(mpi.Arr[0])], mpi.Kb[mpi.Base+uint32(mpi.Arr[len(mpi.Arr)-1])]
	var readRegionSeq []byte
	var start, end int
	if direction == FORWARD {
		start, end = int(mkia.GetPosition()), int(mkib.GetPosition())+int(mkib.GetMapLen())
	} else {
		end, start = int(mkia.GetPosition())+int(mkia.GetMapLen()), int(mkib.GetPosition())
	}
	if start > flankExtendLen {
		start -= flankExtendLen
	} else {
		start = 0
	}
	if end < len(rd.Seq)-flankExtendLen {
		end += flankExtendLen
	} else {
		end = len(rd.Seq)
	}
	readRegionSeq = rd.Seq[start:end]

	edgesSeqArr := make([][]byte, len(bubArr)+len(cyclePathArr)+len(diffLenBubArr))
	// add bubble edges seq
	{
		var strand bool
		for i, eID := range bubArr {
			if eID == epi.Path[idxArr[i]] {
				strand = GetEdgeMappingStrand(eID, mkia.GetStrand(), epi.Path, direction, edgesArr, nodesArr)
			} else {
				e := edgesArr[epi.Path[idxArr[i]]]
				ne := edgesArr[eID]
				if ne.StartNID == e.EndNID && ne.EndNID == e.StartNID {
					strand = !strand
				}
			}
			if strand == PLUS {
				edgesSeqArr[i] = append(edgesSeqArr[i], edgesArr[eID].Ks...)
			} else {
				edgesSeqArr[i] = append(edgesSeqArr[i], GetReverseCompByteArr(edgesArr[eID].Ks)...)
			}
		}
	}
	var twoEdgeCycleIdxArr []int
	// add two edge cycle seq
	{
		var strand bool
		for i, cp := range cyclePathArr {
			var path []uint32
			path = append(path, epi.Path[:cp.Idx]...)
			path = append(path, cp.Path...)
			if cp.Idx < len(epi.Path)-1 {
				path = append(path, epi.Path[cp.Idx+1:]...)
			}
			//e := edgesArr[epi.Path[cp.Idx]]
			le := edgesArr[epi.Path[cp.Idx-1]]
			strand = GetEdgeMappingStrand(le.ID, mkia.GetStrand(), epi.Path, direction, edgesArr, nodesArr)
			startPos, extLen := GetExtendFlank(le, direction, strand, flankExtendLen)
			for _, id := range cp.Path {
				extLen += len(edgesArr[id].Ks) - (kmerlen - 1)
			}
			if len(epi.Path) > cp.Idx+1 {
				ne := edgesArr[epi.Path[cp.Idx+1]]
				if len(ne.Ks)-(kmerlen-1) > flankExtendLen {
					extLen += flankExtendLen
				} else {
					extLen += len(ne.Ks) - (kmerlen - 1)
				}
			}
			if DebugModel {
				fmt.Printf("[BubbleAlignment] cyclePathArr[%d] extLen:%d path: %v\n", i, extLen, path)
			}
			es := GetDirectionEdgesSeqLen(edgesArr, nodesArr, kmerlen, strand, startPos, extLen, path[cp.Idx-1:], le.ID)
			/*if direction == BACKWARD {
				es = GetReverseByteArr(es)
			}*/
			edgesSeqArr[bubNum+i] = es
			twoEdgeCycleIdxArr = append(twoEdgeCycleIdxArr, bubNum+i)
		}
	}

	// add diff path length  Bubble Path seq
	{
		idx := -1
		var strand bool
		for i, bp := range diffLenBubArr {
			le := &edgesArr[epi.Path[bp.Idx-1]]
			if idx < bp.Idx {
				idx = bp.Idx
				strand = GetEdgeMappingStrand(le.ID, mkia.GetStrand(), epi.Path, direction, edgesArr, nodesArr)
			} else if idx == bp.Idx {

			} else {
				log.Fatalf("[BubbleAlignment] DiffLenBubArr: %v, set error\n", diffLenBubArr)
			}
			ne := &edgesArr[bp.Path[0]]
			std := GetNextEdgeStrand2(le, ne, strand, direction)
			//extLen := GetPathSeqLen(bp.Path, edgesArr, kmerlen)
			es := GetDirectionEdgesSeq(edgesArr, nodesArr, kmerlen, std, bp.Path)
			/*if direction == BACKWARD {
				es = GetReverseByteArr(es)
			}*/
			edgesSeqArr[bubNum+len(cyclePathArr)+i] = es
		}
	}
	//fmt.Printf("[BubbleAlignment] readRegionSeq: %v\n", readRegionSeq)
	//fmt.Printf("[BubbleAlignment] edgesSeqArr: %v\n", edgesSeqArr)
	cgArr := AlignmentBlocksBubble(readRegionSeq, edgesSeqArr, SeedLen, WinSize, maxPathInfoArr, bubArr, edgesArr, start)
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
		var eIDArr []uint32
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
	diffLen := 0
	// process two edge cycle arr
	var bestCPArr []CyclePath
	for i := 0; i < len(cyclePathArr); {
		max := 0
		j := i
		for ; j < i+(CycleSize+1); j++ {
			cg := cgArr[bubNum+j]
			score := int(cg.Mch) - int(cg.Mis) - int(cg.Ins) - int(cg.Del)
			if score > max {
				max = score
				bestCPArr = bestCPArr[:0]
				bestCPArr = append(bestCPArr, cyclePathArr[j])
			} else if score == max {
				bestCPArr = append(bestCPArr, cyclePathArr[j])
			}
		}
		idx := IndexEIDAllowDiffLen(epi.Path, cyclePathArr[i].EID, cyclePathArr[i].Idx, diffLen)
		if idx < 0 {
			log.Fatalf("[BubbleAlignment] idx: %v, eID: %v not include path: %v, diffLen: %v\n", idx, cyclePathArr[i].EID, epi.Path, diffLen)
		}
		if len(bestCPArr) == 1 {
			var p0, p1 []uint32
			cp := bestCPArr[0]
			p0 = append(p0, epi.Path[:idx]...)
			p1 = append(p1, epi.Path1[:idx]...)
			p0 = append(p0, cp.Path...)
			tp := make([]uint32, len(cp.Path))
			p1 = append(p1, tp...)
			if idx < len(epi.Path)-1 {
				p0 = append(p0, epi.Path[idx+1:]...)
				p1 = append(p1, epi.Path1[idx+1:]...)
			}
			if len(cp.Path) != len(cyclePathArr[i].Path) {
				diffLen += AbsInt(len(cp.Path) - len(cyclePathArr[i].Path))
			}
			epi.Path = p0
			epi.Path1 = p1
		}
		i = j
	}

	//fmt.Printf("[BubbleAlignment]bub and two cycle path changed Path0: %v\n\tPath1: %v\n", epi.Path, epi.Path1)
	// process diff edge length of bubble
	{
		cgIdx := bubNum + len(cyclePathArr)
		for i := 0; i < len(diffLenBubArr); {
			db := diffLenBubArr[i]
			j := i
			maxScore := math.MinInt32
			var bestDBArr []BubblePath
			for ; j < len(diffLenBubArr); j++ {
				if diffLenBubArr[j].Idx == db.Idx {
					cg := cgArr[cgIdx+j]
					score := int(cg.Mch) - int(cg.Mis) - int(cg.Ins) - int(cg.Del)
					if score > maxScore {
						maxScore = score
						bestDBArr = bestDBArr[:0]
						bestDBArr = append(bestDBArr, diffLenBubArr[j])
					} else if score == maxScore {
						bestDBArr = append(bestDBArr, diffLenBubArr[j])
					}
				} else {
					break
				}
			}
			idx := IndexEIDAllowDiffLen(epi.Path, db.Path[0], db.Idx, diffLen)
			if idx < 0 {
				log.Fatalf("[BubbleAlignment] idx: %v, eID: %v not include path: %v, diffLen: %v\n", idx, db.Path[0], epi.Path, diffLen)
			}
			if len(bestDBArr) == 1 {
				var p0, p1 []uint32
				bp := bestDBArr[0]
				p0 = append(p0, epi.Path[:idx]...)
				p1 = append(p1, epi.Path1[:idx]...)
				p0 = append(p0, bp.Path...)
				tp := make([]uint32, len(bp.Path))
				p1 = append(p1, tp...)
				if idx+len(db.Path) < len(epi.Path) {
					p0 = append(p0, epi.Path[idx+len(db.Path):]...)
					p1 = append(p1, epi.Path1[idx+len(db.Path):]...)
				}
				if len(db.Path) != len(bp.Path) {
					diffLen += AbsInt(len(db.Path) - len(bp.Path))
				}
				epi.Path = p0
				epi.Path1 = p1
			} else if len(bestDBArr) > 1 {
				var p0, p1 []uint32
				bp0 := bestDBArr[0]
				bp1 := bestDBArr[1]
				p0 = append(p0, epi.Path[:idx]...)
				p1 = append(p1, epi.Path1[:idx]...)
				p0 = append(p0, bp0.Path...)
				tp := make([]uint32, len(bp0.Path))
				tp[0] = bp1.Path[0]
				tp[len(tp)-1] = bp1.Path[len(bp1.Path)-1]
				p1 = append(p1, tp...)
				if idx+len(db.Path) < len(epi.Path) {
					p0 = append(p0, epi.Path[idx+len(db.Path):]...)
					p1 = append(p1, epi.Path1[idx+len(db.Path):]...)
				}
				if len(db.Path) != len(bp0.Path) {
					diffLen += AbsInt(len(db.Path) - len(bp0.Path))
				}
				epi.Path = p0
				epi.Path1 = p1
			}
			i = j
			//fmt.Printf("[BubbleAlignment]diffbubblePath idx: %v changed Path0: %v\n\tPath1: %v\n", db.Idx, epi.Path, epi.Path1)
		}
	}
	if DebugModel {
		fmt.Printf("[BubbleAlignment]changed Path0: %v\n\tPath1: %v\n", epi.Path, epi.Path1)
	}
	return epi
}

func GetMaxScoreEPI(stk []ExtendPathInfo) (t ExtendPathInfo, idx int) {
	//max := 0 // max score in the stack
	minMapLen := uint32(math.MaxInt32)
	idx = -1
	// choose minMapLen
	for x, ele := range stk {
		if ele.Score <= 0 {
			continue
		}
		if minMapLen > ele.ExtendLen {
			t = ele
			idx = x
			minMapLen = ele.ExtendLen
		} else if minMapLen == ele.ExtendLen {
			if t.Score < ele.Score {
				t = ele
				idx = x
			}
		}
	}
	/*for x, ele := range stk {
		if ele.Score <= 0 {
			continue
		}
		mapLen := ele.ExtendLen
		fraction := ele.Score * 100 / mapLen
		if fraction > max {
			max = fraction
			t = ele
			idx = x
			minMapLen = ele.ExtendLen
		} else if fraction == max {
			if minMapLen > ele.ExtendLen {
				t = ele
				idx = x
				minMapLen = ele.ExtendLen
			}
		}
	}*/
	return
}

func GetMaxScoreEPI2(stk []ExtendPathInfo) (t ExtendPathInfo, idx int) {
	max := uint32(0) // max score in the stack
	idx = -1
	extendLen := uint32(0)
	// choose max score
	for x := range stk {
		if stk[x].Score <= 0 {
			continue
		}

		if stk[x].Score > max {
			t = stk[x]
			idx = x
			max = t.Score
			extendLen = t.ExtendLen
		} else if stk[x].Score == max {
			if stk[x].ExtendLen < extendLen {
				t = stk[x]
				idx = x
				extendLen = t.ExtendLen
			}
		}
	}
	return
}

func CleanExtendPathInfoArr(stk []ExtendPathInfo) []ExtendPathInfo {
	pos := 0
	for _, ele := range stk {
		if ele.Score <= 0 {
			continue
		}
		stk[pos] = ele
		pos++
	}
	stk = stk[:pos]
	return stk
}

func SumTwoEdgeCycleLen(ne DBGEdge, edgesArr []DBGEdge, nodesArr []DBGNode) (sum int) {
	oeID := GetOtherEdgeInTwoEdgesCycle(edgesArr, nodesArr, ne.ID)
	sum = ne.GetSeqLen() + edgesArr[oeID].GetSeqLen()
	return
}

func GetNoSeedPathLen(epi ExtendPathInfo, edgesArr []DBGEdge, kmerlen int) int {
	fl := int(epi.FlankLen)
	pl := 0 // path length
	extendEPI := epi.ExtendEPI
	for extendEPI.Last != nil {
		extendEPI = *extendEPI.Last
		eID := extendEPI.EID
		sl := edgesArr[eID].GetSeqLen() - (kmerlen - 1)
		if fl >= sl {
			fl -= sl
			pl++
		} else {
			break
		}
	}
	return pl
}

func CheckExtend(nt ExtendPathInfo, MinChainScoreIdentityPercent int) (ok bool) {
	mapLen := nt.ExtendLen - nt.FlankLen
	if mapLen < 300 {
		ok = true
	} else {
		ok = mapLen*uint32(MinChainScoreIdentityPercent)/100 < nt.Score
	}
	return
}

func GetInfoStrand(info uint32) bool {
	if (info & 0x1) > 0 {
		return PLUS
	} else {
		return MINUS
	}
}

const MUSK16 = (1 << 16) - 1
const BASE16 = (1 << 16)

/*func SetKmerNum(info uint32, strand bool) uint32 {
	if strand == PLUS {
		t := (info & MUSK16)
		if t < math.MaxUint16 {
			t++
		} else {
			fmt.Printf("[SetKmerNum] t:%d\n", t)
		}
		info >>= 16
		info = (info << 16) | t
	} else {
		t := (info >> 16)
		if t < math.MaxUint16 {
			t++
		} else {
			fmt.Printf("[SetKmerNum] t:%d\n", t)
		}
		info = (t << 16) | (info & MUSK16)
	}
	return info
}*/

func SetKmerNum(info uint32, strand bool) uint32 {
	if strand {
		return info + 1
	} else {
		return ((((info >> 16) + 1) << 16) | (info & MUSK16))
	}
}

func GetKmerNum(info uint32, strand bool) uint32 {
	if strand == PLUS {
		return (info & MUSK16)
	} else {
		return (info >> 16)
	}
}

func DeleteLowNumMinimersEdge(edgesArr []DBGEdge, edgeMiniKmerNumMap map[uint32]uint32, LongEdgeMinLen, MinChainScoreIdentityPercent int) map[uint32]uint32 {
	//eIDArr := make([]uint32, 0, 1000)
	totalNum, totalMinimer := 0, 0
	newMap := make(map[uint32]uint32, 200)
	for k, v := range edgeMiniKmerNumMap {
		e := edgesArr[k]
		if e.CovD < 1 {
			continue
		}
		/*if e.GetLowFreqMiniMersFlag() > 0 {
			continue
		}*/
		num1, num2 := GetKmerNum(v, PLUS), GetKmerNum(v, MINUS)
		num := MaxInt(int(num1), int(num2))
		el := e.GetSeqLen()
		if el > LongEdgeMinLen {
			if float32(num) >= float32(e.CovD)/(float32(el)/float32(LongEdgeMinLen))/7 {
				newMap[k] = v
				//delete(edgeMiniKmerNumMap, k)
			}
		} else {
			if num >= int(e.CovD)/7 {
				totalNum += num
				totalMinimer += int(e.CovD)
				newMap[k] = v
				//delete(edgeMiniKmerNumMap, k)
			} //else {
			//fmt.Printf("[DeleteLowNumMinimersEdge]el:%d minimerNum: %d mapNum:%d\n", el, e.CovD, num)
			//}
		}
	}

	fmt.Printf("[DeleteLowNumMinimersEdge] avg percent: %v\n", float32(totalNum)/float32(totalMinimer+1))

	return newMap
}

func DeleteLowNumMinimersEdgePercent(edgesArr []DBGEdge, edgeMiniKmerNumMap map[uint32]uint32, LongEdgeMinLen, MinChainScoreIdentityPercent int) map[uint32]uint32 {
	//eIDArr := make([]uint32, 0, 1000)
	totalNum, totalLen := 0, 0
	newMap := make(map[uint32]uint32, 100)
	for k, v := range edgeMiniKmerNumMap {
		e := edgesArr[k]
		num1, num2 := GetKmerNum(v, PLUS), GetKmerNum(v, MINUS)
		num := MaxInt(int(num1), int(num2))
		el := e.GetSeqLen()
		if el > LongEdgeMinLen {
			if num >= LongEdgeMinLen/6/7 {
				newMap[k] = v
				//delete(edgeMiniKmerNumMap, k)
			}
		} else {
			if num >= el/6/7 {
				totalNum += num
				totalLen += el
				newMap[k] = v
				//delete(edgeMiniKmerNumMap, k)
			}
			//fmt.Printf("[DeleteLowNumMinimersEdgePercent]el:%d minimerNum: %d mapNum:%d\n", el, e.CovD, num)
		}
	}
	fmt.Printf("[DeleteLowNumMinimersEdgePercent] avg minimers num: %v\n", float32(totalNum)/float32(totalLen+1))

	return newMap
}

type SeedNumInfo struct {
	//EdgeInfo           uint32 // high 12-bits for EdgeLen, middle 10-bits for highest pos of edge, low 10-bits for lowest pos of edge
	SeedNumP, SeedNumM uint8 // SeedNumP note PLUS strand, SeedNumM note MINUS strand
	EdgeLowFreqSeedNum uint8
}

/*const EdgeLenBits = 12
const MaxEdgeLen = 1<<EdgeLenBits - 1
const EdgePosBits = 10
const MaxEdgePos = 1<<EdgePosBits - 1
const EdgeLenMusk = 1<<(EdgePosBits+EdgePosBits) - 1
const EdgePosMusk = 1<<EdgePosBits - 1
const EdgeHighMusk = (1<<EdgeLenBits-1)<<(EdgePosBits+EdgePosBits) | EdgePosMusk
const PosMusk = MaxEdgeLen << (EdgePosBits + EdgePosBits)

func (sni SeedNumInfo) GetEdgeLen() uint32 {
	return sni.EdgeInfo >> 20
}
func (sni *SeedNumInfo) SetEdgeLen(el uint32) {
	sni.EdgeInfo = (el << 20) | (sni.EdgeInfo & EdgeLenMusk)
}

func (sni SeedNumInfo) GetLowPos() uint32 {
	return sni.EdgeInfo & EdgePosMusk
}

func (sni *SeedNumInfo) SetLowPos(pos uint32) {
	sni.EdgeInfo = ((sni.EdgeInfo >> EdgePosBits) << EdgePosBits) | pos
}

func (sni SeedNumInfo) GetHighPos() uint32 {
	return (sni.EdgeInfo >> EdgePosBits) & EdgePosMusk
}

func (sni *SeedNumInfo) SetHighPos(pos uint32) {
	sni.EdgeInfo = (sni.EdgeInfo & EdgeHighMusk) | (pos << EdgeLenBits)
}
func (sni *SeedNumInfo) ResetPos() {
	sni.EdgeInfo = sni.EdgeInfo & PosMusk
} */

type EdgeStrand struct {
	EdgeID uint32
	Strand bool
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

type EdgeMKINum struct {
	EdgeID uint32
	Num    uint32
}

func Transform2EdgeMKINumArr(esArr []EdgeStrand) []EdgeMKINum {
	enArr := make([]EdgeMKINum, len(esArr))
	for i, es := range esArr {
		enArr[i].EdgeID = es.EdgeID
	}

	return enArr
}

func IndexEdgeMKINumArr(enArr []EdgeMKINum, eID uint32) (idx int) {
	idx = -1
	i := len(enArr) * int(eID) / int(enArr[len(enArr)-1].EdgeID)
	if i >= len(enArr) {
		i = len(enArr) - 1
	}
	if enArr[i].EdgeID == eID {
		idx = i
	} else if enArr[i].EdgeID < eID {
		for j := i + 1; j < len(enArr); j++ {
			if enArr[j].EdgeID > eID {
				break
			} else if enArr[j].EdgeID == eID {
				idx = j
				break
			}
		}
	} else {
		for j := i - 1; j >= 0; j-- {
			if enArr[j].EdgeID < eID {
				break
			} else if enArr[j].EdgeID == eID {
				idx = j
				break
			}
		}
	}

	return idx
}

func DeleteLowSeedNumMinimersEdge(SeedNumCountArr []SeedNumInfo, LongEdgeMinLen uint32, WinSize int, kmerlen, MQ int) ([]EdgeStrand, float32) {
	totalNum := 0
	SeedNumSum := float32(0)
	longEdgeNum := 0
	//EdgeMinLen := uint32((kmerlen-1)*2) - 10
	//minNum := uint8(EdgeMinLen * 2 / WinSize * 4 / 5)
	a := 0
	dif := uint8(20)
	if MQ < 20 {
		dif += (dif - uint8(MQ))
	}
	MinPerCent := float32(0.2)
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
		totalNum++
		//num1, num2 := GetKmerNum(ct, PLUS), GetKmerNum(ct, MINUS)
		//el := sni.GetEdgeLen()
		edgeLowFreqSeedNum := sni.EdgeLowFreqSeedNum

		/*if num > 10000 {
			fmt.Printf("[DeleteLowSeedNumMinimersEdge]num1:%d num2:%d el:%d\n", num1, num2, el)
		}*/
		if num > 20 {
			longEdgeNum++
			var es EdgeStrand
			es.EdgeID = uint32(i)
			if num == numP {
				es.Strand = PLUS
			} else {
				es.Strand = MINUS
			}
			esArr = append(esArr, es)
			//fmt.Printf("[DeleteLowSeedNumMinimersEdge]el:%d SeedNum:%d edgeLowFreqSeedNum:%d\n", el, num, edgeLowFreqSeedNum)
		} else {
			if num > 0 && num > uint8(float32(edgeLowFreqSeedNum)*MinPerCent) {
				var es EdgeStrand
				es.EdgeID = uint32(i)
				if num == numP {
					es.Strand = PLUS
				} else {
					es.Strand = MINUS
				}
				esArr = append(esArr, es)
				SeedNumSum += float32(num) / float32(edgeLowFreqSeedNum)
				a++
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
	}

	//longEdgePercent := uint16(float32(1) / (float32(longEdgeSeedNumSum) / float32(longEdgeNum)) * 0.9 * 100)
	if DebugModel {
		fmt.Printf("[DeleteLowSeedNumMinimersEdge] totalNum:%d longEdgeNum:%d after delete:%d\n", totalNum, longEdgeNum, len(esArr))
	}
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
	var sp float32
	if a == 0 {
		sp = 0.22
	} else {
		sp = SeedNumSum / float32(a)
	}
	return esArr, sp
}

func AddHighFreqKmer(ka []MapingKmerInfo, dbgKmerInfoPac DBGKmerInfoPac, edgesArr []DBGEdge, MinChainScoreIdentityPercent, SeedLen, LongEdgeMinLen, readLen int) []MapingKmerInfo {
	edgeMiniKmerNumMap := make(map[uint32]uint32, 5000)
	//LongEdgeMinLen := 1000
	highFreqNum := 0
	//MinChainScoreIdentityPercent *= 2
	for _, mki := range ka {
		if mki.EID == math.MaxUint32 {
			highFreqNum++
			continue
		}
		if mki.GetStrand() == PLUS {
			edgeMiniKmerNumMap[mki.EID] = SetKmerNum(edgeMiniKmerNumMap[mki.EID], PLUS)
		} else {
			edgeMiniKmerNumMap[mki.EID] = SetKmerNum(edgeMiniKmerNumMap[mki.EID], MINUS)
		}
	}
	beforeDN := len(edgeMiniKmerNumMap)
	edgeMiniKmerNumMap = DeleteLowNumMinimersEdge(edgesArr, edgeMiniKmerNumMap, LongEdgeMinLen, MinChainScoreIdentityPercent)
	fmt.Printf("[AddHighFreqKmer]highFreqNum:%d beforeDN:%d len(edgeMiniKmerNumMap):%d\n", highFreqNum, beforeDN, len(edgeMiniKmerNumMap))

	for _, mki := range ka {
		if mki.EID == math.MaxUint32 {
			eIDArr, ok := dbgKmerInfoPac.LowFreqPercentEdgeMap[mki.GetKmer()]
			if !ok {
				log.Fatalf("[AddHighFreqKmer] mki:%v not include LowFreqPercentEdgeMap\n", PrintMKI(mki))
			}
			for _, eID := range eIDArr {
				v := SetKmerNum(edgeMiniKmerNumMap[eID], PLUS)
				edgeMiniKmerNumMap[eID] = SetKmerNum(v, MINUS)
			}
		}
	}

	beforeDN = len(edgeMiniKmerNumMap)
	edgeMiniKmerNumMap = DeleteLowNumMinimersEdgePercent(edgesArr, edgeMiniKmerNumMap, LongEdgeMinLen, MinChainScoreIdentityPercent)
	fmt.Printf("[AddHighFreqKmer]DeleteLowNumMinimersEdgePercent beforeDN:%d len(edgeMiniKmerNumMap):%d\n", beforeDN, len(edgeMiniKmerNumMap))

	aka := make([]MapingKmerInfo, 0, readLen/3+5)
	for _, mki := range ka {
		if mki.EID == math.MaxUint32 {
			eIDArr, ok := dbgKmerInfoPac.LowFreqPercentEdgeMap[mki.GetKmer()]
			if !ok {
				log.Fatalf("[AddHighFreqKmer] mki:%v not include LowFreqPercentEdgeMap\n", PrintMKI(mki))
			}
			for _, eID := range eIDArr {
				if _, ok1 := edgeMiniKmerNumMap[eID]; ok1 {
					arr := dbgKmerInfoPac.HighFreqMap[mki.GetKmer()]
					for _, idInfo := range arr {
						var nmki MapingKmerInfo
						nmki.SeedInfo = mki.SeedInfo
						nmki.EID = eID
						nmki.SetEPosition(idInfo.GetPosition())
						nmki.SetEStrand(idInfo.GetStrand())
						aka = append(aka, nmki)
					}
				}
			}
		} else {
			if _, ok1 := edgeMiniKmerNumMap[mki.EID]; ok1 {
				aka = append(aka, mki)
			}
		}
	}
	edgeMiniKmerNumMap = nil
	fmt.Printf("[AddHighFreqKmer] len(ka):%d len(aka):%d\n", len(ka), len(aka))
	return aka
}

func IsInMKIArr(ka []MapingKmerInfo, mki MapingKmerInfo) bool {
	ok := false
	for _, t := range ka {
		if t.GetKmer() == mki.GetKmer() {
			ok = true
			break
		}
	}
	return ok
}

func AddHighFreqKmer2(ka []MapingKmerInfo, dbgKmerInfoPac DBGKmerInfoPac, SeedLen uint16, LongEdgeMinLen int, SeedNumCountArr []SeedNumInfo, WinSize, kmerlen, MQ int) ([]MapingKmerInfo, []EdgeStrand, float32) {
	//MinChainScoreIdentityPercent *= 2

	var freqRepeatNum int
	highFreqMap := make(map[uint64][]MapPos)
	var lastMKI MapingKmerInfo
	for _, mki := range ka {
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
			if mki.EID == lastMKI.EID && mki.GetEPosition() == lastMKI.GetEPosition() {
				freqRepeatNum++
				continue
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
		fmt.Printf("[AddHighFreqKmer]len(highFreqMap):%d freqRepeatNum:%d  after delete:%d seedNumPercent:%.2f\n", len(highFreqMap), freqRepeatNum, len(esArr), seedNumPercent)
	}
	//fmt.Printf("[AddHighFreqKmer]esArr:%v\n", esArr)
	//aka = aka[:0]
	// EdgeSeedNumMap
	if len(esArr) == 0 || seedNumPercent < 0.2 {
		ka = ka[:0]
		return ka, esArr, seedNumPercent
	}

	idxKa := 0
	for _, mki := range ka {
		if mki.EID != math.MaxUint32 {
			ok, strand := IsInEdgeStrandArr(esArr, mki.EID)
			if ok && strand == mki.GetMapStrand() {
				ka[idxKa] = mki
				idxKa++
			}
		}
	}
	ka = ka[:idxKa]
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
					nmki.SetMapLen(uint32(SeedLen))
					//fmt.Printf("[AddHighFreqKmer]EPos:%d\n", arr[j].GetPosition())
					//nmki.SeedInfo = uint64(RInfo.Info)
					//nmki.Info = arr[j].Info;
					nmki.EID = es.EdgeID
					ka = append(ka, nmki)
					addHighFreqNum++
				}
			}
		}
	}

	if DebugModel {
		fmt.Printf("[AddHighFreqKmer]addHighFreqNum:%d\n", addHighFreqNum)
	}
	return ka, esArr, seedNumPercent
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

func IsTipEdge(e DBGEdge) bool {
	if e.StartNID == 0 || e.EndNID == 0 {
		return true
	}
	return false
}

func CleanScoreArr(maxScoreArr []ScorePos) {
	for i := range maxScoreArr {
		maxScoreArr[i].Pos = 0
		maxScoreArr[i].Score = 0
	}
}

func FindONTReadsPath(dbgKmerInfoPac DBGKmerInfoPac, cs <-chan ReadsBucketInfo, SeedLen, winSize, GapCost, MaxGapLen, kmerlen int, wc chan<- LongReadMappingInfo, edgesArr []DBGEdge, nodesArr []DBGNode, SeqType uint8, MaxPathLen, MaxBubbleSeqLen, MinChainScoreIdentityPercent int) {
	extendLen := 3 * kmerlen
	TwoEdgeCyleMaxLen := 800
	//MaxPathLen := 10
	//MaxBubbleSeqLen := kmerlen * 3
	//ExtendIdentityPercent := 100
	MaxEnlongationSeqLen := uint32(3 * kmerlen)
	//MaxEnlongationPathLen := 10
	MaxFlankLen := 2000
	MaxStackAllow := 5000
	LongEdgeMinLen := 1000
	SeedNumCountArr := GetEdgeLowFreqSeedNumArr(edgesArr, LongEdgeMinLen)
	var totalReadNum, notFoundSeedEdgeNum, mapOneEdgeNum, ExtendNum int
	bucketInterSecKmerInfoArr := make([][]MapingKmerInfo, BuckReadsNum)
	for i := 0; i < BuckReadsNum; i++ {
		bucketInterSecKmerInfoArr[i] = make([]MapingKmerInfo, 0, 200000)
	}
	aka := make([]MapingKmerInfo, 0, 40000)
	kb := make([]MapingKmerInfo, 0, 1000)
	buf := make([]MapingKmerInfo, 0, 40000)
	//EdgeSeedNumCountArr := make([]EdgeSeedNumInfo, 0, 10000)
	stk := make([]ExtendPathInfo, 0, 1000)
	StkCleanSize := 100
	maxScoreArr := make([]ScorePos, 0, 1000)
	ScoreWinSize := 50
	highQualEPathInfoArr := make([]ExtendPathInfo, 0, 10)
	maxPathInfoArr := make([]MaxPathInfo, 0, 10000)
	for {
		rb := <-cs
		if len(rb.ReadsArr) == 0 {
			var tmp LongReadMappingInfo
			wc <- tmp
			break
		}
		totalReadNum += len(rb.ReadsArr)
		startReadID := rb.ReadsArr[0].ID
		bucketInterSecKmerInfoArr = GetInterSectionKmerInfoB(dbgKmerInfoPac, rb.KmerSortArr, uint32(startReadID), uint16(SeedLen), bucketInterSecKmerInfoArr)
		// loop process every read
		for j, ka := range bucketInterSecKmerInfoArr {
			if len(ka) == 0 {
				continue
			}
			rID := uint32(startReadID) + uint32(j)
			ri := rb.ReadsArr[j]
			rl := len(ri.Seq)
			//blKa := len(ka)
			fmt.Printf("[FindONTReadsPath] read[%v] length: %v\n", rID, rl)
			aka = aka[:0]
			//PrintMKIArr(ka)
			MQ := 20
			bkaL := len(ka)
			var seedNumPercent float32
			var esArr []EdgeStrand
			ka, esArr, seedNumPercent = AddHighFreqKmer2(ka, dbgKmerInfoPac, uint16(SeedLen), LongEdgeMinLen, SeedNumCountArr, winSize, kmerlen, MQ)
			//sort.Sort(MapingKmerInfoArr(ka))
			if DebugModel {
				fmt.Printf("[FindONTReadsPath]before: %d,  after add high freq kmer: %v\n", bkaL, len(ka))
			}
			var enArr []EdgeMKINum
			kb, buf, enArr = GetChainBlocks(ka, buf, esArr)
			//kbEdgeArr := GetEdgeChainBlocks(kb, MaxGapLen)
			/*for y, mki := range kb {
				fmt.Printf("[FindONTReadsPath] kb[%v]EID: %v, ReadPos : %v, EdgePos: %v,  Len: %v\n", y, mki.EID, mki.GetPosition(), mki.GetEPosition(), mki.GetMapLen())
			}*/
			//MinKmerScore := SeedLen * 2
			maxPathInfoArr = GetMaxChainBlockArr(kb, maxPathInfoArr, edgesArr, nodesArr, SeedLen, GapCost, MaxGapLen, kmerlen, SeqType, LongEdgeMinLen, float32(MinChainScoreIdentityPercent)/100, seedNumPercent, os.Stdout)
			//sort.Sort(MapingKmerInfoArr(kb))
			fmt.Printf("[FindONTReadsPath]len(kb):%d len(maxPathInfoArr):%d\n", len(kb), len(maxPathInfoArr))
			//sort.Sort(MaxPathInfoTmpArr(maxPathInfoArr))
			//var cBMRIArr []ChainBlockMapLongReadInfo
			/*if DebugModel {
				for x, mpi := range maxPathInfoArr {
					mka, mkb := mpi.Kb[mpi.Base+int(mpi.Arr[0])], mpi.Kb[mpi.Base+int(mpi.Arr[mpi.GetMapLen()Arr-1])]
					start, end := mka.GetPosition(), mkb.GetPosition()+uint32(mkb.GetMapLen())
					fmt.Printf("[FindONTReadsPath]MPIArr[%v]: EdgeID:%v ReadRegion[%v--%v] eLen:%d Len:%v score:%v Percent:%v\n", x, mka.EID, start, end, edgesArr[mka.EID].GetSeqLen(), end-start, mpi.Score, mpi.Score*100/(int(end)-int(start)))
				}
			}*/

			uniqueMaxPathInfo := GetUniqueMaxPathInfo(maxPathInfoArr, kmerlen, SeedLen, rl, edgesArr, MinChainScoreIdentityPercent)
			if uniqueMaxPathInfo.Score == 0 {
				notFoundSeedEdgeNum++
				fmt.Printf("[FindONTReadsPath] read ID[%v] not found unique boundary mapping region\n", rID)
				continue
			}
			mka, mkb := uniqueMaxPathInfo.Kb[uniqueMaxPathInfo.Base+uint32(uniqueMaxPathInfo.Arr[0])], uniqueMaxPathInfo.Kb[uniqueMaxPathInfo.Base+uint32(uniqueMaxPathInfo.Arr[len(uniqueMaxPathInfo.Arr)-1])]
			start, end := int(mka.GetPosition()), int(mkb.GetPosition())+int(mkb.GetMapLen())
			if int(start) < 2*kmerlen && int(end) > rl-2*kmerlen {
				mapOneEdgeNum++
				fmt.Printf("[FindONTReadsPath] read ID[%v] only map one long edge[%v]\n", rID, mka.EID)
				continue
			}
			//for x := len(maxPathInfoArr) - 1; x >= 0; x-- {
			mpi := uniqueMaxPathInfo
			eID := mpi.EID
			edge := &edgesArr[eID]
			if edge.StartNID == edge.EndNID || (edge.GetUniqueFlag() == 0 && edge.GetBubbleFlag() == 0 && edge.GetBubbleRepeatFlag() == 0 && edge.GetSeqLen() < 2*(kmerlen-1)) {
				fmt.Printf("[FindONTReadsPath] read ID[%v] unique map edge[%v] not allow extend path\n", rID, mka.EID)
				continue
			}
			fmt.Printf("[FindONTReadsPath] uniqueMPI: EdgeID:%v ReadRegion[%v--%v] eLen:%d Len:%v score:%v Percent:%v\n", mka.EID, start, end, edge.GetSeqLen(), end-start, uniqueMaxPathInfo.Score, int(uniqueMaxPathInfo.Score)*100/(int(end)-int(start)))
			// check this region has been processed
			last := mpi.Kb[mpi.Base+uint32(mpi.Arr[len(mpi.Arr)-1])]
			//start, end := int(mpi.Arr[0].GetPosition()), int(last.GetPosition())+int(last.GetMapLen())
			//fmt.Printf("[FindONTReadsPath] Read region[%v---%v], readLen: %v, Edge region[%v---%v], EdgeLen: %v, score: %v\n", mpi.Arr[0].GetPosition(), last.GetPosition()+uint32(last.GetMapLen()), rl, mpi.Arr[0].GetEPosition(), last.GetEPosition(), len(edgesArr[eID].Ks), mpi.Score)
			if edge.GetSeqLen() < 3*kmerlen {
				fmt.Printf("[FindONTReadsPath] Seed Edge length: %d smaller than %d\n", edge.GetSeqLen(), 3*kmerlen)
			}
			//fmt.Printf("[FindONTReadsPath] uniqueMaxPathInfo: \n")
			//PrintAddedMpArr(uniqueMaxPathInfo.Arr)

			/*if mpi.Score < 5*SeedLen || IsInChainBlockMapLongReadInfoArr(cBMRIArr, start, end) {
				continue
			}*/
			/*// get base alignment for seed chain
			readRegionSeq := rb.ReadsArr[j].Seq[start:end]
			edgeRegionSeq := GetEdgeRegionSeq(edgesArr[eID], mpi)
			chainArr, _ := ComputingChainArr(edgeRegionSeq, readRegionSeq, uint32(SeedLen))
			//fmt.Printf(">edge\n%v\n", Transform2Letters(edgeRegionSeq))
			//fmt.Printf(">read\n%v\n", Transform2Letters(readRegionSeq))
			if len(chainArr) == 0 {
				continue
			}
			cg := DPLocalAlign(edgeRegionSeq, readRegionSeq, chainArr)
			avgSeqIdentityPercent := int(cg.Mch) * 100 / len(edgeRegionSeq)
			fmt.Printf("[FindONTReadsPath] Seed EdgeID: %v chainArrLen:%d mpi.Score:%v Qual:%v cg.Mch:%v len(edgeRegionSeq):%v\n", eID, len(chainArr), mpi.Score, avgSeqIdentityPercent, cg.Mch, len(edgeRegionSeq))
			if len(edgeRegionSeq) > 1000 {
				if int(cg.Mch) < len(edgeRegionSeq)*75/100 {
					continue
				}
			} else {
				if int(cg.Mch) < len(edgeRegionSeq)*80/100 {
					continue
				}
			}*/

			var bestEPIF, bestEPIB ExtendPathInfo
			// read FORWARD mapping
			{
				// check if is boundary of long read
				pos := end
				if pos+extendLen < rl {
					anchor := last
					// choose edge  boudnary anchor will miss some better alignment  because we use minimizers
					//var mpArr []MapingKmerInfo
					//mpArr = append(mpArr, anchor)
					var epi ExtendPathInfo
					id := eID
					epi.ExtendEPI.Last = nil
					epi.ExtendEPI.EID = uint32(eID)
					//epi.Path = append(epi.Path, uint32(id))
					epi.SetStrandFlag(anchor.GetStrand())
					var coming uint8
					epi.SetComingFlag(GetEdgeComing(edgesArr[id], FORWARD, epi.GetStrandFlag(), nodesArr, coming))
					if anchor.GetStrand() == PLUS {
						epi.ExtendLen = uint32(len(edgesArr[id].Ks)) - uint32(anchor.GetEPosition())
						epi.FlankLen = uint32(len(edgesArr[id].Ks) - (int(anchor.GetEPosition()) + int(anchor.GetMapLen())))
					} else {
						epi.ExtendLen = anchor.GetEPosition() + uint32(anchor.GetMapLen())
						epi.FlankLen = anchor.GetEPosition()
					}
					epi.MPI = uniqueMaxPathInfo
					//mpiArr := make([]uint8, len(uniqueMaxPathInfo.Arr))
					epi.MPI.Arr = uniqueMaxPathInfo.Arr
					epi.Score = uint32(anchor.GetMapLen())
					//maxScore := epi.Score
					//scoreLen := int(anchor.GetMapLen())
					//maxMapLen := int(anchor.GetMapLen())
					stk = stk[:0]
					stk = append(stk, epi)
					stkAddNum := 1
					addNum := 0
					highQualEPathInfoArr = highQualEPathInfoArr[:0]
					CleanScoreArr(maxScoreArr)
					maxScoreArr = maxScoreArr[:0]
					startPos := int(anchor.GetPosition())
					for {
						// found max score in the stack
						if stkAddNum/StkCleanSize < (stkAddNum+addNum)/StkCleanSize {
							stk = CleanExtendPathInfoArr(stk)
						}
						stkAddNum += addNum
						addNum = 0
						t, idx := GetMaxScoreEPI2(stk)
						if idx >= 0 {
							stk[idx].Score = 0
						}
						if t.Score <= 0 {
							break
						}
						if len(stk) > MaxStackAllow {
							fmt.Fprintf(os.Stderr, "[FindONTReadsPath]readID:%d maxScoreArr:%v\n", ri.ID, maxScoreArr)
							break
						}
						var sa ScorePos
						//sa.ExtendLen = t.ExtendLen

						var etLen int
						te := edgesArr[t.ExtendEPI.EID]
						if te.GetSeqLen() > kmerlen+EdgeMapKmerMin {
							etLen = int(t.ExtendLen)
						}
						mpi1 := &t.MPI
						mkib := mpi1.Kb[mpi1.Base+uint32(mpi1.Arr[len(mpi1.Arr)-1])]
						sa.Score, sa.Pos = uint32(t.Score), uint32(mkib.GetPosition())+mkib.GetMapLen()
						if BiggerThanMaxScoreArr(maxScoreArr, startPos, sa, etLen, FORWARD, ScoreWinSize) == false {
							continue
						}
						// check if has NGS Path
						{
							id := t.ExtendEPI.EID
							le := &edgesArr[id]
							if len(le.NGSPathArr) == 1 {
								path := EdgeFreqArr2Uint32Arr(le.NGSPathArr[0])
								index := IndexUint32(path, id)
								var np []uint32
								if t.GetStrandFlag() == PLUS {
									if len(path)-index >= 3 {
										np = path[index+1:]
									}
								} else {
									if index >= 2 {
										np = GetReverseUint32Arr(path[:index])
									}
								}

								if len(np) >= 2 && edgesArr[np[0]].GetTwoEdgesCycleFlag() == 0 && edgesArr[np[len(np)-1]].GetTwoEdgesCycleFlag() == 0 && edgesArr[np[len(np)-1]].StartNID != edgesArr[np[len(np)-1]].EndNID {
									ok := true
									for _, eID := range np {
										var nt ExtendPathInfo
										ne := &edgesArr[eID]
										strand := GetNextEdgeStrand2(le, ne, t.GetStrandFlag(), FORWARD)
										nt.ExtendEPI.EID = eID
										nt.ExtendEPI.Last = &t.ExtendEPI
										//nt.Path = make([]uint32, len(t.Path)+1)
										//copy(nt.Path, t.Path)
										//nt.Path[len(nt.Path)-1] = eID
										nt.SetStrandFlag(strand)
										enIdx := IndexEdgeMKINumArr(enArr, uint32(eID))
										if enIdx < 0 {
											nt.MPI = t.MPI
											nt.Score = t.Score
											nt.FlankLen = t.FlankLen + uint32(len(ne.Ks)-(kmerlen-1))
										} else {
											kbIdx := enArr[enIdx].Num
											nt.MPI, nt.Score, nt.FlankLen = ExtendPathChain(t.MPI, nt.ExtendEPI, t.FlankLen, kmerlen, FORWARD, *ne, int(kbIdx), strand, t.Score, edgesArr, nodesArr, uint16(SeedLen), maxPathInfoArr, LongEdgeMinLen, MinChainScoreIdentityPercent)
										}
										nt.SetComingFlag(GetEdgeComing(edgesArr[eID], FORWARD, strand, nodesArr, t.GetComingFlag()))
										nt.ExtendLen += t.ExtendLen + uint32(len(ne.Ks)-(kmerlen-1))
										var sp ScorePos
										//sa.ExtendLen = nt.ExtendLen
										var etLen int
										te := edgesArr[nt.ExtendEPI.EID]
										if te.GetSeqLen() > kmerlen+EdgeMapKmerMin {
											etLen = int(nt.ExtendLen)
										}
										mpi1 := &nt.MPI
										mkib := mpi1.Kb[mpi1.Base+uint32(mpi1.Arr[len(mpi1.Arr)-1])]
										sp.Pos = uint32(mkib.GetPosition()) + mkib.GetMapLen()
										sp.Score = uint32(nt.Score)
										if BiggerThanMaxScoreArr(maxScoreArr, start, sp, etLen, FORWARD, ScoreWinSize) {
											if nt.Score > t.Score {
												maxScoreArr, _ = AddedToMaxScoreArrMPI(maxScoreArr, startPos, nt.MPI, int(nt.Score), FORWARD, ScoreWinSize)
											}
										} else {
											ok = false
											break
										}
										le = ne
										t = nt
									}
									if ok {
										if DebugModel {
											fmt.Printf("[FindONTReadsPath] np: %v\n", np)
										}
										mpi1 := &t.MPI
										mkib := mpi1.Kb[mpi1.Base+uint32(mpi1.Arr[len(mpi1.Arr)-1])]
										p := int(mkib.GetPosition()) + int(mkib.GetMapLen())
										if t.FlankLen > MaxEnlongationSeqLen {
											if p > rl-MaxFlankLen {
												highQualEPathInfoArr = append(highQualEPathInfoArr, t)
											}
										} else {
											if p > rl-kmerlen || IsTipEdge(edgesArr[t.ExtendEPI.EID]) {
												highQualEPathInfoArr = append(highQualEPathInfoArr, t)
											} else {
												stk = append(stk, t)
												addNum++
											}
										}
									}
									continue
								}
							}
						}
						// add next near edge info
						//id := t.Path[len(t.Path)-1]
						id := t.ExtendEPI.EID
						e := &edgesArr[id]
						eArr := GetNearDirectionEdgeIDArr(*e, FORWARD, t.GetStrandFlag(), t.GetComingFlag(), nodesArr)
						if DebugModel {
							fmt.Printf("[FindONTReadsPath]idx:%v nextEArr:%v len(t.MPI.Arr):%v Score:%v ExtendLen:%v FlankLen:%v Strand:%v Coming:%v LastEID:%d\n", idx, eArr, len(t.MPI.Arr), t.Score, t.ExtendLen, t.FlankLen, t.GetStrandFlag(), t.GetComingFlag(), t.ExtendEPI.EID)
						}
						for _, eID := range eArr {
							ne := &edgesArr[eID]
							/*if e.GetTwoEdgesCycleFlag() > 0 && SumTwoEdgeCycleLen(ne, edgesArr, nodesArr) < TwoEdgeCyleMaxLen {
								if e.StartNID == ne.EndNID && e.EndNID == ne.StartNID { // after process two edge cycle problem
									continue
								}
							}*/
							var nt ExtendPathInfo
							strand := GetNextEdgeStrand2(e, ne, t.GetStrandFlag(), FORWARD)
							//nt.Path = make([]uint32, len(t.Path)+1)
							//copy(nt.Path, t.Path)
							//nt.Path[len(nt.Path)-1] = eID
							nt.ExtendEPI.EID = eID
							nt.ExtendEPI.Last = &t.ExtendEPI
							nt.SetStrandFlag(strand)
							// extend chains
							//if ne.GetTwoEdgesCycleFlag() > 0 && SumTwoEdgeCycleLen(ne, edgesArr, nodesArr) < TwoEdgeCyleMaxLen {
							//nt.MpArr = append(nt.MpArr, t.MpArr...)
							enIdx := IndexEdgeMKINumArr(enArr, uint32(eID))
							if enIdx < 0 {
								nt.MPI = t.MPI
								nt.Score = t.Score
								nt.FlankLen = t.FlankLen + uint32((len(ne.Ks) - (kmerlen - 1)))
								//nt.FlankLen = t.FlankLen
							} else {
								kbIdx := enArr[enIdx].Num
								nt.MPI, nt.Score, nt.FlankLen = ExtendPathChain(t.MPI, nt.ExtendEPI, t.FlankLen, kmerlen, FORWARD, *ne, int(kbIdx), strand, t.Score, edgesArr, nodesArr, uint16(SeedLen), maxPathInfoArr, LongEdgeMinLen, MinChainScoreIdentityPercent)
							} //else if IsSelfCycleEdge(edgesArr, nodesArr,ne.ID) && ne.GetSeqLen() < 500 {

							//fmt.Printf("[FindONTReadsPath]eID: %v, nt.Score: %v\n", eID, nt.Score)
							nt.SetComingFlag(GetEdgeComing(edgesArr[eID], FORWARD, strand, nodesArr, t.GetComingFlag()))
							nt.ExtendLen += t.ExtendLen + uint32(len(ne.Ks)-(kmerlen-1))
							//if len(nt.MpArr) > len(t.MpArr) {
							//nt.FlankLen = GetMappingDBGFlankLen(nt, FORWARD, len(ne.Ks))
							/*} else {
								nt.FlankLen = t.FlankLen + len(ne.Ks) - (kmerlen - 1)
							}*/
							if DebugModel {
								fmt.Printf("[FindONTReadsPath]eID:%v len(nt.MPI.Arr):%v Score:%v ExtendLen:%v FlankLen:%v Strand:%v Coming:%v\n", eID, len(nt.MPI.Arr), nt.Score, nt.ExtendLen, nt.FlankLen, nt.GetStrandFlag(), nt.GetComingFlag())
							}
							// delete low score extend paths
							var sp ScorePos
							var etLen int
							te := edgesArr[nt.ExtendEPI.EID]
							if te.GetSeqLen() > kmerlen+EdgeMapKmerMin {
								etLen = int(nt.ExtendLen)
							}
							mpi1 := &nt.MPI
							mkib := mpi1.Kb[mpi1.Base+uint32(mpi1.Arr[len(mpi1.Arr)-1])]
							sp.Pos = uint32(mkib.GetPosition()) + mkib.GetMapLen()
							sp.Score = uint32(nt.Score)
							if CheckExtend(nt, MinChainScoreIdentityPercent) && BiggerThanMaxScoreArr(maxScoreArr, startPos, sp, etLen, FORWARD, ScoreWinSize) {
								p := int(sp.Pos)
								if nt.FlankLen > MaxEnlongationSeqLen {
									if p > rl-MaxFlankLen {
										highQualEPathInfoArr = append(highQualEPathInfoArr, nt)
									}
								} else {
									if p > rl-kmerlen || IsTipEdge(*ne) {
										highQualEPathInfoArr = append(highQualEPathInfoArr, nt)
									} else {
										stk = append(stk, nt)
										addNum++
									}
								}
								if nt.Score > t.Score {
									maxScoreArr, _ = AddedToMaxScoreArrMPI(maxScoreArr, startPos, nt.MPI, int(nt.Score), FORWARD, ScoreWinSize)
								}
								if DebugModel {
									fmt.Printf("[FindONTReadsPath]sp:%v maxScoreArr:%v\n", sp, maxScoreArr)
								}
							}
						}
						// process bubble edges, just choose best score edge
						if len(eArr) == 2 && addNum == len(eArr) && edgesArr[eArr[0]].GetBubbleFlag() > 0 {
							var deleteNum int
							stk, deleteNum = MuskBubble(eArr, stk)
							addNum -= deleteNum
						} else if len(eArr) >= 2 && addNum == len(eArr) { // check if diff length bubble
							var deleteNum int
							stk, deleteNum = MuskDiffLenBubble(e, eArr, stk, edgesArr, nodesArr, kmerlen, MaxPathLen, MaxBubbleSeqLen)
							addNum -= deleteNum
						}
						if DebugModel {
							fmt.Printf("[FindONTReadsPath]len(stk):%v addNum:%v\n", len(stk), addNum)
						}
						//stk[idx].Score = 0
						//GetNearEdgeIDArr(nd, eID, coming)
						//GetNeighbourEID(eID, nID, edgesArr, nodesArr, max_insert_size)
						//GetNextDirection(eID, edge, nodesArr)
						//GetNextEID(eID, node)
						//GetNextMappingEID(nd, e, base)
					}

					// process highQualEPathInfoArr
					bestEPIF = GetBestEPathInfo(highQualEPathInfoArr, edgesArr, nodesArr, rb.ReadsArr[j], SeedLen, kmerlen, FORWARD, MaxPathLen, MaxBubbleSeqLen, winSize, maxPathInfoArr, maxScoreArr, anchor, ScoreWinSize)
					if len(bestEPIF.Path) == 0 {
						bestEPIF.Path = GetExtendPathInfoPath(bestEPIF.ExtendEPI)
					}

					if len(bestEPIF.Path) > 1 {
						bestEPIF = BubbleAlignment(bestEPIF, edgesArr, nodesArr, rb.ReadsArr[j], SeedLen, kmerlen, FORWARD, TwoEdgeCyleMaxLen, MaxPathLen, MaxBubbleSeqLen, winSize, maxPathInfoArr)
					}
					if DebugModel {
						fmt.Printf("[FindONTReadsPath] bestEPIF path: %v\n", bestEPIF.Path)
					}
				}
			}
			// read BACKWARD mapping
			{
				// check if is boundary of long read
				pos := start
				if pos-extendLen > 0 {
					anchor := mpi.Kb[mpi.Base+uint32(mpi.Arr[0])]
					// choose edge  boudnary anchor will miss some better alignment  because we use minimizers

					//var mpArr []MapingKmerInfo
					//mpArr = append(mpArr, anchor)
					var epi ExtendPathInfo
					id := eID
					epi.ExtendEPI.Last = nil
					epi.ExtendEPI.EID = uint32(eID)
					//epi.Path = append(epi.Path, uint32(id))
					epi.SetStrandFlag(anchor.GetStrand())
					var coming uint8
					epi.SetComingFlag(GetEdgeComing(edgesArr[id], BACKWARD, epi.GetStrandFlag(), nodesArr, coming))
					if anchor.GetStrand() == PLUS {
						epi.ExtendLen = anchor.GetEPosition() + uint32(anchor.GetMapLen())
						epi.FlankLen = anchor.GetEPosition()
					} else {
						epi.ExtendLen = uint32(len(edgesArr[id].Ks)) - anchor.GetEPosition()
						epi.FlankLen = epi.ExtendLen - uint32(anchor.GetMapLen())
					}
					epi.MPI = uniqueMaxPathInfo
					RevUint8Arr(epi.MPI.Arr, int(epi.MPI.LenArr))
					epi.Score = uint32(anchor.GetMapLen())
					//maxMapLen := int(anchor.GetMapLen())
					stk = stk[:0]
					stk = append(stk, epi)
					stkAddNum := 1
					addNum := 0
					highQualEPathInfoArr = highQualEPathInfoArr[:0]
					CleanScoreArr(maxScoreArr)
					maxScoreArr = maxScoreArr[:0]
					startPos := int(anchor.GetPosition()) + int(anchor.GetMapLen())
					for {
						if stkAddNum/StkCleanSize < (stkAddNum+addNum)/StkCleanSize {
							stk = CleanExtendPathInfoArr(stk)
						}
						stkAddNum += addNum
						addNum = 0
						if len(stk) > MaxStackAllow {
							fmt.Fprintf(os.Stderr, "[FindONTReadsPath]readID:%d maxScoreArr:%v\n", ri.ID, maxScoreArr)
							break
						}
						// found max score in the stack
						t, idx := GetMaxScoreEPI2(stk)
						if idx >= 0 {
							stk[idx].Score = 0
						}
						if t.Score <= 0 {
							break
						}
						var sa ScorePos
						//sa.ExtendLen = t.ExtendLen
						var etLen int
						te := edgesArr[t.ExtendEPI.EID]
						if te.GetSeqLen() > kmerlen+EdgeMapKmerMin {
							etLen = int(t.ExtendLen)
						}
						mpi1 := &t.MPI
						mkib := mpi1.Kb[mpi1.Base+uint32(mpi1.Arr[len(mpi1.Arr)-1])]
						sa.Score, sa.Pos = uint32(t.Score), uint32(mkib.GetPosition())
						if BiggerThanMaxScoreArr(maxScoreArr, startPos, sa, etLen, BACKWARD, ScoreWinSize) == false {
							continue
						}

						// check if has NGS Path
						{
							id := t.ExtendEPI.EID
							le := &edgesArr[id]
							if len(le.NGSPathArr) == 1 {
								path := EdgeFreqArr2Uint32Arr(le.NGSPathArr[0])
								index := IndexUint32(path, id)
								var np []uint32
								if t.GetStrandFlag() == PLUS {
									if index >= 2 {
										np = GetReverseUint32Arr(path[:index])
									}
								} else {
									if len(path)-index >= 3 {
										np = path[index+1:]
									}
								}

								if len(np) >= 2 && edgesArr[np[0]].GetTwoEdgesCycleFlag() == 0 && edgesArr[np[len(np)-1]].GetTwoEdgesCycleFlag() == 0 && edgesArr[np[len(np)-1]].StartNID != edgesArr[np[len(np)-1]].EndNID {
									ok := true
									for _, eID := range np {
										var nt ExtendPathInfo
										ne := &edgesArr[eID]
										strand := GetNextEdgeStrand2(le, ne, t.GetStrandFlag(), FORWARD)
										nt.ExtendEPI.EID = eID
										nt.ExtendEPI.Last = &t.ExtendEPI
										nt.SetStrandFlag(strand)
										enIdx := IndexEdgeMKINumArr(enArr, uint32(eID))
										if enIdx < 0 {
											nt.MPI = t.MPI
											nt.Score = t.Score
											nt.FlankLen = t.FlankLen + uint32(len(ne.Ks)-(kmerlen-1))
										} else {
											kbIdx := enArr[enIdx].Num
											nt.MPI, nt.Score, nt.FlankLen = ExtendPathChain(t.MPI, nt.ExtendEPI, t.FlankLen, kmerlen, BACKWARD, *ne, int(kbIdx), strand, t.Score, edgesArr, nodesArr, uint16(SeedLen), maxPathInfoArr, LongEdgeMinLen, MinChainScoreIdentityPercent)
										}
										nt.SetComingFlag(GetEdgeComing(edgesArr[eID], BACKWARD, strand, nodesArr, t.GetComingFlag()))
										nt.ExtendLen += t.ExtendLen + uint32(len(ne.Ks)-(kmerlen-1))
										var sp ScorePos
										//sp.ExtendLen = nt.ExtendLen
										var etLen int
										te := edgesArr[nt.ExtendEPI.EID]
										if te.GetSeqLen() > kmerlen+EdgeMapKmerMin {
											etLen = int(nt.ExtendLen)
										}
										mpi1 := &nt.MPI
										mkb := mpi1.Kb[mpi1.Base+uint32(mpi1.Arr[len(mpi1.Arr)-1])]
										sp.Pos = uint32(mkb.GetPosition())
										sp.Score = uint32(nt.Score)
										if BiggerThanMaxScoreArr(maxScoreArr, startPos, sp, etLen, BACKWARD, ScoreWinSize) {
											if nt.Score > t.Score {
												maxScoreArr, _ = AddedToMaxScoreArrMPI(maxScoreArr, startPos, nt.MPI, int(nt.Score), BACKWARD, ScoreWinSize)
											}
										} else {
											ok = false
											break
										}
										le = ne
										t = nt
									}
									if ok {
										if DebugModel {
											fmt.Printf("[FindONTReadsPath] np: %v\n", np)
										}
										p := int(mkib.GetPosition())
										if t.FlankLen > MaxEnlongationSeqLen {
											if p < MaxFlankLen {
												highQualEPathInfoArr = append(highQualEPathInfoArr, t)
											}
										} else {
											if p < kmerlen || IsTipEdge(edgesArr[t.ExtendEPI.EID]) {
												highQualEPathInfoArr = append(highQualEPathInfoArr, t)
											} else {
												stk = append(stk, t)
												addNum++
											}
										}
									}
									continue
								}
							}
						}

						// add next near edge info
						id := t.ExtendEPI.EID
						e := &edgesArr[id]
						eArr := GetNearDirectionEdgeIDArr(*e, BACKWARD, t.GetStrandFlag(), t.GetComingFlag(), nodesArr)
						if DebugModel {
							fmt.Printf("[FindONTReadsPath]idx:%v nextEIDArr:%v len(t.MPI.Arr):%v Score:%v ExtendLen:%v FlankLen:%v Strand:%v Coming:%v Path:%v\n", idx, eArr, len(t.MPI.Arr), t.Score, t.ExtendLen, t.FlankLen, t.GetStrandFlag(), t.GetComingFlag(), t.Path)
						}
						for _, eID := range eArr {
							ne := &edgesArr[eID]
							var nt ExtendPathInfo
							strand := GetNextEdgeStrand2(e, ne, t.GetStrandFlag(), FORWARD)
							//nt.Path = make([]uint32, len(t.Path)+1)
							//copy(nt.Path, t.Path)
							//nt.Path[len(nt.Path)-1] = eID
							nt.ExtendEPI.EID = eID
							nt.ExtendEPI.Last = &t.ExtendEPI
							nt.SetStrandFlag(strand)
							enIdx := IndexEdgeMKINumArr(enArr, uint32(eID))
							if enIdx < 0 {
								//nt.MpArr = append(nt.MpArr, t.MpArr...)
								nt.MPI = t.MPI
								nt.Score = t.Score
								nt.FlankLen = t.FlankLen + uint32(len(ne.Ks)-(kmerlen-1))
							} else {
								kbIdx := enArr[enIdx].Num
								nt.MPI, nt.Score, nt.FlankLen = ExtendPathChain(t.MPI, nt.ExtendEPI, t.FlankLen, kmerlen, BACKWARD, *ne, int(kbIdx), strand, t.Score, edgesArr, nodesArr, uint16(SeedLen), maxPathInfoArr, LongEdgeMinLen, MinChainScoreIdentityPercent)
							}
							//PrintAddedMpArr(nt.MpArr[len(t.MpArr):])
							nt.SetComingFlag(GetEdgeComing(edgesArr[eID], BACKWARD, strand, nodesArr, t.GetComingFlag()))
							nt.ExtendLen += t.ExtendLen + uint32(len(ne.Ks)-(kmerlen-1))
							/*if len(nt.MpArr) > len(t.MpArr) {
								nt.FlankLen = GetMappingDBGFlankLen(nt, BACKWARD, len(ne.Ks))
							} else {
								nt.FlankLen = t.FlankLen + len(ne.Ks) - (kmerlen - 1)
							}*/
							if DebugModel {
								fmt.Printf("[FindONTReadsPath]eID:%v len(nt.MPI.Arr):%v Score:%v ExtendLen:%v FlankLen:%v Strand:%v Coming:%v Path:%v\n", eID, len(nt.MPI.Arr), nt.Score, nt.ExtendLen, nt.FlankLen, nt.GetStrandFlag(), nt.GetComingFlag(), nt.Path)
							}
							// delete low score extend paths
							//if nt.Score*100/(nt.ExtendLen-nt.FlankLen) < avgSeqIdentityPercent*9/10 {
							var sp ScorePos
							//sp.ExtendLen = nt.ExtendLen
							var etLen int
							te := edgesArr[nt.ExtendEPI.EID]
							if te.GetSeqLen() > kmerlen+EdgeMapKmerMin {
								etLen = int(nt.ExtendLen)
							}
							mpi1 := &nt.MPI
							mkib := mpi1.Kb[mpi1.Base+uint32(mpi1.Arr[len(mpi1.Arr)-1])]
							sp.Pos = uint32(mkib.GetPosition())
							sp.Score = uint32(nt.Score)
							//fmt.Printf("[FindONTReadsPath]sp: %v, maxScoreArr: %v\n", sp, maxScoreArr)
							if CheckExtend(nt, MinChainScoreIdentityPercent) && BiggerThanMaxScoreArr(maxScoreArr, startPos, sp, etLen, BACKWARD, ScoreWinSize) {
								p := int(sp.Pos)
								if nt.FlankLen > MaxEnlongationSeqLen {
									if p < MaxFlankLen {
										highQualEPathInfoArr = append(highQualEPathInfoArr, nt)
									}
								} else {
									if p < kmerlen || IsTipEdge(*ne) {
										highQualEPathInfoArr = append(highQualEPathInfoArr, nt)
									} else {
										stk = append(stk, nt)
										addNum++
									}
								}
								if nt.Score > t.Score {
									maxScoreArr, _ = AddedToMaxScoreArrMPI(maxScoreArr, startPos, nt.MPI, int(nt.Score), BACKWARD, ScoreWinSize)
								}
								if DebugModel {
									fmt.Printf("[FindONTReadsPath]sp:%v maxScoreArr:%v\n", sp, maxScoreArr)
								}
							}
						}
						// process bubble edges, just choose best score edge
						if len(eArr) == 2 && addNum == len(eArr) && edgesArr[eArr[0]].GetBubbleFlag() > 0 {
							var deleteNum int
							stk, deleteNum = MuskBubble(eArr, stk)
							addNum -= deleteNum
						} else if len(eArr) >= 2 && addNum == len(eArr) { // check if diff length bubble
							var deleteNum int
							stk, deleteNum = MuskDiffLenBubble(e, eArr, stk, edgesArr, nodesArr, kmerlen, MaxPathLen, MaxBubbleSeqLen)
							addNum -= deleteNum
						}
						if DebugModel {
							fmt.Printf("[FindONTReadsPath]len(stk):%v addNum:%v\n", len(stk), addNum)
						}
						//stk[idx].Score = 0
					}
					// process highQualEPathInfoArr
					bestEPIB = GetBestEPathInfo(highQualEPathInfoArr, edgesArr, nodesArr, rb.ReadsArr[j], SeedLen, kmerlen, BACKWARD, MaxPathLen, MaxBubbleSeqLen, winSize, maxPathInfoArr, maxScoreArr, anchor, ScoreWinSize)
					if len(bestEPIB.Path) == 0 {
						bestEPIB.Path = GetExtendPathInfoPath(bestEPIB.ExtendEPI)
					}
					if len(bestEPIB.Path) > 1 {
						bestEPIB = BubbleAlignment(bestEPIB, edgesArr, nodesArr, rb.ReadsArr[j], SeedLen, kmerlen, BACKWARD, TwoEdgeCyleMaxLen, MaxPathLen, MaxBubbleSeqLen, winSize, maxPathInfoArr)
					}
					if DebugModel {
						fmt.Printf("[FindONTReadsPath]bestEPIB path:%v\n", bestEPIB.Path)
					}
				}
			}

			// extend path by DBG
			var longRMInfo LongReadMappingInfo
			//bl, fl := len(bestEPIB.Path), len(bestEPIF.Path)
			//bestEPIB = DeleteUnmapPath(bestEPIB, BACKWARD, edgesArr, kmerlen)
			//bestEPIF = DeleteUnmapPath(bestEPIF, FORWARD, edgesArr, kmerlen)
			/*if len(bestEPIB.Path) < bl || len(bestEPIF.Path) < fl {
				fmt.Printf("[FindONTReadsPath]len(bestEPIB.Path): %v < bl: %v, len(bestEPIF.Path):%v < fl: %v\n", len(bestEPIB.Path), bl, len(bestEPIF.Path), fl)
			}*/
			if DebugModel {
				fmt.Printf("[FindONTReadsPath]bestEPIB MapLen:%v Score:%d Path:%v\n", bestEPIB.ExtendLen-bestEPIB.FlankLen, bestEPIB.Score, bestEPIB.Path)
				fmt.Printf("[FindONTReadsPath]bestEPIF MapLen:%v Score:%d Path:%v\n", bestEPIF.ExtendLen-bestEPIF.FlankLen, bestEPIF.Score, bestEPIF.Path)
			}
			longRMInfo = MergeTwoFlankMapingPath(bestEPIB, bestEPIF, mpi, eID)
			longRMInfo.ID = rID
			longRMInfo.Score = int(bestEPIB.Score + uint32(mpi.Score) + bestEPIF.Score)
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
				fmt.Printf("[FindONTReadsPath]longRMInfo RStart:%d REnd:%d Score:%d Path:%v\n", longRMInfo.RStart, longRMInfo.REnd, longRMInfo.Score, longRMInfo.Path)
			}
			//debug.PrintStack()
			//fmt.Fprintf(os.Stderr, "[FindONTReadsPath]debug.PrintStack():%v\n", debug.PrintStack())
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
			pathArr := make([][]uint32, len(kbArr))
			for m, mp := range kbArr {
				kb := mp.Arr
				var edgesSeq []byte
				var chainArr []Chain
				edgesSeq, pathArr[m], chainArr = GetEdgesSeq(edgesArr, nodesArr, kb, kmerlen)
				readSeq := GetReadSeq(rb.ReadsArr[j].Seq, kb[0].GetPosition(), kb[len(kb)-1].GetPosition()+uint32(kb[len(kb)-1].GetMapLen()))
				cgArr[m] = DPLocalAlign(edgesSeq, readSeq, chainArr)
			}
			maxScore := math.MinInt32
			var count int
			maxArr := make([][]uint32, 0, len(cgArr))
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
				var eID1, eID2 uint32
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

func GetHighRegionLen(maxScoreRegionArr []MaxPathInfo, priorIdxArr []uint16, SeedLen int) (highRegLen int) {
	lastP := 0
	for _, idx := range priorIdxArr {
		mpi := maxScoreRegionArr[idx]
		mka, mkb := mpi.Kb[mpi.Base+uint32(mpi.Arr[0])], mpi.Kb[mpi.Base+uint32(mpi.Arr[len(mpi.Arr)-1])]
		a, b := int(mka.GetPosition()), int(mkb.GetPosition())+int(mkb.GetMapLen())
		if a >= lastP-SeedLen {
			highRegLen += b - a
		} else {
			highRegLen += (b - lastP)
		}
		lastP = b
	}
	return
}

func GetIndexUint16Arr(arr []uint16, c uint16) int {
	idx := -1
	for i, id := range arr {
		if id == c {
			idx = i
			break
		}
	}
	return idx
}

/*
func FindDBGPath(highQMPIArr []MaxPathInfo, priorIdxArr []uint16, anchorIdx int, edgesArr []DBGEdge, nodesArr []DBGNode, kmerlen int, ri ReadInfo) {
	anchorMPI := highQMPIArr[anchorIdx]
	idx := GetIndexUint16Arr(priorIdxArr, anchorIdx)
	mka, mkb := anchorMPI.Kb[anchorMPI.Base+uint32(anchorMPI.Arr[0])], anchorMPI.Kb[anchorMPI.Base+uint32(anchorMPI.Arr[len(anchorMPI.Arr)-1])]
	start, end := int(mka.GetPosition()), int(mkb.GetPosition())+int(mkb.GetMapLen())
	e := edgesArr[anchorMPI.EID]
	tipBoundNum := 0
	path := make([]uint32, 0, 5)
	// BACKWARD
	mpi := anchorMPI
	mpiIdx := anchorIdx
	strand := mka.GetStrand()
	var coming uint8
	path  = append(path, e.ID)
	for start > kmerlen*2 {
		if strand == PLUS {
			if e.StartNID == 0 {
				tipBoundNum++
				break
			}
			nd := nodesArr[e.StartNID]
			coming = GetEdgeComing(e, BACKWARD, strand, nodesArr, coming)
			eArr := GetNearDirectionEdgeIDArr(e, BACKWARD, strand, coming, nodesArr)
			eIDArr := GetShareEID(eArr, highQMPIArr,mpiIdx, priorIdxArr, idx, BACKWARD, edgesArr)

		} else {
			if e.EndNID == 0 {
				tipBoundNum++
				break
			}
		}
	}

	// FORWARD

}*/

func IsInUint16Arr(arr []uint16, m uint16) bool {
	for _, ele := range arr {
		if ele == m {
			return true
		}
	}

	return false
}

func DeletePriorIdxArrPath(priorIdxArr []uint16, pathIdxArr []uint16) []uint16 {
	count := 0
	for _, idx := range priorIdxArr {
		if !IsInUint16Arr(pathIdxArr, idx) {
			priorIdxArr[count] = idx
			count++
		}
	}
	priorIdxArr = priorIdxArr[:count]
	return priorIdxArr
}

func GetRegionSeedInfoArr(seedInfoArr []SeedInfo, start, end, readLen uint64) []SeedInfo {
	idx := len(seedInfoArr) * int(start) / int(readLen)
	if idx >= len(seedInfoArr) {
		idx = len(seedInfoArr) - 1
	}
	var si, ei int
	if seedInfoArr[idx].GetPosition() > start {
		for ; idx >= 0; idx-- {
			if seedInfoArr[idx].GetPosition() < start {
				break
			}
		}
		idx++
		si = idx
		for ; idx < len(seedInfoArr); idx++ {
			if seedInfoArr[idx].GetPosition()+SeedLen > end {
				break
			}
		}
		ei = idx - 1
	} else {
		for ; idx < len(seedInfoArr); idx++ {
			if seedInfoArr[idx].GetPosition() >= start {
				break
			}
		}
		si = idx
		for ; idx < len(seedInfoArr); idx++ {
			if seedInfoArr[idx].GetPosition()+SeedLen > end {
				break
			}
		}
		ei = idx - 1
	}
	if si > ei {
		return nil
	}

	sIA := make([]SeedInfo, ei+1-si)
	copy(sIA, seedInfoArr[si:ei+1])
	return sIA
}

func GetEdgeInterSectionKmerInfo(edgeSeedInfoArr, readSeedInfoArr []SeedInfo, buf []MapingKmerInfo, strand bool, eID uint32) []MapingKmerInfo {
	buf = buf[:0]
	i, j := 0, 0
	rsi := readSeedInfoArr[j]
	for ; i < len(edgeSeedInfoArr) && j < len(readSeedInfoArr); i++ {
		esi := edgeSeedInfoArr[i]
		if esi.GetKmer() < rsi.GetKmer() {
			continue
		} else if esi.GetKmer() > rsi.GetKmer() {
			for j = j + 1; j < len(readSeedInfoArr); j++ {
				rsi = readSeedInfoArr[j]
				if esi.GetKmer() <= rsi.GetKmer() {
					break
				}
			}
		}
		if esi.GetKmer() < rsi.GetKmer() {
			continue
		} else if esi.GetKmer() > rsi.GetKmer() {
			break
		}
		var mk MapingKmerInfo
		mk.EID = eID
		mk.SetEPosition(uint32(esi.GetPosition()))
		mk.SetEStrand(esi.GetStrand())
		mk.SetMapLen(uint32(SeedLen))
		for x := j; x < len(readSeedInfoArr); x++ {
			tsi := readSeedInfoArr[x]
			if tsi.GetKmer() > esi.GetKmer() {
				break
			}
			if strand {
				if esi.GetStrand() != tsi.GetStrand() {
					continue
				}
			} else {
				if esi.GetStrand() == tsi.GetStrand() {
					continue
				}
			}
			mk.SeedInfo = SeedInfo(tsi)
			buf = append(buf, mk)
		}
	}
	return buf
}

func GetEdgeChainBlocks(buf, ka []MapingKmerInfo) []MapingKmerInfo {
	sort.Sort(MapingKmerInfoEPosArr(buf))
	ka = ka[:0]
	for i := 0; i < len(buf); i++ {
		mk := buf[i]
		if mk.GetMapLen() == 0 {
			continue
		}
		mE := mk.GetEPosition()
		last := mk
		lR := uint32(last.GetPosition())
		lE := mE
		strand := last.GetMapStrand()
		j := i + 1
		for ; j < len(buf); j++ {
			nk := buf[j]
			nR := uint32(nk.GetPosition())
			nE := nk.GetEPosition()
			if nE > lE+last.GetMapLen() {
				break
			}

			if nk.GetMapStrand() != strand {
				continue
			}
			//fmt.Printf("[GetChainBlocks] %v\n", PrintMKI(nk))
			if strand == PLUS {
				if nR-lR == nE-lE {
					mk.SetMapLen(nE + nk.GetMapLen() - mE)
					buf[j].SetMapLen(0)
					last = nk
					lR, lE = nR, nE
				}
			} else {
				if nE-lE == (lR+last.GetMapLen())-(nR+nk.GetMapLen()) {
					//fmt.Printf("[GetChainBlocks] EID: %v, mk RPos: %v, EPos: %v,Len: %v, nk RPos: %v, EPos: %v, Len:%v\n", mk.EID, mk.GetPosition(), mk.GetEPosition(), mk.Len, nk.GetPosition(), nk.GetEPosition(), nk.GetMapLen())
					mk.SetPosition(uint64(nR))
					mk.SetMapLen(nE + nk.GetMapLen() - mE)
					//mk.SetEPosition(nE)
					//mk.GetMapLen() += uint16(((nk.GetPosition() + uint32(nk.GetMapLen())) - (mk.GetPosition() + uint32(mk.GetMapLen()))))
					buf[j].SetMapLen(0)
					last = nk
					lR, lE = nR, nE
					//fmt.Printf("[GetChainBlocks] mk RPos: %v, EPos: %v,Len: %v\n", mk.GetPosition(), mk.GetEPosition(), mk.GetMapLen())
				}
			}
		}
		ka = append(ka, mk)
	}

	return ka
}

func GetEdgeMaxPathInfo(emi EdgeMapInfo, readSeedInfoArr []SeedInfo, direction uint8, readLen, startPos int, edgesArr []DBGEdge, buf []MapingKmerInfo, MinChainScoreIdentityPercent float32) *MaxPathInfo {
	eID := emi.EdgeID
	e := &edgesArr[eID]
	el := e.GetSeqLen()
	var start, end int
	if direction == FORWARD {
		start = startPos
		end = start + el*6/5
	} else {
		end = startPos
		start = end - el*6/5
	}
	if start < 0 {
		start = 0
	}
	if end >= readLen {
		end = readLen - 1
	}
	sIA := GetRegionSeedInfoArr(readSeedInfoArr, uint64(start), uint64(end), uint64(readLen))
	sort.Sort(SeedInfoArr(sIA))
	fmt.Printf("[GetEdgeMaxPathInfo]start:%d end:%d len(sIA):%d len(e.SeedInfoArr):%d\n", start, end, len(sIA), len(e.SeedInfoArr))
	buf = GetEdgeInterSectionKmerInfo(e.SeedInfoArr, sIA, buf, emi.Strand, eID)
	fmt.Printf("[GetEdgeMaxPathInfo]len(buf):%d\n", len(buf))
	if len(buf) < len(e.SeedInfoArr)/6 {
		/*fmt.Printf("[GetEdgeMaxPathInfo]e.SeedInfoArr: ")
		for _, si := range e.SeedInfoArr {
			fmt.Printf("%d  ", si.GetKmer())
		}
		fmt.Printf("\n[GetEdgeMaxPathInfo]sIA: ")
		for _, si := range sIA {
			fmt.Printf("%d  ", si.GetKmer())
		}*/
		//fmt.Printf("\n[GetEdgeMaxPathInfo]len(buf):%d\n", len(buf))
		return nil
	}
	var ka []MapingKmerInfo
	if cap(buf) > 2*len(buf) {
		ka = buf[len(buf) : len(buf)*2]
	} else {
		ka = make([]MapingKmerInfo, len(buf))
	}
	ka = GetEdgeChainBlocks(buf, ka)
	PrintMKIArr(ka)
	// ka need sort edge position
	sc, a, b := GetSumMappingKIArrLen(ka, emi.Strand, el)
	fmt.Printf("eID:%d el:%d edgeMinNum:%d len(ka):%d mapping region:%d sc:%d \n", eID, el, edgesArr[eID].CovD, len(ka), b-a, sc)
	if float32(sc) < float32(el)*MinChainScoreIdentityPercent || b-a < el/2 {
		//var mpi MaxPathInfo
		return nil
	}
	bsA := make([]uint8, len(ka))
	bsA = InitStrandInfo(ka, bsA, emi.Strand)
	mpi := GetMaxScoreArr(ka, emi.Strand, SeedLen, 1, 100, 1, el, bsA, true)
	if mpi.LenArr == 0 {
		return nil
	}

	mpi.Kb = ka
	mpi.Base = 0
	start, end, _, startE, endE := GetMPIInfo(&mpi)
	readRegLen := end - start
	edgeRegLen := endE - startE
	min := MinInt(readRegLen, edgeRegLen)
	if float32(mpi.Score) > float32(el)*MinChainScoreIdentityPercent && AbsInt(readRegLen-edgeRegLen) < min/10 && min > el/2 {
		mpi.EID = eID
	} else {
		return nil
	}

	mpi.Kb = make([]MapingKmerInfo, len(ka))
	copy(mpi.Kb, ka)
	PrintMPI(&mpi, edgesArr)
	//PrintMKIArr(mpi.Kb)

	return &mpi
}

func GetMergedEdgesMaxPathInfo(emiArr []EdgeMapInfo, readSeedInfoArr []SeedInfo, direction uint8, readLen, startPos int, edgesArr []DBGEdge, nodesArr []DBGNode, buf []MapingKmerInfo, MinChainScoreIdentityPercent float32, SeedLen, WinSize, kmerlen int, logfp io.Writer) MaxPathInfo {
	var mpi MaxPathInfo
	eID := emiArr[0].EdgeID
	var e DBGEdge
	if len(emiArr) == 1 {
		e = edgesArr[eID]
		if emiArr[0].Strand == MINUS {
			e.Ks = GetReverseCompByteArr(e.Ks)
		}
	} else {
		path := Transform2EdgeFreq(emiArr)
		var nID uint32
		if direction == FORWARD {
			if emiArr[0].Strand {
				nID = edgesArr[emiArr[0].EdgeID].EndNID
			} else {
				nID = edgesArr[emiArr[0].EdgeID].StartNID
			}
		} else {
			if emiArr[0].Strand {
				nID = edgesArr[emiArr[0].EdgeID].StartNID
			} else {
				nID = edgesArr[emiArr[0].EdgeID].EndNID
			}
		}

		e1, e2 := &edgesArr[emiArr[0].EdgeID], &edgesArr[emiArr[1].EdgeID]
		if e1.StartNID == e1.EndNID || nID != GetShareDBGNID(e1, e2) || nID == 0 {
			return mpi
		} else if e1.GetTwoEdgesCycleFlag() > 0 || e2.GetTwoEdgesCycleFlag() > 0 {
			return mpi
		}
		fmt.Fprintf(logfp, "[GetMergedEdgesMaxPathInfo]nID:%d emiArr:%v\n", nID, emiArr)
		e = ConcatEdgePathUtg(path, edgesArr, nodesArr, kmerlen)
		if e.GetSeqLen() == 0 {
			return mpi
		}
		if direction == BACKWARD {
			ReverseCompByteArr(e.Ks)
		}
	}

	el := e.GetSeqLen()
	var start, end int
	if direction == FORWARD {
		start = startPos
		end = start + el*6/5
	} else {
		end = startPos
		start = end - el*6/5
	}
	if start < 0 {
		start = 0
	}
	if end >= readLen {
		end = readLen - 1
	}

	bufSize := WinSize * 10
	kaBuf := make([]KI, bufSize)
	ka := make([]KI, 50)
	//rs := make([]byte, 500)
	ka = GetSeqMiniKmerInfoArr(e.Ks, kaBuf, ka, uint32(e.ID), SeedLen, WinSize)
	sort.Sort(KIArr(ka))
	e.SeedInfoArr = make([]SeedInfo, len(ka))
	for x, ki := range ka {
		e.SeedInfoArr[x] = SeedInfo(ki.SeedInfo)
	}

	sIA := GetRegionSeedInfoArr(readSeedInfoArr, uint64(start), uint64(end), uint64(readLen))
	sort.Sort(SeedInfoArr(sIA))
	fmt.Fprintf(logfp, "[GetMergedEdgesMaxPathInfo]start:%d end:%d len(sIA):%d len(e.SeedInfoArr):%d\n", start, end, len(sIA), len(e.SeedInfoArr))
	buf = GetEdgeInterSectionKmerInfo(e.SeedInfoArr, sIA, buf, PLUS, eID)
	fmt.Fprintf(logfp, "[GetMergedEdgesMaxPathInfo]len(buf):%d\n", len(buf))
	if len(buf) < len(e.SeedInfoArr)/6 {
		/*fmt.Printf("[GetEdgeMaxPathInfo]e.SeedInfoArr: ")
		for _, si := range e.SeedInfoArr {
			fmt.Printf("%d  ", si.GetKmer())
		}
		fmt.Printf("\n[GetEdgeMaxPathInfo]sIA: ")
		for _, si := range sIA {
			fmt.Printf("%d  ", si.GetKmer())
		}*/
		//fmt.Printf("\n[GetEdgeMaxPathInfo]len(buf):%d\n", len(buf))
		return mpi
	}
	var kb []MapingKmerInfo
	if cap(buf) > 2*len(buf) {
		kb = buf[len(buf) : len(buf)*2]
	} else {
		kb = make([]MapingKmerInfo, len(buf))
	}
	//PrintMKIArr(buf)
	kb = GetEdgeChainBlocks(buf, kb)
	//PrintMKIArr(kb)
	//PrintMKI(kb[0])
	//PrintMKI(kb[len(kb)-1])
	// ka need sort edge position
	sc, a, b := GetSumMappingKIArrLen(kb, PLUS, el)
	fmt.Fprintf(logfp, "eID:%d el:%d edgeMinNum:%d len(kb):%d mapping region:%d sc:%d\n", eID, el, edgesArr[eID].CovD, len(kb), b-a, sc)
	if float32(sc) < float32(el)*MinChainScoreIdentityPercent || b-a < el/2 {
		//var mpi MaxPathInfo
		return mpi
	}
	bsA := make([]uint8, len(kb))
	bsA = InitStrandInfo(kb, bsA, PLUS)
	mpi = GetMaxScoreArr(kb, PLUS, SeedLen, 1, 100, 1, el, bsA, true)
	if mpi.LenArr == 0 {
		var nm MaxPathInfo
		return nm
	}

	mpi.Kb = kb
	mpi.Base = 0
	start, end, _, startE, endE := GetMPIInfo(&mpi)
	readRegLen := end - start
	edgeRegLen := endE - startE
	min := MinInt(readRegLen, edgeRegLen)
	if float32(mpi.Score) > float32(el)*MinChainScoreIdentityPercent && AbsInt(readRegLen-edgeRegLen) < min/10 && min > el/2 {
		mpi.EID = eID
	} else {
		var nm MaxPathInfo
		return nm
	}

	mpi.Kb = make([]MapingKmerInfo, len(kb))
	copy(mpi.Kb, kb)
	if DebugModel {
		fmt.Fprintf(logfp, "[GetMergedEdgesMaxPathInfo]%s\n", PrintMPI(&mpi, edgesArr))
	}

	//PrintMKIArr(mpi.Kb)

	return mpi
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

func GetDistance(startE, endE, el int, strand bool, direction uint8) uint16 {
	var distance int
	if direction == FORWARD {
		if strand {
			distance = el - endE
		} else {
			distance = startE
		}
	} else {
		if strand {
			distance = startE
		} else {
			distance = el - endE
		}
	}
	//fmt.Printf("[GetDistance]startE:%d endE:%d el:%d strand:%t direction:%d distance:%d\n", startE, endE, el, strand, direction, distance)
	return uint16(distance)
}

func GetBoundaryMPIDistance(maxPathInfoArr []MaxPathInfo, direction uint8, mpiIdx uint16, distance uint16, kmerlen int, emiArr []EdgeMapInfo) (uint16, uint16) {
	if direction == FORWARD {
		for _, emi := range emiArr {
			//e := edgesArr[emi.EdgeID]
			//el := e.GetSeqLen()
			el := int(emi.EdgeLen)
			if emi.MPI != nil {
				mpiIdx = emi.MPIIdx
				_, _, _, startE, endE := GetMPIInfo(emi.MPI)
				distance = GetDistance(startE, endE, el, emi.Strand, direction)
				continue
			}
			var ok bool
			_, end, _, _, _ := GetMPIInfo(&maxPathInfoArr[mpiIdx])
			endPos := end + int(distance)*11/10 + (el - (kmerlen - 1))
			for i := int(mpiIdx) + 1; i < len(maxPathInfoArr); i++ {
				mpi := &maxPathInfoArr[i]
				start, _, _, startE, endE := GetMPIInfo(mpi)
				if start >= endPos {
					break
				}
				if mpi.EID == emi.EdgeID {
					mpiIdx = uint16(i)
					distance = GetDistance(startE, endE, el, emi.Strand, direction)
					ok = true
					break
				}
			}
			if !ok {
				distance += uint16(el - (kmerlen - 1))
			}
		}
	} else { //BACKWARD
		for _, emi := range emiArr {
			//e := edgesArr[emi.EdgeID]
			el := int(emi.EdgeLen)
			if emi.MPI != nil {
				mpiIdx = emi.MPIIdx
				_, _, _, startE, endE := GetMPIInfo(emi.MPI)
				distance = GetDistance(startE, endE, el, emi.Strand, direction)
				continue
			}
			var ok bool
			start, _, _, _, _ := GetMPIInfo(&maxPathInfoArr[mpiIdx])
			startPos := start - int(distance)*11/10 - (el - (kmerlen - 1))
			for i := int(mpiIdx) - 1; i >= 0; i-- {
				mpi := &maxPathInfoArr[i]
				start, _, _, startE, endE := GetMPIInfo(mpi)
				if start < startPos {
					break
				}
				if mpi.EID == emi.EdgeID {
					mpiIdx = uint16(i)
					distance = GetDistance(startE, endE, el, emi.Strand, direction)
					ok = true
					break
				}
			}
			if !ok {
				distance += uint16(el - (kmerlen - 1))
			}
		}
	}

	return mpiIdx, distance
}

func FindEdgeMPI(emi EdgeMapInfo, maxPathInfoArr []MaxPathInfo, direction uint8, mpiIdx, distance uint16, kmerlen int) (*MaxPathInfo, int) {
	if direction == FORWARD {
		el := int(emi.EdgeLen)
		_, end, _, _, _ := GetMPIInfo(&maxPathInfoArr[mpiIdx])
		endPos := end + int(distance)*11/10 + (el - (kmerlen - 1))
		for i := int(mpiIdx) + 1; i < len(maxPathInfoArr); i++ {
			mpi := &maxPathInfoArr[i]
			start, _, _, _, _ := GetMPIInfo(mpi)
			if start >= endPos {
				break
			}
			if mpi.EID == emi.EdgeID {
				return &maxPathInfoArr[i], i
			}
		}
	} else {
		el := int(emi.EdgeLen)
		start, _, _, _, _ := GetMPIInfo(&maxPathInfoArr[mpiIdx])
		startPos := start - int(distance)*11/10 - (el - (kmerlen - 1))
		for i := int(mpiIdx) - 1; i >= 0; i-- {
			mpi := &maxPathInfoArr[i]
			start, _, _, _, _ := GetMPIInfo(mpi)
			if start < startPos {
				break
			}
			if mpi.EID == emi.EdgeID {
				return &maxPathInfoArr[i], i
			}
		}
	}

	return nil, -1
}

func IsEuqlEdgeFreqArr(path1, path2 []EdgeFreq, rev bool) bool {
	if len(path1) != len(path2) {
		return false
	}
	if len(path1) == 0 {
		return true
	}
	if rev {
		for i, ef := range path1 {
			if ef.ID != path2[len(path2)-1-i].ID {
				return false
			}
		}
	} else {
		for i, ef := range path1 {
			if ef.ID != path2[i].ID {
				return false
			}
		}
	}
	return true
}

func Transform2EdgeFreq(emiPath []EdgeMapInfo) []EdgeFreq {
	path := make([]EdgeFreq, len(emiPath))
	for i, emi := range emiPath {
		path[i].ID = uint32(emi.EdgeID)
	}
	return path
}

func Transform2EdgeMapInfo(path []EdgeFreq, edgesArr []DBGEdge, e *DBGEdge, strand bool, direction uint8) []EdgeMapInfo {
	emiArr := make([]EdgeMapInfo, len(path))
	for i, ef := range path {
		emiArr[i].EdgeID = uint32(ef.ID)
		ne := &edgesArr[ef.ID]
		emiArr[i].EdgeLen = uint32(ne.GetSeqLen())
		emiArr[i].Strand = GetNextEdgeStrand2(e, ne, strand, direction)
		e = ne
	}
	return emiArr
}

func GetNextNGSPath(NGSPathMat [][]EdgeFreq, emiPath []EdgeMapInfo, e *DBGEdge, edgesArr []DBGEdge, strand bool, direction uint8) [][]EdgeMapInfo {
	path := Transform2EdgeFreq(emiPath)
	var emiMat [][]EdgeMapInfo
	for _, p := range NGSPathMat {
		idx := IndexEdgeFreq(p, e.ID)
		if idx == len(p)-1 {
			p = GetReverseEdgeFreqArr2(p)
			idx = 0
		}
		id := p[idx+1].ID
		f1 := IsInComing(e.EdgeIDIncoming, id)
		f2 := IsInComing(e.EdgeIDOutcoming, id)
		if !f1 && !f2 {
			log.Fatalf("[GetNextNGSPath]id:%d f1:%t f2:%t Incoming:%v Outcoming:%v\n", id, f1, f2, e.EdgeIDIncoming, e.EdgeIDOutcoming)
		}

		//fmt.Printf("[GetNextNGSPath]id:%d f1:%t f2:%t Incoming:%v Outcoming:%v\n", id, f1, f2, e.EdgeIDIncoming, e.EdgeIDOutcoming)
		if (direction == FORWARD && strand) || (direction == BACKWARD && !strand) {
			if f2 {
				min := MinInt(idx+1, len(path))
				if len(p)-(idx+1) > 0 && IsEuqlEdgeFreqArr(p[idx+1-min:idx+1], path[len(path)-min:], false) {
					emiArr := make([]EdgeMapInfo, len(p)-(idx+1))
					for i, ef := range p[idx+1:] {
						var emi EdgeMapInfo
						emi.EdgeID = uint32(ef.ID)
						//fmt.Printf("[GetNextNGSPath]added emi:%v\n", emi)
						emiArr[i] = emi
					}
					emiMat = append(emiMat, emiArr)
				}
			} else {
				min := MinInt(len(p)-idx, len(path))
				if idx > 0 && IsEuqlEdgeFreqArr(p[idx:idx+min], path[len(path)-min:], true) {
					emiArr := make([]EdgeMapInfo, idx)
					for i := idx - 1; i >= 0; i-- {
						ef := p[i]
						var emi EdgeMapInfo
						emi.EdgeID = uint32(ef.ID)
						//fmt.Printf("[GetNextNGSPath]added emi:%v\n", emi)
						emiArr[idx-1-i] = emi
					}
					emiMat = append(emiMat, emiArr)
				}
			}
		} else {
			if f1 {
				min := MinInt(idx+1, len(path))
				if len(p)-(idx+1) > 0 && IsEuqlEdgeFreqArr(p[idx+1-min:idx+1], path[len(path)-min:], false) {
					emiArr := make([]EdgeMapInfo, len(p)-(idx+1))
					for i, ef := range p[idx+1:] {
						var emi EdgeMapInfo
						emi.EdgeID = uint32(ef.ID)
						//fmt.Printf("[GetNextNGSPath]added emi:%v\n", emi)
						emiArr[i] = emi
					}
					emiMat = append(emiMat, emiArr)
				}
			} else {
				min := MinInt(len(p)-idx, len(path))
				//fmt.Printf("[GetNextNGSPath]idx:%d min:%d p:%v path:%v\n", idx, min, p, path)
				if idx > 0 && IsEuqlEdgeFreqArr(p[idx:idx+min], path[len(path)-min:], true) {
					emiArr := make([]EdgeMapInfo, idx)
					for i := idx - 1; i >= 0; i-- {
						ef := p[i]
						var emi EdgeMapInfo
						emi.EdgeID = uint32(ef.ID)
						//fmt.Printf("[GetNextNGSPath]added emi:%v\n", emi)
						emiArr[idx-1-i] = emi
					}
					emiMat = append(emiMat, emiArr)
				}
			}
		}
	}

	if len(emiMat) == 0 {
		return emiMat
	}
	// set emi info
	for i := range emiMat {
		lastE := e
		//fmt.Printf("[GetNextNGSPath]emiMat[%d]:%v\n", i, emiMat[i])
		for j, emi := range emiMat[i] {
			ne := &edgesArr[emi.EdgeID]
			strand = GetNextEdgeStrand2(lastE, ne, strand, direction)
			//fmt.Printf("[GetNextNGSPath]id:%d strand:%t Incoming:%v Outcoming:%v\n", emi.EdgeID, strand, ne.EdgeIDIncoming, ne.EdgeIDOutcoming)
			emi.Strand = strand
			emi.EdgeLen = uint32(ne.GetSeqLen())
			emiMat[i][j] = emi
			lastE = ne
		}
	}

	return emiMat
}

func GetContainedMPIEMIArr(emiArr []EdgeMapInfo) (ea []EdgeMapInfo) {
	for _, emi := range emiArr {
		if emi.MPI != nil {
			ea = append(ea, emi)
		}
	}
	return
}

func GetNextEMIArr(edgesArr []DBGEdge, lastEMI EdgeMapInfo, direction uint8, maxPathInfoArr []MaxPathInfo, mpiIdx, distance uint16, kmerlen int, logfp io.Writer) (emiArr []EdgeMapInfo) {
	e := &edgesArr[lastEMI.EdgeID]
	strand := lastEMI.Strand
	var ecoming [BaseTypeNum]uint32
	if direction == FORWARD {
		if strand {
			ecoming = e.EdgeIDOutcoming
		} else {
			ecoming = e.EdgeIDIncoming
		}
	} else {
		if strand {
			ecoming = e.EdgeIDIncoming
		} else {
			ecoming = e.EdgeIDOutcoming
		}
	}
	if DebugModel {
		fmt.Fprintf(logfp, "[GetNextEMIArr]ecoming:%v\n", ecoming)
	}
	var emi EdgeMapInfo
	for _, id := range ecoming {
		if id < 2 {
			continue
		}
		ne := &edgesArr[id]
		emi.EdgeID = uint32(id)
		emi.EdgeLen = uint32(ne.GetSeqLen())
		emi.Strand = GetNextEdgeStrand2(e, ne, strand, direction)
		if IsEdgeNoDB(ne, kmerlen) {
			emi.MPIIdx = mpiIdx
			emi.Distance = distance + uint16(emi.EdgeLen) - uint16(kmerlen-1)
		} else {
			var idx int
			emi.MPI, idx = FindEdgeMPI(emi, maxPathInfoArr, direction, mpiIdx, distance, kmerlen)
			if ne.CovD >= 5 && emi.MPI == nil {
				if IsRepeatEdge(ne) || ne.GetBubbleFlag() > 0 {
					continue
				}
			}
			if emi.MPI != nil {
				emi.MPIIdx = uint16(idx)
				_, _, strand, startE, endE := GetMPIInfo(emi.MPI)
				if strand != emi.Strand {
					fmt.Fprintf(logfp, "[GetNextEMIArr]error!!! strand:%t != emi.Strand:%t\n", strand, emi.Strand)
					emi.MPI = nil
					emi.MPIIdx = mpiIdx
					emi.Distance = distance + uint16(emi.EdgeLen) - uint16(kmerlen-1)
				} else {
					emi.Distance = uint16(GetDistance(startE, endE, int(emi.EdgeLen), strand, direction))
				}
			} else {
				emi.MPIIdx = mpiIdx
				emi.Distance = distance + uint16(emi.EdgeLen) - uint16(kmerlen-1)
			}
		}
		emiArr = append(emiArr, emi)
	}
	return
}

func GetEMIArrPathLen(arr []EdgeMapInfo, kmerlen int) (pl int) {
	for _, emi := range arr {
		pl += int(emi.EdgeLen) - (kmerlen - 1)
	}
	if pl > 0 {
		pl += (kmerlen - 1)
	}
	return
}

func GetNextEMIMat(edgesArr []DBGEdge, lastEMI EdgeMapInfo, direction uint8, MaxPathLen uint32, kmerlen int, maxPathInfoArr []MaxPathInfo, mpiIdx, distance uint16, logfp io.Writer) (emiMat [][]EdgeMapInfo) {
	var emiArrHeap [][]EdgeMapInfo
	ea := GetNextEMIArr(edgesArr, lastEMI, direction, maxPathInfoArr, mpiIdx, distance, kmerlen, logfp)
	for _, emi := range ea {
		if edgesArr[emi.EdgeID].StartNID == edgesArr[emi.EdgeID].EndNID {
			return
		}
		var arr []EdgeMapInfo
		arr = append(arr, emi)
		emiArrHeap = append(emiArrHeap, arr)
	}

	for len(emiArrHeap) > 0 {
		if len(emiArrHeap) >= 10 || len(emiMat) > 20 {
			emiMat = nil
			break
		}
		arr := emiArrHeap[len(emiArrHeap)-1]
		emiArrHeap = emiArrHeap[:len(emiArrHeap)-1]
		//fmt.Fprintf(logfp, "[GetNextEMIMat]len(emiArrHeap):%d arr:%v\n", len(emiArrHeap), arr)

		pl := uint32(GetEMIArrPathLen(arr, kmerlen))
		lm := arr[len(arr)-1]
		ea = GetNextEMIArr(edgesArr, lm, direction, maxPathInfoArr, lm.MPIIdx, lm.Distance, kmerlen, logfp)
		//fmt.Fprintf(logfp, "[GetNextEMIMat]ea:%v\n", ea)
		for _, emi := range ea {
			var ta []EdgeMapInfo
			ta = append(ta, arr...)
			ta = append(ta, emi)
			CheckCyclePath(ta, kmerlen)
			if pl+emi.EdgeLen-uint32(kmerlen-1) >= MaxPathLen {
				emiMat = append(emiMat, ta)
			} else {
				emiArrHeap = append(emiArrHeap, ta)
			}
		}
	}
	return
}

func DeleteEdgeMapInfoArr(mpiArr []MaxPathInfo, emiArr []EdgeMapInfo) []EdgeMapInfo {
	idx := 0
	eaIdx := 0
	for i := range mpiArr {
		eID := mpiArr[i].EID
		for ; eaIdx < len(emiArr); eaIdx++ {
			if emiArr[eaIdx].EdgeID == eID {
				emiArr[idx] = emiArr[eaIdx]
				idx++
				break
			}
		}
	}
	emiArr = emiArr[:idx]
	return emiArr
}

func GetAllPathArr(maxPathInfoArr []MaxPathInfo, direction uint8, mpiIdx, distance int, nea []EdgeMapInfo, edgesArr []DBGEdge) (emiMat [][]EdgeMapInfo) {

	return
}

func GetEdgeType(e *DBGEdge) string {
	if e.GetUniqueFlag() > 0 {
		return "Unique"
	} else if e.GetBubbleFlag() > 0 {
		return "Bubble"
	} else if e.GetBubbleRepeatFlag() > 0 {
		return "BubbleRepeat"
	} else if e.GetTwoEdgesCycleFlag() > 0 {
		return "TwoEdgesCycle"
	} else if e.StartNID == e.EndNID {
		return "SelfCycle"
	}
	return "Repeat"
}

type ScoreInfo struct {
	Score      float32
	Start, End int
	DistanceE  int
}

func GetRevDirection(direction uint8) uint8 {
	if direction == FORWARD {
		direction = BACKWARD
	} else {
		direction = FORWARD
	}
	return direction
}

func GetMaxScorePath(emiMat [][]EdgeMapInfo, edgesArr []DBGEdge, direction uint8, kmerlen int) int {
	idx := -1
	sia := make([]ScoreInfo, len(emiMat))
	min := len(emiMat[0])
	for i := range emiMat {
		if len(emiMat[i]) < min {
			min = len(emiMat[i])
		}
	}
	for j := 0; j < min; j++ {
		var noDBOk bool
		for i, arr := range emiMat {
			e := &edgesArr[arr[j].EdgeID]
			el := e.GetSeqLen()
			if i == 0 {
				noDBOk = IsEdgeNoDB(e, kmerlen)
			} else {
				if IsEdgeNoDB(e, kmerlen) != noDBOk {
					return idx
				}
			}
			if !noDBOk && arr[j].MPI != nil {
				start, end, strand, startE, endE := GetMPIInfo(arr[j].MPI)
				if sia[i].Score == 0 {
					sia[i].Score = arr[j].MPI.Score
					sia[i].Start = start
					sia[i].End = end
					sia[i].DistanceE = int(GetDistance(startE, endE, el, strand, direction))
				} else {
					if direction == FORWARD {
						if start < sia[i].End && end > sia[i].End {
							sia[i].Score += GetMPIRegionScore(arr[j].MPI, uint32(sia[i].End), uint32(end))
						} else {
							sia[i].Score += arr[j].MPI.Score
							distE := sia[i].DistanceE - (kmerlen - 1) + int(GetDistance(startE, endE, el, strand, GetRevDirection(direction)))
							distR := start - sia[i].End
							min := MinInt(distE, distR)
							gapLen := AbsInt(distE - distR)
							sia[i].Score -= (0.1*float32(min) + float32(gapLen))
						}
						sia[i].End = end
						sia[i].DistanceE = int(GetDistance(startE, endE, el, strand, direction))
					} else {
						if end > sia[i].Start && start < sia[i].Start {
							sia[i].Score += GetMPIRegionScore(arr[j].MPI, uint32(start), uint32(sia[i].Start))
						} else {
							sia[i].Score += arr[j].MPI.Score
							distE := sia[i].DistanceE - (kmerlen - 1) + int(GetDistance(startE, endE, el, strand, GetRevDirection(direction)))
							distR := sia[i].Start - end
							min := MinInt(distE, distR)
							gapLen := AbsInt(distE - distR)
							sia[i].Score -= (0.1*float32(min) + float32(gapLen))
						}
						sia[i].Start = start
						sia[i].DistanceE = int(GetDistance(startE, endE, el, strand, direction))
					}
				}
			} else {
				sia[i].DistanceE += el - (kmerlen - 1)
			}
		}
	}

	maxIdx := 0
	maxScore := sia[0].Score
	maxCount := 1
	for i, si := range sia {
		if si.Score > maxScore {
			maxScore = si.Score
			maxIdx = i
			maxCount = 1
		} else if si.Score == maxScore {
			maxCount++
		}
	}
	if maxScore > 0 && maxCount == 1 {
		idx = maxIdx
	}
	return idx
}

func IsInEMIMat(pm [][]EdgeMapInfo, eID uint32) bool {
	var ok bool
	for _, arr := range pm {
		if arr[0].EdgeID == eID {
			ok = true
			break
		}
	}
	return ok
}

func IsInEMIArr(arr []EdgeMapInfo, eID uint32) (emi EdgeMapInfo, ok bool) {
	for _, ele := range arr {
		if ele.EdgeID == eID {
			ok = true
			emi = ele
			break
		}
	}
	return
}

func IndexMaxPathInfoArr(pm [][]EdgeMapInfo, mpi *MaxPathInfo) (idx int) {
	idx = -1
	for i, arr := range pm {
		if arr[0].EdgeID == mpi.EID && arr[0].Distance == uint16(mpi.Score) {
			idx = i
			break
		}
	}
	return idx
}

func CheckNGSEMIMat(pm [][]EdgeMapInfo, emiArr []EdgeMapInfo, edgesArr []DBGEdge, maxPathInfoArr []MaxPathInfo, direction uint8, mpiIdx, distance uint16, kmerlen int) [][]EdgeMapInfo {
	pmIdx := 0
	for _, arr := range pm {
		emi, ok := IsInEMIArr(emiArr, arr[0].EdgeID)
		//fmt.Printf("[CheckNGSEMIMat]i:%d arr:%v ok:%t emi:%v\n", i, arr, ok, emi)
		if !ok {
			continue
		}
		idx1, dist1 := mpiIdx, distance
		arr[0] = emi
		var le *DBGEdge
		var strand bool
		for j, ele := range arr {
			ne := &edgesArr[ele.EdgeID]
			if j == 0 {
				strand = ele.Strand
				if ele.MPI != nil {
					idx1 = ele.MPIIdx
					_, _, _, startE, endE := GetMPIInfo(ele.MPI)
					dist1 = GetDistance(startE, endE, ne.GetSeqLen(), strand, direction)
				}
				le = ne
				continue
			} else {
				strand = GetNextEdgeStrand2(le, ne, strand, direction)
				arr[j].Strand = strand
			}
			if !IsEdgeNoDB(ne, kmerlen) && ne.CovD >= 5 {
				var idx int
				arr[j].MPI, idx = FindEdgeMPI(ele, maxPathInfoArr, direction, idx1, dist1, kmerlen)
				if idx < 0 {
					ok = false
					break
				}
				_, _, strand, startE, endE := GetMPIInfo(arr[j].MPI)
				if strand != arr[j].Strand {
					ok = false
					break
				} else {
					idx1 = uint16(idx)
					dist1 = GetDistance(startE, endE, ne.GetSeqLen(), strand, direction)
				}
			} else {
				dist1 += uint16(ne.GetSeqLen() - (kmerlen - 1))
			}
			arr[j].MPIIdx = idx1
			arr[j].Distance = dist1
			arr[j].EdgeLen = uint32(ne.GetSeqLen())
			le = ne
		}
		if ok {
			pm[pmIdx] = arr
			pmIdx++
		}
	}
	return pm[:pmIdx]
}
func ReversePathMat(pathMat [][][]EdgeMapInfo) [][][]EdgeMapInfo {
	for _, pa := range pathMat {
		for _, arr := range pa {
			al := len(arr)
			for x := 0; x < al/2; x++ {
				arr[x], arr[al-1-x] = arr[al-1-x], arr[x]
			}
		}
	}
	pl := len(pathMat)
	for i := 0; i < pl/2; i++ {
		pathMat[i], pathMat[pl-1-i] = pathMat[pl-1-i], pathMat[i]
	}

	return pathMat
}

func CheckCyclePath(path []EdgeMapInfo, kmerlen int) {
	for i := len(path) - 1; i >= 0; i-- {
		m1 := path[i]
		if m1.MPI == nil {
			continue
		}
		for j := i - 1; j >= 0; j-- {
			m2 := path[j]
			if m2.MPI == nil {
				continue
			}
			if m1.MPIIdx == m2.MPIIdx {
				path[i].MPI = nil
				mpiIdx := path[i-1].MPIIdx
				dist := path[i-1].Distance + uint16(path[i].EdgeLen) - uint16(kmerlen-1)
				path[i].MPIIdx = mpiIdx
				path[i].Distance = dist
				for x := i + 1; x < len(path); x++ {
					if path[x].MPI != nil {
						break
					}
					path[x].MPIIdx = mpiIdx
					dist += uint16(path[x].EdgeLen) - uint16(kmerlen-1)
					path[x].Distance = dist
				}
				break
			}
		}
	}
}

func EqualEdgeMapInfoArr(arr1, arr2 []EdgeMapInfo) bool {
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

func DeleteContainedPath(pm [][]EdgeMapInfo) [][]EdgeMapInfo {
	idx := 0
	for _, arr := range pm {
		var ok bool
		for j := 0; j < idx; j++ {
			min := MinInt(len(arr), len(pm[j]))
			if EqualEdgeMapInfoArr(arr[:min], pm[j][:min]) {
				if len(arr) > len(pm[j]) {
					pm[j] = arr
				}
				ok = true
				break
			}
		}
		if !ok {
			pm[idx] = arr
			idx++
		}
	}

	return pm[:idx]
}

func TransformPathMat2String(pathMat [][][]EdgeMapInfo, maxPathInfoArr []MaxPathInfo, rID int, rl int, anotition string) string {
	//ps = make(string, len(pathMat))
	ps := make([]byte, 0, 80)
	ps = append(ps, []byte(strconv.Itoa(rID))...)
	ps = append(ps, '\t')
	ps = append(ps, []byte(strconv.Itoa(rl))...)
	ps = append(ps, '\t')
	emi := pathMat[0][0][0]
	start, _, _, _, _ := GetMPIInfo(&maxPathInfoArr[emi.MPIIdx])
	path := pathMat[len(pathMat)-1][0]
	idx := path[len(path)-1].MPIIdx
	_, end, _, _, _ := GetMPIInfo(&maxPathInfoArr[idx])
	ps = append(ps, []byte(strconv.Itoa(start))...)
	ps = append(ps, '\t')
	ps = append(ps, []byte(strconv.Itoa(end))...)
	ps = append(ps, '\t')
	if emi.Strand == PLUS {
		ps = append(ps, '+')
	} else {
		ps = append(ps, '-')
	}
	ps = append(ps, '\t')

	// add path
	for _, pa := range pathMat {
		if len(pa) == 1 {
			path := pa[0]
			for _, emi := range path {
				ps = append(ps, []byte(strconv.Itoa(int(emi.EdgeID)))...)
				ps = append(ps, '>')
			}
		} else {
			for _, path := range pa {
				for _, emi := range path {
					ps = append(ps, []byte(strconv.Itoa(int(emi.EdgeID)))...)
					ps = append(ps, '-')
				}
				ps[len(ps)-1] = '|'
			}
			ps[len(ps)-1] = '>'
		}
	}
	ps = append(ps, '\t')
	ps = append(ps, []byte(anotition)...)
	return string(ps)
}

func FindMappingPath(readSeedInfoArr []SeedInfo, maxPathInfoArr []MaxPathInfo, anchorIdx int, priorIdxArr []uint16, edgesArr []DBGEdge, nodesArr []DBGNode, kmerlen, readLen int, buf []MapingKmerInfo, MinChainScoreIdentityPercent float32, WinSize int, logfp io.Writer) [][][][]EdgeMapInfo {
	var mapArr [][][][]EdgeMapInfo
	MaxPathLen := uint32(400)

	for anchorIdx >= 0 {

		var pathMat [][][]EdgeMapInfo
		// sequence read BACKWARD
		{
			mpi := &maxPathInfoArr[anchorIdx]
			_, _, strand, startE, endE := GetMPIInfo(mpi)
			var path []EdgeMapInfo
			var emi EdgeMapInfo
			emi.EdgeID = mpi.EID
			emi.Strand = strand
			emi.EdgeLen = uint32(edgesArr[emi.EdgeID].GetSeqLen())
			emi.MPIIdx = uint16(anchorIdx)
			distance := GetDistance(startE, endE, edgesArr[mpi.EID].GetSeqLen(), strand, BACKWARD)
			emi.Distance = distance
			emi.MPI = mpi
			path = append(path, emi)
			mpiIdx := uint16(anchorIdx)
			//pathMat = append(pathMat, [][]EdgeMapInfo{path})
			for {
				lastEMI := path[len(path)-1]
				CheckCyclePath(path, kmerlen)
				e := &edgesArr[lastEMI.EdgeID]
				start, _, _, _, _ := GetMPIInfo(&maxPathInfoArr[mpiIdx])
				if start-int(distance) <= 0 {
					fmt.Fprintf(logfp, "[FindMappingPath]start-distance:%d <= 0\n", start-int(distance))
					break
				} else if distance > 500 {
					fmt.Fprintf(logfp, "[FindMappingPath]distance:%d >= 500\n", distance)
					break
				}

				if DebugModel {
					max := MaxInt(len(path)-2, 0)
					fmt.Fprintf(logfp, "[FindMappingPath]len(path):%d path:%v mpiIdx:%d distance:%d edgeType:%s\n", len(path), path[max:], mpiIdx, distance, GetEdgeType(e))
				}

				strand = lastEMI.Strand
				emiArr := GetNextEMIArr(edgesArr, lastEMI, BACKWARD, maxPathInfoArr, mpiIdx, distance, kmerlen, logfp)
				if len(emiArr) == 0 {
					fmt.Fprintf(logfp, "[FindMappingPath]Next emiArr len:%d\n", len(emiArr))
					//fmt.Printf("[FindMappingPath]Incoming:%v Outcoming:%v startNd:%d endNd:%d\n", e.EdgeIDIncoming, e.EdgeIDOutcoming, nodesArr[e.StartNID], nodesArr[e.EndNID])
					break
				}
				if DebugModel {
					fmt.Fprintf(logfp, "[FindMappingPath]Next emiArr:%v\n", emiArr)
				}
				//fmt.Printf("[FindMappingPath]NGSPathArr:%v\n", npa)
				var pm [][]EdgeMapInfo
				npa := e.NGSPathArr
				if len(npa) > 0 {
					pm = GetNextNGSPath(npa, path, e, edgesArr, strand, BACKWARD)
					pm = DeleteContainedPath(pm)
					// check pm by mpi
					pm = CheckNGSEMIMat(pm, emiArr, edgesArr, maxPathInfoArr, BACKWARD, mpiIdx, distance, kmerlen)
					if DebugModel {
						fmt.Fprintf(logfp, "[FindMappingPath]NGSPath:%v\n\tpm:%v\n", npa, pm)
					}
					if len(emiArr) == 1 && len(pm) == 1 {
						//ea := Transform2EdgeMapInfo(pm[0], edgesArr, e, strand)
						path = append(path, pm[0]...)
						mpiIdx, distance = GetBoundaryMPIDistance(maxPathInfoArr, BACKWARD, mpiIdx, distance, kmerlen, pm[0])
						continue
					}
				}

				if len(emiArr) == 1 {
					path = append(path, emiArr[0])
					mpiIdx, distance = GetBoundaryMPIDistance(maxPathInfoArr, BACKWARD, mpiIdx, distance, kmerlen, emiArr)
					continue
				}

				// len(emiArr) >= 2
				if start-int(distance) < kmerlen {
					fmt.Fprintf(logfp, "[FindMappingPath]start-distance:%d < kmerlen\n", start-int(distance))
					break
				}

				// found all possible path
				if len(pm) > 0 {
					for _, emi := range emiArr {
						if !IsInEMIMat(pm, emi.EdgeID) {
							var arr []EdgeMapInfo
							arr = append(arr, emi)
							//ne := &edgesArr[emi.EdgeID]
							ea := GetNextEMIArr(edgesArr, emi, BACKWARD, maxPathInfoArr, mpiIdx, distance, kmerlen, logfp)
							if len(ea) == 1 {
								arr = append(arr, ea[0])
							}
							pm = append(pm, arr)
						}
					}
				} else {
					pm = GetNextEMIMat(edgesArr, lastEMI, BACKWARD, MaxPathLen, kmerlen, maxPathInfoArr, mpiIdx, distance, logfp)
				}

				if DebugModel {
					for i, arr := range pm {
						fmt.Fprintf(logfp, "[FindMappingPath]pm[%d]:%v\n", i, arr)
					}
				}

				if len(pm) == 0 {
					fmt.Fprintf(logfp, "[FindMappingPath]len(pm) == 0\n")
					break
				} else if len(pm) == 1 {
					path = append(path, pm[0]...)
					mpiIdx, distance = GetBoundaryMPIDistance(maxPathInfoArr, BACKWARD, mpiIdx, distance, kmerlen, pm[0])
					continue
				} else if len(pm) == 2 {
					e1, e2 := &edgesArr[pm[0][0].EdgeID], &edgesArr[pm[1][0].EdgeID]
					if IsBubble(e1, e2, nodesArr) {
						if len(path) > 0 {
							var pa [][]EdgeMapInfo
							pa = append(pa, path)
							pathMat = append(pathMat, pa)
							var np []EdgeMapInfo
							path = np
							ta := GetNextEMIArr(edgesArr, pm[0][0], BACKWARD, maxPathInfoArr, mpiIdx, distance, kmerlen, logfp)
							if len(ta) == 1 {
								path = append(path, ta[0])
							} else {
								fmt.Fprintf(logfp, "[FindMappingPath]ta:%v\n", ta)
								break
							}
						}
						var p1, p2 []EdgeMapInfo
						p1 = append(p1, pm[0][0])
						p2 = append(p2, pm[1][0])
						var pa [][]EdgeMapInfo
						pa = append(pa, p1)
						pa = append(pa, p2)
						pathMat = append(pathMat, pa)
						if DebugModel {
							fmt.Fprintf(logfp, "[FindMappingPath]added bubble pathArr:%v\n", pa)
						}
						mpiIdx, distance = GetBoundaryMPIDistance(maxPathInfoArr, BACKWARD, mpiIdx, distance, kmerlen, []EdgeMapInfo{pm[0][0], path[0]})
						continue
					}
				}

				// check  self cycle path
				var cycle bool
				for _, arr := range pm {
					eID := arr[0].EdgeID
					if edgesArr[eID].StartNID == edgesArr[eID].EndNID {
						cycle = true
						break
					}
				}
				if cycle {
					break
				}

				// get mpi
				{
					mpiArr := make([]MaxPathInfo, len(pm))
					rStart, _, _, _, _ := GetMPIInfo(&maxPathInfoArr[mpiIdx])
					for i, emiArr := range pm {
						mpiArr[i] = GetMergedEdgesMaxPathInfo(emiArr, readSeedInfoArr, BACKWARD, readLen, rStart-int(distance)*11/10+kmerlen, edgesArr, nodesArr, buf, MinChainScoreIdentityPercent, SeedLen, WinSize, kmerlen, logfp)
						pm[i][0].Distance = uint16(mpiArr[i].Score)
					}
					if DebugModel {
						for i := range mpiArr {
							fmt.Fprintf(logfp, "[FindMappingPath]idx:%d len(pm):%d el:%d mpi:%s\n", i, len(pm), GetEMIArrPathLen(pm[i], kmerlen), PrintMPI(&mpiArr[i], edgesArr))
						}
					}
					mpiArr = FilterMPIArr(mpiArr, kmerlen, readLen, edgesArr, logfp)
					if len(mpiArr) == 1 {
						idx := IndexMaxPathInfoArr(pm, &mpiArr[0])
						path = append(path, pm[idx]...)
						mpiIdx, distance = GetBoundaryMPIDistance(maxPathInfoArr, BACKWARD, mpiIdx, distance, kmerlen, pm[idx])
						continue
					} else {
						fmt.Fprintf(logfp, "[FindMappingPath]len(mpiArr):%d\n", len(mpiArr))
						break
					}
				}

				//idx := GetMaxScorePath(pm, edgesArr, FORWARD, kmerlen)
			}

			if len(path) > 0 {
				pathMat = append(pathMat, [][]EdgeMapInfo{path})
			}
			pathMat = ReversePathMat(pathMat)
			fmt.Fprintf(logfp, "[FindMappingPath]BACKWARD pathMat:%v\n", pathMat)
		}

		// sequence read FORWARD
		{
			/*mpi := &maxPathInfoArr[anchorIdx]
			var path []EdgeMapInfo
			var emi EdgeMapInfo
			emi.EdgeID = mpi.EID
			emi.Strand = strand
			emi.EdgeLen = uint32(edgesArr[emi.EdgeID].GetSeqLen())
			emi.MPI = mpi
			path = append(path, emi) */
			path := pathMat[len(pathMat)-1][0]
			pathMat = pathMat[:len(pathMat)-1]
			fmt.Fprintf(logfp, "[FindMappingPath]FORWARD path:%v\n", path)
			emi := path[len(path)-1]
			_, _, strand, startE, endE := GetMPIInfo(emi.MPI)
			mpiIdx := uint16(anchorIdx)
			distance := GetDistance(startE, endE, edgesArr[emi.EdgeID].GetSeqLen(), strand, FORWARD)
			for {
				lastEMI := path[len(path)-1]
				CheckCyclePath(path, kmerlen)
				e := &edgesArr[lastEMI.EdgeID]
				_, end, _, _, _ := GetMPIInfo(&maxPathInfoArr[mpiIdx])
				//pathMat = append(pathMat, [][]EdgeMapInfo{path})
				if end+int(distance) >= readLen {
					fmt.Fprintf(logfp, "[FindMappingPath]end+distance:%d >= readLen\n", end+int(distance))
					break
				} else if distance > 500 {
					fmt.Fprintf(logfp, "[FindMappingPath]distance:%d >= 500\n", distance)
					break
				}

				if DebugModel {
					max := MaxInt(len(path)-2, 0)
					fmt.Fprintf(logfp, "[FindMappingPath]len(path):%d path:%v mpiIdx:%d distance:%d edgeType:%s\n", len(path), path[max:], mpiIdx, distance, GetEdgeType(e))
				}

				strand = lastEMI.Strand
				emiArr := GetNextEMIArr(edgesArr, lastEMI, FORWARD, maxPathInfoArr, mpiIdx, distance, kmerlen, logfp)
				if len(emiArr) == 0 {
					fmt.Fprintf(logfp, "[FindMappingPath]Next emiArr len:%d\n", len(emiArr))
					//fmt.Printf("[FindMappingPath]Incoming:%v Outcoming:%v startNd:%d endNd:%d\n", e.EdgeIDIncoming, e.EdgeIDOutcoming, nodesArr[e.StartNID], nodesArr[e.EndNID])
					break
				}
				if DebugModel {
					fmt.Fprintf(logfp, "[FindMappingPath]Next emiArr:%v\n", emiArr)
				}
				//fmt.Printf("[FindMappingPath]NGSPathArr:%v\n", npa)
				var pm [][]EdgeMapInfo
				npa := e.NGSPathArr
				if len(npa) > 0 {
					pm = GetNextNGSPath(npa, path, e, edgesArr, strand, FORWARD)
					pm = DeleteContainedPath(pm)
					// check pm by mpi
					pm = CheckNGSEMIMat(pm, emiArr, edgesArr, maxPathInfoArr, FORWARD, mpiIdx, distance, kmerlen)
					if DebugModel {
						fmt.Fprintf(logfp, "[FindMappingPath]NGSPath:%v\n\tpm:%v\n", npa, pm)
					}
					if len(emiArr) == 1 && len(pm) == 1 {
						//ea := Transform2EdgeMapInfo(pm[0], edgesArr, e, strand)
						path = append(path, pm[0]...)
						mpiIdx, distance = GetBoundaryMPIDistance(maxPathInfoArr, FORWARD, mpiIdx, distance, kmerlen, pm[0])
						continue
					}
				}

				if len(emiArr) == 1 {
					path = append(path, emiArr[0])
					mpiIdx, distance = GetBoundaryMPIDistance(maxPathInfoArr, FORWARD, mpiIdx, distance, kmerlen, emiArr)
					continue
				}

				// len(emiArr) >= 2
				if end+int(distance) > readLen-kmerlen {
					fmt.Fprintf(logfp, "[FindMappingPath]end+distance:%d > readLen-kmerlen\n", end+int(distance))
					break
				}

				// found all possible path
				if len(pm) > 0 {
					for _, emi := range emiArr {
						if !IsInEMIMat(pm, emi.EdgeID) {
							var arr []EdgeMapInfo
							arr = append(arr, emi)
							//ne := &edgesArr[emi.EdgeID]
							ea := GetNextEMIArr(edgesArr, emi, FORWARD, maxPathInfoArr, mpiIdx, distance, kmerlen, logfp)
							if len(ea) == 1 {
								arr = append(arr, ea[0])
							}
							pm = append(pm, arr)
						}
					}
				} else {
					pm = GetNextEMIMat(edgesArr, lastEMI, FORWARD, MaxPathLen, kmerlen, maxPathInfoArr, mpiIdx, distance, logfp)
				}

				if DebugModel {
					for i, arr := range pm {
						fmt.Fprintf(logfp, "[FindMappingPath]pm[%d]:%v\n", i, arr)
					}
				}

				if len(pm) == 0 {
					fmt.Fprintf(logfp, "[FindMappingPath]len(pm) == 0\n")
					break
				} else if len(pm) == 1 {
					path = append(path, pm[0]...)
					mpiIdx, distance = GetBoundaryMPIDistance(maxPathInfoArr, FORWARD, mpiIdx, distance, kmerlen, pm[0])
					continue
				} else if len(pm) == 2 {
					e1, e2 := &edgesArr[pm[0][0].EdgeID], &edgesArr[pm[1][0].EdgeID]
					if IsBubble(e1, e2, nodesArr) {
						if len(path) > 0 {
							var pa [][]EdgeMapInfo
							pa = append(pa, path)
							pathMat = append(pathMat, pa)
							var np []EdgeMapInfo
							path = np
							ta := GetNextEMIArr(edgesArr, pm[0][0], FORWARD, maxPathInfoArr, mpiIdx, distance, kmerlen, logfp)
							if len(ta) == 1 {
								path = append(path, ta[0])
							} else {
								fmt.Fprintf(logfp, "[FindMappingPath]ta:%v\n", ta)
								break
							}
						}
						var p1, p2 []EdgeMapInfo
						p1 = append(p1, pm[0][0])
						p2 = append(p2, pm[1][0])
						var pa [][]EdgeMapInfo
						pa = append(pa, p1)
						pa = append(pa, p2)
						pathMat = append(pathMat, pa)
						if DebugModel {
							fmt.Fprintf(logfp, "[FindMappingPath]added bubble pathArr:%v\n", pa)
						}
						mpiIdx, distance = GetBoundaryMPIDistance(maxPathInfoArr, FORWARD, mpiIdx, distance, kmerlen, []EdgeMapInfo{pm[0][0], path[0]})
						continue
					}
				}
				// check  self cycle path
				var cycle bool
				for _, arr := range pm {
					eID := arr[0].EdgeID
					if edgesArr[eID].StartNID == edgesArr[eID].EndNID {
						cycle = true
						break
					}
				}
				if cycle {
					break
				}

				// get mpi
				{
					mpiArr := make([]MaxPathInfo, len(pm))
					_, rEnd, _, _, _ := GetMPIInfo(&maxPathInfoArr[mpiIdx])
					for i, emiArr := range pm {
						mpiArr[i] = GetMergedEdgesMaxPathInfo(emiArr, readSeedInfoArr, FORWARD, readLen, rEnd+int(distance)*9/10-kmerlen, edgesArr, nodesArr, buf, MinChainScoreIdentityPercent, SeedLen, WinSize, kmerlen, logfp)
						pm[i][0].Distance = uint16(mpiArr[i].Score)
					}
					if DebugModel {
						for i := range mpiArr {
							fmt.Fprintf(logfp, "[FindMappingPath]idx:%d len(pm):%d el:%d mpi:%s\n", i, len(pm), GetEMIArrPathLen(pm[i], kmerlen), PrintMPI(&mpiArr[i], edgesArr))
						}
					}
					mpiArr = FilterMPIArr(mpiArr, kmerlen, readLen, edgesArr, logfp)
					if len(mpiArr) == 1 {
						idx := IndexMaxPathInfoArr(pm, &mpiArr[0])
						path = append(path, pm[idx]...)
						mpiIdx, distance = GetBoundaryMPIDistance(maxPathInfoArr, FORWARD, mpiIdx, distance, kmerlen, pm[idx])
						continue
					} else {
						fmt.Fprintf(logfp, "[FindMappingPath]len(mpiArr):%d\n", len(mpiArr))
						break
					}
				}

				//idx := GetMaxScorePath(pm, edgesArr, FORWARD, kmerlen)
			}
		}

		if DebugModel {
			fmt.Fprintf(logfp, "[FindMappingPath]pathMat:%v\n", pathMat)
		}
		// transform pathMat to output string
		//ps := TransformPathMat2String(pathMat)
		if len(pathMat) > 0 {
			mapArr = append(mapArr, pathMat)
		}
		break

		var pathIdxArr []uint16
		// find next anchorIdx
		priorIdxArr = DeletePriorIdxArrPath(priorIdxArr, pathIdxArr)
		anchorIdx = GetAnchorMPIIdx(maxPathInfoArr, priorIdxArr, edgesArr)
	}

	return mapArr
}

func MapingONTReadsChunk(dbgKmerInfoPac DBGKmerInfoPac, cs <-chan ReadsBucketInfo, SeedLen, winSize, GapCost, MaxGapLen, kmerlen int, wc chan<- string, edgesArr []DBGEdge, nodesArr []DBGNode, SeqType uint8, MaxPathLen, MaxBubbleSeqLen int, MinChainScoreIdentityPercent float32, prefix string, processID int) {
	logfn := prefix + "." + strconv.Itoa(processID) + ".log"
	fp, err := os.Create(logfn)
	if err != nil {
		log.Fatalf("[MapingONTReadsChunk] file %s create error, err: %v\n", logfn, err)
	}
	defer fp.Close()
	logfp := fp
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
	bucketInterSecKmerInfoArr := make([][]MapingKmerInfo, BuckReadsNum)
	for i := 0; i < BuckReadsNum; i++ {
		bucketInterSecKmerInfoArr[i] = make([]MapingKmerInfo, 0, 100000)
	}
	//aka := make([]MapingKmerInfo, 0, 50000000)
	//kb := make([]MapingKmerInfo, 0, 300000)
	buf := make([]MapingKmerInfo, 0, 100000)
	//EdgeSeedNumCountArr := make([]EdgeSeedNumInfo, 0, 1200)
	maxPathInfoArr := make([]MaxPathInfo, 0, 1000)
	//highQMPIArr := make([]MaxPathInfo, 0, 3000)
	//stk := make([]ExtendPathInfo, 0, 1000)
	//StkCleanSize := 100
	//maxScoreArr := make([]ScorePos, 0, 1000)
	//ScoreWinSize := 50
	//highQualEPathInfoArr := make([]ExtendPathInfo, 0, 10)
	for {
		rb := <-cs
		if len(rb.ReadsArr) == 0 {
			var ts string
			wc <- ts
			break
		}
		totalReadNum += len(rb.ReadsArr)
		startReadID := rb.ReadsArr[0].ID
		bucketInterSecKmerInfoArr = GetInterSectionKmerInfoB(dbgKmerInfoPac, rb.KmerSortArr, uint32(startReadID), uint16(SeedLen), bucketInterSecKmerInfoArr)
		// loop process every read
		for j, ka := range bucketInterSecKmerInfoArr {
			if len(ka) == 0 {
				continue
			}
			rID := uint32(startReadID) + uint32(j)
			ri := rb.ReadsArr[j]
			rl := len(ri.Seq)
			seedInfoArr := rb.SeedInfoMat[j]
			fmt.Fprintf(logfp, "read[%v] length: %v %s\n", rID, rl, ri.Anotition)
			MQ := 20
			if len(ri.Anotition) > 3 {
				x, err := strconv.Atoi(string(ri.Anotition[3:]))
				if err != nil {
					//MQ = 20
				} else {
					MQ = x
				}
			}

			//aka = aka[:0]
			var seedNumPercent float32
			bkaL := len(ka)
			var esArr []EdgeStrand
			ka, esArr, seedNumPercent = AddHighFreqKmer2(ka, dbgKmerInfoPac, uint16(SeedLen), LongEdgeMinLen, SeedNumCountArr, winSize, kmerlen, MQ)
			//sort.Sort(MapingKmerInfoArr(ka))
			if DebugModel {
				fmt.Fprintf(logfp, "before:%d,  after add high freq kmer:%d\n", bkaL, len(ka))
			}
			if len(ka) < 10 || len(esArr) == 0 {
				fmt.Fprintf(os.Stderr, "[MapingONTReadsChunk]len(ka):%d < 10!!!!\n", len(ka))
				continue
			}
			var enArr []EdgeMKINum
			ka, buf, enArr = GetChainBlocks(ka, buf, esArr)
			if len(enArr) > 1000 {
				fmt.Fprintf(os.Stderr, "[MapingONTReadsChunk]read[%d] len(enArr):%d > 1000\n", ri.ID, len(enArr))
			}
			maxPathInfoArr = GetMaxChainBlockArr(ka, maxPathInfoArr, edgesArr, nodesArr, SeedLen, GapCost, MaxGapLen, kmerlen, SeqType, LongEdgeMinLen, MinChainScoreIdentityPercent, seedNumPercent, logfp)
			fmt.Fprintf(logfp, "len(ka):%d len(enArr):%d len(maxPathInfoArr):%d\n", len(ka), len(enArr), len(maxPathInfoArr))
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
				continue
			} else if len(maxPathInfoArr) == 1 {
				mpi := &maxPathInfoArr[0]
				start, end, _, _, _ := GetMPIInfo(mpi)
				if int(start) < 2*kmerlen && int(end) > rl-2*kmerlen {
					mapOneEdgeNum++
					fmt.Fprintf(logfp, "readID[%d] only map one long edge[%d]\n", rID, mpi.EID)
				}
				continue
			}
			anchorMPI := &maxPathInfoArr[anchorIdx]
			fmt.Fprintf(logfp, "len(mpiArr):%d anchorMPI:%s\n", len(maxPathInfoArr), PrintMPI(anchorMPI, edgesArr))

			// find DBG path
			{
				mapArr := FindMappingPath(seedInfoArr, maxPathInfoArr, anchorIdx, priorIdxArr, edgesArr, nodesArr, kmerlen, rl, buf, MinChainScoreIdentityPercent, winSize, logfp)
				if len(mapArr) > 0 {
					for _, pathMat := range mapArr {
						wc <- TransformPathMat2String(pathMat, maxPathInfoArr, int(ri.ID), rl, string(ri.Anotition))
					}
				}
			}
			/*highRegLen := GetHighRegionLen(maxPathInfoArr, priorIdxArr, SeedLen)
			fmt.Printf("[MapingONTReadsChunk] high region num:%d length:%d percent:%d\n", len(maxPathInfoArr), highRegLen, highRegLen*100/rl)
			if highRegLen >= rl {
				highRegLen = rl - 1
			}*/
			//highRegPercentArr[highRegLen*10/rl]++

			// Find DBG path by long read maping
			//FindDBGPath(highQMPIArr, priorIdxArr, anchorIdx, edgesArr, nodesArr, kmerlen, ri)

			// check this region has been processed
			/*last := mpi.Kb[mpi.Base+uint32(mpi.Arr[len(mpi.Arr)-1])]
					if edge.GetSeqLen() < 3*kmerlen {
						fmt.Printf("[MapingONTReadsChunk] Seed Edge length: %d smaller than %d\n", edge.GetSeqLen(), 3*kmerlen)
					}

					var bestEPIF, bestEPIB ExtendPathInfo
					// read FORWARD mapping
					{
						// check if is boundary of long read
						pos := end
						if pos+extendLen < rl {
							anchor := last
							// choose edge  boudnary anchor will miss some better alignment  because we use minimizers
							//var mpArr []MapingKmerInfo
							//mpArr = append(mpArr, anchor)
							var epi ExtendPathInfo
							id := eID
							epi.ExtendEPI.Last = nil
							epi.ExtendEPI.EID = uint32(eID)
							//epi.Path = append(epi.Path, uint32(id))
							epi.SetStrandFlag(anchor.GetStrand())
							var coming uint8
							epi.SetComingFlag(GetEdgeComing(edgesArr[id], FORWARD, epi.GetStrandFlag(), nodesArr, coming))
							if anchor.GetStrand() == PLUS {
								epi.ExtendLen = uint32(len(edgesArr[id].Ks)) - uint32(anchor.GetEPosition())
								epi.FlankLen = uint32(len(edgesArr[id].Ks) - (int(anchor.GetEPosition()) + int(anchor.GetMapLen())))
							} else {
								epi.ExtendLen = anchor.GetEPosition() + uint32(anchor.GetMapLen())
								epi.FlankLen = anchor.GetEPosition()
							}
							epi.MPI = uniqueMaxPathInfo
							mpiArr := make([]uint16, len(uniqueMaxPathInfo.Arr))
							copy(mpiArr, uniqueMaxPathInfo.Arr)
							epi.MPI.Arr = mpiArr
							epi.Score = uint32(anchor.GetMapLen())
							//maxScore := epi.Score
							//scoreLen := int(anchor.GetMapLen())
							//maxMapLen := int(anchor.GetMapLen())
							stk = stk[:0]
							stk = append(stk, epi)
							stkAddNum := 1
							addNum := 0
							highQualEPathInfoArr = highQualEPathInfoArr[:0]
							CleanScoreArr(maxScoreArr)
							maxScoreArr = maxScoreArr[:0]
							startPos := int(anchor.GetPosition())
							for {
								// found max score in the stack
								if stkAddNum/StkCleanSize < (stkAddNum+addNum)/StkCleanSize {
									stk = CleanExtendPathInfoArr(stk)
								}
								stkAddNum += addNum
								addNum = 0
								t, idx := GetMaxScoreEPI2(stk)
								if idx >= 0 {
									stk[idx].Score = 0
								}
								if t.Score <= 0 {
									break
								}
								if len(stk) > MaxStackAllow {
									fmt.Fprintf(os.Stderr, "[FindONTReadsPath]readID:%d maxScoreArr:%v\n", ri.ID, maxScoreArr)
									break
								}
								var sa ScorePos
								//sa.ExtendLen = t.ExtendLen

								var etLen int
								te := edgesArr[t.ExtendEPI.EID]
								if te.GetSeqLen() > kmerlen+EdgeMapKmerMin {
									etLen = int(t.ExtendLen)
								}
								mpi1 := &t.MPI
								mkib := mpi1.Kb[mpi1.Base+uint32(mpi1.Arr[len(mpi1.Arr)-1])]
								sa.Score, sa.Pos = uint32(t.Score), mkib.GetPosition()+uint32(mkib.GetMapLen())
								if BiggerThanMaxScoreArr(maxScoreArr, startPos, sa, etLen, FORWARD, ScoreWinSize) == false {
									continue
								}
								// check if has NGS Path
								{
									id := t.ExtendEPI.EID
									le := edgesArr[id]
									if len(le.NGSPathArr) == 1 {
										path := le.NGSPathArr[0].IDArr
										index := IndexUint32(path, id)
										var np []uint32
										if t.GetStrandFlag() == PLUS {
											if len(path)-index >= 3 {
												np = path[index+1:]
											}
										} else {
											if index >= 2 {
												np = GetReverseUint32Arr(path[:index])
											}
										}

										if len(np) >= 2 && edgesArr[np[0]].GetTwoEdgesCycleFlag() == 0 && edgesArr[np[len(np)-1]].GetTwoEdgesCycleFlag() == 0 && edgesArr[np[len(np)-1]].StartNID != edgesArr[np[len(np)-1]].EndNID {
											ok := true
											for _, eID := range np {
												var nt ExtendPathInfo
												ne := edgesArr[eID]
												strand := GetNextEdgeStrand(le, ne, nodesArr, t.GetStrandFlag())
												nt.ExtendEPI.EID = eID
												nt.ExtendEPI.Last = &t.ExtendEPI
												//nt.Path = make([]uint32, len(t.Path)+1)
												//copy(nt.Path, t.Path)
												//nt.Path[len(nt.Path)-1] = eID
												nt.SetStrandFlag(strand)
												kbIdx, mok := edgeCountMap[uint32(eID)]
												if !mok {
													nt.MPI = t.MPI
													nt.Score = t.Score
													nt.FlankLen = t.FlankLen + uint32(len(ne.Ks)-(kmerlen-1))
												} else {
													nt.MPI, nt.Score, nt.FlankLen = ExtendPathChain(t.MPI, nt.ExtendEPI, t.FlankLen, kmerlen, FORWARD, ne, int(kbIdx), strand, t.Score, edgesArr, nodesArr, uint16(SeedLen), maxPathInfoArr, LongEdgeMinLen, MinChainScoreIdentityPercent)
												}
												nt.SetComingFlag(GetEdgeComing(edgesArr[eID], FORWARD, strand, nodesArr, t.GetComingFlag()))
												nt.ExtendLen += t.ExtendLen + uint32(len(ne.Ks)-(kmerlen-1))
												var sp ScorePos
												//sa.ExtendLen = nt.ExtendLen
												var etLen int
												te := edgesArr[nt.ExtendEPI.EID]
												if te.GetSeqLen() > kmerlen+EdgeMapKmerMin {
													etLen = int(nt.ExtendLen)
												}
												mpi1 := &nt.MPI
												mkib := mpi1.Kb[mpi1.Base+uint32(mpi1.Arr[len(mpi1.Arr)-1])]
												sp.Pos = mkib.GetPosition() + uint32(mkib.GetMapLen())
												sp.Score = uint32(nt.Score)
												if BiggerThanMaxScoreArr(maxScoreArr, start, sp, etLen, FORWARD, ScoreWinSize) {
													if nt.Score > t.Score {
														maxScoreArr, _ = AddedToMaxScoreArrMPI(maxScoreArr, startPos, nt.MPI, int(nt.Score), FORWARD, ScoreWinSize)
													}
												} else {
													ok = false
													break
												}
												le = ne
												t = nt
											}
											if ok {
												if DebugModel {
													fmt.Printf("[FindONTReadsPath] np: %v\n", np)
												}
												mpi1 := &t.MPI
												mkib := mpi1.Kb[mpi1.Base+uint32(mpi1.Arr[len(mpi1.Arr)-1])]
												p := int(mkib.GetPosition()) + int(mkib.GetMapLen())
												if t.FlankLen > MaxEnlongationSeqLen {
													if p > rl-MaxFlankLen {
														highQualEPathInfoArr = append(highQualEPathInfoArr, t)
													}
												} else {
													if p > rl-kmerlen || IsTipEdge(edgesArr[t.ExtendEPI.EID]) {
														highQualEPathInfoArr = append(highQualEPathInfoArr, t)
													} else {
														stk = append(stk, t)
														addNum++
													}
												}
											}
											continue
										}
									}
								}
								// add next near edge info
								//id := t.Path[len(t.Path)-1]
								id := t.ExtendEPI.EID
								e := edgesArr[id]
								eArr := GetNearDirectionEdgeIDArr(e, FORWARD, t.GetStrandFlag(), t.GetComingFlag(), nodesArr)
								if DebugModel {
									fmt.Printf("[FindONTReadsPath]idx:%v nextEArr:%v len(t.MPI.Arr):%v Score:%v ExtendLen:%v FlankLen:%v Strand:%v Coming:%v LastEID:%d\n", idx, eArr, len(t.MPI.Arr), t.Score, t.ExtendLen, t.FlankLen, t.GetStrandFlag(), t.GetComingFlag(), t.ExtendEPI.EID)
								}
								for _, eID := range eArr {
									ne := edgesArr[eID]
									if e.GetTwoEdgesCycleFlag() > 0 && SumTwoEdgeCycleLen(ne, edgesArr, nodesArr) < TwoEdgeCyleMaxLen {
										if e.StartNID == ne.EndNID && e.EndNID == ne.StartNID { // after process two edge cycle problem
											continue
										}
									}
									var nt ExtendPathInfo
									strand := GetNextEdgeStrand(e, ne, nodesArr, t.GetStrandFlag())
									//nt.Path = make([]uint32, len(t.Path)+1)
									//copy(nt.Path, t.Path)
									//nt.Path[len(nt.Path)-1] = eID
									nt.ExtendEPI.EID = eID
									nt.ExtendEPI.Last = &t.ExtendEPI
									nt.SetStrandFlag(strand)
									// extend chains
									//if ne.GetTwoEdgesCycleFlag() > 0 && SumTwoEdgeCycleLen(ne, edgesArr, nodesArr) < TwoEdgeCyleMaxLen {
									//nt.MpArr = append(nt.MpArr, t.MpArr...)
									kbIdx, mok := edgeCountMap[uint32(eID)]
									if !mok {
										nt.MPI = t.MPI
										nt.Score = t.Score
										nt.FlankLen = t.FlankLen + uint32((len(ne.Ks) - (kmerlen - 1)))
										//nt.FlankLen = t.FlankLen
									} else {
										nt.MPI, nt.Score, nt.FlankLen = ExtendPathChain(t.MPI, nt.ExtendEPI, t.FlankLen, kmerlen, FORWARD, ne, int(kbIdx), strand, t.Score, edgesArr, nodesArr, uint16(SeedLen), maxPathInfoArr, LongEdgeMinLen, MinChainScoreIdentityPercent)
									} //else if IsSelfCycleEdge(edgesArr, nodesArr,ne.ID) && ne.GetSeqLen() < 500 {

									//fmt.Printf("[FindONTReadsPath]eID: %v, nt.Score: %v\n", eID, nt.Score)
									nt.SetComingFlag(GetEdgeComing(edgesArr[eID], FORWARD, strand, nodesArr, t.GetComingFlag()))
									nt.ExtendLen += t.ExtendLen + uint32(len(ne.Ks)-(kmerlen-1))
									//if len(nt.MpArr) > len(t.MpArr) {
									//nt.FlankLen = GetMappingDBGFlankLen(nt, FORWARD, len(ne.Ks))
									} else {
										nt.FlankLen = t.FlankLen + len(ne.Ks) - (kmerlen - 1)
									}
									if DebugModel {
										fmt.Printf("[FindONTReadsPath]eID:%v len(nt.MPI.Arr):%v Score:%v ExtendLen:%v FlankLen:%v Strand:%v Coming:%v\n", eID, len(nt.MPI.Arr), nt.Score, nt.ExtendLen, nt.FlankLen, nt.GetStrandFlag(), nt.GetComingFlag())
									}
									// delete low score extend paths
									var sp ScorePos
									var etLen int
									te := edgesArr[nt.ExtendEPI.EID]
									if te.GetSeqLen() > kmerlen+EdgeMapKmerMin {
										etLen = int(nt.ExtendLen)
									}
									mpi1 := &nt.MPI
									mkib := mpi1.Kb[mpi1.Base+uint32(mpi1.Arr[len(mpi1.Arr)-1])]
									sp.Pos = mkib.GetPosition() + uint32(mkib.GetMapLen())
									sp.Score = uint32(nt.Score)
									if CheckExtend(nt, MinChainScoreIdentityPercent) && BiggerThanMaxScoreArr(maxScoreArr, startPos, sp, etLen, FORWARD, ScoreWinSize) {
										p := int(sp.Pos)
										if nt.FlankLen > MaxEnlongationSeqLen {
											if p > rl-MaxFlankLen {
												highQualEPathInfoArr = append(highQualEPathInfoArr, nt)
											}
										} else {
											if p > rl-kmerlen || IsTipEdge(ne) {
												highQualEPathInfoArr = append(highQualEPathInfoArr, nt)
											} else {
												stk = append(stk, nt)
												addNum++
											}
										}
										if nt.Score > t.Score {
											maxScoreArr, _ = AddedToMaxScoreArrMPI(maxScoreArr, startPos, nt.MPI, int(nt.Score), FORWARD, ScoreWinSize)
										}
										if DebugModel {
											fmt.Printf("[FindONTReadsPath]sp:%v maxScoreArr:%v\n", sp, maxScoreArr)
										}
									}
								}
								// process bubble edges, just choose best score edge
								if len(eArr) == 2 && addNum == len(eArr) && edgesArr[eArr[0]].GetBubbleFlag() > 0 {
									var deleteNum int
									stk, deleteNum = MuskBubble(eArr, stk)
									addNum -= deleteNum
								} else if len(eArr) >= 2 && addNum == len(eArr) { // check if diff length bubble
									var deleteNum int
									stk, deleteNum = MuskDiffLenBubble(e, eArr, stk, edgesArr, nodesArr, kmerlen, MaxPathLen, MaxBubbleSeqLen)
									addNum -= deleteNum
								}
								if DebugModel {
									fmt.Printf("[FindONTReadsPath]len(stk):%v addNum:%v\n", len(stk), addNum)
								}
								//stk[idx].Score = 0
								//GetNearEdgeIDArr(nd, eID, coming)
								//GetNeighbourEID(eID, nID, edgesArr, nodesArr, max_insert_size)
								//GetNextDirection(eID, edge, nodesArr)
								//GetNextEID(eID, node)
								//GetNextMappingEID(nd, e, base)
							}

							// process highQualEPathInfoArr
							bestEPIF = GetBestEPathInfo(highQualEPathInfoArr, edgesArr, nodesArr, rb.ReadsArr[j], SeedLen, kmerlen, FORWARD, MaxPathLen, MaxBubbleSeqLen, winSize, maxPathInfoArr, maxScoreArr, anchor, ScoreWinSize)
							if len(bestEPIF.Path) == 0 {
								bestEPIF.Path = GetExtendPathInfoPath(bestEPIF.ExtendEPI)
							}

							if len(bestEPIF.Path) > 1 {
								bestEPIF = BubbleAlignment(bestEPIF, edgesArr, nodesArr, rb.ReadsArr[j], SeedLen, kmerlen, FORWARD, TwoEdgeCyleMaxLen, MaxPathLen, MaxBubbleSeqLen, winSize, maxPathInfoArr)
							}
							if DebugModel {
								fmt.Printf("[FindONTReadsPath] bestEPIF path: %v\n", bestEPIF.Path)
							}
						}
					}
					// read BACKWARD mapping
					{
						// check if is boundary of long read
						pos := start
						if pos-extendLen > 0 {
							anchor := mpi.Kb[mpi.Base+uint32(mpi.Arr[0])]
							// choose edge  boudnary anchor will miss some better alignment  because we use minimizers

							//var mpArr []MapingKmerInfo
							//mpArr = append(mpArr, anchor)
							var epi ExtendPathInfo
							id := eID
							epi.ExtendEPI.Last = nil
							epi.ExtendEPI.EID = uint32(eID)
							//epi.Path = append(epi.Path, uint32(id))
							epi.SetStrandFlag(anchor.GetStrand())
							var coming uint8
							epi.SetComingFlag(GetEdgeComing(edgesArr[id], BACKWARD, epi.GetStrandFlag(), nodesArr, coming))
							if anchor.GetStrand() == PLUS {
								epi.ExtendLen = anchor.GetEPosition() + uint32(anchor.GetMapLen())
								epi.FlankLen = anchor.GetEPosition()
							} else {
								epi.ExtendLen = uint32(len(edgesArr[id].Ks)) - anchor.GetEPosition()
								epi.FlankLen = epi.ExtendLen - uint32(anchor.GetMapLen())
							}
							epi.MPI = uniqueMaxPathInfo
							epi.MPI.Arr = GetRevUint16Arr(epi.MPI.Arr)
							epi.Score = uint32(anchor.GetMapLen())
							//maxMapLen := int(anchor.GetMapLen())
							stk = stk[:0]
							stk = append(stk, epi)
							stkAddNum := 1
							addNum := 0
							highQualEPathInfoArr = highQualEPathInfoArr[:0]
							CleanScoreArr(maxScoreArr)
							maxScoreArr = maxScoreArr[:0]
							startPos := int(anchor.GetPosition()) + int(anchor.GetMapLen())
							for {
								if stkAddNum/StkCleanSize < (stkAddNum+addNum)/StkCleanSize {
									stk = CleanExtendPathInfoArr(stk)
								}
								stkAddNum += addNum
								addNum = 0
								if len(stk) > MaxStackAllow {
									fmt.Fprintf(os.Stderr, "[FindONTReadsPath]readID:%d maxScoreArr:%v\n", ri.ID, maxScoreArr)
									break
								}
								// found max score in the stack
								t, idx := GetMaxScoreEPI2(stk)
								if idx >= 0 {
									stk[idx].Score = 0
								}
								if t.Score <= 0 {
									break
								}
								var sa ScorePos
								//sa.ExtendLen = t.ExtendLen
								var etLen int
								te := edgesArr[t.ExtendEPI.EID]
								if te.GetSeqLen() > kmerlen+EdgeMapKmerMin {
									etLen = int(t.ExtendLen)
								}
								mpi1 := &t.MPI
								mkib := mpi1.Kb[mpi1.Base+uint32(mpi1.Arr[len(mpi1.Arr)-1])]
								sa.Score, sa.Pos = uint32(t.Score), mkib.GetPosition()
								if BiggerThanMaxScoreArr(maxScoreArr, startPos, sa, etLen, BACKWARD, ScoreWinSize) == false {
									continue
								}

								// check if has NGS Path
								{
									id := t.ExtendEPI.EID
									le := edgesArr[id]
									if len(le.NGSPathArr) == 1 {
										path := le.NGSPathArr[0].IDArr
										index := IndexUint32(path, id)
										var np []uint32
										if t.GetStrandFlag() == PLUS {
											if index >= 2 {
												np = GetReverseUint32Arr(path[:index])
											}
										} else {
											if len(path)-index >= 3 {
												np = path[index+1:]
											}
										}

										if len(np) >= 2 && edgesArr[np[0]].GetTwoEdgesCycleFlag() == 0 && edgesArr[np[len(np)-1]].GetTwoEdgesCycleFlag() == 0 && edgesArr[np[len(np)-1]].StartNID != edgesArr[np[len(np)-1]].EndNID {
											ok := true
											for _, eID := range np {
												var nt ExtendPathInfo
												ne := edgesArr[eID]
												strand := GetNextEdgeStrand(le, ne, nodesArr, t.GetStrandFlag())
												nt.ExtendEPI.EID = eID
												nt.ExtendEPI.Last = &t.ExtendEPI
												nt.SetStrandFlag(strand)
												kbIdx, mok := edgeCountMap[uint32(eID)]
												if !mok {
													nt.MPI = t.MPI
													nt.Score = t.Score
													nt.FlankLen = t.FlankLen + uint32(len(ne.Ks)-(kmerlen-1))
												} else {
													nt.MPI, nt.Score, nt.FlankLen = ExtendPathChain(t.MPI, nt.ExtendEPI, t.FlankLen, kmerlen, BACKWARD, ne, int(kbIdx), strand, t.Score, edgesArr, nodesArr, uint16(SeedLen), maxPathInfoArr, LongEdgeMinLen, MinChainScoreIdentityPercent)
												}
												nt.SetComingFlag(GetEdgeComing(edgesArr[eID], BACKWARD, strand, nodesArr, t.GetComingFlag()))
												nt.ExtendLen += t.ExtendLen + uint32(len(ne.Ks)-(kmerlen-1))
												var sp ScorePos
												//sp.ExtendLen = nt.ExtendLen
												var etLen int
												te := edgesArr[nt.ExtendEPI.EID]
												if te.GetSeqLen() > kmerlen+EdgeMapKmerMin {
													etLen = int(nt.ExtendLen)
												}
												mpi1 := &nt.MPI
												mkb := mpi1.Kb[mpi1.Base+uint32(mpi1.Arr[len(mpi1.Arr)-1])]
												sp.Pos = mkb.GetPosition()
												sp.Score = uint32(nt.Score)
												if BiggerThanMaxScoreArr(maxScoreArr, startPos, sp, etLen, BACKWARD, ScoreWinSize) {
													if nt.Score > t.Score {
														maxScoreArr, _ = AddedToMaxScoreArrMPI(maxScoreArr, startPos, nt.MPI, int(nt.Score), BACKWARD, ScoreWinSize)
													}
												} else {
													ok = false
													break
												}
												le = ne
												t = nt
											}
											if ok {
												if DebugModel {
													fmt.Printf("[FindONTReadsPath] np: %v\n", np)
												}
												p := int(mkib.GetPosition())
												if t.FlankLen > MaxEnlongationSeqLen {
													if p < MaxFlankLen {
														highQualEPathInfoArr = append(highQualEPathInfoArr, t)
													}
												} else {
													if p < kmerlen || IsTipEdge(edgesArr[t.ExtendEPI.EID]) {
														highQualEPathInfoArr = append(highQualEPathInfoArr, t)
													} else {
														stk = append(stk, t)
														addNum++
													}
												}
											}
											continue
										}
									}
								}

								// add next near edge info
								id := t.ExtendEPI.EID
								e := edgesArr[id]
								eArr := GetNearDirectionEdgeIDArr(e, BACKWARD, t.GetStrandFlag(), t.GetComingFlag(), nodesArr)
								if DebugModel {
									fmt.Printf("[FindONTReadsPath]idx:%v nextEIDArr:%v len(t.MPI.Arr):%v Score:%v ExtendLen:%v FlankLen:%v Strand:%v Coming:%v Path:%v\n", idx, eArr, len(t.MPI.Arr), t.Score, t.ExtendLen, t.FlankLen, t.GetStrandFlag(), t.GetComingFlag(), t.Path)
								}
								for _, eID := range eArr {
									ne := edgesArr[eID]
									var nt ExtendPathInfo
									strand := GetNextEdgeStrand(e, ne, nodesArr, t.GetStrandFlag())
									//nt.Path = make([]uint32, len(t.Path)+1)
									//copy(nt.Path, t.Path)
									//nt.Path[len(nt.Path)-1] = eID
									nt.ExtendEPI.EID = eID
									nt.ExtendEPI.Last = &t.ExtendEPI
									nt.SetStrandFlag(strand)
									kbIdx, mok := edgeCountMap[uint32(eID)]
									if !mok {
										//nt.MpArr = append(nt.MpArr, t.MpArr...)
										nt.MPI = t.MPI
										nt.Score = t.Score
										nt.FlankLen = t.FlankLen + uint32(len(ne.Ks)-(kmerlen-1))
									} else {
										nt.MPI, nt.Score, nt.FlankLen = ExtendPathChain(t.MPI, nt.ExtendEPI, t.FlankLen, kmerlen, BACKWARD, ne, int(kbIdx), strand, t.Score, edgesArr, nodesArr, uint16(SeedLen), maxPathInfoArr, LongEdgeMinLen, MinChainScoreIdentityPercent)
									}
									//PrintAddedMpArr(nt.MpArr[len(t.MpArr):])
									nt.SetComingFlag(GetEdgeComing(edgesArr[eID], BACKWARD, strand, nodesArr, t.GetComingFlag()))
									nt.ExtendLen += t.ExtendLen + uint32(len(ne.Ks)-(kmerlen-1))
									if len(nt.MpArr) > len(t.MpArr) {
										nt.FlankLen = GetMappingDBGFlankLen(nt, BACKWARD, len(ne.Ks))
									} else {
										nt.FlankLen = t.FlankLen + len(ne.Ks) - (kmerlen - 1)
									}
									if DebugModel {
										fmt.Printf("[FindONTReadsPath]eID:%v len(nt.MPI.Arr):%v Score:%v ExtendLen:%v FlankLen:%v Strand:%v Coming:%v Path:%v\n", eID, len(nt.MPI.Arr), nt.Score, nt.ExtendLen, nt.FlankLen, nt.GetStrandFlag(), nt.GetComingFlag(), nt.Path)
									}
									// delete low score extend paths
									//if nt.Score*100/(nt.ExtendLen-nt.FlankLen) < avgSeqIdentityPercent*9/10 {
									var sp ScorePos
									//sp.ExtendLen = nt.ExtendLen
									var etLen int
									te := edgesArr[nt.ExtendEPI.EID]
									if te.GetSeqLen() > kmerlen+EdgeMapKmerMin {
										etLen = int(nt.ExtendLen)
									}
									mpi1 := &nt.MPI
									mkib := mpi1.Kb[mpi1.Base+uint32(mpi1.Arr[len(mpi1.Arr)-1])]
									sp.Pos = mkib.GetPosition()
									sp.Score = uint32(nt.Score)
									//fmt.Printf("[FindONTReadsPath]sp: %v, maxScoreArr: %v\n", sp, maxScoreArr)
									if CheckExtend(nt, MinChainScoreIdentityPercent) && BiggerThanMaxScoreArr(maxScoreArr, startPos, sp, etLen, BACKWARD, ScoreWinSize) {
										p := int(sp.Pos)
										if nt.FlankLen > MaxEnlongationSeqLen {
											if p < MaxFlankLen {
												highQualEPathInfoArr = append(highQualEPathInfoArr, nt)
											}
										} else {
											if p < kmerlen || IsTipEdge(ne) {
												highQualEPathInfoArr = append(highQualEPathInfoArr, nt)
											} else {
												stk = append(stk, nt)
												addNum++
											}
										}
										if nt.Score > t.Score {
											maxScoreArr, _ = AddedToMaxScoreArrMPI(maxScoreArr, startPos, nt.MPI, int(nt.Score), BACKWARD, ScoreWinSize)
										}
										if DebugModel {
											fmt.Printf("[FindONTReadsPath]sp:%v maxScoreArr:%v\n", sp, maxScoreArr)
										}
									}
								}
								// process bubble edges, just choose best score edge
								if len(eArr) == 2 && addNum == len(eArr) && edgesArr[eArr[0]].GetBubbleFlag() > 0 {
									var deleteNum int
									stk, deleteNum = MuskBubble(eArr, stk)
									addNum -= deleteNum
								} else if len(eArr) >= 2 && addNum == len(eArr) { // check if diff length bubble
									var deleteNum int
									stk, deleteNum = MuskDiffLenBubble(e, eArr, stk, edgesArr, nodesArr, kmerlen, MaxPathLen, MaxBubbleSeqLen)
									addNum -= deleteNum
								}
								if DebugModel {
									fmt.Printf("[FindONTReadsPath]len(stk):%v addNum:%v\n", len(stk), addNum)
								}
								//stk[idx].Score = 0
							}
							// process highQualEPathInfoArr
							bestEPIB = GetBestEPathInfo(highQualEPathInfoArr, edgesArr, nodesArr, rb.ReadsArr[j], SeedLen, kmerlen, BACKWARD, MaxPathLen, MaxBubbleSeqLen, winSize, maxPathInfoArr, maxScoreArr, anchor, ScoreWinSize)
							if len(bestEPIB.Path) == 0 {
								bestEPIB.Path = GetExtendPathInfoPath(bestEPIB.ExtendEPI)
							}
							if len(bestEPIB.Path) > 1 {
								bestEPIB = BubbleAlignment(bestEPIB, edgesArr, nodesArr, rb.ReadsArr[j], SeedLen, kmerlen, BACKWARD, TwoEdgeCyleMaxLen, MaxPathLen, MaxBubbleSeqLen, winSize, maxPathInfoArr)
							}
							if DebugModel {
								fmt.Printf("[FindONTReadsPath]bestEPIB path:%v\n", bestEPIB.Path)
							}
						}
					}

					// extend path by DBG
					var longRMInfo LongReadMappingInfo
					//bl, fl := len(bestEPIB.Path), len(bestEPIF.Path)
					//bestEPIB = DeleteUnmapPath(bestEPIB, BACKWARD, edgesArr, kmerlen)
					//bestEPIF = DeleteUnmapPath(bestEPIF, FORWARD, edgesArr, kmerlen)
					if len(bestEPIB.Path) < bl || len(bestEPIF.Path) < fl {
						fmt.Printf("[FindONTReadsPath]len(bestEPIB.Path): %v < bl: %v, len(bestEPIF.Path):%v < fl: %v\n", len(bestEPIB.Path), bl, len(bestEPIF.Path), fl)
					}
					if DebugModel {
						fmt.Printf("[FindONTReadsPath]bestEPIB MapLen:%v Score:%d Path:%v\n", bestEPIB.ExtendLen-bestEPIB.FlankLen, bestEPIB.Score, bestEPIB.Path)
						fmt.Printf("[FindONTReadsPath]bestEPIF MapLen:%v Score:%d Path:%v\n", bestEPIF.ExtendLen-bestEPIF.FlankLen, bestEPIF.Score, bestEPIF.Path)
					}
					longRMInfo = MergeTwoFlankMapingPath(bestEPIB, bestEPIF, mpi, eID)
					longRMInfo.ID = rID
					longRMInfo.Score = int(bestEPIB.Score + mpi.Score + bestEPIF.Score)
					longRMInfo.AnchorEdgeID = eID
					longRMInfo.AnchorStrand = last.GetStrand()
					longRMInfo.Qual = uint32(longRMInfo.Score) * 100 / (longRMInfo.REnd - longRMInfo.RStart)
					// add read region to the cBMRIArr
					var cbmr ChainBlockMapLongReadInfo
					cbmr.Start, cbmr.End, cbmr.Score = int(longRMInfo.RStart), int(longRMInfo.REnd), int(longRMInfo.Score)
					cBMRIArr = append(cBMRIArr, cbmr)
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
						fmt.Printf("[FindONTReadsPath]longRMInfo RStart:%d REnd:%d Score:%d Path:%v\n", longRMInfo.RStart, longRMInfo.REnd, longRMInfo.Score, longRMInfo.Path)
					}
				}
			}*/
		}
	}
	fmt.Printf("[MapingONTReadsChunk] total reads num: %v, not found seed edge num: %v, map one edge num: %v, extend read num: %v\n", totalReadNum, notFoundSeedEdgeNum, mapOneEdgeNum, ExtendNum)
	for i, num := range highRegPercentArr {
		fmt.Printf("[MapingONTReadsChunk]highRegPercentArr[%d]:%d\n", i, num)
	}

	for i := 0; i < BuckReadsNum; i++ {
		fmt.Printf("[MapingONTReadsChunk]bucketInterSecKmerInfoArr[%d]cap:%d\n", i, cap(bucketInterSecKmerInfoArr[i]))
	}
	//fmt.Printf("[MapingONTReadsChunk]aka cap:%d\n", cap(aka))
	//fmt.Printf("[MapingONTReadsChunk]ka cap:%d\n", cap(ka))
	fmt.Printf("[MapingONTReadsChunk]buf cap:%d\n", cap(buf))
	//fmt.Printf("[MapingONTReadsChunk]EdgeSeedNumCountArr cap:%d\n", cap())
	fmt.Printf("[MapingONTReadsChunk]maxPathInfoArr cap:%d\n", cap(maxPathInfoArr))

}

func PrintDBGEdges(edgesArr []DBGEdge, prefix string) {
	edgesfn := prefix + ".DBGedges"
	edgesfp, err := os.Create(edgesfn)
	if err != nil {
		log.Fatalf("[PrintTmpDBG] failed to create file: %s, err: %v\n", edgesfn, err)
	}
	defer edgesfp.Close()

	var e *DBGEdge
	for i := range edgesArr {
		e = &edgesArr[i]
		if e.ID < 2 || e.GetDeleteFlag() > 0 {
			continue
		}
		edgesfp.WriteString(fmt.Sprintf("ID:%d len:%d CovD:%d StartNID:%d EndNID:%d Incoming:%v Outcoming:%v\n", e.ID, len(e.Ks), e.CovD, e.StartNID, e.EndNID, e.EdgeIDIncoming, e.EdgeIDOutcoming))
	}

}

func DeconstructDBG(c cli.Command) {
	// check arguments
	opt, suc := checkArgsDDBG(c)
	if suc == false {
		log.Fatalf("[DeconstructDBG] check Arguments error, opt: %v\n", opt)
	}
	DebugModel = opt.Debug
	BuckReadsNum = 100
	EdgeMapKmerMin = 40
	GapCost := 1
	MaxGapLen := 2000
	SeqType := uint8(1)
	MaxPathLen := 10
	MaxBubbleSeqLen := 3 * opt.Kmer
	HighFreqKmerMin := 1000
	HighFreqKmerMax := 100000
	//Kmerlen = opt.Kmer
	fmt.Printf("Arguments: %v\n", c.Flags())

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
	nodesfn := opt.Prefix + ".nodes.MapDBG.Arr.br"
	nodesArr := make([]DBGNode, nSize)
	fc := make(chan int, 1)
	go NodesArrReader(nodesfn, nodesArr, opt.Kmer, fc)
	//if len(nodesArr) != nSize {
	//	log.Fatalf("[DeconstructDBG] len(nodesArr): %v != nodesArr Size: %v in file: %v\n", len(nodesArr), nSize, DBGInfofn)
	//}
	//edgesArr := make([]DBGEdge, eSize)
	edgesfn := opt.Prefix + ".edges.MapDBG.fq.br"
	edgesArr := ReadEdgesFromFile(edgesfn, uint32(eSize))
	<-fc

	go CheckInterConnectivity(edgesArr, nodesArr)
	//PrintTmpDBG(nodesArr, edgesArr, opt.Prefix)

	uniqueNum, bubbleRepeatNum, twoCycleNum, selfCycle, selfCycleSameComingNum, bubbleEdgeNum := SetDBGEdgesUniqueFlag(edgesArr, nodesArr, opt.Kmer)
	AddNodeInfo2DBGEdgeArr(edgesArr, nodesArr)
	runtime.GC()
	fmt.Printf("[DeconstructDBG] unique edge number is: %v, bubbleRepeatNum: %v, twoCycleNum:%v, selfCycle:%v, selfCycleSameComingNum: %v, bubbleEdgeNum: %v\n", uniqueNum, bubbleRepeatNum, twoCycleNum, selfCycle, selfCycleSameComingNum, bubbleEdgeNum)
	t0 := time.Now()
	pathFn := opt.Prefix + ".Path.br"
	mapDBGInfoFn := opt.Prefix + ".ONTMapDBGInfo"
	// no any other align tools for seed edge by ONT reads, use same as daligner,
	// first sort DBG edges by seed kmer, and split ONT reads by bucket, every bucket sort by seed kmer, and found share kmers from edges sort kmersArr
	if _, err := os.Stat(pathFn); err != nil {
		//BucketSize := (1 << 28) // ==2**28
		cfgInfo, err := ParseCfg(opt.CfgFn, false, false)
		if err != nil {
			log.Fatalf("[DeconstructDBG] ParseCfg 'C': %v err :%v\n", opt.CfgFn, err)
		}
		fmt.Printf("[DeconstructDBG] opt: %v\n\tcfgInfo: %v\n", opt, cfgInfo)

		minimerSize := GetCuckoofilterDBGSampleSize(edgesArr, opt.WinSize, opt.Kmer)
		dbgCS := make(chan DBGKmerInfoPac, 1)
		go SortDBGEdgesKmer(dbgCS, edgesArr, SeedLen, opt.WinSize, HighFreqKmerMin, HighFreqKmerMax, int(minimerSize*11/4), opt.Kmer)
		//fmt.Printf("[DeconstructDBG]total found edges minimizers number is: %v\n", len(edgesKmerSortArr))
		//bufSize := int64(minimerSize) / 50
		//bufSize := int64(1000000)
		var processedCPUNum int
		processLongReadCPUNum := opt.NumCPU / 20
		if opt.NumCPU <= 2 {
			processedCPUNum = 1
		} else {
			processedCPUNum = opt.NumCPU - processLongReadCPUNum
		}

		cs := make(chan ReadsBucketInfo, processedCPUNum/4+1)
		wc := make(chan string, processedCPUNum*1000)
		go LoadLongReadsAndSort(cfgInfo, processLongReadCPUNum, cs, SeedLen, opt.WinSize, BuckReadsNum)
		dbgKmerInfoPac := <-dbgCS
		runtime.GC()
		fmt.Printf("[DeconstructDBG] Sort DBGEdge used: %v\n", time.Now().Sub(t0))
		//go PrintDBGEdges(edgesArr, opt.Prefix)
		t0 = time.Now()
		for i := 0; i < processedCPUNum; i++ {
			go MapingONTReadsChunk(dbgKmerInfoPac, cs, SeedLen, opt.WinSize, GapCost, MaxGapLen, opt.Kmer, wc, edgesArr, nodesArr, SeqType, MaxPathLen, MaxBubbleSeqLen, opt.MinChainScoreIdentityPercent, opt.Prefix, i)
		}
		//func FindONTReadsPath(edgesKmerSortArr []KI, cs <-chan ReadsBucketInfo, SeedLen, winSize, GapCost, MaxGapLen, kmerlen int, wc chan<- LongReadMappingInfo, edgesArr []DBGEdge, nodesArr []DBGNode, SeqType uint8) {

		//pathFn := opt.Prefix + ".Path.protobuf"
		WriteLongPathToFile(wc, pathFn, mapDBGInfoFn, processedCPUNum)
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

	t1 := time.Now()
	{
		cs := make(chan MapingInfo, 1000)
		go LoadPathFromFile(pathFn, cs)
		// Simplify using Long Reads Mapping info
		mapRecordSize := LoadONTMapInfoFile(mapDBGInfoFn)
		joinPathArr := SimplifyByLongReadsPath(edgesArr, nodesArr, cs, opt, MaxPathLen, MaxBubbleSeqLen, mapRecordSize)

		//graphfn := opt.Prefix + ".afterLR.dot"
		//GraphvizDBGArr(nodesArr, edgesArr, graphfn)
		DcDBGEdgesfn := opt.Prefix + ".edges.DcDBG.fa"
		ExtractSeq(edgesArr, nodesArr, joinPathArr, DcDBGEdgesfn, opt.Kmer)
		fmt.Printf("[DeconstructDBG] find max edge path and produce sequence used: %v\n", time.Now().Sub(t1))
	}
	fmt.Printf("[DeconstructDBG] total used: %v\n", time.Now().Sub(t0))
	//StoreEdgesToFn(DcDBGEdgesfn, edgesArr)
}
