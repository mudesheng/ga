package deconstructdbg

import (
	"fmt"
	"log"
	"math"
	"os"
	"reflect"
	"sort"

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
		if e.GetUniqueFlag() > 0 {
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

func findMaxPath(sortedEIDIdxArr []IdxLen, edgesArr []constructdbg.DBGEdge, nodesArr []constructdbg.DBGNode, minMapingFreq int) (pathArr []constructdbg.Path) {
	for _, item := range sortedEIDIdxArr {
		e := edgesArr[item.Idx]
		if e.GetDeleteFlag() > 0 || e.GetProcessFlag() > 0 || e.GetUniqueFlag() == 0 || len(e.PathMat) != 1 || len(e.PathMat[0].IDArr) == 0 {
			continue
		}
		fmt.Printf("[findMaxPath] eID: %v, length: %v\n", item.Idx, item.Length)
		maxP := constructdbg.ExtendPath(edgesArr, nodesArr, e, minMapingFreq)
		if len(maxP.IDArr) > 1 {
			pathArr = append(pathArr, maxP)
		}
	}
	return pathArr
}

func constructEdgesPathRelationship(edgesArr []constructdbg.DBGEdge, pathArr [][2][]constructdbg.DBG_MAX_INT, edgesPathRelationArr [][]uint32) {
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

func checkEdgePathArrDirection(edgesArr []constructdbg.DBGEdge, nodesArr []constructdbg.DBGNode, ePathArr [][2][]constructdbg.DBG_MAX_INT, eID constructdbg.DBG_MAX_INT) {
	e := edgesArr[eID]
	var ea1, ea2 []constructdbg.DBG_MAX_INT
	if e.StartNID > 0 {
		ea1 = constructdbg.GetNearEdgeIDArr(nodesArr[e.StartNID], eID)
	}
	if e.EndNID > 0 {
		ea2 = constructdbg.GetNearEdgeIDArr(nodesArr[e.EndNID], eID)
	}

	// found two edge cycle
	if len(ea1) == 1 && len(ea2) == 1 && ea1[0] == ea2[0] {
		edgesArr[eID].ResetUniqueFlag()
		return
	}
	for i, extArr := range ePathArr {
		idx := IndexEID(extArr, eID)
		if idx < len(extArr[0])-1 {
			if constructdbg.IsInDBG_MAX_INTArr(ea1, extArr[0][idx+1]) || constructdbg.IsInDBG_MAX_INTArr(ea1, extArr[1][idx+1]) {
				constructdbg.ReverseDBG_MAX_INTArr(ePathArr[i][0])
				constructdbg.ReverseDBG_MAX_INTArr(ePathArr[i][1])
			}
		} else {
			if constructdbg.IsInDBG_MAX_INTArr(ea2, extArr[0][idx-1]) || constructdbg.IsInDBG_MAX_INTArr(ea2, extArr[1][idx-1]) {
				constructdbg.ReverseDBG_MAX_INTArr(ePathArr[i][0])
				constructdbg.ReverseDBG_MAX_INTArr(ePathArr[i][1])
			}
		}
	}
}

func MergePathArr(edgesArr []constructdbg.DBGEdge, nodesArr []constructdbg.DBGNode, pathArr [][2][]constructdbg.DBG_MAX_INT, edgesPathRelationArr [][]uint32, minMapFreq int) {
	for i, e := range edgesArr {
		if i < 2 || e.GetDeleteFlag() > 0 || len(edgesPathRelationArr[i]) < minMapFreq || e.GetUniqueFlag() == 0 {
			continue
		}

		// check if is a node cycle edge
		if e.StartNID == e.EndNID {
			node := nodesArr[e.StartNID]
			n1, _ := constructdbg.GetEdgeIDComing(node.EdgeIDIncoming)
			n2, _ := constructdbg.GetEdgeIDComing(node.EdgeIDOutcoming)
			if n1 > 1 || n2 > 1 {
				edgesArr[e.ID].ResetUniqueFlag()
			}
			continue
		}

		ePathArr := copyPathArr(pathArr, edgesPathRelationArr[i])
		checkEdgePathArrDirection(edgesArr, nodesArr, ePathArr, e.ID)
		if edgesArr[e.ID].GetUniqueFlag() == 0 {
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
		/*for x, p := range ePathArr {
			fmt.Printf("[MergePathArr]: ePathArr[%v]: %v\n", x, p)
		}*/

		// find consis path
		var path constructdbg.Path
		path.Freq = math.MaxInt64
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
				suc := true
				var id1, id2 constructdbg.DBG_MAX_INT
				var freq1, freq2 int

				for k := 0; k < len(ePathArr); k++ {
					if ePathArr[k][0][j] > 0 {
						eID1 := ePathArr[k][0][j]
						eID2 := ePathArr[k][1][j]
						//fmt.Printf("[MergePathArr]k: %v, eID1:%v, eID2: %v, j: %v id1: %v, id2: %v, freq1: %v, freq2: %v\n", k, eID1, eID2, j, id1, id2, freq1, freq2)
						if id1 == 0 {
							id1 = eID1
							freq1++
						} else if eID1 == id1 {
							freq1++
						} else if id2 == 0 {
							id2 = eID1
							freq2++
						} else if eID1 == id2 {
							freq2++
						} else {
							suc = false
							break
						}
						if eID2 > 0 {
							if eID2 == id1 {
								freq1++
							} else if id2 == 0 {
								id2 = eID2
								freq2++
							} else if eID2 == id2 {
								freq2++
							} else {
								suc = false
								break
							}
						}
					}
				}
				if !suc {
					break
				}
				if freq1 < freq2 {
					id1, id2 = id2, id1
					freq1, freq2 = freq2, freq1
				}
				if freq1 > minMapFreq && float32(freq1)/10 > float32(freq2) {
					path.IDArr = append(path.IDArr, id1)
					if freq1 < path.Freq {
						path.Freq = freq1
					}
				} else {
					break
				}
			}
			//fmt.Printf("[MergePathArr] path: %v\n", path)

			if len(path.IDArr) < 1 {
				break
			}
			if z == 0 {
				path.IDArr = constructdbg.GetReverseDBG_MAX_INTArr(path.IDArr)
			}
		}

		if len(path.IDArr) >= 2 {
			var pm []constructdbg.Path
			pm = append(pm, path)
			edgesArr[i].PathMat = pm
			fmt.Printf("[MergePathArr]: pm: %v\n", pm)
		} else {
			edgesArr[i].PathMat = nil
		}
	}
}

func SimplifyByLongReadsPath(edgesArr []constructdbg.DBGEdge, nodesArr []constructdbg.DBGNode, pathArr [][2][]constructdbg.DBG_MAX_INT, opt Options) {

	sortedEIDIdxArr := sortIdxByUniqueEdgeLen(edgesArr)
	edgesPathRelationArr := make([][]uint32, len(edgesArr))
	constructEdgesPathRelationship(edgesArr, pathArr, edgesPathRelationArr)

	cleanNGSPath(edgesArr)
	MergePathArr(edgesArr, nodesArr, pathArr, edgesPathRelationArr, opt.MinMapFreq)
	//constructdbg.MergePathMat(edgesArr, nodesArr, 2)
	//fmt.Printf("[SimplifyByLongReadsPath] sortedEIDIdxArr: %v\n", sortedEIDIdxArr)
	mergePathArr := findMaxPath(sortedEIDIdxArr, edgesArr, nodesArr, opt.MinMapFreq)
	// constuct map edge ID to the path
	//IDMapPath := constructdbg.ConstructIDMapPath(pathArr)
	//constructdbg.DeleteJoinPathArrEnd(edgesArr, pathArr)
	constructdbg.ReconstructDBG(edgesArr, nodesArr, mergePathArr, opt.Kmer)
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

	opt := Options{gOpt, 0, 0, 0, 0, 0, ""}
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
	constructdbg.Kmerlen = opt.Kmer
	fmt.Printf("Arguments: %v\n", opt)

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

	uniqueNum := constructdbg.SetDBGEdgesUniqueFlag(edgesArr, nodesArr)
	fmt.Printf("[DeconstructDBG] unique edge number is : %v\n", uniqueNum)

	// remap NGS reads to the new samplify DBG
	var copt constructdbg.Options
	copt.Kmer = opt.Kmer
	copt.NumCPU = opt.NumCPU
	copt.WinSize = opt.WinSize
	copt.MaxNGSReadLen = opt.MaxNGSReadLen
	copt.CfgFn = opt.CfgFn
	wrFn := opt.Prefix + ".decdbg.NGSAlignment"
	constructdbg.MapNGS2DBG(copt, nodesArr, edgesArr, wrFn)
	constructdbg.AddPathToDBGEdge(edgesArr, wrFn)
	constructdbg.MergePathMat(edgesArr, nodesArr, opt.MinMapFreq)

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
	SimplifyByLongReadsPath(edgesArr, nodesArr, pathArr, opt)

	graphfn := opt.Prefix + ".afterLR.dot"
	constructdbg.GraphvizDBGArr(nodesArr, edgesArr, graphfn)
	DcDBGEdgesfn := opt.Prefix + ".edges.DcDBG.fq"
	constructdbg.StoreEdgesToFn(DcDBGEdgesfn, edgesArr)
}
