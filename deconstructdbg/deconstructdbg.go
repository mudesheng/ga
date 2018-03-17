package deconstructdbg

import (
	"fmt"
	"log"
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

func WriteLongPathToDBG(wc chan [2][]constructdbg.DBG_MAX_INT, edgesArr []constructdbg.DBGEdge, numCPU int) {
	var finishNum int
	for {
		tmp := <-wc
		fmt.Printf("[WriteLongPathToDBG] path: %v\n", tmp)
		p := tmp[0]
		if len(p) == 0 {
			finishNum++
			if finishNum == numCPU {
				break
			}
		}
		// write path to the DBG
		// all long reads path store this pathArr, and store index of pathArr to the edge PathMat[1]
		if len(p) > 2 {
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
		}
	}

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

func findMaxPath(sortedEIDIdxArr []IdxLen, edgesArr []constructdbg.DBGEdge, nodesArr []constructdbg.DBGNode) (pathArr []constructdbg.Path) {
	for _, item := range sortedEIDIdxArr {
		e := edgesArr[item.Idx]
		if e.GetDeleteFlag() > 0 || e.GetProcessFlag() > 0 || e.GetUniqueFlag() == 0 || len(e.PathMat) == 0 || len(e.PathMat[0].IDArr) == 0 || constructdbg.IsTwoEdgesCyclePath(edgesArr, nodesArr, e.ID) || constructdbg.FreqNumInDBG_MAX_INTArr(edgesArr[e.ID].PathMat[0].IDArr, e.ID) != 1 {
			continue
		}
		fmt.Printf("[findMaxPath] eID: %v, length: %v\n", item.Idx, item.Length)
		maxP := constructdbg.ExtendPath(edgesArr, nodesArr, e)
		if len(maxP.IDArr) > 1 {
			pathArr = append(pathArr, maxP)
		}
	}
	return pathArr
}

func SimplifyByLongReadsPath(edgesArr []constructdbg.DBGEdge, nodesArr []constructdbg.DBGNode, opt Options) {

	//cleanNGSPath(edgesArr)
	constructdbg.MergePathMat(edgesArr, nodesArr, 2)

	sortedEIDIdxArr := sortIdxByUniqueEdgeLen(edgesArr)
	pathArr := findMaxPath(sortedEIDIdxArr, edgesArr, nodesArr)
	// constuct map edge ID to the path
	//IDMapPath := constructdbg.ConstructIDMapPath(pathArr)
	constructdbg.DeleteJoinPathArrEnd(edgesArr, pathArr)
	//pathArr = constructdbg.ReconstructDBG(edgesArr, nodesArr, pathArr, IDMapPath)
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
	if opt.MinMapFreq < 2 && opt.MinMapFreq >= 10 {
		log.Fatalf("[checkArgs] argument 'MinMapFreq': %v must 2 <= MinMapFreq < 10\n", c.Flag("MinMapFreq").String())
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

	WriteLongPathToDBG(wc, edgesArr, opt.NumCPU)

	// Simplify using Long Reads Mapping info
	SimplifyByLongReadsPath(edgesArr, nodesArr, opt)

	graphfn := opt.Prefix + ".afterLR.dot"
	constructdbg.GraphvizDBGArr(nodesArr, edgesArr, graphfn)
	DcDBGEdgesfn := opt.Prefix + ".edges.DcDBG.fq"
	constructdbg.StoreEdgesToFn(DcDBGEdgesfn, edgesArr)
}
