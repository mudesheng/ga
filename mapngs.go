package main

import (
	"encoding/binary"
	"fmt"
	"io"
	"log"
	"math"
	"os"
	"reflect"
	"runtime"
	"sync"
	"time"

	"github.com/jwaldrip/odin/cli"
	"github.com/klauspost/compress/zstd"
)

/*func (e *DBGEdge) InsertPathToEdge(path []uint32, freq int) {

	// check have been added
	added := false
	for i, v := range e.PathMat {
		rp := GetReverseUint32Arr(path)
		if reflect.DeepEqual(v.IDArr, path) || reflect.DeepEqual(v.IDArr, rp) {
			e.PathMat[i].Freq += freq
			added = true
			break
		}
	}
	if added == false {
		var np Path
		np.IDArr = make([]uint32, len(path))
		copy(np.IDArr, path)
		np.Freq = freq
		e.PathMat = append(e.PathMat, np)
	}
}*/

/*func FindConsisPath(pID uint32, e DBGEdge) (consisP Path) {
	var pm []Path
	for _, pa := range e.PathMat {
		p1 := IndexUint32(pa.IDArr, pID)
		p2 := IndexUint32(pa.IDArr, e.ID)
		if p1 >= 0 {
			var npa Path
			if p1 < p2 {
				npa.IDArr = GetReverseUint32Arr(pa.IDArr[:p2+1])
			} else {
				npa.IDArr = pa.IDArr[p2:]
			}
			npa.Freq = pa.Freq
			pm = append(pm, npa)
		}
	}

	// debug code
	for i, p := range pm {
		fmt.Printf("[FindConsisPath] pm[%v]: %v\n", i, p)
	}

	freq := 1
	for k := 0; freq > 0; k++ {
		freq = 0
		for _, p := range pm {
			if len(p.IDArr) <= k {
				continue
			}
			if len(consisP.IDArr) == k {
				consisP.IDArr = append(consisP.IDArr, p.IDArr[k])
			} else {
				if consisP.IDArr[k] != p.IDArr[k] {
					freq = 0
					consisP.IDArr = consisP.IDArr[:len(consisP.IDArr)-1] // remove last element
					break
				}
			}
			freq += p.Freq
		}
		if freq > 0 {
			consisP.Freq = freq
		}
	}
	return consisP
}*/

/*func ExtendPathCDBG(edgesArr []DBGEdge, nodesArr []DBGNode, e DBGEdge, minMappingFreq int, semi bool) (maxP Path) {
	var p Path
	p.IDArr = make([]uint32, len(e.PathMat[0].IDArr))
	copy(p.IDArr, e.PathMat[0].IDArr)
	p.Freq = e.PathMat[0].Freq
	if p.Freq < minMappingFreq {
		return
	}
	idx := IndexUint32(p.IDArr, e.ID)
	mutualArr := []uint32{e.ID}
	fmt.Printf("[ExtendPath] e.ID: %v,e.PathMat[0]: %v\n", e.ID, e.PathMat[0])
	// found left partition path
	for i := 0; i < idx; i++ {
		e2 := edgesArr[p.IDArr[i]]
		//fmt.Printf("[ExtendPath] i: %v, idx: %v, e2.Process: %v, e2.Delete: %v, e2.Unique: %v\n", i, idx, e2.GetProcessFlag(), e2.GetDeleteFlag(), e2.GetUniqueFlag())
		if len(e2.PathMat) != 1 || e2.PathMat[0].Freq < minMappingFreq || e2.GetProcessFlag() > 0 || e2.GetDeleteFlag() > 0 {
			continue
		}
		if semi {
			if e2.GetUniqueFlag() == 0 {
				continue
			}
		} else {
			if e2.GetUniqueFlag() == 0 {
				continue
			}
		}
		if IsInuint32Arr(mutualArr, e2.ID) {
			fmt.Printf("[ExtendPath] e2.ID: %v in the mutualArr: %v\n", e2.ID, mutualArr)
			continue
		}
		p2 := e2.PathMat[0]
		eID1 := mutualArr[len(mutualArr)-1]
		e1 := edgesArr[eID1]
		j1 := IndexUint32(p2.IDArr, e1.ID)
		j2 := IndexUint32(p2.IDArr, e2.ID)
		fmt.Printf("[ExtendPath]BACKWARD e1 ID: %v,  e2 ID: %v, p2: %v\n", e1.ID, e2.ID, p2)
		if j1 < 0 || j2 < 0 {
			continue
		}
		if j1 < j2 {
			p2.IDArr = GetReverseUint32Arr(p2.IDArr)
			j1 = IndexUint32(p2.IDArr, e1.ID)
			//j2 = IndexUint32(p2.IDArr, e2.ID)
		}
		k1 := IndexUint32(p.IDArr, e1.ID)
		//fmt.Printf("[ExtendPath]j1: %v, j2:%v,k1: %v, eID1: %v, eID2: %v, p: %v, p2: %v\n", j1, j2, k1, e1.ID, e2.ID, p, p2)
		if j1 >= k1 && len(p2.IDArr)-j1 <= len(p.IDArr)-k1 {
			if reflect.DeepEqual(p2.IDArr[j1-k1:], p.IDArr[:len(p2.IDArr)-(j1-k1)]) {
				mutualArr = append(mutualArr, e2.ID)
				na := make([]uint32, len(p.IDArr)+j1-k1)
				copy(na[:j1], p2.IDArr[:j1])
				copy(na[j1:], p.IDArr[k1:])
				p.IDArr = na
				if p2.Freq < p.Freq {
					p.Freq = p2.Freq
				}
				idx = IndexUint32(p.IDArr, e2.ID)
				i = -1
			}
		}
		fmt.Printf("[ExtendPath] p: %v\n", p)
	}

	Reverseuint32Arr(mutualArr)
	idx = IndexUint32(p.IDArr, e.ID)
	// find right path
	for i := len(p.IDArr) - 1; i > idx; i-- {
		e2 := edgesArr[p.IDArr[i]]
		if len(e2.PathMat) != 1 || e2.PathMat[0].Freq < minMappingFreq || e2.GetProcessFlag() > 0 || e2.GetDeleteFlag() > 0 {
			continue
		}
		if semi {
			if e2.GetUniqueFlag() == 0 {
				continue
			}
		} else {
			if e2.GetUniqueFlag() == 0 {
				continue
			}
		}
		if IsInuint32Arr(mutualArr, e2.ID) {
			fmt.Printf("[ExtendPath] e2.ID: %v in the mutualArr: %v\n", e2.ID, mutualArr)
			continue
		}
		p2 := e2.PathMat[0]
		eID1 := mutualArr[len(mutualArr)-1]
		e1 := edgesArr[eID1]
		j1 := IndexUint32(p2.IDArr, e1.ID)
		j2 := IndexUint32(p2.IDArr, e2.ID)
		if j1 < 0 || j2 < 0 {
			continue
		}
		if j2 < j1 {
			p2.IDArr = GetReverseUint32Arr(p2.IDArr)
			j1 = IndexUint32(p2.IDArr, e1.ID)
			//j2 = IndexUint32(p2.IDArr, e2.ID)
		}
		k1 := IndexUint32(p.IDArr, e1.ID)
		//fmt.Printf("[ExtendPath] j1: %v, j2: %v, k1: %v, eID1: %v, eID2: %v, p: %v, p2: %v\n", j1, j2, k1, e1.ID, e2.ID, p, p2)
		fmt.Printf("[ExtendPath]FORWARD e1 ID: %v,  e2 ID: %v, p2: %v\n", e1.ID, e2.ID, p2)
		if len(p2.IDArr)-j1 >= len(p.IDArr)-k1 && k1 >= j1 {
			if reflect.DeepEqual(p.IDArr[k1-j1:], p2.IDArr[:len(p.IDArr)-(k1-j1)]) {
				mutualArr = append(mutualArr, e2.ID)
				if len(p2.IDArr)-j1 > len(p.IDArr)-k1 {
					p.IDArr = append(p.IDArr, p2.IDArr[j1+(len(p.IDArr)-k1):]...)
				}
				if p2.Freq < p.Freq {
					p.Freq = p2.Freq
				}
				idx = IndexUint32(p.IDArr, e2.ID)
				i = len(p.IDArr)
			}
		}
	}
	fmt.Printf("[ExtendPath] p: %v\n", p)

	// not allow two cycle edge in the start or end edge
	start, end := 0, len(mutualArr)-1
	for x, eID := range mutualArr {
		if edgesArr[eID].GetTwoEdgesCycleFlag() > 0 {
			start = x + 1
		} else {
			break
		}
	}
	for x := len(mutualArr) - 1; x > start; x-- {
		if edgesArr[mutualArr[x]].GetTwoEdgesCycleFlag() > 0 {
			end = x - 1
		} else {
			break
		}
	}
	if start >= end {
		edgesArr[e.ID].SetProcessFlag()
		return
	}
	mutualArr = mutualArr[start : end+1]

	// set maxP and process DBG
	i1 := IndexUint32(p.IDArr, mutualArr[0])
	i2 := IndexUint32(p.IDArr, mutualArr[len(mutualArr)-1])
	maxP.IDArr = p.IDArr[i1 : i2+1]
	maxP.Freq = p.Freq

	edgesArr[mutualArr[0]].SetProcessFlag()
	for _, id := range mutualArr[1:] {
		edgesArr[id].SetDeleteFlag()
		edgesArr[id].SetProcessFlag()
	}
	fmt.Printf("[ExtendPath] mutualArr: %v\n", mutualArr)
	fmt.Printf("[ExtendPath] maxP: %v\n", maxP)

	return
}*/

func findPathOverlap(jp Path, pathArr []uint32, edgesArr []DBGEdge) (id uint32, num int) {
	fmt.Printf("[findPathOverlap] jpArr: %v\tpathArr: %v\n", jp, pathArr)
	jpArr := jp.IDArr
	if jpArr[0] == pathArr[0] {
		i := 1
		for ; i < len(jpArr) && i < len(pathArr); i++ {
			if jpArr[i] != pathArr[i] {
				break
			}
		}
		if i == len(jpArr) || i == len(pathArr) {
			if edgesArr[jpArr[0]].GetDeleteFlag() == 0 {
				id = jpArr[0]
			} else if edgesArr[jpArr[len(jpArr)-1]].GetDeleteFlag() == 0 {
				id = jpArr[len(jpArr)-1]
			} else {
				for _, item := range jpArr[1 : len(jpArr)-1] {
					if edgesArr[item].GetDeleteFlag() == 0 {
						if id > 0 {
							log.Fatalf("[findPathOverlap] found more than one edge not deleted in the joinPathMat: %v\n", jpArr)
						}
						id = item
					}
				}
			}
			num = i
			return id, num
		}
	} else if jpArr[len(jpArr)-1] == pathArr[0] {
		rp := GetReverseUint32Arr(jpArr)
		i := 1
		fmt.Printf("[findPathOverlap] rp: %v\tpathArr: %v\n", rp, pathArr)
		for ; i < len(rp) && i < len(pathArr); i++ {
			if rp[i] != pathArr[i] {
				break
			}
		}
		if i == len(rp) || i == len(pathArr) {
			if edgesArr[rp[0]].GetDeleteFlag() == 0 {
				id = rp[0]
			} else if edgesArr[rp[len(rp)-1]].GetDeleteFlag() == 0 {
				id = rp[len(rp)-1]
			} else {
				for _, item := range rp[1 : len(rp)-1] {
					if edgesArr[item].GetDeleteFlag() == 0 {
						if id > 0 {
							log.Fatalf("[findPathOverlap] found more than one edge not deleted in the joinPathMat: %v\n", jpArr)
						}
						id = item
					}
				}
			}
			num = i
			return id, num
		}
	}

	return id, num
}

/*func LoadPathFromFile(pathfn string, cs chan<- MapingInfo) {
	bufSize := (1 << 15)
	fncs := make(chan []byte, 2)
	go ReadBrFile(pathfn, fncs, bufSize)
	buf := make([]byte, 0, 1<<16)
	//lb := make([]byte, 100)
	var lb []byte
	pos := 0
	var suc bool
	for {
		lb, pos, suc = GetNextLine(buf, pos)
		if !suc {
			nb, ok := <-fncs
			if !ok {
				break
			}
			if pos < len(buf) {
				copy(buf[0:], buf[pos:len(buf)])
			}
			buf = buf[:len(buf)-pos]
			pos = 0
			buf = append(buf, nb...)
			continue
		}
		var mi MapingInfo
		colList := SplitEIDArr(lb, '\t')
		id, err := strconv.Atoi(string(colList[0]))
		if err != nil {
			log.Fatalf("[LoadPathFromFile]%s transform read id error\n", colList[0])
		}
		mi.ID = uint32(id)
		mi.Anotition = string(colList[6])
		var pathMat [][][]uint32
		var path []uint32
		//flist := strings.Fields(string(lb))
		//eIDArr := make([]uint32, len(colList))
		list := SplitEIDArr(colList[5], '>')
		for _, ele := range list {
			bubbleArr := SplitEIDArr(ele, '|')
			if len(bubbleArr) == 1 {
				bs := bubbleArr[0]
				eID, err := strconv.Atoi(string(bs))
				if err != nil {
					log.Fatalf("[LoadPathFromFile]%s transform eID error\n", bs)
				}
				path = append(path, uint32(eID))
			} else {
				pathMat = append(pathMat, [][]uint32{path})
				var tmp []uint32
				path = tmp
				pathArr := make([][]uint32, len(bubbleArr))
				for x, path := range bubbleArr {
					arr := SplitEIDArr(path, '-')
					for _, bs := range arr {
						//eID, err := ByteArrInt(id)
						eID, err := strconv.Atoi(string(bs))
						if err != nil {
							log.Fatalf("[LoadPathFromFile]%s transform eID error\n", bs)
						}
						pathArr[x] = append(pathArr[x], uint32(eID))
					}
				}
				pathMat = append(pathMat, pathArr)
			}
		}
		//fmt.Printf("[LoadPathFromFile]eIDArr:%v\n", eIDArr)
		mi.PathMat = pathMat
		cs <- mi
	}
	close(cs)
	/*pathfp, err := os.Open(pathfn)
	if err != nil {
		log.Fatalf("[LoadPathFromFile] file %s create error, err: %v\n", pathfn, err)
	}
	defer pathfp.Close()
	cbrofp := cbrotli.NewReaderSize(pathfp, 1<<20)
	defer cbrofp.Close()
	buffp := bufio.NewReaderSize(cbrofp, 1<<20) // 1<<25 == 2**25
	buf, err := ioutil.ReadAll(buffp)
	if err != nil {
		log.Fatalf("[LoadPathFromFile] read file %s failed, err:%v\n", pathfn, err)
	}
	var a PathArr
	err = proto.Unmarshal(buf, &a)
	if err != nil {
		log.Fatalf("[LoadPathFromFile] proto.Unmarshal() err:%v\n", err)
	}
	readPathArr = a

	return
}*/

/*func MappingNGSRead(dk DBGKmer, ri ReadInfo, rpos int, rstrand bool, edgesArr []DBGEdge, nodesArr []DBGNode, kmerlen int) (errorNum, mappingNum int, path []uint32) {
	var strand bool
	if dk.Strand == rstrand {
		strand = PLUS
		dk.Pos += int32(kmerlen)
	} else {
		strand = MINUS
		dk.Pos--
	}
	rpos += kmerlen
	mappingNum += kmerlen
	edgeMapLen := kmerlen
	for i := rpos; i < len(ri.Seq); {
		e := &edgesArr[dk.ID]
		b := len(ri.Seq)
		var j int
		//fmt.Printf("[MappingReadToEdges]i: %v,  strand: %v, dk.Pos: %v, dk.ID: %v, errorNum: %v\n", i, strand, dk.Pos, dk.ID, errorNum)
		if strand == PLUS {
			if len(e.Ks)-int(dk.Pos) < len(ri.Seq)-i {
				b = i + (len(e.Ks) - int(dk.Pos))
			}
			j = int(dk.Pos)
			for ; i < b && j < e.GetSeqLen(); i++ {
				if ri.Seq[i] != e.Ks[j] {
					errorNum++
				} else {
					mappingNum++
					edgeMapLen++
				}
				j++
			}
		} else { // strand == MINUS
			if len(ri.Seq)-i > int(dk.Pos)+1 {
				b = i + (int(dk.Pos) + 1)
			}
			j = int(dk.Pos)
			for ; i < b && j >= 0; i++ {
				//fmt.Printf("[MappingReadToEdges] b: %v,i: %v,  j: %v\n", b, i, j)
				if ri.Seq[i] != BntRev[e.Ks[j]] {
					errorNum++
				} else {
					mappingNum++
					edgeMapLen++
				}
				j--
			}
		}
		path = append(path, e.ID)
		if i >= len(ri.Seq) {
			break
		}
		// find next edge
		{
			eIDArr := GetNextEdgeIDArr(e, nodesArr, strand, FORWARD)
			if len(eIDArr) == 0 {
				break
			} else if len(eIDArr) == 1 {
				dk.ID = eIDArr[0]
				ne := &edgesArr[dk.ID]
				strand = GetNextEdgeStrand(e, ne, nodesArr, strand)
				if strand {
					dk.Pos = int32(kmerlen) - 1
				} else {
					dk.Pos = int32(ne.GetSeqLen()) - int32(kmerlen)
				}
				continue
			}
			ok := false
			for x := 0; x < len(eIDArr); x++ {
				dk.ID = eIDArr[x]
				ne := &edgesArr[dk.ID]
				std := GetNextEdgeStrand(e, ne, nodesArr, strand)
				var eb byte
				if std {
					dk.Pos = int32(kmerlen) - 1
					eb = ne.Ks[dk.Pos]
				} else {
					dk.Pos = int32(ne.GetSeqLen() - kmerlen)
					eb = BntRev[ne.Ks[dk.Pos]]
				}
				if eb == ri.Seq[i] {
					strand = std
					ok = true
					break
				}
			}
			if ok {
				edgeMapLen = 0
				continue
			} else {
				break
			}
		}
	}
	if edgeMapLen < 4 && len(path) > 0 {
		path = path[:len(path)-1]
	}

	return
}*/

func MapNGSReadFindPath(RIPoolChan <-chan RIPool, rbytesPool, RIArrPool *sync.Pool, wc chan<- []uint32, pathArrPool *sync.Pool, nodesArr []DBGNode, edgesArr []DBGEdge, cf CuckooFilterDBGKmer, opt optionsMN) {
	var notFoundSeedNum, highErrorNum, shortMappingNum, allNum, oneMoreLocationNum int
	var totalMapLen int
	ka := make([]DBGKmerSeq, opt.WinSize)
	for i := range ka {
		ka[i].ks = make([]byte, opt.Kmer)
	}
	rSeq := make([]byte, opt.MaxNGSReadLen)
	dbgKArr := make([]DBGKmer, 2)
	var rdbgK DBGKmer
	pa := pathArrPool.Get().([]uint32)
	pa = pa[:0]
	for {
		riPool, ok := <-RIPoolChan
		if !ok {
			var t []uint32
			wc <- t
			break
		}
		allNum += len(riPool.RIArr)
		for _, ri := range riPool.RIArr {
			// found kmer seed position in the DBG edges
			if len(ri.Seq) < opt.Kmer+opt.WinSize*2 {
				continue
			}
			//fmt.Printf("[paraMapNGSAndMerge]ri:%v\n", ri)
			dbgKArr, rdbgK = LocateSeedKmerCF(cf, ri, opt.WinSize, opt.Kmer, edgesArr, ka, rSeq, dbgKArr)
			if len(dbgKArr) == 0 { // not found in the cuckoofilter
				//fmt.Printf("[paraMapNGSAndMerge]read ID:%d rdbgK:%v not found seed!!!\n", ri.ID, rdbgK)
				notFoundSeedNum++
				continue
			} else if len(dbgKArr) > 1 {
				dbgKArr[0] = ReLocation(dbgKArr, rdbgK, ri, edgesArr, opt.Kmer)
				if dbgKArr[0].ID == 0 {
					oneMoreLocationNum++
					//fmt.Fprintf(os.Stderr, "[paraMapNGSAndCorrect]read ID:%d ReLocation failed!!!\n", ri.ID)
					continue
				}
			}
			errorNum, mappingNum, rmi, _ := MappingReadToEdges(dbgKArr[0], rdbgK, ri, edgesArr, opt.Kmer)
			totalMapLen += mappingNum
			if mappingNum < len(ri.Seq)*7/10 {
				shortMappingNum++
				continue
			}
			if errorNum > len(ri.Seq)*1/10 {
				highErrorNum++
				continue
			}
			if len(rmi.PathSeqArr) < 2 {
				continue
			}
			if len(pa)+len(rmi.PathSeqArr)+1 >= cap(pa) {
				wc <- pa
				pa = pathArrPool.Get().([]uint32)
				pa = pa[:0]
			}
			for _, ps := range rmi.PathSeqArr {
				pa = append(pa, ps.ID)
			}
			pa = append(pa, 0)
		}
		rbytesPool.Put(riPool.Cs)
		RIArrPool.Put(riPool.RIArr)
	}
	if len(pa) > 0 {
		wc <- pa
	}
	fmt.Printf("[MapNGSReadFindPath] not found seed read pair number is:%d allNum:%d  percent:%f\n", notFoundSeedNum, allNum, float32(notFoundSeedNum)/float32(allNum))
	fmt.Printf("[MapNGSReadFindPath] high error mapping read pair number is:%d short mapping num:%d avg mapping len:%d\n", highErrorNum, shortMappingNum, totalMapLen/(allNum-notFoundSeedNum+1))
}

func ParaMapReadsFile(fn string, wc chan<- []uint32, pathArrPool *sync.Pool, concurrentNum int, nodesArr []DBGNode, edgesArr []DBGEdge, cf CuckooFilterDBGKmer, opt optionsMN) {
	bufSize := WindowSize
	cs := make(chan []byte, 2)
	rbytesPool := sync.Pool{New: func() interface{} {
		bytes := make([]byte, bufSize)
		return bytes
	}}

	size := bufSize / opt.MinNGSReadLen
	RIArrPool := sync.Pool{New: func() interface{} {
		arr := make([]ReadInfo, size)
		return arr
	}}
	RIPoolChan := make(chan RIPool, concurrentNum)

	for i := 0; i < concurrentNum; i++ {
		go MapNGSReadFindPath(RIPoolChan, &rbytesPool, &RIArrPool, wc, pathArrPool, nodesArr, edgesArr, cf, opt)
	}
	go ReadZstdFile(fn, &rbytesPool, cs)
	GetReadInfoBucket(fn, cs, &RIArrPool, RIPoolChan, false)
	runtime.GC()
}

type IDFreq struct {
	ID   uint32
	Freq uint32
}

func AddIDFreqArr(idFreqArr []IDFreq, eID uint32, freq uint32) []IDFreq {
	ok := false
	for i, idfreq := range idFreqArr {
		if idfreq.ID == eID {
			idFreqArr[i].Freq += idfreq.Freq
			ok = true
			break
		}
	}
	if !ok {
		var t IDFreq
		t.Freq = freq
		t.ID = eID
		idFreqArr = append(idFreqArr, t)
	}
	return idFreqArr
}

type PathFreqExt struct {
	PF         PathFreq
	Start, End int
}

func DeleteLowFreq(idFreqArr []IDFreq, MinMapFreq int) []IDFreq {
	count := 0
	for _, idFreq := range idFreqArr {
		if int(idFreq.Freq) >= MinMapFreq {
			idFreqArr[count] = idFreq
			count++
		}
	}
	idFreqArr = idFreqArr[:count]
	return idFreqArr
}

func GetNumFreqInPathFreqArr(pa []PathFreq, eID uint32) (cn int) {
	cn = 0
	for _, pf := range pa {
		for j, id := range pf.IDArr {
			if id == eID {
				cn += int(pf.FreqArr[j])
			}
		}
	}

	return
}

func GetNumInPathFreqArr(pa []PathFreq, eID uint32) (cn int, idx int) {
	idx = -1
	for i, pf := range pa {
		for _, id := range pf.IDArr {
			if id == eID {
				cn++
				idx = i
				break
			}
		}
	}
	return
}

func CutLowFreq(pf PathFreq, MinMapFreq uint8) (npf PathFreq) {
	low, high := 0, len(pf.IDArr)-1
	for ; low < len(pf.IDArr); low++ {
		if pf.FreqArr[low] >= MinMapFreq {
			break
		}
	}
	for ; high >= low; high-- {
		if pf.FreqArr[high] >= MinMapFreq {
			break
		}
	}
	if low <= high {
		npf.IDArr = pf.IDArr[low : high+1]
	}
	return
}

/*func AdjustEIDPath(e DBGEdge, joinPathArr []Path, IDMapPath map[uint32]uint32) (p Path) {
	ep := e.PathMat[0]
	idx := IndexUint32(ep.IDArr, e.ID)
	p.Freq = ep.Freq
	fmt.Printf("[AdjustEIDPath] e: %v\n", e)
	if x, ok := IDMapPath[e.ID]; ok {
		a1, a2 := ep.IDArr, joinPathArr[x].IDArr
		fmt.Printf("[AdjustEIDPath] ep: %v\n\tjoinPathArr: %v\n", ep, joinPathArr[x])
		if idx < len(a1)-1 && a1[idx] == a2[0] && a1[idx+1] == a2[1] {
			l := len(a2)
			if l > len(a1[idx:]) {
				l = len(a1[idx:])
			}
			if reflect.DeepEqual(a1[idx:idx+l], a2[:l]) {
				p.IDArr = append(p.IDArr, a1[:idx+1]...)
				//p.IDArr = append(p.IDArr, e.ID)
				p.IDArr = append(p.IDArr, a1[idx+l:]...)
			}
		} else if idx < len(a1)-1 && a1[idx] == a2[len(a2)-1] && a1[idx+1] == a2[len(a2)-2] {
			a2 = GetReverseUint32Arr(a2)
			l := len(a2)
			if l > len(a1[idx:]) {
				l = len(a1[idx:])
			}
			if reflect.DeepEqual(a1[idx:idx+l], a2[:l]) {
				p.IDArr = append(p.IDArr, a1[:idx+1]...)
				//p.IDArr = append(p.IDArr, e.ID)
				p.IDArr = append(p.IDArr, a1[idx+l:]...)
			}
		} else if idx > 0 && a1[idx] == a2[0] && a1[idx-1] == a2[1] {
			a2 := GetReverseUint32Arr(a2)
			l := len(a2)
			if l > len(a1[:idx+1]) {
				l = len(a1[:idx+1])
			}
			if reflect.DeepEqual(a1[idx-l+1:idx+1], a2[len(a2)-l:]) {
				p.IDArr = append(p.IDArr, a1[:idx-l+1]...)
				p.IDArr = append(p.IDArr, e.ID)
				p.IDArr = append(p.IDArr, a1[idx+1:]...)
			}
		} else if idx > 0 && a1[idx] == a2[len(a2)-1] && a1[idx-1] == a2[len(a2)-2] {
			l := len(a2)
			if l > len(a1[:idx+1]) {
				l = len(a1[:idx+1])
			}
			if reflect.DeepEqual(a1[idx+1-l:idx+1], a2[len(a2)-l:]) {
				p.IDArr = append(p.IDArr, a1[:idx-l+1]...)
				p.IDArr = append(p.IDArr, e.ID)
				p.IDArr = append(p.IDArr, a1[idx+1:]...)
			}
		}
	}
	if len(p.IDArr) == 0 {
		p.IDArr = ep.IDArr
	}
	fmt.Printf("[AdjustEIDPath] p: %v\n", p)
	return p
}*/

func IsContainBoundaryArr(a1, a2 []uint32) bool {
	if reflect.DeepEqual(a1[:len(a2)], a2) {
		return true
	}
	if reflect.DeepEqual(a1[len(a1)-len(a2):], a2) {
		return true
	}
	ra2 := GetReverseUint32Arr(a2)
	if reflect.DeepEqual(a1[:len(ra2)], ra2) {
		return true
	}
	if reflect.DeepEqual(a1[len(a1)-len(ra2):], ra2) {
		return true
	}
	return false
}

func GetReverseEdgeFreqArr(arr, rArr []EdgeFreq) []EdgeFreq {
	la := len(arr)
	if cap(rArr) < la {
		rArr = make([]EdgeFreq, la)
	} else {
		rArr = rArr[:la]
	}
	for i, ef := range arr {
		rArr[la-1-i] = ef
	}
	return rArr
}

func ReverseEdgeFreqArr(arr []EdgeFreq) []EdgeFreq {
	al := len(arr)
	for i := 0; i < al/2; i++ {
		arr[i], arr[al-1-i] = arr[al-1-i], arr[i]
	}

	return arr
}

func IsUniqueEdge(e *DBGEdge) bool {
	if e.GetUniqueFlag() > 0 || e.GetTwoEdgesCycleFlag() > 0 || e.GetBubbleFlag() > 0 {
		return true
	}
	return false
}

func CountEIDArr(eIDArr []uint32, eID uint32) (count int) {
	for _, id := range eIDArr {
		if id == eID {
			count++
		}
	}
	return
}

func GetExtendUniquePath(e *DBGEdge, edgesArr []DBGEdge) (efArr []EdgeFreq, eIDArr []uint32) {
	max := 1000
	// extend left path
	ep := e.NGSPathArr[0]
	rArr := make([]EdgeFreq, 0, len(ep)*3)
	shareArr := make([]EdgeFreq, 0, 20)
	eIdx := IndexEdgeFreq(ep, e.ID)
	efArr = GetReverseEdgeFreqArr(ep[:eIdx+1], rArr)
	eIDArr = make([]uint32, 0, 5)
	eIDArr = append(eIDArr, e.ID)
	//fmt.Printf("[GetExtendUniquePath]eID:%v, ep:%v\n", e.ID, ep)
	//idx := IndexEdgeFreq(efArr, eIDArr[len(eIDArr)-1])
	for i := 1; i < len(efArr) && i < max; i++ {
		eID := efArr[i].ID
		te := &edgesArr[eID]
		if !IsUniqueEdge(te) || len(te.NGSPathArr) != 1 {
			continue
		}

		if CountEIDArr(eIDArr, te.ID) > 0 {
			break
		}
		if te.GetProcessFlag() > 0 {
			//idx := IndexEdgeFreq(efArr[:i], te.ID)
			for j, id := range eIDArr {
				fmt.Printf("[GetExtendUniquePath]error! eIDArr[%d]:%d NGSPath:%v\n", j, id, edgesArr[id].NGSPathArr)
			}
			fmt.Printf("[GetExtendUniquePath]error! eID:%v has processed,efArr:%v NGSPath:%v\n", te.ID, efArr, te.NGSPathArr)
			break
			//log.Fatalf("[GetExtendUniquePath]error! eID:%v has processed\n", te.ID)
		}
		lastEID := eIDArr[len(eIDArr)-1]
		ep1 := te.NGSPathArr[0]
		idx1 := IndexEdgeFreq(ep1, lastEID)
		idx2 := IndexEdgeFreq(ep1, eID)
		if idx1 >= 0 && idx2 >= 0 {
			idx3 := IndexEdgeFreq(efArr, lastEID)
			if idx1 > idx2 {
				shareArr = GetReverseEdgeFreqArr(ep1, shareArr)
				ep1 = shareArr
				idx1 = len(ep1) - 1 - idx1
			}
			if len(efArr)-idx3 > len(ep1)-idx1 {
				break
			}
			if !EqualEdgeFreqArr(efArr[idx3:], ep1[idx1:idx1+len(efArr)-idx3]) {
				break
				//log.Fatalf("[GetExtendUniquePath]error! i:%d efArr:%v and ep1:%v not extend\n", i, efArr, ep1)
			}
			efArr = append(efArr[:idx3], ep1[idx1:]...)
			//fmt.Printf("[GetExtendUniquePath]eID:%v, ep1:%v\n\tefArr:%v\n", eID, ep1, efArr)
			eIDArr = append(eIDArr, te.ID)
		}
	}
	efArr = ReverseEdgeFreqArr(efArr)
	eIDArr = GetReverseUint32Arr(eIDArr)
	if eIdx < len(ep)-1 {
		efArr = append(efArr, ep[eIdx+1:]...)
	}

	// extend right path
	leftIdx := IndexEdgeFreq(efArr, eIDArr[len(eIDArr)-1])
	for i := leftIdx + 1; i < len(efArr) && i < max; i++ {
		eID := efArr[i].ID
		te := &edgesArr[eID]
		if !IsUniqueEdge(te) || len(te.NGSPathArr) != 1 {
			continue
		}
		if CountEIDArr(eIDArr, te.ID) > 0 {
			break
		}
		if te.GetProcessFlag() > 0 {
			//idx := IndexEdgeFreq(efArr[:i], te.ID)
			fmt.Printf("[GetExtendUniquePath]error! eID:%v has processed,efArr:%v NGSPath:%v \n", te.ID, efArr, te.NGSPathArr)
			break
			//log.Fatalf("[GetExtendUniquePath]error! eID:%v has processed\n", te.ID)
		}
		lastEID := eIDArr[len(eIDArr)-1]
		ep1 := te.NGSPathArr[0]
		idx1 := IndexEdgeFreq(ep1, lastEID)
		idx2 := IndexEdgeFreq(ep1, eID)
		if idx1 >= 0 && idx2 >= 0 {
			idx3 := IndexEdgeFreq(efArr, lastEID)
			if idx1 > idx2 {
				shareArr = GetReverseEdgeFreqArr(ep1, shareArr)
				ep1 = shareArr
				idx1 = len(ep1) - 1 - idx1
			}
			if len(efArr)-idx3 > len(ep1)-idx1 {
				break
			}
			if !EqualEdgeFreqArr(efArr[idx3:], ep1[idx1:idx1+len(efArr)-idx3]) {
				break
				//log.Fatalf("[GetExtendUniquePath]error! i:%d efArr:%v and ep1:%v not extend\n", i, efArr, ep1)
			}
			efArr = append(efArr[:idx3], ep1[idx1:]...)
			//fmt.Printf("[GetExtendUniquePath]eID:%v, ep1:%v efArr:%v\n", te.ID, ep1, efArr)
			eIDArr = append(eIDArr, te.ID)
		}
	}
	//FLANK_ALLOW_FACTORfmt.Printf("[GetExtendUniquePath]efArr:%v\n", efArr)

	if len(efArr) < len(ep) {
		efArr = nil
		eIDArr = nil
	}

	if len(efArr) >= max {
		fmt.Fprintf(os.Stderr, "[GetExtendUniquePath]e.NGSPathArr[0]:%v extend too long\n", e.NGSPathArr[0])
		efArr = nil
		eIDArr = nil
	}
	return
}

func ConcatEdgePath(path []uint32, edgesArr []DBGEdge, nodesArr []DBGNode, kmerlen int) (u Unitig) {
	sl := GetPathSeqLen(path, edgesArr, kmerlen)
	fmt.Printf("[ConcatEdgePath]path:%v seqLen:%d\n", path, sl)
	kn := kmerlen - 1
	u.Ks = make([]byte, sl)
	var pos int
	strand := PLUS
	e := &edgesArr[path[0]]
	nID := GetShareNodeID(e, &edgesArr[path[1]])
	var direction uint8
	if nID == e.StartNID {
		direction = BACKWARD
		pos = sl - e.GetSeqLen()
		copy(u.Ks[pos:], e.Ks)
	} else {
		direction = FORWARD
		copy(u.Ks, e.Ks)
		pos = e.GetSeqLen()
	}
	fmt.Printf("[ConcatEdgePath]direction:%v pos:%d\n\teID:%d %s\n", direction, pos, e.ID, Transform2Char(e.Ks))
	for i := 1; i < len(path); i++ {
		eID := path[i]
		ne := &edgesArr[eID]
		nl := ne.GetSeqLen()
		strand = GetNextEdgeStrand2(e, ne, strand, FORWARD)
		if direction == FORWARD {
			if strand == PLUS {
				copy(u.Ks[pos:], ne.Ks[kn:])
			} else {
				copy(u.Ks[pos:], ne.Ks[:nl-kn])
				ReverseCompByteArr(u.Ks[pos : pos+nl-kn])
			}
			pos += nl - kn
		} else {
			tl := nl - kn
			if strand == PLUS {
				copy(u.Ks[pos-tl:pos], ne.Ks[:tl])
			} else {
				copy(u.Ks[pos-tl:pos], ne.Ks[nl-tl:])
				ReverseCompByteArr(u.Ks[pos-tl : pos])
				ReverseUint8Arr(u.Kq[pos-tl : pos])
			}
			pos -= tl
		}
		if strand == PLUS {
			fmt.Printf("[ConcatEdgePath2]strand:%v pos:%d ne.StartNID:%d ne.EndNID:%d eID:%d %s\n", strand, pos, ne.StartNID, ne.EndNID, ne.ID, Transform2Char(ne.Ks))
		} else {
			fmt.Printf("[ConcatEdgePath2]strand:%v pos:%d ne.StartNID:%d ne.EndNID:%d eID:%d %s\n", strand, pos, ne.StartNID, ne.EndNID, ne.ID, Transform2Char(GetReverseCompByteArr(ne.Ks)))
		}
		e = ne
	}
	fmt.Printf("[ConcatEdgePath]path:%v %s\n", path, Transform2Char(u.Ks))

	if direction == FORWARD {
		if pos != sl {
			fmt.Fprintf(os.Stderr, "[ConcatEdgePath]pos:%d path:%v\n", pos, path)
		}
	} else {
		if pos != 0 {
			fmt.Fprintf(os.Stderr, "[ConcatEdgePath]pos:%d path:%v\n", pos, path)
		}
	}

	return
}

func ConcatEdgePath2(path []uint32, edgesArr []DBGEdge, nodesArr []DBGNode, kmerlen int) (u Unitig) {
	sl := GetPathSeqLen(path, edgesArr, kmerlen)
	fmt.Printf("[ConcatEdgePath2]path:%v seqLen:%d\n", path, sl)
	kn := kmerlen - 1
	u.Ks = make([]byte, sl)
	u.Kq = make([]uint8, sl)
	var pos int
	strand := PLUS
	e := &edgesArr[path[0]]
	nID := GetShareNodeID(e, &edgesArr[path[1]])
	var direction uint8
	if nID == e.StartNID {
		direction = BACKWARD
		pos = sl - e.GetSeqLen()
		copy(u.Ks[pos:], e.Ks)
	} else {
		direction = FORWARD
		copy(u.Ks, e.Ks)
		pos = e.GetSeqLen()
	}
	rSeq := make([]byte, 0, 500)
	fmt.Printf("[ConcatEdgePath2]direction:%v pos:%d\n\teID:%d %s\n", direction, pos, e.ID, Transform2Char(e.Ks))
	for i := 1; i < len(path); i++ {
		eID := path[i]
		ne := &edgesArr[eID]
		nl := ne.GetSeqLen()
		strand = GetNextEdgeStrand2(e, ne, strand, FORWARD)
		ks := ne.Ks
		if strand == MINUS {
			rSeq = GetReverseCompByteArr2(ks, rSeq)
			ks = rSeq
		}
		if direction == FORWARD {
			copy(u.Ks[pos:], ks[kn:])
			pos += nl - kn
		} else {
			tl := nl - kn
			copy(u.Ks[pos-tl:pos], ks[:tl])
			pos -= tl
		}
		fmt.Printf("[ConcatEdgePath2]strand:%v pos:%d ne.StartNID:%d ne.EndNID:%d eID:%d %s\n", strand, pos, ne.StartNID, ne.EndNID, ne.ID, Transform2Char(ks))
		e = ne
	}
	fmt.Printf("[ConcatEdgePath2]path:%v %s\n", path, Transform2Char(u.Ks))

	if direction == FORWARD {
		if pos != sl {
			fmt.Fprintf(os.Stderr, "[ConcatEdgePath2]pos:%d path:%v\n", pos, path)
		}
	} else {
		if pos != 0 {
			fmt.Fprintf(os.Stderr, "[ConcatEdgePath2]pos:%d path:%v\n", pos, path)
		}
	}

	return
}

func IndexNextZeroUint32Arr(arr []uint32, direction uint8) (idx int) {
	idx = -1
	if direction == FORWARD {
		for i, id := range arr {
			if id == 0 {
				idx = i
				return
			}
		}
	} else {
		for i := len(arr) - 1; i >= 0; i-- {
			if arr[i] == 0 {
				idx = i
				return
			}
		}
	}
	return
}

func ReadZstNGSPathFile(zstfn string, arrPool *sync.Pool, cs chan<- []uint32, ngsPathSize int) {
	fp, err1 := os.Open(zstfn)
	if err1 != nil {
		log.Fatalf("[ReadZstNGSPathFile] open file:%s failed..., err: %v\n", zstfn, err1)
	}
	defer fp.Close()
	//size := (1 << 26)
	//brfp1 := cbrotli.NewReader(fp1)
	zr, err2 := zstd.NewReader(fp, zstd.WithDecoderConcurrency(1))
	if err2 != nil {
		log.Fatalf("[ReadZstNGSPathFile] zstd open file:%s failed..., err: %v\n", zstfn, err1)
	}
	defer zr.Close()
	//buffp := brfp
	//buffp := bufio.NewReaderSize(zr, 1<<16)
	var fullCount, emptyCount int
	var err error
	var num int
	var pa [WindowSize / 4]uint32
	InitUint32Slice(&pa)
	for err != io.EOF {
		//buf := make([]byte, size)
		//num, err = buffp.Read(buf)
		err = binary.Read(zr, binary.LittleEndian, &pa)
		buf := arrPool.Get().([]uint32)
		copy(buf, pa[:])
		if err != nil {
			if err != io.EOF {
				log.Fatalf("[ReadZstNGSPathFile] read file:%s, found err: %v\n", zstfn, err)
			} else {
				continue
			}
		}
		if len(cs) == cap(cs) {
			fullCount++
		} else if len(cs) == 0 {
			emptyCount++
		}
		cs <- buf
		//i := IndexNextZeroUint32Arr(buf, BACKWARD)
	}
	close(cs)
	fmt.Printf("[ReadZstNGSPathFile]readNum:%d ngsPathSize:%d\n", num, ngsPathSize)
	fmt.Printf("[ReadZstNGSPathFile]emptyCount:%d fullCount:%d\n", emptyCount, fullCount)
	//time.Sleep(time.Second * 10)
	return
}

func GetProcessFlag(flag *uint8) uint8 {
	return *flag & 0x1
}

func SetProcessFlag(flag *uint8) {
	*flag = *flag | 0x1
}

func ResetProcessFlag(flag *uint8) {
	*flag = *flag & (0xFF - 0x1)
}
func ParaSetProcessFlag(flag *uint8) bool {
	var oc, fc uint8
	for {
		oc = *flag
		if GetProcessFlag(&oc) > 0 {
			continue
		}
		fc = oc
		SetProcessFlag(&fc)
		if CompareAndSwapUint8(flag, oc, fc) {
			return true
		}
		//fmt.Printf("[insert] oc:%d fc:%d\n", oc, fc)
	}
	return false
}

func ParaResetProcessFlag(flag *uint8) bool {
	var oc, fc uint8
	for {
		oc = *flag
		if GetProcessFlag(&oc) == 0 {
			log.Fatalf("[ParaResetProcessFlag] e.Flag:%d Process flag ==0\n", oc)
		}
		fc = oc
		ResetProcessFlag(&fc)
		if CompareAndSwapUint8(flag, oc, fc) {
			return true
		}
		//fmt.Printf("[insert] oc:%d fc:%d\n", oc, fc)
	}
	return false
}

func EqualEdgeFreqArr(arr1, arr2 []EdgeFreq) bool {
	if len(arr1) != len(arr2) || len(arr1) == 0 || len(arr2) == 0 {
		return false
	}

	ok := true
	for i, c := range arr1 {
		if c.ID != arr2[i].ID {
			ok = false
			break
		}
	}
	return ok
}

func EqualEdgeFreqArrPathNoDirection(arr []EdgeFreq, path []uint32) bool {
	if len(arr) != len(path) || len(arr) == 0 || len(path) == 0 {
		return false
	}

	ok := true
	if arr[0].ID == path[0] && arr[len(arr)-1].ID == path[len(path)-1] {
		for i := 1; i < len(arr)-1; i++ {
			if arr[i].ID != path[i] {
				ok = false
				break
			}
		}
	} else if arr[0].ID == path[len(path)-1] && arr[len(arr)-1].ID == path[0] {
		for i := 1; i < len(arr)-1; i++ {
			if arr[i].ID != path[len(path)-1-i] {
				ok = false
				break
			}
		}
	} else {
		ok = false
	}
	return ok
}

func EqualEdgeFreqArrNoDirection(arr1, arr2 []EdgeFreq) bool {
	if len(arr1) != len(arr2) || len(arr1) == 0 || len(arr2) == 0 {
		return false
	}

	ok := true
	if arr1[0].ID == arr2[0].ID && arr1[len(arr1)-1].ID == arr2[len(arr2)-1].ID {
		for i := 1; i < len(arr1)-1; i++ {
			if arr1[i].ID != arr2[i].ID {
				ok = false
				break
			}
		}
	} else if arr1[0].ID == arr2[len(arr2)-1].ID && arr1[len(arr1)-1].ID == arr2[0].ID {
		for i := 1; i < len(arr1)-1; i++ {
			if arr1[i].ID != arr2[len(arr2)-1-i].ID {
				ok = false
				break
			}
		}
	} else {
		ok = false
	}
	return ok
}

func AddEdgeFreqArr(arr1, arr2 []EdgeFreq) bool {
	if len(arr1) != len(arr2) || len(arr1) == 0 || len(arr2) == 0 {
		return false
	}
	ok := true
	if arr1[0].ID == arr2[0].ID && arr1[len(arr1)-1].ID == arr2[len(arr2)-1].ID {
		for i := 0; i < len(arr1); i++ {
			arr1[i].Freq += arr2[i].Freq
		}
	} else if arr1[0].ID == arr2[len(arr2)-1].ID && arr1[len(arr1)-1].ID == arr2[0].ID {
		for i := 0; i < len(arr1); i++ {
			arr1[i].Freq += arr2[len(arr2)-1-i].Freq
		}
	} else {
		ok = false
	}
	return ok
}

func GetPathFreq(pf PathFreq) (freq uint8) {
	freq = 1
	for _, f := range pf.FreqArr {
		if f > freq {
			freq = f
			break
		}
	}
	return
}

func IndexEdgeFreq(path []EdgeFreq, eID uint32) (idx int) {
	idx = -1
	for i, ef := range path {
		if ef.ID == eID {
			idx = i
			break
		}
	}

	return
}

func EdgeFreqArr2Uint32Arr(ep []EdgeFreq) []uint32 {
	arr := make([]uint32, len(ep))
	for i, ef := range ep {
		arr[i] = ef.ID
	}

	return arr
}

func DeleteEID(ea []uint32, eID uint32) []uint32 {
	if len(ea) == 0 {
		return ea
	}
	na := make([]uint32, 0, len(ea))
	for _, id := range ea {
		if id == eID {
			continue
		}
		na = append(na, id)
	}
	if len(na) == len(ea) {
		log.Fatalf("[DeleteEID] eID: %v not include ea: %v\n", eID, ea)
	}
	return na
}

func RePlaceNGSPath2(path, mergePath, revMP []EdgeFreq) []EdgeFreq {
	eID1, eID2 := mergePath[0].ID, mergePath[len(mergePath)-1].ID
	idx1 := IndexEdgeFreq(path, eID1)
	idx2 := IndexEdgeFreq(path, eID2)
	if idx1 >= 0 || idx2 >= 0 {
		//fmt.Printf("[RePlaceNGSPath2]idx1:%d idx2:%d path:%v mergePath:%v\n", idx1, idx2, path, mergePath)
	}
	if idx1 >= 0 && idx2 >= 0 {
		if idx1 < idx2 {
			if EqualEdgeFreqArr(path[idx1:idx2+1], mergePath) {
				if idx2 < len(path)-1 {
					copy(path[idx1+1:], path[idx2+1:])
				}
				path = path[:idx1+1+len(path)-1-idx2]
			}
		} else {
			if EqualEdgeFreqArr(path[idx2:idx1+1], revMP) {
				path[idx2].ID = eID1
				if idx1 < len(path)-1 {
					copy(path[idx2+1:], path[idx1+1:])
				}
				path = path[:idx2+1+len(path)-1-idx1]
			}
		}
	} else if idx1 >= 0 && idx2 < 0 {
		if idx1 < len(path)-1 && len(path)-idx1 <= len(mergePath) && EqualEdgeFreqArr(path[idx1:], mergePath[:len(path)-idx1]) {
			path = path[:idx1+1]
		} else if idx1+1 <= len(revMP) && EqualEdgeFreqArr(revMP[len(revMP)-(idx1+1):], path[:idx1+1]) {
			path = path[idx1:]
		} else if idx1 == len(path)-1 && len(path)-idx1 <= len(mergePath) && EqualEdgeFreqArr(path[idx1:], mergePath[:len(path)-idx1]) {
			path = path[:idx1+1]
		}
	} else if idx1 < 0 && idx2 >= 0 {
		if idx2 > 0 && idx2+1 <= len(mergePath) && EqualEdgeFreqArr(mergePath[len(mergePath)-(idx2+1):], path[:idx2+1]) {
			path[idx2].ID = eID1
			path = path[idx2:]
		} else if len(path)-idx2 <= len(revMP) && EqualEdgeFreqArr(revMP[:len(path)-idx2], path[idx2:]) {
			path[idx2].ID = eID1
			path = path[:idx2+1]
		} else if idx2 == 0 && idx2+1 <= len(mergePath) && EqualEdgeFreqArr(mergePath[len(mergePath)-(idx2+1):], path[:idx2+1]) {
			path[idx2].ID = eID1
			path = path[idx2:]
		}
	}
	if idx1 >= 0 || idx2 >= 0 {
		//fmt.Printf("[RePlaceNGSPath2]idx1:%d idx2:%d path:%v\n", idx1, idx2, path)
		//fmt.Printf("\tchanged path:%v\n", path)
	}
	return path
}

func GetOverlapPath(p1, rP1, p2 []EdgeFreq, idx1, idx1R, idx2, overlapLen int) []EdgeFreq {
	var catArr []EdgeFreq
	ok := true
	if idx2 > overlapLen-1 {
		if idx2 < idx1 && len(p1)-idx1 < len(p2)-idx2 {
			if overlapLen == 1 {
				if idx2+1 < len(p2) && idx1-1 >= 0 && p2[idx2+1].ID == p1[idx1-1].ID {
					ok = false
				}
			}
			if ok && EqualEdgeFreqArr(p1[idx1-idx2:], p2[:idx2+len(p1)-idx1]) {
				catArr = make([]EdgeFreq, len(p1)+len(p2)-(idx2+len(p1)-idx1))
				copy(catArr, p1)
				for i, ef := range p2[:idx2+len(p1)-idx1] {
					if catArr[idx1-idx2+i].Freq+ef.Freq < math.MaxUint8 {
						catArr[idx1-idx2+i].Freq += ef.Freq
					} else {
						catArr[idx1-idx2+i].Freq = math.MaxUint8
					}
				}
				copy(catArr[len(p1):], p2[idx2+len(p1)-idx1:])
			}
		}
		if len(catArr) == 0 && idx2 < idx1R && len(rP1)-idx1R < len(p2)-idx2 {
			if overlapLen == 1 {
				if idx2+1 < len(p2) && idx1R-1 >= 0 && p2[idx2+1].ID == rP1[idx1R-1].ID {
					ok = false
				}
			}
			if ok && EqualEdgeFreqArr(rP1[idx1R-idx2:], p2[:idx2+len(rP1)-idx1R]) {
				catArr = make([]EdgeFreq, len(rP1)+len(p2)-(idx2+len(rP1)-idx1R))
				copy(catArr, rP1)
				for i, ef := range p2[:idx2+len(rP1)-idx1R] {
					if catArr[idx1R-idx2+i].Freq+ef.Freq < math.MaxUint8 {
						catArr[idx1R-idx2+i].Freq += ef.Freq
					} else {
						catArr[idx1R-idx2+i].Freq = math.MaxUint8
					}
				}
				copy(catArr[len(rP1):], p2[idx2+len(rP1)-idx1R:])
			}
		}
	}

	if len(catArr) == 0 && len(p2)-idx2 > overlapLen-1 {
		if idx2 > idx1 && len(p1)-idx1 > len(p2)-idx2 {
			if overlapLen == 1 {
				if idx2-1 >= 0 && idx1+1 < len(p1) && p2[idx2-1].ID == p1[idx1+1].ID {
					ok = false
				}
			}
			if ok && EqualEdgeFreqArr(p2[idx2-idx1:], p1[:idx1+len(p2)-idx2]) {
				catArr = make([]EdgeFreq, len(p2)+len(p1)-(idx1+len(p2)-idx2))
				copy(catArr, p2)
				for i, ef := range p1[:idx1+len(p2)-idx2] {
					if catArr[idx2-idx1+i].Freq+ef.Freq < math.MaxUint8 {
						catArr[idx2-idx1+i].Freq += ef.Freq
					} else {
						catArr[idx2-idx1+i].Freq = math.MaxUint8
					}
				}
				copy(catArr[len(p2):], p1[idx1+len(p2)-idx2:])
			}
		}

		if len(catArr) == 0 && idx2 > idx1R && len(rP1)-idx1R > len(p2)-idx2 {
			if overlapLen == 1 {
				if idx2-1 >= 0 && idx1R+1 < len(rP1) && p2[idx2-1].ID == rP1[idx1R+1].ID {
					ok = false
				}
			}
			if ok && EqualEdgeFreqArr(p2[idx2-idx1R:], rP1[:idx1R+len(p2)-idx2]) {
				catArr = make([]EdgeFreq, len(p2)+len(rP1)-(idx1R+len(p2)-idx2))
				copy(catArr, p2)
				for i, ef := range rP1[:idx1R+len(p2)-idx2] {
					if catArr[idx2-idx1R+i].Freq+ef.Freq < math.MaxUint8 {
						catArr[idx2-idx1R+i].Freq += ef.Freq
					} else {
						catArr[idx2-idx1R+i].Freq = math.MaxUint8
					}
				}
				copy(catArr[len(p2):], rP1[idx1R+len(p2)-idx2:])
			}
		}
	}

	return catArr
}

/*func AddPathToEdge2(e *DBGEdge, path, rpath []EdgeFreq) *DBGEdge {
	ok := false
	eID := e.ID
	idx := IndexEdgeFreq(path, eID)
	ridx := len(path) - 1 - idx
	//ridx := IndexEdgeFreq(rpath, eID)
	//fmt.Printf("[AddPathToEdge2]eID:%d, NGSPathArr:%v, path:%v\n", e.ID, e.NGSPathArr, path)
	for i, np := range e.NGSPathArr {
		nidx := IndexEdgeFreq(np, eID)
		min := Min(len(np), len(path))
		if len(np) >= len(path) {
			var p []EdgeFreq
			var x int
			if (idx > 0 && nidx > 0 && path[idx-1].ID == np[nidx-1].ID) || (idx < len(path)-1 && nidx < len(np)-1 && path[idx+1].ID == np[nidx+1].ID) {
				p = path
				x = idx
			} else if (ridx > 0 && nidx > 0 && rpath[ridx-1].ID == np[nidx-1].ID) || (ridx < len(rpath)-1 && nidx < len(np)-1 && rpath[ridx+1].ID == np[nidx+1].ID) {
				p = rpath
				x = ridx
			} else {
				continue
			}
			if nidx >= x && len(np)-nidx >= len(p)-x && EqualEdgeFreqArr(np[nidx-x:nidx-x+min], p) {
				for j, ef := range p {
					if int(np[nidx-x+j].Freq)+int(ef.Freq) < math.MaxUint8 {
						np[nidx-x+j].Freq += ef.Freq
					} else {
						np[nidx-x+j].Freq = math.MaxUint8
					}
				}
				e.NGSPathArr[i] = np
				ok = true
				break
			}
		} else {
			var p []EdgeFreq
			var x int
			if (idx > 0 && nidx > 0 && path[idx-1].ID == np[nidx-1].ID) || (idx < len(path)-1 && nidx < len(np)-1 && path[idx+1].ID == np[nidx+1].ID) {
				p = path
				x = idx
			} else if (ridx > 0 && nidx > 0 && rpath[ridx-1].ID == np[nidx-1].ID) || (ridx < len(rpath)-1 && nidx < len(np)-1 && rpath[ridx+1].ID == np[nidx+1].ID) {
				p = rpath
				x = ridx
			} else {
				continue
			}
			if nidx <= x && len(np)-nidx <= len(p)-x && EqualEdgeFreqArr(np, p[x-nidx:x-nidx+min]) {
				efArr := make([]EdgeFreq, len(p))
				copy(efArr, p)
				for j, ef := range np {
					if int(efArr[x-nidx+j].Freq)+int(ef.Freq) < math.MaxUint8 {
						efArr[x-nidx+j].Freq += ef.Freq
					} else {
						efArr[x-nidx+j].Freq = math.MaxUint8
					}
				}
				e.NGSPathArr[i] = efArr
				ok = true
				break
			}
		}
	}

	if !ok {
		np := make([]EdgeFreq, len(path))
		copy(np, path)
		e.NGSPathArr = append(e.NGSPathArr, np)
	}
	//fmt.Printf("[AddPathToEdge2]eID:%d, e.NGSPathArr:%v\n", e.ID, e.NGSPathArr)

	return e
}*/

func LookUpNearNGSPathArr(arr [][]EdgeFreq, eID uint32) bool {
	for i := range arr {
		if arr[i][1].ID == eID {
			if arr[i][1].Freq < math.MaxUint8 {
				arr[i][0].Freq++
				arr[i][1].Freq++
			}
			return true
		}
	}
	return false
}

func AddFreq(arr []EdgeFreq) {
	for i := range arr {
		if arr[i].Freq < math.MaxUint8 {
			arr[i].Freq++
		}
	}
}
func IsNeedNGSPath(e *DBGEdge, kmerlen int) bool {
	return len(e.Ks) > kmerlen*3/2 || e.GetBubbleFlag() > 0 || e.GetUniqueFlag() > 0
}

func AddPathToEdge(e *DBGEdge, path []uint32, idx, kmerlen int) {
	added := false
	if IsNeedNGSPath(e, kmerlen) {
		for i, p := range e.NGSPathArr {
			if len(path) != len(p) {
				continue
			}
			if EqualEdgeFreqArrPathNoDirection(p, path) {
				AddFreq(e.NGSPathArr[i])
				added = true
				break
			}
		}
		if added {
			return
		}
		// need complex alignment
		for i, p := range e.NGSPathArr {
			if len(path) >= len(p) {
				continue
			}
			if EqualEdgeFreqArrPathNoDirection(p[:len(path)], path) {
				AddFreq(e.NGSPathArr[i][:len(path)])
				added = true
				break
			} else if EqualEdgeFreqArrPathNoDirection(p[len(p)-len(path):], path) {
				AddFreq(e.NGSPathArr[i][len(p)-len(path):])
				added = true
				break
			}
		}
		if added {
			return
		}
		np := make([]EdgeFreq, len(path))
		for j := range path {
			np[j].ID = path[j]
			np[j].Freq = 1
		}
		e.NGSPathArr = append(e.NGSPathArr, np)
	} else { // just need added near edgeID
		if idx > 0 {
			if !LookUpNearNGSPathArr(e.NGSPathArr, path[idx-1]) {
				np := make([]EdgeFreq, 2)
				np[0].ID, np[1].ID = path[idx], path[idx-1]
				np[0].Freq, np[1].Freq = 1, 1
				e.NGSPathArr = append(e.NGSPathArr, np)
			}
		}
		if idx < len(path)-1 {
			if !LookUpNearNGSPathArr(e.NGSPathArr, path[idx+1]) {
				np := make([]EdgeFreq, 2)
				np[0].ID, np[1].ID = path[idx], path[idx+1]
				np[0].Freq, np[1].Freq = 1, 1
				e.NGSPathArr = append(e.NGSPathArr, np)
			}
		}
	}
	//fmt.Printf("[AddPathToEdge2]eID:%d, e.NGSPathArr:%v\n", e.ID, e.NGSPathArr)
}

const MaxPathSize = 9

func paraAddNGSPathToDBG(cs <-chan []uint32, arrPool *sync.Pool, edgesArr []DBGEdge, opt optionsMN) {
	var path []uint32
	var e *DBGEdge
	for {
		arr, ok := <-cs
		if !ok {
			break
		}

		for i := 0; i < len(arr)-1; {
			if arr[i] == 0 {
				break
			}
			idx := IndexNextZeroUint32Arr(arr[i:], FORWARD)
			if idx < 0 {
				log.Fatalf("[AddNGSPathToDBG]idx < 0\n")
			}
			path = arr[i : i+idx]
			for j, id := range path {
				e = &edgesArr[id]
				if len(e.NGSPathArr) >= MaxPathSize {
					continue
				}
				ParaSetProcessFlag(&e.Flag)
				AddPathToEdge(e, path, j, opt.Kmer)
				ParaResetProcessFlag(&e.Flag)
			}
			i += idx + 1
		}
		arrPool.Put(arr)
	}
}

// merge contained NGSPath
func MergeNGSPathArr(edgesArr []DBGEdge) {
	maxCount, mergeNum := 0, 0
	for i := range edgesArr {
		e := &edgesArr[i]
		if i < 2 || e.GetDeleteFlag() > 0 {
			continue
		}
		if len(e.NGSPathArr) == MaxPathSize {
			//e.NGSPathArr = nil
			//fmt.Printf("[MergeNGSPathArr]max eID:%d el:%d NGSPathArr:%v\n", e.ID, len(e.Ks), e.NGSPathArr)
			maxCount++
			continue
		}
		if len(e.NGSPathArr) <= 1 {
			continue
			//fmt.Fprintf(os.Stderr, "[AddNGSPathToDBG]e.NGSPathArr:%v is nil\n", e.NGSPathArr)
		}
		for j := len(e.NGSPathArr) - 1; j > 0; j-- {
			if e.NGSPathArr[j] == nil {
				continue
			}
			for k := 0; k < j; k++ {
				if e.NGSPathArr[k] == nil {
					continue
				}
				if len(e.NGSPathArr[j]) < len(e.NGSPathArr[k]) {
					x := len(e.NGSPathArr[k]) - len(e.NGSPathArr[j])
					for l := 0; l <= x; l++ {
						jl := len(e.NGSPathArr[j])
						if EqualEdgeFreqArrNoDirection(e.NGSPathArr[k][l:l+jl], e.NGSPathArr[j]) {
							AddEdgeFreqArr(e.NGSPathArr[k][l:l+jl], e.NGSPathArr[j])
							e.NGSPathArr[j] = nil
							break
						}
					}
				} else {
					x := len(e.NGSPathArr[j]) - len(e.NGSPathArr[k])
					for l := 0; l <= x; l++ {
						kl := len(e.NGSPathArr[k])
						if EqualEdgeFreqArrNoDirection(e.NGSPathArr[j][l:l+kl], e.NGSPathArr[k]) {
							AddEdgeFreqArr(e.NGSPathArr[j][l:l+kl], e.NGSPathArr[j])
							e.NGSPathArr[j] = nil
							break
						}
					}
				}
			}
		}
		ngsi := 0
		for j := range e.NGSPathArr {
			if e.NGSPathArr[j] != nil {
				e.NGSPathArr[ngsi] = e.NGSPathArr[j]
				ngsi++
			} else {
				mergeNum++
			}
		}
		e.NGSPathArr = e.NGSPathArr[:ngsi]
		//fmt.Printf("[MergeNGSPathArr]eID:%d el:%d NGSPathArr:%v\n", e.ID, len(e.Ks), e.NGSPathArr)
	}
	fmt.Printf("[MergeNGSPathArr]maxCount:%d mergeNum:%d\n", maxCount, mergeNum)
}

func DeleteLowFreqEdge(edgesArr []DBGEdge, nodesArr []DBGNode) {
	deleteEdgeNum := 0
	for i := range edgesArr {
		e := &edgesArr[i]
		if e.ID < 2 || e.GetDeleteFlag() > 0 {
			continue
		}
		if len(e.NGSPathArr) == MaxPathSize {
			e.NGSPathArr = nil
			continue
		}
		if len(e.NGSPathArr) == 0 {
			fmt.Printf("[DeleteLowFreqEdge]eID:%d el:%d Flag:%b\n", e.ID, len(e.Ks), e.Flag)
			edgesArr[i].SetDeleteFlag()
			deleteEdgeNum++
			if e.StartNID > 0 {
				if !SubstituteEdgeID(nodesArr, e.StartNID, e.ID, 0) {
					fmt.Fprintf(os.Stderr, "[DeleteLowFreqEdge]v: %v\ne.ID:%d substitute by 0 failed\n", nodesArr[e.StartNID], e.ID)
				}
			}
			if e.EndNID > 0 {
				if !SubstituteEdgeID(nodesArr, e.EndNID, e.ID, 0) {
					fmt.Fprintf(os.Stderr, "[DeleteLowFreqEdge]v:%v\ne.ID:%d substitute by 0 failed\n", nodesArr[e.EndNID], e.ID)
				}
			}
		}
	}
	fmt.Printf("[DeleteLowFreqEdge]deleteEdgeNum:%d\n", deleteEdgeNum)
}

// clean low freq path
func DeleteLowFreqNGSPath(edgesArr []DBGEdge, MinMapFreq int) {
	deleteNum := 0
	for i := range edgesArr {
		e := &edgesArr[i]
		if e.ID < 2 || e.GetDeleteFlag() > 0 {
			continue
		}
		if len(e.NGSPathArr) == 0 {
			continue
		}

		idx := 0
		for j, path := range e.NGSPathArr {
			idx1, idx2 := -1, -1
			for k, ef := range path {
				if int(ef.Freq) >= MinMapFreq {
					idx1 = k
					break
				}
			}
			for k := len(path) - 1; k >= 0; k-- {
				if int(path[k].Freq) >= MinMapFreq {
					idx2 = k
					break
				}
			}

			if idx2+1-idx1 >= 3 {
				e.NGSPathArr[idx] = e.NGSPathArr[j][idx1 : idx2+1]
				idx++
			} else {
				deleteNum++
			}
		}
		e.NGSPathArr = e.NGSPathArr[:idx]
	}
	fmt.Printf("[DeleteLowFreqNGSPath]deleteNum:%d\n", deleteNum)

}

// cat left and rigth edge path or merge consistence path
func CatLeftRigthPath(edgesArr []DBGEdge) {
	catNum := 0
	for i := range edgesArr {
		e := &edgesArr[i]
		if e.ID < 2 || e.GetDeleteFlag() > 0 || len(e.NGSPathArr) != 2 {
			continue
		}
		ok := true
		for j, p := range e.NGSPathArr {
			if p[0].ID != e.ID {
				if p[len(p)-1].ID != e.ID {
					ok = false
					break
				} else {
					e.NGSPathArr[j] = ReverseEdgeFreqArr(e.NGSPathArr[j])
				}
			}
		}
		if !ok {
			continue
		}
		p0, p1 := e.NGSPathArr[0], e.NGSPathArr[1]
		if e.GetTwoEdgesCycleFlag() > 0 {
			if p0[1].ID == p1[1].ID {
				e1 := &edgesArr[p0[1].ID]
				ea1 := GetComingEdgeArr(e1.EdgeIDIncoming)
				ea2 := GetComingEdgeArr(e1.EdgeIDOutcoming)
				if len(ea1) == 2 && len(ea2) == 2 && ((IsInuint32Arr(ea1, p0[2].ID) && IsInuint32Arr(ea2, p1[2].ID)) || (IsInuint32Arr(ea1, p1[2].ID) && IsInuint32Arr(ea2, p0[2].ID))) {

				} else {
					ok = false
				}
			}
		} else {
			if p0[1].ID == e.ID || p1[1].ID == e.ID {
				ok = false
			} else {
				if (IsInuint32Slice(e.EdgeIDIncoming, p0[1].ID) && IsInuint32Slice(e.EdgeIDOutcoming, p1[1].ID)) || (IsInuint32Slice(e.EdgeIDIncoming, p1[1].ID) && IsInuint32Slice(e.EdgeIDOutcoming, p0[1].ID)) {

				} else {
					ok = false
				}
			}
		}
		if ok {
			ea := make([]EdgeFreq, 0, len(p0)+len(p1)-1)
			for k := len(p0) - 1; k >= 0; k-- {
				ea = append(ea, p0[k])
			}
			for k := 1; k < len(p1); k++ {
				ea = append(ea, p1[k])
			}
			//fmt.Printf("[CatLeftRigthPath]p0:%v p1:%v catArr:%v\n", p0, p1, ea)
			e.NGSPathArr[0] = ea
			e.NGSPathArr = e.NGSPathArr[:1]
			catNum++
		}
	}
	fmt.Printf("[CatLeftRigthPath]catNum:%d\n", catNum)
}

func CountEdgeFreqArr(path []EdgeFreq, eID uint32) (count int) {
	for _, ef := range path {
		if ef.ID == eID {
			count++
		}
	}
	return
}

// merge Unique path
func MergeUniquePath(edgesArr []DBGEdge, nodesArr []DBGNode, kmerlen int) {
	mergeNum := 0
	for i := range edgesArr {
		e := &edgesArr[i]
		if len(e.NGSPathArr) > 0 {
			fmt.Fprintf(os.Stderr, "[MergeUniquePath]eID:%d NGSPathArr:%v\n", e.ID, e.NGSPathArr)
		}
		if e.ID < 2 || e.GetDeleteFlag() > 0 || len(e.NGSPathArr) != 1 || e.GetProcessFlag() > 0 {
			continue
		}
		if e.GetUniqueFlag() == 0 && e.GetBubbleFlag() == 0 {
			continue
		}
		ep := e.NGSPathArr[0]
		// check cycle path
		if CountEdgeFreqArr(ep, ep[0].ID) > 1 || CountEdgeFreqArr(ep, ep[len(ep)-1].ID) > 1 {
			continue
		}
		efArr, eIDArr := GetExtendUniquePath(e, edgesArr)
		if len(eIDArr) < 2 || edgesArr[eIDArr[0]].GetTwoEdgesCycleFlag() > 0 {
			continue
		}
		//fmt.Printf("[MergeUniquePath]efArr:%v UniqueArr:%v\n", efArr, eIDArr)
		// check cycle path
		{
			cycleOk := false
			for _, eID := range eIDArr {
				if CountEdgeFreqArr(efArr, eID) > 1 {
					fmt.Fprintf(os.Stderr, "[MergeUniquePath]found cycle path eID:%d in path:%v\n", eID, efArr)
					cycleOk = true
					break
				}
			}
			if cycleOk {
				continue
			}
		}
		// set process flag
		for _, eID := range eIDArr {
			edgesArr[eID].SetProcessFlag()
		}
		mergeNum++
		idx1, idx2 := IndexEdgeFreq(efArr, eIDArr[0]), IndexEdgeFreq(efArr, eIDArr[len(eIDArr)-1])

		e, e1 := &edgesArr[efArr[idx1].ID], &edgesArr[efArr[idx1+1].ID]
		nID := GetShareNodeID(e, e1)
		te := ConcatEdgePathUtg(efArr[idx1:idx2+1], edgesArr, nodesArr, kmerlen)
		// reset DBG
		if !SubstituteEdgeID(nodesArr, nID, e.ID, 0) {
			fmt.Printf("[MergeUniquePath]e.StartNID:%d e.EndNID:%d e1.StartNID:%d e1.EndNID:%d\n", e.StartNID, e.EndNID, e1.StartNID, e1.EndNID)
			log.Fatalf("[MergeUniquePath]nID:%d v: %v\ne.ID:%d substitute by 0 failed\n", nID, nodesArr[nID], e.ID)
		}
		if te.EndNID > 0 {
			if !SubstituteEdgeID(nodesArr, te.EndNID, eIDArr[len(eIDArr)-1], e.ID) {
				log.Fatalf("[MergeUniquePath]v: %v\ne.ID:%d substitute by :%d failed\n", nodesArr[te.EndNID], eIDArr[len(eIDArr)-1], e.ID)
			}
		}
		for _, eID := range eIDArr[1:] {
			ne := &edgesArr[eID]
			if ne.StartNID > 0 && ne.StartNID != te.EndNID && !SubstituteEdgeID(nodesArr, ne.StartNID, ne.ID, 0) {
				log.Fatalf("[MergeUniquePath]v: %v\ne.ID:%d substitute by 0 failed\n", nodesArr[ne.StartNID], ne.ID)
			}
			if ne.EndNID > 0 && ne.EndNID != te.EndNID && !SubstituteEdgeID(nodesArr, ne.EndNID, ne.ID, 0) {
				log.Fatalf("[MergeUniquePath]v: %v\ne.ID:%d substitute by 0 failed\n", nodesArr[ne.EndNID], ne.ID)
			}
			edgesArr[eID].SetDeleteFlag()
		}
		edgesArr[te.ID] = te

		// reset NGSPathArr
		{
			eID1, eID2 := eIDArr[0], eIDArr[1]
			for _, ef := range efArr[idx1+1 : idx2] {
				//eID := ef.ID
				idx := IndexUint32(eIDArr, ef.ID)
				if idx >= 0 && idx+1 < len(eIDArr) {
					eID1 = eID2
					eID2 = eIDArr[idx+1]
					continue
				}
				pa := edgesArr[ef.ID].NGSPathArr
				k := 0
				for _, ep := range pa {
					idx5 := IndexEdgeFreq(ep, eID1)
					idx6 := IndexEdgeFreq(ep, eID2)
					if idx5 < 0 && idx6 < 0 {
						pa[k] = ep
						k++
					}
				}
				edgesArr[ef.ID].NGSPathArr = pa[:k]
			}
		}

		fmt.Printf("[MergeUniquePath]merged path:%v UniqueArr:%v\n", efArr[idx1:idx2+1], eIDArr)
		mergePath := efArr[idx1 : idx2+1]
		revMP := GetReverseEdgeFreqArr2(mergePath)
		for _, ef := range efArr {
			id := ef.ID
			if IsInuint32Arr(eIDArr[1:], id) {
				continue
			}
			for j, tmp := range edgesArr[id].NGSPathArr {
				edgesArr[id].NGSPathArr[j] = RePlaceNGSPath2(tmp, mergePath, revMP)
				//if tl != len(edgesArr[id].NGSPathArr[j]) {
				//	fmt.Printf("[AddNGSPathToDBG]eID:%d j:%d path:%v\n\t\tchanged:%v\n", id, j, td, edgesArr[id].NGSPathArr[j])
				//}
			}
		}
	}
	fmt.Printf("[MergeUniquePath]mergeNum:%d\n", mergeNum)
}

func AddNGSPathToDBG(NGSPathfn string, nodesArr []DBGNode, edgesArr []DBGEdge, opt optionsMN, ngsPathSize int) {
	bufSize := WindowSize
	cs := make(chan []uint32, opt.NumCPU)
	arrPool := sync.Pool{New: func() interface{} {
		arr := make([]uint32, bufSize/4)
		return arr
	}}
	fmt.Printf("[AddNGSPathToDBG] begin processe file:%s...\n", NGSPathfn)
	for i := 0; i < opt.NumCPU-1; i++ {
		go paraAddNGSPathToDBG(cs, &arrPool, edgesArr, opt)
	}
	ReadZstNGSPathFile(NGSPathfn, &arrPool, cs, ngsPathSize)
	time.Sleep(time.Millisecond)

	// merge contained NGSPath
	MergeNGSPathArr(edgesArr)
	// remove low freq edge
	//DeleteLowFreqEdge(edgesArr, nodesArr)
	DeleteLowFreqNGSPath(edgesArr, opt.MinMapFreq)
	// merge contained NGSPath
	MergeNGSPathArr(edgesArr)
	// cat left and rigth edge path or merge consistence path
	CatLeftRigthPath(edgesArr)
	// merge Unique path
	MergeUniquePath(edgesArr, nodesArr, opt.Kmer)
}

/*
func AddNGSPathToDBG2(NGSPathfn string, nodesArr []DBGNode, edgesArr []DBGEdge, opt optionsMN) {
bufSize := 10000
minAvgEdgeFreq := opt.MinMapFreq
cs := make(chan []uint32, bufSize)
go LoadNGSPathFromFile(NGSPathfn, cs)

MaxPathSize := 16
emptyChannelNum, fullChannelNum := 0, 0
//AllowMaxPathNum := 8
//MinNGSPathNum := 1
pArr := make([]EdgeFreq, 10)
rPath := make([]EdgeFreq, 10)
processedNum := 0
for {
	eIDArr, ok := <-cs
	if !ok {
		break
	}
	if len(cs) == 0 {
		emptyChannelNum++
	} else if len(cs) == cap(cs) {
		fullChannelNum++
	}
	if len(eIDArr) < 3 {
		continue
	}
	pArr = Transform2EdgeFreqArr(eIDArr, pArr)
	rPath = GetReverseEdgeFreqArr(pArr, rPath)
	processedNum++
	for _, ef := range pArr {
		e := &edgesArr[ef.ID]
		if len(e.NGSPathArr) >= MaxPathSize {
			continue
		}
		//fmt.Printf("[AddNGSPathToDBG]eID:%d el:%d NGSPathArr:%v\n", eID, e.GetSeqLen(), e.NGSPathArr)
		//AddPathToEdge(edgesArr, nodesArr, e, path, rArr, MaxPathSize)
		AddPathToEdge2(e, pArr, rPath)
		//fmt.Printf("after addPath:%v, eID:%d el:%d NGSPathArr:%v\n", path, eID, e.GetSeqLen(), edgesArr[eID].NGSPathArr)
	}
}
//fmt.Printf("[AddNGSPathToDBG]AddPathToEdge2 eID:299779 NGSPathArr:%v\n", edgesArr[299779].NGSPathArr)
fmt.Printf("[AddNGSPathToDBG] processed path:%d emptyChannelNum:%d fullChannelNum:%d\n", processedNum, emptyChannelNum, fullChannelNum)
bigMaxPathSizeNum := 0
repeatSinglePathNum := 0
singlePathNum := 0
noPathNum := 0
manyPathNum := 0
cyclePathNum := 0
// clean low freq path
for i, e := range edgesArr {
	if e.ID < 2 || e.GetDeleteFlag() > 0 {
		continue
	}
	if len(e.NGSPathArr) == 0 {
		noPathNum++
		continue
	} else if len(e.NGSPathArr) >= MaxPathSize {
		bigMaxPathSizeNum++
		if e.GetUniqueFlag() > 0 {
			edgesArr[i].ResetUniqueFlag()
		}
		edgesArr[i].NGSPathArr = nil
		continue
	}

	pathArr := e.NGSPathArr
	idx := -1
	// process MaxPathSize case
	//fmt.Printf("[AddNGSPathToDBG]eID:%d, NGSPathArr:%v\n", i, e.NGSPathArr)

	for _, path := range e.NGSPathArr {
		idx1, idx2 := -1, -1
		for j, ef := range path {
			if int(ef.Freq) >= minAvgEdgeFreq {
				idx1 = j
				break
			}
		}
		for j := len(path) - 1; j >= 0; j-- {
			if int(path[j].Freq) >= minAvgEdgeFreq {
				idx2 = j
				break
			}
		}
		//if idx2+1-idx1 >= 3 && GetAvgEdgeFreq(path[idx1:idx2+1]) >= minAvgEdgeFreq {
		if idx2+1-idx1 >= 3 {
			if CountEdgeFreqArr(path, e.ID) > 1 {
				cyclePathNum++
				idx = -1
				break
			}
			//fmt.Printf("[AddNGSPathToDBG]idx1:%d idx2:%d len(path):%d\n", idx1, idx2, len(path))
			// check conatined element
			contained := false
			for j := 0; j <= idx; j++ {
				p0 := pathArr[j]
				p1 := path[idx1 : idx2+1]
				if len(p0) < len(p1) {
					p0, p1 = p1, p0
				}
				x1, x2 := IndexEdgeFreq(p0, p1[0].ID), IndexEdgeFreq(p0, p1[len(p1)-1].ID)
				if x1 >= 0 && x2 >= 0 {
					if x1 > x2 {
						p1 = ReverseEdgeFreqArr(p1)
						x1, x2 = x2, x1
					}
					if EqualEdgeFreqArr(p0[x1:x2+1], p1) {
						for h, ef := range p1 {
							if int(p0[x1+h].Freq)+int(ef.Freq) < math.MaxUint8 {
								p0[x1+h].Freq += ef.Freq
							} else {
								p0[x1+h].Freq = math.MaxUint8
							}
						}
						pathArr[j] = p0
						contained = true
						//fmt.Printf("[AddNGSPathToDBG]contained p0:%v p1:%v\n", p0, p1)
						break
					}
				}
			}
			if !contained {
				idx++
				pathArr[idx] = path[idx1 : idx2+1]
			}
		}
	}
	edgesArr[i].NGSPathArr = pathArr[:idx+1]
}

// cat left and rigth edge path or merge consistence path

catNum := 0
//rPath := make([]EdgeFreq, 20)
for i := range edgesArr {
	e := &edgesArr[i]
	if e.ID < 2 || e.GetDeleteFlag() > 0 || len(e.NGSPathArr) < 2 {
		continue
	}

	overlapLen := 2
	if e.GetTwoEdgesCycleFlag() > 0 {
		overlapLen = 3
	} else if e.GetUniqueFlag() > 0 || e.GetBubbleFlag() > 0 {
		overlapLen = 1
	}

	//fmt.Printf("[AddNGSPathToDBG]eID:%d overlapLen:%d el:%d\n", i, overlapLen, e.GetSeqLen())
	//for x, p := range e.NGSPathArr {
	//	fmt.Printf("Arr[%d]:%v\n", x, p)
	//}

	deleteFlagArr := make([]bool, len(e.NGSPathArr))
	for x, px := range e.NGSPathArr {
		if deleteFlagArr[x] {
			continue
		}
		count := 0
		overIdx := -1
		var ca []EdgeFreq
		ix := IndexEdgeFreq(px, e.ID)
		rPath := GetReverseEdgeFreqArr(px, rPath)
		ixR := len(px) - 1 - ix
		for y, py := range e.NGSPathArr {
			if x == y || deleteFlagArr[y] {
				continue
			}
			iy := IndexEdgeFreq(py, e.ID)
			catArr := GetOverlapPath(px, rPath, py, ix, ixR, iy, overlapLen)
			if len(catArr) > 0 {
				count++
				if count > 1 {
					break
				}
				overIdx = y
				ca = catArr
			}
		}
		if count == 1 {
			e.NGSPathArr[x] = ca
			deleteFlagArr[overIdx] = true
		}
	}

	idx := 0
	for x, p := range e.NGSPathArr {
		if deleteFlagArr[x] {
			continue
		}
		e.NGSPathArr[idx] = p
		idx++
	}
	if idx < len(e.NGSPathArr) {
		catNum += len(e.NGSPathArr) - idx
		//for x, p := range e.NGSPathArr[:idx] {
		//	fmt.Printf("catArr[%d]:%v\n", x, p)
		//}
	}
	e.NGSPathArr = e.NGSPathArr[:idx]
}

fmt.Printf("[AddNGSPathToDBG]catNum:%d bigMaxPathSizeNum:%d, singlePathNum:%d,repeatSinglePathNum:%d, manyPathNum:%d, noPathNum:%d cyclePathNum:%d\n", catNum, bigMaxPathSizeNum, singlePathNum, repeatSinglePathNum, manyPathNum, noPathNum, cyclePathNum)

// merge Unique path
mergeNum := 0
for i := range edgesArr {
	e := &edgesArr[i]
	if e.ID < 2 || e.GetDeleteFlag() > 0 || len(e.NGSPathArr) != 1 || e.GetProcessFlag() > 0 {
		continue
	}
	if e.GetUniqueFlag() == 0 && e.GetTwoEdgesCycleFlag() == 0 && e.GetBubbleFlag() == 0 {
		continue
	}

	ep := e.NGSPathArr[0]
	avgFreq := GetAvgEdgeFreq(ep)
	if avgFreq < minAvgEdgeFreq*2 {
		continue
	}
	// check cycle path
	if CountEdgeFreqArr(ep, ep[0].ID) > 1 || CountEdgeFreqArr(ep, ep[len(ep)-1].ID) > 1 {
		continue
	}
	//ep := e.NGSPathArr[0]
	efArr, eIDArr := GetExtendUniquePath(e, edgesArr, minAvgEdgeFreq)
	if len(eIDArr) < 2 {
		continue
	}
	if edgesArr[eIDArr[0]].GetTwoEdgesCycleFlag() > 0 || edgesArr[eIDArr[len(eIDArr)-1]].GetTwoEdgesCycleFlag() > 0 {
		continue
	}
	fmt.Printf("[AddNGSPathToDBG]efArr:%v UniqueArr:%v\n", efArr, eIDArr)
	//if efArr[idx1].Freq < minAvgEdgeFreq*2 || efArr[idx1]
	// check cycle path
	{
		cycleOk := false
		for _, eID := range eIDArr {
			if CountEdgeFreqArr(efArr, eID) > 1 {
				fmt.Fprintf(os.Stderr, "[AddNGSPathToDBG]found cycle path eID:%d in path:%v\n", eID, efArr)
				cycleOk = true
				break
			}
		}
		if cycleOk {
			continue
		}
	}
	// set process flag
	for _, eID := range eIDArr {
		edgesArr[eID].SetProcessFlag()
	}
	mergeNum++
	idx1, idx2 := IndexEdgeFreq(efArr, eIDArr[0]), IndexEdgeFreq(efArr, eIDArr[len(eIDArr)-1])

	e, e1 := &edgesArr[efArr[idx1].ID], &edgesArr[efArr[idx1+1].ID]
	nID := GetShareNodeID(e, e1)
	te := ConcatEdgePathUtg(efArr[idx1:idx2+1], nID, edgesArr, nodesArr, opt.Kmer)
	// reset DBG

	if !SubstituteEdgeID(nodesArr, nID, e.ID, 0) {
		fmt.Printf("[AddNGSPathToDBG]e.StartNID:%d e.EndNID:%d e1.StartNID:%d e1.EndNID:%d\n", e.StartNID, e.EndNID, e1.StartNID, e1.EndNID)
		log.Fatalf("[AddNGSPathToDBG]nID:%d v: %v\ne.ID:%d substitute by 0 failed\n", nID, nodesArr[nID], e.ID)
	}

	if te.EndNID > 0 {
		if !SubstituteEdgeID(nodesArr, te.EndNID, eIDArr[len(eIDArr)-1], e.ID) {
			log.Fatalf("[AddNGSPathToDBG]v: %v\ne.ID:%d substitute by :%d failed\n", nodesArr[te.EndNID], eIDArr[len(eIDArr)-1], e.ID)
		}
	}

	for _, eID := range eIDArr[1:] {
		ne := &edgesArr[eID]
		if ne.StartNID > 0 && ne.StartNID != te.EndNID && !SubstituteEdgeID(nodesArr, ne.StartNID, ne.ID, 0) {
			log.Fatalf("[AddNGSPathToDBG]v: %v\ne.ID:%d substitute by 0 failed\n", nodesArr[ne.StartNID], ne.ID)
		}
		if ne.EndNID > 0 && ne.EndNID != te.EndNID && !SubstituteEdgeID(nodesArr, ne.EndNID, ne.ID, 0) {
			log.Fatalf("[AddNGSPathToDBG]v: %v\ne.ID:%d substitute by 0 failed\n", nodesArr[ne.EndNID], ne.ID)
		}
		edgesArr[eID].SetDeleteFlag()
	}
	edgesArr[te.ID] = te
	//nPath := make([]EdgeFreq, idx1+1+len(efArr)-1-idx2)
	//copy(nPath[:idx1+1], efArr[:idx1+1])
	//copy(nPath[idx1+1:], efArr[idx2:])
	//edgesArr[eID].NGSPathArr[0] = nPath

	// reset EdgeID if path
	{

		eID1, eID2 := eIDArr[0], eIDArr[1]
		for _, ef := range efArr[idx1+1 : idx2] {
			//eID := ef.ID
			idx := IndexUint32(eIDArr, ef.ID)
			if idx >= 0 && idx+1 < len(eIDArr) {
				eID1 = eID2
				eID2 = eIDArr[idx+1]
				continue
			}
			pa := edgesArr[ef.ID].NGSPathArr
			nArr := pa[:0]
			for _, ep := range pa {
				idx5 := IndexEdgeFreq(ep, eID1)
				idx6 := IndexEdgeFreq(ep, eID2)
				if idx5 < 0 && idx6 < 0 {
					nArr = append(nArr, ep)
				}
			}
			edgesArr[ef.ID].NGSPathArr = nArr
		}
	}

	fmt.Printf("[AddNGSPathToDBG]merged path:%v UniqueArr:%v\n", efArr[idx1:idx2+1], eIDArr)
	mergePath := efArr[idx1 : idx2+1]
	revMP := GetReverseEdgeFreqArr2(mergePath)
	for _, ef := range efArr {
		id := ef.ID
		if IsInuint32Arr(eIDArr[1:], id) {
			continue
		}
		for j, tmp := range edgesArr[id].NGSPathArr {
			edgesArr[id].NGSPathArr[j] = RePlaceNGSPath2(tmp, mergePath, revMP)
			//if tl != len(edgesArr[id].NGSPathArr[j]) {
			//	fmt.Printf("[AddNGSPathToDBG]eID:%d j:%d path:%v\n\t\tchanged:%v\n", id, j, td, edgesArr[id].NGSPathArr[j])
			//}
		}
	}

}
fmt.Printf("[AddNGSPathToDBG]mergeNum:%d\n", mergeNum)

/*
	// merge path by edge
	tooManyPathNum := 0
	noPathNum := 0
	morePathNum := 0
	totalPathNum := 0
	repeatNum, otherNum := 0, 0
	for i, e := range edgesArr {
		if e.ID < 2 || e.GetDeleteFlag() > 0 {
			continue
		}
		eID := e.ID
		pa := e.NGSPathArr
		if len(pa) <= 2 {
			cn := GetNumFreqInPathFreqArr(pa, eID)
			min := 0
			if e.StartNID > 0 {
				min += MinNGSPathNum
			}
			if e.EndNID > 0 {
				min += MinNGSPathNum
			}
			if min == 2 {
				min = 1
			}
			if cn <= min && e.GetSeqLen() < opt.Kmer*2+20 {
				edgesArr[i].SetDeleteFlag()
				noPathNum++
				if e.GetUniqueFlag() == 0 && e.GetBubbleFlag() == 0 && e.GetBubbleRepeatFlag() == 0 {
					repeatNum++
				} else {
					otherNum++
				}
				fmt.Fprintf(os.Stderr, "[AddNGSPathToDBG] eID:%d el:%d no NGSPath delete\n", eID, e.GetSeqLen())
				edgesArr[eID].NGSPathArr = nil
				continue
			}
			if len(pa) == 0 {
				continue
			}
		}
		if len(pa) > MaxPathSize {
			//fmt.Printf("[AddNGSPathToDBG]el:%d len(pa):%d pa:%v\n", e.GetSeqLen(), len(pa), pa)
			tooManyPathNum++
			edgesArr[eID].NGSPathArr = nil
			continue
		}


		fmt.Printf("[AddNGSPathToDBG]eID:%d el:%d pa:%v\n", eID, e.GetSeqLen(), pa)

		if len(pa) == 1 {
			npf := CutLowFreq(pa[0], uint8(opt.MinMapFreq))
			if len(npf.IDArr) < 3 {
				edgesArr[eID].NGSPathArr = nil
			} else {
				edgesArr[eID].NGSPathArr[0] = npf
			}
			fmt.Printf("[AddNGSPathToDBG]eID:%d el:%d Unique:%v|Bubble:%v NGSPathArr:%v\n", eID, e.GetSeqLen(), e.GetUniqueFlag() > 0, e.GetBubbleFlag() > 0, edgesArr[eID].NGSPathArr)
			continue
		}

		leftMax, rightMax := 0, 0
		for _, pf := range pa {
			idx := IndexUint32(pf.IDArr, eID)
			if idx > leftMax {
				leftMax = idx
			}
			if len(pf.IDArr)-idx > rightMax {
				rightMax = len(pf.IDArr) - idx
			}
		}
		al := leftMax + rightMax
		pArr := make([]PathFreqExt, len(pa))
		for j, pf := range pa {
			var pfe PathFreqExt
			pfe.PF.IDArr = make([]uint32, al)
			pfe.PF.FreqArr = make([]uint8, al)
			idx := IndexUint32(pf.IDArr, eID)
			copy(pfe.PF.IDArr[leftMax-idx:], pf.IDArr)
			copy(pfe.PF.FreqArr[leftMax-idx:], pf.FreqArr)
			pfe.Start = leftMax - idx
			pfe.End = pfe.Start + len(pf.IDArr)
			pArr[j] = pfe
		}

		stk := make([]PathFreqExt, 0, 10)
		var pfe PathFreqExt
		pfe.PF.IDArr = make([]uint32, al)
		pfe.PF.FreqArr = make([]uint8, al)
		pfe.PF.IDArr[leftMax] = eID
		pfe.PF.FreqArr[leftMax] = 20
		pfe.Start = leftMax - 1
		pfe.End = leftMax + 1
		stk = append(stk, pfe)

		rstk := make([]PathFreqExt, 0, 10)

		for len(stk) > 0 {
			p1 := stk[len(stk)-1]
			stk = stk[:len(stk)-1]
			j := p1.Start
			for ; j >= 0; j-- {
				var idFreqArr []IDFreq
				for _, pfe := range pArr {
					if pfe.Start > j {
						continue
					}
					if Equaluint32Arr(p1.PF.IDArr[j+1:p1.End], pfe.PF.IDArr[j+1:p1.End]) {
						idFreqArr = AddIDFreqArr(idFreqArr, pfe.PF.IDArr[j], int(pfe.PF.FreqArr[j]))
					}
				}
				idFreqArr = DeleteLowFreq(idFreqArr, opt.MinMapFreq)
				if len(idFreqArr) == 0 {
					p1.Start = j
					rstk = append(rstk, p1)
					break
				} else if len(idFreqArr) == 1 {
					p1.PF.IDArr[j] = idFreqArr[0].ID
					if idFreqArr[0].Freq > math.MaxUint8 {
						p1.PF.FreqArr[j] = math.MaxUint8
					} else {
						p1.PF.FreqArr[j] = uint8(idFreqArr[0].Freq)
					}
				} else {
					for x, idfreq := range idFreqArr {
						var npf PathFreqExt
						if x == len(idFreqArr)-1 {
							npf = p1
						} else {
							npf.PF.IDArr = make([]uint32, al)
							npf.PF.FreqArr = make([]uint8, al)
							copy(npf.PF.IDArr, p1.PF.IDArr)
							copy(npf.PF.FreqArr, p1.PF.FreqArr)
						}
						npf.PF.IDArr[j] = idfreq.ID
						if idfreq.Freq > math.MaxUint8 {
							npf.PF.FreqArr[j] = math.MaxUint8
						} else {
							npf.PF.FreqArr[j] = uint8(idfreq.Freq)
						}
						npf.Start = j - 1
						npf.End = p1.End
						stk = append(stk, npf)
					}
					break
				}
			}
			if j < 0 {
				p1.Start = j
				rstk = append(rstk, p1)
			}
		}

		var startFlag bool
		if len(rstk) == 0 {
			startFlag = true
		}

		var rpfe PathFreqExt
		rpfe.PF.IDArr = make([]uint32, al)
		rpfe.PF.FreqArr = make([]uint8, al)
		rpfe.PF.IDArr[leftMax] = eID
		rpfe.PF.FreqArr[leftMax] = 20
		rpfe.Start = leftMax - 1
		rpfe.End = leftMax + 1
		rstk = append(rstk, rpfe)

		// found right part path
		for len(rstk) > 0 {
			p1 := rstk[len(rstk)-1]
			rstk = rstk[:len(rstk)-1]
			j := p1.End
			for ; j < al; j++ {
				var idFreqArr []IDFreq
				for _, pfe := range pArr {
					if pfe.End <= j {
						continue
					}
					var start int
					if startFlag {
						start = MaxInt(p1.Start, pfe.Start)
					} else {
						if p1.Start == leftMax-1 {
							start = p1.Start
						} else {
							start = MinInt(p1.Start, pfe.Start)
						}
					}
					if Equaluint32Arr(p1.PF.IDArr[start+1:p1.End], pfe.PF.IDArr[start+1:p1.End]) {
						idFreqArr = AddIDFreqArr(idFreqArr, pfe.PF.IDArr[j], int(pfe.PF.FreqArr[j]))
					}
				}
				idFreqArr = DeleteLowFreq(idFreqArr, opt.MinMapFreq)
				//fmt.Printf("[AddNGSPathToDBG]idFreqArr:%v p1:%v\n", idFreqArr, p1)
				if len(idFreqArr) == 0 {
					p1.End = j
					stk = append(stk, p1)
					break
				} else if len(idFreqArr) == 1 {
					p1.PF.IDArr[j] = idFreqArr[0].ID
					if idFreqArr[0].Freq > math.MaxUint8 {
						p1.PF.FreqArr[j] = math.MaxUint8
					} else {
						p1.PF.FreqArr[j] = uint8(idFreqArr[0].Freq)
					}
				} else {
					for x, idfreq := range idFreqArr {
						var npf PathFreqExt
						if x == len(idFreqArr)-1 {
							npf = p1
						} else {
							npf.PF.IDArr = make([]uint32, al)
							npf.PF.FreqArr = make([]uint8, al)
							copy(npf.PF.IDArr, p1.PF.IDArr)
							copy(npf.PF.FreqArr, p1.PF.FreqArr)
						}
						npf.PF.IDArr[j] = idfreq.ID
						if idfreq.Freq > math.MaxUint8 {
							npf.PF.FreqArr[j] = math.MaxUint8
						} else {
							npf.PF.FreqArr[j] = uint8(idfreq.Freq)
						}
						npf.Start = p1.Start
						npf.End = j + 1
						rstk = append(rstk, npf)
					}
					break
				}
			}
			if j >= al {
				p1.End = j
				stk = append(stk, p1)
			}
		}

		if len(stk) == 2 {
			t0, t1 := stk[0], stk[1]
			l0, l1 := t0.End-(t0.Start+1), t1.End-(t1.Start+1)
			ok := false
			if t0.End == leftMax+1 && t1.Start == leftMax-1 {
				ok = true
			} else if t1.End == leftMax+1 && t0.Start == leftMax-1 {
				t0, t1 = t1, t0
				l0, l1 = l1, l0
				ok = true
			}
			if ok {
				if l0 >= 3 || l1 >= 3 {
					var pf PathFreq
					pf.IDArr = make([]uint32, l0+l1-1)
					copy(pf.IDArr, t0.PF.IDArr[t0.Start+1:t0.End])
					copy(pf.IDArr[l0:], t1.PF.IDArr[t1.Start+1+1:])
					edgesArr[i].NGSPathArr[0] = pf
					edgesArr[i].NGSPathArr = edgesArr[i].NGSPathArr[:1]
					fmt.Printf("[AddNGSPathToDBG]eID:%d el:%d Unique:%v|Bubble:%v NGSPathArr:%v\n", eID, e.GetSeqLen(), e.GetUniqueFlag() > 0, e.GetBubbleFlag() > 0, edgesArr[eID].NGSPathArr)
				} else {
					edgesArr[i].NGSPathArr = nil
				}
				//fmt.Printf("[AddNGSPathToDBG]eID:%v NGSPathArr:%v\n", eID, edgesArr[i].NGSPathArr)
				continue
			}
		}
		if len(stk) <= AllowMaxPathNum {
			//fmt.Printf("[AddNGSPathToDBG]eID:%d el:%d stk:%v\n", eID, e.GetSeqLen(), stk)
			edgesArr[i].NGSPathArr = make([]PathFreq, 0, len(stk))
			for _, pfe := range stk {
				var pf PathFreq
				pf.IDArr = pfe.PF.IDArr[pfe.Start+1 : pfe.End]
				edgesArr[i].NGSPathArr = append(edgesArr[i].NGSPathArr, pf)
				//if len(pf.IDArr) >= 3 {
				//}
			}
			//fmt.Printf("[AddNGSPathToDBG]eID:%d el:%d Unique:%v|Bubble:%v NGSPathArr:%v\n", eID, e.GetSeqLen(), e.GetUniqueFlag() > 0, e.GetBubbleFlag() > 0, edgesArr[i].NGSPathArr)
			if len(edgesArr[i].NGSPathArr) > 1 {
				count := 0
				//pa := edgesArr[i].NGSPathArr
				c0 := GetNumInPath0(edgesArr[i].NGSPathArr, eID)
				for x, pf := range edgesArr[i].NGSPathArr {
					if pf.IDArr[0] != eID {
						edgesArr[i].NGSPathArr[count] = pf
						count++
						continue
					}
					ok := true
					for y, pf2 := range edgesArr[i].NGSPathArr {
						if y == x {
							continue
						}
						idx2 := IndexUint32(pf2.IDArr, eID)
						if len(pf2.IDArr)-idx2 >= len(pf.IDArr) && Equaluint32Arr(pf.IDArr, pf2.IDArr[idx2:idx2+len(pf.IDArr)]) {
							ok = false
							break
						} else if c0 == 1 && len(pf2.IDArr)-idx2 < len(pf.IDArr) && Equaluint32Arr(pf.IDArr[:len(pf2.IDArr)-idx2], pf2.IDArr[idx2:]) {
							edgesArr[i].NGSPathArr[y].IDArr = append(pf2.IDArr, pf.IDArr[len(pf2.IDArr)-idx2:]...)
						}
					}
					if ok {
						edgesArr[i].NGSPathArr[count] = pf
						count++
					}
				}
				edgesArr[i].NGSPathArr = edgesArr[i].NGSPathArr[:count]
			}
			if len(edgesArr[i].NGSPathArr) > 0 {
				fmt.Printf("[AddNGSPathToDBG]eID:%d el:%d Unique:%v|Bubble:%v NGSPathArr:%v\n", eID, e.GetSeqLen(), e.GetUniqueFlag() > 0, e.GetBubbleFlag() > 0, edgesArr[eID].NGSPathArr)
			}
			//fmt.Printf("[AddNGSPathToDBG]eID:%v NGSPathArr:%v\n", eID, stk)
		} else {
			morePathNum++
			totalPathNum += len(stk)
			edgesArr[i].NGSPathArr = nil
		}
	}

	fmt.Printf("[AddNGSPathToDBG]repeatNum:%d otherNum:%d\n", repeatNum, otherNum)
	fmt.Printf("[AddNGSPathToDBG]tooManyPathNum:%d not found any NGS path edge num:%d morePathNum:%d avg Path Num:%d\n", tooManyPathNum, noPathNum, morePathNum, totalPathNum/(morePathNum+1))

	// cat path just unique path extend
	noPathNum = 0
	catEdgeNum := 0
	loopMax := 1
	for i := 0; i < loopMax; i++ {
		for _, e := range edgesArr {
			if e.ID < 2 || e.GetDeleteFlag() > 0 {
				continue
			}
			eID := e.ID
			pa := e.NGSPathArr
			if len(pa) == 0 {
				noPathNum++
				continue
			}
			//pl := len(pa)
			EqualFlag := false
			for j, pf := range pa {
				idxArr, ix2Arr := GetUniquePathIdxArr(pa, pf, eID, edgesArr, EqualFlag)
				// process unique path
				if len(idxArr) < 1 {
					continue
				}
				l1 := len(pf.IDArr)
				idx := IndexUint32(pf.IDArr, eID)
				fmt.Printf("[AddNGSPathToDBG]eID:%v pf:%v idxArr:%v\n", eID, pf, idxArr)
				var npf PathFreq
				if idxArr[0] < idx {
					x := idxArr[0]
					e2 := edgesArr[pf.IDArr[x]]
					pf2 := e2.NGSPathArr[ix2Arr[0]]
					idx2 := IndexUint32(pf2.IDArr, eID)
					l2 := len(pf2.IDArr)
					y := IndexUint32(pf2.IDArr, e2.ID)
					if y > idx2 {
						pf2.IDArr = GetReverseUint32Arr(pf2.IDArr)
						y = l2 - 1 - y
					}
					lmin := MinInt(x, y)
					rmin := MinInt(l1-x, l2-y)
					if !Equaluint32Arr(pf.IDArr[x-lmin:x+rmin], pf2.IDArr[y-lmin:y+rmin]) {
						log.Fatalf("[AddNGSPathToDBG]pf:%v pf2:%v not found consis path", pf.IDArr, pf2.IDArr)
					}
					npf.IDArr = make([]uint32, y+l1-x)
					copy(npf.IDArr[:y], pf2.IDArr[:y])
					copy(npf.IDArr[y:], pf.IDArr[x:])
				}

				if idxArr[len(idxArr)-1] > idx {
					x := idxArr[len(idxArr)-1]
					e2 := edgesArr[pf.IDArr[x]]
					pf2 := e2.NGSPathArr[ix2Arr[len(ix2Arr)-1]]
					idx2 := IndexUint32(pf2.IDArr, eID)
					l2 := len(pf2.IDArr)
					y := IndexUint32(pf2.IDArr, e2.ID)
					if y < idx2 {
						pf2.IDArr = GetReverseUint32Arr(pf2.IDArr)
						y = l2 - 1 - y
					}
					//lmin := MinInt(x, y)
					//rmin := MinInt(l1-x, l2-y)
					if !Equaluint32Arr(pf.IDArr[x-y:], pf2.IDArr[:y+l1-x]) {
						log.Fatalf("[AddNGSPathToDBG]pf:%v pf2:%v not found consis path", pf.IDArr, pf2.IDArr)
					}

					if l1-x > l2-y {
						log.Fatalf("[AddNGSPathToDBG]pf:%v pf2:%v l1-x:%d > l2-y:%d\n", pf.IDArr, pf2.IDArr, x, y)
					}
					if len(npf.IDArr) == 0 {
						npf.IDArr = append(npf.IDArr, pf.IDArr...)
					}
					npf.IDArr = append(npf.IDArr, pf2.IDArr[y+l1-x:]...)
				}
				fmt.Printf("[AddNGSPathToDBG]eID:%v npf:%v\n", eID, npf)
				if len(npf.IDArr) <= len(pf.IDArr) {
					continue
				}
				catEdgeNum++
				var rnpf PathFreq
				rnpf.IDArr = GetReverseUint32Arr(npf.IDArr)
				for j, t := range idxArr {
					p2 := ix2Arr[j]
					id := pf.IDArr[t]
					e2 := edgesArr[id]
					pf2 := e2.NGSPathArr[p2]
					idx2 := IndexUint32(pf2.IDArr, eID)
					fmt.Printf("[AddNGSPathToDBG]e2ID:%v e2 pf:%v\n", id, pf2)
					y := IndexUint32(pf2.IDArr, id)
					if t < idx {
						if y < idx2 {
							edgesArr[id].NGSPathArr[p2] = npf
						} else {
							edgesArr[id].NGSPathArr[p2] = rnpf
						}
					} else {
						if y > idx2 {
							edgesArr[id].NGSPathArr[p2] = npf
						} else {
							edgesArr[id].NGSPathArr[p2] = rnpf
						}
					}
				}
				edgesArr[eID].NGSPathArr[j] = npf
			}
		}
		fmt.Printf("[AddNGSPathToDBG]loop:%d catEdgeNum:%d\n", i, catEdgeNum)
	}

	// merge unique edge path
	pathMap := make(map[uint32][]uint32)
	eArr := make([]DBGEdge, 0, 10000)
	// merge len(e.NGSPathArr) == 1
	for _, e := range edgesArr {
		if e.ID < 2 || e.GetDeleteFlag() > 0 || len(e.NGSPathArr) != 1 {
			continue
		}
		eID := e.ID
		pa := e.NGSPathArr
		if len(pa) == 0 {
			noPathNum++
			continue
		}
		EqualFlag := true
		fmt.Printf("[AddNGSPathToDBG]eID:%d pa:%v\n", eID, pa)
		for j, pf := range edgesArr[eID].NGSPathArr {
			var npf PathFreq
			npf.IDArr = make([]uint32, len(pf.IDArr))
			copy(npf.IDArr, pf.IDArr)
			pf = npf
			idxArr, _ := GetUniquePathIdxArr(pa, pf, eID, edgesArr, EqualFlag)
			if len(idxArr) < 1 {
				continue
			}
			fmt.Printf("[AddNGSPathToDBG]j:%d pf:%v idxArr:%v\n", j, pf, idxArr)
			//idx := IndexUint32(pf.IDArr, eID)
			//if len(edgesArr[eID].NGSPathArr) == 1 {
			//	idxArr, ix2Arr = AddToIdxArr(idxArr, ix2Arr, idx, j)
			//}
			idx1, idx2 := GetMaxUniqueRegion(idxArr, pf, edgesArr)
			fmt.Printf("[AddNGSPathToDBG]idxArr:%v, idx1:%d idx2:%d\n", idxArr, idx1, idx2)
			if idx1 >= 0 && idx2 >= 0 && idx1 < idx2 {
				start, end := idxArr[idx1], idxArr[idx2]
				path := make([]uint32, end+1-start)
				copy(path, pf.IDArr[start:end+1])
				if len(path) < 3 {
					continue
				}
				if path[0] == path[len(path)-1] {
					continue
				}
				id1, id2 := pf.IDArr[idxArr[idx1]], pf.IDArr[idxArr[idx2]]
				v1, ok1 := pathMap[id1]
				v2, ok2 := pathMap[id2]
				if ok1 != ok2 {
					fmt.Fprintf(os.Stderr, "[AddNGSPathToDBG]path:%v ok1:%v p1:%v ok2:%v p2:%v\n", path, ok1, v1, ok2, v2)
					continue
				}
				if ok1 || ok2 {
					continue
				}
				pathMap[id1] = path
				pathMap[id2] = path

				for x := idx1 + 1; x < idx2; x++ {
					idx := idxArr[x]
					id := pf.IDArr[idx]
					if len(edgesArr[id].NGSPathArr) != 1 {
						continue
					}
					v, ok := pathMap[id]
					if ok {
						if v[0] == path[0] && Equaluint32Arr(v, path) {

						} else if v[0] == path[len(path)-1] && v[len(v)-1] == path[0] {

						} else {
							fmt.Printf("[AddNGSPathToDBG]hasCat pathMap[%d]:%v path:%v \n", id, pathMap[id], path)
						}
					} else {
						pathMap[id] = path
					}
				}
				e := edgesArr[id1]
				fmt.Printf("[AddNGSPathToDBG]path:%v idxArr:%v\n", path, idxArr)
				utg := ConcatEdgePath2(path, edgesArr, nodesArr, opt.Kmer)
				e.= utg
				var eP Path
				eP.IDArr = path
				e.PathMat = append(e.PathMat, eP)
				eArr = append(eArr, e)
			}
		}
	}

	fmt.Printf("[AddNGSPathToDBG]len(eArr):%d\n", len(eArr))
	// merge len(e.NGSPathArr) != 1
	for _, e := range edgesArr {
		if e.ID < 2 || e.GetDeleteFlag() > 0 || len(e.NGSPathArr) == 1 {
			continue
		}
		eID := e.ID
		pa := e.NGSPathArr
		if len(pa) == 0 {
			noPathNum++
			continue
		}
		EqualFlag := true
		fmt.Printf("[AddNGSPathToDBG]eID:%d pa:%v\n", eID, pa)
		for j, pf := range edgesArr[eID].NGSPathArr {
			var npf PathFreq
			npf.IDArr = make([]uint32, len(pf.IDArr))
			copy(npf.IDArr, pf.IDArr)
			pf = npf
			idxArr, _ := GetUniquePathIdxArr(pa, pf, eID, edgesArr, EqualFlag)
			if len(idxArr) < 1 {
				continue
			}
			fmt.Printf("[AddNGSPathToDBG]j:%d pf:%v idxArr:%v\n", j, pf, idxArr)
			//idx := IndexUint32(pf.IDArr, eID)
			//if len(edgesArr[eID].NGSPathArr) == 1 {
			//	idxArr, ix2Arr = AddToIdxArr(idxArr, ix2Arr, idx, j)
			//}
			idx1, idx2 := GetMaxUniqueRegion(idxArr, pf, edgesArr)
			fmt.Printf("[AddNGSPathToDBG]idxArr:%v, idx1:%d idx2:%d\n", idxArr, idx1, idx2)
			if idx1 >= 0 && idx2 >= 0 && idx1 < idx2 {
				start, end := idxArr[idx1], idxArr[idx2]
				path := make([]uint32, end+1-start)
				copy(path, pf.IDArr[start:end+1])
				if len(path) < 3 {
					continue
				}
				if path[0] == path[len(path)-1] {
					continue
				}
				id1, id2 := pf.IDArr[idxArr[idx1]], pf.IDArr[idxArr[idx2]]
				v1, ok1 := pathMap[id1]
				v2, ok2 := pathMap[id2]
				if ok1 != ok2 {
					fmt.Fprintf(os.Stderr, "[AddNGSPathToDBG]path:%v ok1:%v p1:%v ok2:%v p2:%v\n", path, ok1, v1, ok2, v2)
					continue
				}
				if ok1 || ok2 {
					continue
				}
				pathMap[id1] = path
				pathMap[id2] = path

				for x := idx1 + 1; x < idx2; x++ {
					idx := idxArr[x]
					id := pf.IDArr[idx]
					if len(edgesArr[id].NGSPathArr) != 1 {
						continue
					}
					v, ok := pathMap[id]
					if ok {
						if v[0] == path[0] && Equaluint32Arr(v, path) {

						} else if v[0] == path[len(path)-1] && v[len(v)-1] == path[0] {

						} else {
							fmt.Printf("[AddNGSPathToDBG]hasCat pathMap[%d]:%v path:%v \n", id, pathMap[id], path)
						}
					} else {
						pathMap[id] = path
					}
				}
				e := edgesArr[id1]
				fmt.Printf("[AddNGSPathToDBG]path:%v idxArr:%v\n", path, idxArr)
				utg := ConcatEdgePath2(path, edgesArr, nodesArr, opt.Kmer)
				e.= utg
				var eP Path
				eP.IDArr = path
				e.PathMat = append(e.PathMat, eP)
				eArr = append(eArr, e)
			}
		}
	}

	fmt.Printf("[AddNGSPathToDBG]len(eArr):%d\n", len(eArr))
	mergePathNum := len(eArr)
	// change utg
	for _, e0 := range eArr {
		eID := e0.ID
		path := e0.PathMat[0].IDArr
		e0.PathMat = nil
		edgesArr[eID].= e0.Utg
		edgesArr[eID].ResetDeleteFlag()
		for x := 1; x < len(path); x++ {
			id := path[x]
			if len(edgesArr[id].NGSPathArr) == 1 {
				v, ok := pathMap[id]
				if ok && v[0] == id {
				} else {
					edgesArr[id].SetDeleteFlag()
				}
			}
		}
		fmt.Printf("[AddNGSPathToDBG]eID:%d path:%v\n", eID, path)
		e, e1, elast2, e2 := edgesArr[path[0]], edgesArr[path[1]], edgesArr[path[len(path)-2]], edgesArr[path[len(path)-1]]
		nID := GetShareNodeID(e, e1)
		nID2 := GetShareNodeID(elast2, e2)
		if e.EndNID == nID {
			if e.EndNID > 0 && !SubstituteEdgeID(nodesArr, e.EndNID, e.ID, 0) {
				log.Fatalf("[AddNGSPathToDBG]v: %v\ne.ID:%d substitute by 0 failed\n", nodesArr[e.EndNID], e.ID)
			}
			if e2.StartNID == nID2 {
				if e2.EndNID > 0 && !SubstituteEdgeID(nodesArr, e2.EndNID, e2.ID, e.ID) {
					log.Fatalf("[AddNGSPathToDBG]v: %v\ne2.ID:%d substitute by e.ID:%d  failed\n", nodesArr[e2.EndNID], e2.ID, e.ID)
				}
				edgesArr[eID].EndNID = e2.EndNID
			} else {
				if e2.StartNID > 0 && !SubstituteEdgeID(nodesArr, e2.StartNID, e2.ID, e.ID) {
					log.Fatalf("[AddNGSPathToDBG]v: %v\ne2.ID:%d substitute by e.ID:%d  failed\n", nodesArr[e2.StartNID], e2.ID, e.ID)
				}
				edgesArr[eID].EndNID = e2.StartNID
			}
		} else {
			if e.StartNID > 0 && !SubstituteEdgeID(nodesArr, e.StartNID, e.ID, 0) {
				log.Fatalf("[AddNGSPathToDBG]v: %v\ne.ID:%d substitute by 0 failed\n", nodesArr[e.StartNID], e.ID)
			}
			if e2.StartNID == nID2 {
				if e2.EndNID > 0 && !SubstituteEdgeID(nodesArr, e2.EndNID, e2.ID, e.ID) {
					log.Fatalf("[AddNGSPathToDBG]v: %v\ne2.ID:%d substitute by e.ID:%d  failed\n", nodesArr[e2.EndNID], e2.ID, e.ID)
				}
				edgesArr[eID].StartNID = e2.EndNID
			} else {
				if e2.StartNID > 0 && !SubstituteEdgeID(nodesArr, e2.StartNID, e2.ID, e.ID) {
					log.Fatalf("[AddNGSPathToDBG]v: %v\ne2.ID:%d substitute by e.ID:%d  failed\n", nodesArr[e2.StartNID], e2.ID, e.ID)
				}
				edgesArr[eID].StartNID = e2.StartNID
			}
		}
	}

	for i, e := range edgesArr {
		if e.ID < 2 || e.GetDeleteFlag() > 0 || len(e.NGSPathArr) == 0 {
			continue
		}

		for x, pf := range e.NGSPathArr {
			for _, id := range pf.IDArr {
				v, ok := pathMap[id]
				if ok {
					if v[0] != id {
						id = v[0]
						v = GetReverseUint32Arr(v)
					}
					fmt.Printf("[AddNGSPathToDBG]v:%v pf:%v id:%v\n", v, pf, id)
					edgesArr[i].NGSPathArr[x] = RePlaceNGSPath(pf, v, id)
					if len(pf.IDArr) != len(edgesArr[i].NGSPathArr[x].IDArr) {
						fmt.Printf("[AddNGSPathToDBG]pf:%v, RePlaceNGSPath:%v\n", pf, edgesArr[i].NGSPathArr[x])
					}
				}
			}
		}
	}

	fmt.Printf("[AddNGSPathToDBG]mergePathNum:%d\n", mergePathNum)

	numArr := make([]int, AllowMaxPathNum+1)
	for _, e := range edgesArr {
		if e.ID < 2 || e.GetDeleteFlag() > 0 {
			continue
		}
		pa := e.NGSPathArr
		pl := len(pa)
		if pl > AllowMaxPathNum {
			pl = AllowMaxPathNum
		}
		numArr[pl]++
	}

	for i, num := range numArr {
		fmt.Printf("[AddNGSPathToDBG]numArr[%d]:%d\n", i, num)

	}
	return
}
*/

func InitUint32Slice(pa *[WindowSize / 4]uint32) {
	for i := range pa {
		pa[i] = 0
	}
}

func WriteNGSPath(pathfn string, wc <-chan []uint32, pathArrPool *sync.Pool, writeNumC chan<- int) {
	pathfp, err := os.Create(pathfn)
	if err != nil {
		log.Fatalf("[WriteNGSPath] file %s create error, err: %v\n", pathfn, err)
	}
	//sbuf := make([]byte, 0, 100)
	defer pathfp.Close()
	zfp, err1 := zstd.NewWriter(pathfp, zstd.WithEncoderCRC(false), zstd.WithEncoderConcurrency(1), zstd.WithEncoderLevel(1))
	defer zfp.Close()
	if err1 != nil {
		log.Fatalf("[WriteNGSPath]open write file:%s err:%v\n", pathfn, err1)
	}
	//var finishNum int
	readNum := 0
	var pa [WindowSize / 4]uint32
	InitUint32Slice(&pa)
	//var pathArr PathArr
	for {
		arr, ok := <-wc
		if !ok {
			break
		}
		copy(pa[:], arr)
		readNum += len(arr)
		binary.Write(zfp, binary.LittleEndian, pa)
		pathArrPool.Put(arr)
		InitUint32Slice(&pa)
		//pathArr.Arr = append(pathArr.Arr, &rp)
	}
	//time.Sleep(time.Second * 10)
	fmt.Printf("[WriteNGSPath] write NGSPath num:%d to file:%s\n", readNum, pathfn)
	writeNumC <- readNum
}

type MapingInfo struct {
	ID        uint32
	PathMat   [][][]uint32
	Anotition string
}

type optionsMN struct {
	ArgsOpt
	TipMaxLen     int
	WinSize       int
	MaxNGSReadLen int
	MinNGSReadLen int
	MinMapFreq    int
}

func checkArgsMN(c cli.Command) (opt optionsMN, suc bool) {
	var tmp int
	var err error
	gOpt, suc := CheckGlobalArgs(c.Parent())
	if suc == false {
		log.Fatalf("[checkArgsMN] check global Arguments error, opt: %v\n", gOpt)
	}
	opt.ArgsOpt = gOpt
	opt.MaxNGSReadLen = c.Flag("MaxNGSReadLen").Get().(int)
	opt.MinNGSReadLen = c.Flag("MinNGSReadLen").Get().(int)
	opt.WinSize = c.Flag("WinSize").Get().(int)
	if opt.WinSize < 1 || opt.WinSize > 20 {
		log.Fatalf("[checkArgsMN] argument 'WinSize':%v, must bewteen [3~20] set error: %v\n", c.Flag("WinSize"), err)
	}
	tmp = c.Flag("TipMaxLen").Get().(int)
	if tmp == 0 {
		opt.TipMaxLen = opt.MaxNGSReadLen
	} else {
		if tmp < 0 || tmp < opt.Kmer*2 || tmp > opt.MaxNGSReadLen {
			log.Fatalf("[checkArgsMN] argument 'tipMaxLen': %v must between [%v~%v]\n", tmp, opt.Kmer*2, opt.MaxNGSReadLen)
		}
		opt.TipMaxLen = tmp
	}
	opt.MinMapFreq = c.Flag("MinMapFreq").Get().(int)
	if opt.MinMapFreq < 1 || opt.MinMapFreq > MaxC {
		log.Fatalf("[checkArgsMN]the argument 'MinMapFreq': %v must [%v ~ %v]\n", opt.MinMapFreq, 1, MaxC)
	}
	suc = true
	return opt, suc
}

func MapingNGSFindDBGPath(opt optionsMN, nodesArr []DBGNode, edgesArr []DBGEdge) int {
	// construct cuckoofilter of DBG sample
	cfSize := GetCuckoofilterDBGSampleSize(edgesArr, opt.WinSize, opt.Kmer)
	fmt.Printf("[MapingNGSFindDBGPath]cfSize:%d\n", cfSize)
	cf := MakeCuckooFilterDBGKmer(uint64(cfSize * 2))
	fmt.Printf("[MapingNGSFindDBGPath]cf.BucketPow:%d\n", cf.BucketPow)
	ConstructCFDBGMinimizers(&cf, edgesArr, opt.Kmer, opt.WinSize)
	fmt.Printf("[MapingNGSFindDBGPath]construct Smaple of DBG edges cuckoofilter number is:%d\n", cf.Count)

	cfgInfo, err := ParseCfg(opt.CfgFn, false, false)
	if err != nil {
		log.Fatal("[ParseCfg] found err")
	}
	fmt.Printf("[MapingNGSFindDBGPath]cfgInfo: %v\n", cfgInfo)

	runtime.GOMAXPROCS(opt.NumCPU)

	writeNumC := make(chan int, 1)
	wc := make(chan []uint32, opt.NumCPU)
	//defer close(wc)
	bufSize := WindowSize
	pathArrPool := sync.Pool{New: func() interface{} {
		arr := make([]uint32, bufSize/4)
		return arr
	}}

	pathfn := opt.Prefix + ".NGSPath.zst"
	go WriteNGSPath(pathfn, wc, &pathArrPool, writeNumC)
	for _, lib := range cfgInfo.Libs {
		if lib.AsmFlag != AllState && lib.SeqProfile != 1 {
			continue
		}
		for i := 0; i < len(lib.FnName); i++ {
			ParaMapReadsFile(lib.FnName[i], wc, &pathArrPool, opt.NumCPU-1, nodesArr, edgesArr, cf, opt)
		}
	}
	close(wc)
	return <-writeNumC
}

func NGSPathSizeWriter(fn string, size int) {
	fp, err := os.Create(fn)
	if err != nil {
		log.Fatalf("[NGSPathSizeWriter] file:%s create error:%v\n", fn, err)
	}
	defer fp.Close()
	fmt.Fprintf(fp, "NGSPath size:%d\n", size)
}

func NGSPathSizeReader(fn string) (size int) {
	fp, err := os.Open(fn)
	if err != nil {
		log.Fatalf("[NGSPathSizeReader] file:%s Open error:%v\n", fn, err)
	}
	defer fp.Close()
	_, err = fmt.Fscanf(fp, "NGSPath size:%d\n", &size)
	if err != nil {
		log.Fatalf("[NGSPathSizeReader]file:%s parse error:%v\n", fn, err)
	}
	return
}

func MapNGS(c cli.Command) {

	t0 := time.Now()
	// check agruments

	opt, suc := checkArgsMN(c)
	if suc == false {
		log.Fatalf("[MapNGS] check Arguments error, opt:%v\n", opt)
	}
	fmt.Printf("Arguments: %v\n", opt)
	/*DebugModel := true
	if DebugModel {
		profileFn := opt.Prefix + ".MapNGS.prof"
		cpuprofilefp, err := os.Create(profileFn)
		if err != nil {
			log.Fatalf("[MapNGS] open cpuprofile file: %v failed\n", profileFn)
		}
		pprof.StartCPUProfile(cpuprofilefp)
		defer pprof.StopCPUProfile()
	}*/

	// read files and construt DBG
	DBGInfofn := opt.Prefix + ".smfy.DBGInfo"
	eSize, nSize := DBGInfoReader(DBGInfofn)
	nodesfn := opt.Prefix + ".nodes.smfy.Arr"
	nodesArr := make([]DBGNode, nSize)
	fc := make(chan int, 1)
	go NodesArrReader(nodesfn, nodesArr, opt.Kmer, fc)
	edgesfn := opt.Prefix + ".edges.smfy.fa.zst"
	//LoadEdgesfqFromFn(edgesfn, edgesArr, true)
	edgesArr := ReadEdgesFromFile(edgesfn, uint32(eSize), opt.Kmer)
	// check NodesArrReader() finished
	<-fc

	CheckInterConnectivity(edgesArr, nodesArr)
	AddNodeInfo2DBGEdgeArr(edgesArr, nodesArr)

	uniqueNum, bubbleRepeatNum, twoCycleNum, selfCycle, selfCycleSameComingNum, bubbleEdgeNum := SetDBGEdgesUniqueFlag(edgesArr, nodesArr, opt.Kmer)
	fmt.Printf("[MapNGS] unique edge number is:%d bubbleRepeatNum:%d twoCycleNum:%d selfCycle:%d selfCycleSameComingNum:%d bubbleEdgeNum:%d\n", uniqueNum, bubbleRepeatNum, twoCycleNum, selfCycle, selfCycleSameComingNum, bubbleEdgeNum)
	//PrintTmpDBG(nodesArr, edgesArr, opt.Prefix)

	// map Illumina reads to the DBG and find reads map path
	NGSPathfn := opt.Prefix + ".NGSPath.zst"
	NGSPathSizefn := opt.Prefix + ".NGSPathSize"
	if _, err := os.Stat(NGSPathfn); err != nil {
		size := MapingNGSFindDBGPath(opt, nodesArr, edgesArr)
		NGSPathSizeWriter(NGSPathSizefn, size)
		t2 := time.Now()
		fmt.Printf("[MapNGS]MapingNGSFindDBGPath:%v\n", t2.Sub(t0))
		runtime.GC()
	}
	ngsPathSize := NGSPathSizeReader(NGSPathSizefn)
	AddNGSPathToDBG(NGSPathfn, nodesArr, edgesArr, opt, ngsPathSize)
	var optSF optionsSF
	optSF.ArgsOpt, optSF.TipMaxLen, optSF.MinMapFreq = opt.ArgsOpt, opt.TipMaxLen, opt.MinMapFreq
	SmfyDBG(nodesArr, edgesArr, optSF)

	//edgesArr = SetEdgeID(nodesArr, edgesArr)
	//nodesArr = SetNodeID(nodesArr, edgesArr)

	CheckDBGSelfCycle(nodesArr, edgesArr, opt.Kmer)
	CheckInterConnectivity(edgesArr, nodesArr)

	//store DBG
	mapDBGNodesfn := opt.Prefix + ".nodes.MapDBG.Arr"
	fc2 := make(chan int, 1)
	go NodesArrWriter(nodesArr, mapDBGNodesfn, fc2)

	mapDBGEdgesfn := opt.Prefix + ".edges.MapDBG.fa.zst"
	StoreEdgesToFn(mapDBGEdgesfn, edgesArr)
	<-fc2
	//mappingEdgefn := opt.Prefix + ".edges.mapping.fa"
	// StoreMappingEdgesToFn(mappingEdgefn, edgesArr, opt.MaxMapEdgeLen)
	//	adpaterEdgesfn := prefix + ".edges.adapter.fq"
	//	StoreEdgesToFn(adpaterEdgesfn, edgesArr, true)

	DBGInfofn = opt.Prefix + ".MapDBG.DBGInfo"
	DBGInfoWriter(DBGInfofn, len(edgesArr), len(nodesArr))
	/*wrFn := opt.Prefix + ".smfy.NGSAlignment"
	MapNGS2DBG(opt, nodesArr, edgesArr, wrFn)
	//CheckInterConnectivity(edgesArr, nodesArr)
	SimplifyByNGS(opt, nodesArr, edgesArr, wrFn)
	SmfyDBG(nodesArr, edgesArr, opt)
	*/
	//CheckInterConnectivity(edgesArr, nodesArr)
	// simplify DBG
	//IDMapPath := ConstructIDMapPath(joinPathArr)
	//AdjustPathMat(edgesArr, nodesArr, joinPathArr, IDMapPath)

	t2 := time.Now()
	fmt.Printf("[MapNGS]MapNGS total used:%v\n", t2.Sub(t0))
}
