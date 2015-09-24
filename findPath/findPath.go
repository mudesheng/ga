package findPath

import (
	"bufio"
	// "compress/gzip"
	// "container/list"
	// "encoding/binary"
	// "encoding/gob"
	"container/list"
	"fmt"
	"ga/bnt"
	// "ga/constructcf"
	"ga/constructdbg"
	// "ga/cuckoofilter"
	"io"
	"log"
	"os"
	// "runtime"
	"math"
	"strconv"

	"github.com/awalterschulze/gographviz"
	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/sam"
	"github.com/jwaldrip/odin/cli"
	"strings"
)

var Kmerlen int

const (
	MAX_ALIGN_GAP_ALLOW     = 400
	MAX_SHORT_READ_LEN      = 450
	MAX_LONG_READ_LEN       = 50000
	MIN_PATH_LEN            = 3
	FLANK_ALLOW_FACTOR      = 8
	MAX_FLANK_ALLOW_LEN     = 50
	ALLOW_FACTOR            = 10
	MIN_REF_LEN_ALLOW_MERGE = 1000
	QUY_OVL_FACTOR          = 2
	REF_OVL_FACTOR          = 2 // the min overlap of two connective edge
	MAX_INSERT_SIZE         = 700
	MIN_EDGE_LEN_EXTEND     = 500
	MIN_ALLOW_PATH_FREQ     = 3
)

type PathCrossInfo struct {
	EdgeID    constructdbg.DBG_MAX_INT
	Node      constructdbg.DBGNode
	RemainLen int
}

type LA struct {
	RefID   constructdbg.DBG_MAX_INT
	QuyID   constructdbg.DBG_MAX_INT // query ID
	AlgnLen int
	Diff    int
	Idty    float64
	RefB    int // the alignment start position of reference
	RefE    int // the alignment end position of reference
	RefLen  int // the reference length
	QuyB    int // the alignment start position of Query
	QuyE    int // the alignment end position of Query
	QuyLen  int // the query length
}

type MapInfo struct {
	Start, End int // read start and end mapping position
	RefID      constructdbg.DBG_MAX_INT
	RefStart   int
}

type FastMapRecord struct {
	ReadID string
	Rlen   int32
	Minfo  []MapInfo
}

type GapRegion struct {
	Begin, End int
}

func GetSamRecord(bamfn string, rc chan []sam.Record, numCPU int) {
	fp, err := os.Open(bamfn)
	if err != nil {
		log.Fatalf("[GetSamRecord] open file: %s failed, err: %v\n", bamfn, err)
	}
	defer fp.Close()
	bamfp, err := bam.NewReader(fp, 0)
	if err != nil {
		log.Fatalf("[GetSamRecord] create bam.NewReader err: %v\n", err)
	}
	defer bamfp.Close()
	var rArr []sam.Record
	// var cigar sam.Cigar
	var NM = sam.Tag{'N', 'M'}
	var AS = sam.Tag{'A', 'S'}
	for {
		r, err := bamfp.Read()
		if err != nil {
			break
		}
		if len(r.Cigar) < 2 || len(r.Cigar) > 3 {
			continue
		} else {
			v := r.AuxFields.Get(NM).Value()
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
			}
			// fmt.Printf("[GetSamRecord]AS:%d, cigar: %v\n", asi, r.Cigar.String())
			if asi < Kmerlen {
				continue
			}
		}
		// Debug
		// if r.Cigar.String() == "400M" {
		// 	continue
		// }
		//Debug
		if len(rArr) > 0 {
			if rArr[0].Name != r.Name {
				if len(rArr) >= 3 {
					rc <- rArr
					rArr = nil
				} else {
					rArr = rArr[0:0]
				}
			}
		}
		rArr = append(rArr, *r)
	}
	if len(rArr) >= 3 {
		rc <- rArr
	}
	// send terminal signals
	for i := 0; i < numCPU; i++ {
		var nilArr []sam.Record
		rc <- nilArr
	}
}

func GetFastMapRecord(fastmapfn string, rc chan [2]FastMapRecord, numCPU int) {
	fp, err := os.Open(fastmapfn)
	if err != nil {
		log.Fatalf("[GetFastMapRecord] open file: %s failed, err: %v\n", fastmapfn, err)
	}
	defer fp.Close()
	buffp := bufio.NewReader(fp)
	// if err != nill {
	// 	log.Fatalf("[GetFastMapRecord] create bufio.NewReader err : %v\n", err)
	// }
	// bamfp, err := bam.NewReader(fp, 0)
	// if err != nil {
	// 	log.Fatalf("[GetSamRecord] create bam.NewReader err: %v\n", err)
	// }
	// defer bamfp.Close()
	// var rArr []sam.Record
	// // var cigar sam.Cigar
	// var NM = sam.Tag{'N', 'M'}
	// var AS = sam.Tag{'A', 'S'}
	var pairinfo [2]FastMapRecord
	var SQ string
	for {
		r, err := buffp.ReadString('\n')
		if err != nil {
			break
		}
		fields := strings.Split(strings.TrimSpace(r), "\t")
		if fields[0] == "//" {
			continue
		} else if fields[0] == "SQ" {
			if SQ == "" || fields[1][:len(fields[1])-2] != SQ[:len(SQ)-2] {
				// write last pair
				if len(pairinfo[0].Minfo) > 0 && len(pairinfo[1].Minfo) > 0 && (len(pairinfo[0].Minfo)+len(pairinfo[0].Minfo) > 2 || pairinfo[0].Minfo[0].RefID != pairinfo[1].Minfo[0].RefID) {
					rc <- pairinfo
				}
				if fields[1][len(fields[1])-1:] != "1" {
					log.Fatalf("[GetFastMapRecord] reads pair mapping error\n")
				}
				pairinfo[0].ReadID = fields[1]
				rlen, err := strconv.Atoi(fields[2])
				if err != nil {
					log.Fatalf("[GetFastMapRecord] Atoi err: %v\n", err)
				}
				pairinfo[0].Rlen = int32(rlen)
				pairinfo[0].Minfo = nil
			} else { // same pair
				pairinfo[1].ReadID = fields[1]
				rlen, err := strconv.Atoi(fields[2])
				if err != nil {
					log.Fatalf("[GetFastMapRecord] Atoi err: %v\n", err)
				}
				pairinfo[1].Rlen = int32(rlen)
				pairinfo[1].Minfo = nil
			}
			SQ = fields[1]
		} else if fields[0] == "EM" {
			if fields[3] != "1" {
				// log.Fatalf("[GetFastMapRecord] EM[3]: %s != 1\n%s", fields[3], r)
				fmt.Printf("[GetFastMapRecord] EM[3]: %s != 1\n%s", fields[3], r)
			}

			var info MapInfo
			info.Start, err = strconv.Atoi(fields[1])
			if err != nil {
				log.Fatalf("[GetFastMapRecord] Atoi err: %v\n", err)
			}
			info.End, err = strconv.Atoi(fields[2])
			if err != nil {
				log.Fatalf("[GetFastMapRecord] Atoi err: %v\n", err)
			}
			refinfo := strings.Split(fields[4], ":")
			tmp, err := strconv.Atoi(refinfo[0])
			if err != nil {
				log.Fatalf("[GetFastMapRecord] Atoi err: %v\n", err)
			}
			info.RefID = constructdbg.DBG_MAX_INT(tmp)
			info.RefStart, err = strconv.Atoi(refinfo[1])
			if err != nil {
				log.Fatalf("[GetFastMapRecord] Atoi err: %v\n", err)
			}
			if info.RefStart < 0 {
				info.RefStart = -info.RefStart
			}
			if SQ == pairinfo[0].ReadID {
				pairinfo[0].Minfo = append(pairinfo[0].Minfo, info)
			} else {
				pairinfo[1].Minfo = append(pairinfo[1].Minfo, info)
			}
		} else {
			log.Fatalf("[GetFastMapRecord] readline: %s, err\n", r)
		}

	}
	if len(pairinfo[0].Minfo) > 0 && len(pairinfo[1].Minfo) > 0 && (len(pairinfo[0].Minfo)+len(pairinfo[0].Minfo) > 2 || pairinfo[0].Minfo[0].RefID != pairinfo[1].Minfo[0].RefID) {
		rc <- pairinfo
	}
	// send terminal signals
	for i := 0; i < numCPU; i++ {
		var nilArr [2]FastMapRecord
		rc <- nilArr
	}
}

func AccumulateCigar(cigar sam.Cigar) (Mnum, Inum, Dnum int) {
	for _, co := range cigar {
		if co.Type() == sam.CigarInsertion {
			Inum += co.Len()
		} else if co.Type() == sam.CigarDeletion {
			Dnum += co.Len()
		} else if co.Type() == sam.CigarMatch {
			Mnum += co.Len()
		}
	}

	return
}

func findNextEdge(node constructdbg.DBGNode, eID constructdbg.DBG_MAX_INT, rArr []sam.Record) (neID constructdbg.DBG_MAX_INT, recd sam.Record) {

	var direction uint8
	// find direction
	for i := 0; i < bnt.BaseTypeNum; i++ {
		if node.EdgeIDIncoming[i] == eID {
			direction = constructdbg.FORWARD
			break
		}
		if node.EdgeIDOutcoming[i] == eID {
			direction = constructdbg.BACKWARD
			break
		}
	}
	// get Reference ID arr
	var refIDArr []constructdbg.DBG_MAX_INT
	for _, v := range rArr {
		id, err := strconv.Atoi(v.Ref.Name())
		if err != nil {
			log.Fatalf("[findNextEdge] ref Name not numberical err: %v\n", err)
		}
		refIDArr = append(refIDArr, constructdbg.DBG_MAX_INT(id))
	}

	num := 0
	var teID constructdbg.DBG_MAX_INT
	for i := 0; i < bnt.BaseTypeNum; i++ {
		id := node.EdgeIDIncoming[i]
		if direction == constructdbg.FORWARD {
			id = node.EdgeIDOutcoming[i]
		}
		for j, v := range refIDArr {
			if id == v {
				teID = id
				recd = rArr[j]
				num++
			}
		}
		if num == 1 {
			neID = teID
		}
	}

	return
}

func Cigar2String(cigar sam.Cigar) (cs string) {
	for _, v := range cigar {
		cs += strconv.Itoa(v.Len()) + "---"
	}
	return
}

func isInEdgesArr(arr []constructdbg.DBG_MAX_INT, e constructdbg.DBG_MAX_INT) bool {
	for _, v := range arr {
		if v == e {
			return true
		}
	}
	return false
}

func paraFindShortMappingPath(rc chan [2]FastMapRecord, wc chan []constructdbg.DBG_MAX_INT, edgesArr []constructdbg.DBGEdge, nodesArr []constructdbg.DBGNode) {
	// totalM, totalI, totalD := 0, 0, 0
	for {
		rArr := <-rc
		if len(rArr[0].Minfo) == 0 {
			var guardPath []constructdbg.DBG_MAX_INT
			wc <- guardPath
			break
		}

		// write read1~read2 rrelative path
		for _, p1 := range rArr[0].Minfo {
			var path []constructdbg.DBG_MAX_INT
			path = append(path, p1.RefID)
			if len(rArr[1].Minfo) == 2 && path[0] == rArr[1].Minfo[len(rArr[1].Minfo)-1].RefID {
				continue
			}
			for j := len(rArr[1].Minfo) - 1; j >= 0; j-- {
				path = append(path, rArr[1].Minfo[j].RefID)
			}

			wc <- path
		}

		// write read2~read1 relative path
		for _, p2 := range rArr[1].Minfo {
			var path []constructdbg.DBG_MAX_INT
			path = append(path, p2.RefID)
			if len(rArr[0].Minfo) == 2 && path[0] == rArr[0].Minfo[len(rArr[0].Minfo)-1].RefID {
				continue
			}
			for j := len(rArr[0].Minfo) - 1; j >= 0; j-- {
				path = append(path, rArr[0].Minfo[j].RefID)
			}

			wc <- path
		}

		// if len(rArr) > 5 {
		// 	for i := 0; i < len(rArr); i++ {
		// 		fmt.Printf("[paraFindShortMappingPath]refName: %v, refLen: %v, pos:%v, cigar: %v\n", rArr[i].Ref.Name(), rArr[i].Ref.Len(), rArr[i].Start(), rArr[i].Cigar.String())
		// 	}
		// 	fmt.Printf("[paraFindShortMappingPath]\n")
		// }

		// for _, v := range rArr[:1] {
		// 	var NM, MD, SA sam.Tag
		// 	NM[0] = 'N'
		// 	NM[1] = 'M'
		// 	MD[0] = 'M'
		// 	MD[1] = 'D'
		// 	SA[0] = 'S'
		// 	SA[1] = 'A'
		// 	nm := v.AuxFields.Get(NM).Value()
		// 	var nmv int
		// 	switch nm.(type) {
		// 	case uint8:
		// 		nmv = int(nm.(uint8))
		// 	case uint16:
		// 		nmv = int(nm.(uint16))
		// 	}
		// 	// nmv, err := strconv.Atoi(nm.String()[5:])
		// 	// if err != nil {
		// 	// 	log.Fatalf("[paraFindPacbioMappingPath] convert 'NM' err: %v\n", nmv)
		// 	// }
		// 	// fmt.Printf("[paraFindPacbioMappingPath]nm type: %v\n", v.AuxFields.Get(NM).Type())
		// 	// if nm != nil {
		// 	// 	// snm = nm.Value()
		// 	// } else {
		// 	// 	log.Fatalf("[paraFindPacbioMappingPath] nm: %v\n", v.Ref)
		// 	// }
		// 	Mnum, Inum, Dnum := AccumulateCigar(v.Cigar)
		// 	if Mnum < 1000 {
		// 		continue
		// 	}
		// 	fmt.Printf("ref\tpos\tmapQ\tNM\tCigar\n")
		// 	totalM += Mnum
		// 	totalI += Inum
		// 	totalD += Dnum
		// 	totalMis += (int(nmv) - Inum - Dnum)
		// 	// as := Cigar2String(acgr)
		// 	fmt.Printf("%s\t%d\t%v\t%v\n", v.Ref.Name(), v.Pos, v.MapQ, v.Cigar)
		// 	// sav := v.AuxFields.Get(SA)
		// 	// // fmt.Printf("%v\n", sav)
		// 	// if sav != nil {
		// 	// 	sa := v.AuxFields.Get(SA).Value().(string)
		// 	// 	// fmt.Printf("%v\n", sa)
		// 	// 	// sav := strings.Split(sa, ":")
		// 	// 	// fmt.Printf("%v\n", sav)
		// 	// 	saArr := strings.Split(sa[:len(sa)-1], ";")
		// 	// 	for _, e := range saArr {
		// 	// 		arr := strings.Split(e, ",")
		// 	// 		// fmt.Printf("%v\n", arr)
		// 	// 		// cgr, _ := sam.ParseCigar([]byte(arr[3]))
		// 	// 		// cr := AccumulateCigar(cgr)
		// 	// 		// as := Cigar2String(cr)
		// 	// 		fmt.Printf("%s\t%s\t%s\t%s\t%s\n", arr[0], arr[1], arr[4], arr[5], arr[3])
		// 	// 	}

		// 	// }
		// }
	}

	// fmt.Printf("[paraFindPacbioMappingPath] totalM: %d, totalMis: %d, totalI: %d, totalD: %d\n", totalM, totalMis, totalI, totalD)
}

func IsCyclePath(path []constructdbg.DBG_MAX_INT) bool {
	cycle := false
	for i := 0; i < len(path); i++ {
		for j := 0; j < i; j++ {
			if path[i] == path[j] {
				cycle = true
				break
			}
		}
	}
	return cycle
}

func IsNodeCyclePath(path []constructdbg.DBG_MAX_INT, edgesArr []constructdbg.DBGEdge) bool {
	// node array
	var narr []constructdbg.DBG_MAX_INT
	for i := 1; i < len(path); i++ {
		e1 := edgesArr[path[i-1]]
		e2 := edgesArr[path[i]]
		if e1.StartNID == e2.StartNID || e1.StartNID == e2.EndNID {
			if e1.EndNID == e2.StartNID || e1.EndNID == e2.EndNID {
				return true
			} else {
				narr = append(narr, e1.StartNID)
			}
		} else if e1.EndNID == e2.StartNID || e1.EndNID == e2.EndNID {
			narr = append(narr, e1.EndNID)
		} else {
			log.Fatalf("[IsNodeCyclePath] not found Identity Node e1: %v\ne2: %v\n", e1, e2)
		}
	}
	var start, end constructdbg.DBG_MAX_INT
	// fmt.Printf("[IsNodeCyclePath] path: %v, narr: %v\n", path, narr)
	if edgesArr[path[0]].StartNID == narr[0] {
		start = edgesArr[path[0]].EndNID
	} else {
		start = edgesArr[path[0]].StartNID
	}
	if edgesArr[path[len(path)-1]].StartNID == narr[len(narr)-1] {
		end = edgesArr[path[len(path)-1]].EndNID
	} else {
		end = edgesArr[path[len(path)-1]].StartNID
	}
	if start > 0 && start == end {
		return true
	}
	if start > 0 {
		var tmp []constructdbg.DBG_MAX_INT
		tmp = append(tmp, start)
		narr = append(tmp, narr...)
	}
	if end > 0 {
		narr = append(narr, end)
	}

	// check cycle
	return IsCyclePath(narr)

}

func WriteShortPathToDBG(wc chan []constructdbg.DBG_MAX_INT, edgesArr []constructdbg.DBGEdge, numCPU int) {
	terminalNum := 0
	totalPathNum, totalPathLen := 0, 0
	for {
		path := <-wc
		if len(path) == 0 {
			terminalNum++
			if terminalNum == numCPU {
				break
			} else {
				continue
			}

		}

		// write path to the DBG edgesArr
		{
			if IsCyclePath(path) {
				continue
			}
			edgesArr[path[0]].InsertPathToEdge(path[1:], 1)
			// path = constructdbg.ReverseDBG_MAX_INTArr(path)
			// edgesArr[path[0]].InsertPathToEdge(path, 1)
		}
		totalPathNum++
		totalPathLen += len(path)
		// if len(path) > 5 {
		// 	fmt.Printf("[WriteShortPathToDBG] pathArr: %v\n", path)
		// }
	}

	fmt.Printf("[WriteShortPathToDBG] total path num: %d, total path length: %d\n", totalPathNum, totalPathLen)
}

func GetNextEID(eID constructdbg.DBG_MAX_INT, node constructdbg.DBGNode) (neID constructdbg.DBG_MAX_INT) {
	var direction uint8
	for i := 0; i < bnt.BaseTypeNum; i++ {
		if node.EdgeIDIncoming[i] == eID {
			direction = constructdbg.BACKWARD
			break
		}
		if node.EdgeIDOutcoming[i] == eID {
			direction = constructdbg.FORWARD
			break
		}
	}
	if direction != constructdbg.BACKWARD && direction != constructdbg.FORWARD {
		log.Fatalf("[GetNextEID] direction not set\n")
	}
	num := 0
	for i := 0; i < bnt.BaseTypeNum; i++ {
		if direction == constructdbg.FORWARD {
			if node.EdgeIDIncoming[i] > 0 {
				num++
				neID = node.EdgeIDIncoming[i]
			}
		} else {
			if node.EdgeIDOutcoming[i] > 0 {
				num++
				neID = node.EdgeIDOutcoming[i]
			}
		}
	}
	if num != 1 {
		log.Fatalf("[GetNextEID] found %d edges\n", num)
	}
	return neID
}

func GetNextEArr(eID constructdbg.DBG_MAX_INT, node constructdbg.DBGNode) (eIDArr []constructdbg.DBG_MAX_INT) {
	var direction uint8
	for i := 0; i < bnt.BaseTypeNum; i++ {
		if node.EdgeIDIncoming[i] == eID {
			direction = constructdbg.BACKWARD
			break
		}
		if node.EdgeIDOutcoming[i] == eID {
			direction = constructdbg.FORWARD
			break
		}
	}
	if direction != constructdbg.BACKWARD && direction != constructdbg.FORWARD {
		log.Fatalf("[GetNextEID] direction not set\n")
	}
	for i := 0; i < bnt.BaseTypeNum; i++ {
		if direction == constructdbg.FORWARD {
			if node.EdgeIDIncoming[i] > 0 {
				eIDArr = append(eIDArr, node.EdgeIDIncoming[i])
			}
		} else {
			if node.EdgeIDOutcoming[i] > 0 {
				eIDArr = append(eIDArr, node.EdgeIDOutcoming[i])
			}
		}
	}
	return eIDArr
}

func IsDirectionUniqueEdge(edge constructdbg.DBGEdge, node constructdbg.DBGNode) bool {
	var direction uint8
	in, out := 0, 0
	for i := 0; i < bnt.BaseTypeNum; i++ {
		if node.EdgeIDIncoming[i] > 0 {
			in++
			if node.EdgeIDIncoming[i] == edge.ID {
				direction = constructdbg.BACKWARD
			}
		}
		if node.EdgeIDOutcoming[i] > 0 {
			out++
			if node.EdgeIDOutcoming[i] == edge.ID {
				direction = constructdbg.FORWARD
			}
		}
	}
	if direction == constructdbg.FORWARD {
		if in == 1 {
			return true
		} else {
			return false
		}
	} else if direction == constructdbg.BACKWARD {
		if out == 1 {
			return true
		} else {
			return false
		}
	} else {
		return false
		// log.Fatalf("[IsDirectionUniqueEdge] direction set error\nnode:%v\nedge:%v\n", node, edge)
	}

	return false
}

func IsUniqueEdge(edge constructdbg.DBGEdge, nodesArr []constructdbg.DBGNode) bool {
	if IsDirectionUniqueEdge(edge, nodesArr[edge.StartNID]) == false {
		return false
	}
	if IsDirectionUniqueEdge(edge, nodesArr[edge.EndNID]) == false {
		return false
	}
	return true
}

func MergePath(pathMat []constructdbg.Path, beID constructdbg.DBG_MAX_INT) (uniquePath constructdbg.Path, num int) {
	for _, p := range pathMat {
		if p.IDArr[1] == beID {
			if len(uniquePath.IDArr) == 0 {
				uniquePath = p
				num++
			} else {
				for j, eID := range p.IDArr {
					if j < len(uniquePath.IDArr) {
						if eID != uniquePath.IDArr[j] {
							num++
							break
						}
					} else {
						break
					}
				}
				if num == 1 {
					if len(uniquePath.IDArr) < len(p.IDArr) {
						uniquePath = p
					}
				} else {
					break
				}
			}
		}
	}
	return
}

func IndexEID(arr []constructdbg.DBG_MAX_INT, eID constructdbg.DBG_MAX_INT) int {
	for i, id := range arr {
		if id == eID {
			return i
		}
	}
	return -1
}

func IndexUniqueEdge(consisPathArr []constructdbg.DBG_MAX_INT, edgesArr []constructdbg.DBGEdge, direction uint8) int {
	if direction == constructdbg.FORWARD {
		for i, eID := range consisPathArr {
			if edgesArr[eID].GetUniqueFlag() > 0 {
				return i
			}
		}
	} else if direction == constructdbg.BACKWARD {
		for i := len(consisPathArr) - 1; i >= 0; i-- {
			if edgesArr[consisPathArr[i]].GetUniqueFlag() > 0 {
				return i
			}
		}
	} else {
		log.Fatalf("[IndexUniqueEdge] direction set error\n")
	}

	return -1
}

func GetReverseDBG_MAX_INTArr(arr []constructdbg.DBG_MAX_INT) []constructdbg.DBG_MAX_INT {
	la := len(arr)
	rarr := make([]constructdbg.DBG_MAX_INT, la)
	for i := la - 1; i >= 0; i-- {
		rarr[la-i-1] = arr[i]
	}

	return rarr

}

func findConsistenceAndMergePath(consisPathArr1, consisPathArr2 []constructdbg.DBG_MAX_INT, eID1, eID2 constructdbg.DBG_MAX_INT) (mergePathArr []constructdbg.DBG_MAX_INT) {
	idx1 := IndexEID(consisPathArr1, eID1)
	idx2 := IndexEID(consisPathArr1, eID2)

	idx3 := IndexEID(consisPathArr2, eID1)
	idx4 := IndexEID(consisPathArr2, eID2)
	if idx1 < 0 || idx2 < 0 || idx3 < 0 || idx4 < 0 {
		log.Fatalf("[findConsistenceAndMergePath] eID not found in the consisPathArr\n")
	}

	if idx1 > idx2 {
		if idx3 < idx4 {
			log.Fatalf("[findConsistenceAndMergePath] two consisPathArr not consistence arr1: %v, arr2:%v\n", consisPathArr1, consisPathArr2)
		}
		idx1, idx2 = idx2, idx1
		idx3, idx4 = idx4, idx3
	} else {
		if idx3 > idx4 {
			log.Fatalf("[findConsistenceAndMergePath] two consisPathArr not consistence arr1: %v, arr2:%v\n", consisPathArr1, consisPathArr2)
		}
	}

	var consisFlag bool = true
	if idx1 <= idx3 {
		i := 0
		j := idx3 - idx1
		for i < len(consisPathArr1) && j < len(consisPathArr2) {
			if consisPathArr1[i] != consisPathArr2[j] {
				consisFlag = false
				break
			}
			i++
			j++
		}
		if consisFlag == false {
			return mergePathArr
		} else {
			if i < len(consisPathArr1) {
				mergePathArr = append(consisPathArr2, consisPathArr1[i:]...)
			} else {
				mergePathArr = consisPathArr2
			}
		}
	} else {
		i := idx1 - idx3
		j := 0
		for i < len(consisPathArr1) && j < len(consisPathArr2) {
			if consisPathArr1[i] != consisPathArr2[j] {
				consisFlag = false
				break
			}
			i++
			j++
		}
		if consisFlag == false {
			return mergePathArr
		} else {
			if j < len(consisPathArr2) {
				mergePathArr = append(consisPathArr1, consisPathArr2[j:]...)
			} else {
				mergePathArr = consisPathArr1
			}
		}
	}
	return mergePathArr
}

func SubstituteEdgeID(nodesArr []constructdbg.DBGNode, nID, srcEID, dstEID constructdbg.DBG_MAX_INT) bool {
	var success bool
	for i := 0; i < bnt.BaseTypeNum; i++ {
		if nodesArr[nID].EdgeIDIncoming[i] == srcEID {
			nodesArr[nID].EdgeIDIncoming[i] = dstEID
			success = true
			break
		}
		if nodesArr[nID].EdgeIDOutcoming[i] == srcEID {
			nodesArr[nID].EdgeIDOutcoming[i] = dstEID
			success = true
			break
		}
	}
	return success
}

func GetDirection(node constructdbg.DBGNode, eID constructdbg.DBG_MAX_INT) (direction uint8) {
	if eID == 0 {
		log.Fatalf("[GetDirection] not allow found eid:%d\n", eID)
	}
	for i := 0; i < bnt.BaseTypeNum; i++ {
		if eID == node.EdgeIDIncoming[i] {
			direction = constructdbg.BACKWARD
			break
		} else if eID == node.EdgeIDOutcoming[i] {
			direction = constructdbg.FORWARD
			break
		}
	}
	if direction != constructdbg.BACKWARD && direction != constructdbg.FORWARD {
		log.Fatalf("[GetDirection] direction error\n")
	}

	return direction
}

func connectNodeID(e1, e2 constructdbg.DBGEdge) constructdbg.DBG_MAX_INT {
	if e1.StartNID == e2.StartNID || e1.StartNID == e2.EndNID {
		return e1.StartNID
	} else if e1.EndNID == e2.StartNID || e1.EndNID == e2.EndNID {
		return e1.EndNID
	} else {
		log.Fatalf("[connectNodeID] Not found connect node ID\n")
	}
	return 0
}

func ResetProcessFlag(edgesArr []constructdbg.DBGEdge) {
	for i, _ := range edgesArr {
		edgesArr[i].ResetProcessFlag()
	}
}

func DeleteEdgeID(nodesArr []constructdbg.DBGNode, nID, eID constructdbg.DBG_MAX_INT) (success bool) {
	for i := 0; i < bnt.BaseTypeNum; i++ {
		if nodesArr[nID].EdgeIDIncoming[i] == eID {
			nodesArr[nID].EdgeIDIncoming[i] = 0
			success = true
			break
		}
		if nodesArr[nID].EdgeIDOutcoming[i] == eID {
			nodesArr[nID].EdgeIDOutcoming[i] = 0
			success = true
			break
		}
	}

	return success
}

func IsInDBG_MAX_INTArr(arr []constructdbg.DBG_MAX_INT, eID constructdbg.DBG_MAX_INT) bool {

	for _, id := range arr {
		if id == eID {
			return true
		}
	}

	return false
}

func CheckConsistenceArr(arr1, arr2 []constructdbg.DBG_MAX_INT, eID constructdbg.DBG_MAX_INT) (mp []constructdbg.DBG_MAX_INT, ok bool) {
	m1 := IndexEID(arr1, eID)
	m2 := IndexEID(arr2, eID)
	if m1 < 0 || m2 < 0 {
		log.Fatalf("[IsConsistenceArr] Not found center position\n")
	}
	i, j := m1-1, m2-1
	for i >= 0 && j >= 0 {
		if arr1[i] != arr2[j] {
			return
		}
		i--
		j--
	}
	if i == 0 {
		mp = arr2[:m2+1]
	} else {
		mp = arr1[:m1+1]
	}
	i, j = m1+1, m2+1
	for i < len(arr1) && j < len(arr2) {
		if arr1[i] != arr2[j] {
			return
		}
		i++
		j++
	}
	if i == len(arr1) {
		mp = append(mp, arr2[m2+1:]...)
	} else {
		mp = append(mp, arr1[m1+1:]...)
	}
	ok = true
	return mp, ok
}

func FoundAllPath(edge constructdbg.DBGEdge, eIDArr []constructdbg.DBG_MAX_INT) (path constructdbg.Path, num int) {
	for _, eID := range eIDArr {
		var mergePath constructdbg.Path
		for _, p := range edge.PathMat {
			if IsInDBG_MAX_INTArr(p.IDArr, eID) {
				if len(mergePath.IDArr) == 0 {
					num++
					mergePath = p
				} else {
					if mp, ok := CheckConsistenceArr(mergePath.IDArr, p.IDArr, eID); ok {

						mergePath.Freq += p.Freq
						mergePath.IDArr = mp
					} else {
						num++
						break
					}
				}
			}
		}
		if len(mergePath.IDArr) > 0 {
			path = mergePath
		}
		if num > 1 {
			break
		}
	}

	return path, num
}

func GetUniqueOne(src, dst []constructdbg.DBG_MAX_INT) (unique constructdbg.DBG_MAX_INT) {
	for _, sID := range src {
		for _, dID := range dst {
			if sID == dID {
				unique = sID
				break
			}
		}
		if unique != 0 {
			break
		}
	}

	if unique == 0 {
		log.Fatalf("[GetUniqueOne] not found unique ID\n")
	}

	return unique
}

func FindMaxUnqiuePath(edgesArr []constructdbg.DBGEdge, nodesArr []constructdbg.DBGNode) {

	consistenceA := make([][]constructdbg.DBG_MAX_INT, len(edgesArr))

	for _, e := range edgesArr {
		if e.ID == 0 || len(e.Utg.Ks) < MIN_EDGE_LEN_EXTEND {
			continue
		}

		var right, left []constructdbg.DBG_MAX_INT
		// extend right
		if e.EndNID > 0 && IsDirectionUniqueEdge(e, nodesArr[e.EndNID]) {
			right = append(right, e.ID)
			reID := GetNextEID(e.ID, nodesArr[e.EndNID])
			if IsInDBG_MAX_INTArr(right, reID) == false {
				right = append(right, reID)
				ne := edgesArr[reID]
				var node constructdbg.DBGNode
				if ne.StartNID == e.EndNID {
					node = nodesArr[ne.EndNID]
				} else {
					node = nodesArr[ne.StartNID]
				}
				for node.ID > 0 {
					eIDArr := GetNextEArr(ne.ID, node)
					if len(eIDArr) == 0 {
						break
					} else if len(eIDArr) == 1 {
						if IsInDBG_MAX_INTArr(right, eIDArr[0]) {
							break
						}
						right = append(right, eIDArr[0])
						ne = edgesArr[eIDArr[0]]
						if ne.StartNID == node.ID {
							node = nodesArr[ne.EndNID]
						} else {
							node = nodesArr[ne.StartNID]
						}
					} else { // len(eIDArr) > 1
						var addedFlag bool
						backLen := len(edgesArr[right[len(right)-1]].Utg.Ks) - Kmerlen
						for i := len(right) - 2; i >= 0 && backLen < MAX_INSERT_SIZE; i-- {
							edge := edgesArr[right[i]]
							path, num := FoundAllPath(edge, eIDArr)
							if num == 1 {
								if path.Freq >= MIN_ALLOW_PATH_FREQ {
									uniqueID := GetUniqueOne(path.IDArr, eIDArr)
									cp := IndexEID(path.IDArr, uniqueID)
									if cp > 0 {
										if _, ok := CheckConsistenceArr(right, path.IDArr, path.IDArr[cp-1]); ok == false {
											break

										}
									}
									if IsInDBG_MAX_INTArr(right, uniqueID) {
										break
									}
									right = append(right, uniqueID)
									addedFlag = true
									ne = edgesArr[uniqueID]
									if ne.StartNID == node.ID {
										node = nodesArr[ne.EndNID]
									} else {
										node = nodesArr[ne.StartNID]
									}
									break
								}
							}
							backLen += len(edge.Utg.Ks) - Kmerlen
						}
						if addedFlag == false {
							break
						}
					}
				}
			}
		}
		// extend left
		if e.StartNID > 0 && IsDirectionUniqueEdge(e, nodesArr[e.StartNID]) {
			left = append(left, e.ID)
			leID := GetNextEID(e.ID, nodesArr[e.StartNID])
			if IsInDBG_MAX_INTArr(left, leID) == false {
				left = append(left, leID)
				ne := edgesArr[leID]
				var node constructdbg.DBGNode
				if ne.StartNID == e.StartNID {
					node = nodesArr[ne.EndNID]
				} else {
					node = nodesArr[ne.StartNID]
				}

				for node.ID > 0 {
					eIDArr := GetNextEArr(ne.ID, node)
					if len(eIDArr) == 0 {
						break
					} else if len(eIDArr) == 1 {
						if IsInDBG_MAX_INTArr(left, eIDArr[0]) {
							break
						}
						left = append(left, eIDArr[0])
						ne = edgesArr[eIDArr[0]]
						if ne.StartNID == node.ID {
							node = nodesArr[ne.EndNID]
						} else {
							node = nodesArr[ne.StartNID]
						}
					} else { // len(eIDArr) > 1
						var addedFlag bool
						backLen := len(edgesArr[left[len(left)-1]].Utg.Ks) - Kmerlen
						for i := len(left) - 2; i >= 0 && backLen < MAX_INSERT_SIZE; i-- {
							edge := edgesArr[left[i]]
							path, num := FoundAllPath(edge, eIDArr)
							if num == 1 {
								if path.Freq >= MIN_ALLOW_PATH_FREQ {
									uniqueID := GetUniqueOne(path.IDArr, eIDArr)
									cp := IndexEID(path.IDArr, uniqueID)
									if cp > 0 {
										if _, ok := CheckConsistenceArr(left, path.IDArr, path.IDArr[cp-1]); ok == false {
											break
										}
									}
									if IsInDBG_MAX_INTArr(left, uniqueID) {
										break
									}
									left = append(left, uniqueID)
									addedFlag = true
									ne = edgesArr[uniqueID]
									if ne.StartNID == node.ID {
										node = nodesArr[ne.EndNID]
									} else {
										node = nodesArr[ne.StartNID]
									}
									break
								}
							}
							backLen += len(edge.Utg.Ks) - Kmerlen
						}
						if addedFlag == false {
							break
						}
					}
				}
			}
		}

		// merge left and right extend Array
		var ma []constructdbg.DBG_MAX_INT
		if len(left) >= MIN_PATH_LEN {
			ma = GetReverseDBG_MAX_INTArr(left)
		}
		if len(right) >= MIN_PATH_LEN {
			if len(ma) > 0 {
				ma = append(ma, right[1:]...)
			} else {
				ma = right
			}
		}
		if len(ma) >= MIN_PATH_LEN {
			consistenceA[e.ID] = ma
			fmt.Printf("[FindMaxUnqiuePath] e.ID: %d, ma: %v\n", e.ID, ma)
		}
	}

	/*// merge all cross edge paths to the PathMat
	for _, e := range edgesArr {
		if e.ID == 0 || len(e.Utg.Ks) >= MAX_SHORT_READ_LEN {
			continue
		}
		if IsUniqueEdge(e, nodesArr) == false {
			continue
		}
		var eID constructdbg.DBG_MAX_INT
		var node constructdbg.DBGNode
		leID := GetNextEID(e.ID, nodesArr[e.StartNID])
		_, lnum := MergePath(e.PathMat, leID)
		reID := GetNextEID(e.ID, nodesArr[e.EndNID])
		_, rnum := MergePath(e.PathMat, reID)
		if lnum != 1 || rnum != 1 {
			continue
		}
		eID = reID
		node = nodesArr[e.EndNID]

		remainLen := MAX_SHORT_READ_LEN - len(e.Utg.Ks)
		stk := list.New()
		var pci PathCrossInfo
		pci.EdgeID, pci.Node, pci.RemainLen = eID, node, remainLen
		stk.PushBack(pci)
		for stk.Len() > 0 {
			ele := stk.Back()
			stk.Remove(ele)
			pci = ele.Value.(PathCrossInfo)
			if pci.RemainLen <= 0 {
				log.Fatalf("[FindMaxUnqiuePath] RemainLen smaller than zero\n")
			}
			ne := edgesArr[pci.EdgeID]
			for _, path := range ne.PathMat {
				j := IndexEID(path.IDArr, e.ID)
				if MIN_PATH_LEN-1 <= j && j < len(path.IDArr)-1 {
					slc := GetReverseDBG_MAX_INTArr(path.IDArr[:j+1])
					edgesArr[e.ID].InsertPathToEdge(slc, path.Freq)
					if len(path.IDArr)-j >= MIN_PATH_LEN {
						edgesArr[e.ID].InsertPathToEdge(path.IDArr[j:], path.Freq)
					}
				}
			}

			pci.RemainLen -= (len(ne.Utg.Ks) - Kmerlen + 1)
			if pci.RemainLen <= 0 {
				continue
			}

			// add to the stack next edges
			if ne.StartNID == pci.Node.ID {
				node = nodesArr[ne.EndNID]
			} else {
				node = nodesArr[ne.StartNID]
			}
			direction := GetDirection(node, pci.EdgeID)
			for j := 0; j < bnt.BaseTypeNum; j++ {
				if direction == constructdbg.FORWARD && node.EdgeIDIncoming[j] > 0 {
					var npci PathCrossInfo
					npci.EdgeID = node.EdgeIDIncoming[j]
					npci.Node = node
					npci.RemainLen = pci.RemainLen
					stk.PushBack(npci)
				} else if direction == constructdbg.BACKWARD && node.EdgeIDOutcoming[j] > 0 {
					var npci PathCrossInfo
					npci.EdgeID = node.EdgeIDOutcoming[j]
					npci.Node = node
					npci.RemainLen = pci.RemainLen
					stk.PushBack(npci)
				}
			}
		}
	}

	// fmt.Printf("[FindMaxUnqiuePath]<DeBug>edgesArr[382].PathMat: %v\n", edgesArr[382].PathMat)
	// find edge max unique path
	consistenceA := make([][]constructdbg.DBG_MAX_INT, len(edgesArr))
	{
		for _, e := range edgesArr {
			if e.ID == 0 {
				continue
			}
			if IsUniqueEdge(e, nodesArr) == false {
				continue
			}

			beID := GetNextEID(e.ID, nodesArr[e.StartNID])
			buniquePath, bnum := MergePath(e.PathMat, beID)
			feID := GetNextEID(e.ID, nodesArr[e.EndNID])
			funiquePath, fnum := MergePath(e.PathMat, feID)
			if bnum > 1 || fnum > 1 {
				continue
			}

			if bnum == 1 {
				tp := GetReverseDBG_MAX_INTArr(buniquePath.IDArr)
				consistenceA[e.ID] = tp
				edgesArr[e.ID].SetUniqueFlag()
			}
			if fnum == 1 {
				if len(consistenceA[e.ID]) > 0 {
					consistenceA[e.ID] = append(consistenceA[e.ID], funiquePath.IDArr[1:]...)
					if IsCyclePath(consistenceA[e.ID]) {
						consistenceA[e.ID] = nil
					}

				} else {
					consistenceA[e.ID] = funiquePath.IDArr
					edgesArr[e.ID].SetUniqueFlag()
				}
			}
		}
	}
	// fmt.Printf("[FindMaxUnqiuePath]edgesArr[353].PathMat: %v\n", consistenceA[353])

	// extend edges max unique path
	for i, e := range edgesArr {
		if e.GetProcessFlag() > 0 {
			continue
		}
		edgesArr[i].SetProcessFlag()

		if len(consistenceA[i]) < MIN_PATH_LEN {
			continue
		}
		// Debug
		if len(consistenceA[i]) >= 5 {
			fmt.Printf("[FindMaxUnqiuePath] consisteance[%d] Path: %v\n", i, consistenceA[i])
		}
		j := IndexEID(consistenceA[i], e.ID)
		if j < 0 {
			log.Fatalf("[FindMaxUnqiuePath] not index Edge ID: %d\n", i)
		}
		// extend left region
		if j >= MIN_PATH_LEN-1 {
			extendPathLen := j
			se := e
			for x := extendPathLen - 1; x >= 0; x-- {
				ne := edgesArr[consistenceA[i][x]]
				if len(consistenceA[ne.ID]) < MIN_PATH_LEN || IsUniqueEdge(ne, nodesArr) == false {
					continue
				}
				y := IndexEID(consistenceA[ne.ID], ne.ID)
				z := IndexEID(consistenceA[ne.ID], se.ID)
				if y >= 0 && z >= 0 {
					if y >= z {
						consistenceA[ne.ID] = GetReverseDBG_MAX_INTArr(consistenceA[ne.ID])
					}
					consistenceMergePath := findConsistenceAndMergePath(consistenceA[i], consistenceA[ne.ID], se.ID, ne.ID)
					if len(consistenceMergePath) >= len(consistenceA[i]) {
						if IsCyclePath(consistenceMergePath) {
							break
						}
						consistenceA[i] = consistenceMergePath
						x = IndexEID(consistenceA[i], ne.ID)
						edgesArr[ne.ID].SetProcessFlag()
						consistenceA[ne.ID] = nil
						se = ne
					} else {
						consistenceA[ne.ID] = nil
					}
				}
			}
		}
		if len(consistenceA[i])-j >= MIN_PATH_LEN { // extend right region
			se := e
			for x := j + 1; x < len(consistenceA[i]); x++ {
				ne := edgesArr[consistenceA[i][x]]
				if len(consistenceA[ne.ID]) < MIN_PATH_LEN || IsUniqueEdge(ne, nodesArr) {
					continue
				}
				y := IndexEID(consistenceA[ne.ID], ne.ID)
				z := IndexEID(consistenceA[ne.ID], se.ID)
				if y >= 0 && z >= 0 {
					if z >= y {
						consistenceA[ne.ID] = GetReverseDBG_MAX_INTArr(consistenceA[ne.ID])
					}
					consistenceMergePath := findConsistenceAndMergePath(consistenceA[i], consistenceA[ne.ID], se.ID, ne.ID)
					if len(consistenceMergePath) >= len(consistenceA[i]) {
						if IsCyclePath(consistenceMergePath) {
							break
						}
						consistenceA[i] = consistenceMergePath
						x = IndexEID(consistenceA[i], ne.ID)
						edgesArr[ne.ID].SetProcessFlag()
						consistenceA[ne.ID] = nil
						se = ne
					} else {
						consistenceA[ne.ID] = nil
					}
				}
			}
		}
	}

	// merge unique path
	ResetProcessFlag(edgesArr)
	{
		for i, _ := range edgesArr {
			if edgesArr[i].GetProcessFlag() > 0 || edgesArr[i].GetDeleteFlag() > 0 || len(consistenceA[i]) < MIN_PATH_LEN {
				continue
			}
			x := IndexUniqueEdge(consistenceA[i], edgesArr, constructdbg.FORWARD)
			y := IndexUniqueEdge(consistenceA[i], edgesArr, constructdbg.BACKWARD)
			if x < y && y-x+1 >= MIN_PATH_LEN {
				// tmp1 := consistenceA[i][0]
				// tmp2 := consistenceA[i][len(consistenceA[i])-1]
				// fmt.Printf("[FindMaxUnqiuePath]i: %d, merge path: %v, x:%d, y:%d\npathMat[%d]: %v\npathMat[%d]: %v\n", i, consistenceA[i], x, y, tmp1, edgesArr[tmp1].PathMat, tmp2, edgesArr[tmp2].PathMat)
				// e2 := edgesArr[consistenceA[i][x+1]]
				// if e1.StartNID == e2.StartNID || e1.StartNID == e2.EndNID {
				// 	constructdbg.RCEdge(edgesArr, e1.ID)
				// }
				if IsNodeCyclePath(consistenceA[i][x:y+1], edgesArr) {
					continue
				}
				e1 := edgesArr[consistenceA[i][x]]
				fmt.Printf("[FindMaxUnqiuePath] MergePath: %v\n", consistenceA[i][x:y+1])
				edgesArr[e1.ID].SetProcessFlag()
				for j := x + 1; j <= y; j++ {
					e2 := edgesArr[consistenceA[i][j]]
					edgesArr[e2.ID].SetProcessFlag()
					nID := connectNodeID(e1, e2)
					fmt.Printf("[FindMaxUnqiuePath] node: %v\ne1:%v\ne2:%v\n", nodesArr[nID], e1, e2)
					if GetDirection(nodesArr[nID], e2.ID) == constructdbg.FORWARD {
						if e1.StartNID == nID {
							constructdbg.RCEdge(edgesArr, e1.ID)
						}
						if e2.EndNID == nID {
							constructdbg.RCEdge(edgesArr, e2.ID)
						}
						e1 = edgesArr[e1.ID]
						e2 = edgesArr[e2.ID]
						DeleteEdgeID(nodesArr, e1.EndNID, e1.ID)
						constructdbg.ConcatEdges(edgesArr, e1.ID, e2.ID, e1.ID)
						// if SubstituteEdgeID(nodesArr, e1.EndNID, e1.ID, 0) == false {
						// 	log.Fatalf("[FindMaxUnqiuePath]1 SubstituteEdgeID failed")
						// }
						// edgesArr[e1.ID].EndNID = e2.EndNID
						if e2.GetUniqueFlag() > 0 {
							edgesArr[e2.ID].SetDeleteFlag()
							DeleteEdgeID(nodesArr, e2.StartNID, e2.ID)
							if SubstituteEdgeID(nodesArr, e2.EndNID, e2.ID, e1.ID) == false {
								log.Fatalf("[FindMaxUnqiuePath]1 SubstituteEdgeID failed")
							}
						}
						e1 = edgesArr[e1.ID]
					} else { // == constructdbg.FORWARD
						if e1.EndNID == nID {
							constructdbg.RCEdge(edgesArr, e1.ID)
						}
						if e2.StartNID == nID {
							constructdbg.RCEdge(edgesArr, e2.ID)
						}
						e1 = edgesArr[e1.ID]
						e2 = edgesArr[e2.ID]
						DeleteEdgeID(nodesArr, e1.StartNID, e1.ID)
						constructdbg.ConcatEdges(edgesArr, e2.ID, e1.ID, e1.ID)
						if e2.GetUniqueFlag() > 0 {
							edgesArr[e2.ID].SetDeleteFlag()
							DeleteEdgeID(nodesArr, e2.EndNID, e2.ID)
							if SubstituteEdgeID(nodesArr, e2.StartNID, e2.ID, e1.ID) == false {
								log.Fatalf("[FindMaxUnqiuePath]1 SubstituteEdgeID failed")
							}
						}
						e1 = edgesArr[e1.ID]
					}
				}
			}
		}
	} */

	/*for i, e := range edgesArr {
		// edge has been processed and  Unique
		if e.GetProcessFlag() > 0 {
			continue
		}
		if IsUniqueEdge(e, nodesArr) == false {
			continue
		}
		// constructdbg.FORWARD extend
		feID := GetNextEID(e.ID, nodesArr[e.EndNID])
		uniquePath, num := MergePath(e.PathMat, feID)
		if num == 1 && uniquePath.Freq >= MIN_PATH_FREQ {
			for i := 1; i < len(uniquePath.Path); i++ {
				edge := edgesArr[uniquePath.Path[i]]
				// check path consistence
				if IsConsistencePath(uniquePath.Path, i, edge.PathMat) == false {
					uniquePath.Path = uniquePath.Path[:i]
					break
				}
				if IsUniqueEdge(edge, nodesArr) {
					eid := GetNextEID(edge.ID, nodesArr[edge.EndNID])
					if eid == uniquePath.Path[i-1] {
						eid = GetNextEID(edge.ID, nodesArr[edge.StartNID])
					}
					up, n := MergePath(edge.PathMat)
				}
			}
		}
	}*/
}

func FindMaxLongReadsUnqiuePath(edgesArr []constructdbg.DBGEdge, nodesArr []constructdbg.DBGNode) {
	// merge all cross edge paths to the PathMat
	for _, e := range edgesArr {
		if e.ID == 0 || len(e.Utg.Ks) >= MAX_LONG_READ_LEN {
			continue
		}
		if IsUniqueEdge(e, nodesArr) == false {
			continue
		}
		// leID := GetNextEID(e.ID, nodesArr[e.StartNID])
		// _, lnum := MergePath(e.PathMat, leID)
		reID := GetNextEID(e.ID, nodesArr[e.EndNID])
		// _, rnum := MergePath(e.PathMat, reID)
		// if lnum != 1 || rnum != 1 {
		// 	continue
		// }
		eID := reID
		node := nodesArr[e.EndNID]

		remainLen := MAX_LONG_READ_LEN - len(e.Utg.Ks)
		stk := list.New()
		var pci PathCrossInfo
		pci.EdgeID, pci.Node, pci.RemainLen = eID, node, remainLen
		stk.PushBack(pci)
		for stk.Len() > 0 {
			ele := stk.Back()
			stk.Remove(ele)
			pci = ele.Value.(PathCrossInfo)
			if pci.RemainLen <= 0 {
				log.Fatalf("[FindMaxUnqiuePath] RemainLen smaller than zero\n")
			}
			ne := edgesArr[pci.EdgeID]
			for _, path := range ne.PathMat {
				j := IndexEID(path.IDArr, e.ID)
				if MIN_PATH_LEN-1 <= j && j < len(path.IDArr)-1 {
					slc := GetReverseDBG_MAX_INTArr(path.IDArr[:j+1])
					edgesArr[e.ID].InsertPathToEdge(slc, path.Freq)
					if len(path.IDArr)-j >= MIN_PATH_LEN {
						edgesArr[e.ID].InsertPathToEdge(path.IDArr[j:], path.Freq)
					}
				}
			}

			pci.RemainLen -= (len(ne.Utg.Ks) - Kmerlen + 1)
			if pci.RemainLen <= 0 {
				continue
			}

			// add to the stack next edges
			if ne.StartNID == pci.Node.ID {
				node = nodesArr[ne.EndNID]
			} else {
				node = nodesArr[ne.StartNID]
			}
			direction := GetDirection(node, pci.EdgeID)
			for j := 0; j < bnt.BaseTypeNum; j++ {
				if direction == constructdbg.FORWARD && node.EdgeIDIncoming[j] > 0 {
					var npci PathCrossInfo
					npci.EdgeID = node.EdgeIDIncoming[j]
					npci.Node = node
					npci.RemainLen = pci.RemainLen
					stk.PushBack(npci)
				} else if direction == constructdbg.BACKWARD && node.EdgeIDOutcoming[j] > 0 {
					var npci PathCrossInfo
					npci.EdgeID = node.EdgeIDOutcoming[j]
					npci.Node = node
					npci.RemainLen = pci.RemainLen
					stk.PushBack(npci)
				}
			}
		}
	}

	// fmt.Printf("[FindMaxUnqiuePath]<DeBug>edgesArr[382].PathMat: %v\n", edgesArr[382].PathMat)
	// find edge max unique path
	consistenceA := make([][]constructdbg.DBG_MAX_INT, len(edgesArr))
	{
		for _, e := range edgesArr {
			if e.ID == 0 {
				continue
			}
			if IsUniqueEdge(e, nodesArr) == false {
				continue
			}

			beID := GetNextEID(e.ID, nodesArr[e.StartNID])
			buniquePath, bnum := MergePath(e.PathMat, beID)
			feID := GetNextEID(e.ID, nodesArr[e.EndNID])
			funiquePath, fnum := MergePath(e.PathMat, feID)
			if bnum > 1 || fnum > 1 {
				continue
			}

			if bnum == 1 {
				tp := GetReverseDBG_MAX_INTArr(buniquePath.IDArr)
				consistenceA[e.ID] = tp
				edgesArr[e.ID].SetUniqueFlag()
			}
			if fnum == 1 {
				if len(consistenceA[e.ID]) > 0 {
					consistenceA[e.ID] = append(consistenceA[e.ID], funiquePath.IDArr[1:]...)
					if IsCyclePath(consistenceA[e.ID]) {
						consistenceA[e.ID] = nil
					}

				} else {
					consistenceA[e.ID] = funiquePath.IDArr
					edgesArr[e.ID].SetUniqueFlag()
				}
			}
		}
	}
	// fmt.Printf("[FindMaxUnqiuePath]edgesArr[353].PathMat: %v\n", consistenceA[353])

	// extend edges max unique path
	for i, e := range edgesArr {
		if e.GetProcessFlag() > 0 {
			continue
		}
		edgesArr[i].SetProcessFlag()

		if len(consistenceA[i]) < MIN_PATH_LEN {
			continue
		}
		// Debug
		if len(consistenceA[i]) >= 5 {
			fmt.Printf("[FindMaxUnqiuePath] consisteance[%d] Path: %v\n", i, consistenceA[i])
		}
		j := IndexEID(consistenceA[i], e.ID)
		if j < 0 {
			log.Fatalf("[FindMaxUnqiuePath] not index Edge ID: %d\n", i)
		}
		// extend left region
		if j >= MIN_PATH_LEN-1 {
			extendPathLen := j
			se := e
			for x := extendPathLen - 1; x >= 0; x-- {
				ne := edgesArr[consistenceA[i][x]]
				if len(consistenceA[ne.ID]) < MIN_PATH_LEN || IsUniqueEdge(ne, nodesArr) == false {
					continue
				}
				y := IndexEID(consistenceA[ne.ID], ne.ID)
				z := IndexEID(consistenceA[ne.ID], se.ID)
				if y >= 0 && z >= 0 {
					if y >= z {
						consistenceA[ne.ID] = GetReverseDBG_MAX_INTArr(consistenceA[ne.ID])
					}
					consistenceMergePath := findConsistenceAndMergePath(consistenceA[i], consistenceA[ne.ID], se.ID, ne.ID)
					if len(consistenceMergePath) >= len(consistenceA[i]) {
						if IsCyclePath(consistenceMergePath) {
							break
						}
						consistenceA[i] = consistenceMergePath
						x = IndexEID(consistenceA[i], ne.ID)
						edgesArr[ne.ID].SetProcessFlag()
						consistenceA[ne.ID] = nil
						se = ne
					} else {
						consistenceA[ne.ID] = nil
					}
				}
			}
		}
		if len(consistenceA[i])-j >= MIN_PATH_LEN { // extend right region
			se := e
			for x := j + 1; x < len(consistenceA[i]); x++ {
				ne := edgesArr[consistenceA[i][x]]
				if len(consistenceA[ne.ID]) < MIN_PATH_LEN || IsUniqueEdge(ne, nodesArr) {
					continue
				}
				y := IndexEID(consistenceA[ne.ID], ne.ID)
				z := IndexEID(consistenceA[ne.ID], se.ID)
				if y >= 0 && z >= 0 {
					if z >= y {
						consistenceA[ne.ID] = GetReverseDBG_MAX_INTArr(consistenceA[ne.ID])
					}
					consistenceMergePath := findConsistenceAndMergePath(consistenceA[i], consistenceA[ne.ID], se.ID, ne.ID)
					if len(consistenceMergePath) >= len(consistenceA[i]) {
						if IsCyclePath(consistenceMergePath) {
							break
						}
						consistenceA[i] = consistenceMergePath
						x = IndexEID(consistenceA[i], ne.ID)
						edgesArr[ne.ID].SetProcessFlag()
						consistenceA[ne.ID] = nil
						se = ne
					} else {
						consistenceA[ne.ID] = nil
					}
				}
			}
		}
	}

	// merge unique path
	ResetProcessFlag(edgesArr)
	{
		for i, _ := range edgesArr {
			if edgesArr[i].GetProcessFlag() > 0 || edgesArr[i].GetDeleteFlag() > 0 || len(consistenceA[i]) < MIN_PATH_LEN {
				continue
			}
			x := IndexUniqueEdge(consistenceA[i], edgesArr, constructdbg.FORWARD)
			y := IndexUniqueEdge(consistenceA[i], edgesArr, constructdbg.BACKWARD)
			if x < y && y-x+1 >= MIN_PATH_LEN {
				// tmp1 := consistenceA[i][0]
				// tmp2 := consistenceA[i][len(consistenceA[i])-1]
				// fmt.Printf("[FindMaxUnqiuePath]i: %d, merge path: %v, x:%d, y:%d\npathMat[%d]: %v\npathMat[%d]: %v\n", i, consistenceA[i], x, y, tmp1, edgesArr[tmp1].PathMat, tmp2, edgesArr[tmp2].PathMat)
				// e2 := edgesArr[consistenceA[i][x+1]]
				// if e1.StartNID == e2.StartNID || e1.StartNID == e2.EndNID {
				// 	constructdbg.RCEdge(edgesArr, e1.ID)
				// }
				if IsNodeCyclePath(consistenceA[i][x:y+1], edgesArr) {
					continue
				}
				e1 := edgesArr[consistenceA[i][x]]
				fmt.Printf("[FindMaxUnqiuePath] MergePath: %v\n", consistenceA[i][x:y+1])
				edgesArr[e1.ID].SetProcessFlag()
				for j := x + 1; j <= y; j++ {
					e2 := edgesArr[consistenceA[i][j]]
					edgesArr[e2.ID].SetProcessFlag()
					nID := connectNodeID(e1, e2)
					fmt.Printf("[FindMaxUnqiuePath] node: %v\ne1:%v\ne2:%v\n", nodesArr[nID], e1, e2)
					if GetDirection(nodesArr[nID], e2.ID) == constructdbg.FORWARD {
						if e1.StartNID == nID {
							constructdbg.RCEdge(edgesArr, e1.ID)
						}
						if e2.EndNID == nID {
							constructdbg.RCEdge(edgesArr, e2.ID)
						}
						e1 = edgesArr[e1.ID]
						e2 = edgesArr[e2.ID]
						DeleteEdgeID(nodesArr, e1.EndNID, e1.ID)
						constructdbg.ConcatEdges(edgesArr, e1.ID, e2.ID, e1.ID)
						// if SubstituteEdgeID(nodesArr, e1.EndNID, e1.ID, 0) == false {
						// 	log.Fatalf("[FindMaxUnqiuePath]1 SubstituteEdgeID failed")
						// }
						// edgesArr[e1.ID].EndNID = e2.EndNID
						if e2.GetUniqueFlag() > 0 {
							edgesArr[e2.ID].SetDeleteFlag()
							DeleteEdgeID(nodesArr, e2.StartNID, e2.ID)
							if SubstituteEdgeID(nodesArr, e2.EndNID, e2.ID, e1.ID) == false {
								log.Fatalf("[FindMaxUnqiuePath]1 SubstituteEdgeID failed")
							}
						}
						e1 = edgesArr[e1.ID]
					} else { // == constructdbg.FORWARD
						if e1.EndNID == nID {
							constructdbg.RCEdge(edgesArr, e1.ID)
						}
						if e2.StartNID == nID {
							constructdbg.RCEdge(edgesArr, e2.ID)
						}
						e1 = edgesArr[e1.ID]
						e2 = edgesArr[e2.ID]
						DeleteEdgeID(nodesArr, e1.StartNID, e1.ID)
						constructdbg.ConcatEdges(edgesArr, e2.ID, e1.ID, e1.ID)
						if e2.GetUniqueFlag() > 0 {
							edgesArr[e2.ID].SetDeleteFlag()
							DeleteEdgeID(nodesArr, e2.EndNID, e2.ID)
							if SubstituteEdgeID(nodesArr, e2.StartNID, e2.ID, e1.ID) == false {
								log.Fatalf("[FindMaxUnqiuePath]1 SubstituteEdgeID failed")
							}
						}
						e1 = edgesArr[e1.ID]
					}
				}
			}
		}
	}

	/*for i, e := range edgesArr {
		// edge has been processed and  Unique
		if e.GetProcessFlag() > 0 {
			continue
		}
		if IsUniqueEdge(e, nodesArr) == false {
			continue
		}
		// constructdbg.FORWARD extend
		feID := GetNextEID(e.ID, nodesArr[e.EndNID])
		uniquePath, num := MergePath(e.PathMat, feID)
		if num == 1 && uniquePath.Freq >= MIN_PATH_FREQ {
			for i := 1; i < len(uniquePath.Path); i++ {
				edge := edgesArr[uniquePath.Path[i]]
				// check path consistence
				if IsConsistencePath(uniquePath.Path, i, edge.PathMat) == false {
					uniquePath.Path = uniquePath.Path[:i]
					break
				}
				if IsUniqueEdge(edge, nodesArr) {
					eid := GetNextEID(edge.ID, nodesArr[edge.EndNID])
					if eid == uniquePath.Path[i-1] {
						eid = GetNextEID(edge.ID, nodesArr[edge.StartNID])
					}
					up, n := MergePath(edge.PathMat)
				}
			}
		}
	}*/
}

func IsComingInNode(node constructdbg.DBGNode, eID constructdbg.DBG_MAX_INT) bool {
	for i := 0; i < bnt.BaseTypeNum; i++ {
		if node.EdgeIDIncoming[i] == eID {
			return true
		}
		if node.EdgeIDOutcoming[i] == eID {
			return true
		}
	}
	return false
}

func CleanDBG(edgesArr []constructdbg.DBGEdge, nodesArr []constructdbg.DBGNode) {
	deleteNum := 0
	for i, e := range edgesArr {
		if e.GetDeleteFlag() > 0 {
			edgesArr[i] = edgesArr[0]
			deleteNum++
		} else {
			if e.StartNID > 0 {
				if IsComingInNode(nodesArr[e.StartNID], e.ID) == false {
					log.Fatalf("[CleanDBG] edge ID not include Node coming arr\nnode: %v\nedge: %v\n", nodesArr[e.StartNID], e)
					// fmt.Printf("[CleanDBG] edge ID not include Node coming arr\nnode: %v\nedge: %v\n", nodesArr[e.StartNID], e)
				}
			}
			if e.EndNID > 0 {
				if IsComingInNode(nodesArr[e.EndNID], e.ID) == false {
					log.Fatalf("[CleanDBG] edge ID not include Node coming arr\nnode: %v\nedge: %v\n", nodesArr[e.StartNID], e)
					// fmt.Printf("[CleanDBG] edge ID not include Node coming arr\nnode: %v\nedge: %v\n", nodesArr[e.EndNID], e)
				}
			}
			edgesArr[i].Flag = 0
		}
	}

	// check and clean node
	for i, n := range nodesArr {
		if n.ID > 0 {
			for j := 0; j < bnt.BaseTypeNum; j++ {
				if n.EdgeIDIncoming[j] > 0 {
					eID := n.EdgeIDIncoming[j]
					if edgesArr[eID].StartNID != n.ID && edgesArr[eID].EndNID != n.ID {
						nodesArr[i].EdgeIDIncoming[j] = 0
					}
				}
				if n.EdgeIDOutcoming[j] > 0 {
					eID := n.EdgeIDOutcoming[j]
					if edgesArr[eID].StartNID != n.ID && edgesArr[eID].EndNID != n.ID {
						nodesArr[i].EdgeIDOutcoming[j] = 0
					}
				}
			}
		}
	}

	fmt.Printf("[CleanDBG] delete edges number is : %d\n", deleteNum)
}

func GraphvizDBG(nodesArr []constructdbg.DBGNode, edgesArr []constructdbg.DBGEdge, graphfn string) {
	// create a new graph
	g := gographviz.NewGraph()
	g.SetName("G")
	g.SetDir(true)
	g.SetStrict(false)
	for _, v := range nodesArr {
		attr := gographviz.NewAttrs()
		attr.Add("color", "Green")
		attr.Add("shape", "record")
		var labels string
		labels = "{" + strconv.Itoa(int(v.EdgeIDIncoming[0])) + "|" + strconv.Itoa(int(v.EdgeIDIncoming[1])) + "|" + strconv.Itoa(int(v.EdgeIDIncoming[2])) + "|" + strconv.Itoa(int(v.EdgeIDIncoming[3])) + "}|" + strconv.Itoa(int(v.ID)) + "|{" + strconv.Itoa(int(v.EdgeIDOutcoming[0])) + "|" + strconv.Itoa(int(v.EdgeIDOutcoming[1])) + "|" + strconv.Itoa(int(v.EdgeIDOutcoming[2])) + "|" + strconv.Itoa(int(v.EdgeIDOutcoming[3])) + "}"
		attr.Add("label", labels)
		g.AddNode("G", strconv.Itoa(int(v.ID)), attr)
	}
	g.AddNode("G", "0", nil)

	for i := 1; i < len(edgesArr); i++ {
		e := edgesArr[i]
		if e.ID == 0 || e.GetDeleteFlag() > 0 {
			continue
		}
		attr := gographviz.NewAttrs()
		attr.Add("color", "Blue")
		labels := "ID:" + strconv.Itoa(int(e.ID)) + " len:" + strconv.Itoa(len(e.Utg.Ks))
		attr.Add("label", labels)
		g.AddEdge(strconv.Itoa(int(e.StartNID)), strconv.Itoa(int(e.EndNID)), true, attr)
	}
	// output := graph.String()
	gfp, err := os.Create(graphfn)
	if err != nil {
		log.Fatalf("[GraphvizDBG] Create file: %s failed, err: %v\n", graphfn, err)
	}
	defer gfp.Close()
	gfp.WriteString(g.String())
}

func FSpath(c cli.Command) {
	k := c.Parent().Flag("K").String()
	var err error = nil
	Kmerlen, err = strconv.Atoi(k)
	constructdbg.Kmerlen = Kmerlen
	if err != nil {
		log.Fatalf("[Fpath] argument: %s set error: %v\n", k, err)
	}
	numCPU, err := strconv.Atoi(c.Parent().Flag("t").String())
	if err != nil {
		log.Fatalf("[Fpath] argument: 't' set error: %v\n", err)
	}
	prefix := c.Parent().Flag("p").String()
	// read nodes file and transform to array mode, for more quickly access
	smfyNodesfn := prefix + ".nodes.smfy.mmap"
	nodeMap := constructdbg.NodeMapMmapReader(smfyNodesfn)
	nodesStatfn := prefix + ".nodes.stat"
	nodesSize := constructdbg.NodesStatReader(nodesStatfn)
	nodesArr := make([]constructdbg.DBGNode, nodesSize)
	constructdbg.NodeMap2NodeArr(nodeMap, nodesArr)
	nodeMap = nil
	// Restore edges info
	edgesStatfn := prefix + ".edges.stat"
	edgesSize := constructdbg.EdgesStatReader(edgesStatfn)
	edgesArr := make([]constructdbg.DBGEdge, edgesSize)
	edgesfn := prefix + ".edges.smfy.fq"
	constructdbg.LoadEdgesfqFromFn(edgesfn, edgesArr)

	bamfn := prefix + ".bam"
	rc := make(chan []sam.Record, numCPU*2)
	wc := make(chan []constructdbg.DBG_MAX_INT, numCPU*2)
	// defer rc.Close()
	go GetSamRecord(bamfn, rc, numCPU)

	for i := 0; i < numCPU; i++ {
		// go paraFindShortMappingPath(rc, wc, edgesArr, nodesArr)
	}

	WriteShortPathToDBG(wc, edgesArr, numCPU)

	// Find Max unique path   and merge neighbour edges
	// fmt.Printf("[FSpath] edgesArr[76]: %v\n", edgesArr[76])
	FindMaxUnqiuePath(edgesArr, nodesArr)
	CleanDBG(edgesArr, nodesArr)
	// simplify DBG
	// SmfyDBG(edgesArr, nodesArr)
	graphfn := prefix + ".ShortPath.dot"
	GraphvizDBG(nodesArr, edgesArr, graphfn)
	// Write to files
	edgesfn = prefix + ".edges.ShortPath.fq"
	constructdbg.StoreEdgesToFn(edgesfn, edgesArr, true)
	// constructdbg.StoreEdgesToFn(edgesfn, edgesArr, false)
	nodesfn := prefix + ".nodes.ShortPath.Arr"
	constructdbg.NodesArrWriter(nodesArr, nodesfn)
}

func Convert2LA(fields []string, RefIDMapArr []constructdbg.DBG_MAX_INT) (la LA) {
	id, err := strconv.Atoi(fields[0])
	if err != nil {
		log.Fatalf("[Convert2LA] convert %s error\n", fields[0])
	}
	// fmt.Printf("[Convert2LA] id: %d\n", id)
	la.RefID = RefIDMapArr[id]
	id, err = strconv.Atoi(fields[1])
	if err != nil {
		log.Fatalf("[Convert2LA] convert %s error\n", fields[1])
	}
	la.QuyID = constructdbg.DBG_MAX_INT(id)
	la.AlgnLen, err = strconv.Atoi(fields[2])
	if err != nil {
		log.Fatalf("[Convert2LA] convert %s error\n", fields[2])
	}
	if la.AlgnLen < 0 {
		la.AlgnLen = 0 - la.AlgnLen
	}
	la.Idty, err = strconv.ParseFloat(fields[3], 64)
	if err != nil {
		log.Fatalf("[Convert2LA] convert %s error\n", fields[3])
	}
	la.Idty /= 100
	la.RefB, err = strconv.Atoi(fields[5])
	if err != nil {
		log.Fatalf("[Convert2LA] convert %s error\n", fields[5])
	}
	la.RefE, err = strconv.Atoi(fields[6])
	if err != nil {
		log.Fatalf("[Convert2LA] convert %s error\n", fields[6])
	}
	la.RefLen, err = strconv.Atoi(fields[7])
	if err != nil {
		log.Fatalf("[Convert2LA] convert %s error\n", fields[7])
	}
	la.QuyB, err = strconv.Atoi(fields[9])
	if err != nil {
		log.Fatalf("[Convert2LA] convert %s error\n", fields[9])
	}
	la.QuyE, err = strconv.Atoi(fields[10])
	if err != nil {
		log.Fatalf("[Convert2LA] convert %s error\n", fields[10])
	}
	la.QuyLen, err = strconv.Atoi(fields[11])
	if err != nil {
		log.Fatalf("[Convert2LA] convert %s error\n", fields[11])
	}

	return la
}

func InsertIDArr(RefIDArr []constructdbg.DBG_MAX_INT, ID constructdbg.DBG_MAX_INT) []constructdbg.DBG_MAX_INT {
	for _, e := range RefIDArr {
		if ID == e {
			return RefIDArr
		}
	}
	RefIDArr = append(RefIDArr, ID)
	return RefIDArr
}

func Round(n float64) float64 {
	if n < 0 {
		return math.Ceil(n - 0.5)

	}
	return math.Floor(n + 0.5)
}

func GetLARecord(lafn string, RefIDMapArr []constructdbg.DBG_MAX_INT, lac chan []LA, numCPU int) {
	fp, err := os.Open(lafn)
	if err != nil {
		log.Fatalf("[GetSamRecord] open file: %s failed, err: %v\n", lafn, err)
	}
	defer fp.Close()
	lafp := bufio.NewReader(fp)
	var laArr []LA
	var RefIDArr []constructdbg.DBG_MAX_INT
	parseNum := 0
	totalNum := 0
	for {
		line, err := lafp.ReadString('\n')
		if err != nil {
			if err == io.EOF {
				break
			} else {
				log.Fatalf("[GetLARecord] read %s file error\n", lafn)
			}
		}
		fields := strings.Fields(line)
		// fmt.Printf("[GetLARecord]fields:%v, err: %v\n", fields, err)
		la := Convert2LA(fields, RefIDMapArr)
		// fmt.Printf("[GetLARecord]la:%v\n", la)
		if len(laArr) > 0 {
			if laArr[0].QuyID != la.QuyID {
				if len(RefIDArr) >= MIN_PATH_LEN {
					lac <- laArr
					parseNum++
				}
				laArr = nil
				RefIDArr = RefIDArr[:0]
				totalNum++
			}
		}
		la.Diff = int(Round(float64((la.QuyE-la.QuyB)+(la.RefE-la.RefB)) * (1 - la.Idty) / 2))
		laArr = append(laArr, la)
		RefIDArr = InsertIDArr(RefIDArr, la.RefID)
	}
	if len(laArr) >= MIN_PATH_LEN && len(RefIDArr) >= MIN_PATH_LEN {
		lac <- laArr
		parseNum++
	}
	fmt.Printf("[GetLARecord] total num is : %d\tparse num is : %d\n", totalNum, parseNum)
	// send terminal signals
	for i := 0; i < numCPU; i++ {
		var nilArr []LA
		lac <- nilArr
	}
}

func MaxIndexEID(laArr []LA, start int) int {
	end := -1
	for i := start; i < len(laArr); i++ {
		if laArr[i].RefID == laArr[start].RefID {
			end = i
		} else {
			return end
		}
	}

	return end
}

func findMaxLengthMappingLA(laArr []LA) (maxLa LA) {
	maxLa = laArr[0]
	for i := 1; i < len(laArr); i++ {
		if laArr[i].AlgnLen > maxLa.AlgnLen {
			maxLa = laArr[i]
		}
	}

	return maxLa
}

func findNearLA(laArr []LA, src LA, direction uint8) (la LA) {
	for _, e := range laArr {
		if direction == constructdbg.FORWARD {
			if e.RefB >= src.RefE {
				if la.RefID == 0 {
					la = e
				} else {
					if e.RefB < la.RefB {
						la = e
					} else if e.RefB == la.RefB && e.RefE > la.RefE {
						la = e
					}
				}
			}
		} else {
			if e.RefE <= src.RefB {
				if la.RefID == 0 {
					la = e
				} else {
					if e.RefE > la.RefE {
						la = e
					} else if e.RefE == la.RefE && e.RefB < la.RefB {
						la = e
					}
				}
			}
		}
	}

	return la
}

func IsInGapRegArr(gapRegArr []GapRegion, gp GapRegion) bool {
	for _, g := range gapRegArr {
		if gp.Begin <= g.Begin && gp.End >= g.End {
			if ((gp.End - gp.Begin) - (g.End - g.Begin)) <= (g.End - g.Begin) {
				return true
			} else {
				return false
			}
		} else {
			if gp.End < g.Begin {
				return false
			}
		}
	}
	return false
}

func MergeLA(laArr []LA, gapRegArr []GapRegion, initLa LA) (la LA) {
	// extend left partition
	for {
		left := findNearLA(laArr, initLa, constructdbg.BACKWARD)
		if left.RefID == 0 {
			break
		}
		var gp GapRegion
		if left.QuyE < initLa.QuyB {
			gp.Begin, gp.End = left.QuyE, initLa.QuyB
		} else {
			gp.Begin, gp.End = initLa.QuyE, left.QuyB
		}
		if IsInGapRegArr(gapRegArr, gp) == false {
			break
		}

		if initLa.QuyB > left.QuyE {
			initLa.Diff += left.Diff + (initLa.QuyB - left.QuyE)
			initLa.QuyB = left.QuyB
		} else {
			initLa.Diff += left.Diff + (left.QuyB - initLa.QuyE)
			initLa.QuyE = left.QuyE
		}
		initLa.AlgnLen += left.AlgnLen
		initLa.RefB = left.RefB
	}

	// extend right parition
	for {
		right := findNearLA(laArr, initLa, constructdbg.FORWARD)
		if right.RefID == 0 {
			break
		}
		// initLa.QuyE must small than right.QuyB
		var gp GapRegion
		if right.QuyB > initLa.QuyE {
			gp.Begin, gp.End = initLa.QuyE, right.QuyB
		} else {
			gp.Begin, gp.End = right.QuyE, initLa.QuyB
		}
		if IsInGapRegArr(gapRegArr, gp) == false {
			break
		}

		if initLa.QuyE < right.QuyB {
			initLa.Diff += right.Diff + (right.QuyB - initLa.QuyE)
			initLa.QuyE = right.QuyE
		} else {
			initLa.Diff += right.Diff + (initLa.QuyB - right.QuyE)
			initLa.QuyB = right.QuyB
		}
		initLa.AlgnLen += right.AlgnLen
		initLa.RefE = right.RefE
	}
	la = initLa
	la.Idty = 1 - float64(la.Diff)/float64((la.QuyE-la.QuyB))
	return la
}

func GetMaxLenLA(la1, la2 LA) LA {
	n1 := la1.QuyE - la1.QuyB - la1.Diff
	n2 := la2.QuyE - la2.QuyB - la2.Diff
	if n1 > n2 {
		return la1
	} else {
		return la2
	}
}

func GetRobustLA(la1, la2 LA) LA {
	n1 := la1.QuyE - la1.QuyB - la1.Diff
	n2 := la2.QuyE - la2.QuyB - la2.Diff
	if n1 > n2 {
		if (n1 - n2) < n1/10 {
			if la2.Idty > la1.Idty {
				fmt.Printf("[GetRobustLA] la1: %v\nla2: %v\n", la1, la2)
			}
		}
		return la1
	} else {
		if (n2 - n1) < n2/10 {
			if la1.Idty > la2.Idty {
				fmt.Printf("[GetRobustLA] la1: %v\nla2: %v\n", la1, la2)
			}
		}
		return la2
	}
}

func GetMaxMappingLength(laArr []LA, gapRegArr []GapRegion) (la LA) {
	la = laArr[0]
	for i := 0; i < len(laArr); i++ {
		j := MaxIndexEID(laArr, i)
		maxLa := findMaxLengthMappingLA(laArr[i : j+1])
		if j-i > 0 {
			maxLa = MergeLA(laArr[i:j+1], gapRegArr, maxLa)
		}
		la = GetMaxLenLA(la, maxLa)
		// fmt.Printf("[paraFindLongMappingPath] la: %v\n", la)
	}

	return la
}

func Min(x, y int) int {
	if x < y {
		return x
	}
	return y
}

func GetLARegion(laArr []LA, eID constructdbg.DBG_MAX_INT) []LA {
	b, e := -1, -1
	for i, la := range laArr {
		if la.RefID == eID {
			if b == -1 {
				b, e = i, i+1
			} else {
				e = i + 1
			}
		} else {
			if b >= 0 {
				break
			}
		}
	}
	// if b == -1 {
	// 	log.Fatalf("[GetLARegion]b:%d, e:%d, not found edgeID:%d\n", b, e, eID)
	// }
	if b == -1 {
		return laArr[0:0]
	} else {
		return laArr[b:e]
	}
}

func GetMinDistanceFocus(la LA, focus int) (min int) {
	if la.QuyB > focus {
		min = la.QuyB - focus
	} else if la.QuyE < focus {
		min = focus - la.QuyE
	} else {
		min = Min(focus-la.QuyB, la.QuyE-focus)
	}

	return min
}

func GetNearstLA(focus int, edgeLaArr []LA) (la LA) {
	la = edgeLaArr[0]
	min := GetMinDistanceFocus(la, focus)
	for i := 1; i < len(edgeLaArr); i++ {
		n := GetMinDistanceFocus(edgeLaArr[i], focus)
		if n < min {
			la = edgeLaArr[i]
			min = n
		}
	}
	return la
}

func findEdgeLA(eID constructdbg.DBG_MAX_INT, laArr []LA, gapRegArr []GapRegion) (la LA) {
	edgeLaArr := GetLARegion(laArr, eID)
	if len(edgeLaArr) == 1 {
		la = edgeLaArr[0]
	} else if len(edgeLaArr) > 1 {
		maxLa := findMaxLengthMappingLA(edgeLaArr)
		la = MergeLA(edgeLaArr, gapRegArr, maxLa)
	}
	// if la.RefE-la.RefB < Kmerlen {
	// 	var tmp LA
	// 	la = tmp
	// }
	return la
}

func findNextProEdge(node constructdbg.DBGNode, eID constructdbg.DBG_MAX_INT, laArr []LA, gapRegArr []GapRegion, la LA) (nLa LA, num int) {
	direction := GetDirection(node, eID)
	for i := 0; i < bnt.BaseTypeNum; i++ {
		if direction == constructdbg.FORWARD {
			if node.EdgeIDIncoming[i] > 0 {
				proLa := findEdgeLA(node.EdgeIDIncoming[i], laArr, gapRegArr)
				if proLa.RefID == 0 {
					break
				}
				if proLa.QuyE < la.QuyE {
					if proLa.QuyE-la.QuyB < Kmerlen/QUY_OVL_FACTOR {
						break
					}
				} else {
					if la.QuyE-proLa.QuyB < Kmerlen/QUY_OVL_FACTOR {
						break
					}
				}
				num++
				if nLa.RefID == 0 {
					nLa = proLa
				} else {
					fmt.Printf("[findNextProEdge] proLa:%v\nnLa:%v\n", proLa, nLa)
					if proLa.AlgnLen == nLa.AlgnLen {
						if proLa.Idty == nLa.Idty {
							var la LA
							nLa = la
							break
						} else if proLa.Idty > nLa.Idty {
							nLa = proLa
						}
					} else {
						if proLa.AlgnLen > nLa.AlgnLen {
							if proLa.Idty >= nLa.Idty {
								nLa = proLa
							} else {
								if proLa.AlgnLen > 2*nLa.AlgnLen {
									nLa = proLa
								} else {
									var la LA
									nLa = la
									break
								}
							}
						} else {
							if proLa.Idty > nLa.Idty && proLa.AlgnLen*2 > nLa.AlgnLen {
								var la LA
								nLa = la
								break
							}
						}
					}
				}
			}
		} else {
			if node.EdgeIDOutcoming[i] > 0 {
				proLa := findEdgeLA(node.EdgeIDOutcoming[i], laArr, gapRegArr)
				if proLa.RefID == 0 {
					break
				}
				if proLa.QuyE > la.QuyE {
					if la.QuyE-proLa.QuyB < Kmerlen/QUY_OVL_FACTOR {
						break
					}
				} else {
					if proLa.QuyE-la.QuyB < Kmerlen/QUY_OVL_FACTOR {
						break
					}
				}
				num++
				if nLa.RefID == 0 {
					nLa = proLa
				} else {
					fmt.Printf("[findNextProEdge] proLa:%v\nnLa:%v\n", proLa, nLa)
					if proLa.AlgnLen == nLa.AlgnLen {
						if proLa.Idty == nLa.Idty {
							var la LA
							nLa = la
							break
						} else if proLa.Idty > nLa.Idty {
							nLa = proLa
						}
					} else {
						if proLa.AlgnLen > nLa.AlgnLen {
							if proLa.Idty >= nLa.Idty {
								nLa = proLa
							} else {
								if proLa.AlgnLen > 2*nLa.AlgnLen {
									nLa = proLa
								} else {
									var la LA
									nLa = la
									break
								}
							}
						} else {
							if proLa.Idty > nLa.Idty && proLa.AlgnLen*2 > nLa.AlgnLen {
								var la LA
								nLa = la
								break
							}
						}
					}
				}
			}
		}
	}
	return nLa, num
}

func GetLAArrGapRegion(laArr []LA) []GapRegion {
	var gra1, gra2 []GapRegion
	var gp GapRegion
	gp.Begin, gp.End = 0, laArr[0].QuyLen
	gra1 = append(gra1, gp)
	gra2 = append(gra2, gp)
	// cycle fill gap region when encounter alignment
	for _, la := range laArr {
		gra2 = gra2[0:0]
		// fmt.Printf("[GetLAArrGapRegion]la:%v\n", la)
		for _, g := range gra1 {
			if g.End <= la.QuyB || g.Begin >= la.QuyE {
				gra2 = append(gra2, g)
			} else if g.Begin <= la.QuyB && la.QuyB < g.End {
				if g.Begin < la.QuyB {
					gp.Begin, gp.End = g.Begin, la.QuyB
					gra2 = append(gra2, gp)
				}
				if la.QuyE < g.End {
					gp.Begin, gp.End = la.QuyE, g.End
					gra2 = append(gra2, gp)
				} else {
					la.QuyB = g.End
				}

			} else if g.Begin < la.QuyE && la.QuyE <= g.End {
				if la.QuyE < g.End {
					gp.Begin, gp.End = la.QuyE, g.End
					gra2 = append(gra2, gp)
				}
			}
		}
		gra1, gra2 = gra2, gra1
		// fmt.Printf("[GetLAArrGapRegion]gra1:%v\n", gra1)
	}

	return gra1
}

func paraFindLongMappingPath(lac chan []LA, wc chan []constructdbg.DBG_MAX_INT, edgesArr []constructdbg.DBGEdge, nodesArr []constructdbg.DBGNode, RefIDMapArr []constructdbg.DBG_MAX_INT) {
	for {
		laArr := <-lac
		if len(laArr) == 0 {
			var guardPath []constructdbg.DBG_MAX_INT
			wc <- guardPath
			break
		}

		for i, te := range laArr {
			fmt.Printf("[paraFindLongMappingPath] laArr[%d]: %v\n", i, te)
		}
		gapRegArr := GetLAArrGapRegion(laArr)
		fmt.Printf("[paraFindLongMappingPath] gapRegArr: %v\n", gapRegArr)
		// var path []constructdbg.DBG_MAX_INT
		maxLa := GetMaxMappingLength(laArr, gapRegArr)
		fmt.Printf("[paraFindLongMappingPath] maxLa: %v\n", maxLa)
		edge := edgesArr[maxLa.RefID]
		// fmt.Printf("[paraFindLongMappingPath] RefLen: %d, edges Length: %d\n", maxLa.RefLen, len(edge.Utg.Ks))
		var left, right []constructdbg.DBG_MAX_INT
		// if maxLa.RefB < Min(len(edge.Utg.Ks)/FLANK_ALLOW_FACTOR, MAX_FLANK_ALLOW_LEN)+len(constructdbg.AdpaterSeq) {
		flankLen := Kmerlen/REF_OVL_FACTOR + len(constructdbg.AdpaterSeq)
		if maxLa.RefB < flankLen {
			la := maxLa
			eID := la.RefID
			nID := edgesArr[eID].StartNID

			fmt.Printf("[paraFindLongMappingPath] extend left region\n")
			for nID > 0 {
				//neID, b, e, rl
				nla, foundNum := findNextProEdge(nodesArr[nID], eID, laArr, gapRegArr, la)
				if nla.RefID == 0 {
					break
				}
				// fmt.Printf("[paraFindLongMappingPath] next La: %v\n", nla)
				nedge := edgesArr[nla.RefID]
				// flankLen := Min(len(nedge.Utg.Ks)/FLANK_ALLOW_FACTOR, MAX_FLANK_ALLOW_LEN) + len(constructdbg.AdpaterSeq)
				if edgesArr[nla.RefID].EndNID == nID {
					if foundNum > 1 && nla.RefE < len(nedge.Utg.Ks)-flankLen {
						break
					}
					if nla.RefID != eID && isInEdgesArr(left, nla.RefID) == false {
						left = append(left, nla.RefID)
						if foundNum > 1 && nla.RefB > flankLen {
							break
						}
						la = nla
						eID = nla.RefID
						nID = edgesArr[nla.RefID].StartNID
					} else {
						break
					}
				} else {
					if foundNum > 1 && nla.RefB > flankLen {
						break
					}
					if nla.RefID != eID && isInEdgesArr(left, nla.RefID) == false {
						left = append(left, nla.RefID)
						if foundNum > 1 && nla.RefE < len(nedge.Utg.Ks)-flankLen {
							break
						}
						la = nla
						eID = nla.RefID
						nID = edgesArr[nla.RefID].EndNID
					} else {
						break
					}
				}
			}
		}
		// extend rigth path
		if maxLa.RefE > len(edge.Utg.Ks)-flankLen {
			la := maxLa
			eID := la.RefID
			nID := edgesArr[eID].EndNID
			fmt.Printf("[paraFindLongMappingPath] extend right region\n")
			for nID > 0 {
				// neID, b, e, rl
				nla, foundNum := findNextProEdge(nodesArr[nID], eID, laArr, gapRegArr, la)
				if nla.RefID == 0 {
					break
				}
				// fmt.Printf("[paraFindLongMappingPath] next La: %v\n", nla)
				nedge := edgesArr[nla.RefID]
				if edgesArr[nla.RefID].StartNID == nID {
					if foundNum > 1 && nla.RefB > flankLen {
						break
					}
					if nla.RefID != eID && isInEdgesArr(right, nla.RefID) == false {
						right = append(right, nla.RefID)
						if foundNum > 1 && nla.RefE < len(nedge.Utg.Ks)-flankLen {
							break
						}
						eID = nla.RefID
						la = nla
						nID = edgesArr[nla.RefID].EndNID
					} else {
						break
					}
				} else {
					if foundNum > 1 && nla.RefE < len(nedge.Utg.Ks)-flankLen {
						break
					}
					if nla.RefID != eID && isInEdgesArr(right, nla.RefID) == false {
						right = append(right, nla.RefID)
						if foundNum > 1 && nla.RefB > flankLen {
							break
						}
						eID = nla.RefID
						la = nla
						nID = edgesArr[nla.RefID].StartNID
					} else {
						break
					}
				}
			}
		}

		// cat left and  right path
		var catArr []constructdbg.DBG_MAX_INT
		for i := len(left) - 1; i >= 0; i-- {
			catArr = append(catArr, left[i])
		}
		catArr = append(catArr, edge.ID)
		catArr = append(catArr, right...)
		fmt.Printf("[paraFindLongMappingPath] catArr: %v\n", catArr)
		if len(catArr) >= MIN_PATH_LEN {
			wc <- catArr
		}
	}
}

func GetOrderID(edgesArr []constructdbg.DBGEdge) (ReadIDMapArr []constructdbg.DBG_MAX_INT) {
	for _, e := range edgesArr {
		if e.ID > 0 {
			ReadIDMapArr = append(ReadIDMapArr, e.ID)
		}
	}
	return ReadIDMapArr
}
func WriteLongPathToDBG(wc chan []constructdbg.DBG_MAX_INT, edgesArr []constructdbg.DBGEdge, numCPU int) {
	WriteShortPathToDBG(wc, edgesArr, numCPU)
}

func checkDBG(edgesArr []constructdbg.DBGEdge, nodesArr []constructdbg.DBGNode) {

	// print pathMat
	for _, e := range edgesArr {
		if e.ID == 0 {
			continue
		}
		fmt.Printf("[checkDBG] edge ID : %d\n", e.ID)
		for _, p := range e.PathMat {
			fmt.Printf("[checkDBG] path: %v\n", p)
		}
	}
}

func FLpath(c cli.Command) {
	k := c.Parent().Flag("K").String()
	var err error = nil
	Kmerlen, err = strconv.Atoi(k)
	constructdbg.Kmerlen = Kmerlen
	if err != nil {
		log.Fatalf("[Fpath] argument: %s set error: %v\n", k, err)
	}
	numCPU, err := strconv.Atoi(c.Parent().Flag("t").String())
	if err != nil {
		log.Fatalf("[Fpath] argument: 't' set error: %v\n", err)
	}
	prefix := c.Parent().Flag("p").String()
	// read nodesArr file
	spNodesfn := prefix + ".nodes.ShortPath.Arr"
	nodesArr := constructdbg.NodesArrReader(spNodesfn)
	// Restore edges info
	edgesStatfn := prefix + ".edges.stat"
	edgesSize := constructdbg.EdgesStatReader(edgesStatfn)
	edgesArr := make([]constructdbg.DBGEdge, edgesSize)
	edgesfn := prefix + ".edges.ShortPath.fq"
	constructdbg.LoadEdgesfqFromFn(edgesfn, edgesArr)

	LongReadPathfn := prefix + ".LA"
	lac := make(chan []LA, numCPU*2)
	wc := make(chan []constructdbg.DBG_MAX_INT, numCPU*2)
	RefIDMapArr := GetOrderID(edgesArr)
	// defer rc.Close()
	go GetLARecord(LongReadPathfn, RefIDMapArr, lac, numCPU)

	for i := 0; i < numCPU; i++ {
		go paraFindLongMappingPath(lac, wc, edgesArr, nodesArr, RefIDMapArr)
	}

	WriteLongPathToDBG(wc, edgesArr, numCPU)

	// Find Max unique path   and merge neighbour edges
	// fmt.Printf("[FSpath] edgesArr[76]: %v\n", edgesArr[76])
	// SmfyDBG(edgesArr, nodesArr)
	checkDBG(edgesArr, nodesArr)
	FindMaxUnqiuePath(edgesArr, nodesArr)
	// simplify DBG
	graphfn := prefix + ".LongPath.dot"
	GraphvizDBG(nodesArr, edgesArr, graphfn)
	// Write to files
	edgesfn = prefix + ".edges.LongPath.fq"
	constructdbg.StoreEdgesToFn(edgesfn, edgesArr, false)
	// constructdbg.StoreEdgesToFn(edgesfn, edgesArr, false)
	nodesfn := prefix + ".nodes.LongPath.Arr"
	constructdbg.NodesArrWriter(nodesArr, nodesfn)
}

func Fpath(c cli.Command) {
	k := c.Parent().Flag("K").String()
	var err error = nil
	Kmerlen, err = strconv.Atoi(k)
	constructdbg.Kmerlen = Kmerlen
	if err != nil {
		log.Fatalf("[Fpath] argument: %s set error: %v\n", k, err)
	}
	numCPU, err := strconv.Atoi(c.Parent().Flag("t").String())
	if err != nil {
		log.Fatalf("[Fpath] argument: 't' set error: %v\n", err)
	}
	prefix := c.Parent().Flag("p").String()
	// read nodes file and transform to array mode, for more quickly access
	smfyNodesfn := prefix + ".nodes.smfy.mmap"
	nodeMap := constructdbg.NodeMapMmapReader(smfyNodesfn)
	nodesStatfn := prefix + ".nodes.stat"
	nodesSize := constructdbg.NodesStatReader(nodesStatfn)
	nodesArr := make([]constructdbg.DBGNode, nodesSize)
	constructdbg.NodeMap2NodeArr(nodeMap, nodesArr)
	nodeMap = nil
	// Restore edges info
	edgesStatfn := prefix + ".edges.stat"
	edgesSize := constructdbg.EdgesStatReader(edgesStatfn)
	edgesArr := make([]constructdbg.DBGEdge, edgesSize)
	edgesfn := prefix + ".edges.smfy.fq"
	constructdbg.LoadEdgesfqFromFn(edgesfn, edgesArr)

	fastmapfn := prefix + ".fastmap"
	rc := make(chan [2]FastMapRecord, numCPU*2)
	wc := make(chan []constructdbg.DBG_MAX_INT, numCPU*2)
	// defer rc.Close()
	go GetFastMapRecord(fastmapfn, rc, numCPU)

	for i := 0; i < numCPU; i++ {
		go paraFindShortMappingPath(rc, wc, edgesArr, nodesArr)
	}

	WriteShortPathToDBG(wc, edgesArr, numCPU)

	// Find Max unique path   and merge neighbour edges
	// fmt.Printf("[FSpath] edgesArr[76]: %v\n", edgesArr[76])
	FindMaxUnqiuePath(edgesArr, nodesArr)
	// CleanDBG(edgesArr, nodesArr)
	// // simplify DBG
	// // SmfyDBG(edgesArr, nodesArr)
	// graphfn := prefix + ".ShortPath.dot"
	// GraphvizDBG(nodesArr, edgesArr, graphfn)
	// // Write to files
	// edgesfn = prefix + ".edges.ShortPath.fq"
	// constructdbg.StoreEdgesToFn(edgesfn, edgesArr, true)
	// // constructdbg.StoreEdgesToFn(edgesfn, edgesArr, false)
	// nodesfn := prefix + ".nodes.ShortPath.Arr"
	// constructdbg.NodesArrWriter(nodesArr, nodesfn)
}
