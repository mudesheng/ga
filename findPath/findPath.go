package findPath

import (
	"bufio"
	"sort"
	// "compress/gzip"
	// "encoding/binary"
	// "encoding/gob"
	"container/list"
	"fmt"
	"ga/bnt"
	"reflect"
	// "ga/constructcf"
	"ga/constructdbg"
	// "ga/cuckoofilter"
	"io"
	"log"
	"os"
	// "unsafe"
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
	MAX_LONG_READ_LEN       = 20000
	MIN_LONG_REG_LEN        = 5000
	TAG_ZOOM_SIZE           = 100
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
	MATCH_SCORE             = 1
)

type PathCrossInfo struct {
	EdgeID    constructdbg.DBG_MAX_INT
	NID       constructdbg.DBG_MAX_INT
	RemainLen int
}

type IDLen struct {
	EID  constructdbg.DBG_MAX_INT
	ELen int32
}

type LA struct {
	RefID   constructdbg.DBG_MAX_INT
	QuyID   constructdbg.DBG_MAX_INT // query ID
	AlgnLen int
	Diff    int
	Minus   bool // flag if have been aligned Minus strand
	Flag    bool // Flag used set flag state by function
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

func GetLastMapRecord(lastfn string, rc chan [2]FastMapRecord, numCPU int) {
	fp, err := os.Open(lastfn)
	if err != nil {
		log.Fatalf("[GetLastRecord] open file: %s failed, err: %v\n", lastfn, err)
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
	for {
		r, err := buffp.ReadString('\n')
		if err != nil {
			break
		}
		// fmt.Printf("[GetLastMapRecord] r[0]:%c\n", r[0])
		if r[0] == '#' {
			if len(pairinfo[0].Minfo) > 0 && len(pairinfo[1].Minfo) > 0 && (len(pairinfo[0].Minfo)+len(pairinfo[1].Minfo) > 2 || pairinfo[0].Minfo[0].RefID != pairinfo[1].Minfo[0].RefID) {
				// fmt.Printf("[GetLastMapRecord] pairinfo:%v\n", pairinfo)
				rc <- pairinfo
			}
			pairinfo[0].Minfo = nil
			pairinfo[1].Minfo = nil
		} else {
			fields := strings.Split(strings.TrimSpace(r), "\t")
			// fmt.Printf("[GetLastMapRecord] fields:%v\n", fields)
			alnlen, err1 := strconv.Atoi(fields[3])
			score, err2 := strconv.Atoi(fields[0])
			if err1 != nil || err2 != nil {
				log.Fatalf("[GetLastMapRecord] r: %s , err: %v\n", r, err)
			}
			if alnlen < Kmerlen || score < MATCH_SCORE*alnlen {
				continue
			}
			// fmt.Printf("[GetLastMapRecord] alnlen:%d\tscore:%d\n", alnlen, score)
			var idx int
			if fields[6][len(fields[6])-1:] == "1" {
				idx = 0
			} else {
				idx = 1
			}
			pairinfo[idx].ReadID = fields[6]
			tmp, err3 := strconv.Atoi(fields[10])
			if err3 != nil {
				log.Fatalf("[GetLastMapRecord] r: %s , err: %v\n", r, err)
			}
			pairinfo[idx].Rlen = int32(tmp)
			var info MapInfo
			tmp, err3 = strconv.Atoi(fields[7])
			tmp1, err4 := strconv.Atoi(fields[8])
			if err3 != nil || err4 != nil {
				log.Fatalf("[GetLastMapRecord] r: %s , err3: %v, err4: %v\n", r, err3, err4)
			}
			// fmt.Printf("[GetLastMapRecord] fields[9]:%v\n", fields[9])
			// info.Start = tmp
			// info.End = info.Start + tmp1
			if fields[9] == "+" {
				info.Start = tmp
				info.End = info.Start + tmp1
			} else { // minus strand
				info.Start = int(pairinfo[idx].Rlen) - (tmp1 + tmp)
				info.End = int(pairinfo[idx].Rlen) - tmp
			}
			tmp, err3 = strconv.Atoi(fields[1])
			if err3 != nil {
				log.Fatalf("[GetLastMapRecord] r: %s , err: %v\n", r, err)
			}
			info.RefID = constructdbg.DBG_MAX_INT(tmp)
			tmp, err3 = strconv.Atoi(fields[2])
			if err3 != nil {
				log.Fatalf("[GetLastMapRecord] r: %s , err: %v\n", r, err)
			}
			info.RefStart = tmp
			pairinfo[idx].Minfo = append(pairinfo[idx].Minfo, info)
		}
	}
	// send terminal signals
	for i := 0; i < numCPU; i++ {
		var nilArr [2]FastMapRecord
		rc <- nilArr
	}
}

func GetLastMapRecordOne(lastfn string, rc chan FastMapRecord, numCPU int) {
	fp, err := os.Open(lastfn)
	if err != nil {
		log.Fatalf("[GetLastRecord] open file: %s failed, err: %v\n", lastfn, err)
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
	for {
		r, err := buffp.ReadString('\n')
		if err != nil {
			break
		}
		var info FastMapRecord
		// fmt.Printf("[GetLastMapRecord] r[0]:%c\n", r[0])
		if r[0] != '#' {
			fields := strings.Split(strings.TrimSpace(r), "\t")
			// fmt.Printf("[GetLastMapRecord] fields:%v\n", fields)
			alnlen, err1 := strconv.Atoi(fields[8])
			score, err2 := strconv.Atoi(fields[0])
			if err1 != nil || err2 != nil {
				log.Fatalf("[GetLastMapRecord] r: %s , err1: %v, er2: %v\n", r, err1, err2)
			}
			tmp, err3 := strconv.Atoi(fields[10])
			if err3 != nil {
				log.Fatalf("[GetLastMapRecord] r: %s , err3: %v\n", r, err3)
			}
			info.Rlen = int32(tmp)
			if alnlen < Kmerlen || score < MATCH_SCORE*int(info.Rlen) {
				continue
			}
			// fmt.Printf("[GetLastMapRecord] alnlen:%d\tscore:%d\n", alnlen, score)
			info.ReadID = fields[6]
			// if int(info.Minfo[0].Ref) > score*MATCH_SCORE {
			// 	continue
			// }
			var minfo MapInfo
			tmp, err3 = strconv.Atoi(fields[2])
			tmp1, err4 := strconv.Atoi(fields[3])
			if err3 != nil || err4 != nil {
				log.Fatalf("[GetLastMapRecord] r: %s , err3: %v, err4: %v\n", r, err3, err4)
			}

			// fmt.Printf("[GetLastMapRecord] fields[9]:%v\n", fields[9])
			minfo.Start = tmp
			minfo.End = minfo.Start + tmp1
			minfo.RefStart = tmp
			// if fields[9] == "+" {
			// 	minfo.Start = tmp
			// 	minfo.End = minfo.Start + tmp1
			// } else { // minus strand
			// 	minfo.Start = int(info.Rlen) - (tmp1 + tmp)
			// 	minfo.End = int(info.Rlen) - tmp
			// }
			// tmp, err3 = strconv.Atoi(fields[1])
			// if err3 != nil {
			// 	log.Fatalf("[GetLastMapRecord] r: %s , err3: %v\n", r, err3)
			// }
			// minfo.RefID = constructdbg.DBG_MAX_INT(tmp)
			info.Minfo = append(info.Minfo, minfo)
			rc <- info
		}
	}
	// send terminal signals
	for i := 0; i < numCPU; i++ {
		var nilArr FastMapRecord
		rc <- nilArr
	}
}

func computeCoverageSmfyEdge(prefix string) {
	numCPU := 1
	rc := make(chan FastMapRecord, numCPU*2)
	lastfn := prefix + ".mapRef.last"
	coverfn := prefix + ".edgeCoverGap"
	cfp, err := os.Create(coverfn)
	if err != nil {
		log.Fatalf("[GetLastRecord] open file: %s failed, err: %v\n", coverfn, err)
	}
	defer cfp.Close()
	go GetLastMapRecordOne(lastfn, rc, numCPU)
	refLen := 5000000
	high := 1000
	coverArr := make([]int, refLen)

	for {
		r := <-rc
		// fmt.Fprintf(os.Stderr, "eID: %v\n", r.ReadID)
		if len(r.Minfo) == 0 {
			break
		}
		for i := 0; i < int(r.Minfo[0].End-r.Minfo[0].Start); i++ {
			if i < Kmerlen-1 || i > int(r.Rlen)-Kmerlen {
				coverArr[r.Minfo[0].Start+i] += 1
			} else {
				coverArr[r.Minfo[0].Start+i] += 1
			}
		}
	}
	gapLen := 0
	notNode := false
	for _, v := range coverArr {
		// fmt.Fprintf(cfp, "%d\n", v)
		if v < 2 {
			notNode = true
			gapLen++
		} else if v < high {
			gapLen++
		} else {
			if notNode {
				if gapLen > 0 {
					fmt.Fprintf(cfp, "%d\n", gapLen)
				}
			} else {
				if gapLen > 0 && gapLen < Kmerlen-1 {
					fmt.Fprintf(cfp, "%d\n", gapLen)
				}
			}
			gapLen = 0
			notNode = false
		}
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

type MapInfoSlice []MapInfo

func (p MapInfoSlice) Len() int {
	return len(p)
}

func (p MapInfoSlice) Less(i, j int) bool {
	return p[i].Start < p[j].Start
}

func (p MapInfoSlice) Swap(i, j int) {
	p[i], p[j] = p[j], p[i]
}
func bubbleSort(data sort.Interface) {
	r := data.Len() - 1
	for i := 0; i < r; i++ {
		for j := r; j > i; j-- {
			if data.Less(j, j-1) {
				data.Swap(j, j-1)
			}
		}
	}
}

/*func AscendSortMinfo(Minfo []MapInfo) (sm []MapInfo) {
	for i := 0; i < len(Minfo); i++ {
		var info MapInfo
		info.Start = math.MaxInt32
		for j := 0; j < len(Minfo); j++ {
			if info.Start > Minfo[j] {
				info = Minfo[j]
			}
		}
	}
} */

func IsConnectEdges(miArr []MapInfo, edgesArr []constructdbg.DBGEdge, nodesArr []constructdbg.DBGNode) bool {
	if len(miArr) < 2 {
		return true
	}
	e1 := edgesArr[miArr[0].RefID]
	e2 := edgesArr[miArr[1].RefID]
	var nID constructdbg.DBG_MAX_INT
	if e1.StartNID == e2.StartNID || e1.EndNID == e2.StartNID {
		nID = e2.StartNID
	} else if e1.StartNID == e2.EndNID || e2.EndNID == e2.EndNID {
		nID = e2.EndNID
	} else {
		return false
	}

	for i := 2; i < len(miArr); i++ {
		te := edgesArr[miArr[i].RefID]
		if te.StartNID == nID {
			nID = te.EndNID
		} else if te.EndNID == nID {
			nID = te.StartNID
		} else {
			return false
		}
	}

	return true

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
		// sort by read alignment start position
		// if len(rArr[0].Minfo) > 2 {
		// 	fmt.Printf("[paraFindPacbioMappingPath]rArr[0]: %v\n", rArr[0])
		// }
		bubbleSort(MapInfoSlice(rArr[0].Minfo))
		// if len(rArr[0].Minfo) > 2 {
		// 	fmt.Printf("[paraFindPacbioMappingPath]rArr[0]: %v\n", rArr[0])
		// }
		bubbleSort(MapInfoSlice(rArr[1].Minfo))
		// if len(rArr[0].Minfo) > 1 && rArr[0].Minfo[0].Start == 0 && rArr[0].Minfo[1].Start == 0 {
		// fmt.Printf("[paraFindPacbioMappingPath]rArr[0].Minfo: %v\n", rArr[0].Minfo)
		// }
		// if len(rArr[1].Minfo) > 1 && rArr[1].Minfo[0].Start == 0 && rArr[1].Minfo[1].Start == 0 {
		// fmt.Printf("[paraFindPacbioMappingPath]rArr[1].Minfo: %v\n", rArr[1].Minfo)
		// }
		if IsConnectEdges(rArr[0].Minfo, edgesArr, nodesArr) == false || IsConnectEdges(rArr[1].Minfo, edgesArr, nodesArr) == false {
			continue
		}

		// write read1~read2 rrelative path
		for _, p1 := range rArr[0].Minfo {
			var path []constructdbg.DBG_MAX_INT
			path = append(path, p1.RefID)
			if len(rArr[1].Minfo) <= 2 && path[0] == rArr[1].Minfo[len(rArr[1].Minfo)-1].RefID {
				continue
			}
			j := len(rArr[1].Minfo) - 1
			if path[0] != rArr[1].Minfo[j].RefID {
				path = append(path, rArr[1].Minfo[j].RefID)
			}
			for j--; j >= 0; j-- {
				path = append(path, rArr[1].Minfo[j].RefID)
			}

			wc <- path
		}

		// write read2~read1 relative path
		for _, p2 := range rArr[1].Minfo {
			var path []constructdbg.DBG_MAX_INT
			path = append(path, p2.RefID)
			if len(rArr[0].Minfo) <= 2 && path[0] == rArr[0].Minfo[len(rArr[0].Minfo)-1].RefID {
				continue
			}
			j := len(rArr[0].Minfo) - 1
			if path[0] != rArr[0].Minfo[j].RefID {
				path = append(path, rArr[0].Minfo[j].RefID)
			}
			for j--; j >= 0; j-- {
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
		// if len(path) > 3 {
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
		log.Fatalf("[GetNextEArr] direction not set, eID:%d, node:%v\n", eID, node)
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

func GetOtherEArr(node constructdbg.DBGNode, eID constructdbg.DBG_MAX_INT) (eIDArr []constructdbg.DBG_MAX_INT) {
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
		log.Fatalf("[GetOtherEArrEID] direction not set\n")
	}
	for i := 0; i < bnt.BaseTypeNum; i++ {
		if direction == constructdbg.FORWARD {
			if node.EdgeIDOutcoming[i] > 0 && node.EdgeIDOutcoming[i] != eID {
				eIDArr = append(eIDArr, node.EdgeIDOutcoming[i])
			}
		} else {
			if node.EdgeIDIncoming[i] > 0 && node.EdgeIDIncoming[i] != eID {
				eIDArr = append(eIDArr, node.EdgeIDIncoming[i])
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

/*func FindMaxUnqiuePath(edgesArr []constructdbg.DBGEdge, nodesArr []constructdbg.DBGNode) [][]constructdbg.DBG_MAX_INT {

consistenceA := make([][]constructdbg.DBG_MAX_INT, len(edgesArr))

for _, e := range edgesArr {
	if e.ID == 0 || len(e.Utg.Ks) < Kmerlen*3/2 { // Kmerlen * 1.5
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

return consistenceA

// merge all cross edge paths to the PathMat
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
	}
} */

func IsInPathCrossInfoArr(arr []PathCrossInfo, p PathCrossInfo) bool {
	for _, ele := range arr {
		if ele == p {
			return true
		}
	}

	return false
}

func GetNeighbourEID(eID constructdbg.DBG_MAX_INT, nID constructdbg.DBG_MAX_INT, edgesArr []constructdbg.DBGEdge, nodesArr []constructdbg.DBGNode, max_insert_size int) (neighbourEIDArr []constructdbg.DBG_MAX_INT) {

	var stk []PathCrossInfo
	if nID > 0 {
		var pci PathCrossInfo
		pci.EdgeID, pci.NID, pci.RemainLen = eID, nID, max_insert_size
		stk = append(stk, pci)
	}

	for len(stk) > 0 {
		pci := stk[len(stk)-1]
		stk = stk[:len(stk)-1]
		// fmt.Printf("[GetNeighbourEID] pci: %v\n", pci)
		eIDArr := GetNextEArr(pci.EdgeID, nodesArr[pci.NID])
		for _, eID := range eIDArr {
			if IsInDBG_MAX_INTArr(neighbourEIDArr, eID) {
				continue
			}
			neighbourEIDArr = append(neighbourEIDArr, eID)
			if edgesArr[eID].StartNID == pci.NID {
				nID = edgesArr[eID].EndNID
			} else {
				nID = edgesArr[eID].StartNID
			}
			rlen := pci.RemainLen - (len(edgesArr[eID].Utg.Ks) - Kmerlen + 1)
			if nID > 0 && rlen > 0 {
				var p PathCrossInfo
				p.EdgeID = eID
				p.NID = nID
				p.RemainLen = rlen
				if IsInPathCrossInfoArr(stk, p) == false {
					stk = append(stk, p)
				}
			}
		}

	}

	return neighbourEIDArr
}

func IsIntersection(rightEIDArr, leftEIDArr []constructdbg.DBG_MAX_INT) bool {
	for _, eID := range rightEIDArr {
		if IsInDBG_MAX_INTArr(leftEIDArr, eID) {
			return true
		}
	}
	return false
}

func IsContained(subArr, allArr []constructdbg.DBG_MAX_INT) bool {
	for _, id := range subArr {
		if IsInDBG_MAX_INTArr(allArr, id) == false {
			return false
		}
	}

	return true
}

func IsInPathMat(pathMat []constructdbg.Path, eID constructdbg.DBG_MAX_INT) bool {
	for _, p := range pathMat {
		if IsInDBG_MAX_INTArr(p.IDArr, eID) {
			return true
		}
	}
	return false
}

// if has any other edges share same start(or end) node and start direction
func IsParaEdge(edge constructdbg.DBGEdge, nodesArr []constructdbg.DBGNode, pathMat []constructdbg.Path) bool {
	for i := 0; i < 2; i++ {
		var nID constructdbg.DBG_MAX_INT
		if i == 0 {
			nID = edge.StartNID
		} else {
			nID = edge.EndNID
		}
		if nID == 0 {
			continue
		}
		eIDArr := GetOtherEArr(nodesArr[nID], edge.ID)
		for _, eID := range eIDArr {
			if IsInPathMat(pathMat, eID) {
				return true
			}
		}
	}

	return false
}

func GetInterNodeID(edgeSrc, edgeDst constructdbg.DBGEdge) constructdbg.DBG_MAX_INT {
	if edgeSrc.StartNID > 0 {
		if edgeSrc.StartNID == edgeDst.StartNID || edgeSrc.StartNID == edgeDst.EndNID {
			return edgeSrc.StartNID
		}
	}
	if edgeSrc.EndNID > 0 {
		if edgeSrc.EndNID == edgeDst.StartNID || edgeSrc.EndNID == edgeDst.EndNID {
			return edgeSrc.EndNID
		}
	}

	// log.Fatalf("[GetInterNodeID] not found shared node\n")

	return 0
}

func GetMergePathArr(EpathMat []constructdbg.Path, edgesArr []constructdbg.DBGEdge, nodesArr []constructdbg.DBGNode) (mergePathArr []constructdbg.Path) {
	var pathMat []constructdbg.Path
	for _, pi := range EpathMat {
		var p constructdbg.Path
		p.IDArr = make([]constructdbg.DBG_MAX_INT, len(pi.IDArr))
		copy(p.IDArr, pi.IDArr)
		p.Freq = pi.Freq
		pathMat = append(pathMat, p)
	}

	for i := 0; i < len(pathMat); i++ {
		if pathMat[i].Freq == 0 {
			continue
		}
		p := pathMat[i]
		// check unique path
		var np constructdbg.Path
		np.IDArr = make([]constructdbg.DBG_MAX_INT, len(p.IDArr))
		copy(np.IDArr, p.IDArr)
		uniqueFlag := true
		for _, eID := range np.IDArr {
			if IsParaEdge(edgesArr[eID], nodesArr, pathMat) {
				uniqueFlag = false
			}
		}

		// extend right flank
		if uniqueFlag {
			eID := np.IDArr[len(np.IDArr)-1]
			var nID constructdbg.DBG_MAX_INT
			if len(np.IDArr) == 1 {
				nID = edgesArr[eID].EndNID
			} else {
				tmp := GetInterNodeID(edgesArr[np.IDArr[len(np.IDArr)-2]], edgesArr[eID])
				if tmp == 0 {
					log.Fatalf("[GetMergePathArr] not found shared node\n")
				}
				if tmp == edgesArr[eID].StartNID {
					nID = edgesArr[eID].EndNID
				} else {
					nID = edgesArr[eID].StartNID
				}
			}
			for nID > 0 && eID > 0 {
				// fmt.Printf("[GetMergePathArr] extend right nID: %d, eID:%v\n", nID, eID)
				eIDArr := GetNextEArr(eID, nodesArr[nID])
				num := 0
				for _, ne := range eIDArr {
					if IsInPathMat(pathMat, ne) {
						eID = ne
						np.IDArr = append(np.IDArr, ne)
						if IsCyclePath(np.IDArr) {
							num += 2
							break

						} else {
							num++
						}
					}
				}
				if num == 0 {
					break
				} else if num > 1 {
					uniqueFlag = false
					break
				} else {
					tmp := GetInterNodeID(edgesArr[np.IDArr[len(np.IDArr)-2]], edgesArr[eID])
					if tmp == 0 {
						log.Fatalf("[GetMergePathArr] not found shared node\n")
					}
					if tmp == edgesArr[eID].StartNID {
						nID = edgesArr[eID].EndNID
					} else {
						nID = edgesArr[eID].StartNID
					}
				}
			}
		}
		if uniqueFlag {
			// extend left flank
			eID := np.IDArr[0]
			var nID constructdbg.DBG_MAX_INT
			if len(np.IDArr) == 1 {
				nID = edgesArr[eID].StartNID
			} else {
				tmp := GetInterNodeID(edgesArr[np.IDArr[1]], edgesArr[eID])
				if tmp == 0 {
					log.Fatalf("[GetMergePathArr] not found shared node\n")
				}
				if tmp == edgesArr[eID].StartNID {
					nID = edgesArr[eID].EndNID
				} else {
					nID = edgesArr[eID].StartNID
				}
			}
			for nID > 0 && eID > 0 {
				// fmt.Printf("[GetMergePathArr] extend left nID: %d, eID:%v\n", nID, eID)
				eIDArr := GetNextEArr(eID, nodesArr[nID])
				num := 0
				for _, ne := range eIDArr {
					if IsInPathMat(pathMat, ne) {
						eID = ne
						var na []constructdbg.DBG_MAX_INT
						na = append(na, ne)
						na = append(na, np.IDArr...)
						np.IDArr = na
						if IsCyclePath(np.IDArr) {
							num += 2
							break
						} else {
							num++
						}
					}
				}
				if num == 0 {
					break
				} else if num > 1 {
					uniqueFlag = false
					break
				} else {
					tmp := GetInterNodeID(edgesArr[np.IDArr[1]], edgesArr[eID])
					if tmp == 0 {
						log.Fatalf("[GetMergePathArr] not found shared node\n")
					}
					if tmp == edgesArr[eID].StartNID {
						nID = edgesArr[eID].EndNID
					} else {
						nID = edgesArr[eID].StartNID
					}
				}
			}

			if uniqueFlag {
				mergePathArr = append(mergePathArr, np)
			}
		}

		// flag have processed path
		for j := i + 1; j < len(pathMat); j++ {
			if IsIntersection(pathMat[j].IDArr, np.IDArr) {
				pathMat[j].Freq = 0
			}
		}
	}

	return mergePathArr
}

func FindMaxUnqiuePath(edgesArr []constructdbg.DBGEdge, nodesArr []constructdbg.DBGNode) [][2][]constructdbg.DBG_MAX_INT {

	consistenceA := make([][2][]constructdbg.DBG_MAX_INT, len(edgesArr))

	for _, e := range edgesArr {
		if e.ID == 0 || len(e.Utg.Ks) < Kmerlen*3/2 { // Kmerlen * 1.5
			continue
		}

		if (e.StartNID == 0 || IsDirectionUniqueEdge(e, nodesArr[e.StartNID]) == false) && (e.EndNID == 0 || IsDirectionUniqueEdge(e, nodesArr[e.EndNID]) == false) {
			continue
		}

		mergePathArr := GetMergePathArr(e.PathMat, edgesArr, nodesArr)
		// fmt.Printf("[FindMaxUnqiuePath] eID: %d, pathMat:%v\n", e.ID, e.PathMat)
		if 0 < len(mergePathArr) && len(mergePathArr) <= 2 {
			rightEIDArr := GetNeighbourEID(e.ID, e.EndNID, edgesArr, nodesArr, MAX_INSERT_SIZE)
			leftEIDArr := GetNeighbourEID(e.ID, e.StartNID, edgesArr, nodesArr, MAX_INSERT_SIZE)
			if IsIntersection(rightEIDArr, leftEIDArr) == false {
				// fmt.Printf("[FindMaxUnqiuePath] mergePathArr:%v\n", mergePathArr)
				// fmt.Printf("[FindMaxUnqiuePath] rightEIDArr:%v\nleftEIDArr:%v\n", rightEIDArr, leftEIDArr)
				var flag uint8
				for _, pi := range mergePathArr {
					var found bool
					if flag != constructdbg.FORWARD {
						if IsContained(pi.IDArr, rightEIDArr) {
							consistenceA[e.ID][0] = pi.IDArr
							// fmt.Printf("[FindMaxUnqiuePath] pi:%v\n", pi)
							found = true
							flag = constructdbg.FORWARD
						}
					}

					if found == false && flag != constructdbg.BACKWARD {
						if IsContained(pi.IDArr, leftEIDArr) {
							consistenceA[e.ID][1] = pi.IDArr
							// fmt.Printf("[FindMaxUnqiuePath] pi:%v\n", pi)
							found = true
							flag = constructdbg.BACKWARD
						}
					}

					// if found == false {
					// 	log.Fatalf("[FindMaxUnqiuePath] not found consistence path\n")
					// }
				}
			}
		}
	}

	return consistenceA

	/*var right, left []constructdbg.DBG_MAX_INT
	// extend right
	if e.EndNID > 0 && IsDirectionUniqueEdge(e, nodesArr[e.EndNID]) {
		// right = append(right, e.ID)
		reID := GetNextEID(e.ID, nodesArr[e.EndNID])
		focusIdx := GetFocusEID(nodesArr[e.EndNID], reID, edgesArr, nodesArr, pathMat)

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
	} */

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
		pci.EdgeID, pci.NID, pci.RemainLen = eID, node.ID, remainLen
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
			if ne.StartNID == pci.NID {
				node = nodesArr[ne.EndNID]
			} else {
				node = nodesArr[ne.StartNID]
			}
			direction := GetDirection(node, pci.EdgeID)
			for j := 0; j < bnt.BaseTypeNum; j++ {
				if direction == constructdbg.FORWARD && node.EdgeIDIncoming[j] > 0 {
					var npci PathCrossInfo
					npci.EdgeID = node.EdgeIDIncoming[j]
					npci.NID = node.ID
					npci.RemainLen = pci.RemainLen
					stk.PushBack(npci)
				} else if direction == constructdbg.BACKWARD && node.EdgeIDOutcoming[j] > 0 {
					var npci PathCrossInfo
					npci.EdgeID = node.EdgeIDOutcoming[j]
					npci.NID = node.ID
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
	constructdbg.LoadEdgesfqFromFn(edgesfn, edgesArr, false)

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

	id, err = strconv.Atoi(fields[8])
	if err != nil {
		log.Fatalf("[Convert2LA] convert %s error\n", fields[7])
	}
	if id > 0 {
		la.Minus = true
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
		if la.AlgnLen < Kmerlen*2/3 {
			continue
		}
		// fmt.Printf("[GetLARecord]la:%v\n", la)
		if len(laArr) > 0 {
			if laArr[0].QuyID != la.QuyID {
				// if len(RefIDArr) >= MIN_PATH_LEN {
				lac <- laArr
				parseNum++
				// }
				laArr = nil
				RefIDArr = RefIDArr[:0]
				totalNum++
			}
		}
		la.Diff = int(Round(float64((la.QuyE-la.QuyB)+(la.RefE-la.RefB)) * (1 - la.Idty) / 2))
		la.Idty = 1 - float64(la.Diff)/float64(la.AlgnLen)
		// fmt.Fprintf(os.Stderr, "[GetLARecord] %v\n", la)
		laArr = append(laArr, la)
		RefIDArr = InsertIDArr(RefIDArr, la.RefID)
	}
	if len(laArr) >= MIN_PATH_LEN && len(RefIDArr) >= MIN_PATH_LEN {
		lac <- laArr
		parseNum++
	}
	fmt.Printf("[GetLARecord] total num is : %d\tparse LARecord num is : %d\n", totalNum, parseNum)
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
		if e.RefID != src.RefID {
			continue
		}
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

func IsCrossGapRegArr(gapRegArr []GapRegion, gp GapRegion) bool {
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

func IsInGapRegArr(gapRegArr []GapRegion, gp GapRegion) bool {
	for _, g := range gapRegArr {
		if gp.Begin >= g.Begin && gp.End <= g.End {
			return true
		} else if gp.Begin < g.Begin {
			return false
		}
	}
	return false
}

func MergeLA(laArr []LA, gapRegArr []GapRegion, initLa LA) (la LA) {
	// extend left partition
	for {
		left := findNearLA(laArr, initLa, constructdbg.BACKWARD)
		// fmt.Fprintf(os.Stderr, "[MergeLA]initLa: %v\nleft: %v\n", initLa, left)
		if left.RefID == 0 || left.Minus != initLa.Minus {
			break
		}
		var gp GapRegion
		if left.QuyE < initLa.QuyB {
			gp.Begin, gp.End = left.QuyE, initLa.QuyB
		} else {
			gp.Begin, gp.End = initLa.QuyE, left.QuyB
		}
		if IsCrossGapRegArr(gapRegArr, gp) == false {
			break
		}

		if initLa.QuyB > left.QuyE {
			// initLa.Diff += left.Diff + (initLa.QuyB - left.QuyE)
			initLa.QuyB = left.QuyB
		} else {
			// initLa.Diff += left.Diff + (left.QuyB - initLa.QuyE)
			initLa.QuyE = left.QuyE
		}
		initLa.AlgnLen += left.AlgnLen
		initLa.Diff += left.Diff
		initLa.RefB = left.RefB
	}

	// extend right parition
	for {
		right := findNearLA(laArr, initLa, constructdbg.FORWARD)
		// fmt.Fprintf(os.Stderr, "[MergeLA]right: %v\n", right)
		if right.RefID == 0 || right.Minus != initLa.Minus {
			break
		}
		// initLa.QuyE must small than right.QuyB
		var gp GapRegion
		if right.QuyB > initLa.QuyE {
			gp.Begin, gp.End = initLa.QuyE, right.QuyB
		} else {
			gp.Begin, gp.End = right.QuyE, initLa.QuyB
		}
		if IsCrossGapRegArr(gapRegArr, gp) == false {
			break
		}

		if initLa.QuyE < right.QuyB {
			// initLa.Diff += right.Diff + (right.QuyB - initLa.QuyE)
			initLa.QuyE = right.QuyE
		} else {
			// initLa.Diff += right.Diff + (initLa.QuyB - right.QuyE)
			initLa.QuyB = right.QuyB
		}
		initLa.AlgnLen += right.AlgnLen
		initLa.Diff += right.Diff
		initLa.RefE = right.RefE
	}
	la = initLa
	la.Idty = 1 - float64(la.Diff)/float64(la.AlgnLen)
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

type LAArr []LA

func (arr LAArr) Len() int {
	return len(arr)
}

func (arr LAArr) Less(i, j int) bool {
	return arr[i].AlgnLen < arr[j].AlgnLen
}

func (arr LAArr) Swap(i, j int) {
	arr[i], arr[j] = arr[j], arr[i]
}

func IsOverlapLA(la1, la2 LA) (ol bool) {
	if la1.QuyB < la2.QuyB {
		if la2.QuyB < la1.QuyE-2*Kmerlen {
			ol = true
		}
	} else {
		if la1.QuyB < la2.QuyE-2*Kmerlen {
			ol = true
		}
	}
	return ol
}

func GetMostProMax(maxArr []LA) (la LA) {
	sort.Sort(LAArr(maxArr))
	// fmt.Printf("[GetMostProMax]maxArr: %v", maxArr)
	for i := len(maxArr) - 1; i >= 0; i-- {
		var ol bool
		for j := len(maxArr) - 1; j >= 0; j-- {
			if i == j {
				continue
			}
			if IsOverlapLA(maxArr[i], maxArr[j]) && maxArr[i].Idty < maxArr[j].Idty {
				if maxArr[i].AlgnLen-maxArr[i].Diff < maxArr[j].AlgnLen-maxArr[j].Diff {
					ol = true
					break
				}
				if maxArr[j].Idty-maxArr[i].Idty < 0.05 {
					ol = true
					break
				}
			}
		}
		if !ol {
			la = maxArr[i]
			break
		}
	}

	return la

}

func GetSectionNum(la LA, uniqueRegTagArr []uint8) (num int) {
	for i := la.QuyB + TAG_ZOOM_SIZE; i < la.QuyE-TAG_ZOOM_SIZE; i += TAG_ZOOM_SIZE {
		if uniqueRegTagArr[i/TAG_ZOOM_SIZE] == 1 {
			num++
		}
	}
	return num
}

func GetMaxMappingLength(laArr []LA, gapRegArr []GapRegion, uniqueRegTagArr []uint8) (la LA) {
	// la = laArr[0]
	// var maxArr []LA
	for i := 0; i < len(laArr); i++ {
		// j := MaxIndexEID(laArr[i:], 0) + i
		var gp GapRegion
		gp.Begin, gp.End = laArr[i].QuyB, laArr[i].QuyE
		if laArr[i].Flag == true || IsInGapRegArr(gapRegArr, gp) {
			continue
		}
		// maxLa := findMaxLengthMappingLA(laArr[i : j+1])
		maxLa := MergeLA(laArr, gapRegArr, laArr[i])
		// fmt.Fprintf(os.Stderr, "[GetMaxMappingLength]i: %d, j: %d\n", i, j)
		if GetSectionNum(maxLa, uniqueRegTagArr) > 1 {
			la = maxLa
			break
		}
		// if la.RefID > 0 {
		// 	la = GetMaxLenLA(la, maxLa)
		// } else {
		// 	la = maxLa
		// }
		// fmt.Printf("[paraFindLongMappingPath] la: %v\n", la)
	}

	// la = GetMostProMax(maxArr)

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

func GetNearstLA(edgeLaArr []LA, QuyPos int, QuyDirection uint8) (la LA) {
	min := math.MaxInt32
	if QuyDirection == constructdbg.FORWARD {
		for i := 0; i < len(edgeLaArr); i++ {
			if edgeLaArr[i].QuyE > QuyPos {
				if edgeLaArr[i].QuyB-QuyPos < min {
					la = edgeLaArr[i]
					min = edgeLaArr[i].QuyB - QuyPos
				}
			}

		}
	} else {
		for i := 0; i < len(edgeLaArr); i++ {
			if edgeLaArr[i].QuyB < QuyPos {
				if QuyPos-edgeLaArr[i].QuyE < min {
					la = edgeLaArr[i]
					min = QuyPos - edgeLaArr[i].QuyE
				}
			}
		}
	}
	return la
}

func findEdgeLA(eID constructdbg.DBG_MAX_INT, laArr []LA, gapRegArr []GapRegion, QuyPos int, QuyDirection uint8) (la LA) {
	edgeLaArr := GetLARegion(laArr, eID)
	if len(edgeLaArr) == 1 {
		tmp := edgeLaArr[0]
		if QuyDirection == constructdbg.FORWARD {
			if tmp.QuyE >= QuyPos {
				la = tmp
			}
		} else {
			if tmp.QuyB <= QuyPos {
				la = tmp
			}
		}
	} else if len(edgeLaArr) > 1 {
		nearLa := GetNearstLA(edgeLaArr, QuyPos, QuyDirection)
		la = MergeLA(edgeLaArr, gapRegArr, nearLa)
	}
	return la
}

func IsOverlapQuy(la LA, QuyPos int, QuyDirection uint8) bool {
	if QuyDirection == constructdbg.FORWARD {
		if la.QuyB < QuyPos {
			return true
		} else {
			return false
		}
	} else {
		if la.QuyE > QuyPos {
			return true
		} else {
			return false
		}
	}
}

func IsBoundary(la LA) (b bool) {
	if la.QuyB < Kmerlen {
		b = true
	}
	if la.QuyE == la.QuyLen {
		b = true
	}

	return b
}

func findNextProEdge(eIDArr []constructdbg.DBG_MAX_INT, laArr []LA, gapRegArr []GapRegion, QuyPos int, QuyDirection uint8) (nLaArr []LA) {
	// direction := GetDirection(node, eID)
	// for i := 0; i < bnt.BaseTypeNum; i++ {
	// 	var nEID constructdbg.DBG_MAX_INT
	// 	if direction == constructdbg.FORWARD {
	// 		nEID = node.EdgeIDIncoming[i]
	// 	} else {
	// 		nEID = node.EdgeIDOutcoming[i]
	// 	}
	// 	if nEID == 0 {
	// 		continue
	// 	}
	for _, eID := range eIDArr {
		proLa := findEdgeLA(eID, laArr, gapRegArr, QuyPos, QuyDirection)
		// fmt.Printf("[findNextProEdge] eID: %v\tQuyPos: %v\tQuyDirection: %v\nproLa:%v\n", eID, QuyPos, QuyDirection, proLa)
		if proLa.RefID == 0 {
			continue
		}

		if len(nLaArr) == 0 {
			nLaArr = append(nLaArr, proLa)
		} else {
			// fmt.Printf("[findNextProEdge] proLa:%v\nnLaArr:%v\n", proLa, nLaArr)
			if proLa.AlgnLen == nLaArr[0].AlgnLen {
				if proLa.Idty == nLaArr[0].Idty {
					nLaArr = append(nLaArr, proLa)
				} else if proLa.Idty > nLaArr[0].Idty {
					nLaArr[0] = proLa
				}
			} else {
				var choosed bool
				if proLa.AlgnLen > nLaArr[0].AlgnLen {
					if proLa.Idty >= nLaArr[0].Idty {
						nLaArr[0] = proLa
						choosed = true
					}
				} else {
					if proLa.Idty <= nLaArr[0].Idty {
						choosed = true
					}
				}
				if choosed == false {
					if IsOverlapQuy(proLa, QuyPos, QuyDirection) {
						if IsOverlapQuy(nLaArr[0], QuyPos, QuyDirection) == false {
							nLaArr[0] = proLa
							choosed = true
						}
					} else {
						if IsOverlapQuy(nLaArr[0], QuyPos, QuyDirection) {
							choosed = true
						}
					}
				}
				if !choosed {
					proBoundary := IsBoundary(proLa)
					nlaBoundary := IsBoundary(nLaArr[0])
					if proBoundary && !nlaBoundary {
						nLaArr[0] = proLa
						choosed = true
					} else if !proBoundary && nlaBoundary {
						choosed = true
					} else if proBoundary && nlaBoundary {
						nLaArr = append(nLaArr, proLa)
						choosed = true
					}
				}
				if !choosed {

					proLenP := float64(proLa.RefE-proLa.RefB) / float64(proLa.RefLen-2*len(constructdbg.AdpaterSeq))
					nlaLenP := float64(nLaArr[0].RefE-nLaArr[0].RefB) / float64(nLaArr[0].RefLen-2*len(constructdbg.AdpaterSeq))
					if math.Abs(proLenP-nlaLenP) > 0.3 {
						if proLenP > 0.9 {
							nLaArr[0] = proLa
							choosed = true
						} else if nlaLenP > 0.9 {
							choosed = true
						}
					}
				}

				if !choosed {
					nLaArr = append(nLaArr, proLa)
				}
			}
		}
	}
	/*} else {
			if node.EdgeIDOutcoming[i] > 0 {
				proLa := findEdgeLA(node.EdgeIDOutcoming[i], laArr, gapRegArr)
				if proLa.RefID == 0 {
					continue
				}
				if proLa.QuyE > la.QuyE {
					if la.QuyE-proLa.QuyB < Kmerlen/QUY_OVL_FACTOR {
						continue
					}
				} else {
					if proLa.QuyE-la.QuyB < Kmerlen/QUY_OVL_FACTOR {
						continue
					}
				}
				// nLaArr = append(nLaArr, proLa)
				if len(nLaArr) == 0 {
					nLaArr = append(nLaArr, proLa)
				} else {
					fmt.Printf("[findNextProEdge] proLa:%v\nnLaArr:%v\n", proLa, nLaArr)
					if proLa.AlgnLen == nLaArr[0].AlgnLen {
						if proLa.Idty == nLaArr[0].Idty {
							nLaArr = append(nLaArr, proLa)
						} else if proLa.Idty > nLaArr[0].Idty {
							nLaArr[0] = proLa
						}
					} else {
						if proLa.AlgnLen > nLaArr[0].AlgnLen {
							if proLa.Idty >= nLaArr[0].Idty {
								nLaArr[0] = proLa
							} else {
								if proLa.AlgnLen > 2*nLaArr[0].AlgnLen {
									nLa = proLa
								} else {
									var la LA
									nLa = la
									break
								}
							}
						} else {
							if proLa.Idty > nLaArr[0].Idty && proLa.AlgnLen*2 > nLaArr[0].AlgnLen {
								var la LA
								nLa = la
								break
							}
						}
					}
				}
			}
		}
	} */
	return nLaArr
}

func findNextProPath(subp [2][]constructdbg.DBG_MAX_INT, laArr []LA, gapRegArr []GapRegion, QuyPos int, QuyDirection uint8) int {
	idx := -1
	var QuyStart, QuyEnd, mapLen, errorNum [2]int
	for i := 0; i < 2; i++ {
		for _, id := range subp[i] {
			proLa := findEdgeLA(id, laArr, gapRegArr, QuyPos, QuyDirection)
			if proLa.RefID == 0 {
				continue
			}
			if QuyDirection == constructdbg.FORWARD {
				if QuyStart[i] == 0 {
					QuyStart[i] = proLa.QuyB
				}
				QuyEnd[i] = proLa.QuyE
				QuyPos = proLa.QuyE
			} else {
				if QuyStart[i] == 0 {
					QuyStart[i] = proLa.QuyE
				}
				QuyEnd[i] = proLa.QuyB
				QuyPos = proLa.QuyB
			}
			mapLen[i] += proLa.AlgnLen
			errorNum[i] += proLa.Diff
		}
	}

	// var bl int
	// if QuyEnd > QuyStart {
	//   bl = QuyEnd -  QuyStart
	// } else {
	//   bl = QuyStart - QuyEnd
	// }

	if mapLen[0] < 2*Kmerlen-1 || mapLen[1] < 2*Kmerlen-1 {
		return idx
	}

	if mapLen[0] > mapLen[1] {
		if errorNum[0] < errorNum[1] {
			idx = 0
		}
	} else if mapLen[0] < mapLen[1] {
		if errorNum[0] > errorNum[1] {
			idx = 1
		}
	} else {
		if errorNum[0] > errorNum[1] {
			idx = 1
		} else if errorNum[0] < errorNum[1] {
			idx = 0
		}
	}

	return idx
}

func SumLen(eIDArr []constructdbg.DBG_MAX_INT, edgesArr []constructdbg.DBGEdge) (sl int) {
	sl = Kmerlen - 1
	for _, eID := range eIDArr {
		sl += len(edgesArr[eID].Utg.Ks) - (Kmerlen - 1)
	}

	return sl
}

func IsSubPath(srcArr, subArr []constructdbg.DBG_MAX_INT) (is bool) {
	for i, eID := range srcArr {
		if i+len(subArr) > len(srcArr) {
			break
		}
		if eID == subArr[0] && reflect.DeepEqual(srcArr[i:i+len(subArr)], subArr) {
			is = true
			break
		}
	}

	return is
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

func ComputeRegLen(la LA, gapRegArr []GapRegion, maxGapLenAllow int) (regLen int) {
	regLen = la.QuyLen
	for i, v := range gapRegArr {
		if la.QuyB < v.Begin {
			if i > 0 {
				for j := i - 1; j >= 0; j-- {
					if gapRegArr[j].End-gapRegArr[j].Begin > maxGapLenAllow {
						regLen -= gapRegArr[j].End
						break
					}
				}
			} else {
				regLen = v.Begin
				break
			}
			// fmt.Printf("[ComputeRegLen] regLen: %v\n", regLen)
			for j := i; j < len(gapRegArr); j++ {
				if gapRegArr[j].End-gapRegArr[j].Begin > maxGapLenAllow {
					regLen -= (la.QuyLen - gapRegArr[j].Begin)
					break
				}
			}
			break
		}
		if i == len(gapRegArr)-1 && la.QuyB > v.End {
			regLen -= v.End
			break
		}
	}

	return regLen
}

func AddGapRegion(gapRegArr []GapRegion, mapReg GapRegion) (newArr []GapRegion) {
	if len(gapRegArr) == 0 {
		mapReg.Begin = 0
		mapReg.End = 1000000
		newArr = append(newArr, mapReg)
	} else {
		i := 0
		for ; i < len(gapRegArr); i++ {
			if mapReg.Begin <= gapRegArr[i].End {
				var mergeReg GapRegion
				if i == 0 {
					mergeReg.Begin = 0
				} else {
					mergeReg.Begin = gapRegArr[i-1].Begin
				}
				j := i
				for ; j < len(gapRegArr); j++ {
					if mapReg.End <= gapRegArr[j].End {
						break
					}
				}
				if j == len(gapRegArr) {
					mergeReg.End = 1000000
				} else {
					mergeReg.End = gapRegArr[j].End
				}
				if i > 0 {
					newArr = append(newArr, gapRegArr[:i-1]...)
				}
				newArr = append(newArr, mergeReg)
				if j < len(gapRegArr)-1 {
					newArr = append(newArr, gapRegArr[j+1:]...)
				}
				break
			}
		}
		// fmt.Printf("[AddGapRegion] i: %v\n", i)
		if i == len(gapRegArr) {
			gapRegArr[i-1].End = 1000000
			newArr = gapRegArr
		}
	}
	return newArr
}

func MaxInt(n, m int) (max int) {
	max = n
	if m > max {
		max = m
	}

	return max
}

func AddScore(laArr []LA) (ql, ml int, sc float64) {
	for _, e := range laArr {
		if e.RefE > Kmerlen-1 {
			ml += e.RefE - MaxInt(Kmerlen-1, e.RefB)
		}
		ql := e.QuyE - e.QuyB
		sc += float64(ql) * e.Idty
	}

	return ql, ml, sc
}

func GetNextShortPath(consistenceA [][]constructdbg.DBG_MAX_INT, i, j int) (ni, nj int) {
	idx := IndexEID(consistenceA[i], constructdbg.DBG_MAX_INT(i))
	if idx < j {
		if j+1 < len(consistenceA[i]) {
			ni, nj = i, j+1
		}
	} else {
		if j-1 >= 0 {
			ni, nj = i, j-1
		}
	}
	return ni, nj
}

func IsInComing(eIDcoming [4]constructdbg.DBG_MAX_INT, eID constructdbg.DBG_MAX_INT) bool {
	for _, id := range eIDcoming {
		if id == eID {
			return true
		}
	}
	return false
}

func GetNextDirection(eID constructdbg.DBG_MAX_INT, edge constructdbg.DBGEdge, nodesArr []constructdbg.DBGNode) (direction uint8) {
	if edge.StartNID > 0 {
		n1 := nodesArr[edge.StartNID]
		if GetDirection(n1, edge.ID) == constructdbg.FORWARD {
			if IsInComing(n1.EdgeIDIncoming, eID) {
				return constructdbg.BACKWARD
			}
		} else {
			if IsInComing(n1.EdgeIDOutcoming, eID) {
				return constructdbg.BACKWARD
			}
		}
	}

	if edge.EndNID > 0 {
		n2 := nodesArr[edge.EndNID]
		if GetDirection(n2, edge.ID) == constructdbg.FORWARD {
			if IsInComing(n2.EdgeIDIncoming, eID) {
				return constructdbg.FORWARD
			}
		} else {
			if IsInComing(n2.EdgeIDOutcoming, eID) {
				return constructdbg.BACKWARD
			}
		}
	}

	// not reachable here
	log.Fatalf("[GetNextDirection] direction not set\n")
	return direction
}

func GetNextConsisPath(edgesArr []constructdbg.DBGEdge, nodesArr []constructdbg.DBGNode, left []constructdbg.DBG_MAX_INT, consistenceA [][]constructdbg.DBG_MAX_INT) (n, m int) {
	for i := len(left) - 1; i > 0; i-- {
		if len(consistenceA[left[i]]) > 0 {
			forA := consistenceA[left[i]]
			direction := GetNextDirection(left[i-1], edgesArr[left[i]], nodesArr)
			if direction == constructdbg.FORWARD {
				forA = GetReverseDBG_MAX_INTArr(forA)
			}
			_, ok := CheckConsistenceArr(left, forA, left[i])
			if ok == false {
				log.Fatalf("[GetNextConsisPath] left and forA not consistence\n")
			}
			idx := IndexEID(forA, left[i])
			if len(left)-i < len(forA)-idx {
				n = int(left[i])
				m = idx + (len(left) - i)
				if direction == constructdbg.FORWARD {
					m = len(forA) - 1 - m
				}
			}
			break
		}
	}

	return n, m
}

func GetNextConsistenceAPosition(consistenceA [][]constructdbg.DBG_MAX_INT, eID constructdbg.DBG_MAX_INT, direction uint8) (n, m int) {
	idx := IndexEID(consistenceA[eID], eID)
	if direction == constructdbg.BACKWARD {
		if idx > 0 {
			n = int(eID)
			m = int(consistenceA[n][idx-1])
		}
	} else {
		if idx < len(consistenceA[eID])-1 {
			n = int(eID)
			m = int(consistenceA[n][idx+1])
		}
	}

	return n, m
}

func IsInLAArr(laArr []LA, eID constructdbg.DBG_MAX_INT) bool {
	for _, la := range laArr {
		if la.RefID == eID {
			return true
		}
	}

	return false
}

func GetRemainArr(csArr []constructdbg.DBG_MAX_INT, eID, lasteID constructdbg.DBG_MAX_INT) (remainArr []constructdbg.DBG_MAX_INT) {
	idx := IndexEID(csArr, eID)
	if 0 < idx && idx < len(csArr)-1 {
		if lasteID == csArr[idx-1] {
			if len(csArr)-(idx+1) > 0 {
				remainArr = make([]constructdbg.DBG_MAX_INT, len(csArr)-(idx+1))
				copy(remainArr, csArr[idx+1:])
			}
		} else if lasteID == csArr[idx+1] {
			remainArr = GetReverseDBG_MAX_INTArr(csArr[0:idx])
		} else {
			log.Fatalf("[GetRemainArr] eID:%d not in csArr:%v\n", lasteID, csArr)
		}
	} else {
		if idx == 0 {
			if lasteID != csArr[idx+1] {
				remainArr = make([]constructdbg.DBG_MAX_INT, len(csArr)-(1))
				copy(remainArr, csArr[1:])
			}
		} else { // idx == len(csArr)-1
			if lasteID != csArr[len(csArr)-2] {
				remainArr = GetReverseDBG_MAX_INTArr(csArr[0 : len(csArr)-1])
			}
		}
	}

	return remainArr
}

/*func GetNeedSortEID(edgesArr []constructdbg.DBGEdge, rc chan constructdbg.DBG_MAX_INT, numCPU int) {
	for _, e := range edgesArr {
		if e.ID > 0 && len(e.PathMat) > 0 {
			rc <- e.ID
		}
	}
	// fmt.Printf("[GetNeedSortEID] rc: %v\n", rc)
	for i := 0; i < numCPU; i++ {
		rc <- 0
	}
}*/

type PathM []constructdbg.Path

func (pm PathM) Len() int {
	return len(pm)
}

func (pm PathM) Less(i, j int) bool {
	return pm[i].Freq < pm[j].Freq
}
func (pm PathM) Swap(i, j int) {
	pm[i], pm[j] = pm[j], pm[i]
}

/*func paraSortPathMat(edgesArr []constructdbg.DBGEdge, rc chan constructdbg.DBG_MAX_INT) {

	for {
		eID := <-rc
		if eID == 0 {
			break
		}
		// fmt.Printf("[paraSortPathMat] pathMat: %v\n", edgesArr[eID].PathMat)
		sort.Sort(sort.Reverse(PathM(edgesArr[eID].PathMat)))
		fmt.Printf("[paraSortPathMat] pathMat: %v\n", edgesArr[eID].PathMat)
	}
	// fmt.Printf("[paraSortPathMat] edgesArr: %v\n", edgesArr)
	fmt.Printf("[paraSortPathMat] edgesArr \n")
	// for i, _ := range edgesArr {
	// 	fmt.Printf("[Fpath] pathMat: %v\n", edgesArr[i].PathMat)
	// }
}*/

func SortPathMat(edgesArr []constructdbg.DBGEdge) {

	for i, _ := range edgesArr {
		if len(edgesArr[i].PathMat) > 1 {
			sort.Sort(sort.Reverse(PathM(edgesArr[i].PathMat)))
		}
	}
	// for i, _ := range edgesArr {
	// 	fmt.Printf("[Fpath] pathMat: %v\n", edgesArr[i].PathMat)
	// }
}

func GetNextPathLA(eIDArr []constructdbg.DBG_MAX_INT, laArr []LA) (remainLaArr []LA) {
	for _, eID := range eIDArr {
		edgeLaArr := GetLARegion(laArr, eID)
		if len(edgeLaArr) > 0 {
			la := findMaxLengthMappingLA(edgeLaArr)
			remainLaArr = append(remainLaArr, la)
		} else {
			break
		}
	}

	return remainLaArr
}

type PathFreqArr struct {
	id   []constructdbg.DBG_MAX_INT
	freq []int
}

/*func findConsistenceExtend(pathMat []constructdbg.Path, edgesArr []constructdbg.DBGEdge, nodesArr []constructdbg.DBGNode) (consis []PathFreqArr) {

	for _, path := range pathMat {
		if len(consis) == 0 {
			var pfa PathFreqArr
			pfa.id = append(pfa.id, path.IDArr...)
			pfa.freq = make([]int, len(pfa.id))
			for i := 0; i < len(pfa.freq); i++ {
				pfa.freq[i] = path.Freq
			}
			consis = append(consis, pfa)
		} else {
			for j := 0; j < len(consis); j++ {
				n := IntersectionNum(consis[j].id, path.IDArr)
				if n == 0 {

				} else {

				}

			}

		}
	}
} */

func GetInterSection(arr1, arr2 []constructdbg.DBG_MAX_INT) (interSecArr []constructdbg.DBG_MAX_INT) {
	for _, e1 := range arr1 {
		if IsInDBG_MAX_INTArr(arr2, e1) {
			interSecArr = append(interSecArr, e1)
		}
	}
	return interSecArr
}

func CheckAndCutLocalCycle(eIDArr []constructdbg.DBG_MAX_INT) (cutArr []constructdbg.DBG_MAX_INT, cuted bool) {
	var cycleArr []constructdbg.DBG_MAX_INT
	cutArr = eIDArr
	i := len(eIDArr) - 2
	for i >= 0 {
		if eIDArr[i] == eIDArr[len(eIDArr)-1] {
			cycleArr = eIDArr[i+1:]
			break
		}
		i--
	}
	i++
	// fmt.Printf("[CheckAndCutLocalCycle] cycleArr: %v\n", cycleArr)
	for i -= len(cycleArr); i >= 0; i -= len(cycleArr) {
		if reflect.DeepEqual(eIDArr[i:i+len(cycleArr)], cycleArr) {
			cutArr = eIDArr[:i]
			cuted = true
		} else {
			break
		}
	}

	return cutArr, cuted
}

func IsOverlapGapReg(g, gp GapRegion) (ok bool) {
	if g.Begin < gp.Begin {
		if gp.Begin < g.End {
			ok = true
		}
	} else {
		if g.Begin < gp.End {
			ok = true
		}
	}

	return ok
}

func IsOverlapGapRegArr(gapRegArr []GapRegion, QuyPos int, QuyDirection uint8, extLen int) (ok bool) {
	var gp GapRegion
	if QuyDirection == constructdbg.FORWARD {
		gp.Begin = QuyPos
		gp.End = QuyPos + extLen
	} else {
		gp.Begin = QuyPos - extLen
		gp.End = QuyPos
	}

	for _, g := range gapRegArr {
		if g.Begin > gp.End {
			break
		}
		if IsOverlapGapReg(g, gp) {
			ok = true
			break
		}
	}

	return ok
}

func ExtendGap(la LA, gapRegArr []GapRegion) LA {
	extLa := la
	if la.RefB > len(constructdbg.AdpaterSeq)+Kmerlen/REF_OVL_FACTOR {
		var gp GapRegion
		if la.Minus {
			gp.Begin = la.QuyE
			gp.End = la.QuyE + la.RefB - len(constructdbg.AdpaterSeq)
		} else {
			gp.Begin = la.QuyB - (la.RefB - len(constructdbg.AdpaterSeq))
			gp.End = la.QuyB
		}

		if IsCrossGapRegArr(gapRegArr, gp) {
			extLa.RefB = len(constructdbg.AdpaterSeq)
		}
	}
	if la.RefE < la.RefLen-len(constructdbg.AdpaterSeq)-Kmerlen/REF_OVL_FACTOR {
		var gp GapRegion
		if la.Minus {
			gp.Begin = la.QuyB - (la.RefLen - la.RefE - len(constructdbg.AdpaterSeq))
			gp.End = la.QuyB
		} else {
			gp.Begin = la.QuyE
			gp.End = la.QuyE + la.RefLen - la.RefE - len(constructdbg.AdpaterSeq)
		}
		if IsCrossGapRegArr(gapRegArr, gp) {
			extLa.RefE = la.RefLen - len(constructdbg.AdpaterSeq)
		}
	}

	return extLa
}

func IsHasLongReg(gapRegArr []GapRegion, minLongRegLen, readLen int) bool {
	for i := 0; i < len(gapRegArr); i++ {
		var begin int
		if i == 0 {
			if gapRegArr[i].Begin-minLongRegLen > 0 {
				return true
			}
		} else {
			begin = gapRegArr[i-1].End
		}

		for ; i < len(gapRegArr); i++ {
			if gapRegArr[i].End-gapRegArr[i].Begin > MAX_ALIGN_GAP_ALLOW {
				break
			}
		}
		if i == len(gapRegArr) {
			if readLen-begin-minLongRegLen > 0 {
				return true
			}
		} else {
			if gapRegArr[i].Begin-begin-minLongRegLen > 0 {
				return true
			}
		}
	}

	return false
}

func SetLAFlag(laArr []LA, eID constructdbg.DBG_MAX_INT) {
	for i, _ := range laArr {
		if laArr[i].RefID == eID {
			laArr[i].Flag = true
		}
	}
}

type LAReg struct {
	Start, End int
	LaIndex    constructdbg.DBG_MAX_INT
}

func MaskRepeatMapping(laArr []LA) []uint8 {
	for i := 1; i < len(laArr); i++ {
		v := laArr[i]
		for j, e := range laArr[:i] {
			if e.Flag {
				continue
			}
			if e.QuyB <= v.QuyB && v.QuyE <= e.QuyE {
				if v.Idty < e.Idty || (v.RefID == e.RefID && e.RefB <= v.RefB && v.RefE <= e.RefE) {
					laArr[i].Flag = true
				}
				break
			} else if v.QuyB <= e.QuyB && e.QuyE <= v.QuyE {
				if e.Idty < v.Idty || (v.RefID == e.RefID && v.RefB <= e.RefB && e.RefE <= v.RefE) {
					laArr[j].Flag = true
				}
				break
			}
		}
	}

	// get unique region
	uniqueRegTagArr := make([]uint8, (laArr[0].QuyLen+TAG_ZOOM_SIZE-1)/TAG_ZOOM_SIZE)
	for _, v := range laArr {
		if v.Flag {
			continue
		}
		for j := v.QuyB + TAG_ZOOM_SIZE; j < v.QuyE-TAG_ZOOM_SIZE; j += TAG_ZOOM_SIZE {
			uniqueRegTagArr[j/TAG_ZOOM_SIZE]++
		}
	}

	return uniqueRegTagArr
}

func paraFindLongMappingPath(lac chan []LA, wc chan []constructdbg.DBG_MAX_INT, edgesArr []constructdbg.DBGEdge, nodesArr []constructdbg.DBGNode, RefIDMapArr []constructdbg.DBG_MAX_INT, consistenceA [][2][]constructdbg.DBG_MAX_INT) {
	for {
		laArr := <-lac
		if len(laArr) == 0 {
			var guardPath []constructdbg.DBG_MAX_INT
			wc <- guardPath
			break
		}

		sort.Sort(sort.Reverse(LAArr(laArr)))
		gapRegArr := GetLAArrGapRegion(laArr)
		// found repeat region and masked pseudo mappming record, get Unique region
		uniqueRegTagArr := MaskRepeatMapping(laArr)
		fmt.Printf("[paraFindLongMappingPath] uniqueRegTagArr: %v\n", uniqueRegTagArr)
		lowCount := 0
		for _, t := range uniqueRegTagArr {
			if t <= 2 {
				lowCount++
			}
		}
		fmt.Printf("[paraFindLongMappingPath] len(uniqueRegTagArr): %v\tlowCount: %v\n", len(uniqueRegTagArr), lowCount)
		// MaskRepeatMapping(laArr)
		for i, te := range laArr {
			fmt.Printf("[paraFindLongMappingPath] laArr[%d]: %v\n", i, te)
		}
		for _, gr := range gapRegArr {
			if gr.Begin != 0 {
				// fmt.Fprintf(os.Stderr, "\t[%d,\t%d],\n", gr.Begin, gr.End)
				fmt.Fprintf(os.Stderr, "%d\t%d\n", gr.Begin, gr.End)
			}
		}
		// found all regions map to the DBG, when found one region and added to the mapRegArr used tagged
		for {
			fmt.Printf("[paraFindLongMappingPath] gapRegArr: %v\n", gapRegArr)
			if IsHasLongReg(gapRegArr, MIN_LONG_REG_LEN, laArr[0].QuyLen) == false {
				break
			}
			// var path []constructdbg.DBG_MAX_INT
			maxLa := GetMaxMappingLength(laArr, gapRegArr, uniqueRegTagArr)
			if maxLa.RefID == 0 {
				// for i, te := range laArr {
				// 	fmt.Printf("[paraFindLongMappingPath] laArr[%d]: %v\n", i, te)
				// }
				// log.Fatalf("[paraFindLongMappingPath] not found proper max la\n")
				break
			}
			SetLAFlag(laArr, maxLa.RefID)
			tmp := ComputeRegLen(maxLa, gapRegArr, MAX_ALIGN_GAP_ALLOW)
			fmt.Printf("[paraFindLongMappingPath] maxLa: %v\n", maxLa)
			fmt.Printf("[paraFindLongMappingPath] ComputeRegLen: %v\n", tmp)
			if ComputeRegLen(maxLa, gapRegArr, MAX_ALIGN_GAP_ALLOW) < MIN_LONG_REG_LEN {
				continue
			}
			maxLa = ExtendGap(maxLa, gapRegArr)
			var mapReg GapRegion
			mapReg.Begin, mapReg.End = maxLa.QuyB, maxLa.QuyE

			edge := edgesArr[maxLa.RefID]
			// fmt.Printf("[paraFindLongMappingPath] RefLen: %d, edges Length: %d\n", maxLa.RefLen, len(edge.Utg.Ks))
			var extArr [2][]constructdbg.DBG_MAX_INT
			// if maxLa.RefB < Min(len(edge.Utg.Ks)/FLANK_ALLOW_FACTOR, MAX_FLANK_ALLOW_LEN)+len(constructdbg.AdpaterSeq) {
			flankLen := Kmerlen/REF_OVL_FACTOR + len(constructdbg.AdpaterSeq)
			for i := 0; i < 2; i++ {
				if (i == 0 && maxLa.RefB > flankLen) || (i == 1 && maxLa.RefE < len(edge.Utg.Ks)-flankLen) {
					continue
				}
				la := maxLa
				eID := la.RefID
				var nID constructdbg.DBG_MAX_INT
				var QuyPos int
				var QuyDirection uint8
				if i == 0 {
					nID = edgesArr[eID].StartNID
					if la.Minus {
						QuyPos = la.QuyE
						QuyDirection = constructdbg.FORWARD
					} else {
						QuyPos = la.QuyB
						QuyDirection = constructdbg.BACKWARD
					}
				} else {
					nID = edgesArr[eID].EndNID
					if la.Minus {
						QuyPos = la.QuyB
						QuyDirection = constructdbg.BACKWARD
					} else {
						QuyPos = la.QuyE
						QuyDirection = constructdbg.FORWARD
					}
				}
				fmt.Printf("[paraFindLongMappingPath] extend %d cycle\n", i)
				// found consistence path of short path
				// n, m := -1, -1 // the position of consistenceA
				// var direction uint8
				// if i == 0 {
				// 	direction = constructdbg.BACKWARD
				// } else {
				// 	direction = constructdbg.FORWARD
				// }
				// if len(consistenceA[eID]) > 0 {
				// 	n, m = GetNextConsistenceAPosition(consistenceA, eID, direction)
				// }
				continueNoMapEdgeNum := 0
				for nID > 0 {
					//neID, b, e, rl
					// fmt.Printf("[paraFindLongMappingPath]extArr: %v\n", extArr[i])
					// fmt.Printf("[paraFindLongMappingPath]nID: %v\n", nID)
					eIDArr := GetNextEArr(eID, nodesArr[nID])
					if len(eIDArr) == 0 {
						break
					}
					nlaArr := findNextProEdge(eIDArr, laArr, gapRegArr, QuyPos, QuyDirection)
					fmt.Printf("[paraFindLongMappingPath]nlaArr: %v\n", nlaArr)
					if len(nlaArr) > 0 {
						if nlaArr[0].QuyB < mapReg.Begin {
							mapReg.Begin = nlaArr[0].QuyB
						}
						if nlaArr[0].QuyE > mapReg.End {
							mapReg.End = nlaArr[0].QuyE
						}
						continueNoMapEdgeNum = 0
					} else {
						continueNoMapEdgeNum++
					}
					flag := true // flag if need short read path assist
					if len(nlaArr) == 0 {
						if len(eIDArr) == 1 {
							extArr[i] = append(extArr[i], eIDArr[0])
							flag = false
						} else if len(eIDArr) == 2 {
							if IsOverlapGapRegArr(gapRegArr, QuyPos, QuyDirection, len(edgesArr[eIDArr[0]].Utg.Ks)-Kmerlen) {
								e1 := edgesArr[eIDArr[0]]
								e2 := edgesArr[eIDArr[1]]
								if (e1.StartNID == e2.StartNID && e1.EndNID == e2.EndNID) || (e1.StartNID == e2.EndNID && e1.EndNID == e1.StartNID) {
									extArr[i] = append(extArr[i], 0)
									flag = false
									eID = e1.ID
								} else {
									var arr [2][]constructdbg.DBG_MAX_INT
									for j := 0; j < 2; j++ {
										edge := edgesArr[eIDArr[j]]
										var nodeID constructdbg.DBG_MAX_INT
										if edge.StartNID == nID {
											nodeID = edge.EndNID
										} else {
											nodeID = edge.StartNID
										}
										if nodeID > 0 {
											arr[j] = GetNextEArr(edge.ID, nodesArr[nodeID])
										}
									}
									nextLAArr := findNextProEdge(append(arr[0], arr[1]...), laArr, gapRegArr, QuyPos, QuyDirection)
									if len(nextLAArr) == 1 {
										id := nextLAArr[0].RefID
										if IsInDBG_MAX_INTArr(arr[0], id) {
											extArr[i] = append(extArr[i], e1.ID)
											extArr[i] = append(extArr[i], id)
											if nID == e1.StartNID {
												nID = e1.EndNID
											} else {
												nID = e1.StartNID
											}
										} else {
											extArr[i] = append(extArr[i], e2.ID)
											extArr[i] = append(extArr[i], id)
											if nID == e2.StartNID {
												nID = e2.EndNID
											} else {
												nID = e2.StartNID
											}
										}
										flag = false
									}

								}
							}
						}
						if flag && QuyPos > 500 && QuyPos < maxLa.QuyLen-Kmerlen-30 {
							fmt.Printf("[paraFindLongMappingPath] QuyPos: %d\nextArr: %v\neIDArr: %v\n", QuyPos, extArr[i], eIDArr)
							// for i, te := range laArr {
							// 	fmt.Printf("[paraFindLongMappingPath] laArr[%d]: %v\n", i, te)
							// }
						}
						// break
					} else if len(nlaArr) == 1 {
						extArr[i] = append(extArr[i], nlaArr[0].RefID)
						if QuyDirection == constructdbg.FORWARD {
							QuyPos = nlaArr[0].QuyE
						} else {
							QuyPos = nlaArr[0].QuyB
						}
						flag = false
					} else if len(nlaArr) == 2 {
						if SumLen(extArr[i], edgesArr) > MAX_INSERT_SIZE {
							var neID constructdbg.DBG_MAX_INT
							var nosub bool
							num := 0
							idx := -1
							for j := 0; j < 2; j++ {
								id := nlaArr[j].RefID
								var subp []constructdbg.DBG_MAX_INT
								if edgesArr[id].StartNID == nID {
									subp = consistenceA[id][0]
								} else {
									subp = consistenceA[id][1]
								}
								// fmt.Printf("[paraFindLongMappingPath] before subp: %v\n", subp)
								if len(subp) == 0 {
									nosub = true
									break
								}
								if IsSubPath(extArr[i], GetReverseDBG_MAX_INTArr(subp)) {
									neID = id
									idx = j
									num++
								}
							}
							if nosub == false && num == 1 {
								extArr[i] = append(extArr[i], neID)
								if QuyDirection == constructdbg.FORWARD {
									QuyPos = nlaArr[idx].QuyE
								} else {
									QuyPos = nlaArr[idx].QuyB
								}
								flag = false
							}
						}

						if flag && len(edgesArr[nlaArr[0].RefID].Utg.Ks) == len(edgesArr[nlaArr[1].RefID].Utg.Ks) {
							var subp [2][]constructdbg.DBG_MAX_INT
							var noSub bool
							for j := 0; j < 2; j++ {
								id := nlaArr[j].RefID
								if edgesArr[id].StartNID == nID {
									subp[j] = consistenceA[id][1]
								} else {
									subp[j] = consistenceA[id][0]
								}
								// fmt.Printf("[paraFindLongMappingPath] after subp: %v\n", subp[j])
								if len(subp[j]) == 0 {
									noSub = true
									break
								}
							}
							if noSub == false {
								idx := findNextProPath(subp, laArr, gapRegArr, QuyPos, QuyDirection)
								if idx >= 0 {
									extArr[i] = append(extArr[i], nlaArr[idx].RefID)
									if QuyDirection == constructdbg.FORWARD {
										QuyPos = nlaArr[idx].QuyE
									} else {
										QuyPos = nlaArr[idx].QuyB
									}
									flag = false
								}

							}
						}

						if flag {
							e1 := edgesArr[nlaArr[0].RefID]
							e2 := edgesArr[nlaArr[1].RefID]
							if (e1.StartNID == e2.StartNID && e1.EndNID == e2.EndNID) || (e1.StartNID == e2.EndNID && e1.EndNID == e2.StartNID) {
								extArr[i] = append(extArr[i], 0)
								if QuyDirection == constructdbg.FORWARD {
									QuyPos = nlaArr[0].QuyE
								} else {
									QuyPos = nlaArr[0].QuyB
								}
								flag = false
								eID = e1.ID
							}
						}
					}

					if flag || continueNoMapEdgeNum > 1 {
						// if len(nlaArr) > 0 && (nlaArr[0].QuyB > 500 || nlaArr[0].QuyE < nlaArr[0].QuyLen-Kmerlen-30) {
						// 	fmt.Printf("[paraFindLongMappingPath] nlaArr: %v\n", nlaArr)
						// 	for i, te := range laArr {
						// 		fmt.Printf("[paraFindLongMappingPath] laArr[%d]: %v\n", i, te)
						// 	}
						// }
						break
					}
					// fmt.Printf("[paraFindLongMappingPath] len(extArr[i]):%d, last eID: %d\n", len(extArr[i]), extArr[i][len(extArr[i])-1])
					// fmt.Printf("[paraFindLongMappingPath] len(extArr[i]):%d, last eID: %v\n", len(extArr[i]), extArr[i])
					if len(extArr[i]) > 100 {
						log.Fatalf("[paraFindLongMappingPath] len(extArr[i]) > 100\nextArr[%d]: %v\n", i, extArr[i])
					}

					/*if flag && len(nlaArr) == 2 {

					}

					if flag && len(edgesArr[extArr[i][len(extArr[i])-1]].Utg.Ks) < MAX_INSERT_SIZE {
						remainLen := MAX_INSERT_SIZE
						for j := len(extArr[i]) - 1; j >= 0 && remainLen > 0; j-- {
							leID := extArr[i][j]
							if len(consistenceA[leID][0]) > 0 {
								inter := GetInterSection(eIDArr, consistenceA[leID][0])
								if len(inter) == 1 {
									extArr[i] = append(extArr[i], inter[0])
									flag = false
									break
								}
							}

							if len(consistenceA[leID][1]) > 0 {
								inter := GetInterSection(eIDArr, consistenceA[leID][1])
								if len(inter) == 1 {
									extArr[i] = append(extArr[i], inter[0])
									flag = false
									break
								}
							}

							remainLen -= (len(edgesArr[leID].Utg.Ks) - Kmerlen + 1)
						}

					}

					if flag {
						break
					}
					if n > 0 && m >= 0 {
						if IsInLAArr(nlaArr, consistenceA[n][m]) {
							extArr[i] = append(extArr[i], consistenceA[n][m])
						} else {
							log.Fatalf("[paraFindLongMappingPath] edgeID:%d not in nlaArr:%v\n", consistenceA[n][m], nlaArr)
						}
					} else {
						if len(nlaArr) > 2 {
							log.Fatalf("[paraFindLongMappingPath] nlaArr:%v\n", nlaArr)
						} else {
							// var eIDArr [2]constructdbg.DBG_MAX_INT
							// e[0], e[1] := edgesArr[nlaArr[0].RefID], edgesArr[nlaArr[1].RefID]
							// fmt.Printf("[paraFindLongMappingPath] nlaArr:%v\nextArr[%d]:%v\te1 pathMat: %v\te2 pathMat:%v\n", nlaArr, i, extArr[i], e1.PathMat, e2.PathMat)
							var laRArr [2][]LA // region laArr
							var remainArr [2][]constructdbg.DBG_MAX_INT
							var maplen, quylen [2]int
							var score [2]float64
							for j := 0; j < 2; j++ {
								if nlaArr[j].QuyB > Kmerlen && nlaArr[j].QuyE < nlaArr[j].QuyLen-Kmerlen/QUY_OVL_FACTOR {
									consis := findConsistenceExtend(edgesArr[nlaArr[j].RefID].PathMat)
									if len(consis) != 2 {
										continue
									}
									if IsSubPath(extArr[i], consis[0]) {
										remainArr[j] = consis[1]
									} else if IsSubPath(extArr[i], consise1[1]) {
										remainArr[j] = consis[0]
									}
									if len(remainArr[j]) > 0 {
										laRArr[j] = GetNextPathLA(remainArr[j], laArr)
									}
									maplen[j] = nlaArr[j].RefE - nlaArr[j].RefB
									quylen[j] = nlaArr[j].QuyE - nlaArr[j].QuyB
									score[j] = float64(quylen[j]) * nlaArr[j].Idty
									if len(laRArr[j]) > 0 {
										ql, ml, sc := AddScore(laRArr[j])
										quylen[j] += ql
										maplen[j] += ml
										score[j] += sc
									}
								}
							}

							// check
							if maplen[0] < maplen[1] {
								if score[0]/float64(quylen[0]) < score[1]/float64(quylen[1]) {
									extArr[i] = append(extArr[i], nlaArr[1].RefID)
									extArr[i] = append(extArr[i], remainArr[1]...)
									la = laRArr[1][len(laRArr[1])-1]
								} else {
									fmt.Printf("[paraFindLongMappingPath] not found most probability path laArr1:%v\nlaArr2:%v\n", laRArr[0], laRArr[1])
									break
								}
							} else {
								if score[0]/float64(quylen[0]) > score[1]/float64(quylen[1]) {
									extArr[i] = append(extArr[i], nlaArr[0].RefID)
									extArr[i] = append(extArr[i], remainArr[0]...)
									la = laRArr[0][len(laRArr[0])-1]
								} else {
									fmt.Printf("[paraFindLongMappingPath] not found most probability path laArr1:%v\nlaArr2:%v\n", laRArr[0], laRArr[1])
									break
								}
							}
						}
					}*/

					if extArr[i][len(extArr[i])-1] > 0 {
						eID = extArr[i][len(extArr[i])-1]
					}
					// check if has cycle path
					var OK bool
					if extArr[i], OK = CheckAndCutLocalCycle(extArr[i]); OK {
						break
					}
					// if IsInDBG_MAX_INTArr(extArr[i][:len(extArr[i])-1], eID) {
					// 	break
					// }
					// refresh next info
					if edgesArr[eID].StartNID == nID {
						nID = edgesArr[eID].EndNID
					} else {
						nID = edgesArr[eID].StartNID
					}

					// get next n and m
					// if n > 0 && m >= 0 {
					// 	n, m = GetNextShortPath(consistenceA, n, m)
					// }
					// if n == 0 {
					// 	n, m = GetNextConsisPath(edgesArr, nodesArr, extArr[i], consistenceA)
					// }
				}
			}

			// cat left and  right path
			var catArr []constructdbg.DBG_MAX_INT
			for i := len(extArr[0]) - 1; i >= 0; i-- {
				catArr = append(catArr, extArr[0][i])
			}
			catArr = append(catArr, edge.ID)
			catArr = append(catArr, extArr[1]...)
			fmt.Printf("[paraFindLongMappingPath] catArr: %v\n", catArr)
			fmt.Printf("[paraFindLongMappingPath] mapReg: %v\n", mapReg)
			gapRegArr = AddGapRegion(gapRegArr, mapReg)
			if len(catArr) >= MIN_PATH_LEN {
				wc <- catArr
			}
			if (len(extArr[0]) > 0 && extArr[0][len(extArr[0])-1] == 0) || (len(extArr[1]) > 0 && extArr[1][len(extArr[1])-1] == 0) {
				for i, te := range laArr {
					fmt.Printf("[paraFindLongMappingPath] laArr[%d]: %v\n", i, te)
				}
			}
		}

		/*if nla.RefID == 0 {
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
		      } */

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

func WriteLongPathToDBG(wc chan []constructdbg.DBG_MAX_INT, edgesArr []constructdbg.DBGEdge, numCPU int) [][]constructdbg.Path {
	LongPathMatArr := make([][]constructdbg.Path, len(edgesArr))
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
			// if IsCyclePath(path) {
			// 	continue
			// }

			// edgesArr[path[0]].InsertPathToEdge(path, 1)
		}
		if len(path) > 3 {
			if path[0] == constructdbg.DBG_MAX_INT(0) {
				path = path[1:]
			}
			if path[len(path)-1] == constructdbg.DBG_MAX_INT(0) {
				path = path[:len(path)-1]
			}
		}
		if len(path) > 3 {
			LongPathMatArr[path[0]] = constructdbg.InsertPathToEdge(LongPathMatArr[path[0]], path[1:], 1)
			path = constructdbg.ReverseDBG_MAX_INTArr(path)
			LongPathMatArr[path[0]] = constructdbg.InsertPathToEdge(LongPathMatArr[path[0]], path[1:], 1)
			// fmt.Printf("[WriteLongPathToDBG] LongPathMatArr[%d]: %v\n", path[0], LongPathMatArr[path[0]])
			// fmt.Fprintf(os.Stderr, "[WriteLongPathToDBG] pathArr: %v\n", path)
			totalPathNum++
			totalPathLen += len(path)
		}
	}

	fmt.Printf("[WriteLongPathToDBG] total path num: %d, total path length: %d\n", totalPathNum, totalPathLen)

	return LongPathMatArr
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

}

/*
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
} */

func GetEIDLenArr(edgesArr []constructdbg.DBGEdge, minLen int) (idLenArr []IDLen) {
	for _, edge := range edgesArr {
		if len(edge.Utg.Ks) < minLen {
			continue
		}
		var idLen IDLen
		idLen.EID = edge.ID
		idLen.ELen = int32(len(edge.Utg.Ks))
		idLenArr = append(idLenArr, idLen)
	}

	return idLenArr
}

type IDLenA []IDLen

func (arr IDLenA) Len() int {
	return len(arr)
}

func (arr IDLenA) Less(i, j int) bool {
	return arr[i].ELen < arr[j].ELen
}

func (arr IDLenA) Swap(i, j int) {
	arr[i], arr[j] = arr[j], arr[i]
}

func GetSameDirectionPath(pathArr []constructdbg.Path, eID constructdbg.DBG_MAX_INT) (sdPathArr []constructdbg.Path) {
	for _, p := range pathArr {
		if p.IDArr[0] == eID {
			sdPathArr = append(sdPathArr, p)
		}
	}
	return sdPathArr
}

func AllowZeroDeepEqual(arr1, arr2 []constructdbg.DBG_MAX_INT) bool {
	if arr1[0] != arr2[0] {
		return false
	}
	for i, _ := range arr1 {
		if arr1[i] != 0 && arr2[i] != 0 && arr1[i] != arr2[i] {
			return false
		}
	}

	return true
}

func GetOverlapArr(arr []constructdbg.DBG_MAX_INT, pathArr []constructdbg.Path) (overlapMat []constructdbg.Path) {
	for _, p := range pathArr {
		if len(p.IDArr) > len(arr) && AllowZeroDeepEqual(arr, p.IDArr[:len(arr)]) {
			overlapMat = append(overlapMat, p)
		}
	}

	return overlapMat
}

func GetTotalFreq(pathArr []constructdbg.Path, idx int, eID constructdbg.DBG_MAX_INT) (freqSum int) {
	// fmt.Printf("[GetTotalFreq] idx: %v\n", idx)
	for _, p := range pathArr {
		if p.IDArr[idx] == eID {
			freqSum += p.Freq
		}
	}

	return freqSum
}

func CheckSignificantDiff(scoreArr []int, eIDArr []constructdbg.DBG_MAX_INT) (eID constructdbg.DBG_MAX_INT, ok bool) {
	tmp := make([]int, len(scoreArr))
	copy(tmp, scoreArr)
	sort.Ints(tmp)
	tl := len(tmp)
	if tmp[tl-1] > tmp[tl-2]*2 {
		for i, _ := range scoreArr {
			if scoreArr[i] == tmp[tl-1] {
				eID = eIDArr[i]
				break
			}
		}
		ok = true
	}

	return eID, ok

}

func GetIndexIntSlice(arr []int, ele int) int {
	idx := -1
	for i, id := range arr {
		if id == ele {
			idx = i
			break
		}
	}

	return idx
}

func GetMinCutPathArr(pathArr []constructdbg.Path, minPathFreq int) (consisPathArr []constructdbg.DBG_MAX_INT) {
	for i := 0; ; i++ {
		var eIDArr []constructdbg.DBG_MAX_INT
		var freqArr []int
		for _, p := range pathArr {
			if len(p.IDArr) > i {
				if idx := GetIndexDBG_MAX_INTArr(eIDArr, p.IDArr[i]); idx >= 0 {
					freqArr[idx] += p.Freq
				} else if p.IDArr[i] != 0 {
					eIDArr = append(eIDArr, p.IDArr[i])
					freqArr = append(freqArr, p.Freq)
				}
			}
		}
		if len(eIDArr) == 0 {
			break
		} else if len(eIDArr) == 1 {
			if freqArr[0] >= minPathFreq {
				consisPathArr = append(consisPathArr, eIDArr[0])
			} else {
				break
			}
		} else {
			copyArr := make([]int, len(freqArr))
			copy(copyArr, freqArr)
			sort.Ints(copyArr)
			max, sec := copyArr[len(copyArr)-1], copyArr[len(copyArr)-2]
			if max >= minPathFreq && sec <= minPathFreq && max >= 2*sec {
				idx := GetIndexIntSlice(freqArr, max)
				consisPathArr = append(consisPathArr, eIDArr[idx])
			} else {
				if max >= minPathFreq {
					consisPathArr = nil
				}
				break
			}
		}
	}

	return consisPathArr
}

/*func GetMinCutPathArr(pathArr []constructdbg.Path, minPathFreq int) (consisPathArr []constructdbg.Path) {
	var freqArr [][]int
	for _, p := range pathArr {
		var added bool
		for j, cp := range consisPathArr {
			if len(cp.IDArr) <= len(p.IDArr) {
				if reflect.DeepEqual(cp.IDArr, p.IDArr[:len(cp.IDArr)]) {
					consisPathArr[j] = p
					nf := make([]int, len(p.IDArr))
					copy(nf[:len(cp.IDArr)], freqArr[j])
					nf[len(p.IDArr)-1] += p.Freq
					freqArr[j] = nf
					added = true
					break
				}
			} else {
				if reflect.DeepEqual(p.IDArr, cp.IDArr[:len(p.IDArr)]) {
					freqArr[j][len(p.IDArr)-1] += p.Freq
					added = true
					break
				}
			}
		}
		if added == false {
			consisPathArr = append(consisPathArr, p)
			f := make([]int, len(p.IDArr))
			f[len(f)-1] = p.Freq
			freqArr = append(freqArr, f)
		}
	}

	for _, fa := range freqArr {
		sum := 0
		for j := len(fa) - 1; j >= 0; j-- {
			if fa[j] > 0 {
				sum += fa[j]
			}
			fa[j] = sum
		}
	}
	if len(consisPathArr) == 1 {
		if freqArr[0][0] >= minPathFreq {
			i := 0
			for ; i < len(freqArr[0]); i++ {
				if freqArr[0][i] < minPathFreq {
					break
				}
			}
			consisPathArr[0].IDArr = consisPathArr[0].IDArr[:i]
			consisPathArr[0].Freq = freqArr[0][0]
		} else {
			consisPathArr = consisPathArr[0:0]
		}
	} else if len(consisPathArr) > 1 {
		num := 0
		totalFreq := 0
		for _, f := range freqArr {
			if f[0] >= minPathFreq {
				num++
			}
			totalFreq += f[0]
		}
		if num == 0 {
			if totalFreq > minPathFreq {
				var mergeArr []constructdbg.DBG_MAX_INT
				var stop bool
				j := 0
				for ; j < len(consisPathArr[0].IDArr); j++ {
					freq := freqArr[0][j]
					for x := 1; x < len(consisPathArr); x++ {
						if len(consisPathArr[x].IDArr) <= j {
							stop = true
							break
						}
						if consisPathArr[0].IDArr[j] != consisPathArr[x].IDArr[j] {
							stop = true
							break
						}
					}
					if stop {
						break
					}
					if freq >= minPathFreq {
						mergeArr = append(mergeArr, consisPathArr[0].IDArr[j])
					}
				}
				consisPathArr = consisPathArr[:1]
				consisPathArr[0].IDArr = mergeArr
			} else {
				consisPathArr = consisPathArr[0:0]
			}
		} else if num == 1 {
			for j := 0; j < len(freqArr); j++ {
				if freqArr[j][0] >= minPathFreq {
					x := 0
					for ; x < len(freqArr[j]); x++ {
						if freqArr[j][x] < minPathFreq {
							break
						}
					}
					consisPathArr[j].IDArr = consisPathArr[j].IDArr[:x]
					consisPathArr = consisPathArr[j : j+1]
					break
				}
			}
		} else {
			var mergePath []constructdbg.DBG_MAX_INT
			j := 0
			for {
				idty := true
				var eID constructdbg.DBG_MAX_INT
				for x := 0; x < len(consisPathArr); x++ {
					if j < len(consisPathArr[x].IDArr) && freqArr[x][j] >= minPathFreq {
						if eID == 0 {
							eID = consisPathArr[x].IDArr[j]
						} else {
							if consisPathArr[x].IDArr[j] != eID {
								idty = false
								break
							}
						}
					}
				}
				if eID == 0 {
					break
				}
				if idty {
					mergePath = append(mergePath, eID)
				} else {
					mergePath = nil
					break
				}
			}

			if len(mergePath) > 0 {
				consisPathArr[0].IDArr = mergePath
				consisPathArr = consisPathArr[:1]
			}
		}
	}

	return consisPathArr
} */

func GetIndexDBG_MAX_INTArr(arr []constructdbg.DBG_MAX_INT, ele constructdbg.DBG_MAX_INT) int {
	idx := -1
	for i, id := range arr {
		if id == ele {
			idx = i
			break
		}
	}
	return idx
}

func GetCrossPathMat(neighbourEIDArr []constructdbg.DBG_MAX_INT, LongPathMatArr [][]constructdbg.Path, eID constructdbg.DBG_MAX_INT) (crossPathA []constructdbg.Path) {

	for _, id := range neighbourEIDArr {
		for _, p := range LongPathMatArr[id] {
			idx := GetIndexDBG_MAX_INTArr(p.IDArr, eID)
			if idx >= 0 {
				crossPathA = append(crossPathA, p)
			}
		}
	}

	return crossPathA
}

func GetFlankPath(crossPathA []constructdbg.Path, eID, neID constructdbg.DBG_MAX_INT) (fpArr []constructdbg.Path) {
	for _, p := range crossPathA {
		idx := GetIndexDBG_MAX_INTArr(p.IDArr, eID)
		if idx >= 0 {
			var np constructdbg.Path
			if idx > 0 {
				if p.IDArr[idx-1] == neID {
					np.IDArr = GetReverseDBG_MAX_INTArr(p.IDArr[:idx])
					np.Freq = p.Freq
					fpArr = append(fpArr, np)
				}
			}
			if idx < len(p.IDArr)-1 {
				if p.IDArr[idx+1] == neID {
					np.IDArr = p.IDArr[idx+1:]
					np.Freq = p.Freq
					fpArr = append(fpArr, np)
				}
			}
		}
	}

	return fpArr
}

/*func ParseUniquePath(edgesArr []constructdbg.DBGEdge, nodesArr []constructdbg.DBGNode, LongPathMatArr [][]constructdbg.Path) (mergeMat [][]constructdbg.DBG_MAX_INT) {
	IDLenArr := GetEIDLenArr(edgesArr, 2*Kmerlen-10)
	sort.Sort(sort.Reverse(IDLenA(IDLenArr)))

	// traverse edge from longest to shorter
	for _, IDLen := range IDLenArr {
		edge := edgesArr[IDLen.EID]
		if edge.ID == 0 || edge.GetProcessFlag() > 0 {
			continue
		}
		if IsUniqueEdge(edge, nodesArr) == false {
			continue
		}
		fmt.Printf("[ParseUniquePath] edge[%d] length %d\n", edge.ID, len(edge.Utg.Ks))

		// get crossPathInfo
		var crossPathA [2][]constructdbg.Path
		if edge.StartNID > 0 {
			neighbourEIDArr := GetNeighbourEID(edge.ID, edge.StartNID, edgesArr, nodesArr, MAX_LONG_READ_LEN)
			fmt.Printf("[ParseUniquePath] neighbourEIDArr:%v\n", neighbourEIDArr)
			crossPathA = GetCrossPathMat(neighbourEIDArr, LongPathMatArr, edge.ID)
			cs0 := GetMinCutPathArr(crossPathA[0], 2)
			cs1 := GetMinCutPathArr(crossPathA[1], 2)
			if len(crossPathA[0]) > 0 && len(cs0) == 0 {
				continue
			}
			if len(crossPathA[1]) > 0 && len(cs1) == 0 {
				continue
			}
			fmt.Printf("[ParseUniquePath] crossPathA:%v\n", crossPathA)

		}
		// extend two flank
		var pathArr [2][]constructdbg.DBG_MAX_INT
		for j := 0; j < 2; j++ {
			// var d uint8 // direction
			var neID constructdbg.DBG_MAX_INT
			if j == 0 {
				// d = constructdbg.FORWARD
				if edge.EndNID == 0 {
					continue
				}
				neID = GetNextEID(edge.ID, nodesArr[edge.EndNID])

			} else {
				// d = constructdbg.BACKWARD
				if edge.StartNID == 0 {
					continue
				}
				neID = GetNextEID(edge.ID, nodesArr[edge.StartNID])
			}
			// if len(LongPathMatArr[edge.ID]) > 0 {
			// 	fmt.Printf("[ParseUniquePath] LongPathMatArr: %v\n", LongPathMatArr[edge.ID])
			// }
			dPathArr := GetSameDirectionPath(LongPathMatArr[edge.ID], neID)
			dPathArr = append(dPathArr, crossPathA[j]...)
			for x := 0; x < len(dPathArr); x++ {
				fmt.Printf("[ParseUniquePath] dPathArr[%d]: %v\n", x, dPathArr[x])
			}
			// consisPathMat := GetMinCutPathArr(dPathArr, MIN_ALLOW_PATH_FREQ)
			consisPathArr := GetMinCutPathArr(dPathArr, 2)
			fmt.Printf("[ParseUniquePath] consisPathArr: %v\n", consisPathArr)
			if len(consisPathArr) < MIN_PATH_LEN {
				continue
			}
			pathArr[j] = consisPathArr
			for {
				leID := pathArr[j][len(pathArr[j])-1]
				nID := GetInterNodeID(edgesArr[pathArr[j][len(pathArr[j])-2]], edgesArr[leID])
				if nID == edgesArr[leID].StartNID {
					nID = edgesArr[leID].EndNID
				} else {
					nID = edgesArr[leID].StartNID
				}
				if nID == 0 {
					break
				}
				eIDArr := GetNextEArr(leID, nodesArr[nID])
				if len(eIDArr) == 0 {
					break
				} else if len(eIDArr) == 1 {
					pathArr[j] = append(pathArr[j], eIDArr[0])
				} else {
					minLen := 1000
					scoreArr := make([]int, len(eIDArr))
					extLen := len(edgesArr[leID].Utg.Ks) - (Kmerlen - 1)
					for x := len(pathArr[j]) - 2; x >= 0 && extLen < MAX_LONG_READ_LEN; x-- {
						eID := pathArr[j][x]
						if extLen > minLen {
							overlapPathArr := GetOverlapArr(pathArr[j][x+1:], LongPathMatArr[eID])
							for y, id := range eIDArr {
								scoreArr[y] += GetTotalFreq(overlapPathArr, len(pathArr[j])-(x+1)-1, id) * extLen
							}
						}
						extLen += len(edgesArr[eID].Utg.Ks) - (Kmerlen - 1)

					}
					if id, OK := CheckSignificantDiff(scoreArr, eIDArr); OK {
						pathArr[j] = append(pathArr[j], id)
					} else {
						break
					}
				}
			}
		}
		// merge two flank

		catArr := GetReverseDBG_MAX_INTArr(pathArr[1])
		catArr = append(catArr, edge.ID)
		catArr = append(catArr, pathArr[0]...)
		for _, id := range catArr {
			edgesArr[id].SetProcessFlag()
		}
		if len(catArr) >= MIN_PATH_LEN {
			mergeMat = append(mergeMat, catArr)
		}
		fmt.Printf("[ParseUniquePath] catArr: %v\n", catArr)

	}

	return mergeMat
} */

func ParseUniquePath(edgesArr []constructdbg.DBGEdge, nodesArr []constructdbg.DBGNode, LongPathMatArr [][]constructdbg.Path) (mergeMat [][]constructdbg.DBG_MAX_INT) {
	IDLenArr := GetEIDLenArr(edgesArr, 2*Kmerlen-10)
	sort.Sort(sort.Reverse(IDLenA(IDLenArr)))
	minDepth := 2

	// traverse edge from longest to shorter
	for _, IDLen := range IDLenArr {
		edge := edgesArr[IDLen.EID]
		if edge.ID == 0 || edge.GetProcessFlag() > 0 {
			continue
		}
		if IsUniqueEdge(edge, nodesArr) == false {
			continue
		}
		fmt.Printf("[ParseUniquePath] edge[%d] length %d\n", edge.ID, len(edge.Utg.Ks))

		// get crossPathInfo
		var crossPathA []constructdbg.Path
		if edge.StartNID > 0 {
			neighbourEIDArr := GetNeighbourEID(edge.ID, edge.StartNID, edgesArr, nodesArr, MAX_LONG_READ_LEN-len(edge.Utg.Ks))
			// fmt.Printf("[ParseUniquePath] neighbourEIDArr:%v\n", neighbourEIDArr)
			crossPathA = GetCrossPathMat(neighbourEIDArr, LongPathMatArr, edge.ID)
			fmt.Printf("[ParseUniquePath] crossPathA:%v\n", crossPathA)
		}
		// extend two flank
		var pathArr [2][]constructdbg.DBG_MAX_INT
		for j := 0; j < 2; j++ {
			// var d uint8 // direction
			var neID constructdbg.DBG_MAX_INT
			if j == 0 {
				fmt.Printf("[ParseUniquePath] direction : FORWARD\n")
				// d = constructdbg.FORWARD
				if edge.EndNID == 0 {
					continue
				}
				neID = GetNextEID(edge.ID, nodesArr[edge.EndNID])

			} else {
				fmt.Printf("[ParseUniquePath] direction : BACKWARD\n")
				// d = constructdbg.BACKWARD
				if edge.StartNID == 0 {
					continue
				}
				neID = GetNextEID(edge.ID, nodesArr[edge.StartNID])
			}
			// if len(LongPathMatArr[edge.ID]) > 0 {
			// 	fmt.Printf("[ParseUniquePath] LongPathMatArr: %v\n", LongPathMatArr[edge.ID])
			// }
			dPathArr := GetSameDirectionPath(LongPathMatArr[edge.ID], neID)
			dPathArr = append(dPathArr, GetFlankPath(crossPathA, edge.ID, neID)...)
			for x := 0; x < len(dPathArr); x++ {
				fmt.Printf("[ParseUniquePath] dPathArr[%d]: %v\n", x, dPathArr[x])
			}
			// consisPathMat := GetMinCutPathArr(dPathArr, MIN_ALLOW_PATH_FREQ)
			consisPathArr := GetMinCutPathArr(dPathArr, minDepth)
			fmt.Printf("[ParseUniquePath] consisPathArr: %v\n", consisPathArr)
			if len(consisPathArr) < MIN_PATH_LEN {
				continue
			}
			pathArr[j] = consisPathArr
			y := 0
			var pileArr [][]constructdbg.Path
			pileArr = append(pileArr, dPathArr)
			for x := 0; x < len(pathArr[j]); x++ {
				var overlapPathArr []constructdbg.Path
				if x == len(pathArr[j])-1 {
					for _, p := range LongPathMatArr[pathArr[j][x]] {
						if p.IDArr[0] != pathArr[j][x-1] {
							overlapPathArr = append(overlapPathArr, p)
						}
					}
				} else {
					overlapPathArr = GetOverlapArr(pathArr[j][x+1:], LongPathMatArr[pathArr[j][x]])
				}
				for _, opa := range overlapPathArr {
					fmt.Printf("[ParseUniquePath] overlap: %v\n", opa)
				}
				pileArr = append(pileArr, overlapPathArr)
				if len(overlapPathArr) == 0 {
					continue
				}
				for {
					var pa []constructdbg.Path
					for z := y; z <= x+1; z++ {
						for w := 0; w < len(pileArr[z]); w++ {
							if len(pileArr[z][w].IDArr) > len(pathArr[j])-z {
								var p constructdbg.Path
								p.IDArr = append(p.IDArr, pileArr[z][w].IDArr[len(pathArr[j])-z])
								p.Freq = pileArr[z][w].Freq
								if p.IDArr[0] > 0 {
									pa = constructdbg.InsertPathToEdge(pa, p.IDArr, p.Freq)
								}
							}
						}
						if len(pa) == 0 {
							y++
						}
					}
					if len(pa) == 0 {
						break
					}
					fmt.Printf("[ParseUniquePath] x: %d, pa: %v\n", x, pa)
					if len(pa) == 1 {
						if pa[0].Freq >= minDepth {
							pathArr[j] = append(pathArr[j], pa[0].IDArr[0])
						} else {
							break
						}
					} else if len(pa) == 2 {
						var max, min int
						if pa[0].Freq < pa[1].Freq {
							max = pa[1].Freq
							min = pa[0].Freq
						} else {
							max = pa[0].Freq
							min = pa[1].Freq
						}
						if max > minDepth && min <= minDepth && max >= 3*min {
							if max == pa[0].Freq {
								pathArr[j] = append(pathArr[j], pa[0].IDArr[0])
							} else {
								pathArr[j] = append(pathArr[j], pa[1].IDArr[0])
							}
						} else {
							break
						}
					} else {
						break
						// log.Fatalf("[ParseUniquePath] len(pa) == %d case not processed\n", len(pa))
					}
				}
				// fmt.Printf("[ParseUniquePath] pathArr[%d]: %v\n", j, pathArr[j])

			}
			/*for {
				fmt.Printf("[ParseUniquePath]pathArr[%d]: %v\n", j, pathArr[j])
				leID := pathArr[j][len(pathArr[j])-1]
				nID := GetInterNodeID(edgesArr[pathArr[j][len(pathArr[j])-2]], edgesArr[leID])
				if nID == 0 {
					pathArr[j] = pathArr[j][:len(pathArr[j])-1]
					break
				}
				if nID == edgesArr[leID].StartNID {
					nID = edgesArr[leID].EndNID
				} else {
					nID = edgesArr[leID].StartNID
				}
				if nID == 0 {
					break
				}
				eIDArr := GetNextEArr(leID, nodesArr[nID])
				if len(eIDArr) == 0 {
					break
				} else if len(eIDArr) == 1 {
					pathArr[j] = append(pathArr[j], eIDArr[0])
				} else {
					minLen := 1000
					scoreArr := make([]int, len(eIDArr))
					extLen := len(edgesArr[leID].Utg.Ks) - (Kmerlen - 1)
					for x := len(pathArr[j]) - 2; x >= 0 && extLen < MAX_LONG_READ_LEN; x-- {
						eID := pathArr[j][x]
						if extLen > minLen {
							// for _, opa := range overlapPathArr {
							// 	fmt.Printf("[ParseUniquePath] overlap: %v\n", opa)
							// }
							overlapPathArr := GetOverlapArr(pathArr[j][x+1:], LongPathMatArr[eID])
							for _, opa := range overlapPathArr {
								fmt.Printf("[ParseUniquePath] overlap: %v\n", opa)
							}
							for y, id := range eIDArr {
								scoreArr[y] += GetTotalFreq(overlapPathArr, len(pathArr[j])-(x+1), id) * extLen
							}
						}
						extLen += len(edgesArr[eID].Utg.Ks) - (Kmerlen - 1)

					}
					if id, OK := CheckSignificantDiff(scoreArr, eIDArr); OK {
						pathArr[j] = append(pathArr[j], id)
					} else {
						fmt.Printf("[ParseUniquePath] eIDArr: %v\nscoreArr: %v\n", eIDArr, scoreArr)
						break
					}
					fmt.Printf("[ParseUniquePath] pathArr: %v\n", pathArr[j])
				}
			} */
		}
		// merge two flank

		catArr := GetReverseDBG_MAX_INTArr(pathArr[1])
		catArr = append(catArr, edge.ID)
		catArr = append(catArr, pathArr[0]...)
		for _, id := range catArr {
			edgesArr[id].SetProcessFlag()
		}
		if len(catArr) >= MIN_PATH_LEN {
			mergeMat = append(mergeMat, catArr)
		}
		fmt.Printf("[ParseUniquePath] catArr: %v\n", catArr)

	}

	return mergeMat
}

func IsBubbleEdge(edge constructdbg.DBGEdge, nodesArr []constructdbg.DBGNode, edgesArr []constructdbg.DBGEdge) (b bool) {
	if edge.StartNID > 0 {
		eIDArr := GetOtherEArr(nodesArr[edge.StartNID], edge.ID)
		if len(eIDArr) == 1 {
			e2 := edgesArr[eIDArr[0]]
			if e2.StartNID == edge.StartNID && e2.EndNID == edge.EndNID {
				b = true
			} else if e2.EndNID == edge.StartNID && e2.StartNID == edge.EndNID {
				b = true
			}
		}
	}

	return b
}

func StaticsLongPath(LongPathMatArr [][]constructdbg.Path, edgesArr []constructdbg.DBGEdge, nodesArr []constructdbg.DBGNode) {
	freqNum := make([]int, len(edgesArr))
	num := 0
	totalBaseLen := 0
	for _, pathMat := range LongPathMatArr {
		if len(pathMat) > 0 {
			for _, p := range pathMat {
				initLen := Kmerlen - 1
				for _, eID := range p.IDArr {
					initLen += len(edgesArr[eID].Utg.Ks) - (Kmerlen - 1)
					freqNum[eID] += p.Freq
				}
				totalBaseLen += initLen * p.Freq
				num += p.Freq
			}
		}
	}
	fmt.Printf("[StaticsLongPath]LongPathMatArr has path num: %d, avgBaseLen: %d\n", num, totalBaseLen/num)

	// statics freq
	max := 256
	bubbleNum := 0
	freqdistri := make([]int, max)
	for i, freq := range freqNum {
		if freq >= max {
			freqdistri[max-1]++
		} else if freq > 0 {
			freqdistri[freq]++
		} else {
			freqdistri[freq]++
			if IsBubbleEdge(edgesArr[i], nodesArr, edgesArr) {
				bubbleNum++
				// arr := GetOtherEArr(nodesArr[edgesArr[i].StartNID], edgesArr[i].ID)
				// fmt.Fprintf(os.Stderr, "%v\n", freqNum[arr[0]])
			} else {
				if edgesArr[i].ID > 0 {
					// fmt.Fprintf(os.Stderr, "eID: %d, len: %v\n", i, len(edgesArr[i].Utg.Ks))

				}
			}
		}
	}
	fmt.Printf("[StaticsLongPath] bubbleNum: %d\n", bubbleNum)
	// for _, depth := range freqdistri {
	// 	fmt.Fprintf(os.Stderr, "%d\n", depth)
	// }
}

func StaticsMergeMat(matArr [][]constructdbg.DBG_MAX_INT, edgesArr []constructdbg.DBGEdge) {

	num := 0
	totalBaseLen := 0
	for _, p := range matArr {
		initLen := Kmerlen - 1
		for _, eID := range p {
			initLen += len(edgesArr[eID].Utg.Ks) - (Kmerlen - 1)
		}
		// fmt.Fprintf(os.Stderr, "%d\n", initLen)
		totalBaseLen += initLen
		num++
	}
	fmt.Printf("[StaticsMergeMat]MergeMat has path num: %d, avgBaseLen: %d\n", num, totalBaseLen/num)
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
	nodesStatfn := prefix + ".nodes.stat"
	nodesSize := constructdbg.NodesStatReader(nodesStatfn)
	nodesArr := make([]constructdbg.DBGNode, nodesSize)
	smfyNodesfn := prefix + ".nodes.smfy.Arr"
	nodesArr = constructdbg.NodesArrReader(smfyNodesfn)
	// constructdbg.NodeMap2NodeArr(nodeMap, nodesArr)

	// Restore edges info
	edgesStatfn := prefix + ".edges.stat"
	edgesSize := constructdbg.EdgesStatReader(edgesStatfn)
	edgesArr := make([]constructdbg.DBGEdge, edgesSize)
	edgesfn := prefix + ".edges.smfy.fq"
	constructdbg.LoadEdgesfqFromFn(edgesfn, edgesArr, true)

	// get coverage of smfy edge
	computeCoverageSmfyEdge(prefix)

	// get short read mapping path
	lastfn := prefix + ".last"
	rc := make(chan [2]FastMapRecord, numCPU*2)
	wc := make(chan []constructdbg.DBG_MAX_INT, numCPU*2)
	// defer rc.Close()
	// go GetFastMapRecord(fastmapfn, rc, numCPU)
	go GetLastMapRecord(lastfn, rc, numCPU)

	for i := 0; i < numCPU; i++ {
		go paraFindShortMappingPath(rc, wc, edgesArr, nodesArr)
	}

	WriteShortPathToDBG(wc, edgesArr, numCPU)

	SortPathMat(edgesArr)

	// for _, e := range edgesArr {
	// 	fmt.Printf("[Fpath] pathMat: %v\n", e.PathMat)
	// }

	// Find Max unique path and merge neighbour edges
	// fmt.Printf("[FSpath] edgesArr[76]: %v\n", edgesArr[76])
	consistenceA := FindMaxUnqiuePath(edgesArr, nodesArr)

	// get long read mapping path
	LongReadPathfn := prefix + ".LA"
	lac := make(chan []LA, numCPU*2)
	wc = make(chan []constructdbg.DBG_MAX_INT, numCPU*2)
	RefIDMapArr := GetOrderID(edgesArr)
	// for i, eID := range RefIDMapArr {
	// 	fmt.Printf("[Fpath]idx: %d  eID: %d\n", i, eID)
	// }
	go GetLARecord(LongReadPathfn, RefIDMapArr, lac, numCPU)

	for i := 0; i < numCPU; i++ {
		go paraFindLongMappingPath(lac, wc, edgesArr, nodesArr, RefIDMapArr, consistenceA)
	}

	LongPathMatArr := WriteLongPathToDBG(wc, edgesArr, numCPU)
	StaticsLongPath(LongPathMatArr, edgesArr, nodesArr)
	// fmt.Printf("[Fpath]sizeof(LongPathMatArr): %d\n", unsafe.Sizeof(LongPathMatArr))
	mergeMat := ParseUniquePath(edgesArr, nodesArr, LongPathMatArr)
	StaticsMergeMat(mergeMat, edgesArr)
	// graphfn := prefix + ".LongPath.dot"
	// GraphvizDBG(nodesArr, edgesArr, graphfn)
	// Write to files
	edgesfn = prefix + ".edges.LongPath.fq"
	constructdbg.StoreEdgesToFn(edgesfn, edgesArr, false)
	// constructdbg.StoreEdgesToFn(edgesfn, edgesArr, false)
	nodesfn := prefix + ".nodes.LongPath.Arr"
	constructdbg.NodesArrWriter(nodesArr, nodesfn)
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
