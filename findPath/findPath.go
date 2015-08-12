package findPath

import (
	// "bufio"
	// "compress/gzip"
	// "container/list"
	// "encoding/binary"
	// "encoding/gob"
	"fmt"
	// "ga/bnt"
	// "ga/constructcf"
	"ga/constructdbg"
	// "ga/cuckoofilter"
	// "io"
	"log"
	"os"
	// "runtime"
	"strconv"

	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/sam"
	"github.com/jwaldrip/odin/cli"
	// "strings"
)

var Kmerlen int

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
	for {
		r, err := bamfp.Read()
		if err != nil {
			break
		}
		if len(r.Cigar) <= 1 {
			continue
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

func Cigar2String(cigar sam.Cigar) (cs string) {
	for _, v := range cigar {
		cs += strconv.Itoa(v.Len()) + "---"
	}
	return
}

func paraFindShortMappingPath(rc chan []sam.Record, wc chan []constructdbg.DBG_MAX_INT, edgesArr []constructdbg.DBGEdge, nodesArr []constructdbg.DBGNode) {
	totalM, totalI, totalD := 0, 0, 0
	for {
		rArr := <-rc
		if len(rArr) == 0 {
			var guardPath []constructdbg.DBG_MAX_INT
			wc <- guardPath
			break
		}

		fmt.Printf("[paraFindShortMappingPath]ReadID: %s\n", rArr[0].Name)
		for _, v := range rArr[:1] {
			var NM, MD, SA sam.Tag
			NM[0] = 'N'
			NM[1] = 'M'
			MD[0] = 'M'
			MD[1] = 'D'
			SA[0] = 'S'
			SA[1] = 'A'
			nm := v.AuxFields.Get(NM).Value()
			var nmv int
			switch nm.(type) {
			case uint8:
				nmv = int(nm.(uint8))
			case uint16:
				nmv = int(nm.(uint16))
			}
			// nmv, err := strconv.Atoi(nm.String()[5:])
			// if err != nil {
			// 	log.Fatalf("[paraFindPacbioMappingPath] convert 'NM' err: %v\n", nmv)
			// }
			// fmt.Printf("[paraFindPacbioMappingPath]nm type: %v\n", v.AuxFields.Get(NM).Type())
			// if nm != nil {
			// 	// snm = nm.Value()
			// } else {
			// 	log.Fatalf("[paraFindPacbioMappingPath] nm: %v\n", v.Ref)
			// }
			Mnum, Inum, Dnum := AccumulateCigar(v.Cigar)
			if Mnum < 1000 {
				continue
			}
			fmt.Printf("ref\tpos\tmapQ\tNM\tCigar\n")
			totalM += Mnum
			totalI += Inum
			totalD += Dnum
			totalMis += (int(nmv) - Inum - Dnum)
			// as := Cigar2String(acgr)
			fmt.Printf("%s\t%d\t%v\t%v\n", v.Ref.Name(), v.Pos, v.MapQ, v.Cigar)
			// sav := v.AuxFields.Get(SA)
			// // fmt.Printf("%v\n", sav)
			// if sav != nil {
			// 	sa := v.AuxFields.Get(SA).Value().(string)
			// 	// fmt.Printf("%v\n", sa)
			// 	// sav := strings.Split(sa, ":")
			// 	// fmt.Printf("%v\n", sav)
			// 	saArr := strings.Split(sa[:len(sa)-1], ";")
			// 	for _, e := range saArr {
			// 		arr := strings.Split(e, ",")
			// 		// fmt.Printf("%v\n", arr)
			// 		// cgr, _ := sam.ParseCigar([]byte(arr[3]))
			// 		// cr := AccumulateCigar(cgr)
			// 		// as := Cigar2String(cr)
			// 		fmt.Printf("%s\t%s\t%s\t%s\t%s\n", arr[0], arr[1], arr[4], arr[5], arr[3])
			// 	}

			// }
		}
	}

	fmt.Printf("[paraFindPacbioMappingPath] totalM: %d, totalMis: %d, totalI: %d, totalD: %d\n", totalM, totalMis, totalI, totalD)
}

func WriteShortPathToDBG(wc chan []constructdbg.DBG_MAX_INT, edgesArr []constructdbg.DBGEdge, nodesArr []constructdbg.DBGNode, numCPU int) {
	terminalNum := 0
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
	}
}

func FSpath(c cli.Command) {
	k := c.Parent().Flag("K").String()
	var err error = nil
	Kmerlen, err = strconv.Atoi(k)
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
		go paraFindShortMappingPath(rc, wc, edgesArr, nodesArr)
	}

	WriteShortPathToDBG(wc, edgesArr, nodesArr, numCPU)
}
