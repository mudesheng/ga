package preprocess

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
	"strconv"
	"strings"

	//"github.com/google/brotli/cbrotli"
	"github.com/jwaldrip/odin/cli"
	"github.com/mudesheng/ga/cbrotli"
	"github.com/mudesheng/ga/constructcf"
	"github.com/mudesheng/ga/constructdbg"
	"github.com/mudesheng/ga/utils"
)

type Options struct {
	utils.ArgsOpt
	CFSize        int64
	TipMaxLen     int
	WinSize       int
	MaxNGSReadLen int
	Correct       bool
}

func checkArgs(c cli.Command) (opt Options, suc bool) {
	tmp, err := strconv.Atoi(c.Flag("S").String())
	if err != nil {
		log.Fatalf("[checkArgs] argument 'S': %v set error: %v\n", c.Flag("S"), err)
	}
	if tmp < 1024*1024 {
		log.Fatalf("the argument 'S': %v must bigger than 1024 * 1024\n", c.Flag("S"))
	}
	opt.CFSize = int64(tmp)
	tmp, err = strconv.Atoi(c.Flag("MaxNGSReadLen").String())
	if err != nil {
		log.Fatalf("[checkArgs] argument 'MaxNGSReadLen': %v set error: %v\n", c.Flag("MaxNGSReadLen"), err)
	}
	opt.MaxNGSReadLen = tmp
	cor := c.Flag("Correct").Get().(bool)
	if cor == true {
		opt.Correct = true
	} else {
		log.Fatalf("[checkArgs] argument 'Correct': %v set error, must set 'true'\n", c.Flag("Correct"))
	}
	opt.MaxNGSReadLen = tmp
	tmp, err = strconv.Atoi(c.Flag("WinSize").String())
	if err != nil {
		log.Fatalf("[checkArgs] argument 'MaxNGSReadLen': %v set error: %v\n", c.Flag("MaxNGSReadLen"), err)
	}
	if tmp < 3 || tmp > 20 {
		log.Fatalf("[checkArgs] argument 'WinSize': %v, must bewteen [3~20] set error: %v\n", c.Flag("WinSize"), err)
	}
	opt.WinSize = tmp
	if 81 < opt.Kmer && opt.Kmer < 149 && opt.Kmer%2 == 0 {
		log.Fatalf("the argument 'K': %v must between [81~149] and tmp must been even\n", c.Parent().Flag("K"))
	}

	tmp, err = strconv.Atoi(c.Flag("tipMaxLen").String())
	if err != nil {
		log.Fatalf("[checkArgs] argument 'tipMaxLen': %v set error: %v\n", c.Flag("tipMaxLen"), err)
	}
	if tmp == 0 {
		opt.TipMaxLen = opt.MaxNGSReadLen
	} else {
		if tmp < 0 || tmp < opt.Kmer*2 || tmp > opt.MaxNGSReadLen {
			log.Fatalf("[checkArgs] argument 'tipMaxLen': %v must between [%v~%v]\n", tmp, opt.Kmer*2, opt.MaxNGSReadLen)
		}
		opt.TipMaxLen = tmp
	}
	suc = true
	return opt, suc
}

func LoadNGSReads(brfn1, brfn2 string, cs chan<- [2]constructcf.ReadInfo, kmerlen int) {
	fp1, err1 := os.Open(brfn1)
	if err1 != nil {
		log.Fatalf("[LoadNGSReads] open file: %v failed..., err: %v\n", brfn1, err1)
	}
	defer fp1.Close()
	fp2, err2 := os.Open(brfn2)
	if err2 != nil {
		log.Fatalf("[LoadNGSReads] open file: %v failed..., err: %v\n", brfn2, err2)
	}
	defer fp2.Close()

	//brfp1 := cbrotli.NewReader(fp1)
	brfp1 := cbrotli.NewReaderSize(fp1, 1<<25)
	defer brfp1.Close()
	brfp2 := cbrotli.NewReaderSize(fp2, 1<<25)
	defer brfp2.Close()
	// brfp, err := brotli.NewReader(fp, nil)
	buffp1 := bufio.NewReader(brfp1) // 1<<24 == 2**24
	buffp2 := bufio.NewReader(brfp2)
	format1 := constructcf.GetReadsFileFormat(brfn1)
	format2 := constructcf.GetReadsFileFormat(brfn2)
	var count int
	var EOF error
	for EOF != io.EOF {
		ri1, err1 := constructcf.GetReadFileRecord(buffp1, format1, false)
		ri2, err2 := constructcf.GetReadFileRecord(buffp2, format2, false)
		if ri1.ID == 0 || ri2.ID == 0 {
			if err1 == err2 && err1 == io.EOF {
				break
			} else {
				log.Fatalf("[LoadNGSReads] file : %v not consis with file : %v\n", brfn1, brfn2)
			}
		}
		if ri1.ID != ri2.ID {
			log.Fatalf("[LoadNGSReads] read1 ID : %v != read2 ID: %v\n", ri1.ID, ri2.ID)
		}
		var pairRI [2]constructcf.ReadInfo
		pairRI[0].ID = ri1.ID
		pairRI[1].ID = ri2.ID
		if len(ri1.Seq) < kmerlen+20 || len(ri2.Seq) < kmerlen+20 {
			continue
		}
		pairRI[0].Seq = ri1.Seq
		pairRI[1].Seq = ri2.Seq
		cs <- pairRI
		count++
	}
	fmt.Printf("[LoadNGSReads] processed %d pair reads from pair files: %s , %s\n", count, brfn1, brfn2)
	close(cs)
}


type ReadMapInfo struct {
	ID           int64
	StartP, EndP int
	//NID          constructdbg.DBG_MAX_INT
	Seq     []byte
	PathSeqArr    []constructdbg.PathSeq
	Strands []bool
}

func GetExtendPathArr(path []constructdbg.DBG_MAX_INT, nID constructdbg.DBG_MAX_INT, edgesArr []constructdbg.DBGEdge, nodesArr []constructdbg.DBGNode, kmerlen, extLen int) (epArr [][]constructdbg.DBG_MAX_INT) {
	type pathExtend struct {
		Path      []constructdbg.DBG_MAX_INT
		NID       constructdbg.DBG_MAX_INT
		ExtendLen int
	}
	var pathQueue []pathExtend
	var pe pathExtend
	pe.Path, pe.NID, pe.ExtendLen = path, nID, extLen
	pathQueue = append(pathQueue, pe)
	for len(pathQueue) > 0 {
		lpe := pathQueue[len(pathQueue)-1]
		pathQueue = pathQueue[:len(pathQueue)-1]
		if lpe.ExtendLen >= extLen {
			epArr = append(epArr, lpe.Path)
			continue
		}

		// add next edge
		ea := constructdbg.GetNearEdgeIDArr(nodesArr[lpe.NID], lpe.Path[len(lpe.Path)-1])
		for _, eID := range ea {
			ne := edgesArr[eID]
			var npe pathExtend
			npe.Path = append(npe.Path, lpe.Path...)
			npe.Path = append(npe.Path, eID)
			if ne.StartNID == lpe.NID {
				npe.NID = ne.EndNID
			} else {
				npe.NID = ne.StartNID
			}
			npe.ExtendLen = lpe.ExtendLen + len(ne.Utg.Ks) - (kmerlen - 1)
			pathQueue = append(pathQueue, npe)
		}
	}

	return
}

func ReversePathSeqArr(psArr []PathSeq) {
	la := len(psArr)
	for i := 0; i < len(psArr)/2; i++ {
		tmp := psArr[i]
		psArr[i].ID, psArr[i].NID, psArr[i].Start, psArr[i].End = psArr[la-i-1].ID, psArr[la-i-1].NID, psArr[la-i-1].End, psArr[la-i-1].Start
		psArr[la-i-1].ID, psArr[la-i-1].NID psArr[la-i-1].Start, psArr[la-i-1].End = tmp.ID, tmp.NID, tmp.End, tmp.Start
	}
	if len(psArr)%2 == 1 {
		i := len(psArr)/2
		psArr[i].Start, psArr[i].End = psArr[i].End, psArr[i].Start
	}
}


func MergePathSeqArr(pathSeqArr1, pathSeqArr2 []PathSeq) (mergePathSeqArr []PathSeq) {
	// reverse pathSeqArr2
	ReversePathSeqArr(pathSeqArr2)
	var sharePathLen int
	for i1, i2 := len(pathSeqArr1)-1, len(pathSeqArr2)-1; i1 >= 0 && i2 >= 0; i2-- {
		if pathSeqArr1[i1].ID == pathSeqArr2[i2].ID && pathSeqArr1[i1].NID == pathSeqArr2[i2].NID {
			if i2 > 0 {
				var spl int
				for j1, j2 := i1-1, i2-1; j1 >= 0 && j2 >= 0; j1, j2 = j1-1, j2-1 {
					if pathSeqArr1[j1].ID == pathSeqArr2[j2].ID && pathSeqArr1[j1].NID == pathSeqArr2[j2].NID {
						spl++
					} else { break }
				}
				if spl == i2 {
					sharePathLen = spl+1
					break
				}
			} else {
				sharePathLen = 1
				break
			}
		}
	}
	if sharePathLen < 1 {
		return
	}

	for i1,i2 := 0,sharePathLen-1; i1 <= len(pathSeqArr1) && i2 < len(pathSeqArr2) && i1 <= i2; {
		if i1 < len(pathSeqArr1) -1 {
				mergePathSeqArr = append(mergePathSeqArr, pathSeqArr1[i1])
				i1++
		} else if i1 == len(pathSeqArr1) -1 {
			var ps PathSeq
			ps.ID = pathSeqArr1[i1].ID
				ps.Start = pathSeqArr1[i1].Start
				ps.End = pathSeqArr2[i2].End
			mergePathSeqArr = append(mergePathSeqArr, ps)
			i1++; i2++
		} else {
			mergePathSeqArr = append(mergePathSeqArr, pathSeqArr1[i2])
			i2++
		}
	}
}


func GetPathSeq(m ReadMapInfo, edgesArr []constructdbg.DBGEdge, nodesArr []constructdbg.DBGNode, kmerlen int) (seq []byte) {
	//fmt.Printf("[GetPathSeq] m: %v\n", m)
	if len(m.PathSeqArr) == 1 {
		if m.StartP > m.EndP {
			seq = constructdbg.GetReverseCompByteArr(edgesArr[m.PathSeqArr[0].ID].Utg.Ks[m.EndP : m.StartP])
		} else {
			seq = edgesArr[m.PathSeqArr[0]].Utg.Ks[m.StartP : m.EndP]
		}
		return seq
	}

	var strand bool // true denote "+" DNA strand map, false denote  "-" DNA strand map
	var coming bool // true denote EdgeOutcoming, false denote EdgeIncoming
	for i := 1; i < len(m.PathSeqArr); i++ {
		ps := m.PathSeqArr[i-1]
		ps1 := m.PathSeqArr[i]
		e := edgesArr[ps.ID]
		e1 := edgesArr[ps1.ID]
			if i == 1 {
			if e.StartNID == e.EndNID { return } // test if self cycle
			if constructdbg.IsInComing(nodesArr[ps.NID].EdgeIDIncoming, e.ID) && constructdbg.IsInComing(nodesArr[ps.NID].EdgeIDOutcoming, e1.ID) {
				coming = true
			} else {
				coming = false
			}
				if e.EndNID == ps.NID && ps.End == len(e.Utg.Ks) {
					strand = true
					seq = e.Utg.Ks[ps.Start:ps.End]
				} else {
					strand = false
					seq = constructdbg.GetReverseCompByteArr(e.Utg.Ks[ps.End+1: ps.Start+1])
				}
		}
		if e1.StartNID == e1.EndNID {
				fmt.Printf("[GetPathSeq] v: %v\ne: %v\ne: %v\n", nodesArr[ps.NID], e, e1)
				if coming {
					if strand {
						seq = append(seq, e1.Utg.Ks[kmerlen-1:ps1.End])
					} else {
						seq = append(seq, constructdbg.GetReverseCompByteArr(e1.Utg.Ks[ps1.End+1:ps1.Start+1-(kmerlen-1)]))
					}
				} else {
					if strand {
						seq = append(seq, constructdbg.GetReverseCompByteArr(e1.Utg.Ks[ps1.End+1:ps1.Start+1-(kmerlen-1)]))
					} else {
						seq = append(seq, e1.Utg.Ks[kmerlen-1:ps1.End])
					}
					strand = !strand
				}
			coming = !coming
		} else {
			if (e.EndNID == ps.NID && e1.StartNID == ps.NID) || (e.StartNID == ps.NID && e1.EndNID == ps.NID) {
				if strand {
					seq = append(seq, e1.Utg.Ks[kmerlen-1:ps1.End])
				} else {
					seq = append(seq, constructdbg.GetReverseCompByteArr(e1.Utg.Ks[ps1.End+1:ps1.Start+1-(kmerlen-1)]))
				}
			} else {
				if strand {
					seq = append(seq, constructdbg.GetReverseCompByteArr(e1.Utg.Ks[ps1.End+1:ps1.Start+1-(kmerlen-1)]))
				} else {
					seq = append(seq, e1.Utg.Ks[kmerlen-1:ps1.End])
				}
				strand = !strand
			}
			if constructdbg.IsInComing(nodesArr[ps.NID].EdgeIDIncoming, e.ID) && constructdbg.IsInComing(nodesArr[ps.NID].EdgeIDOutcoming, e1.ID) {
				coming = true
			} else {
				coming = false
			}
		}

	/*strand := m.Strands[0]
	e0, e1 := edgesArr[m.Path[0]], edgesArr[m.Path[1]]
	//if m.NID == 0 {
	//	log.Fatalf("[GetPathSeq] m.NID: %v not set\n", m.NID)
	//}
	//nID := m.NID
	var nID constructdbg.DBG_MAX_INT
	if strand == constructdbg.PLUS {
		seq = e0.Utg.Ks[m.StartP+1:]
		nID = e0.EndNID
	} else {
		seq = constructdbg.GetReverseCompByteArr(e0.Utg.Ks[:m.StartP])
		nID = e0.StartNID
	}
	for _, eID := range m.Path[1 : len(m.Path)-1] {
		e1 = edgesArr[eID]
		if e1.StartNID == e1.EndNID {
			if strand == constructdbg.PLUS {
				if reflect.DeepEqual(e0.Utg.Ks[len(e0.Utg.Ks)-kmerlen+1:], e1.Utg.Ks[:kmerlen-1]) {

				} else {
					strand = !strand
				}
			} else {
				if reflect.DeepEqual(e0.Utg.Ks[:kmerlen-1], e1.Utg.Ks[len(e1.Utg.Ks)-kmerlen+1:]) {

				} else {
					strand = !strand
				}
			}
		} else {
			if e0.EndNID == nID {
				if e1.StartNID == nID {
					nID = e1.EndNID
				} else {
					strand = !strand
					nID = e1.StartNID
				}
			} else {
				if e1.EndNID == nID {
					nID = e1.StartNID
				} else {
					strand = !strand
					nID = e1.EndNID
				}
			}
		}

		var es []byte
		if strand == constructdbg.PLUS {
			es = e1.Utg.Ks[kmerlen-1:]
		} else {
			es = constructdbg.GetReverseCompByteArr(e1.Utg.Ks[:len(e1.Utg.Ks)-(kmerlen-1)])
		}
		seq = append(seq, es...)
		e0 = e1
	}

	e1 = edgesArr[m.Path[len(m.Path)-1]]
	if e1.StartNID == e1.EndNID {
		if strand == constructdbg.PLUS {
			if reflect.DeepEqual(e0.Utg.Ks[len(e0.Utg.Ks)-kmerlen+1:], e1.Utg.Ks[:kmerlen-1]) {

			} else {
				strand = !strand
			}
		} else {
			if reflect.DeepEqual(e0.Utg.Ks[:kmerlen-1], e1.Utg.Ks[len(e1.Utg.Ks)-kmerlen+1:]) {

			} else {
				strand = !strand
			}
		}
	} else {
		if e0.EndNID == nID {
			if e1.StartNID == nID {
			} else {
				strand = !strand
			}
		} else {
			if e1.EndNID == nID {
			} else {
				strand = !strand
			}
		}
	}
	//fmt.Printf("[GetPathSeq] len(e.Utg.Ks): %v, e: %v\n", len(e1.Utg.Ks), e1)
	var es []byte
	if strand == constructdbg.PLUS {
		es = e1.Utg.Ks[kmerlen-1 : m.EndP]
	} else {
		es = constructdbg.GetReverseCompByteArr(e1.Utg.Ks[m.EndP+1 : len(e1.Utg.Ks)-(kmerlen-1)])
	}
	seq = append(seq, es...) */

	return
}

func WriteChanPairReads(riArr [2]ReadMapInfo, edgesArr []constructdbg.DBGEdge, kmerlen int, wc chan<- constructcf.ReadInfo) {
	for x := 0; x < 2; x++ {
		var mr constructcf.ReadInfo
		//mr.Seq = GetPathSeq(riArr[x], edgesArr, kmerlen)
		mr.Seq = riArr[x].Seq
		mr.ID = riArr[x].ID
		mr.Anotition = fmt.Sprintf("%v:", riArr[x].StartP)
		for _, eID := range riArr[x].Path {
			mr.Anotition += fmt.Sprintf("%v:", eID)
		}
		mr.Anotition += fmt.Sprintf("%v", riArr[x].EndP)
		wc <- mr
	}
}

func GetReverseBoolArr(arr []bool) []bool {
	revArr := make([]bool, len(arr))
	for i, t := range arr {
		revArr[len(arr)-1-i] = t
	}

	return revArr
}

{
	if l0 == 1 && l1 == 1 && riArr[0].Path[0] == riArr[1].Path[0] {

		m.StartP = riArr[0].StartP
		m.Path = append(m.Path, riArr[0].Path[0])
		m.Strands = append(m.Strands, riArr[0].Strands[0])
		m.EndP = riArr[1].StartP
	} else {

	}
	// check if path has cycle
	e := edgesArr[riArr[0].Path[0]]
	merged := true
	if e.StartNID == e.EndNID {
		if len(e.Utg.Ks) > MaxPairLen-2*SD {
			if riArr[0].StartP < riArr[0].EndP {
				if riArr[1].StartP > riArr[1].EndP && riArr[1].EndP > riArr[0].StartP {

				} else {
					merged = false
				}
			} else {
				if riArr[1].StartP < riArr[1].EndP && riArr[0].EndP > riArr[1].StartP {

				} else {
					merged = false
				}
			}
		} else {
			merged = false
		}
	}
	/*if merged == false {
		WriteChanPairReads(riArr, edgesArr, cf.Kmerlen, wc)
		continue
	}*/


	// check if path has single or two edges cycle
	if l1 > l0 {
		riArr[0], riArr[1] = riArr[1], riArr[0]
		l0, l1 = l1, l0
	}
	if l0 > 1 {
		merged := true
		e1, e2 := edgesArr[riArr[0].Path[l0-2]], edgesArr[riArr[0].Path[l0-1]]
		if e1.ID == e2.ID {
			merged = false
		} else if constructdbg.IsBubble(e1.ID, e2.ID, edgesArr) {
			merged = false
		}

		if merged == false {
			WriteChanPairReads(riArr, edgesArr, cf.Kmerlen, wc)
			continue
		}
	}
	extLen := cf.Kmerlen - (len(pairRI[0].Seq) + len(pairRI[1].Seq) - MaxPairLen) + 50
	var nID constructdbg.DBG_MAX_INT
	pl := len(riArr[0].Path)
	if pl > 1 {
		e2 := edgesArr[riArr[0].Path[pl-1]]
		for _, eID := range riArr[0].Path[1 : len(riArr[0].Path)-1] {
			if edgesArr[eID].StartNID == nID {
				nID = edgesArr[eID].EndNID
			} else {
				nID = edgesArr[eID].StartNID
			}
		}
		//nID = constructdbg.GetInterNodeID(e1, e2)
		if nID == e2.StartNID {
			nID = e2.EndNID
			extLen -= (len(e2.Utg.Ks) - riArr[0].EndP)
		} else {
			nID = e2.StartNID
			extLen -= (riArr[0].EndP)
		}
	} else {
		e := edgesArr[riArr[0].Path[0]]
		if riArr[0].StartP < riArr[0].EndP {
			nID = e.EndNID
			extLen -= (len(e.Utg.Ks) - riArr[0].EndP)
		} else {
			nID = e.StartNID
			extLen -= (riArr[0].EndP)
		}
	}
	var rev1 ReadMapInfo
	rev1.Path = constructdbg.GetReverseDBG_MAX_INTArr(riArr[1].Path)
	rev1.StartP, rev1.EndP = riArr[1].EndP, riArr[1].StartP
	var linkPathArr [][]constructdbg.DBG_MAX_INT
	if extLen <= 0 {
		min := constructdbg.Min(len(riArr[0].Path), len(rev1.Path))
		for j := len(riArr[0].Path) - min; j < len(riArr[0].Path); j++ {
			if reflect.DeepEqual(riArr[0].Path[j:], rev1.Path[:len(riArr[0].Path)-j]) {
				var mergeP []constructdbg.DBG_MAX_INT
				mergeP = append(mergeP, riArr[0].Path...)
				mergeP = append(mergeP, rev1.Path[len(riArr[0].Path)-j:]...)
				linkPathArr = append(linkPathArr, mergeP)
				break
			}
		}
		//log.Fatalf("[paraMapNGSAndMerge] extLen: %v < 0, riArr: %v\n", extLen, riArr)
	} else {
		epArr := GetExtendPathArr(riArr[0].Path, nID, edgesArr, nodesArr, cf.Kmerlen, extLen)
		for _, ep := range epArr {
			idx := constructdbg.IndexEID(rev1.Path, ep[len(ep)-1])
			if len(ep)-1 < idx || idx < 0 {
				continue
			}
			//fmt.Printf("[paraMapNGSAndMerge]ep: %v,  riArr[1]: %v\n", ep, rev1)
			if reflect.DeepEqual(ep[len(ep)-1-idx:], rev1.Path[:idx+1]) {
				ep = append(ep, rev1.Path[idx+1:]...)
				linkPathArr = append(linkPathArr, ep)
			}
		}
	}
	if len(linkPathArr) == 1 {
		m.StartP = riArr[0].StartP
		m.Path = append(m.Path, linkPathArr[0]...)
		m.Strands = append(m.Strands, riArr[0].Strands...)
		m.EndP = rev1.EndP
	}
}

func paraMapNGSAndMerge(cs <-chan [2]constructcf.ReadInfo, wc chan<- constructcf.ReadInfo, nodesArr []constructdbg.DBGNode, edgesArr []constructdbg.DBGEdge, cf constructdbg.CuckooFilter, winSize, MaxPairLen, SD int) {
	var notFoundSeedNum, notPerfectNum, allNum, notMergeNum int
	for {
		var mR constructcf.ReadInfo
		pairRI, ok := <-cs
		if !ok {
			wc <- mR
			break
		}

		allNum++
		var riArr [2]ReadMapInfo
		for j := 0; j < 2; j++ {
			// found kmer seed position in the DBG edges
			//fmt.Printf("[paraMapNGSAndMerge] pairRI[%v]: %v\n", j, pairRI[j])
			dbgK, pos, strand := constructdbg.LocateSeedKmerCF(cf, pairRI[j], winSize, edgesArr)
			if dbgK.GetCount() == 0 { // not found in the cuckoofilter
				//fmt.Printf("[paraMapNGSAndMerge] read ID: %v not found seed!!!\n", pairRI[j].ID)
				notFoundSeedNum++
				break
			}
			// extend seed map to the edges
			// map the start partition of read sequence
			//fmt.Printf("[paraMapNGSAndMerge] dbgK: %v, pos: %v, strand: %v\n", dbgK, pos, strand)
			errorNum1, aB := constructdbg.MappingReadToEdgesBackWard(dbgK, pairRI[j], int(pos), strand, edgesArr, nodesArr, cf.Kmerlen, true)
			//func MappingReadToEdgesBackWard(dk DBGKmer, ri ReadInfo, rpos int, rstrand bool, edgesArr []DBGEdge, nodesArr []DBGNode, kmerlen int, correct bool) (errorNum int, ai AlignInfo) {
			aB.Paths = constructdbg.ReverseDBG_MAX_INTArr(aB.Paths)
			aB.Strands = GetReverseBoolArr(aB.Strands)
			// map the end partition of read sequence
			errorNum2, aF := constructdbg.MappingReadToEdgesForWard(dbgK, pairRI[j], int(pos)+cf.Kmerlen, strand, edgesArr, nodesArr, cf.Kmerlen, true)
			aB.ID, aF.ID = pairRI[j].ID, pairRI[j].ID
			//fmt.Printf("[paraMapNGSAndMerge] aB: %v\n\taF: %v\n", aB, aF)
			if aB.EndPos < -1 || aF.EndPos < -1 {
				fmt.Printf("[paraMapNGSAndMerge] not found end boundary\n")
				break
			}
			if errorNum1+errorNum2 > len(pairRI[j].Seq)*3/100 {
				fmt.Printf("[paraMapNGSAndMerge] errorNum: %v > %v\n", errorNum1+errorNum2, len(pairRI[j].Seq)*3/100)
				break
			}
			if len(aB.Paths) > 0 && len(aF.Paths) > 0 && aB.Paths[len(aB.Paths)-1] == aF.Paths[0] {
				riArr[j].PathSeqArr = append(aB.Paths, aF.Paths[1:]...)
				riArr[j].Strands = append(aB.Strands, aF.Strands[1:]...)
				riArr[j].Seq = append(aB.Seq, aF.Seq[cf.Kmerlen:]...)
				riArr[j].ID = pairRI[j].ID
				riArr[j].StartP = aB.EndPos
				riArr[j].EndP = aF.EndPos
			} else {
				log.Fatalf("[paraMapNGSAndMerge] aB: %v and aF: %v not consis\n", aB, aF)
			}
			//fmt.Printf("[paraMapNGSAndMerge] riArr[%v]: %v\n", j, riArr[j])
		}

		l0, l1 := len(riArr[0].Path), len(riArr[1].Path)
		if l0 < 1 ||  l1 < 1 {
			notPerfectNum++
			continue
		}

		// find link Path
		//var ri constructcf.ReadInfo
		//ri.ID = pairRI[0].ID
		// get link path sequence and pass to wc
		var m ReadMapInfo
		m.ID = riArr[0].ID
		m.PathSeqArr = MergePathSeqArr(riArr[0].PathSeqArr, riArr[1].PathSeqArr)
		if len(m.PathSeqArr) < 1 {
			notMergeNum++
			continue
		}
		m.Seq = GetPathSeq(m.PathSeqArr, edgesArr, nodesArr, cf.Kmerlen)
		m.Anotition = fmt.Sprintf("Path-- ")
		for _, pathSeq := range m.PathSeqArr {
			m.Anotition += fmt.Sprintf("%v:", pathSeq.ID)
		}
		wc <- m
	}
	fmt.Printf("[paraMapNGSAndMerge] not found seed read pair number is : %v,allNum: %v,  percent: %v\n", notFoundSeedNum,allNum, float32(notFoundSeedNum)/float32(allNum))
	fmt.Printf("[paraMapNGSAndMerge] too more error mapping read pair number is : %v, notMergeNum: %v\n", notPerfectNum-notFoundSeedNum, notMergeNum)
}

func writeCorrectReads(browfn string, wc <-chan constructcf.ReadInfo, numCPU int) (readNum int) {
	fp, err := os.Create(browfn)
	if err != nil {
		log.Fatalf("[writeCorrectReads] failed to create file: %s, err: %v\n", browfn, err)
	}
	defer fp.Close()
	cbrofp := cbrotli.NewWriter(fp, cbrotli.WriterOptions{Quality: 1})
	//cbrofp := cbrotli.NewWriter(fp, cbrotli.WriterOptions{Quality: 1})
	defer cbrofp.Close()
	//buffp := bufio.NewWriter(cbrofp)
	buffp := bufio.NewWriterSize(cbrofp, 1<<25) // 1<<24 == 2**24
	//gzfp := gzip.NewWriter(fp)
	//defer gzfp.Close()
	var finishNum int
	for {
		ri := <-wc
		if len(ri.Seq) == 0 {
			finishNum++
			if finishNum == numCPU {
				break
			} else {
				continue
			}
		}
		//seq := constructdbg.Transform2Letters(ri.Seq).String()
		seq := constructdbg.Transform2Char(ri.Seq)
		s := fmt.Sprintf(">%v\tpath:%s\n%s\n", ri.ID, ri.Anotition, string(seq))
		//cbrofp.Write([]byte(s))
		buffp.WriteString(s)
		readNum++
	}
	if err := buffp.Flush(); err != nil {
		log.Fatalf("[writeCorrectReads] failed to flush file: %s, err: %v\n", browfn, err)
	}

	if err := cbrofp.Flush(); err != nil {
		log.Fatalf("[writeCorrectReads] failed to flush file: %s, err: %v\n", browfn, err)
	}

	return
}

func paraProcessReadsFile(fn1, fn2 string, concurrentNum int, nodesArr []constructdbg.DBGNode, edgesArr []constructdbg.DBGEdge, cf constructdbg.CuckooFilter, opt Options, MaxPairLen, InsertSD int, processT chan int) {
	idx1 := strings.LastIndex(fn1, "1")
	idx2 := strings.LastIndex(fn2, "2")
	if !(idx1 > 0 && idx2 > 0 && fn1[:idx1] == fn2[:idx2]) {
		log.Fatalf("[paraProcessReadsFile] fn1: %v, fn2: %v, must used suffix *[1|2].[fasta|fa|fq|fastq].br\n", fn1, fn2)
	}

	bufSize := 60000
	cs := make(chan [2]constructcf.ReadInfo, bufSize)
	wc := make(chan constructcf.ReadInfo, bufSize)
	go LoadNGSReads(fn1, fn2, cs, opt.Kmer)
	for j := 0; j < concurrentNum; j++ {
		go paraMapNGSAndMerge(cs, wc, nodesArr, edgesArr, cf, opt.WinSize, MaxPairLen, InsertSD)
	}
	// write function
	brwfn := fn1[:idx1] + ".Correct.fa.br"
	writeNum := writeCorrectReads(brwfn, wc, concurrentNum)
	fmt.Printf("[paraProcessReadsFile] write correct reads num: %d to file: %s\n", writeNum, brwfn)
	processT <- 1
}

// MappingNGSAndCorrect() function parallel Map NGS reads to the DBG edges, then merge path output long single read
func MappingNGSAndMerge(opt Options, nodesArr []constructdbg.DBGNode, edgesArr []constructdbg.DBGEdge) {
	// construct cuckoofilter of DBG sample
	// use all edge seq, so set MaxNGSReadLen to MaxInt
	cfSize := constructdbg.GetCuckoofilterDBGSampleSize(edgesArr, int64(opt.WinSize), int64(math.MaxInt32), int64(opt.Kmer))
	fmt.Printf("[MappingNGSAndCorrect] cfSize: %v\n", cfSize)
	cf := constructdbg.MakeCuckooFilter(uint64(cfSize*5), opt.Kmer)
	fmt.Printf("[MappingNGSAndCorrect] cf.numItems: %v\n", cf.NumItems)
	count := constructdbg.ConstructCFDBGMinimizers(cf, edgesArr, opt.WinSize, math.MaxInt32)
	fmt.Printf("[MappingNGSAndCorrect]construct Smaple of DBG edges cuckoofilter number is : %v\n", count)

	cfgInfo, err := constructcf.ParseCfg(opt.CfgFn, opt.Correct)
	if err != nil {
		log.Fatal("[constructcf.ParseCfg] found err")
	}
	fmt.Printf("[MappingNGSAndCorrect] cfgInfo: %v\n", cfgInfo)

	runtime.GOMAXPROCS(opt.NumCPU + 2)

	concurrentNum := 3
	totalNumT := opt.NumCPU/(concurrentNum+1) + 1
	processT := make(chan int, totalNumT)
	for i := 0; i < totalNumT; i++ {
		processT <- 1
	}
	for _, lib := range cfgInfo.Libs {
		if lib.AsmFlag != constructcf.AllState && lib.SeqProfile != 1 {
			continue
		}
		MaxPairLen := lib.InsertSize + lib.InsertSD

		for i := 0; i < len(lib.FnName)-1; i += 2 {
			<-processT
			go paraProcessReadsFile(lib.FnName[i], lib.FnName[i+1], concurrentNum, nodesArr, edgesArr, cf, opt, MaxPairLen, lib.InsertSD, processT)
		}
	}

	for i := 0; i < totalNumT; i++ {
		<-processT
	}
}

func Correct(c cli.Command) {
	gOpt, suc := utils.CheckGlobalArgs(c.Parent())
	if suc == false {
		log.Fatalf("[Correct] check global Arguments error, opt: %v\n", gOpt)
	}
	opt := Options{gOpt, 0, 0, 0, 0, true}
	tmp, suc := checkArgs(c)
	if suc == false {
		log.Fatalf("[Correct] check Arguments error, opt: %v\n", tmp)
	}
	opt.CFSize = tmp.CFSize
	opt.MaxNGSReadLen = tmp.MaxNGSReadLen
	opt.TipMaxLen = tmp.TipMaxLen
	opt.WinSize = tmp.WinSize
	opt.Correct = tmp.Correct
	cfgInfo, err := constructcf.ParseCfg(opt.CfgFn, opt.Correct)
	if err != nil {
		log.Fatalf("[Correct] ParseCfg 'C': %v err :%v\n", opt.CfgFn, err)
	}
	fmt.Printf("[Correct] opt: %v\n\tcfgInfo: %v\n", opt, cfgInfo)

	// construct cuckoofilter and construct DBG
	constructcf.CCF(c)
	constructdbg.CDBG(c)

	profileFn := opt.Prefix + ".preprocess.prof"
	cpuprofilefp, err := os.Create(profileFn)
	if err != nil {
		log.Fatalf("[Correct] open cpuprofile file: %v failed\n", profileFn)
	}
	pprof.StartCPUProfile(cpuprofilefp)
	defer pprof.StopCPUProfile()
	// smfy DBG
	// read nodes file and transform to array mode for more quckly access
	nodesfn := opt.Prefix + ".nodes.mmap"
	nodeMap := constructdbg.NodeMapMmapReader(nodesfn)
	DBGStatfn := opt.Prefix + ".DBG.stat"
	nodesSize, edgesSize := constructdbg.DBGStatReader(DBGStatfn)
	//nodesSize := len(nodeMap)
	fmt.Printf("[Correct] len(nodeMap): %v, length of edge array: %v\n", nodesSize, edgesSize)
	// read edges file
	edgesfn := opt.Prefix + ".edges.fq"
	edgesArr := constructdbg.ReadEdgesFromFile(edgesfn, edgesSize)

	nodesArr := make([]constructdbg.DBGNode, nodesSize)
	constructdbg.NodeMap2NodeArr(nodeMap, nodesArr)
	nodeMap = nil // nodeMap any more used
	var copt constructdbg.Options
	copt.CfgFn, copt.Kmer, copt.MaxNGSReadLen, copt.NumCPU, copt.Prefix, copt.TipMaxLen, copt.WinSize = opt.CfgFn, opt.Kmer, opt.MaxNGSReadLen, opt.NumCPU, opt.Prefix, opt.TipMaxLen, opt.WinSize
	gfn := opt.Prefix + ".beforeSmfyDBG.dot"
	constructdbg.GraphvizDBGArr(nodesArr, edgesArr, gfn)
	constructdbg.SmfyDBG(nodesArr, edgesArr, copt)
	constructdbg.CheckDBGSelfCycle(nodesArr, edgesArr)

	constructdbg.PrintTmpDBG(nodesArr, edgesArr, opt.Prefix)
	gfn1 := opt.Prefix + ".afterSmfyDBG.dot"
	constructdbg.GraphvizDBGArr(nodesArr, edgesArr, gfn1)

	// Mapping NGS to DBG and Correct
	MappingNGSAndMerge(opt, nodesArr, edgesArr)
}
