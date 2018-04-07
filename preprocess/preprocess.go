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

func LoadNGSReads(brfn1, brfn2 string, cs chan<- [2]constructdbg.ReadInfo, kmerlen int) {
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
	var blockLineNum1, blockLineNum2 int
	if format1 == "fa" {
		blockLineNum1 = 2
	} else {
		blockLineNum1 = 4
	}
	if format2 == "fa" {
		blockLineNum2 = 2
	} else {
		blockLineNum2 = 4
	}
	var count int
	var EOF error
	for EOF != io.EOF {
		var b1, b2 [][]byte
		b1 = make([][]byte, blockLineNum1)
		b2 = make([][]byte, blockLineNum2)
		var i1, i2 int
		var error1, error2 error
		for ; i1 < blockLineNum1; i1++ {
			b1[i1], error1 = buffp1.ReadBytes('\n')
			if error1 != nil {
				break
			}
		}
		for ; i2 < blockLineNum2; i2++ {
			b2[i2], error2 = buffp2.ReadBytes('\n')
			if error2 != nil {
				break
			}
		}
		//fmt.Printf("[paraLoadNGSReads] b1: %v\n\tb2: %v\n", b1, b2)
		if error1 != nil || error2 != nil {
			if error1 == io.EOF || error2 == io.EOF {
				if !(error1 == io.EOF && error2 == io.EOF) {
					log.Fatalf("[LoadNGSReads] file : %v not consis with file : %v\n", brfn1, brfn2)
				}
				if i1 == 0 && i2 == 0 {
					break
				}
				if i1 != blockLineNum1-1 || i2 != blockLineNum2-1 {
					log.Fatalf("[LoadNGSReads] file1: %v, file2: %v,i1: %v, i2: %v, found not unbroken record\n", brfn1, brfn2, i1, i2)
				}
				EOF = io.EOF
			} else {
				if error1 != nil {
					log.Fatalf("[LoadNGSReads] file : %v encounter err: %v\n", brfn1, error1)
				}
				if error2 != nil {
					log.Fatalf("[LoadNGSReads] file : %v encounter err: %v\n", brfn2, error2)
				}
			}
		}
		id1, err := strconv.Atoi(string(b1[0][1 : len(b1[0])-1]))
		if err != nil {
			log.Fatalf("[LoadNGSReads] load fn: '%v' file, read ID: %v not digits, please convert to digits...\n", brfn1, b1[0][:len(b1[0])-1])
		}
		id2, err := strconv.Atoi(string(b2[0][1 : len(b2[0])-1]))
		if err != nil {
			log.Fatalf("[LoadNGSReads] load fn: '%v' file, read ID: %v not digits, please convert to digits...\n", brfn2, b2[0][:len(b2[0])-1])
		}
		if id1 != id2 {
			log.Fatalf("[LoadNGSReads] read1 ID : %v != read2 ID: %v\n", id1, id2)
		}
		var pairRI [2]constructdbg.ReadInfo
		pairRI[0].ID = int64(id1)
		pairRI[1].ID = int64(id2)
		b1[1], b2[1] = b1[1][:len(b1[1])-1], b2[1][:len(b2[1])-1]
		if len(b1[1]) < kmerlen+20 || len(b2[1]) < kmerlen+20 {
			continue
		}
		pairRI[0].Seq = constructdbg.Transform2BntByte(b1[1])
		pairRI[1].Seq = constructdbg.Transform2BntByte(b2[1])
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
	Path    []constructdbg.DBG_MAX_INT
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

func GetPathSeq(m ReadMapInfo, edgesArr []constructdbg.DBGEdge, kmerlen int) (seq []byte) {
	//fmt.Printf("[GetPathSeq] m: %v\n", m)
	if len(m.Path) == 1 {
		if m.StartP > m.EndP {
			seq = constructdbg.GetReverseCompByteArr(edgesArr[m.Path[0]].Utg.Ks[m.EndP+1 : m.StartP])
		} else {
			seq = edgesArr[m.Path[0]].Utg.Ks[m.StartP+1 : m.EndP]
		}
		return seq
	}
	strand := m.Strands[0]
	e0, e1 := edgesArr[m.Path[0]], edgesArr[m.Path[1]]
	/*if m.NID == 0 {
		log.Fatalf("[GetPathSeq] m.NID: %v not set\n", m.NID)
	}*/
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
	seq = append(seq, es...)

	return
}

func WriteChanPairReads(riArr [2]ReadMapInfo, edgesArr []constructdbg.DBGEdge, kmerlen int, wc chan<- constructdbg.ReadInfo) {
	for x := 0; x < 2; x++ {
		var mr constructdbg.ReadInfo
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

func paraMapNGSAndMerge(cs <-chan [2]constructdbg.ReadInfo, wc chan<- constructdbg.ReadInfo, nodesArr []constructdbg.DBGNode, edgesArr []constructdbg.DBGEdge, cf constructdbg.CuckooFilter, winSize, MaxPairLen, SD int) {
	var notFoundSeedNum, notPerfectNum, allNum int
	for {
		var mR constructdbg.ReadInfo
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
				riArr[j].Path = append(aB.Paths, aF.Paths[1:]...)
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

		if len(riArr[0].Path) == 0 || len(riArr[1].Path) == 0 {
			notPerfectNum++
			continue
		}
		l0, l1 := len(riArr[0].Path), len(riArr[1].Path)
		// find link Path
		var ri constructdbg.ReadInfo
		ri.ID = pairRI[0].ID
		var m ReadMapInfo
		if l0 == 1 && l1 == 1 && riArr[0].Path[0] == riArr[1].Path[0] {
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
			if merged == false {
				WriteChanPairReads(riArr, edgesArr, cf.Kmerlen, wc)
				continue
			}
			m.StartP = riArr[0].StartP
			m.Path = append(m.Path, riArr[0].Path[0])
			m.Strands = append(m.Strands, riArr[0].Strands[0])
			m.EndP = riArr[1].StartP
		} else {
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

		// get link path sequence and pass to wc
		if len(m.Path) > 0 {
			var mr constructdbg.ReadInfo
			mr.Seq = GetPathSeq(m, edgesArr, cf.Kmerlen)
			mr.ID = riArr[0].ID
			mr.Anotition = fmt.Sprintf("%v:", m.StartP)
			for _, eID := range m.Path {
				mr.Anotition += fmt.Sprintf("%v:", eID)
			}
			mr.Anotition += fmt.Sprintf("%v", m.EndP)
			wc <- mr
		} else {
			WriteChanPairReads(riArr, edgesArr, cf.Kmerlen, wc)
		}
	}
	fmt.Printf("[paraMapNGSAndMerge] not found seed read pair number is : %v, percent: %v\n", notFoundSeedNum, float32(notFoundSeedNum)/float32(allNum))
	fmt.Printf("[paraMapNGSAndMerge] too more error mapping read pair number is : %v\n", notPerfectNum-notFoundSeedNum)
}

func writeCorrectReads(browfn string, wc <-chan constructdbg.ReadInfo, numCPU int) (readNum int) {
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
	cs := make(chan [2]constructdbg.ReadInfo, bufSize)
	wc := make(chan constructdbg.ReadInfo, bufSize)
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

	cfgInfo, err := constructcf.ParseCfg(opt.CfgFn)
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
	opt := Options{gOpt, 0, 0, 0, 0}
	tmp, suc := checkArgs(c)
	if suc == false {
		log.Fatalf("[Correct] check Arguments error, opt: %v\n", tmp)
	}
	opt.CFSize = tmp.CFSize
	opt.MaxNGSReadLen = tmp.MaxNGSReadLen
	opt.TipMaxLen = tmp.TipMaxLen
	opt.WinSize = tmp.WinSize
	cfgInfo, err := constructcf.ParseCfg(opt.CfgFn)
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
	constructdbg.SmfyDBG(nodesArr, edgesArr, copt)

	gfn1 := opt.Prefix + ".afterSmfyDBG.dot"
	constructdbg.GraphvizDBGArr(nodesArr, edgesArr, gfn1)

	// Mapping NGS to DBG and Correct
	MappingNGSAndMerge(opt, nodesArr, edgesArr)
}
