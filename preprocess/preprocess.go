package preprocess

import (
	"bufio"
	"fmt"
	"log"
	"math"
	"os"
	"reflect"
	"runtime"
	"strconv"
	"strings"
	"time"
	"unsafe"
	//"github.com/google/brotli/cbrotli"
	"github.com/jwaldrip/odin/cli"
	"github.com/mudesheng/ga/bnt"
	//"github.com/mudesheng/ga/cbrotli"
	"github.com/google/brotli/go/cbrotli"
	"github.com/mudesheng/ga/constructcf"
	"github.com/mudesheng/ga/constructdbg"
	"github.com/mudesheng/ga/cuckoofilter"
	"github.com/mudesheng/ga/utils"
)

type Options struct {
	utils.ArgsOpt
	CFSize        int64
	TipMaxLen     int
	WinSize       int
	MaxNGSReadLen int
	Correct       bool
	Merge         bool
	MinKmerFreq   int
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
	merge := c.Flag("Merge").Get().(bool)
	opt.Merge = merge

	tmp, err = strconv.Atoi(c.Flag("WinSize").String())
	if err != nil {
		log.Fatalf("[checkArgs] argument 'MaxNGSReadLen': %v set error: %v\n", c.Flag("MaxNGSReadLen"), err)
	}
	if tmp < 1 || tmp > 20 {
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
	minKmerFreq := c.Flag("MinKmerFreq").Get().(int)
	if minKmerFreq < 1 || minKmerFreq > cuckoofilter.MAX_C {
		log.Fatalf("[checkArgs]the argument 'MinKmerFreq': %v must [%v ~ %v]\n", minKmerFreq, 1, cuckoofilter.MAX_C)
	}
	opt.MinKmerFreq = minKmerFreq
	suc = true
	return opt, suc
}

func LoadNGSReads(brfn string, cs chan<- constructcf.ReadInfo, kmerlen int) {
	bufSize := (1 << 20)
	//size := 8000
	fncs := make(chan []byte, 10)
	//recordChan := make(chan constructcf.ReadInfo, size)
	format := constructcf.GetReadsFileFormat(brfn)

	go constructcf.ReadBrFile(brfn, fncs, bufSize)
	count := constructcf.GetReadFileRecord(fncs, cs, format, bufSize, true)

	fmt.Printf("[LoadNGSReads] processed %d reads from files: %s\n", count, brfn)
	//time.Sleep(time.Second * 10)
}

func LoadNGSReadsPair(brfn1, brfn2 string, cs chan<- [2]constructcf.ReadInfo, kmerlen int, filterLowReadLen bool) {
	bufSize := (1 << 20)
	size := 8000
	fncs1 := make(chan []byte, 10)
	recordChan1 := make(chan constructcf.ReadInfo, size)
	format1 := constructcf.GetReadsFileFormat(brfn1)

	go constructcf.ReadBrFile(brfn1, fncs1, bufSize)
	go constructcf.GetReadFileRecord(fncs1, recordChan1, format1, bufSize, true)

	format2 := constructcf.GetReadsFileFormat(brfn2)
	fncs2 := make(chan []byte, 10)
	go constructcf.ReadBrFile(brfn2, fncs2, bufSize)
	recordChan2 := make(chan constructcf.ReadInfo, size)
	go constructcf.GetReadFileRecord(fncs2, recordChan2, format2, bufSize, true)

	var count int
	for {
		ri1, ok1 := <-recordChan1
		ri2, ok2 := <-recordChan2

		if ri1.ID != ri2.ID || ok1 != ok2 {
			log.Fatalf("[LoadNGSReadsPair] read1 ID : %v != read2 ID: %v\n", ri1.ID, ri2.ID)
		}
		if !ok1 {
			break
		}
		if ri1.ID == 0 {
			log.Fatalf("[LoadNGSReadsPair] read1 ID : %v\n", ri1.ID)
		}
		if filterLowReadLen {
			if len(ri1.Seq) < 140 || len(ri2.Seq) < 140 {
				continue
			}
		}
		var pairRI [2]constructcf.ReadInfo
		pairRI[0].ID = ri1.ID
		pairRI[1].ID = ri2.ID
		/*if len(ri1.Seq) < kmerlen || len(ri2.Seq) < kmerlen {
			continue
		}*/
		pairRI[0].Seq = ri1.Seq
		pairRI[1].Seq = ri2.Seq

		cs <- pairRI
		count++
	}
	fmt.Printf("[LoadNGSReadsPair] processed %d pair reads from pair files: %s , %s\n", count, brfn1, brfn2)
	close(cs)
	//time.Sleep(time.Second * 10)
}

/*func GetExtendPathArr(path []constructdbg.DBG_MAX_INT, nID constructdbg.DBG_MAX_INT, edgesArr []constructdbg.DBGEdge, nodesArr []constructdbg.DBGNode, kmerlen, extLen int) (epArr [][]constructdbg.DBG_MAX_INT) {
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
}*/

func ReversePathSeqArr(psArr []constructdbg.PathSeq) {
	la := len(psArr)
	for i := 0; i < len(psArr)/2; i++ {
		tmp := psArr[i]

		psArr[i].ID, psArr[i].Strand = psArr[la-i-1].ID, !psArr[la-i-1].Strand
		psArr[la-i-1].ID, psArr[la-i-1].Strand = tmp.ID, !tmp.Strand
	}
	if len(psArr)%2 == 1 {
		i := len(psArr) / 2
		psArr[i].Strand = !psArr[i].Strand
	}
}

func MappingReadToEdgesAndAddFreq(dk constructdbg.DBGKmer, ri constructcf.ReadInfo, rpos int, rstrand bool, edgesArr []constructdbg.DBGEdge, nodesArr []constructdbg.DBGNode, kmerlen int) (mappingNum int) {
	var strand bool
	if dk.Strand == rstrand {
		strand = constructdbg.PLUS
		dk.Pos += (int32(kmerlen) - 1)
	} else {
		strand = constructdbg.MINUS
	}
	rpos += kmerlen - 1
	mappingNum += kmerlen
	e := &edgesArr[dk.ID]
	for i := rpos; i < len(ri.Seq); i++ {
		edgesArr[dk.ID].Utg.Kq[dk.Pos]++
		if i == len(ri.Seq)-1 {
			break
		}
		if strand == constructdbg.PLUS {
			if dk.Pos < int32(e.GetSeqLen())-1 {
				eb := e.Utg.Ks[dk.Pos+1]
				if eb == ri.Seq[i+1] {
					dk.Pos++
					mappingNum++
				} else {
					break
				}
			} else {
				if e.EndNID == 0 {
					break
				}
				// found next edge
				eIDArr := constructdbg.GetNextEdgeIDArr(e, nodesArr, strand, constructdbg.FORWARD)
				if len(eIDArr) == 0 {
					break
				}
				ok := false
				for j := 0; j < len(eIDArr); j++ {
					ne := &edgesArr[eIDArr[j]]
					strand = constructdbg.GetNextEdgeStrand(e, ne, nodesArr, strand)
					e = ne
					dk.ID = ne.ID
					var eb byte
					if strand {
						dk.Pos = int32(kmerlen) - 1
						eb = e.Utg.Ks[dk.Pos]
					} else {
						dk.Pos = int32(e.GetSeqLen() - kmerlen)
						eb = bnt.BntRev[e.Utg.Ks[dk.Pos]]
					}
					if eb == ri.Seq[i+1] {
						//fmt.Printf("[MappingReadToEdgesAndAddFreq]el: %v i: %v, dk.ID: %v\n", e.GetSeqLen(), i, dk.ID)
						mappingNum++
						ok = true
						break
					}
				}
				if ok == false {
					break
				}
			}
		} else {
			if dk.Pos > 0 {
				eb := bnt.BntRev[e.Utg.Ks[dk.Pos-1]]
				if eb == ri.Seq[i+1] {
					dk.Pos--
					mappingNum++
				} else {
					break
				}
			} else {
				if e.StartNID == 0 {
					break
				}
				// found next edge
				eIDArr := constructdbg.GetNextEdgeIDArr(e, nodesArr, strand, constructdbg.FORWARD)
				if len(eIDArr) == 0 {
					break
				}
				ok := false
				for j := 0; j < len(eIDArr); j++ {
					ne := &edgesArr[eIDArr[j]]
					strand = constructdbg.GetNextEdgeStrand(e, ne, nodesArr, strand)
					e = ne
					dk.ID = ne.ID
					var eb byte
					if strand {
						dk.Pos = int32(kmerlen) - 1
						eb = e.Utg.Ks[dk.Pos]
					} else {
						dk.Pos = int32(e.GetSeqLen() - kmerlen)
						eb = bnt.BntRev[e.Utg.Ks[dk.Pos]]
					}
					if eb == ri.Seq[i+1] {
						//fmt.Printf("[MappingReadToEdgesAndAddFreq]el: %v i: %v, dk.ID: %v\n", e.GetSeqLen(), i, dk.ID)
						mappingNum++
						ok = true
						break
					}
				}
				if ok == false {
					break
				}
			}
		}
	}
	return
}

type PathSeqInfo struct {
	Seq     []byte
	Path    []constructdbg.DBG_MAX_INT
	Strand  bool
	NeedLen int
	EndPos  int
}

func GetMinSeqLen(psiArr []PathSeqInfo) int {
	min := math.MaxInt32
	for _, psi := range psiArr {
		if len(psi.Seq) < min {
			min = len(psi.Seq)
		}
	}
	return min
}

func MappingReadToEdges(dk constructdbg.DBGKmer, ri constructcf.ReadInfo, rpos int, rstrand bool, edgesArr []constructdbg.DBGEdge, nodesArr []constructdbg.DBGNode, kmerlen int) (int, int, constructdbg.ReadMapInfo, constructcf.ReadInfo) {
	var errorNum, mappingNum int
	var rmi constructdbg.ReadMapInfo
	rmi.PathSeqArr = make([]constructdbg.PathSeq, 0, 10)
	var strand bool
	var ps constructdbg.PathSeq
	if dk.Strand == rstrand {
		rmi.StartP = int(dk.Pos)
		dk.Pos += int32(kmerlen)
		strand = constructdbg.PLUS
	} else {
		rmi.StartP = int(dk.Pos) + kmerlen - 1
		dk.Pos--
		strand = constructdbg.MINUS
	}
	rpos += kmerlen
	mappingNum += kmerlen

	//fmt.Printf("[MappingReadToEdges]len(ri.Seq): %v\n", len(ri.Seq))
	for i := rpos; i < len(ri.Seq); {
		e := &edgesArr[dk.ID]
		b := len(ri.Seq)
		var j int
		//fmt.Printf("[MappingReadToEdges]i: %v,  strand: %v, dk.Pos: %v, dk.ID: %v, errorNum: %v\n", i, strand, dk.Pos, dk.ID, errorNum)
		if strand == constructdbg.PLUS {
			if len(e.Utg.Ks)-int(dk.Pos) < len(ri.Seq)-i {
				b = i + (len(e.Utg.Ks) - int(dk.Pos))
			}
			j = int(dk.Pos)
			for ; i < b && j < e.GetSeqLen(); i++ {
				if ri.Seq[i] != e.Utg.Ks[j] {
					/*if e.Utg.Ks[j] > 3 {
						fmt.Printf("[MappingReadToEdges] b: %v,i: %v,  j: %v, es:%v, eID:%d\n\te: %v", b, i, j, e.Utg.Ks[j], e.ID, e)
					}*/
					errorNum++
					ri.Seq[i] = e.Utg.Ks[j]
				} else {
					mappingNum++
				}
				j++
			}
		} else { // strand == MINUS
			if len(ri.Seq)-i > int(dk.Pos)+1 {
				b = i + (int(dk.Pos) + 1)
			}
			j = int(dk.Pos)
			for ; i < b && j >= 0; i++ {
				/*if e.Utg.Ks[j] > 3 {
					fmt.Printf("[MappingReadToEdges] b: %v,i: %v,  j: %v, es:%v, eID:%d\n\te: %v", b, i, j, e.Utg.Ks[j], e.ID, e)
				}*/
				if ri.Seq[i] != bnt.BntRev[e.Utg.Ks[j]] {
					errorNum++
					ri.Seq[i] = bnt.BntRev[e.Utg.Ks[j]]
				} else {
					mappingNum++
				}
				j--
			}
		}

		/*if errorNum > 0 {
			fmt.Printf("[MappingReadToEdges]not perfect end i: %v,edge ID: %v,len(e.Utg.Ks): %v,  dk.Pos: %v, pos: %v, b: %v\n", i, dk.ID, len(e.Utg.Ks), dk.Pos, rpos, b)
			break
		}*/

		ps.Strand = strand
		ps.ID = e.ID
		rmi.PathSeqArr = append(rmi.PathSeqArr, ps)
		//fmt.Printf("[MappingReadToEdges]ps: %v, strand: %v, i : %v\n", ps, strand, i)

		// check if move to the end of read
		if strand {
			rmi.EndP = j - 1
		} else {
			rmi.EndP = j + 1
		}
		//rmi.EndP = j
		if i >= len(ri.Seq) {
			break
		}
		// find next edge
		{
			rpos = i
			eIDArr := constructdbg.GetNextEdgeIDArr(e, nodesArr, strand, constructdbg.FORWARD)
			if len(eIDArr) == 0 {
				ri.Seq = ri.Seq[:i]
				//ri.Qual = ri.Qual[:i]
				break
			} else if len(eIDArr) == 1 {
				dk.ID = eIDArr[0]
				ne := &edgesArr[dk.ID]
				strand = constructdbg.GetNextEdgeStrand(e, ne, nodesArr, strand)
				if strand {
					dk.Pos = int32(kmerlen) - 1
				} else {
					dk.Pos = int32(ne.GetSeqLen() - kmerlen)
				}
				continue
			}
			ok := false
			for x := 0; x < len(eIDArr); x++ {
				dk.ID = eIDArr[x]
				ne := &edgesArr[dk.ID]
				std := constructdbg.GetNextEdgeStrand(e, ne, nodesArr, strand)
				var eb byte
				if std {
					dk.Pos = int32(kmerlen) - 1
					eb = ne.Utg.Ks[dk.Pos]
				} else {
					dk.Pos = int32(ne.GetSeqLen() - kmerlen)
					eb = bnt.BntRev[ne.Utg.Ks[dk.Pos]]
				}
				if eb == ri.Seq[i] {
					mappingNum++
					strand = std
					ok = true
					break
				}
			}
			if ok {
				continue
			}

			// need found all  paths seq for alignment
			{
				MaxPathLen := 4
				lps := rmi.PathSeqArr[len(rmi.PathSeqArr)-1]
				var psi PathSeqInfo
				psi.NeedLen = len(ri.Seq) - i
				psi.Path = append(psi.Path, lps.ID)
				psi.Strand = lps.Strand
				psiArr := make([]PathSeqInfo, 0, 20)
				stk := make([]PathSeqInfo, 0, 20)
				stk = append(stk, psi)
				for len(stk) > 0 {
					ps := stk[len(stk)-1]
					//fmt.Printf("[MappingReadToEdges]ps: %v\n", ps)
					stk = stk[:len(stk)-1]
					eID := ps.Path[len(ps.Path)-1]
					eIDArr := constructdbg.GetNextEdgeIDArr(&edgesArr[eID], nodesArr, ps.Strand, constructdbg.FORWARD)
					for x := 0; x < len(eIDArr); x++ {
						var nps PathSeqInfo
						id := eIDArr[x]
						ne := &edgesArr[id]
						std := constructdbg.GetNextEdgeStrand(&edgesArr[eID], ne, nodesArr, ps.Strand)
						nps.Strand = std
						nps.Path = make([]constructdbg.DBG_MAX_INT, len(ps.Path)+1)
						copy(nps.Path, ps.Path)
						nps.Path[len(nps.Path)-1] = id
						nl := len(ne.Utg.Ks) - (kmerlen - 1)
						if std {
							if nl >= ps.NeedLen {
								nl = ps.NeedLen
							}
							nps.EndPos = (kmerlen - 1) + nl
							nps.NeedLen = ps.NeedLen - nl
							//fmt.Printf("[MappingReadToEdges] nl: %v, el: %v\n", nl, ne.GetSeqLen())
							nps.Seq = make([]byte, len(ps.Seq)+nl)
							copy(nps.Seq, ps.Seq)
							copy(nps.Seq[len(ps.Seq):], ne.Utg.Ks[kmerlen-1:kmerlen-1+nl])
						} else {
							if nl >= ps.NeedLen {
								nl = ps.NeedLen
							}
							nps.EndPos = len(ne.Utg.Ks) - (kmerlen - 1) - nl - 1
							nps.NeedLen = ps.NeedLen - nl
							nps.Seq = make([]byte, len(ps.Seq)+nl)
							copy(nps.Seq, ps.Seq)
							copy(nps.Seq[len(ps.Seq):], constructdbg.GetReverseCompByteArr(ne.Utg.Ks[len(ne.Utg.Ks)-(kmerlen-1)-nl:len(ne.Utg.Ks)-(kmerlen-1)]))
							//fmt.Printf("[MappingReadToEdges] nl: %v, el: %v\n", nl, ne.GetSeqLen())
						}
						if nps.NeedLen == 0 || len(nps.Path) >= MaxPathLen {
							psiArr = append(psiArr, nps)
						} else {
							stk = append(stk, nps)
						}
					}
				}

				// porcess psiArr
				correct := false
				if len(psiArr) > 0 {
					minLen := GetMinSeqLen(psiArr)
					mchArr := make([]int, len(psiArr))
					maxScore := math.MinInt32
					maxCount := 0
					idx := -1
					//fmt.Printf("[MappingReadToEdges] ri.Seq[%v:]: %v\n", rpos, ri.Seq[rpos:])
					for x, psi := range psiArr {
						//fmt.Printf("[MappingReadToEdges] psiArr[%v]: %v\n", x, psi.Seq)
						for y := 0; y < minLen; y++ {
							if ri.Seq[rpos+y] == psi.Seq[y] {
								mchArr[x]++
							}
						}
						if mchArr[x] == 0 {
							continue
						}
						if mchArr[x] > maxScore {
							maxScore = mchArr[x]
							idx = x
							maxCount = 1
						} else if mchArr[x] == maxScore {
							maxCount++
						}
					}
					if maxCount == 1 {
						psi := psiArr[idx]
						en, mn := 0, 0
						for y := 0; y < len(psi.Seq); y++ {
							if ri.Seq[rpos+y] == psi.Seq[y] {
								mn++
							} else {
								en++
							}
						}
						if (len(psi.Seq) > 10 && en >= len(psi.Seq)/3) || (en >= len(psi.Seq)/2) {

						} else {
							errorNum += en
							mappingNum += mn
							//fmt.Printf("[MappingReadToEdges]i: %v, Len: %v, maxScore: %v, len(psi.Seq): %v\n", i, minLen, maxScore, len(psi.Seq))
							copy(ri.Seq[rpos:rpos+len(psi.Seq)], psi.Seq)
							rmi.EndP = psi.EndPos
							std := rmi.PathSeqArr[len(rmi.PathSeqArr)-1].Strand
							le := &edgesArr[psi.Path[0]]
							for x := 1; x < len(psi.Path); x++ {
								var ps constructdbg.PathSeq
								ps.ID = psi.Path[x]
								ne := &edgesArr[ps.ID]
								std = constructdbg.GetNextEdgeStrand(le, ne, nodesArr, std)
								ps.Strand = std
								rmi.PathSeqArr = append(rmi.PathSeqArr, ps)
								le = ne
							}
							rpos += len(psi.Seq)
							i = rpos
							dk.ID = psi.Path[len(psi.Path)-1]
							dk.Pos = int32(psi.EndPos)
							strand = psi.Strand
							correct = true
						}
					}
				}
				if !correct {
					ri.Seq = ri.Seq[:rpos]
					break
					//ri.Qual = ri.Qual[:rpos]
				}
			}
		}
	}
	return errorNum, mappingNum, rmi, ri
}

type Item struct {
	Strand    bool
	ExtendLen int
	EIDArr    []constructdbg.DBG_MAX_INT
}

func GetPairReadsLinkPathArr(e1, e2 constructdbg.DBGEdge, strand bool, maxAllowLen int, nodesArr []constructdbg.DBGNode, edgesArr []constructdbg.DBGEdge, kmerlen int) (extPathArr [][]constructdbg.DBG_MAX_INT) {
	var stack []Item
	var t Item
	t.EIDArr = append(t.EIDArr, e1.ID)
	t.Strand = strand
	stack = append(stack, t)

	// process stack
	for len(stack) > 0 {
		t := stack[len(stack)-1]
		stack = stack[:len(stack)-1]
		e := &edgesArr[t.EIDArr[len(t.EIDArr)-1]]
		std := t.Strand

		ea := constructdbg.GetNextEdgeIDArr(e, nodesArr, std, constructdbg.FORWARD)
		for _, id := range ea {
			var nt Item
			nt.EIDArr = make([]constructdbg.DBG_MAX_INT, len(t.EIDArr)+1)
			copy(nt.EIDArr, t.EIDArr)
			nt.EIDArr[len(nt.EIDArr)-1] = id
			ne := edgesArr[id]
			if ne.StartNID == ne.EndNID || ne.GetTwoEdgesCycleFlag() > 0 {
				extPathArr = nil
				return
			}
			if ne.ID == e2.ID {
				extPathArr = append(extPathArr, nt.EIDArr)
				continue
			}
			nt.Strand = std
			nt.ExtendLen = t.ExtendLen + (len(e.Utg.Ks) - (kmerlen - 1))
			if nt.ExtendLen >= maxAllowLen {
				continue
			}
			stack = append(stack, nt)
		}
	}

	return
}

func MergePathSeqArr(rmi1, rmi2 constructdbg.ReadMapInfo, edgesArr []constructdbg.DBGEdge, nodesArr []constructdbg.DBGNode, kmerlen int) (mergeRMI constructdbg.ReadMapInfo) {
	pathSeqArr1 := rmi1.PathSeqArr
	pathSeqArr2 := make([]constructdbg.PathSeq, len(rmi2.PathSeqArr))
	copy(pathSeqArr2, rmi2.PathSeqArr)
	lastT1 := pathSeqArr1[len(pathSeqArr1)-1]
	//e1, e2 := edgesArr[lastT1.ID], edgesArr[lastT2.ID]
	//strand1, strand2 := lastT1.Strand, lastT2.Strand
	mergeRMI.StartP = rmi1.StartP
	mergeRMI.EndP = rmi2.StartP
	// reverse pathSeqArr2
	ReversePathSeqArr(pathSeqArr2)
	//fmt.Printf("[MergePathSeqArr]rmi1: %v\n\t\trmi1: %v\n", rmi1, rmi2)
	//fmt.Printf("[MergePathSeqArr]psArr1: %v\n\t\tpsArr2: %v\n", pathSeqArr1, pathSeqArr2)
	var sharePathLen int
	for i := len(pathSeqArr2) - 1; i >= 0; i-- {
		id := pathSeqArr2[i].ID
		if lastT1.ID == id && lastT1.Strand == pathSeqArr2[i].Strand {
			if i <= len(pathSeqArr1)-1 {
				if reflect.DeepEqual(pathSeqArr1[len(pathSeqArr1)-1-i:], pathSeqArr2[:i+1]) {
					sharePathLen = i + 1
					break
				}
			}
		}
	}

	if sharePathLen < 1 {
		lastT2 := pathSeqArr2[0]
		e1, e2 := edgesArr[lastT1.ID], edgesArr[lastT2.ID]
		if e1.StartNID == e1.EndNID {
			return
		}
		//fmt.Printf("[MergePathSeqArr]lastT1: %v\tlastT2: %v\n", lastT1, lastT2)
		concat := false
		if lastT1.Strand == constructdbg.PLUS {
			if lastT2.Strand == constructdbg.PLUS {
				if e1.EndNID == e2.StartNID {
					concat = true
				}
			} else {
				if e1.EndNID == e2.EndNID {
					concat = true
				}
			}
		} else {
			if lastT2.Strand == constructdbg.PLUS {
				if e1.StartNID == e2.StartNID {
					concat = true
				}
			} else {
				if e1.StartNID == e2.EndNID {
					concat = true
				}
			}
		}

		if !concat {
			return
		}
		// maybe can find unique link path
		/*pairReadLen := len(pairRI[0].Seq) + len(pairRI[1].Seq)
		maxAllowLen := insertSize - pairReadLen + kmerlen + SD - flank1 - flank2
		if maxAllowLen <= 0 {
			return
		}
		extPathArr := GetPairReadsLinkPathArr(e1, e2, strand1, maxAllowLen, nodesArr, edgesArr, kmerlen)
		if len(extPathArr) == 1 {
			//fmt.Printf("[MergePathSeqArr]found link path: %v\n", extPathArr[0])
			// add forward read path
			mergeRMI.PathSeqArr = append(mergeRMI.PathSeqArr, pathSeqArr1...)
			// add link path
			strand := strand1
			for x := 1; x < len(extPathArr[0]); x++ {
				var t constructdbg.PathSeq
				id := extPathArr[0][x]
				e := edgesArr[extPathArr[0][x-1]]
				ne := edgesArr[id]
				t.ID = id
				strand = constructdbg.GetNextEdgeStrand(e, ne, nodesArr, strand)
				t.Strand = strand
				//fmt.Printf("[MergePathSeqArr]t: %v\n", t)
				mergeRMI.PathSeqArr = append(mergeRMI.PathSeqArr, t)
			}
			if strand != pathSeqArr2[0].Strand {
				mergeRMI.PathSeqArr = nil
				fmt.Printf("[MergePathSeqArr] deduce strand: %v not consis with BACKWARD read strand: %v\n", strand, pathSeqArr2[0].Strand)
				return
			}
			// add BACKWARD read path
			mergeRMI.PathSeqArr = append(mergeRMI.PathSeqArr, pathSeqArr2[1:]...)
		} */
	}

	mergeRMI.PathSeqArr = append(mergeRMI.PathSeqArr, pathSeqArr1...)
	mergeRMI.PathSeqArr = append(mergeRMI.PathSeqArr, pathSeqArr2[sharePathLen:]...)
	// check not a pair end reads
	if len(mergeRMI.PathSeqArr) == 1 {
		if pathSeqArr1[0].Strand == constructdbg.PLUS {
			if mergeRMI.StartP > mergeRMI.EndP {
				fmt.Fprintf(os.Stderr, "[MergePathSeqArr]Postion error, p1: %v, p2: %v, mergeRMI: %v\n", pathSeqArr1, pathSeqArr2, mergeRMI)
				mergeRMI.PathSeqArr = nil
			}
		} else {
			if mergeRMI.StartP < mergeRMI.EndP {
				fmt.Fprintf(os.Stderr, "[MergePathSeqArr]Postion error, p1: %v, p2: %v, mergeRMI: %v\n", pathSeqArr1, pathSeqArr2, mergeRMI)
				mergeRMI.PathSeqArr = nil
			}
		}
	}
	return
}

func GetPathSeq(m constructdbg.ReadMapInfo, edgesArr []constructdbg.DBGEdge, nodesArr []constructdbg.DBGNode, kmerlen int, riArr [2]constructdbg.ReadMapInfo) (seq []byte) {
	//fmt.Printf("[GetPathSeq] m: %v\n", m)
	strand := m.PathSeqArr[0].Strand
	if len(m.PathSeqArr) == 1 {
		ps := m.PathSeqArr[0]
		if strand == constructdbg.PLUS {
			seq = make([]byte, m.EndP-m.StartP)
			copy(seq, edgesArr[ps.ID].Utg.Ks[m.StartP:m.EndP])
		} else {
			/*if m.EndP < 0 || m.StartP > len(edgesArr[ps.ID].Utg.Ks) {
				fmt.Printf("[GetPathSeq] error set m: %v, len(e.Utg.Ks): %v\n", m, len(edgesArr[ps.ID].Utg.Ks))
				seq = nil
				return
			}*/
			//fmt.Printf("[GetPathSeq] m: %v, len(e.Utg.Ks): %v\n", m, len(edgesArr[ps.ID].Utg.Ks))
			seq = constructdbg.GetReverseCompByteArr(edgesArr[ps.ID].Utg.Ks[m.EndP:m.StartP])
		}
		return seq
	}

	seq = make([]byte, 0, 500)
	//var strand bool // true denote "+" DNA strand map, false denote  "-" DNA strand map
	//var coming bool // true denote EdgeOutcoming, false denote EdgeIncoming
	for i := 0; i < len(m.PathSeqArr); i++ {
		ps := m.PathSeqArr[i]
		e := &edgesArr[ps.ID]
		if strand != ps.Strand {
			fmt.Fprintf(os.Stderr, "[GetPathSeq]m.PathSeqArr[%d]: %v, strand: %v, Cycle:%v\n", i, m.PathSeqArr[i], strand, edgesArr[m.PathSeqArr[i].ID].StartNID == edgesArr[m.PathSeqArr[i].ID].EndNID)
		}
		if i == 0 {
			if strand == constructdbg.PLUS {
				seq = append(seq, e.Utg.Ks[m.StartP:]...)
			} else {
				seq = append(seq, constructdbg.GetReverseCompByteArr(e.Utg.Ks[:m.StartP])...)
			}
		} else if i == len(m.PathSeqArr)-1 {
			//fmt.Printf("[GetPathSeq]len(e.Utg.Ks): %v, strand: %v, m.EndP: %v\n", len(e.Utg.Ks), strand, m.EndP)
			if strand == constructdbg.PLUS {
				if m.EndP < kmerlen-1 {
					fmt.Fprintf(os.Stderr, "[GetPathSeq]m.EndP: %v set error, len(e.Utg.Ks): %v, strand: %v, new m.EndP: %v\n\tm: %v\n\triArr[0]: %v\n\triArr[1]: %v\n", m.EndP, len(e.Utg.Ks), strand, kmerlen-1, m, riArr[0], riArr[1])
					m.EndP = kmerlen - 1
				} else if m.EndP > len(e.Utg.Ks) {
					fmt.Fprintf(os.Stderr, "[GetPathSeq]m.EndP: %v set error, len(e.Utg.Ks): %v, strand: %v, new m.EndP: %v\n\tm: %v\n\triArr[0]: %v\n\triArr[1]: %v\n", m.EndP, len(e.Utg.Ks), strand, len(e.Utg.Ks), m, riArr[0], riArr[1])
					m.EndP = len(e.Utg.Ks)
				}
				seq = append(seq, e.Utg.Ks[kmerlen-1:m.EndP]...)
			} else {
				if m.EndP < 0 {
					fmt.Fprintf(os.Stderr, "[GetPathSeq]m.EndP: %v set error, len(e.Utg.Ks): %v, strand: %v, new m.EndP: %v\n\tm: %v\n\triArr[0]: %v\n\triArr[1]: %v\n", m.EndP, len(e.Utg.Ks), strand, 0, m, riArr[0], riArr[1])
					m.EndP = 0
				} else if m.EndP > len(e.Utg.Ks)-(kmerlen-1) {
					fmt.Fprintf(os.Stderr, "[GetPathSeq]m.EndP: %v set error, len(e.Utg.Ks): %v, strand: %v, new m.EndP: %v\n\tm: %v\n\triArr[0]: %v\n\triArr[1]: %v\n", m.EndP, len(e.Utg.Ks), strand, len(e.Utg.Ks)-(kmerlen-1), m, riArr[0], riArr[1])
					m.EndP = len(e.Utg.Ks) - (kmerlen - 1)
				}
				seq = append(seq, constructdbg.GetReverseCompByteArr(e.Utg.Ks[m.EndP:len(e.Utg.Ks)-(kmerlen-1)])...)
			}
		} else {
			if strand == constructdbg.PLUS {
				seq = append(seq, e.Utg.Ks[kmerlen-1:]...)
			} else {
				seq = append(seq, constructdbg.GetReverseCompByteArr(e.Utg.Ks[:len(e.Utg.Ks)-(kmerlen-1)])...)
			}
		}
		if i < len(m.PathSeqArr)-1 {
			ne := &edgesArr[m.PathSeqArr[i+1].ID]
			strand = constructdbg.GetNextEdgeStrand(e, ne, nodesArr, strand)
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

	return seq
}

/*func WriteChanPairReads(riArr [2]constructdbg.ReadMapInfo, edgesArr []constructdbg.DBGEdge, kmerlen int, wc chan<- constructcf.ReadInfo) {
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
} */

func GetReverseBoolArr(arr []bool) []bool {
	revArr := make([]bool, len(arr))
	for i, t := range arr {
		revArr[len(arr)-1-i] = t
	}

	return revArr
}

/*{
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
	//if merged == false {
	//	WriteChanPairReads(riArr, edgesArr, cf.Kmerlen, wc)
	//	continue
	//}


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
	var rev1 constructdbg.ReadMapInfo
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
} */

func CountNumInPathSeqArr(eID constructdbg.DBG_MAX_INT, arr []constructdbg.PathSeq) (count int) {
	for _, item := range arr {
		if eID == item.ID {
			count++
		}
	}
	return count
}

func paraMapNGSAndCorrect(cs <-chan constructcf.ReadInfo, wc chan<- constructcf.ReadInfo, nodesArr []constructdbg.DBGNode, edgesArr []constructdbg.DBGEdge, cf constructdbg.CuckooFilter, winSize, kmerlen int) {
	var notFoundSeedNum, shortMappingNum, highErrorNum, allNum int
	MaxStepNum := 10
	ka := make([]constructdbg.KmerInfo, MaxStepNum*winSize)
	rSeq := make([]byte, 300)
	for {
		ri, ok := <-cs
		if !ok {
			var t constructcf.ReadInfo
			wc <- t
			break
		}

		if len(ri.Seq) < 140 {
			continue
		}

		allNum++
		// found kmer seed position in the DBG edges
		dbgK, pos, strand := constructdbg.LocateSeedKmerCF(cf, ri, winSize, edgesArr, MaxStepNum, ka, rSeq)
		if dbgK.GetCount() == 0 { // not found in the cuckoofilter
			//fmt.Printf("[paraMapNGSAndMerge] read ID: %v not found seed!!!\n", pairRI[j].ID)
			notFoundSeedNum++
			wc <- ri
			continue
		}
		//fmt.Printf("[paraMapNGSAndMerge] pairRI[%v]: %v\n", j, pairRI[j])

		//fmt.Printf("[paraMapNGSAndMerge]readID: %v, dbgK: %v, pos: %v, strand: %v\n", pairRI[0].ID, dbgK, pos, strand)
		errorNum, mappingNum, _, correctRI := MappingReadToEdges(dbgK, ri, int(pos), strand, edgesArr, nodesArr, cf.Kmerlen)
		// extend seed map to the edges
		// map the start partition of read sequence
		//errorNum1, aB := constructdbg.MappingReadToEdgesBackWard(dbgK, pairRI[j], int(pos), strand, edgesArr, nodesArr, cf.Kmerlen, true)
		//func MappingReadToEdgesBackWard(dk DBGKmer, ri ReadInfo, rpos int, rstrand bool, edgesArr []DBGEdge, nodesArr []DBGNode, kmerlen int, correct bool) (errorNum int, ai AlignInfo) {
		//aB.Paths = constructdbg.ReverseDBG_MAX_INTArr(aB.Paths)
		//aB.Strands = GetReverseBoolArr(aB.Strands)
		// map the end partition of read sequence
		//errorNum2, aF := constructdbg.MappingReadToEdgesForWard(dbgK, pairRI[j], int(pos)+cf.Kmerlen, strand, edgesArr, nodesArr, cf.Kmerlen, true)
		//aB.ID, aF.ID = pairRI[j].ID, pairRI[j].ID
		//fmt.Printf("[paraMapNGSAndMerge] aB: %v\n\taF: %v\n", aB, aF)
		//if aB.EndPos < -1 || aF.EndPos < -1 {
		//	fmt.Printf("[paraMapNGSAndMerge] not found end boundary\n")
		//	break
		//}

		if mappingNum < (len(ri.Seq)-pos)*7/10 {
			//fmt.Printf("[paraMapNGSAndMerge] mappingNum: %v < %v\n", mappingNum, len(pairRI[j].Seq)*8/10)
			shortMappingNum++
			wc <- ri
			continue
		}

		if errorNum > (len(ri.Seq)-pos)*1/6 {
			//fmt.Printf("[paraMapNGSAndMerge] errorNum: %v > %v\n", errorNum[j], (len(pairRI[j].Seq)-pos)*8/100)
			highErrorNum++
			wc <- ri
			continue
		}

		var mR constructcf.ReadInfo
		mR.Seq = correctRI.Seq
		mR.ID = ri.ID
		//mR.Anotition = ""
		mR.Anotition = append(mR.Anotition, fmt.Sprintf("EN:%v", errorNum)...)
		//fmt.Printf("[paraMapNGSAndMerge]len(mR.Seq): %v,  mR.Seq: %v\n", len(mR.Seq), mR.Seq)
		wc <- mR
	}
	fmt.Printf("[paraMapNGSAndCorrect] not found seed read number is : %v,allNum: %v,  percent: %v\n", notFoundSeedNum, allNum, float32(notFoundSeedNum)/float32(allNum))
	fmt.Printf("[paraMapNGSAndCorrect] high error mapping read number is : %v  short mapping number is: %v\n", highErrorNum, shortMappingNum)
}

func paraMapNGSAndCorrectPair(cs <-chan [2]constructcf.ReadInfo, wc chan<- [2]constructcf.ReadInfo, nodesArr []constructdbg.DBGNode, edgesArr []constructdbg.DBGEdge, cf constructdbg.CuckooFilter, winSize, kmerlen int) {
	var notFoundSeedNum, shortMappingNum, highErrorNum, allNum int
	MaxStepNum := 10
	ka := make([]constructdbg.KmerInfo, MaxStepNum*winSize)
	rSeq := make([]byte, 300)
	for {
		pairRI, ok := <-cs
		if !ok {
			var t [2]constructcf.ReadInfo
			wc <- t
			break
		}

		if len(pairRI[0].Seq) < 140 || len(pairRI[1].Seq) < 140 {
			continue
		}

		var wr [2]constructcf.ReadInfo
		allNum++
		// read1
		{
			ok := true
			// found kmer seed position in the DBG edges
			dbgK, pos, strand := constructdbg.LocateSeedKmerCF(cf, pairRI[0], winSize, edgesArr, MaxStepNum, ka, rSeq)
			if dbgK.GetCount() == 0 { // not found in the cuckoofilter
				//fmt.Printf("[paraMapNGSAndMerge] read ID: %v not found seed!!!\n", pairRI[j].ID)
				notFoundSeedNum++
				wr[0] = pairRI[0]
				ok = false
			}
			//fmt.Printf("[paraMapNGSAndMerge] pairRI[%v]: %v\n", j, pairRI[j])

			//fmt.Printf("[paraMapNGSAndMerge]readID: %v, dbgK: %v, pos: %v, strand: %v\n", pairRI[0].ID, dbgK, pos, strand)
			if ok {
				errorNum, mappingNum, _, correctRI := MappingReadToEdges(dbgK, pairRI[0], int(pos), strand, edgesArr, nodesArr, cf.Kmerlen)
				// extend seed map to the edges
				// map the start partition of read sequence
				//errorNum1, aB := constructdbg.MappingReadToEdgesBackWard(dbgK, pairRI[j], int(pos), strand, edgesArr, nodesArr, cf.Kmerlen, true)
				//func MappingReadToEdgesBackWard(dk DBGKmer, ri ReadInfo, rpos int, rstrand bool, edgesArr []DBGEdge, nodesArr []DBGNode, kmerlen int, correct bool) (errorNum int, ai AlignInfo) {
				//aB.Paths = constructdbg.ReverseDBG_MAX_INTArr(aB.Paths)
				//aB.Strands = GetReverseBoolArr(aB.Strands)
				// map the end partition of read sequence
				//errorNum2, aF := constructdbg.MappingReadToEdgesForWard(dbgK, pairRI[j], int(pos)+cf.Kmerlen, strand, edgesArr, nodesArr, cf.Kmerlen, true)
				//aB.ID, aF.ID = pairRI[j].ID, pairRI[j].ID
				//fmt.Printf("[paraMapNGSAndMerge] aB: %v\n\taF: %v\n", aB, aF)
				//if aB.EndPos < -1 || aF.EndPos < -1 {
				//	fmt.Printf("[paraMapNGSAndMerge] not found end boundary\n")
				//	break
				//}

				if mappingNum < (len(pairRI[0].Seq)-pos)*7/10 {
					//fmt.Printf("[paraMapNGSAndMerge] mappingNum: %v < %v\n", mappingNum, len(pairRI[j].Seq)*8/10)
					shortMappingNum++
					wr[0] = pairRI[0]
					ok = false
				}

				if ok && errorNum > (len(pairRI[0].Seq)-pos)*1/6 {
					//fmt.Printf("[paraMapNGSAndMerge] errorNum: %v > %v\n", errorNum[j], (len(pairRI[j].Seq)-pos)*8/100)
					highErrorNum++
					wr[0] = pairRI[0]
					ok = false
				}

				if ok {
					var mR constructcf.ReadInfo
					mR.Seq = correctRI.Seq
					mR.ID = pairRI[0].ID
					//mR.Anotition = ""
					mR.Anotition = append(mR.Anotition, fmt.Sprintf("EN:%v", errorNum)...)
					//fmt.Printf("[paraMapNGSAndMerge]len(mR.Seq): %v,  mR.Seq: %v\n", len(mR.Seq), mR.Seq)
					wr[0] = mR
				}
			}
		}

		// read2
		{
			ok := true
			// found kmer seed position in the DBG edges
			dbgK, pos, strand := constructdbg.LocateSeedKmerCF(cf, pairRI[1], winSize, edgesArr, MaxStepNum, ka, rSeq)
			if dbgK.GetCount() == 0 { // not found in the cuckoofilter
				//fmt.Printf("[paraMapNGSAndMerge] read ID: %v not found seed!!!\n", pairRI[j].ID)
				notFoundSeedNum++
				wr[1] = pairRI[1]
				ok = false
			}
			//fmt.Printf("[paraMapNGSAndMerge] pairRI[%v]: %v\n", j, pairRI[j])

			//fmt.Printf("[paraMapNGSAndMerge]readID: %v, dbgK: %v, pos: %v, strand: %v\n", pairRI[0].ID, dbgK, pos, strand)
			if ok {
				errorNum, mappingNum, _, correctRI := MappingReadToEdges(dbgK, pairRI[1], int(pos), strand, edgesArr, nodesArr, cf.Kmerlen)
				// extend seed map to the edges
				// map the start partition of read sequence
				//errorNum1, aB := constructdbg.MappingReadToEdgesBackWard(dbgK, pairRI[j], int(pos), strand, edgesArr, nodesArr, cf.Kmerlen, true)
				//func MappingReadToEdgesBackWard(dk DBGKmer, ri ReadInfo, rpos int, rstrand bool, edgesArr []DBGEdge, nodesArr []DBGNode, kmerlen int, correct bool) (errorNum int, ai AlignInfo) {
				//aB.Paths = constructdbg.ReverseDBG_MAX_INTArr(aB.Paths)
				//aB.Strands = GetReverseBoolArr(aB.Strands)
				// map the end partition of read sequence
				//errorNum2, aF := constructdbg.MappingReadToEdgesForWard(dbgK, pairRI[j], int(pos)+cf.Kmerlen, strand, edgesArr, nodesArr, cf.Kmerlen, true)
				//aB.ID, aF.ID = pairRI[j].ID, pairRI[j].ID
				//fmt.Printf("[paraMapNGSAndMerge] aB: %v\n\taF: %v\n", aB, aF)
				//if aB.EndPos < -1 || aF.EndPos < -1 {
				//	fmt.Printf("[paraMapNGSAndMerge] not found end boundary\n")
				//	break
				//}

				if mappingNum < (len(pairRI[1].Seq)-pos)*7/10 {
					//fmt.Printf("[paraMapNGSAndMerge] mappingNum: %v < %v\n", mappingNum, len(pairRI[j].Seq)*8/10)
					shortMappingNum++
					wr[1] = pairRI[1]
					ok = false
				}

				if ok && errorNum > (len(pairRI[1].Seq)-pos)*1/6 {
					//fmt.Printf("[paraMapNGSAndMerge] errorNum: %v > %v\n", errorNum[j], (len(pairRI[j].Seq)-pos)*8/100)
					highErrorNum++
					wr[1] = pairRI[1]
					ok = false
				}

				if ok {
					var mR constructcf.ReadInfo
					mR.Seq = correctRI.Seq
					mR.ID = pairRI[1].ID
					//mR.Anotition = ""
					mR.Anotition = append(mR.Anotition, fmt.Sprintf("EN:%v", errorNum)...)
					//mR.Anotition = fmt.Sprintf("EN:%v", errorNum) // EN note errorNum
					//fmt.Printf("[paraMapNGSAndMerge]len(mR.Seq): %v,  mR.Seq: %v\n", len(mR.Seq), mR.Seq)
					wr[1] = mR
				}
			}
		}
		wc <- wr
	}
	fmt.Printf("[paraMapNGSAndCorrectPair] not found seed read number is : %v,allNum: %v,  percent: %v\n", notFoundSeedNum, allNum, float32(notFoundSeedNum)/float32(allNum))
	fmt.Printf("[paraMapNGSAndCorrectPair] high error mapping read number is : %v  short mapping number is: %v\n", highErrorNum, shortMappingNum)
}

func paraMapNGSAndMerge(cs <-chan [2]constructcf.ReadInfo, wc chan<- constructcf.ReadInfo, nodesArr []constructdbg.DBGNode, edgesArr []constructdbg.DBGEdge, cf constructdbg.CuckooFilter, winSize, kmerlen int) {
	var notFoundSeedNum, shortMappingNum, highErrorNum, allNum, singleNum, mergeNum, unMergeNum int
	MaxStepNum := 10
	ka := make([]constructdbg.KmerInfo, MaxStepNum*winSize)
	rSeq := make([]byte, 300)
	for {
		pairRI, ok := <-cs
		if !ok {
			var t constructcf.ReadInfo
			wc <- t
			break
		}

		if len(pairRI[0].Seq) < 140 || len(pairRI[1].Seq) < 140 {
			continue
		}

		allNum += 2
		var riArr [2]constructdbg.ReadMapInfo
		var errorNum [2]int
		needMerge := true
		for j := 0; j < 2; j++ {
			// found kmer seed position in the DBG edges
			dbgK, pos, strand := constructdbg.LocateSeedKmerCF(cf, pairRI[j], winSize, edgesArr, MaxStepNum, ka, rSeq)
			if dbgK.GetCount() == 0 { // not found in the cuckoofilter
				//fmt.Printf("[paraMapNGSAndMerge] read ID: %v not found seed!!!\n", pairRI[j].ID)
				notFoundSeedNum++
				needMerge = false
				continue
			}
			//fmt.Printf("[paraMapNGSAndMerge] pairRI[%v]: %v\n", j, pairRI[j])

			//fmt.Printf("[paraMapNGSAndMerge]readID: %v, dbgK: %v, pos: %v, strand: %v\n", pairRI[0].ID, dbgK, pos, strand)
			var mappingNum int
			errorNum[j], mappingNum, riArr[j], pairRI[j] = MappingReadToEdges(dbgK, pairRI[j], int(pos), strand, edgesArr, nodesArr, cf.Kmerlen)
			// extend seed map to the edges
			// map the start partition of read sequence
			//errorNum1, aB := constructdbg.MappingReadToEdgesBackWard(dbgK, pairRI[j], int(pos), strand, edgesArr, nodesArr, cf.Kmerlen, true)
			//func MappingReadToEdgesBackWard(dk DBGKmer, ri ReadInfo, rpos int, rstrand bool, edgesArr []DBGEdge, nodesArr []DBGNode, kmerlen int, correct bool) (errorNum int, ai AlignInfo) {
			//aB.Paths = constructdbg.ReverseDBG_MAX_INTArr(aB.Paths)
			//aB.Strands = GetReverseBoolArr(aB.Strands)
			// map the end partition of read sequence
			//errorNum2, aF := constructdbg.MappingReadToEdgesForWard(dbgK, pairRI[j], int(pos)+cf.Kmerlen, strand, edgesArr, nodesArr, cf.Kmerlen, true)
			//aB.ID, aF.ID = pairRI[j].ID, pairRI[j].ID
			//fmt.Printf("[paraMapNGSAndMerge] aB: %v\n\taF: %v\n", aB, aF)
			//if aB.EndPos < -1 || aF.EndPos < -1 {
			//	fmt.Printf("[paraMapNGSAndMerge] not found end boundary\n")
			//	break
			//}

			if mappingNum < (len(pairRI[j].Seq)-pos)*7/10 {
				needMerge = false
				//fmt.Printf("[paraMapNGSAndMerge] mappingNum: %v < %v\n", mappingNum, len(pairRI[j].Seq)*8/10)
				riArr[j].PathSeqArr = nil
				shortMappingNum++
				//continue
			}

			if errorNum[j] > (len(pairRI[j].Seq)-pos)*1/6 {
				needMerge = false
				//fmt.Printf("[paraMapNGSAndMerge] errorNum: %v > %v\n", errorNum[j], (len(pairRI[j].Seq)-pos)*8/100)
				riArr[j].PathSeqArr = nil
				highErrorNum++
				//break
			}
			if errorNum[j] > (len(pairRI[j].Seq)-pos)*1/10 && len(riArr[j].PathSeqArr) > 1 {
				le := edgesArr[riArr[j].PathSeqArr[len(riArr[j].PathSeqArr)-1].ID]
				if le.GetSeqLen()-kmerlen < 5 {
					riArr[j].PathSeqArr = riArr[j].PathSeqArr[:len(riArr[j].PathSeqArr)-1]
				}
			}
			//fmt.Printf("[paraMapNGSAndMerge] riArr[%v]: %v\n", j, riArr[j])
		}

		/*l0, l1 := len(riArr[0].PathSeqArr), len(riArr[1].PathSeqArr)
		if l0 < 1 || l1 < 1 {
			continue
		}*/

		// find link Path
		//var ri constructcf.ReadInfo
		//ri.ID = pairRI[0].ID
		//fmt.Printf("[paraMapNGSAndMerge] riArr[0]: %v\n\triArr[1]: %v\n", riArr[0], riArr[1])
		// get link path sequence and pass to wc
		// program unknow how many times repeat cycle
		var m constructdbg.ReadMapInfo
		if needMerge && len(riArr[0].PathSeqArr) > 0 && len(riArr[1].PathSeqArr) > 0 {
			lasteID1 := riArr[0].PathSeqArr[len(riArr[0].PathSeqArr)-1]
			lasteID2 := riArr[1].PathSeqArr[len(riArr[1].PathSeqArr)-1]
			laste1, laste2 := edgesArr[lasteID1.ID], edgesArr[lasteID2.ID]
			if laste1.GetTwoEdgesCycleFlag() == 0 && laste1.StartNID != laste1.EndNID {
				m = MergePathSeqArr(riArr[0], riArr[1], edgesArr, nodesArr, kmerlen)
				m.ID = int64(pairRI[0].ID)
			} else if laste2.GetTwoEdgesCycleFlag() == 0 && laste2.StartNID != laste2.EndNID {
				m = MergePathSeqArr(riArr[1], riArr[0], edgesArr, nodesArr, kmerlen)
				m.ID = int64(pairRI[0].ID)
			} else {
				if CountNumInPathSeqArr(lasteID1.ID, riArr[0].PathSeqArr) == 1 {
					m = MergePathSeqArr(riArr[0], riArr[1], edgesArr, nodesArr, kmerlen)
					m.ID = int64(pairRI[0].ID)
				} else if CountNumInPathSeqArr(lasteID2.ID, riArr[1].PathSeqArr) == 1 {
					m = MergePathSeqArr(riArr[1], riArr[0], edgesArr, nodesArr, kmerlen)
					m.ID = int64(pairRI[0].ID)
				} else {
					m = MergePathSeqArr(riArr[0], riArr[1], edgesArr, nodesArr, kmerlen)
					m.ID = int64(pairRI[0].ID)
				}
			}
			if len(m.PathSeqArr) == 0 {
				//fmt.Printf("[paraMapNGSAndMerge] unMerge pair riArr[0]: %v\n\t\triArr[1]: %v\n", riArr[0], riArr[1])
				unMergeNum++
			}
		}

		if len(m.PathSeqArr) < 1 {
			for j := 0; j < 2; j++ {
				var mR constructcf.ReadInfo
				if len(riArr[j].PathSeqArr) < 1 {
					continue
				}
				//mR.Seq = GetPathSeq(riArr[j], edgesArr, nodesArr, cf.Kmerlen)
				mR.Seq = pairRI[j].Seq
				if len(mR.Seq) < 140 {
					continue
				}
				singleNum++
				mR.ID = pairRI[j].ID
				//mR.Anotition = ""
				/*for _, pathSeq := range riArr[j].PathSeqArr {
					mR.Anotition += fmt.Sprintf("%v-", pathSeq.ID)
				}*/
				mR.Anotition = append(mR.Anotition, fmt.Sprintf("EN:%v", errorNum)...)
				//mR.Anotition = fmt.Sprintf("EN:%d", errorNum[j])
				//fmt.Printf("[paraMapNGSAndMerge]len(mR.Seq): %v,  mR.Seq: %v\n", len(mR.Seq), mR.Seq)
				wc <- mR
			}
		} else {
			var mR constructcf.ReadInfo
			//fmt.Printf("[paraMapNGSAndMerge] merge ReadMapInfo: %v\n", m)
			seq := GetPathSeq(m, edgesArr, nodesArr, cf.Kmerlen, riArr)
			if len(seq) <= len(pairRI[0].Seq) || len(seq) < 200 {
				continue
			}
			mR.Seq = seq
			mergeNum++
			mR.ID = uint32(m.ID)
			/*for _, pathSeq := range m.PathSeqArr {
				mR.Anotition += fmt.Sprintf("%v-", pathSeq.ID)
			}*/
			mR.Anotition = append(mR.Anotition, fmt.Sprintf("EN:%v", errorNum)...)
			//mR.Anotition = fmt.Sprintf("EN:%d len:%d", errorNum[0]+errorNum[1], len(mR.Seq))
			//fmt.Printf("[paraMapNGSAndMerge]len(mR.Seq): %v,  mR.Seq: %v\n", len(mR.Seq), mR.Seq)
			wc <- mR
		}
	}
	fmt.Printf("[paraMapNGSAndMerge] not found seed read pair number is : %v,allNum: %v,  percent: %v\n", notFoundSeedNum, allNum, float32(notFoundSeedNum)/float32(allNum))
	fmt.Printf("[paraMapNGSAndMerge] high error mapping read pair number is : %v\tshort mapping number is: %v\tsingle mapping number: %v\t merge number: %v, unMerge number: %v\n", highErrorNum, shortMappingNum, singleNum, mergeNum*2, unMergeNum)
}

func paraMapNGSAndAddFreq(cs <-chan constructcf.ReadInfo, nodesArr []constructdbg.DBGNode, edgesArr []constructdbg.DBGEdge, cf constructdbg.CuckooFilter, winSize, kmerlen int) {
	var notFoundSeedNum, allNum int
	var totalMappingLen int
	MaxStepNum := 10
	ka := make([]constructdbg.KmerInfo, MaxStepNum*winSize)
	rSeq := make([]byte, 300)
	for {
		ri, ok := <-cs
		if !ok {
			break
		}

		allNum++

		// found kmer seed position in the DBG edges
		dbgK, pos, strand := constructdbg.LocateSeedKmerCF(cf, ri, winSize, edgesArr, MaxStepNum, ka, rSeq)
		if dbgK.GetCount() == 0 { // not found in the cuckoofilter
			//fmt.Printf("[paraMapNGSAndMerge] read ID: %v not found seed!!!\n", pairRI[j].ID)
			notFoundSeedNum++
			continue
		}

		mappingNum := MappingReadToEdgesAndAddFreq(dbgK, ri, int(pos), strand, edgesArr, nodesArr, cf.Kmerlen)
		totalMappingLen += mappingNum
	}
	if allNum > 0 && allNum > notFoundSeedNum {
		fmt.Printf("[paraMapNGSAndAddFreq] not found seed read pair number is : %v,allNum: %v,  percent: %v, avgMapping seq len: %v\n", notFoundSeedNum, allNum, float32(notFoundSeedNum)/float32(allNum), totalMappingLen/(allNum-notFoundSeedNum))
	} else {
		fmt.Printf("[paraMapNGSAndAddFreq] not found seed read pair number is : %v,allNum: %v,  \n", notFoundSeedNum, allNum)
	}
	return
}

func writeCorrectReads(browfn string, wc <-chan constructcf.ReadInfo, numCPU, bufSize int) (readNum int) {
	t0 := time.Now()
	fp, err := os.Create(browfn)
	if err != nil {
		log.Fatalf("[writeCorrectReads] failed to create file: %s, err: %v\n", browfn, err)
	}
	defer fp.Close()
	//cbrofp := cbrotli.NewWriter(fp, cbrotli.WriterOptions{Quality: 1})
	cbrofp := cbrotli.NewWriter(fp, cbrotli.WriterOptions{Quality: 1})
	defer cbrofp.Close()
	buffp := bufio.NewWriterSize(cbrofp, 1<<20) // 1<<24 == 2**24
	//buffp := bufio.NewWriter(cbrofp)
	//gzfp := gzip.NewWriter(fp)
	//defer gzfp.Close()
	var finishNum int
	for {
		ri := <-wc
		if ri.ID == 0 {
			finishNum++
			if finishNum == numCPU {
				break
			} else {
				continue
			}
		}
		//seq := constructdbg.Transform2Letters(ri.Seq).String()
		seq := constructdbg.Transform2Char2(ri.Seq)
		str := (*string)(unsafe.Pointer(&seq))
		//s := fmt.Sprintf(">%v\t%s\n%s\n", ri.ID, ri.Anotition, *str)
		fmt.Fprintf(buffp, ">%v\t%s\n%s\n", ri.ID, ri.Anotition, *str)
		//cbrofp.Write([]byte(s))
		//buffp.WriteString(s)
		readNum++
		//if len(wc) > bufSize-100 {
		//	fmt.Printf("[writeCorrectReads] write correct read delay, len(wc): %v\n", len(wc))
		//}
	}
	fmt.Printf("[writeCorrectReads] write correct read num: %v to file: %v, used time: %v\n", readNum, browfn, time.Now().Sub(t0))
	if err := buffp.Flush(); err != nil {
		log.Fatalf("[writeCorrectReads] failed to flush file: %s, err: %v\n", browfn, err)
	}

	if err := cbrofp.Flush(); err != nil {
		log.Fatalf("[writeCorrectReads] failed to flush file: %s, err: %v\n", browfn, err)
	}

	return
}

func paraCorrectReadsFile(fn string, concurrentNum int, nodesArr []constructdbg.DBGNode, edgesArr []constructdbg.DBGEdge, cf constructdbg.CuckooFilter, opt Options, processT chan int) {
	/*idx := strings.LastIndex(fn, "1")
	if idx < 0 {
		idx = strings.LastIndex(fn, "2")
	}*/

	bufSize := 8000
	cs := make(chan constructcf.ReadInfo, bufSize)
	wc := make(chan constructcf.ReadInfo, bufSize)
	//filterLowReadLen := true
	fmt.Printf("[paraCorrectReadsFile]begin process files  fn: %v\n", fn)
	go LoadNGSReads(fn, cs, opt.Kmer)
	for j := 0; j < concurrentNum; j++ {
		go paraMapNGSAndCorrect(cs, wc, nodesArr, edgesArr, cf, opt.WinSize, opt.Kmer)
	}
	// write function
	farr := strings.Split(fn, ".")
	brwfn := strings.Join(farr[:len(farr)-3], ".") + ".Correct." + strings.Join(farr[len(farr)-3:], ".")
	writeNum := writeCorrectReads(brwfn, wc, concurrentNum, bufSize)
	close(wc)
	//time.Sleep(time.Second * 10)
	fmt.Printf("[paraCorrectReadsFile] write correct reads num: %d to file: %s\n", writeNum, brwfn)
	processT <- 1
}

func SplitWrite(wc <-chan [2]constructcf.ReadInfo, wc1, wc2 chan<- constructcf.ReadInfo, numCPU int) {
	count := 0
	for {
		pri := <-wc
		wc1 <- pri[0]
		wc2 <- pri[1]
		if pri[0].ID == 0 {
			count++
		}
		if count == numCPU {
			break
		}
		continue
	}
	//close(wc1)
	//close(wc2)
}

func paraCorrectReadsFilePair(fn1, fn2 string, concurrentNum int, nodesArr []constructdbg.DBGNode, edgesArr []constructdbg.DBGEdge, cf constructdbg.CuckooFilter, opt Options, processT chan int) {
	/*idx := strings.LastIndex(fn, "1")
	if idx < 0 {
		idx = strings.LastIndex(fn, "2")
	}*/

	bufSize := 8000
	cs := make(chan [2]constructcf.ReadInfo, bufSize)
	wc := make(chan [2]constructcf.ReadInfo, bufSize)
	//filterLowReadLen := true
	fmt.Printf("[paraCorrectReadsFilePair]begin process files  fn1: %v, fn2:%v\n", fn1, fn2)
	go LoadNGSReadsPair(fn1, fn2, cs, opt.Kmer, false)
	for j := 0; j < concurrentNum; j++ {
		go paraMapNGSAndCorrectPair(cs, wc, nodesArr, edgesArr, cf, opt.WinSize, opt.Kmer)
	}
	// write function
	wc1 := make(chan constructcf.ReadInfo, bufSize)
	wc2 := make(chan constructcf.ReadInfo, bufSize)
	go SplitWrite(wc, wc1, wc2, concurrentNum)
	farr := strings.Split(fn1, ".")
	brwfn1 := strings.Join(farr[:len(farr)-3], ".") + ".Correct." + strings.Join(farr[len(farr)-3:], ".")
	go writeCorrectReads(brwfn1, wc1, concurrentNum, bufSize)
	farr = strings.Split(fn2, ".")
	brwfn2 := strings.Join(farr[:len(farr)-3], ".") + ".Correct." + strings.Join(farr[len(farr)-3:], ".")
	writeNum := writeCorrectReads(brwfn2, wc2, concurrentNum, bufSize)
	close(wc)
	close(wc1)
	close(wc2)
	for len(wc1) > 0 {
		time.Sleep(time.Second)
	}
	time.Sleep(time.Second)
	fmt.Printf("[paraCorrectReadsFilePair] write correct reads num: %d to file: %s,%s\n", writeNum, brwfn1, brwfn2)
	processT <- 1
}

func paraProcessReadsFile(fn1, fn2 string, concurrentNum int, nodesArr []constructdbg.DBGNode, edgesArr []constructdbg.DBGEdge, cf constructdbg.CuckooFilter, opt Options, processT chan int) {
	idx1 := strings.LastIndex(fn1, "1")
	idx2 := strings.LastIndex(fn2, "2")
	if idx1 < 0 || idx2 < 0 {
		log.Fatalf("[paraProcessReadsFile] fn1: %v,fn2: %v, must used suffix *[1|2].[fasta|fa|fq|fastq].br\n", fn1, fn2)
	}
	if fn1[:idx1] != fn2[:idx2] {
		log.Fatalf("[paraProcessReadsFile] fn1: %v, fn2: %v, must used suffix *[1|2].[fasta|fa|fq|fastq].br\n", fn1, fn2)
	}

	bufSize := 8000
	cs := make(chan [2]constructcf.ReadInfo, bufSize)
	wc := make(chan constructcf.ReadInfo, bufSize)
	filterLowReadLen := true
	fmt.Printf("[paraProcessReadsFile]begin process files  fn1: %v, fn2: %v\n", fn1, fn2)
	go LoadNGSReadsPair(fn1, fn2, cs, opt.Kmer, filterLowReadLen)
	for j := 0; j < concurrentNum; j++ {
		go paraMapNGSAndMerge(cs, wc, nodesArr, edgesArr, cf, opt.WinSize, opt.Kmer)
	}
	// write function
	idx := strings.LastIndex(fn1, "Correct")
	if idx < 0 {
		log.Fatalf("[paraProcessReadsFile] fn1: %v,fn2: %v, must used suffix *.Correct.[1|2].[fasta|fa|fq|fastq].br\n", fn1, fn2)
	}
	brwfn := fn1[:idx] + "Merge.fa.br"
	writeNum := writeCorrectReads(brwfn, wc, concurrentNum, bufSize)
	close(wc)
	//time.Sleep(time.Second * 10)
	fmt.Printf("[paraProcessReadsFile] write correct reads num: %d to file: %s\n", writeNum, brwfn)
	processT <- 1
}

func paraMapReadsFile(fn string, concurrentNum int, nodesArr []constructdbg.DBGNode, edgesArr []constructdbg.DBGEdge, cf constructdbg.CuckooFilter, opt Options, processT chan int) {

	bufSize := 20000
	cs := make(chan constructcf.ReadInfo, bufSize)
	for j := 0; j < concurrentNum; j++ {
		go paraMapNGSAndAddFreq(cs, nodesArr, edgesArr, cf, opt.WinSize, opt.Kmer)
	}
	fmt.Printf("[paraMapReadsFile]begin process files  fn: %v\n", fn)
	LoadNGSReads(fn, cs, opt.Kmer)
	for len(cs) > 0 {
		time.Sleep(time.Second)
	}
	processT <- 1
}

func MappingNGSAndCorrect(opt Options, nodesArr []constructdbg.DBGNode, edgesArr []constructdbg.DBGEdge) {
	// construct cuckoofilter of DBG sample
	// use all edge seq, so set MaxNGSReadLen to MaxInt
	cfSize := constructdbg.GetCuckoofilterDBGSampleSize(edgesArr, int64(opt.WinSize), int64(opt.Kmer))
	fmt.Printf("[MappingNGSAndCorrect] cfSize: %v\n", cfSize)
	cf := constructdbg.MakeCuckooFilter(uint64(cfSize*2), opt.Kmer)
	fmt.Printf("[MappingNGSAndCorrect] cf.numItems: %v\n", cf.NumItems)
	count := constructdbg.ConstructCFDBGMinimizers(cf, edgesArr, opt.WinSize)
	fmt.Printf("[MappingNGSAndCorrect]construct Smaple of DBG edges cuckoofilter number is : %v\n", count)

	cfgInfo, err := constructcf.ParseCfg(opt.CfgFn, opt.Merge, opt.Correct)
	if err != nil {
		log.Fatal("[constructcf.ParseCfg] found err")
	}
	fmt.Printf("[MappingNGSAndCorrect] cfgInfo: %v\n", cfgInfo)

	go runtime.GC()
	runtime.GOMAXPROCS(opt.NumCPU + 10)

	concurrentNum := opt.Kmer / 15
	if concurrentNum < 13 {
		concurrentNum = 13
	}
	totalNumT := (opt.NumCPU + concurrentNum - 1) / concurrentNum

	if opt.NumCPU == 2 {
		concurrentNum = 1
		totalNumT = 1
	}
	//concurrentNum := 1
	// test
	//totalNumT := 1
	//concurrentNum := 6
	processT := make(chan int, totalNumT)
	for i := 0; i < totalNumT; i++ {
		processT <- 1
	}
	for _, lib := range cfgInfo.Libs {
		if lib.AsmFlag != constructcf.AllState && lib.SeqProfile != 1 {
			continue
		}
		//MaxPairLen := lib.InsertSize + lib.InsertSD

		for i := 0; i < len(lib.FnName)-1; i += 2 {
			<-processT
			go paraCorrectReadsFilePair(lib.FnName[i], lib.FnName[i+1], concurrentNum, nodesArr, edgesArr, cf, opt, processT)
		}
		for i := 0; i < len(lib.SingleFn); i++ {
			<-processT
			go paraCorrectReadsFile(lib.SingleFn[i], concurrentNum, nodesArr, edgesArr, cf, opt, processT)
		}
	}

	for i := 0; i < totalNumT; i++ {
		<-processT
	}
	//time.Sleep(time.Second)
}

// MappingNGSAndCorrect() function parallel Map NGS reads to the DBG edges, then merge path output long single read
func MappingNGSAndMerge(opt Options, nodesArr []constructdbg.DBGNode, edgesArr []constructdbg.DBGEdge) {
	// construct cuckoofilter of DBG sample
	// use all edge seq, so set MaxNGSReadLen to MaxInt
	cfSize := constructdbg.GetCuckoofilterDBGSampleSize(edgesArr, int64(opt.WinSize), int64(opt.Kmer))
	fmt.Printf("[MappingNGSAndCorrect] cfSize: %v\n", cfSize)
	cf := constructdbg.MakeCuckooFilter(uint64(cfSize*2), opt.Kmer)
	fmt.Printf("[MappingNGSAndCorrect] cf.numItems: %v\n", cf.NumItems)
	count := constructdbg.ConstructCFDBGMinimizers(cf, edgesArr, opt.WinSize)
	fmt.Printf("[MappingNGSAndCorrect]construct Smaple of DBG edges cuckoofilter number is : %v\n", count)

	cfgInfo, err := constructcf.ParseCfg(opt.CfgFn, opt.Merge, opt.Correct)
	if err != nil {
		log.Fatal("[constructcf.ParseCfg] found err")
	}
	fmt.Printf("[MappingNGSAndCorrect] cfgInfo: %v\n", cfgInfo)

	go runtime.GC()
	runtime.GOMAXPROCS(opt.NumCPU + 10)

	concurrentNum := opt.Kmer / 10
	if concurrentNum < 14 {
		concurrentNum = 14
	}
	totalNumT := (opt.NumCPU + concurrentNum - 1) / concurrentNum

	if opt.NumCPU == 2 {
		concurrentNum = 1
		totalNumT = 1
	}
	//concurrentNum := 1
	// test
	//totalNumT := 1
	//concurrentNum := 6
	processT := make(chan int, totalNumT)
	for i := 0; i < totalNumT; i++ {
		processT <- 1
	}
	for _, lib := range cfgInfo.Libs {
		if lib.AsmFlag != constructcf.AllState && lib.SeqProfile != 1 {
			continue
		}
		//MaxPairLen := lib.InsertSize + lib.InsertSD

		for i := 0; i < len(lib.FnName)-1; i += 2 {
			<-processT
			go paraProcessReadsFile(lib.FnName[i], lib.FnName[i+1], concurrentNum, nodesArr, edgesArr, cf, opt, processT)
		}
		/*for i := 0; i < len(lib.SingleFn); i++ {
			<-processT
			go paraProcessReadsFile(lib.SingleFn[i], "", concurrentNum, nodesArr, edgesArr, cf, opt, lib.InsertSize, lib.InsertSD, processT, true)
		}*/
	}

	for i := 0; i < totalNumT; i++ {
		<-processT
	}
	time.Sleep(time.Second)
}

func PrintCF(cf constructdbg.CuckooFilter) {
	for i, bt := range cf.Hash {
		for j := 0; j < constructdbg.BucketSize; j++ {
			if bt.Bkt[j].ID > 0 {
				fmt.Printf("[PrintCF]cf.Hash[%v].Bkt[%v]: %v, count: %v\n", i, j, bt.Bkt[j], bt.Bkt[j].GetCount())
			}
		}
	}
}

// MappingNGSAndCorrect() function parallel Map NGS reads to the DBG edges, then merge path output long single read
func MappingNGSAndAddFreq(opt Options, nodesArr []constructdbg.DBGNode, edgesArr []constructdbg.DBGEdge) {
	// construct cuckoofilter of DBG sample
	cfSize := constructdbg.GetCuckoofilterDBGSampleSize(edgesArr, int64(opt.WinSize), int64(opt.Kmer))
	fmt.Printf("[MappingNGSAndAddFreq] cfSize: %v\n", cfSize)
	cf := constructdbg.MakeCuckooFilter(uint64(cfSize*2), opt.Kmer)
	fmt.Printf("[MappingNGSAndAddFreq] cf.numItems: %v\n", cf.NumItems)
	count := constructdbg.ConstructCFDBGMinimizers(cf, edgesArr, opt.WinSize)
	fmt.Printf("[MappingNGSAndAddFreq]construct Smaple of DBG edges cuckoofilter number is : %v\n", count)

	cfgInfo, err := constructcf.ParseCfg(opt.CfgFn, opt.Merge, opt.Correct)
	if err != nil {
		log.Fatal("[constructcf.ParseCfg] found err")
	}
	fmt.Printf("[MappingNGSAndAddFreq] cfgInfo: %v\n", cfgInfo)

	runtime.GOMAXPROCS(opt.NumCPU + 10)

	//PrintCF(cf)

	concurrentNum := opt.Kmer / 15
	if concurrentNum < 7 {
		concurrentNum = 7
	}
	//concurrentNum := 1
	totalNumT := (opt.NumCPU + concurrentNum - 1) / concurrentNum
	// test
	//totalNumT := 1
	//concurrentNum := 6
	processT := make(chan int, totalNumT)
	for i := 0; i < totalNumT; i++ {
		processT <- 1
	}
	for _, lib := range cfgInfo.Libs {
		if lib.AsmFlag != constructcf.AllState && lib.SeqProfile != 1 {
			continue
		}

		for i := 0; i < len(lib.FnName); i++ {
			<-processT
			go paraMapReadsFile(lib.FnName[i], concurrentNum, nodesArr, edgesArr, cf, opt, processT)
		}
		for i := 0; i < len(lib.SingleFn); i++ {
			<-processT
			go paraMapReadsFile(lib.SingleFn[i], concurrentNum, nodesArr, edgesArr, cf, opt, processT)
		}
	}

	for i := 0; i < totalNumT; i++ {
		<-processT
	}
	time.Sleep(time.Second)
}

func Correct(c cli.Command) {
	gOpt, suc := utils.CheckGlobalArgs(c.Parent())
	if suc == false {
		log.Fatalf("[Correct] check global Arguments error, opt: %v\n", gOpt)
	}
	opt := Options{gOpt, 0, 0, 0, 0, true, true, 0}
	tmp, suc := checkArgs(c)
	if suc == false {
		log.Fatalf("[Correct] check Arguments error, opt: %v\n", tmp)
	}
	//debug.PrintStack()
	opt.CFSize = tmp.CFSize
	opt.MaxNGSReadLen = tmp.MaxNGSReadLen
	opt.TipMaxLen = tmp.TipMaxLen
	opt.WinSize = tmp.WinSize
	opt.Correct = tmp.Correct
	opt.Merge = tmp.Merge
	opt.MinKmerFreq = tmp.MinKmerFreq
	cfgInfo, err := constructcf.ParseCfg(opt.CfgFn, opt.Merge, opt.Correct)
	if err != nil {
		log.Fatalf("[Correct] ParseCfg 'C': %v err :%v\n", opt.CfgFn, err)
	}
	fmt.Printf("[Correct] opt: %v\n\tcfgInfo: %v\n", opt, cfgInfo)

	/*profileFn := opt.Prefix + ".preprocess.prof"
	cpuprofilefp, err := os.Create(profileFn)
	if err != nil {
		log.Fatalf("[Correct] open cpuprofile file: %v failed\n", profileFn)
	}
	pprof.StartCPUProfile(cpuprofilefp)
	defer pprof.StopCPUProfile()*/
	//runtime.Stack(buf, all)
	//runtime.GoroutineProfile(p)
	// construct cuckoofilter and construct DBG
	nodesfn := opt.Prefix + ".nodes.Arr.br"
	if _, err := os.Stat(nodesfn); err != nil {
		cf := constructcf.CCFFunc(c)
		constructdbg.CDBGFunc(c, cf)
	}
	// smfy DBG
	// read nodes file and transform to array mode for more quckly access
	DBGStatfn := opt.Prefix + ".DBG.stat"
	nodesSize, edgesSize := constructdbg.DBGStatReader(DBGStatfn)
	//nodesSize := len(nodeMap)
	fmt.Printf("[Correct] len(nodesArr): %v, length of edge array: %v\n", nodesSize, edgesSize)
	nodesArr := make([]constructdbg.DBGNode, nodesSize)
	fc := make(chan int, 1)
	constructdbg.NodesArrReader(nodesfn, nodesArr, opt.Kmer, fc)
	// read edges file
	edgesbrfn := opt.Prefix + ".edges.fq.br"
	edgesArr := constructdbg.ReadEdgesFromFile(edgesbrfn, edgesSize)

	var copt constructdbg.Options
	copt.CfgFn, copt.Kmer, copt.MaxNGSReadLen, copt.NumCPU, copt.Prefix, copt.TipMaxLen, copt.WinSize = opt.CfgFn, opt.Kmer, opt.MaxNGSReadLen, opt.NumCPU, opt.Prefix, opt.TipMaxLen, opt.WinSize
	copt.MinMapFreq = opt.MinKmerFreq
	//gfn := opt.Prefix + ".beforeSmfyDBG.dot"
	//constructdbg.GraphvizDBGArr(nodesArr, edgesArr, gfn)
	t0 := time.Now()
	//constructdbg.PrintTmpDBG(nodesArr, edgesArr, opt.Prefix+".beforeSmfy")
	//constructdbg.CheckInterConnectivity(edgesArr, nodesArr)
	for i := 0; i < 5; i++ {
		fmt.Printf("[Correct] %v cycle Smfy DBG....\n", i)
		constructdbg.SmfyDBG(nodesArr, edgesArr, copt, false)
	}
	constructdbg.CheckDBGSelfCycle(nodesArr, edgesArr, copt.Kmer)
	//constructdbg.PrintTmpDBG(nodesArr, edgesArr, opt.Prefix)
	constructdbg.CheckInterConnectivity(edgesArr, nodesArr)

	uniqueNum, bubbleRepeatNum, twoEdgeCycleNum, selfCycleNum, selfCycleSameComingNum, bubbleEdgeNum := constructdbg.SetDBGEdgesUniqueFlag(edgesArr, nodesArr, opt.Kmer)
	fmt.Printf("[Correct] unique edge number is : %v, bubbleRepeatNum: %v, twoEdgeCycleNum: %v,selfCycleNum: %v, selfCycleSameComingNum: %v,bubbleEdgeNum: %v\n", uniqueNum, bubbleRepeatNum, twoEdgeCycleNum, selfCycleNum, selfCycleSameComingNum, bubbleEdgeNum)
	fmt.Printf("[Correct] Smfy DBG used : %v\n", time.Now().Sub(t0))
	t0 = time.Now()
	//gfn1 := opt.Prefix + ".afterSmfyDBG.dot"
	//constructdbg.GraphvizDBGArr(nodesArr, edgesArr, gfn1)

	// Mapping NGS to DBG, Delete low freq edge and Correct
	//copt.MinMapFreq = 5
	// write edges fq
	//ppfqfn := opt.Prefix + ".edges.pp.smfy.fq.br"
	//constructdbg.StoreEdgesToFn(ppfqfn, edgesArr)
	go runtime.GC()
	if opt.Merge == false {
		MappingNGSAndAddFreq(opt, nodesArr, edgesArr)
		fmt.Printf("[Correct] Mapping  NGS and add edge mapping freq info used : %v\n", time.Now().Sub(t0))
		t0 = time.Now()
		DeleteLowFreqEdgeFlag := true
		go runtime.GC()
		constructdbg.SmfyDBG(nodesArr, edgesArr, copt, DeleteLowFreqEdgeFlag)
		MappingNGSAndCorrect(opt, nodesArr, edgesArr)
		fmt.Printf("[Correct] Mapping  NGS and Correct base used : %v\n", time.Now().Sub(t0))
	} else {
		MappingNGSAndMerge(opt, nodesArr, edgesArr)
		fmt.Printf("[Correct] Mapping and Merge NGS reads used : %v\n", time.Now().Sub(t0))
	}

	/*graphfn := opt.Prefix + ".afterSmfy.dot"
	go constructdbg.GraphvizDBGArr(nodesArr, edgesArr, graphfn)*/
	//go runtime.GC()
}
