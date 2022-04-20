package main

import (
	"fmt"
	"log"
	"math"
	"os"
	"runtime"
	"strings"
	"sync"
	"time"

	"github.com/jwaldrip/odin/cli"
	"github.com/klauspost/compress/zstd"
)

type optionsPP struct {
	ArgsOpt
	CFSize        int
	TipMaxLen     int
	WinSize       int
	MaxNGSReadLen int
	MinNGSReadLen int
	Correct       bool
	Merge         bool
	MinKmerFreq   int
}

func checkArgsPP(c cli.Command) (opt optionsPP, suc bool) {
	gOpt, suc := CheckGlobalArgs(c.Parent())
	if suc == false {
		log.Fatalf("[Correct] check global Arguments error, opt: %v\n", gOpt)
	}
	opt.ArgsOpt = gOpt
	opt.CFSize = c.Flag("S").Get().(int)
	if opt.CFSize < 1024*1024 {
		log.Fatalf("the argument 'S':%s must bigger than 1024 * 1024\n", c.Flag("S"))
	}
	opt.MaxNGSReadLen = c.Flag("MaxNGSReadLen").Get().(int)
	opt.MinNGSReadLen = c.Flag("MinNGSReadLen").Get().(int)
	cor := c.Flag("Correct").Get().(bool)
	if cor == true {
		opt.Correct = true
	} else {
		//log.Fatalf("[checkArgs] argument 'Correct': %v set error, must set 'true'\n", c.Flag("Correct"))
	}
	merge := c.Flag("Merge").Get().(bool)
	opt.Merge = merge
	opt.WinSize = c.Flag("WinSize").Get().(int)
	if opt.WinSize < 1 || opt.WinSize > 20 {
		log.Fatalf("[checkArgs] argument 'WinSize':%s, must bewteen [3~20]\n", c.Flag("WinSize"))
	}
	if 81 < opt.Kmer && opt.Kmer < 149 && opt.Kmer%2 == 0 {
		log.Fatalf("the argument 'K': %v must between [81~149] and tmp must been even\n", c.Parent().Flag("K"))
	}

	opt.TipMaxLen = c.Flag("TipMaxLen").Get().(int)

	if opt.TipMaxLen == 0 {
		opt.TipMaxLen = opt.MaxNGSReadLen
	} else {
		if opt.TipMaxLen < 0 || opt.TipMaxLen < opt.Kmer*2 || opt.TipMaxLen > opt.MaxNGSReadLen {
			log.Fatalf("[checkArgs] argument 'TipMaxLen':%d must between [%d~%d]\n", opt.TipMaxLen, opt.Kmer*2, opt.MaxNGSReadLen)
		}
	}
	opt.MinKmerFreq = c.Flag("MinKmerFreq").Get().(int)
	if opt.MinKmerFreq < 1 || opt.MinKmerFreq > MaxC {
		log.Fatalf("[checkArgs]the argument 'MinKmerFreq':%d must [%d~%d]\n", opt.MinKmerFreq, 1, MaxC)
	}
	suc = true
	return opt, suc
}

/*func LoadNGSReads(brfn string, cs chan<- ReadInfo, kmerlen int) {
	bufSize := (1 << 20)
	//size := 8000
	fncs := make(chan []byte, 10)
	//recordChan := make(chan ReadInfo, size)
	format := GetReadsFileFormat(brfn)

	go ReadBrFile(brfn, fncs, bufSize)
	count := GetReadFileRecord(fncs, cs, format, bufSize, true)

	fmt.Printf("[LoadNGSReads] processed %d reads from files: %s\n", count, brfn)
	//time.Sleep(time.Second * 10)
}*/

/*func LoadNGSReadsPair(brfn1, brfn2 string, cs chan<- [2]ReadInfo, kmerlen int, filterLowReadLen bool) {
	bufSize := (1 << 20)
	size := 8000
	fncs1 := make(chan []byte, 10)
	recordChan1 := make(chan ReadInfo, size)
	format1 := GetReadsFileFormat(brfn1)

	go ReadBrFile(brfn1, fncs1, bufSize)
	go GetReadFileRecord(fncs1, recordChan1, format1, bufSize, true)

	format2 := GetReadsFileFormat(brfn2)
	fncs2 := make(chan []byte, 10)
	go ReadBrFile(brfn2, fncs2, bufSize)
	recordChan2 := make(chan ReadInfo, size)
	go GetReadFileRecord(fncs2, recordChan2, format2, bufSize, true)

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
		var pairRI [2]ReadInfo
		pairRI[0].ID = ri1.ID
		pairRI[1].ID = ri2.ID
		//if len(ri1.Seq) < kmerlen || len(ri2.Seq) < kmerlen {
		//	continue
		//}
		pairRI[0].Seq = ri1.Seq
		pairRI[1].Seq = ri2.Seq

		cs <- pairRI
		count++
	}
	fmt.Printf("[LoadNGSReadsPair] processed %d pair reads from pair files: %s , %s\n", count, brfn1, brfn2)
	close(cs)
	//time.Sleep(time.Second * 10)
}*/

/*func GetExtendPathArr(path []uint32, nID uint32, edgesArr []DBGEdge, nodesArr []DBGNode, kmerlen, extLen int) (epArr [][]uint32) {
	type pathExtend struct {
		Path      []uint32
		NID       uint32
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
		ea := GetNearEdgeIDArr(nodesArr[lpe.NID], lpe.Path[len(lpe.Path)-1])
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
			npe.ExtendLen = lpe.ExtendLen + len(ne.Ks) - (kmerlen - 1)
			pathQueue = append(pathQueue, npe)
		}
	}

	return
}*/

func ReversePathSeqArr(psArr []PathSeq) {
	la := len(psArr)
	for i := 0; i < la/2; i++ {
		tmp := psArr[i]
		psArr[i].ID, psArr[i].Strand = psArr[la-i-1].ID, !psArr[la-i-1].Strand
		psArr[la-i-1].ID, psArr[la-i-1].Strand = tmp.ID, !tmp.Strand
	}
	if la%2 == 1 {
		i := la / 2
		psArr[i].Strand = !psArr[i].Strand
	}
}

func EqualPathSeqArr(a1, a2 []PathSeq) (ok bool) {
	if len(a1) != len(a2) {
		return
	}
	ok = true
	for i := 0; i < len(a1); i++ {
		if a1[i].ID != a2[i].ID || a1[i].Strand != a2[i].Strand {
			ok = false
			break
		}
	}
	return
}

/*func MappingReadToEdgesAndAddFreq(dk DBGKmer, ri ReadInfo, rpos int, rstrand bool, edgesArr []DBGEdge, nodesArr []DBGNode, kmerlen int) (mappingNum int) {
	var strand bool
	if dk.Strand == rstrand {
		strand = PLUS
		dk.Pos += (int32(kmerlen) - 1)
	} else {
		strand = MINUS
	}
	rpos += kmerlen - 1
	mappingNum += kmerlen
	e := &edgesArr[dk.ID]
	for i := rpos; i < len(ri.Seq); i++ {
		edgesArr[dk.ID].Kq[dk.Pos]++
		if i == len(ri.Seq)-1 {
			break
		}
		if strand == PLUS {
			if dk.Pos < int32(e.GetSeqLen())-1 {
				eb := e.Ks[dk.Pos+1]
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
				eIDArr := GetNextEdgeIDArr(e, nodesArr, strand, FORWARD)
				if len(eIDArr) == 0 {
					break
				}
				ok := false
				for j := 0; j < len(eIDArr); j++ {
					ne := &edgesArr[eIDArr[j]]
					strand = GetNextEdgeStrand(e, ne, nodesArr, strand)
					e = ne
					dk.ID = ne.ID
					var eb byte
					if strand {
						dk.Pos = int32(kmerlen) - 1
						eb = e.Ks[dk.Pos]
					} else {
						dk.Pos = int32(e.GetSeqLen() - kmerlen)
						eb = BntRev[e.Ks[dk.Pos]]
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
				eb := BntRev[e.Ks[dk.Pos-1]]
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
				eIDArr := GetNextEdgeIDArr(e, nodesArr, strand, FORWARD)
				if len(eIDArr) == 0 {
					break
				}
				ok := false
				for j := 0; j < len(eIDArr); j++ {
					ne := &edgesArr[eIDArr[j]]
					strand = GetNextEdgeStrand(e, ne, nodesArr, strand)
					e = ne
					dk.ID = ne.ID
					var eb byte
					if strand {
						dk.Pos = int32(kmerlen) - 1
						eb = e.Ks[dk.Pos]
					} else {
						dk.Pos = int32(e.GetSeqLen() - kmerlen)
						eb = BntRev[e.Ks[dk.Pos]]
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
}*/

type PathSeqInfo struct {
	Seq     []byte
	Path    []uint32
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

func MappingReadToEdges(dk, rdk DBGKmer, ri ReadInfo, edgesArr []DBGEdge, kmerlen int) (int, int, ReadMapInfo, ReadInfo) {
	var errorNum, mappingNum int
	var rmi ReadMapInfo
	rmi.PathSeqArr = make([]PathSeq, 0, 5)
	var strand bool
	var ps PathSeq
	if dk.GetStrand() == rdk.GetStrand() {
		rmi.StartP = int(dk.GetPos())
		dk.SetPos(uint64(rmi.StartP) + uint64(kmerlen))
		strand = PLUS
	} else {
		rmi.StartP = int(dk.GetPos()) + kmerlen - 1
		//dk.Pos-- // need check if Pos == 0
		strand = MINUS
	}
	rdk.SetPos(rdk.GetPos() + uint64(kmerlen))
	mappingNum += int(kmerlen)
	eIDArr := make([]uint32, 0, BaseTypeNum)

	//fmt.Printf("[MappingReadToEdges]len(ri.Seq): %v\n", len(ri.Seq))
	for i := int(rdk.GetPos()); i < len(ri.Seq); {
		e := &edgesArr[dk.GetID()]
		b := len(ri.Seq)
		var j int
		//fmt.Printf("[MappingReadToEdges]i: %v,  strand: %v, dk.Pos: %v, dk.ID: %v, errorNum: %v\n", i, strand, dk.Pos, dk.ID, errorNum)
		if strand == PLUS {
			if len(e.Ks)-int(dk.GetPos()) < len(ri.Seq)-i {
				b = i + (len(e.Ks) - int(dk.GetPos()))
			}
			j = int(dk.GetPos())
			for ; i < b && j < len(e.Ks); i++ {
				if ri.Seq[i] != e.Ks[j] {
					/*if e.Ks[j] > 3 {
						fmt.Printf("[MappingReadToEdges] b: %v,i: %v,  j: %v, es:%v, eID:%d\n\te: %v", b, i, j, e.Ks[j], e.ID, e)
					}*/
					errorNum++
					ri.Seq[i] = e.Ks[j]
				} else {
					mappingNum++
				}
				j++
			}
		} else { // strand == MINUS
			if len(ri.Seq)-i > int(dk.GetPos()) {
				b = i + int(dk.GetPos())
			}
			j = int(dk.GetPos()) - 1
			for ; i < b && j >= 0; i++ {
				/*if e.Ks[j] > 3 {
					fmt.Printf("[MappingReadToEdges] b: %v,i: %v,  j: %v, es:%v, eID:%d\n\te: %v", b, i, j, e.Ks[j], e.ID, e)
				}*/
				if ri.Seq[i] != BntRev[e.Ks[j]] {
					errorNum++
					ri.Seq[i] = BntRev[e.Ks[j]]
				} else {
					mappingNum++
				}
				j--
			}
		}

		/*if errorNum > 0 {
			fmt.Printf("[MappingReadToEdges]not perfect end i: %v,edge ID: %v,len(e.Ks): %v,  dk.Pos: %v, pos: %v, b: %v\n", i, dk.ID, len(e.Ks), dk.Pos, rpos, b)
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
			rdk.SetPos(uint64(i))
			eIDArr = GetNextEdgeIDArr(e, strand, FORWARD, eIDArr)
			if len(eIDArr) == 0 {
				ri.Seq = ri.Seq[:i]
				//ri.Qual = ri.Qual[:i]
				errorNum += len(ri.Seq) - i
				//fmt.Printf("[MappingReadToEdges]len(eIDArr):%d\n", len(eIDArr))
				break
			} else if len(eIDArr) == 1 {
				dk.SetID(uint64(eIDArr[0]))
				ne := &edgesArr[dk.GetID()]
				strand = GetNextEdgeStrand2(e, ne, strand, FORWARD)
				if strand {
					dk.SetPos(uint64(kmerlen) - 1)
				} else {
					dk.SetPos(uint64(ne.GetSeqLen() - (kmerlen - 1)))
				}
				//fmt.Printf("[MappingReadToEdges]i:%d, dk.Pos:%d errorNum:%d\n", i, dk.Pos, errorNum)
				continue
			}

			// alignment all eIDArr and choose best score
			seqArr := make([][]byte, len(eIDArr))
			strandArr := make([]bool, len(eIDArr))
			for x := 0; x < len(eIDArr); x++ {
				id := eIDArr[x]
				ne := &edgesArr[id]
				std := GetNextEdgeStrand2(e, ne, strand, FORWARD)
				strandArr[x] = std
				maxLen := len(ri.Seq) - i
				if ne.GetSeqLen()-(kmerlen-1) <= maxLen {
					maxLen = ne.GetSeqLen() - (kmerlen - 1)
				}
				if std {
					seqArr[x] = ne.Ks[kmerlen-1 : kmerlen-1+maxLen]
				} else {
					seqArr[x] = ne.Ks[ne.GetSeqLen()-(kmerlen-1)-maxLen : ne.GetSeqLen()-(kmerlen-1)]
				}
			}

			minLen := GetMinLenByteArr(seqArr)
			var maxScore, nextMaxScore int
			maxIdx := -1
			for x, a := range seqArr {
				s := CountEqualByte(a, ri.Seq[i:i+minLen], strandArr[x])
				if s > maxScore {
					maxScore = s
					maxIdx = x
					nextMaxScore = 0
				} else if s == maxScore {
					nextMaxScore = s
				}
			}

			if maxScore > 0 && nextMaxScore == 0 {
				dk.SetID(uint64(eIDArr[maxIdx]))
				strand = strandArr[maxIdx]
				if strand {
					dk.SetPos(uint64(kmerlen) - 1)
				} else {
					dk.SetPos(uint64(edgesArr[dk.GetID()].GetSeqLen() - (kmerlen - 1)))
				}
				//fmt.Printf("[MappingReadToEdges]maxScore:%d i:%d, dk.Pos:%d errorNum:%d\n", maxScore, i, dk.Pos, errorNum)
			} else {
				/*fmt.Printf("[MappingReadToEdges]minLen:%d maxScore:%d nextMaxScore:%d\n\tri.Seq[%d:%d]:%v\n", minLen, maxScore, nextMaxScore, i, i+minLen, ri.Seq[i:i+minLen])
				for x, a := range seqArr {
					if strandArr[x] {
						fmt.Printf("std:%v seq:%v\n", strandArr[x], a[:minLen])
					} else {
						fmt.Printf("std:%v seq:%v\n", strandArr[x], GetReverseCompByteArr(a[len(a)-minLen:]))
					}
				}*/
				break
			}

			// need found all  paths seq for alignment
			/*{
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
					eIDArr := GetNextEdgeIDArr(&edgesArr[eID], nodesArr, ps.Strand, FORWARD)
					for x := 0; x < len(eIDArr); x++ {
						var nps PathSeqInfo
						id := eIDArr[x]
						ne := &edgesArr[id]
						std := GetNextEdgeStrand(&edgesArr[eID], ne, nodesArr, ps.Strand)
						nps.Strand = std
						nps.Path = make([]uint32, len(ps.Path)+1)
						copy(nps.Path, ps.Path)
						nps.Path[len(nps.Path)-1] = id
						nl := len(ne.Ks) - (kmerlen - 1)
						if std {
							if nl >= ps.NeedLen {
								nl = ps.NeedLen
							}
							nps.EndPos = (kmerlen - 1) + nl
							nps.NeedLen = ps.NeedLen - nl
							//fmt.Printf("[MappingReadToEdges] nl: %v, el: %v\n", nl, ne.GetSeqLen())
							nps.Seq = make([]byte, len(ps.Seq)+nl)
							copy(nps.Seq, ps.Seq)
							copy(nps.Seq[len(ps.Seq):], ne.Ks[kmerlen-1:kmerlen-1+nl])
						} else {
							if nl >= ps.NeedLen {
								nl = ps.NeedLen
							}
							nps.EndPos = len(ne.Ks) - (kmerlen - 1) - nl - 1
							nps.NeedLen = ps.NeedLen - nl
							nps.Seq = make([]byte, len(ps.Seq)+nl)
							copy(nps.Seq, ps.Seq)
							copy(nps.Seq[len(ps.Seq):], GetReverseCompByteArr(ne.Ks[len(ne.Ks)-(kmerlen-1)-nl:len(ne.Ks)-(kmerlen-1)]))
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
								var ps PathSeq
								ps.ID = psi.Path[x]
								ne := &edgesArr[ps.ID]
								std = GetNextEdgeStrand(le, ne, nodesArr, std)
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
			}*/
		}
	}
	return errorNum, mappingNum, rmi, ri
}

type Item struct {
	Strand    bool
	ExtendLen int
	EIDArr    []uint32
}

func GetPairReadsLinkPathArr(e1, e2 DBGEdge, strand bool, maxAllowLen int, nodesArr []DBGNode, edgesArr []DBGEdge, kmerlen int) (extPathArr [][]uint32) {
	var stack []Item
	var t Item
	t.EIDArr = append(t.EIDArr, e1.ID)
	t.Strand = strand
	stack = append(stack, t)
	ea := make([]uint32, 0, BaseTypeNum)

	// process stack
	for len(stack) > 0 {
		t := stack[len(stack)-1]
		stack = stack[:len(stack)-1]
		e := &edgesArr[t.EIDArr[len(t.EIDArr)-1]]
		std := t.Strand

		ea = GetNextEdgeIDArr(e, std, FORWARD, ea)
		for _, id := range ea {
			var nt Item
			nt.EIDArr = make([]uint32, len(t.EIDArr)+1)
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
			nt.ExtendLen = t.ExtendLen + (len(e.Ks) - (kmerlen - 1))
			if nt.ExtendLen >= maxAllowLen {
				continue
			}
			stack = append(stack, nt)
		}
	}

	return
}

func MergePathSeqArr(rmi1, rmi2 ReadMapInfo, edgesArr []DBGEdge, kmerlen int) (mergeRMI ReadMapInfo) {
	psa1 := rmi1.PathSeqArr
	psa2 := make([]PathSeq, len(rmi2.PathSeqArr))
	copy(psa2, rmi2.PathSeqArr)
	ReversePathSeqArr(psa2)
	//e1, e2 := edgesArr[lastT1.ID], edgesArr[lastT2.ID]
	//strand1, strand2 := lastT1.Strand, lastT2.Strand
	mergeRMI.StartP = rmi1.StartP
	mergeRMI.EndP = rmi2.StartP
	// reverse pathSeqArr2
	//fmt.Printf("[MergePathSeqArr]rmi1: %v\n\t\trmi1: %v\n", rmi1, rmi2)
	//fmt.Printf("[MergePathSeqArr]psArr1: %v\n\t\tpsArr2: %v\n", pathSeqArr1, pathSeqArr2)
	var sharePathLen int
	min := MinInt(len(psa1), len(psa2))
	for i := len(psa1) - min; i < len(psa1); i++ {
		if EqualPathSeqArr(psa1[i:], psa2[:len(psa1)-i]) {
			sharePathLen = len(psa1) - i
			break
		}
	}

	if sharePathLen > 0 {
		mergeRMI.PathSeqArr = append(mergeRMI.PathSeqArr, psa1...)
		mergeRMI.PathSeqArr = append(mergeRMI.PathSeqArr, psa2[sharePathLen:]...)
		// check not a pair end reads
		if len(mergeRMI.PathSeqArr) == 1 {
			if mergeRMI.PathSeqArr[0].Strand == PLUS {
				if mergeRMI.StartP > mergeRMI.EndP {
					fmt.Fprintf(os.Stderr, "[MergePathSeqArr]Postion error, p1:%v p2:%v mergeRMI:%v\n", psa1, psa2, mergeRMI)
					mergeRMI.PathSeqArr = nil
				}
			} else {
				if mergeRMI.StartP < mergeRMI.EndP {
					fmt.Fprintf(os.Stderr, "[MergePathSeqArr]Postion error, p1:%v p2:%v mergeRMI:%v\n", psa1, psa2, mergeRMI)
					mergeRMI.PathSeqArr = nil
				}
			}
		}
	} else {
		l1, l2 := psa1[len(psa1)-1], psa2[0]
		e1, e2 := &edgesArr[l1.ID], &edgesArr[l2.ID]
		var coming1, coming2 []uint32
		if l1.Strand == PLUS {
			coming1 = GetComingEdgeArr(e1.EdgeIDOutcoming)
		} else {
			coming1 = GetComingEdgeArr(e1.EdgeIDIncoming)
		}
		if l2.Strand == PLUS {
			coming2 = GetComingEdgeArr(e2.EdgeIDIncoming)
		} else {
			coming2 = GetComingEdgeArr(e2.EdgeIDOutcoming)
		}
		if len(coming1) == 1 && len(coming2) == 1 && coming1[0] == coming2[0] {
			e := &edgesArr[coming1[0]]
			st1, st2 := GetNextEdgeStrand2(e1, e, l1.Strand, FORWARD), GetNextEdgeStrand2(e2, e, l2.Strand, BACKWARD)
			if st1 == st2 {
				var ps PathSeq
				ps.ID = e.ID
				ps.Strand = st1
				mergeRMI.PathSeqArr = append(mergeRMI.PathSeqArr, psa1...)
				mergeRMI.PathSeqArr = append(mergeRMI.PathSeqArr, ps)
				mergeRMI.PathSeqArr = append(mergeRMI.PathSeqArr, psa2...)
			}
		}
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
			var t PathSeq
			id := extPathArr[0][x]
			e := edgesArr[extPathArr[0][x-1]]
			ne := edgesArr[id]
			t.ID = id
			strand = GetNextEdgeStrand(e, ne, nodesArr, strand)
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
	return
}

func GetPathSeq(m ReadMapInfo, edgesArr []DBGEdge, kmerlen int, riArr [2]ReadMapInfo, seq []byte) []byte {
	//fmt.Printf("[GetPathSeq] m: %v\n", m)
	strand := m.PathSeqArr[0].Strand
	seq = seq[:0]
	if len(m.PathSeqArr) == 1 {
		ps := m.PathSeqArr[0]
		if strand == PLUS {
			//seq = make([]byte, m.EndP-m.StartP)
			seq = append(seq, edgesArr[ps.ID].Ks[m.StartP:m.EndP+1]...)
		} else {
			/*if m.EndP < 0 || m.StartP > len(edgesArr[ps.ID].Ks) {
				fmt.Printf("[GetPathSeq] error set m: %v, len(e.Ks): %v\n", m, len(edgesArr[ps.ID].Ks))
				seq = nil
				return
			}*/
			//fmt.Printf("[GetPathSeq] m: %v, len(e.Ks): %v\n", m, len(edgesArr[ps.ID].Ks))
			seq = append(seq, edgesArr[ps.ID].Ks[m.EndP:m.StartP+1]...)
			ReverseCompByteArr(seq)
		}
		return seq
	}

	//var strand bool // true denote "+" DNA strand map, false denote  "-" DNA strand map
	//var coming bool // true denote EdgeOutcoming, false denote EdgeIncoming
	for i := 0; i < len(m.PathSeqArr); i++ {
		ps := m.PathSeqArr[i]
		e := &edgesArr[ps.ID]
		if strand != ps.Strand {
			fmt.Fprintf(os.Stderr, "[GetPathSeq]m.PathSeqArr[%d]:%v strand:%v Cycle:%v\n", i, m.PathSeqArr[i], strand, edgesArr[m.PathSeqArr[i].ID].StartNID == edgesArr[m.PathSeqArr[i].ID].EndNID)
		}
		if i == 0 {
			if strand == PLUS {
				seq = append(seq, e.Ks[m.StartP:]...)
			} else {
				seq = append(seq, GetReverseCompByteArr(e.Ks[:m.StartP+1])...)
			}
		} else if i == len(m.PathSeqArr)-1 {
			//fmt.Printf("[GetPathSeq]len(e.Ks): %v, strand: %v, m.EndP: %v\n", len(e.Ks), strand, m.EndP)
			if strand == PLUS {
				if m.EndP < kmerlen-1 {
					fmt.Fprintf(os.Stderr, "[GetPathSeq]m.EndP: %v set error, len(e.Ks): %v, strand: %v, new m.EndP: %v\n\tm: %v\n\triArr[0]: %v\n\triArr[1]: %v\n", m.EndP, len(e.Ks), strand, kmerlen-1, m, riArr[0], riArr[1])
					m.EndP = kmerlen - 1
				} else if m.EndP > len(e.Ks) {
					fmt.Fprintf(os.Stderr, "[GetPathSeq]m.EndP: %v set error, len(e.Ks): %v, strand: %v, new m.EndP: %v\n\tm: %v\n\triArr[0]: %v\n\triArr[1]: %v\n", m.EndP, len(e.Ks), strand, len(e.Ks), m, riArr[0], riArr[1])
					m.EndP = len(e.Ks)
				}
				seq = append(seq, e.Ks[kmerlen-1:m.EndP+1]...)
			} else {
				if m.EndP < 0 {
					fmt.Fprintf(os.Stderr, "[GetPathSeq]m.EndP: %v set error, len(e.Ks): %v, strand: %v, new m.EndP: %v\n\tm: %v\n\triArr[0]: %v\n\triArr[1]: %v\n", m.EndP, len(e.Ks), strand, 0, m, riArr[0], riArr[1])
					m.EndP = 0
				} else if m.EndP > len(e.Ks)-(kmerlen-1) {
					fmt.Fprintf(os.Stderr, "[GetPathSeq]m.EndP: %v set error, len(e.Ks): %v, strand: %v, new m.EndP: %v\n\tm: %v\n\triArr[0]: %v\n\triArr[1]: %v\n", m.EndP, len(e.Ks), strand, len(e.Ks)-(kmerlen-1), m, riArr[0], riArr[1])
					m.EndP = len(e.Ks) - (kmerlen - 1)
				}
				seq = append(seq, GetReverseCompByteArr(e.Ks[m.EndP:len(e.Ks)-(kmerlen-1)])...)
			}
		} else {
			if strand == PLUS {
				seq = append(seq, e.Ks[kmerlen-1:]...)
			} else {
				seq = append(seq, GetReverseCompByteArr(e.Ks[:len(e.Ks)-(kmerlen-1)])...)
			}
		}
		if i < len(m.PathSeqArr)-1 {
			ne := &edgesArr[m.PathSeqArr[i+1].ID]
			strand = GetNextEdgeStrand2(e, ne, strand, FORWARD)
		}
	}

	/*strand := m.Strands[0]
	e0, e1 := edgesArr[m.Path[0]], edgesArr[m.Path[1]]
	//if m.NID == 0 {
	//	log.Fatalf("[GetPathSeq] m.NID: %v not set\n", m.NID)
	//}
	//nID := m.NID
	var nID uint32
	if strand == PLUS {
		seq = e0.Ks[m.StartP+1:]
		nID = e0.EndNID
	} else {
		seq = GetReverseCompByteArr(e0.Ks[:m.StartP])
		nID = e0.StartNID
	}
	for _, eID := range m.Path[1 : len(m.Path)-1] {
		e1 = edgesArr[eID]
		if e1.StartNID == e1.EndNID {
			if strand == PLUS {
				if reflect.DeepEqual(e0.Ks[len(e0.Ks)-kmerlen+1:], e1.Ks[:kmerlen-1]) {

				} else {
					strand = !strand
				}
			} else {
				if reflect.DeepEqual(e0.Ks[:kmerlen-1], e1.Ks[len(e1.Ks)-kmerlen+1:]) {

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
		if strand == PLUS {
			es = e1.Ks[kmerlen-1:]
		} else {
			es = GetReverseCompByteArr(e1.Ks[:len(e1.Ks)-(kmerlen-1)])
		}
		seq = append(seq, es...)
		e0 = e1
	}

	e1 = edgesArr[m.Path[len(m.Path)-1]]
	if e1.StartNID == e1.EndNID {
		if strand == PLUS {
			if reflect.DeepEqual(e0.Ks[len(e0.Ks)-kmerlen+1:], e1.Ks[:kmerlen-1]) {

			} else {
				strand = !strand
			}
		} else {
			if reflect.DeepEqual(e0.Ks[:kmerlen-1], e1.Ks[len(e1.Ks)-kmerlen+1:]) {

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
	//fmt.Printf("[GetPathSeq] len(e.Ks): %v, e: %v\n", len(e1.Ks), e1)
	var es []byte
	if strand == PLUS {
		es = e1.Ks[kmerlen-1 : m.EndP]
	} else {
		es = GetReverseCompByteArr(e1.Ks[m.EndP+1 : len(e1.Ks)-(kmerlen-1)])
	}
	seq = append(seq, es...) */

	return seq
}

/*func WriteChanPairReads(riArr [2]ReadMapInfo, edgesArr []DBGEdge, kmerlen int, wc chan<- ReadInfo) {
	for x := 0; x < 2; x++ {
		var mr ReadInfo
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
		if len(e.Ks) > MaxPairLen-2*SD {
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
		} else if IsBubble(e1.ID, e2.ID, edgesArr) {
			merged = false
		}

		if merged == false {
			WriteChanPairReads(riArr, edgesArr, cf.Kmerlen, wc)
			continue
		}
	}
	extLen := cf.Kmerlen - (len(pairRI[0].Seq) + len(pairRI[1].Seq) - MaxPairLen) + 50
	var nID uint32
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
		//nID = GetInterNodeID(e1, e2)
		if nID == e2.StartNID {
			nID = e2.EndNID
			extLen -= (len(e2.Ks) - riArr[0].EndP)
		} else {
			nID = e2.StartNID
			extLen -= (riArr[0].EndP)
		}
	} else {
		e := edgesArr[riArr[0].Path[0]]
		if riArr[0].StartP < riArr[0].EndP {
			nID = e.EndNID
			extLen -= (len(e.Ks) - riArr[0].EndP)
		} else {
			nID = e.StartNID
			extLen -= (riArr[0].EndP)
		}
	}
	var rev1 ReadMapInfo
	rev1.Path = GetReverseuint32Arr(riArr[1].Path)
	rev1.StartP, rev1.EndP = riArr[1].EndP, riArr[1].StartP
	var linkPathArr [][]uint32
	if extLen <= 0 {
		min := Min(len(riArr[0].Path), len(rev1.Path))
		for j := len(riArr[0].Path) - min; j < len(riArr[0].Path); j++ {
			if reflect.DeepEqual(riArr[0].Path[j:], rev1.Path[:len(riArr[0].Path)-j]) {
				var mergeP []uint32
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
			idx := IndexEID(rev1.Path, ep[len(ep)-1])
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

func CountNumInPathSeqArr(eID uint32, arr []PathSeq) (count int) {
	for _, item := range arr {
		if eID == item.ID {
			count++
		}
	}
	return count
}

func paraMapNGSAndCorrect(cs <-chan RIPool, wc chan<- RIPool, nodesArr []DBGNode, edgesArr []DBGEdge, cf CuckooFilterDBGKmer, opt optionsPP) {
	var notFoundSeedNum, shortMappingNum, highErrorNum, allNum int
	var oneMoreLocationNum int
	//MaxStepNum := opt.WinSize
	ka := make([]DBGKmerSeq, opt.WinSize)
	rSeq := make([]byte, opt.MaxNGSReadLen)
	dbgKArr := make([]DBGKmer, 2)
	var rdbgK DBGKmer
	for {
		riPool, ok := <-cs
		if !ok {
			var t RIPool
			wc <- t
			break
		}

		for i, ri := range riPool.RIArr {
			// found kmer seed position in the DBG edges
			if len(ri.Seq) < opt.Kmer+opt.WinSize*2 {
				continue
			}
			dbgKArr, rdbgK = LocateSeedKmerCF(cf, ri, opt.WinSize, opt.Kmer, edgesArr, ka, rSeq, dbgKArr)
			if len(dbgKArr) == 0 { // not found in the cuckoofilter
				//fmt.Printf("[paraMapNGSAndMerge] read ID: %v not found seed!!!\n", pairRI[j].ID)
				notFoundSeedNum++
				continue
			} else if len(dbgKArr) > 1 {
				dbgKArr[0] = ReLocation(dbgKArr, rdbgK, ri, edgesArr, opt.Kmer)
				if dbgKArr[0].GetID() == 0 {
					oneMoreLocationNum++
					//fmt.Fprintf(os.Stderr, "[paraMapNGSAndCorrect]read ID:%d ReLocation failed!!!\n", ri.ID)
					continue
				}
			}

			//fmt.Printf("[paraMapNGSAndMerge] pairRI[%v]: %v\n", j, pairRI[j])

			//dk := dbgKArr[0]
			//fmt.Printf("[paraMapNGSAndMerge]readID: %v, dbgK: %v, pos: %v, strand: %v\n", pairRI[0].ID, dbgK, pos, strand)
			errorNum, mappingNum, _, correctRI := MappingReadToEdges(dbgKArr[0], rdbgK, ri, edgesArr, opt.Kmer)
			// extend seed map to the edges

			if mappingNum < len(ri.Seq)*7/10 {
				//fmt.Printf("[paraMapNGSAndMerge]errorNum:%d mappingNum:%d rmi:%v\n", errorNum, mappingNum, rmi)
				shortMappingNum++
				continue
			}

			if errorNum > len(ri.Seq)*1/10 {
				//fmt.Printf("[paraMapNGSAndMerge]errorNum:%d mappingNum:%d rmi:%v\n", errorNum, mappingNum, rmi)
				//fmt.Printf("[paraMapNGSAndMerge] errorNum: %v > %v\n", errorNum[j], (len(pairRI[j].Seq)-pos)*8/100)
				highErrorNum++
				continue
			}

			correctRI.ID = ri.ID
			correctRI.Qual = ri.Qual
			correctRI.Anotition = append(correctRI.Anotition, fmt.Sprintf("EN:%d", errorNum)...)
			riPool.RIArr[i] = correctRI
			//fmt.Printf("[paraMapNGSAndMerge]len(mR.Seq): %v,  mR.Seq: %v\n", len(mR.Seq), mR.Seq)
		}

		wc <- riPool

		allNum += len(riPool.RIArr)

		// map the start partition of read sequence
		//errorNum1, aB := MappingReadToEdgesBackWard(dbgK, pairRI[j], int(pos), strand, edgesArr, nodesArr, cf.Kmerlen, true)
		//func MappingReadToEdgesBackWard(dk DBGKmer, ri ReadInfo, rpos int, rstrand bool, edgesArr []DBGEdge, nodesArr []DBGNode, kmerlen int, correct bool) (errorNum int, ai AlignInfo) {
		//aB.Paths = Reverseuint32Arr(aB.Paths)
		//aB.Strands = GetReverseBoolArr(aB.Strands)
		// map the end partition of read sequence
		//errorNum2, aF := MappingReadToEdgesForWard(dbgK, pairRI[j], int(pos)+cf.Kmerlen, strand, edgesArr, nodesArr, cf.Kmerlen, true)
		//aB.ID, aF.ID = pairRI[j].ID, pairRI[j].ID
		//fmt.Printf("[paraMapNGSAndMerge] aB: %v\n\taF: %v\n", aB, aF)
		//if aB.EndPos < -1 || aF.EndPos < -1 {
		//	fmt.Printf("[paraMapNGSAndMerge] not found end boundary\n")
		//	break
		//}
	}
	fmt.Printf("[paraMapNGSAndCorrect] not found seed read number is:%d allNum:%d  percent:%f\n", notFoundSeedNum, allNum, float32(notFoundSeedNum)/float32(allNum))
	fmt.Printf("[paraMapNGSAndCorrect] high error mapping read number is:%d  short mapping number is:%d oneMoreLocationNum:%d\n", highErrorNum, shortMappingNum, oneMoreLocationNum)
}

/*func paraMapNGSAndCorrectPair(cs <-chan RIPool, wc chan<- RIPool, nodesArr []DBGNode, edgesArr []DBGEdge, cf CuckooFilterDBGKmer, winSize, kmerlen int) {
	var notFoundSeedNum, shortMappingNum, highErrorNum, allNum int
	MaxStepNum := 10
	ka := make([]KmerInfo, MaxStepNum*winSize)
	rSeq := make([]byte, 300)
	for {
		pairRI, ok := <-cs
		if !ok {
			var t [2]ReadInfo
			wc <- t
			break
		}

		if len(pairRI[0].Seq) < 140 || len(pairRI[1].Seq) < 140 {
			continue
		}

		var wr [2]ReadInfo
		allNum++
		// read1
		{
			ok := true
			// found kmer seed position in the DBG edges
			dbgK, pos, strand := LocateSeedKmerCF(cf, pairRI[0], winSize, edgesArr, MaxStepNum, ka, rSeq)
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
				//errorNum1, aB := MappingReadToEdgesBackWard(dbgK, pairRI[j], int(pos), strand, edgesArr, nodesArr, cf.Kmerlen, true)
				//func MappingReadToEdgesBackWard(dk DBGKmer, ri ReadInfo, rpos int, rstrand bool, edgesArr []DBGEdge, nodesArr []DBGNode, kmerlen int, correct bool) (errorNum int, ai AlignInfo) {
				//aB.Paths = Reverseuint32Arr(aB.Paths)
				//aB.Strands = GetReverseBoolArr(aB.Strands)
				// map the end partition of read sequence
				//errorNum2, aF := MappingReadToEdgesForWard(dbgK, pairRI[j], int(pos)+cf.Kmerlen, strand, edgesArr, nodesArr, cf.Kmerlen, true)
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
					var mR ReadInfo
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
			dbgK, pos, strand := LocateSeedKmerCF(cf, pairRI[1], winSize, edgesArr, MaxStepNum, ka, rSeq)
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
				//errorNum1, aB := MappingReadToEdgesBackWard(dbgK, pairRI[j], int(pos), strand, edgesArr, nodesArr, cf.Kmerlen, true)
				//func MappingReadToEdgesBackWard(dk DBGKmer, ri ReadInfo, rpos int, rstrand bool, edgesArr []DBGEdge, nodesArr []DBGNode, kmerlen int, correct bool) (errorNum int, ai AlignInfo) {
				//aB.Paths = Reverseuint32Arr(aB.Paths)
				//aB.Strands = GetReverseBoolArr(aB.Strands)
				// map the end partition of read sequence
				//errorNum2, aF := MappingReadToEdgesForWard(dbgK, pairRI[j], int(pos)+cf.Kmerlen, strand, edgesArr, nodesArr, cf.Kmerlen, true)
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
					var mR ReadInfo
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
}*/

func ReLocation(da []DBGKmer, ra DBGKmer, ri ReadInfo, edgesArr []DBGEdge, kmerlen int) DBGKmer {
	//seq1 := make([]byte, kmerlen)
	//seq2 := make([]byte, kmerlen)
	var seq1, seq2 []byte
	if ra.GetStrand() {
		seq1 = ri.Seq[ra.GetPos() : int(ra.GetPos())+kmerlen]
	} else {
		seq1 = GetReverseCompByteArr(ri.Seq[ra.GetPos() : int(ra.GetPos())+kmerlen])
	}
	count := 0
	//fmt.Fprintf(os.Stderr, "[ReLocation]seq1:%s\n", Transform2Char(ri.Seq[ra.Pos:int(ra.Pos)+kmerlen]))
	var dbgK DBGKmer
	for _, dk := range da {
		//fmt.Fprintf(os.Stderr, "[ReLocation]eID:%d len(e.ks):%d seq2:%s\n", dk.ID, len(edgesArr[dk.ID].Ks), Transform2Char(edgesArr[dk.ID].Ks[dk.Pos:int(dk.Pos)+kmerlen]))
		//fmt.Fprintf(os.Stderr, "[ReLocation]seq2:%s\n", Transform2Char(edgesArr[dk.ID].Ks[dk.Pos:int(dk.Pos)+kmerlen]))
		if dk.GetStrand() {
			seq2 = edgesArr[dk.GetID()].Ks[dk.GetPos() : int(dk.GetPos())+kmerlen]
		} else {
			seq2 = GetReverseCompByteArr(edgesArr[dk.GetID()].Ks[dk.GetPos() : int(dk.GetPos())+kmerlen])
		}
		if EqualByteArr(seq1, seq2) {
			count++
			dbgK = dk
		}
	}
	if count >= 2 {
		if da[0].GetID() == da[1].GetID() {
			if ra.GetStrand() == da[0].GetStrand() {
				dbgK = da[0]
			} else if ra.GetStrand() == da[1].GetStrand() {
				dbgK = da[1]
			} else {
				dbgK.SetID(0)
			}
		} else {
			fmt.Fprintf(os.Stderr, "[ReLocation]count:%d >1 DBGKmer arr[0]:%s arr[1]:%s\n", count, GetDkString(da[0]), GetDkString(da[1]))
		}
	} else if count == 0 {
		fmt.Fprintf(os.Stderr, "[ReLocation]count:%d ==0 DBGKmer arr:%v\n", count, da)
		dbgK.SetID(0)
	}
	return dbgK
}

func paraMapNGSAndMerge(cs <-chan RIPool, wc chan<- RIPool, nodesArr []DBGNode, edgesArr []DBGEdge, cf CuckooFilterDBGKmer, opt optionsPP, MergeSeqArrPool *sync.Pool) {
	var notFoundSeedNum, shortMappingNum, highErrorNum, allNum, needMergeNum, mergeNum, unMergeNum, shortMergeNum, oneMoreLocationNum int
	ka := make([]DBGKmerSeq, opt.WinSize)
	for i := range ka {
		ka[i].ks = make([]byte, opt.Kmer)
	}
	rSeq := make([]byte, opt.MaxNGSReadLen)
	var da [2][]DBGKmer
	da[0], da[1] = make([]DBGKmer, 2), make([]DBGKmer, 2)
	var ri ReadInfo
	for {
		riPool, ok := <-cs
		if !ok {
			var t RIPool
			wc <- t
			break
		}
		riPool.MergeSeqArr = MergeSeqArrPool.Get().([][]byte)
		mIdx := 0
		allNum += len(riPool.RIArr)

		for i := 0; i < len(riPool.RIArr)-1; i += 2 {
			if riPool.RIArr[i].ID != riPool.RIArr[i+1].ID {
				log.Fatalf("[paraMapNGSAndMerge]riPool.RIArr[%d].ID:%d riPool.RIArr[%d+1].ID:%d\n", i, riPool.RIArr[i].ID, i, riPool.RIArr[i+1].ID)
			}
			if len(riPool.RIArr[i].Seq)+len(riPool.RIArr[i+1].Seq) < 300 {
				continue
			}
			var ra [2]DBGKmer
			var cRI [2]ReadInfo
			var riArr [2]ReadMapInfo
			var errorNum, mappingNum int
			needMerge := true
			for j := i; j < i+2; j++ {
				// found kmer seed position in the DBG edges
				ri = riPool.RIArr[j]
				da[j-i], ra[j-i] = LocateSeedKmerCF(cf, ri, opt.WinSize, opt.Kmer, edgesArr, ka, rSeq, da[j-i])
				if len(da[j-i]) == 0 { // not found in the cuckoofilter
					//fmt.Printf("[paraMapNGSAndMerge] read ID: %v not found seed!!!\n", pairRI[j].ID)
					notFoundSeedNum++
					needMerge = false
					riPool.RIArr[j].Seq = nil
					continue
				} else if len(da[j-i]) > 1 {
					da[j-i][0] = ReLocation(da[j-i], ra[j-i], ri, edgesArr, opt.Kmer)
					if da[j-i][0].GetID() == 0 {
						oneMoreLocationNum++
						needMerge = false
						//fmt.Fprintf(os.Stderr, "[paraMapNGSAndMerge]read ID:%d ReLocation failed!!!\n", ri.ID)
						continue
					}
				}

				dk := da[j-i][0]
				//fmt.Printf("[paraMapNGSAndMerge]readID: %v, dbgK: %v, pos: %v, strand: %v\n", pairRI[0].ID, dbgK, pos, strand)
				// extend seed map to the edges
				errorNum, mappingNum, riArr[j-i], cRI[j-i] = MappingReadToEdges(dk, ra[j-i], ri, edgesArr, opt.Kmer)
				if mappingNum < len(ri.Seq)*7/10 {
					needMerge = false
					//fmt.Printf("[paraMapNGSAndMerge] mappingNum: %v < %v\n", mappingNum, len(pairRI[j].Seq)*8/10)
					riArr[j-i].PathSeqArr = nil
					shortMappingNum++
					break
				}

				if errorNum > len(ri.Seq)*1/10 {
					needMerge = false
					//fmt.Printf("[paraMapNGSAndMerge] errorNum: %v > %v\n", errorNum[j], (len(pairRI[j].Seq)-pos)*8/100)
					riArr[j-i].PathSeqArr = nil
					highErrorNum++
					break
				}
			}
			// check need merge path
			if !needMerge {
				continue
			}
			needMergeNum++
			// get link path sequence
			// program unknow how many times repeat selfcycle
			//var m ReadMapInfo
			m := MergePathSeqArr(riArr[0], riArr[1], edgesArr, opt.Kmer)
			m.ID = riPool.RIArr[i].ID
			if len(m.PathSeqArr) == 0 {
				//fmt.Printf("[paraMapNGSAndMerge] unMerge pair riArr[0]: %v\n\t\triArr[1]: %v\n", riArr[0], riArr[1])
				unMergeNum++
			} else if len(m.PathSeqArr) > 0 {
				//fmt.Printf("[paraMapNGSAndMerge] merge ReadMapInfo: %v\n", m)
				seq := GetPathSeq(m, edgesArr, opt.Kmer, riArr, riPool.MergeSeqArr[mIdx])
				if len(seq) < opt.MaxNGSReadLen {
					shortMergeNum++
					continue
				}
				mIdx++
				xl := int(ra[0].GetPos() + ra[1].GetPos())
				if len(seq)+xl > cap(seq) {
					tmp := make([]byte, len(seq)+xl)
					copy(tmp[ra[0].GetPos():], seq)
					seq = tmp
				} else {
					seq = seq[:len(seq)+xl]
					copy(seq[ra[0].GetPos():], seq[:len(seq)-xl])
				}
				copy(seq[:ra[0].GetPos()], riPool.RIArr[i].Seq[:ra[0].GetPos()])
				copy(seq[len(seq)-int(ra[1].GetPos()):], riPool.RIArr[i+1].Seq[:ra[1].GetPos()])
				ReverseCompByteArr(seq[len(seq)-int(ra[1].GetPos()):])
				var mR ReadInfo
				mR.Seq = seq
				mR.ID = m.ID
				mR.Anotition = append(mR.Anotition, fmt.Sprintf("SL:%d", len(mR.Seq))...)
				riPool.RIArr[i] = mR
				riPool.RIArr[i+1].Seq = nil
				mergeNum++
				/*for _, pathSeq := range m.PathSeqArr {
				mR.Anotition += fmt.Sprintf("%v-", pathSeq.ID)
				}*/
				//mR.Anotition = fmt.Sprintf("EN:%d len:%d", errorNum[0]+errorNum[1], len(mR.Seq))
				//fmt.Printf("[paraMapNGSAndMerge]len(mR.Seq): %v,  mR.Seq: %v\n", len(mR.Seq), mR.Seq)
			}
		}
		wc <- riPool
		/*ps1 := riArr[0].PathSeqArr[len(riArr[0].PathSeqArr)-1]
		ps2 := riArr[1].PathSeqArr[len(riArr[1].PathSeqArr)-1]
		e1, e2 := &edgesArr[ps1.ID], &edgesArr[ps2.ID]
		if e1.GetTwoEdgesCycleFlag() > 0 && e2.GetTwoEdgesCycleFlag() > 0 {
			if e1.ID == e2.ID {

			} else if IsTwoEdgesCycle(e1, e2) {

			}
		} else if e1.StartNID == e1.EndNID && e2.StartNID == e2.EndNID {
			if e1.ID == e2.ID {

			}
		} else {

		}*/
		/*if laste1.GetTwoEdgesCycleFlag() == 0 && laste1.StartNID != laste1.EndNID {
			m = MergePathSeqArr(riArr[0], riArr[1], edgesArr, nodesArr, opt.Kmer)
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
		}*/
	}
	fmt.Printf("[paraMapNGSAndMerge] not found seed read pair number is:%d oneMoreLocationNum:%d allNum:%d  percent:%f\n", notFoundSeedNum, oneMoreLocationNum, allNum, float32(notFoundSeedNum+oneMoreLocationNum)*2/float32(allNum))
	fmt.Printf("[paraMapNGSAndMerge] high error mapping read pair number is:%d short mapping number is:%d need merge number:%d merge number:%d unMerge number:%d shortMergeNum:%d\n", highErrorNum, shortMappingNum, needMergeNum, mergeNum, unMergeNum, shortMergeNum)
}

func paraMapNGSAndAddFreq(cs <-chan ReadInfo, nodesArr []DBGNode, edgesArr []DBGEdge, cf CuckooFilterDBGKmer, winSize, kmerlen int) {
	var notFoundSeedNum, allNum int
	var totalMappingLen int
	ka := make([]DBGKmerSeq, winSize)
	rSeq := make([]byte, 300)
	dbgKArr := make([]DBGKmer, 2)
	for {
		ri, ok := <-cs
		if !ok {
			break
		}

		allNum++

		// found kmer seed position in the DBG edges
		dbgKArr, _ = LocateSeedKmerCF(cf, ri, winSize, kmerlen, edgesArr, ka, rSeq, dbgKArr)
		if len(dbgKArr) == 0 { // not found in the cuckoofilter
			//fmt.Printf("[paraMapNGSAndMerge] read ID: %v not found seed!!!\n", pairRI[j].ID)
			notFoundSeedNum++
			continue
		}

		//mappingNum := MappingReadToEdgesAndAddFreq(dbgKArr[0], ri, int(pos), strand, edgesArr, nodesArr, cf.Kmerlen)
		//totalMappingLen += mappingNum
	}
	if allNum > 0 && allNum > notFoundSeedNum {
		fmt.Printf("[paraMapNGSAndAddFreq] not found seed read pair number is : %v,allNum: %v,  percent: %v, avgMapping seq len: %v\n", notFoundSeedNum, allNum, float32(notFoundSeedNum)/float32(allNum), totalMappingLen/(allNum-notFoundSeedNum))
	} else {
		fmt.Printf("[paraMapNGSAndAddFreq] not found seed read pair number is : %v,allNum: %v,  \n", notFoundSeedNum, allNum)
	}
	return
}

func writeCorrectReads(fn string, wc <-chan RIPool, rbytesPool, riArrPool, MergeSeqArrPool *sync.Pool, numCPU, bufSize int) {
	t0 := time.Now()
	fp, err := os.Create(fn)
	if err != nil {
		log.Fatalf("[writeCorrectReads] failed to create file:%s, err:%v\n", fn, err)
	}
	defer fp.Close()
	//cbrofp := cbrotli.NewWriter(fp, cbrotli.WriteroptionsPP{Quality: 1})
	zfp, err1 := zstd.NewWriter(fp, zstd.WithEncoderCRC(false), zstd.WithEncoderConcurrency(1), zstd.WithEncoderLevel(1))
	defer zfp.Close()
	if err1 != nil {
		log.Fatalf("[writeCorrectReads]open write file:%s err: %v\n", fn, err1)
	}
	//buffp := bufio.NewWriterSize(cbrofp, 1<<20) // 1<<24 == 2**24
	//buffp := bufio.NewWriter(cbrofp)
	//gzfp := gzip.NewWriter(fp)
	//defer gzfp.Close()
	var readNum, finishNum int
	for {
		riPool := <-wc
		if len(riPool.RIArr) == 0 {
			finishNum++
			if finishNum == numCPU {
				break
			} else {
				continue
			}
		}
		//seq := Transform2Letters(ri.Seq).String()
		for _, ri := range riPool.RIArr {
			if len(ri.Seq) == 0 {
				continue
			}
			Transform2Char2(ri.Seq)
			//str := (*string)(unsafe.Pointer(&seq))
			//s := fmt.Sprintf(">%v\t%s\n%s\n", ri.ID, ri.Anotition, *str)
			fmt.Fprintf(zfp, ">%v\t%s\n%s\n", ri.ID, ri.Anotition, ri.Seq)
		}
		readNum += len(riPool.RIArr)

		rbytesPool.Put(riPool.Cs)
		riArrPool.Put(riPool.RIArr)
		if len(riPool.MergeSeqArr) > 0 {
			MergeSeqArrPool.Put(riPool.MergeSeqArr)
		}
		//if len(wc) > bufSize-100 {
		//	fmt.Printf("[writeCorrectReads] write correct read delay, len(wc): %v\n", len(wc))
		//}
	}
	fmt.Printf("[writeCorrectReads] write correct read num: %v to file: %v, used time: %v\n", readNum, fn, time.Now().Sub(t0))
	return
}

func SplitWrite(wc <-chan [2]ReadInfo, wc1, wc2 chan<- ReadInfo, numCPU int) {
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

func paraCorrectReadsFile(fn string, concurrentNum int, nodesArr []DBGNode, edgesArr []DBGEdge, cf CuckooFilterDBGKmer, opt optionsPP) {
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
	go ReadZstdFile(fn, &rbytesPool, cs)
	go GetReadInfoBucket(fn, cs, &RIArrPool, RIPoolChan, true)

	wc := make(chan RIPool, concurrentNum)
	defer close(wc)
	//filterLowReadLen := true
	fmt.Printf("[paraCorrectReadsFile]concurrentNum:%d begin process files:%s\n", concurrentNum, fn)
	//go LoadNGSReadsPair(fn1, fn2, cs, opt.Kmer, false)
	for j := 0; j < concurrentNum; j++ {
		go paraMapNGSAndCorrect(RIPoolChan, wc, nodesArr, edgesArr, cf, opt)
	}
	// write function
	farr := strings.Split(fn, ".")
	wfn := strings.Join(farr[:len(farr)-2], ".") + ".Correct." + "fa.zst"
	writeCorrectReads(wfn, wc, &rbytesPool, &RIArrPool, nil, concurrentNum, bufSize)
	//fmt.Printf("[paraCorrectReadsFile] write correct reads num:%d to file:%s\n", writeNum, wfn)
}

func paraProcessReadsFile(fn string, concurrentNum int, nodesArr []DBGNode, edgesArr []DBGEdge, cf CuckooFilterDBGKmer, opt optionsPP) {
	fmt.Printf("[paraProcessReadsFile]begin process files  fn:%s\n", fn)
	pair := true
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
	MergeSeqArrPool := sync.Pool{New: func() interface{} {
		arr := make([][]byte, size/2)
		for i := range arr {
			arr[i] = make([]byte, opt.MaxNGSReadLen*2)
		}
		return arr
	}}

	RIPoolChan := make(chan RIPool, concurrentNum)
	go ReadZstdFile(fn, &rbytesPool, cs)
	go GetReadInfoBucket(fn, cs, &RIArrPool, RIPoolChan, pair)

	wc := make(chan RIPool, concurrentNum)
	defer close(wc)

	for j := 0; j < concurrentNum; j++ {
		go paraMapNGSAndMerge(RIPoolChan, wc, nodesArr, edgesArr, cf, opt, &MergeSeqArrPool)
	}
	// write function
	idx := strings.LastIndex(fn, "Correct")
	if idx < 0 {
		log.Fatalf("[paraProcessReadsFile] fn:%s, must used suffix *.Correct.[fasta|fa|fq|fastq].zst\n", fn)
	}
	wfn := fn[:idx] + "Merged.fa.zst"
	writeCorrectReads(wfn, wc, &rbytesPool, &RIArrPool, &MergeSeqArrPool, concurrentNum, bufSize)
	//time.Sleep(time.Second * 10)
	//fmt.Printf("[paraProcessReadsFile] write correct reads num: %d to file: %s\n", writeNum, wfn)
}

func paraMapReadsFile(fn string, concurrentNum int, nodesArr []DBGNode, edgesArr []DBGEdge, cf CuckooFilterDBGKmer, opt optionsPP) {
	bufSize := 20000
	cs := make(chan ReadInfo, bufSize)
	for j := 0; j < concurrentNum; j++ {
		go paraMapNGSAndAddFreq(cs, nodesArr, edgesArr, cf, opt.WinSize, opt.Kmer)
	}
	fmt.Printf("[paraMapReadsFile]begin process files  fn: %v\n", fn)
	//LoadNGSReads(fn, cs, opt.Kmer)
}

func MappingNGSAndCorrect(opt optionsPP, nodesArr []DBGNode, edgesArr []DBGEdge) {
	// construct cuckoofilter of DBG sample
	// use all edge seq, so set MaxNGSReadLen to MaxInt
	cfSize := GetCuckoofilterDBGSampleSize(edgesArr, opt.WinSize, opt.Kmer)
	fmt.Printf("[MappingNGSAndCorrect]cfSize:%d\n", cfSize)
	cf := MakeCuckooFilterDBGKmer(uint64(cfSize * 2))
	ConstructCFDBGMinimizers(&cf, edgesArr, opt.Kmer, opt.WinSize)
	fmt.Printf("[MappingNGSAndCorrect]construct Smaple of DBG edges cuckoofilter number is:%d\n", cf.Count)

	cfgInfo, err := ParseCfg(opt.CfgFn, opt.Merge, opt.Correct)
	if err != nil {
		log.Fatal("[ParseCfg] found err")
	}
	fmt.Printf("[MappingNGSAndCorrect] cfgInfo: %v\n", cfgInfo)

	go runtime.GC()
	runtime.GOMAXPROCS(opt.NumCPU)
	for _, lib := range cfgInfo.Libs {
		if lib.AsmFlag != AllState && lib.SeqProfile != 1 {
			continue
		}
		//MaxPairLen := lib.InsertSize + lib.InsertSD

		for i := 0; i < len(lib.FnName); i++ {
			paraCorrectReadsFile(lib.FnName[i], opt.NumCPU-1, nodesArr, edgesArr, cf, opt)
		}
	}
	//time.Sleep(time.Second)
}

// MappingNGSAndCorrect() function parallel Map NGS reads to the DBG edges, then merge path output long single read
func MappingNGSAndMerge(opt optionsPP, nodesArr []DBGNode, edgesArr []DBGEdge) {
	// construct cuckoofilter of DBG sample
	// use all edge seq, so set MaxNGSReadLen to MaxInt
	cfSize := GetCuckoofilterDBGSampleSize(edgesArr, opt.WinSize, opt.Kmer)
	fmt.Printf("[MappingNGSAndMerge] cfSize:%d\n", cfSize)
	cf := MakeCuckooFilterDBGKmer(uint64(cfSize * 2))
	fmt.Printf("[MappingNGSAndMerge]len(cf.Hash):%d\n", len(cf.Hash)*BucketSize)
	ConstructCFDBGMinimizers(&cf, edgesArr, opt.Kmer, opt.WinSize)
	fmt.Printf("[MappingNGSAndMerge]construct Smaple of DBG edges cuckoofilter number is:%d\n", cf.Count)

	cfgInfo, err := ParseCfg(opt.CfgFn, opt.Merge, opt.Correct)
	if err != nil {
		log.Fatal("[ParseCfg] found err")
	}
	fmt.Printf("[MappingNGSAndMerge] cfgInfo: %v\n", cfgInfo)

	go runtime.GC()
	runtime.GOMAXPROCS(opt.NumCPU)
	for _, lib := range cfgInfo.Libs {
		if lib.AsmFlag != AllState && lib.SeqProfile != 1 {
			continue
		}
		//MaxPairLen := lib.InsertSize + lib.InsertSD
		for i := 0; i < len(lib.FnName); i++ {
			paraProcessReadsFile(lib.FnName[i], opt.NumCPU-1, nodesArr, edgesArr, cf, opt)
		}
	}
}

/*func PrintCF(cf CuckooFilter) {
	for i, bt := range cf.Hash {
		for j := 0; j < BucketSize; j++ {
			if bt.Bkt[j].ID > 0 {
				fmt.Printf("[PrintCF]cf.Hash[%v].Bkt[%v]: %v, count: %v\n", i, j, bt.Bkt[j], bt.Bkt[j].GetCount())
			}
		}
	}
}*/

// MappingNGSAndCorrect() function parallel Map NGS reads to the DBG edges, then merge path output long single read
/*func MappingNGSAndAddFreq(opt optionsPP, nodesArr []DBGNode, edgesArr []DBGEdge) {
	// construct cuckoofilter of DBG sample
	cfSize := GetCuckoofilterDBGSampleSize(edgesArr, int64(opt.WinSize), int64(opt.Kmer))
	fmt.Printf("[MappingNGSAndAddFreq] cfSize: %v\n", cfSize)
	cf := MakeCuckooFilter(uint64(cfSize*2), opt.Kmer)
	fmt.Printf("[MappingNGSAndAddFreq] cf.numItems: %v\n", cf.NumItems)
	count := ConstructCFDBGMinimizers(cf, edgesArr, opt.WinSize)
	fmt.Printf("[MappingNGSAndAddFreq]construct Smaple of DBG edges cuckoofilter number is : %v\n", count)

	cfgInfo, err := ParseCfg(opt.CfgFn, opt.Merge, opt.Correct)
	if err != nil {
		log.Fatal("[ParseCfg] found err")
	}
	fmt.Printf("[MappingNGSAndAddFreq] cfgInfo: %v\n", cfgInfo)

	runtime.GOMAXPROCS(opt.NumCPU)

	//PrintCF(cf)

	// test
	//totalNumT := 1
	//concurrentNum := 6
	for _, lib := range cfgInfo.Libs {
		if lib.AsmFlag != AllState && lib.SeqProfile != 1 {
			continue
		}
		for i := 0; i < len(lib.FnName); i++ {
			paraMapReadsFile(lib.FnName[i], opt.NumCPU, nodesArr, edgesArr, cf, opt)
		}
	}
}*/

func Correct(c cli.Command) {

	//opt := optionsPP{gOpt, 0, 0, 0, 0, true, true, 0}
	opt, suc := checkArgsPP(c)
	if suc == false {
		log.Fatalf("[Correct] check Arguments error, opt:%v\n", opt)
	}
	//debug.PrintStack()
	cfgInfo, err := ParseCfg(opt.CfgFn, opt.Merge, opt.Correct)
	if err != nil {
		log.Fatalf("[Correct] ParseCfg 'C': %v err :%v\n", opt.CfgFn, err)
	}
	fmt.Printf("[Correct] opt:%+v\n\tcfgInfo:%+v\n", opt, cfgInfo)

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
	uniqKmerfn := opt.Prefix + ".uniqkmerseq.zst"
	if _, err := os.Stat(uniqKmerfn); err != nil {
		CCF(c)
	}
	nodesfn := opt.Prefix + ".nodes.Arr"
	if _, err := os.Stat(nodesfn); err != nil {
		CDBG(c)
	}
	// smfy DBG
	// read nodes file and transform to array mode for more quckly access
	DBGStatfn := opt.Prefix + ".DBG.stat"
	nodesSize, edgesSize := DBGStatReader(DBGStatfn)
	//nodesSize := len(nodeMap)
	fmt.Printf("[Correct] len(nodesArr):%d length of edge array:%d\n", nodesSize, edgesSize)
	nodesArr := make([]DBGNode, nodesSize)
	fc := make(chan int, 1)
	defer close(fc)
	go NodesArrReader(nodesfn, nodesArr, opt.Kmer, fc)
	// read edges file
	edgesfn := opt.Prefix + ".edges.fa.zst"
	edgesArr := ReadEdgesFromFile(edgesfn, edgesSize, opt.Kmer)
	<-fc

	var copt optionsSF
	copt.ArgsOpt = opt.ArgsOpt
	copt.MaxNGSReadLen, copt.TipMaxLen, copt.WinSize, copt.MinMapFreq = opt.MaxNGSReadLen, opt.TipMaxLen, opt.WinSize, opt.MinKmerFreq
	//gfn := opt.Prefix + ".beforeSmfyDBG.dot"
	//GraphvizDBGArr(nodesArr, edgesArr, gfn)
	t0 := time.Now()
	//PrintTmpNodesArr(nodesArr, opt.Prefix+".beforeSmfy")
	//PrintTmpDBG(nodesArr, edgesArr, opt.Prefix+".beforeSmfy")
	//CheckInterConnectivity(edgesArr, nodesArr)
	for i := 0; i < 1; i++ {
		fmt.Printf("[Correct] %v cycle Smfy DBG....\n", i)
		SmfyDBG(nodesArr, edgesArr, copt)
	}
	//edgesArr = SetEdgeID(nodesArr, edgesArr)
	//nodesArr = SetNodeID(nodesArr, edgesArr)

	CheckDBGSelfCycle(nodesArr, edgesArr, copt.Kmer)
	AddNodeInfo2DBGEdgeArr(edgesArr, nodesArr)
	//PrintTmpDBG(nodesArr, edgesArr, opt.Prefix)
	CheckInterConnectivity(edgesArr, nodesArr)

	uniqueNum, bubbleRepeatNum, twoEdgeCycleNum, selfCycleNum, selfCycleSameComingNum, bubbleEdgeNum := SetDBGEdgesUniqueFlag(edgesArr, nodesArr, opt.Kmer)
	fmt.Printf("[Correct] unique edge number is : %v, bubbleRepeatNum: %v, twoEdgeCycleNum: %v,selfCycleNum: %v, selfCycleSameComingNum: %v,bubbleEdgeNum: %v\n", uniqueNum, bubbleRepeatNum, twoEdgeCycleNum, selfCycleNum, selfCycleSameComingNum, bubbleEdgeNum)
	fmt.Printf("[Correct] Smfy DBG used : %v\n", time.Now().Sub(t0))
	t0 = time.Now()
	//gfn1 := opt.Prefix + ".afterSmfyDBG.dot"
	//GraphvizDBGArr(nodesArr, edgesArr, gfn1)

	// Mapping NGS to DBG, Delete low freq edge and Correct
	//copt.MinMapFreq = 5
	// write edges fq
	//ppfqfn := opt.Prefix + ".edges.pp.smfy.fq.br"
	//StoreEdgesToFn(ppfqfn, edgesArr)
	runtime.GC()
	if opt.Merge == false {
		//MappingNGSAndAddFreq(opt, nodesArr, edgesArr)
		//fmt.Printf("[Correct] Mapping  NGS and add edge mapping freq info used : %v\n", time.Now().Sub(t0))
		//t0 = time.Now()
		//DeleteLowFreqEdgeFlag := true
		//go runtime.GC()
		//SmfyDBG(nodesArr, edgesArr, copt, DeleteLowFreqEdgeFlag)
		MappingNGSAndCorrect(opt, nodesArr, edgesArr)
		fmt.Printf("[Correct] Mapping  NGS and Correct base used : %v\n", time.Now().Sub(t0))
	} else {
		MappingNGSAndMerge(opt, nodesArr, edgesArr)
		fmt.Printf("[Correct] Mapping and Merge NGS reads used : %v\n", time.Now().Sub(t0))
	}

	/*graphfn := opt.Prefix + ".afterSmfy.dot"
	go GraphvizDBGArr(nodesArr, edgesArr, graphfn)*/
	//go runtime.GC()
}
