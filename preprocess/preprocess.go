package preprocess

import (
	"compress/gzip"
	"fmt"
	"io"
	"log"
	"math"
	"os"
	"reflect"
	"runtime"
	"strconv"
	"strings"

	"github.com/biogo/biogo/alphabet"
	"github.com/biogo/biogo/io/seqio/fastq"
	"github.com/biogo/biogo/seq/linear"
	"github.com/jwaldrip/odin/cli"
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
	tmp, err := strconv.Atoi(c.Flag("s").String())
	if err != nil {
		log.Fatalf("[checkArgs] argument 's': %v set error: %v\n", c.Flag("s"), err)
	}
	if tmp < 1024*1024 {
		log.Fatalf("the argument 's': %v must bigger than 1024 * 1024\n", c.Flag("s"))
	}
	opt.CFSize = int64(tmp)
	tmp, err = strconv.Atoi(c.Flag("MaxNGSReadLen").String())
	if err != nil {
		log.Fatalf("[checkArgs] argument 'MaxNGSReadLen': %v set error: %v\n", c.Flag("MaxNGSReadLen"), err)
	}
	opt.MaxNGSReadLen = tmp
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

func LoadNGSReads(fn1, fn2 string, cs chan<- [2]constructdbg.ReadInfo, kmerlen int) {
	fp1, err1 := os.Open(fn1)
	if err1 != nil {
		log.Fatalf("[paraLoadNGSReads] open file: %v failed..., err: %v\n", fn1, err1)
	}
	fp2, err2 := os.Open(fn2)
	if err2 != nil {
		log.Fatalf("[paraLoadNGSReads] open file: %v failed..., err: %v\n", fn2, err2)
	}

	gzfp1, err3 := gzip.NewReader(fp1)
	if err3 != nil {
		log.Fatalf("[paraLoadNGSReads]gz read file: %v failed..., err: %v\n", fn1, err3)
	}
	gzfp2, err4 := gzip.NewReader(fp2)
	if err4 != nil {
		log.Fatalf("[paraLoadNGSReads]gz read file: %v failed..., err: %v\n", fn2, err4)
	}
	// brfp, err := brotli.NewReader(fp, nil)
	defer gzfp1.Close()
	defer gzfp2.Close()
	defer fp1.Close()
	defer fp2.Close()
	// defer fp.Close()
	// buffp := bufio.NewReader(gzfp)
	fqfp1 := fastq.NewReader(gzfp1, linear.NewQSeq("", nil, alphabet.DNA, alphabet.Sanger))
	fqfp2 := fastq.NewReader(gzfp2, linear.NewQSeq("", nil, alphabet.DNA, alphabet.Sanger))
	s1, err1 := fqfp1.Read()
	s2, err2 := fqfp2.Read()
	var count int
	for err1 == nil && err2 == nil {
		var ri1, ri2 constructdbg.ReadInfo
		fq1 := s1.(*linear.QSeq)
		fq2 := s2.(*linear.QSeq)
		if len(fq1.Seq) < kmerlen+10 || len(fq2.Seq) < kmerlen+10 { // read length is short for found DBG paths
			continue
		}
		id1, err := strconv.Atoi(fq1.Name())
		if err != nil {
			log.Fatalf("[paraLoadNGSReads] load fn: '%v' file, read ID: %v not digits, please convert to digits...\n", fn1, fq1.Name())
		}
		id2, err := strconv.Atoi(fq2.Name())
		if err != nil {
			log.Fatalf("[paraLoadNGSReads] load fn: '%v' file, read ID: %v not digits, please convert to digits...\n", fn2, fq2.Name())
		}
		if id1 != id2 {
			log.Fatalf("[paraLoadNGSReads] read1 ID : %v != read2 ID: %v\n", id1, id2)
		}
		ri1.ID = int64(id1)
		ri2.ID = int64(id2)
		ri1.Seq = constructdbg.Transform2Unitig(fq1.Seq, false).Ks
		ri2.Seq = constructdbg.Transform2Unitig(fq2.Seq, false).Ks
		var pairRI [2]constructdbg.ReadInfo
		pairRI[0], pairRI[1] = ri1, ri2
		cs <- pairRI
		count++
		s1, err1 = fqfp1.Read()
		s2, err2 = fqfp2.Read()
	}
	if err1 != io.EOF {
		log.Fatalf("[LoadNGSReads] Failed to read file %v, err: %v\n", fn2, err1)
	}
	if err2 != io.EOF {
		log.Fatalf("[LoadNGSReads] Failed to read file %v, err: %v\n", fn2, err2)
	}

	close(cs)
}

type ReadMapInfo struct {
	ID           int64
	StartP, EndP int
	Path         []constructdbg.DBG_MAX_INT
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
	if len(m.Path) == 1 {
		seq = edgesArr[m.Path[0]].Utg.Ks[m.StartP+1 : m.EndP]
		if m.StartP > m.EndP {
			seq = constructdbg.GetReverseCompByteArr(seq)
		}
		return seq
	}
	var strand bool
	e0, e1 := edgesArr[m.Path[0]], edgesArr[m.Path[1]]
	nID := constructdbg.GetLinkNodeID(e0, e1)
	if nID == e0.EndNID {
		strand = constructdbg.PLUS
		seq = e0.Utg.Ks[m.StartP+1:]
	} else {
		strand = constructdbg.MINUS
		seq = constructdbg.GetReverseCompByteArr(e0.Utg.Ks[:m.StartP])
	}
	for _, eID := range m.Path[1 : len(m.Path)-1] {
		e := edgesArr[eID]
		if e.EndNID == nID {
			strand = !strand
			nID = e.StartNID
		} else {
			nID = e.EndNID
		}
		var es []byte
		if strand == constructdbg.PLUS {
			es = e.Utg.Ks[kmerlen-1:]
		} else {
			es = constructdbg.GetReverseCompByteArr(e.Utg.Ks[:len(e.Utg.Ks)-(kmerlen-1)])
		}
		seq = append(seq, es...)
	}
	e := edgesArr[m.Path[len(m.Path)-1]]
	if e.EndNID == nID {
		strand = !strand
		nID = e.StartNID
	} else {
		nID = e.EndNID
	}
	var es []byte
	if strand == constructdbg.PLUS {
		es = e.Utg.Ks[kmerlen-1 : m.EndP]
	} else {
		es = constructdbg.GetReverseCompByteArr(e.Utg.Ks[m.EndP+1 : len(e.Utg.Ks)-(kmerlen-1)])
	}
	seq = append(seq, es...)

	return
}

func paraMapNGSAndMerge(cs <-chan [2]constructdbg.ReadInfo, wc chan<- constructdbg.ReadInfo, nodesArr []constructdbg.DBGNode, edgesArr []constructdbg.DBGEdge, cf constructdbg.CuckooFilter, winSize, MaxPairLen int) {
	var notFoundSeedNum, notPerfectNum int
	for {
		var mR constructdbg.ReadInfo
		pairRI, ok := <-cs
		if !ok {
			wc <- mR
			break
		}

		var riArr [2]ReadMapInfo
		for j := 0; j < 2; j++ {
			// found kmer seed position in the DBG edges
			dbgK, pos, strand := constructdbg.LocateSeedKmerCF(cf, pairRI[j], winSize, edgesArr)
			if dbgK.GetCount() == 0 { // not found in the cuckoofilter
				notFoundSeedNum++
				break
			}
			// extend seed map to the edges
			// map the start partition of read sequence
			errorNum1, ar := constructdbg.MappingReadToEdgesBackWard(dbgK, pairRI[j], int(pos), strand, edgesArr, nodesArr, cf, false)
			ar.Paths = constructdbg.ReverseDBG_MAX_INTArr(ar.Paths)
			// map the end partition of read sequence
			errorNum2, al := constructdbg.MappingReadToEdgesForWard(dbgK, pairRI[j], int(pos)+cf.Kmerlen, strand, edgesArr, nodesArr, cf, false)
			if errorNum1+errorNum2 > len(pairRI[j].Seq)*3/100 {
				notPerfectNum++
				fmt.Printf("[paraMapNGSAndMerge] errorNum: %v > %v\n", errorNum1+errorNum2, len(pairRI[j].Seq)*3/100)
			}
			if len(ar.Paths) > 0 && len(al.Paths) > 0 && ar.Paths[len(ar.Paths)-1] == al.Paths[0] {
				riArr[j].Path = append(ar.Paths, al.Paths[1:]...)
				riArr[j].ID = pairRI[j].ID
				riArr[j].StartP = ar.EndPos
				riArr[j].EndP = al.EndPos
			} else {
				log.Fatalf("[paraMapNGS2DBG] ar: %v and al: %v not consis\n", ar, al)
			}
		}

		// find link Path
		var ri constructdbg.ReadInfo
		ri.ID = pairRI[0].ID
		var m ReadMapInfo
		if len(riArr[0].Path) == 1 && len(riArr[1].Path) == 1 && riArr[0].Path[0] == riArr[1].Path[0] {
			m.StartP = riArr[0].StartP
			m.Path = append(m.Path, riArr[0].Path[0])
			m.EndP = riArr[1].StartP
		} else {
			extLen := cf.Kmerlen - (len(pairRI[0].Seq) + len(pairRI[1].Seq) - MaxPairLen) + 50
			var nID constructdbg.DBG_MAX_INT
			pl := len(riArr[0].Path)
			if pl > 1 {
				e1, e2 := edgesArr[riArr[0].Path[pl-2]], edgesArr[riArr[0].Path[pl-1]]
				nID = constructdbg.GetInterNodeID(e1, e2)
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
			if extLen < 0 {
				log.Fatalf("[paraMapNGS2DBG] extLen: %v < 0, riArr: %v\n", extLen, riArr)
			}
			epArr := GetExtendPathArr(riArr[0].Path, nID, edgesArr, nodesArr, cf.Kmerlen, extLen)
			riArr[1].Path = constructdbg.GetReverseDBG_MAX_INTArr(riArr[1].Path)
			var linkPathArr [][]constructdbg.DBG_MAX_INT
			for _, ep := range epArr {
				idx := constructdbg.IndexEID(riArr[1].Path, ep[len(ep)-1])
				if idx > len(ep)-1 {
					continue
				}
				if reflect.DeepEqual(ep[len(ep)-1-idx:], riArr[1].Path[:idx+1]) {
					ep = append(ep, riArr[1].Path[idx+1:]...)
					linkPathArr = append(linkPathArr, ep)
				}
			}
			if len(linkPathArr) == 1 {
				m.StartP = riArr[0].StartP
				m.Path = append(m.Path, linkPathArr[0]...)
				m.EndP = riArr[1].StartP
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
		}
	}
	fmt.Printf("[paraMapNGS2DBG] not found seed read pair number is : %v\n", notFoundSeedNum)
	fmt.Printf("[paraMapNGS2DBG] too more error mapping read pair number is : %v\n", notPerfectNum)
}

func writeCorrectReads(gzwfn string, wc <-chan constructdbg.ReadInfo, numCPU int) (readNum int) {
	fp, err := os.Create(gzwfn)
	if err != nil {
		log.Fatalf("[writeCorrectReads] failed to create file: %s, err: %v\n", gzwfn, err)
	}
	gzfp := gzip.NewWriter(fp)
	defer gzfp.Close()
	defer fp.Close()
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
		s := fmt.Sprintf(">%v\tpath:%s\n%s\n", ri.ID, ri.Anotition, string(ri.Seq))
		gzfp.Write([]byte(s))
		readNum++
	}

	if err := gzfp.Flush(); err != nil {
		log.Fatalf("[writeCorrectReads] failed to flush file: %s, err: %v\n", gzwfn, err)
	}

	return
}

// MappingNGSAndCorrect() function parallel Map NGS reads to the DBG edges, then merge path output long single read
func MappingNGSAndMerge(opt Options, nodesArr []constructdbg.DBGNode, edgesArr []constructdbg.DBGEdge) {
	// construct cuckoofilter of DBG sample
	// use all edge seq, so set MaxNGSReadLen to MaxInt
	cfSize := constructdbg.GetCuckoofilterDBGSampleSize(edgesArr, opt.WinSize, math.MaxInt64, opt.Kmer)
	fmt.Printf("[MappingNGSAndCorrect] cfSize: %v\n", cfSize)
	cf := constructdbg.MakeCuckooFilter(uint64(cfSize*5), opt.Kmer)
	fmt.Printf("[MappingNGSAndCorrect] cf.numItems: %v\n", cf.NumItems)
	count := constructdbg.ConstructCFDBGMinimizers(cf, edgesArr, opt.WinSize, math.MaxInt64)
	fmt.Printf("[MappingNGSAndCorrect]construct Smaple of DBG edges cuckoofilter number is : %v\n", count)

	cfgInfo, err := constructcf.ParseCfg(opt.CfgFn)
	if err != nil {
		log.Fatal("[constructcf.ParseCfg] found err")
	}
	fmt.Printf("[MappingNGSAndCorrect] cfgInfo: %v\n", cfgInfo)

	runtime.GOMAXPROCS(opt.NumCPU + 2)

	for _, lib := range cfgInfo.Libs {
		if lib.AsmFlag != constructcf.AllState && lib.SeqProfile != 1 {
			continue
		}
		MaxPairLen := lib.InsertSize + lib.InsertSD

		for i := 0; i < len(lib.FnName)-1; i += 2 {
			fn1, fn2 := lib.FnName[i], lib.FnName[i+1]
			idx1 := strings.LastIndex(fn1, "1")
			idx2 := strings.LastIndex(fn2, "2")
			if idx1 > 0 && idx2 > 0 && fn1[:idx1] == fn2[:idx2] {

			} else {
				log.Fatalf("[MappingNGSAndCorrect] fn1: %v, fn2: %v, must used suffix *[1|2].[fasta|fa|fq|fastq].gz\n", fn1, fn2)
			}

			bufSize := opt.NumCPU
			cs := make(chan [2]constructdbg.ReadInfo, bufSize)
			wc := make(chan constructdbg.ReadInfo, bufSize)
			go LoadNGSReads(fn1, fn2, cs, opt.Kmer)
			for j := 0; j < opt.NumCPU; j++ {
				go paraMapNGSAndMerge(cs, wc, nodesArr, edgesArr, cf, opt.WinSize, MaxPairLen)
			}
			// write function
			gzwfn := fn1[:idx1] + ".Correct.fa.gz"
			writeCorrectReads(gzwfn, wc, opt.NumCPU)
		}
	}

}

func Correct(c cli.Command) {
	constructcf.CCF(c)
	constructdbg.CDBG(c)
	gOpt, suc := utils.CheckGlobalArgs(c.Parent())
	if suc == false {
		log.Fatalf("[Smfy] check global Arguments error, opt: %v\n", gOpt)
	}
	opt := Options{gOpt, 0, 0, 0, 0}
	tmp, suc := checkArgs(c)
	if suc == false {
		log.Fatalf("[Smfy] check Arguments error, opt: %v\n", tmp)
	}
	opt.CFSize = tmp.CFSize
	opt.MaxNGSReadLen = tmp.MaxNGSReadLen
	opt.TipMaxLen = tmp.TipMaxLen
	opt.WinSize = tmp.WinSize
	cfgInfo, err := constructcf.ParseCfg(opt.CfgFn)
	if err != nil {
		log.Fatalf("[CCF] ParseCfg 'C': %v err :%v\n", opt.CfgFn, err)
	}
	fmt.Println(cfgInfo)

	// smfy DBG
	// read nodes file and transform to array mode for more quckly access
	nodesfn := opt.Prefix + ".nodes.mmap"
	nodeMap := constructdbg.NodeMapMmapReader(nodesfn)
	DBGStatfn := opt.Prefix + ".DBG.stat"
	nodesSize, edgesSize := constructdbg.DBGStatReader(DBGStatfn)
	//nodesSize := len(nodeMap)
	fmt.Printf("[Smfy] len(nodeMap): %v, length of edge array: %v\n", nodesSize, edgesSize)
	// read edges file
	edgesfn := opt.Prefix + ".edges"
	edgesArr := constructdbg.ReadEdgesFromFile(edgesfn, edgesSize)
	gfn1 := opt.Prefix + ".beforeSmfyDBG.dot"
	constructdbg.GraphvizDBG(nodeMap, edgesArr, gfn1)

	nodesArr := make([]constructdbg.DBGNode, nodesSize)
	constructdbg.NodeMap2NodeArr(nodeMap, nodesArr)
	nodeMap = nil // nodeMap any more used
	var copt constructdbg.Options
	copt.CfgFn, copt.Kmer, copt.MaxNGSReadLen, copt.NumCPU, copt.Prefix, copt.TipMaxLen, copt.WinSize = opt.CfgFn, opt.Kmer, opt.MaxNGSReadLen, opt.NumCPU, opt.Prefix, opt.TipMaxLen, opt.WinSize
	constructdbg.SmfyDBG(nodesArr, edgesArr, copt)

	// Mapping NGS to DBG and Correct
	MappingNGSAndMerge(opt, nodesArr, edgesArr)
}
