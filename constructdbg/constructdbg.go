package constructdbg

import (
	"bufio"
	"bytes"
	"encoding/binary"
	"encoding/gob"
	"fmt"
	"io"
	"log"
	"math"
	"os"
	"reflect"
	"runtime"
	"runtime/pprof"
	"sort"
	"strconv"
	"strings"
	"sync"
	"time"

	"github.com/mudesheng/ga/bnt"
	//"github.com/mudesheng/ga/cbrotli"
	"github.com/google/brotli/go/cbrotli"
	"github.com/mudesheng/ga/constructcf"
	"github.com/mudesheng/ga/cuckoofilter"
	"github.com/mudesheng/ga/utils"
	// "time"

	"github.com/awalterschulze/gographviz"
	"github.com/biogo/biogo/alphabet"
	"github.com/biogo/biogo/io/seqio/fasta"
	"github.com/biogo/biogo/io/seqio/fastq"
	"github.com/biogo/biogo/seq/linear"
	"github.com/jwaldrip/odin/cli"
)

type DBG_MAX_INT uint32 // use marker DBG node ID and DBG edge

type DBGNode struct {
	ID              DBG_MAX_INT
	EdgeIDIncoming  [bnt.BaseTypeNum]DBG_MAX_INT // the ID of EdgeID inclming
	EdgeIDOutcoming [bnt.BaseTypeNum]DBG_MAX_INT
	Seq             []uint64 // kmer node seq, use 2bit schema
	SubGID          uint8    // Sub graph ID, most allow 2**8 -1 subgraphs
	Flag            uint8    // from low~high, 1:Process, 2:
}

const (
	DIN             uint8 = 1 // note DBGNode and DBGEdge relationship
	DOUT            uint8 = 2 // note DBGNode and DBGEdge relationship
	PLUS, MINUS           = true, false
	NODEMAP_KEY_LEN       = 7
)

type NodeInfoByEdge struct {
	n1, n2 DBGNode
	c1, c2 uint8
	i1, i2 uint
	Flag   bool
	edgeID DBG_MAX_INT
}

func (n *DBGNode) GetProcessFlag() uint8 {
	return n.Flag & 0x1
}

func (n *DBGNode) SetProcessFlag() {
	// n.Flag = n.Flag & 0xFE
	n.Flag |= 0x1
}

func (n *DBGNode) ResetProcessFlag() {
	n.Flag &^= 0x1
}

func (n *DBGNode) GetDeleteFlag() uint8 {
	return n.Flag & 0x2
}

func (n *DBGNode) SetDeleteFlag() {
	// n.Flag = n.Flag & 0xFE
	n.Flag |= 0x2
}

type Unitig struct {
	Ks []byte
	Kq []uint8
}

type Path struct {
	IDArr  []DBG_MAX_INT
	AltArr []DBG_MAX_INT
	Freq   int
}

type EdgeFreq struct {
	ID   DBG_MAX_INT
	Freq uint8
}

type PathFreq struct {
	IDArr   []DBG_MAX_INT
	FreqArr []uint8
}

type Path1 struct {
	IDArr []DBG_MAX_INT
	//AltArr []DBG_MAX_INT
	Freq int
}

const SeedLen = 21
const SeqPositionBitNum = 64 - 1 - SeedLen*bnt.NumBitsInBase
const MaxAllowSeqLen = (1 << SeqPositionBitNum) - 1
const KmerMask = (1 << (SeqPositionBitNum + 1)) - 1
const PositionMask = (1 << SeqPositionBitNum) - 1
const KmerReSetPosition = (((1 << (SeedLen * bnt.NumBitsInBase)) - 1) << (SeqPositionBitNum + 1)) | 0x1
const KmerResetStrand = ((1 << (SeedLen*bnt.NumBitsInBase + SeqPositionBitNum)) - 1) << 1

type SeedInfo uint64

/*type KmerInfo struct {
	Kmer uint64 // kmerBase|Position|strand [SeedLen*bnt.NumBitsInBase:SeqPositionBitNum:1] just store 2-bits base, allow max 32 kmer base
	ID   uint32 // edgeID or reads ID
	//Info uint32 // the first high 31-bit for Position,and last bit for strand info
}*/

func (k SeedInfo) GetKmer() uint64 {
	return uint64(k >> (SeqPositionBitNum + 1))
}
func (k *SeedInfo) SetKmer(kmer uint64) {
	//fmt.Printf("[SetKmer]kmer:%b\n", kmer)
	*k = (*k << (SeqPositionBitNum + 1)) | (*k & KmerMask)
	//fmt.Printf("[SetKmer]k.Kmer:%b\n", k.GetKmer())
}

func (k SeedInfo) GetPosition() uint64 {
	return uint64((k >> 1) & PositionMask)
}
func (k *SeedInfo) SetPosition(p uint64) {
	if p >= (1 << SeqPositionBitNum) {
		log.Fatalf("[SetPosition]Pos:%d must < %d\n", p, 1<<SeqPositionBitNum)
	}
	//fmt.Printf("[SetPosition]Pos:%b\n", p)
	*k = SeedInfo((p << 1) | (uint64(*k) & KmerReSetPosition))
	//fmt.Printf("[SetPosition]k.Kmer:%b\n", k.Kmer)
}

func (k SeedInfo) GetStrand() bool {
	if (k & 0x1) > 0 {
		return true
	} else {
		return false
	}
}

func (k *SeedInfo) SetStrand(s bool) {
	if s == PLUS {
		*k |= 0x1
	} else {
		//k.Kmer = k.Kmer & KmerResetStrand
		*k >>= 1
		*k <<= 1
	}
}

type SeedInfoArr []SeedInfo

func (arr SeedInfoArr) Len() int {
	return len(arr)
}

func (arr SeedInfoArr) Less(i, j int) bool {
	return arr[i].GetKmer() < arr[j].GetKmer()
}
func (arr SeedInfoArr) Swap(i, j int) {
	arr[i], arr[j] = arr[j], arr[i]
}

type DBGEdge struct {
	ID                DBG_MAX_INT
	StartNID          DBG_MAX_INT                  // start node ID
	EndNID            DBG_MAX_INT                  // end node ID
	EdgeIDIncoming    [bnt.BaseTypeNum]DBG_MAX_INT // incoming edgeIDArr
	EdgeIDOutcoming   [bnt.BaseTypeNum]DBG_MAX_INT // outcoming edgeIDArr
	CovD              uint16                       // Coverage Depth or number of  low freq minimers
	Flag              uint8                        //
	Utg               Unitig
	NGSPathArr        [][]EdgeFreq
	PathMat, PathMat1 []Path // read Path matrix
	GFAPath           []*[][][]uint32
	SeedInfoArr       []SeedInfo // store seq seed array
}

/*
func EqualDBG_MAX_INT(path1, path2 []DBG_MAX_INT) bool {
	if len(path1) != len(path2) {
		return false
	}

	for i, v := range path1 {
		if v != path2[i] {
			return false
		}
	}
	return true
} */

func (e *DBGEdge) InsertPathToEdge(path []DBG_MAX_INT, freq int) {

	// check have been added
	added := false
	for i, v := range e.PathMat {
		rp := GetReverseDBG_MAX_INTArr(path)
		if reflect.DeepEqual(v.IDArr, path) || reflect.DeepEqual(v.IDArr, rp) {
			e.PathMat[i].Freq += freq
			added = true
			break
		}
	}
	if added == false {
		var np Path
		np.IDArr = make([]DBG_MAX_INT, len(path))
		copy(np.IDArr, path)
		np.Freq = freq
		e.PathMat = append(e.PathMat, np)
	}
}

/*func InsertPathToEdge(pathMat []Path, path []DBG_MAX_INT, freq int) []Path {

	// check have been added
	added := false
	for i, v := range pathMat {
		if EqualDBG_MAX_INT(v.IDArr, path) {
			pathMat[i].Freq += freq
			added = true
			break
		}
	}
	if added == false {
		var np Path
		np.IDArr = make([]DBG_MAX_INT, len(path))
		copy(np.IDArr, path)
		np.Freq = freq
		pathMat = append(pathMat, np)
	}

	return pathMat
}*/

func (e *DBGEdge) GetProcessFlag() uint8 {
	return e.Flag & 0x1
}

func (e *DBGEdge) GetDeleteFlag() uint8 {
	return e.Flag & 0x2
}

func (e *DBGEdge) GetUniqueFlag() uint8 {
	return e.Flag & 0x4
}

func (e *DBGEdge) GetBubbleRepeatFlag() uint8 {
	return e.Flag & 0x8
}

func (e *DBGEdge) GetTwoEdgesCycleFlag() uint8 {
	return e.Flag & 0x10
}

func (e *DBGEdge) GetMergedFlag() uint8 {
	return e.Flag & 0x20
}

func (e *DBGEdge) GetBubbleFlag() uint8 {
	return e.Flag & 0x40
}

func (e *DBGEdge) GetLowFreqMiniMersFlag() uint8 {
	return e.Flag & 0x80
}

func (e *DBGEdge) SetBubbleFlag() {
	e.Flag = e.Flag | 0x40
}

func (e *DBGEdge) SetBubbleRepeatFlag() {
	e.Flag = e.Flag | 0x8
}

func (e *DBGEdge) SetTwoEdgesCycleFlag() {
	e.Flag = e.Flag | 0x10
}

func (e *DBGEdge) SetMergedFlag() {
	e.Flag = e.Flag | 0x20
}
func (e *DBGEdge) SetLowFreqMiniMersFlag() {
	e.Flag = e.Flag | 0x80
}

func (e *DBGEdge) ResetBubbleFlag() {
	e.Flag = e.Flag & (0xFF - 0x40)
}

func (e *DBGEdge) ResetTwoEdgesCycleFlag() {
	e.Flag = e.Flag & (0xFF - 0x10)
}

func (e *DBGEdge) ResetMergedFlag() {
	e.Flag = e.Flag & (0xFF - 0x20)
}
func (e *DBGEdge) ResetLowFreqMiniMersFlag() {
	e.Flag = e.Flag & (0xFF - 0x80)
}

func (e *DBGEdge) SetUniqueFlag() {
	e.Flag = e.Flag | 0x4
}

func (e *DBGEdge) ResetUniqueFlag() {
	e.Flag = e.Flag & (0xFF - 0x4)
}

func (e *DBGEdge) ResetBubbleRepeatFlag() {
	e.Flag = e.Flag & (0xFF - 0x8)
}

func (e *DBGEdge) SetProcessFlag() {
	e.Flag = e.Flag | 0x1
}

func (e *DBGEdge) ResetProcessFlag() {
	e.Flag = e.Flag & (0xFF - 0x1)
}

func (e *DBGEdge) SetDeleteFlag() {
	e.Flag = e.Flag | 0x2
}

func (e *DBGEdge) ResetDeleteFlag() {
	e.Flag = e.Flag & (0xFF - 0x2)
}

func (e *DBGEdge) GetSeqLen() int {
	return len(e.Utg.Ks)
}

//var MIN_KMER_COUNT uint16 = 3
var BACKWARD uint8 = 0
var FORWARD uint8 = 1

//var Kmerlen int = -1

func ReverseDBG_MAX_INTArr(path []DBG_MAX_INT) []DBG_MAX_INT {
	al := len(path)
	dl := al / 2
	for i := 0; i < dl; i++ {
		path[i], path[al-1-i] = path[al-1-i], path[i]
	}

	return path
}

func GetKmerRecord(fncs <-chan []byte, cs chan<- constructcf.KmerBntBucket, bufSize, kmerlen int, kmerNumC chan<- int) {
	var readNum int
	buf := make([]byte, 0, bufSize+(1<<10))
	pos := 0
	nb, ok := <-fncs
	if !ok {
		//var b []byte
		//kmerChan <- b
		fmt.Printf("[GetKmerRecord] not found any data\n")
		//return readNum
	}
	var buck constructcf.KmerBntBucket
	size := 10000
	count := 0
	kM := make([]byte, kmerlen*size)
	buf = append(buf, nb...)
	//fmt.Printf("[GetKmerRecord] nb num is: %d, buf len: %d\n", len(nb), len(buf))
	for {
		//b := make([]byte, kmerlen)
		if len(buf)-pos < kmerlen {
			nb, _ := <-fncs
			if len(nb) == 0 {
				break
			}
			copy(buf[0:], buf[pos:len(buf)])
			buf = buf[:len(buf)-pos]
			pos = 0
			buf = append(buf, nb...)
			//fmt.Printf("[GetKmerRecord] nb num is: %d, buf len: %d\n", len(nb), len(buf))
		}

		if len(buf)-pos < kmerlen {
			continue
		}
		if count >= size {
			kM = make([]byte, kmerlen*size)
			count = 0
		}
		b := kM[count*kmerlen : (count+1)*kmerlen]
		copy(b, buf[pos:pos+kmerlen])
		pos += kmerlen
		var rb constructcf.KmerBnt
		rb.Seq = b
		rb.Len = kmerlen
		if buck.Count >= constructcf.ReadSeqSize {
			cs <- buck
			var nb constructcf.KmerBntBucket
			buck = nb
		}
		buck.KmerBntBuf[buck.Count] = rb
		buck.Count++
		readNum++
		count++
	}
	if pos%kmerlen != 0 {
		log.Fatalf("[GetKmerRecord] process kmer file pos:%v err\n", pos)
	}

	if buck.Count > 0 {
		cs <- buck
	}
	fmt.Printf("[GetKmerRecord] read kmer Record num is: %v\n", readNum)
	//time.Sleep(time.Second * 10)
	//close(kmerChan)
	kmerNumC <- readNum
	//return readNum
}

func ParaReadUniqKmer(prefix string, cs chan constructcf.KmerBntBucket, kmerlen int) {
	bufSize := (1 << 20)
	//size := 40000
	//kmerChan := make(chan []byte, size)
	kmerNumC := make(chan int, 1)
	for i := 0; i < constructcf.UniqueKmerFileNum; i++ {
		brfn := prefix + "." + strconv.Itoa(i) + ".uniqkmerseq.br"
		fncs := make(chan []byte, 10)
		go constructcf.ReadBrFile(brfn, fncs, bufSize)
		go GetKmerRecord(fncs, cs, bufSize, kmerlen, kmerNumC)
	}

	kmerNum := 0
	for i := 0; i < constructcf.UniqueKmerFileNum; i++ {
		kmerNum += <-kmerNumC
	}
	fmt.Printf("[readUniqKmer] total read kmer number is : %d\n", kmerNum)
	// send read finish signal
	//time.Sleep(time.Second * 10)
	close(cs)
}

/*func UniDirectExtend(kb, rkb constructcf.KmerBnt, cf cuckoofilter.CuckooFilter, min_kmer_count uint16, direction uint8) (ec int) {
	kmerlen := cf.Kmerlen
	//var nBnt constructcf.ReadBnt
	//nBnt.Seq = make([]byte, kmerlen)
	if direction == FORWARD {
		//copy(nBnt.Seq[:kmerlen-1], nb.Seq)
		for i := 0; i < bnt.BaseTypeNum; i++ {
			b := uint64(i)
			//nBnt.Seq[kmerlen-1] = b
			ks := constructcf.GetNextKmer(kb, b, kmerlen)
			rs := constructcf.GetPreviousKmer(rkb, uint64(bnt.BntRev[b]), kmerlen)
			if ks.BiggerThan(rs) {
				ks = rs
			}
			if count := cf.GetCountAllowZero(ks.Seq); count >= min_kmer_count {
				ec++
			}
		}
	} else { // direction == BACKWARD
		//copy(nBnt.Seq[1:], nb.Seq)
		for i := 0; i < bnt.BaseTypeNum; i++ {
			b := uint64(i)
			//nBnt.Seq[0] = b
			ks := constructcf.GetPreviousKmer(kb, b, kmerlen)
			rs := constructcf.GetNextKmer(rkb, uint64(bnt.BntRev[b]), kmerlen)
			if ks.BiggerThan(rs) {
				ks = rs
			}
			if count := cf.GetCountAllowZero(ks.Seq); count >= min_kmer_count {
				ec++
			}
		}
	}

	return
}*/

func ExtendNodeKmer(nkb, rkb constructcf.KmerBnt, cf cuckoofilter.CuckooFilter, min_kmer_count uint16, m []byte) (nd DBGNode) {
	kmerlen := cf.Kmerlen
	var kb2, rb2 []byte
	if len(m) < 2*kmerlen {
		log.Fatalf("[ExtendNodeKmer] len(m):%d < 2 *kmerlen:%d\n", len(m), 2*kmerlen)
	}
	kb2, rb2 = m[:kmerlen:kmerlen], m[kmerlen:2*kmerlen]
	//rb2 = make([]byte, kmerlen)
	//var nBnt constructcf.ReadBnt
	//nBnt.Seq = make([]byte, kmerlen+1)
	//nd.Seq = constructcf.GetReadBntKmer(nkb.Seq, 0, kmerlen-1)
	//copy(nBnt.Seq[1:], nodeBnt.Seq)
	//rkb := constructcf.ReverseComplet(nkb)
	//nd.Seq = nkb.Seq
	//fmt.Printf("[ExtendNodeKmer]nkb: %v\n\trkb: %v\n", nkb, rkb)
	for i := 0; i < bnt.BaseTypeNum; i++ {
		bi := byte(i)
		//nBnt.Seq[0] = bi
		{
			//ks := constructcf.GetPreviousKmer(nkb, bi, kmerlen)
			//rs := constructcf.GetNextKmer(rkb, uint64(bnt.BntRev[bi]), kmerlen)
			//kb2 = constructcf.NoAllocGetPreviousKmer(nkb, kb2, bi, kmerlen)
			//rb2 = constructcf.NoAllocGetNextKmer(rkb, rb2, uint64(bnt.BntRev[bi]), kmerlen)
			kb2[0] = bi
			copy(kb2[1:], nkb.Seq)
			copy(rb2[:kmerlen-1], rkb.Seq)
			rb2[kmerlen-1] = bnt.BntRev[bi]
			min := kb2
			if constructcf.BiggerThan(kb2, rb2) {
				min = rb2
			}
			count := cf.GetCountAllowZero(min)
			//fmt.Printf("[ExtendNodeKmer]min: %v, count: %v\n", min, count)
			if count >= min_kmer_count {
				//var nb constructcf.ReadBnt
				//nb.Seq = append(nb.Seq, nBnt.Seq[:kmerlen-1]...)
				//nb.Length = len(nb.Seq)
				/*if ec := UniDirectExtend(ks, rs, cf, min_kmer_count, BACKWARD); ec > 0 {
					if direction == FORWARD {
						baseBnt = uint8(bi)
						baseCount = count
					}
					//fmt.Printf("[ExtendNodeKmer] BACKWARD bi: %v, ks.Seq: %v\n", bi, ks.Seq)
				}*/
				nd.EdgeIDIncoming[bi] = 1
			}
		}

		{
			//nBnt.Seq[cf.Kmerlen] = bi
			//kb2 = constructcf.NoAllocGetNextKmer(nkb, kb2, bi, kmerlen)
			//rb2 = constructcf.NoAllocGetPreviousKmer(rkb, rb2, uint64(bnt.BntRev[bi]), kmerlen)
			copy(kb2, nkb.Seq)
			kb2[kmerlen-1] = bi
			rb2[0] = bnt.BntRev[bi]
			copy(rb2[1:], rkb.Seq)
			min := kb2
			if BiggerThan(kb2, rb2) {
				min = rb2
			}
			count := cf.GetCountAllowZero(min)
			//fmt.Printf("[ExtendNodeKmer]min: %v, count: %v\n", min, count)
			if count >= min_kmer_count {
				//fmt.Printf("[ExtendNodeKmer]count: %v\n", count)
				//var nb constructcf.ReadBnt
				//nb.Seq = append(nb.Seq, nBnt.Seq[2:]...)
				//nb.Length = len(nb.Seq)
				/*if ec := UniDirectExtend(ks, rs, cf, min_kmer_count, FORWARD); ec > 0 {
					if direction == BACKWARD {
						baseBnt = uint8(i)
						baseCount = count
					}
					//fmt.Printf("[ExtendNodeKmer] FORWARD bi: %v, ks.Seq: %v\n", bi, ks.Seq)
				}*/
				nd.EdgeIDOutcoming[bi] = 1
			}
		}
	}

	return
}

func ChangeEdgeIDComing(nd DBGNode) DBGNode {
	for i := 0; i < bnt.BaseTypeNum; i++ {
		if nd.EdgeIDIncoming[i] == math.MaxUint32 {
			nd.EdgeIDIncoming[i] = 1
		}
		if nd.EdgeIDOutcoming[i] == math.MaxUint32 {
			nd.EdgeIDOutcoming[i] = 1
		}
	}
	return nd
}

func GetMinDBGNode(nd DBGNode, kmerlen int) (minN DBGNode) {
	var nkb constructcf.KmerBnt
	nkb.Len = kmerlen - 1
	tmp := make([]uint64, len(nd.Seq))
	copy(tmp, nd.Seq)
	rs := constructcf.ReverseCompletBnt(nd.Seq, kmerlen-1)
	//fmt.Printf("[GetMinDBGNode] node: %v\nRC node: %v\n", nodeBnt, rnode)
	if constructcf.BiggerThanBnt(nd.Seq, rs) {
		minN.Seq = rs
		for i := 0; i < bnt.BaseTypeNum; i++ {
			minN.EdgeIDIncoming[i] = nd.EdgeIDOutcoming[bnt.BaseTypeNum-1-i]
			minN.EdgeIDOutcoming[i] = nd.EdgeIDIncoming[bnt.BaseTypeNum-1-i]
		}
	} else {
		minN.Seq = nd.Seq
		minN.EdgeIDIncoming = nd.EdgeIDIncoming
		minN.EdgeIDOutcoming = nd.EdgeIDOutcoming
	}
	return minN
}

func paraLookupComplexNode(cs chan constructcf.KmerBntBucket, wc chan DBGNode, cf cuckoofilter.CuckooFilter, MinKmerFreq uint16) {
	// var wrsb constructcf.ReadSeqBucket
	var rkb constructcf.KmerBnt
	rkb.Len = cf.Kmerlen
	rkb.Seq = make([]byte, rkb.Len)
	m := make([]byte, cf.Kmerlen*2)
	for {
		kbBucket, _ := <-cs
		if kbBucket.Count == 0 {
			var t DBGNode
			wc <- t
			break
		}
		// if rsb.Count < constructcf.ReadSeqSize {
		// 	fmt.Printf("rsb.ReadBuf length is : %d\n", len(rsb.ReadBuf))
		// }
		// if found kmer count is 1 , this kmer will be ignore, and skip this branch
		for j := 0; j < kbBucket.Count; j++ {
			kb := kbBucket.KmerBntBuf[j]
			rkb = constructcf.GetReverseComplet2(kb, rkb)
			//fmt.Printf("[paraLookupComplexNode] kb : %v\n", kb)
			{ // check fisrt node of kmer
				var nb, rb constructcf.KmerBnt
				nb.Len = cf.Kmerlen - 1
				rb.Len = cf.Kmerlen - 1
				//nkb.Seq = make([]byte, nkb.Len)
				nb.Seq = kb.Seq[:cf.Kmerlen-1]
				rb.Seq = rkb.Seq[1:]
				//nkb, _ := constructcf.DeleteLastBaseKmer(kb)
				//fmt.Printf("[paraLookupComplexNode] nkb : %v\n", nkb)
				//rkb := constructcf.ReverseComplet(nkb)
				//fmt.Printf("[paraLookupComplexNode] rkb : %v\n", rkb)
				//extkb := constructcf.ExtendKmerBnt2Byte(kb)
				//extnkb := constructcf.ExtendKmerBnt2Byte(nkb)
				//extrkb := constructcf.ExtendKmerBnt2Byte(rkb)
				//fmt.Printf("[paraLookupComplexNode]first kb: %v\n\tbase: %v,nkb: %v\n\trkb: %v\n", constructcf.ExtendKmerBnt2Byte(kb), base, constructcf.ExtendKmerBnt2Byte(nkb), constructcf.ExtendKmerBnt2Byte(rkb))
				//nodeBnt.Length = len(nodeBnt.Seq)
				nd := ExtendNodeKmer(nb, rb, cf, MinKmerFreq, m)
				var leftcount, rightcount int
				for i := 0; i < bnt.BaseTypeNum; i++ {
					if nd.EdgeIDIncoming[i] == 1 {
						leftcount++
					}
					if nd.EdgeIDOutcoming[i] == 1 {
						rightcount++
					}
				}
				if leftcount > 1 || rightcount > 1 {
					var wd DBGNode
					//wd.Seq = make([]uint64, len(nkb.Seq))
					if constructcf.BiggerThan(nb.Seq, rb.Seq) {
						nb = rb
						for i := 0; i < bnt.BaseTypeNum; i++ {
							wd.EdgeIDIncoming[i] = nd.EdgeIDOutcoming[bnt.BaseTypeNum-1-i]
							wd.EdgeIDOutcoming[i] = nd.EdgeIDIncoming[bnt.BaseTypeNum-1-i]
						}
					} else {
						wd.EdgeIDIncoming = nd.EdgeIDIncoming
						wd.EdgeIDOutcoming = nd.EdgeIDOutcoming
					}
					wd.Seq = constructcf.GetReadBntKmer(nb.Seq, 0, nb.Len)
					//copy(wd.Seq, nkb.Seq)
					//tn := GetMinDBGNode(nd, cf.Kmerlen)
					//fmt.Printf("[paraLookupComplexNode] nd: %v\n", wd)
					//fmt.Printf("[paraLookupComplexNode] node: %v\n", wd)
					wc <- wd
				}
			}

			{ // check second node of kmer
				var nb, rb constructcf.KmerBnt
				nb.Len = cf.Kmerlen - 1
				rb.Len = cf.Kmerlen - 1
				nb.Seq = kb.Seq[1:]
				rb.Seq = rkb.Seq[:cf.Kmerlen-1]
				//fmt.Printf("[paraLookupComplexNode]before nb: %v\n", nb)
				//rb := constructcf.GetReverseComplet(nkb)
				//nkb, _ := constructcf.DeleteFirstBaseKmer(kb)
				//fmt.Printf("[paraLookupComplexNode]second kb: %v\nbase: %v,nkb: %v\n", constructcf.ExtendKmerBnt2Byte(kb), base, constructcf.ExtendKmerBnt2Byte(nkb))
				//rkb := constructcf.ReverseComplet(nkb)
				nd := ExtendNodeKmer(nb, rb, cf, MinKmerFreq, m)
				var leftcount, rightcount int
				for i := 0; i < bnt.BaseTypeNum; i++ {
					if nd.EdgeIDIncoming[i] == 1 {
						leftcount++
					}
					if nd.EdgeIDOutcoming[i] == 1 {
						rightcount++
					}
				}
				if leftcount > 1 || rightcount > 1 {
					var wd DBGNode
					//wd.Seq = make([]uint64, len(nkb.Seq))
					if constructcf.BiggerThan(nb.Seq, rb.Seq) {
						nb = rb
						for i := 0; i < bnt.BaseTypeNum; i++ {
							wd.EdgeIDIncoming[i] = nd.EdgeIDOutcoming[bnt.BaseTypeNum-1-i]
							wd.EdgeIDOutcoming[i] = nd.EdgeIDIncoming[bnt.BaseTypeNum-1-i]
						}
					} else {
						wd.EdgeIDIncoming = nd.EdgeIDIncoming
						wd.EdgeIDOutcoming = nd.EdgeIDOutcoming
					}
					//fmt.Printf("[paraLookupComplexNode] nb: %v\n", nb)
					wd.Seq = constructcf.GetReadBntKmer(nb.Seq, 0, nb.Len)
					//copy(wd.Seq, nkb.Seq)
					wc <- wd
				}
			}
		}
	}
}

func GetNodeRecord(fncs <-chan []byte, nodeChan chan<- DBGNode, bufSize, NBntUint64Len int) int {
	var readNum int
	buf := make([]byte, 0, bufSize+(1<<12))
	pos := 0
	nb, ok := <-fncs
	if !ok {
		log.Fatalf("[GetNodeRecord] not found any data\n")
	}
	buf = append(buf, nb...)
	//fmt.Printf("[GetKmerRecord] nb num is: %d, buf len: %d\n", len(nb), len(buf))
	//var t DBGNodeInfo
	sizeDBGNode := NBntUint64Len*8 + bnt.BaseTypeNum*2*4
	//fmt.Printf("[GetNodeRecord]sizeDBGNodeInfo : %v\n", sizeDBGNode)
	for {
		//b := make([]byte, kmerlen)

		if len(buf)-pos < sizeDBGNode {
			nb, _ := <-fncs
			if len(nb) == 0 {
				break
			}
			copy(buf[0:], buf[pos:len(buf)])
			buf = buf[:len(buf)-pos]
			pos = 0
			buf = append(buf, nb...)
			//fmt.Printf("[GetKmerRecord] nb num is: %d, buf len: %d\n", len(nb), len(buf))
		}
		if len(buf)-pos < sizeDBGNode {
			continue
		}
		//bio := bytes.Buffer(buf[pos:])
		bnr := bytes.NewReader(buf[pos:])
		var node DBGNode
		node.Seq = make([]uint64, NBntUint64Len)
		err := binary.Read(bnr, binary.LittleEndian, node.Seq)
		if err != nil {
			log.Fatalf("[constructNodeMap] err: %v\n", err)
		}
		if err := binary.Read(bnr, binary.LittleEndian, &node.EdgeIDIncoming); err != nil {
			log.Fatalf("[constructNodeMap] err: %v\n", err)
		}
		if err := binary.Read(bnr, binary.LittleEndian, &node.EdgeIDOutcoming); err != nil {
			log.Fatalf("[constructNodeMap] err: %v\n", err)
		}
		pos += sizeDBGNode
		nodeChan <- node
		readNum++
	}
	if pos%sizeDBGNode != 0 {
		log.Fatalf("[GetNodeRecord] process node file pos:%v err\n", pos)
	}
	fmt.Printf("[GetNodeRecord] read node Record num is: %v\n", readNum)
	close(nodeChan)
	//time.Sleep(time.Second * 10)
	return readNum
}

func constructNodeMap(complexKmerbrfn string, nodeMap map[[NODEMAP_KEY_LEN]uint64]DBGNode, NBntUint64Len int) DBG_MAX_INT {
	nodeID := DBG_MAX_INT(2)
	bufSize := (1 << 20)
	size := 40000
	fncs := make(chan []byte, 10)
	nodeChan := make(chan DBGNode, size)

	go constructcf.ReadBrFile(complexKmerbrfn, fncs, bufSize)
	go GetNodeRecord(fncs, nodeChan, bufSize, NBntUint64Len)
	readnodeNum := 0
	for {
		node, _ := <-nodeChan
		if len(node.Seq) == 0 {
			break
		}
		// fmt.Fprintf(os.Stderr, "[constructNodeMap] node: %v\n", node)
		// n, err := ckfp.Read(b)
		// if n != NBntByteLen {
		// 	log.Fatalf("[constructNodeMap] read node seq err: n(%d) != NBntByteLen(%d\n)", n, NBntByteLen)
		// }
		readnodeNum++
		node.ID = nodeID
		var key [NODEMAP_KEY_LEN]uint64
		copy(key[:], node.Seq)
		//fmt.Printf("[constructNodeMap] key: %v\n\tnode: %v\n", key, node)
		if _, ok := nodeMap[key]; ok == false {
			nodeMap[key] = node
			nodeID++
			//fmt.Printf("[constructNodeMap] node: %v\n", node)
		} else {
			//fmt.Printf("[constructNodeMap] repeat node: %v\n", node)
		}
	}

	fmt.Printf("[constructNodeMap] read node number is : %d\n", readnodeNum)

	return nodeID
}

func AddNodeToNodeMap(node DBGNode, nodeMap map[[NODEMAP_KEY_LEN]uint64]DBGNode, nodeID DBG_MAX_INT) DBG_MAX_INT {
	/*if node.Flag != 1 {
		log.Fatalf("[AddNodeToNodeMap] found node.Flag: %v != 1\n", node.Flag)
	}*/
	var key [NODEMAP_KEY_LEN]uint64
	copy(key[:], node.Seq)
	if _, ok := nodeMap[key]; ok == false {
		node.ID = nodeID
		nodeMap[key] = node
		nodeID++
	} else {
		log.Fatalf("[AddNodeToNodeMap] node: %v has been exist in the nodeMap\n", node)
	}

	return nodeID
}

func AddNewDBGNode(narr []DBGNode, anc <-chan DBGNode, finish <-chan bool) {
loop:
	for {
		select {
		case nd := <-anc:
			narr = append(narr, nd)
		case <-finish:
			for nd := range anc {
				narr = append(narr, nd)
			}
			break loop
		}
	}
}

func CollectAddedDBGNode(anc chan DBGNode, nodeMap map[[NODEMAP_KEY_LEN]uint64]DBGNode, nc chan<- DBGNode, nodeID *DBG_MAX_INT, readNodeMapFinishedC <-chan int) {
	var narr []DBGNode
	var addedNum int
loop:
	for {
		select {
		case nd := <-anc:
			narr = append(narr, nd)
		case <-readNodeMapFinishedC:
			break loop
		}
	}

	narrFlag := make(chan bool, 1)
	j := 0
	for {
		lj := j
		nl := len(narr)
		//fmt.Printf("[CollectAddedDBGNode]before select j:%d, len(narr): %v\n", j, len(narr))
		if nl > 0 && j < nl {
			narrFlag <- true
		}
		select {
		case nd := <-anc:
			narr = append(narr, nd)
		case <-narrFlag:
			{
				var key [NODEMAP_KEY_LEN]uint64
				copy(key[:], narr[j].Seq)
				muRW.RLock()
				_, ok := nodeMap[key]
				muRW.RUnlock()
				if !ok {
					narr[j].ID = *nodeID
					*nodeID++
					muRW.Lock()
					nodeMap[key] = narr[j]
					muRW.Unlock()
					fmt.Fprintf(os.Stderr, "[CollectAddedDBGNode] new added node: %v\n", narr[j])
					nc <- narr[j]
					addedNum++
				}
				j++
			}
		default:
			time.Sleep(time.Microsecond)
		}
		//fmt.Printf("[CollectAddedDBGNode]after select j:%d, len(narr): %v\n", j, len(narr))
		//if lj == j && len(narr) > 0 && j < len(narr) {
		if lj == j && nl > 0 && j < nl {
			<-narrFlag
		}
		if j == len(narr) && len(nc) == 0 && len(anc) == 0 {
			time.Sleep(time.Second)
			if j == len(narr) && len(nc) == 0 && len(anc) == 0 {
				break
			}
		}
	}
	close(anc)
	close(nc)
	time.Sleep(time.Second * 5)
	fmt.Printf("[CollectAddedDBGNode] added node number is: %d, total node number is: %d\n", addedNum, *nodeID)
}

// ChangeNodeMap add new Node to the nodeMap and check node edge has been output
/*func ChangeNodeMap(nodeMap map[[NODEMAP_KEY_LEN]uint64]DBGNode, anc chan<- DBGNode, finishedC <-chan int, nIEC <-chan NodeInfoByEdge, flagNIEC chan<- NodeInfoByEdge, Kmerlen int, nodeID DBG_MAX_INT) (nID DBG_MAX_INT, edgeID DBG_MAX_INT) {
	oldNodeID := nodeID
	edgeID = DBG_MAX_INT(2)
loop:
	for {
		select {
		case <-finishedC:
			break loop
		case nie := <-nIEC:
			var key [NODEMAP_KEY_LEN]uint64
			copy(key[:], nie.n1.Seq)
			v1 := nodeMap[key]
			if nie.c1 == DIN {
				if v1.EdgeIDIncoming[nie.i1] == 1 {
					v1.EdgeIDIncoming[nie.i1] = edgeID
				} else {
					nie.Flag = false
					fmt.Printf("[ChangeNodeMap] add edge: start node ID: %v, end node ID: %v, same as edgeID: %v\n", nie.n1.ID, nie.n2.ID, v1.EdgeIDIncoming[nie.i1])
					flagNIEC <- nie
					continue loop
				}
			} else {
				if v1.EdgeIDOutcoming[nie.i1] == 1 {
					v1.EdgeIDOutcoming[nie.i1] = edgeID
				} else {
					nie.Flag = false
					fmt.Printf("[ChangeNodeMap] add edge: start node ID: %v, end node ID: %v, same as edgeID: %v\n", nie.n1.ID, nie.n2.ID, v1.EdgeIDOutcoming[nie.i1])
					flagNIEC <- nie
					continue loop
				}
			}
			nie.edgeID = edgeID
			nie.n1 = v1
			nie.Flag = true
			copy(key[:], v1.Seq)
			nodeMap[key] = v1
			nd := nie.n2
			if len(nd.Seq) > 0 {
				tn := GetMinDBGNode(nd, Kmerlen)
				copy(key[:], tn.Seq)
				v2, ok := nodeMap[key]
				if !ok { // this is new node
					v2 = tn
				}
				if reflect.DeepEqual(v2.Seq, nd.Seq) {
					if nie.c2 == DIN {
						v2.EdgeIDIncoming[nie.i2] = edgeID
					} else {
						v2.EdgeIDOutcoming[nie.i2] = edgeID
					}
				} else {
					b := bnt.BntRev[nie.i2]
					if nie.c2 == DIN {
						v2.EdgeIDOutcoming[b] = edgeID
					} else {
						v2.EdgeIDIncoming[b] = edgeID
					}
				}
				copy(key[:], v2.Seq)
				if ok {
					nodeMap[key] = v2
					nie.n2 = v2
				} else {
					v2.Flag = 1
					nodeID = AddNodeToNodeMap(v2, nodeMap, nodeID)
					//fmt.Printf("[ChangeNodeMap] v2: %v\nnodeID: %v\n", v2, nodeID-1)
					t := nodeMap[key]
					nie.n2 = t
					anc <- t
				}
			}
			edgeID++
			flagNIEC <- nie
		}
	}

	fmt.Printf("[ChangeNodeMap] added nodes number is : %d\n", nodeID-oldNodeID)
	nID = nodeID
	return
}*/

var muRW sync.RWMutex

// ReadDBGNodeToChan  read DBG nodeMap and simultaneously add new node to the nodeMap
func ReadDBGNodeToChan(nodeArr []DBGNode, nodeMap map[[NODEMAP_KEY_LEN]uint64]DBGNode, nc chan<- DBGNode, readNodeMapFinished chan<- int) {
	for _, value := range nodeArr {
		if len(value.Seq) > 0 && value.Flag == 0 {
			var key [NODEMAP_KEY_LEN]uint64
			copy(key[:], value.Seq)
			muRW.RLock()
			value = nodeMap[key]
			muRW.RUnlock()
			nc <- value
		}
		//value.Flag = uint8(1)
		//nodeMap[string(value.Seq)] = value
		//}
	}

	// notice Function ChangeNodeMap() has finished read nodeMap
	time.Sleep(time.Second)
	readNodeMapFinished <- 1
	close(readNodeMapFinished)
}

func ReverseByteArr(ba []byte) {
	lba := len(ba)
	for i := 0; i < lba/2; i++ {
		ba[i], ba[lba-1-i] = ba[lba-1-i], ba[i]
	}
}

func GetReverseByteArr(ba []byte) (na []byte) {
	na = make([]byte, len(ba))
	for i := 0; i < len(ba); i++ {
		na[len(ba)-1-i] = ba[i]
	}
	return
}

func GetReverseByteArr2(ba, rArr []byte) []byte {
	if cap(rArr) < len(ba) {
		rArr = make([]byte, len(ba))
	} else {
		rArr = rArr[:len(ba)]
	}
	for i := 0; i < len(ba); i++ {
		rArr[len(ba)-1-i] = ba[i]
	}
	return rArr
}

func GetCompByteArr(seq []byte) (cseq []byte) {
	cseq = make([]byte, len(seq))
	for i, b := range seq {
		cseq[i] = bnt.BntRev[b]
	}
	return
}

func ReverseCompByteArr(seq []byte) {
	ls := len(seq)
	for i := 0; i < ls/2; i++ {
		seq[i], seq[ls-1-i] = bnt.BntRev[seq[ls-1-i]], bnt.BntRev[seq[i]]
	}
	if ls%2 != 0 {
		seq[ls/2] = bnt.BntRev[seq[ls/2]]
	}
}

func GetReverseCompByteArr(seq []byte) []byte {
	sl := len(seq)
	rv := make([]byte, sl)
	for i := 0; i < len(rv); i++ {
		rv[i] = bnt.BntRev[seq[sl-1-i]]
	}

	return rv
}

func GetReverseCompByteArr2(seq []byte, rSeq []byte) []byte {
	sl := len(seq)
	var rv []byte
	if cap(rSeq) < sl {
		rv = make([]byte, sl)
	} else {
		rv = rSeq[:sl]
	}
	for i := 0; i < len(seq); i++ {
		rv[i] = bnt.BntRev[seq[sl-1-i]]
	}

	return rv
}

func ReverseUint8Arr(ua []uint8) {
	lua := len(ua)
	for i := 0; i < lua/2; i++ {
		ua[i], ua[lua-1-i] = ua[lua-1-i], ua[i]
	}
}

func GetEdges(cf cuckoofilter.CuckooFilter, kb, rkb constructcf.KmerBnt, count uint8, direction uint8, MIN_KMER_COUNT uint16) (edge DBGEdge, nd DBGNode) {
	var kb2, rb2 constructcf.KmerBnt
	kb2.Len, rb2.Len = kb.Len-1, rkb.Len-1
	kb2.Seq = make([]byte, cf.Kmerlen-1)
	rb2.Seq = make([]byte, cf.Kmerlen-1)
	//tb.Seq = make([]uint64, (kmerlen+bnt.NumBaseInUint64-1)/bnt.NumBaseInUint64)
	//seq := constructcf.ExtendKmerBnt2Byte(kb)
	if direction == FORWARD {
		edge.Utg.Ks = append(edge.Utg.Ks, kb.Seq...)
		edge.Utg.Kq = make([]uint8, cf.Kmerlen-1)
		edge.Utg.Kq = append(edge.Utg.Kq, count)
		bi := kb.Seq[0]
		var nkb, nr constructcf.KmerBnt
		nkb.Len, nr.Len = kb.Len-1, rkb.Len-1
		nkb.Seq = make([]byte, cf.Kmerlen-1)
		nr.Seq = make([]byte, cf.Kmerlen-1)
		copy(nkb.Seq, kb.Seq[1:])
		copy(nr.Seq, rkb.Seq[:rkb.Len-1])
		m := make([]byte, cf.Kmerlen*2)
		for {
			node := ExtendNodeKmer(nkb, nr, cf, MIN_KMER_COUNT, m)
			var leftcount, rightcount int
			var baseBnt byte
			for i := 0; i < bnt.BaseTypeNum; i++ {
				if node.EdgeIDIncoming[i] == 1 {
					leftcount++
				}
				if node.EdgeIDOutcoming[i] == 1 {
					baseBnt = uint8(i)
					rightcount++
				}
			}
			if leftcount == 1 && rightcount == 1 {
				edge.Utg.Ks = append(edge.Utg.Ks, baseBnt)
				edge.Utg.Kq = append(edge.Utg.Kq, uint8(MIN_KMER_COUNT))
				copy(kb2.Seq, nkb.Seq[1:])
				kb2.Seq[kb2.Len-1] = baseBnt
				rb2.Seq[0] = bnt.BntRev[baseBnt]
				copy(rb2.Seq[1:], nr.Seq[:nr.Len-1])
				//nkb = constructcf.GetNextKmer(nkb, uint64(baseBnt), cf.Kmerlen-1)
				//kb2 = constructcf.NoAllocGetNextKmer(nkb, kb2, uint64(baseBnt), cf.Kmerlen-1)
				//rb2 = constructcf.NoAllocGetPreviousKmer(nr, rb2, uint64(bnt.BntRev[baseBnt]), cf.Kmerlen-1)
				nkb, kb2 = kb2, nkb
				nr, rb2 = rb2, nr
				bi = edge.Utg.Ks[len(edge.Utg.Ks)-cf.Kmerlen]
			} else {
				if leftcount > 1 || rightcount > 1 {
					nd = node
					nd.Seq = constructcf.GetReadBntKmer(nkb.Seq, 0, nkb.Len)
					//nd.Flag = leftcount + rightcount
					nd.EdgeIDIncoming[bi] = math.MaxUint32
				}
				break
			}
		}
	} else {
		seq := make([]byte, kb.Len)
		copy(seq, kb.Seq)
		ReverseByteArr(seq)
		edge.Utg.Ks = append(edge.Utg.Ks, seq...)
		edge.Utg.Kq = make([]uint8, cf.Kmerlen-1)
		edge.Utg.Kq = append(edge.Utg.Kq, count)
		bi := seq[0]
		var nkb, nr constructcf.KmerBnt
		nkb.Len, nr.Len = kb.Len-1, rkb.Len-1
		nkb.Seq = make([]byte, nkb.Len)
		copy(nkb.Seq, kb.Seq[:kb.Len-1])
		nr.Seq = make([]byte, nr.Len)
		copy(nr.Seq, rkb.Seq[1:])
		m := make([]byte, cf.Kmerlen*2)
		for {
			node := ExtendNodeKmer(nkb, nr, cf, MIN_KMER_COUNT, m)
			var leftcount, rightcount int
			var baseBnt byte
			for i := 0; i < bnt.BaseTypeNum; i++ {
				if node.EdgeIDIncoming[i] == 1 {
					baseBnt = uint8(i)
					leftcount++
				}
				if node.EdgeIDOutcoming[i] == 1 {
					rightcount++
				}
			}
			if leftcount == 1 && rightcount == 1 {
				edge.Utg.Ks = append(edge.Utg.Ks, baseBnt)
				edge.Utg.Kq = append(edge.Utg.Kq, uint8(MIN_KMER_COUNT))
				kb2.Seq[0] = baseBnt
				copy(kb2.Seq[1:], nkb.Seq[:nkb.Len-1])
				copy(rb2.Seq[:rb2.Len-1], nr.Seq[1:])
				rb2.Seq[rb2.Len-1] = bnt.BntRev[baseBnt]
				//nkb = constructcf.GetPreviousKmer(nkb, uint64(baseBnt), cf.Kmerlen-1)
				//nr = constructcf.GetNextKmer(nr, uint64(bnt.BntRev[baseBnt]), cf.Kmerlen-1)
				//kb2 = constructcf.NoAllocGetPreviousKmer(nkb, kb2, uint64(baseBnt), cf.Kmerlen-1)
				//rb2 = constructcf.NoAllocGetNextKmer(nr, rb2, uint64(bnt.BntRev[baseBnt]), cf.Kmerlen-1)
				nkb, kb2 = kb2, nkb
				nr, rb2 = rb2, nr
				bi = edge.Utg.Ks[len(edge.Utg.Ks)-cf.Kmerlen]
			} else {
				if leftcount > 1 || rightcount > 1 {
					nd = node
					nd.Seq = constructcf.GetReadBntKmer(nkb.Seq, 0, nkb.Len)
					nd.EdgeIDOutcoming[bi] = math.MaxUint32
				}
				// Reverse the seq and quality count
				ReverseByteArr(edge.Utg.Ks)
				ReverseUint8Arr(edge.Utg.Kq)
				break
			}
		}
	}

	return
}

//var mu sync.Mutex
//var mapRWMu sync.RWMutex

type EdgeNode struct {
	Edge         DBGEdge
	NodeS, NodeE DBGNode // NodeS note Start Node, NodeE note End Node
}

func paraGenerateDBGEdges(nc <-chan DBGNode, cf cuckoofilter.CuckooFilter, wc chan<- EdgeNode, MinKmerFreq uint16, MaxNGSReadLen int) {
	for {
		node, _ := <-nc
		if node.ID == 0 {
			var en EdgeNode
			wc <- en
			break
		}
		// read edge seq from cuckoofilter
		var nb constructcf.KmerBnt
		nb.Seq = constructcf.ExtendKmerBnt2Byte(node.Seq, cf.Kmerlen-1)
		nb.Len = cf.Kmerlen - 1
		rnb := constructcf.GetReverseComplet(nb)
		rnb.Len = cf.Kmerlen - 1
		//extRSeq := constructcf.ExtendKmerBnt2Byte(kb)
		// leftcount, rightcount, _, _ := ExtendNodeKmer(extRBnt, cf, MIN_KMER_COUNT, FORWARD)
		//fmt.Printf("[paraGenerateDBGEdges] nb: %v\n\trnb:%v\n", nb, rnb)
		for i := uint(0); i < bnt.BaseTypeNum; i++ {
			bi := byte(i)
			if node.EdgeIDIncoming[i] == 1 {
				var kb, rb constructcf.KmerBnt
				kb.Len, rb.Len = cf.Kmerlen, cf.Kmerlen
				kb.Seq = make([]byte, kb.Len)
				rb.Seq = make([]byte, rb.Len)
				kb.Seq[0] = bi
				copy(kb.Seq[1:], nb.Seq)
				copy(rb.Seq[:rb.Len-1], rnb.Seq)
				rb.Seq[rb.Len-1] = bnt.BntRev[bi]
				min := kb.Seq
				if constructcf.BiggerThan(kb.Seq, rb.Seq) {
					min = rb.Seq
				}
				count := cf.GetCountAllowZero(min)
				//fmt.Printf("[paraGenerateDBGEdges]min: %v, count: %v\n", min, count)
				if count < MinKmerFreq {
					log.Fatalf("[paraGenerateDBGEdges] found count[%v]: < [%v], node: %v", count, MinKmerFreq, node)
				}
				// get edge sequence
				edge, nd := GetEdges(cf, kb, rb, uint8(count), BACKWARD, MinKmerFreq)
				//fmt.Printf("[paraGenerateDBGEdges]Incoming i:%v, edge: %v\n\tnd: %v\n", i, edge, nd)
				//writedEdge := false
				if len(nd.Seq) > 0 || len(edge.Utg.Ks) > MaxNGSReadLen {
					//edge.EndNID = node.ID
					var en EdgeNode
					en.NodeE.Seq = constructcf.GetReadBntKmer(edge.Utg.Ks, len(edge.Utg.Ks)-(cf.Kmerlen-1), cf.Kmerlen-1)
					en.NodeE.EdgeIDIncoming[edge.Utg.Ks[len(edge.Utg.Ks)-cf.Kmerlen]] = math.MaxUint32
					en.Edge = edge
					if len(nd.Seq) > 0 {
						en.NodeS.Seq = constructcf.GetReadBntKmer(edge.Utg.Ks, 0, cf.Kmerlen-1)
						if reflect.DeepEqual(en.NodeS.Seq, nd.Seq) {
							en.NodeS = nd
						} else {
							for j := 0; j < bnt.BaseTypeNum; j++ {
								en.NodeS.EdgeIDIncoming[j] = nd.EdgeIDOutcoming[bnt.BaseTypeNum-1-j]
								en.NodeS.EdgeIDOutcoming[j] = nd.EdgeIDIncoming[bnt.BaseTypeNum-1-j]
							}
						}
						en.NodeS.EdgeIDOutcoming[edge.Utg.Ks[cf.Kmerlen-1]] = math.MaxUint32
						//fmt.Printf("[paraGenerateDBGEdges]nd: %v\n\ten.NodeS: %v\n", nd, en.NodeS)
					}
					wc <- en
				}
			}

			if node.EdgeIDOutcoming[i] == 1 {
				var kb, rb constructcf.KmerBnt
				kb.Len, rb.Len = cf.Kmerlen, cf.Kmerlen
				kb.Seq = make([]byte, kb.Len)
				rb.Seq = make([]byte, rb.Len)
				copy(kb.Seq[:kb.Len-1], nb.Seq)
				kb.Seq[kb.Len-1] = bi
				rb.Seq[0] = bnt.BntRev[bi]
				copy(rb.Seq[1:], rnb.Seq)
				min := kb.Seq
				if BiggerThan(kb.Seq, rb.Seq) {
					min = rb.Seq
				}
				count := cf.GetCountAllowZero(min)
				//fmt.Printf("[paraGenerateDBGEdges]min: %v, count: %v\n", min, count)
				if count < MinKmerFreq {
					log.Fatalf("[paraGenerateDBGEdges] found count[%v]: < [%v], node: %v", count, MinKmerFreq, node)
				}
				edge, nd := GetEdges(cf, kb, rb, uint8(count), FORWARD, MinKmerFreq)
				//fmt.Printf("[paraGenerateDBGEdges]Outcoming i:%v, edge: %v\n\tnd: %v\n", i, edge, nd)
				// writedEdge := false
				if len(nd.Seq) > 0 || len(edge.Utg.Ks) > MaxNGSReadLen {
					//edge.StartNID = node.ID
					var en EdgeNode
					en.Edge = edge
					en.NodeS.Seq = constructcf.GetReadBntKmer(edge.Utg.Ks, 0, cf.Kmerlen-1)
					en.NodeS.EdgeIDOutcoming[edge.Utg.Ks[cf.Kmerlen-1]] = math.MaxUint32
					if len(nd.Seq) > 0 {
						en.NodeE.Seq = constructcf.GetReadBntKmer(edge.Utg.Ks, len(edge.Utg.Ks)-(cf.Kmerlen-1), cf.Kmerlen-1)
						if reflect.DeepEqual(en.NodeE.Seq, nd.Seq) {
							en.NodeE = nd
						} else {
							for j := 0; j < bnt.BaseTypeNum; j++ {
								en.NodeE.EdgeIDIncoming[j] = nd.EdgeIDOutcoming[bnt.BaseTypeNum-1-j]
								en.NodeE.EdgeIDOutcoming[j] = nd.EdgeIDIncoming[bnt.BaseTypeNum-1-j]
							}
						}
						en.NodeE.EdgeIDIncoming[edge.Utg.Ks[len(edge.Utg.Ks)-cf.Kmerlen]] = math.MaxUint32
						//fmt.Printf("[paraGenerateDBGEdges]nd: %v\n\ten.NodeE: %v\n", nd, en.NodeE)
					}
					wc <- en
				}
			}
		}
	}
}

// WritefqRecord write one record of fastq to the file
func WritefqRecord(edgesbuffp io.Writer, ei DBGEdge) {
	fmt.Fprintf(edgesbuffp, "@%d\t%d\t%d\n", ei.ID, ei.StartNID, ei.EndNID)
	// fmt.Fprintf(edgesgzfp, "%s\n", ei.Utg.Ks)
	for i := 0; i < len(ei.Utg.Ks); i++ {
		if ei.Utg.Ks[i] > 3 || ei.Utg.Ks[i] < 0 {
			log.Fatalf("[WriteEdgesToFn] not correct base: %d:%d\n", i, ei.Utg.Ks[i])
		}
		fmt.Fprintf(edgesbuffp, "%c", bnt.BitNtCharUp[ei.Utg.Ks[i]])
	}
	fmt.Fprintf(edgesbuffp, "\n+\n")
	// write quality to the file
	for i := 0; i < len(ei.Utg.Kq); i++ {
		fmt.Fprintf(edgesbuffp, "%c", ei.Utg.Kq[i]+33)
	}
	if len(ei.Utg.Kq) != len(ei.Utg.Ks) {
		log.Fatalf("[WriteEdgesToFn] len(ei.Utg.Kq):%d != len(ei.Utg.Ks):%d\n", len(ei.Utg.Kq), len(ei.Utg.Ks))
	}
	fmt.Fprintf(edgesbuffp, "\n")
	// fmt.Printf("[WriteEdgesToFn] edge: %v\n", ei)
}

// WriteEdgesToFn write edges seq to the file
func WriteEdgesToFn(edgesbrfn string, wc <-chan EdgeNode, numCPU int, nodeMap map[[NODEMAP_KEY_LEN]uint64]DBGNode, anc chan<- DBGNode, kmerlen int) (edgeID DBG_MAX_INT) {
	//oldNodeID := nodeID
	edgeID = DBG_MAX_INT(2)
	edgesNum := 0
	edgesfp, err := os.Create(edgesbrfn)
	if err != nil {
		log.Fatalf("[WriteEdgesToFn] Open file %s failed: %v\n", edgesbrfn, err)
	}
	defer edgesfp.Close()
	brfp := cbrotli.NewWriter(edgesfp, cbrotli.WriterOptions{Quality: 1})
	defer brfp.Close()
	edgesbuffp := bufio.NewWriterSize(brfp, 1<<20)

	finishNum := 0
	for {
		en := <-wc
		ei := en.Edge
		if len(ei.Utg.Ks) == 0 {
			if ei.StartNID != 0 || ei.EndNID != 0 {
				log.Fatalf("[WriteEdgesToFn] err edge: %v\n", ei)
			}
			//fmt.Printf("[WriteEdgesToFn]  edge: %v\n", ei)
			finishNum++
			if finishNum == numCPU {
				break
			}
			continue
		}
		//fmt.Printf("[WriteEdgesToFn] edgeID: %v, en: %v\n", edgeID, en)
		//fmt.Printf("[WriteEdgesToFn] edgeID: %v,len(en.Edge.Utg.Ks): %v,  en.NodeS: %v, en.NodeE: %v\n", edgeID, len(en.Edge.Utg.Ks), en.NodeS, en.NodeE)
		// set edge's node info
		{
			//muRW.Lock()
			var keyS, keyE [NODEMAP_KEY_LEN]uint64
			var vS, vE, tnS, tnE DBGNode
			var okS, okE bool
			if len(en.NodeS.Seq) > 0 {
				tnS = GetMinDBGNode(en.NodeS, kmerlen)
				copy(keyS[:], tnS.Seq)
				muRW.RLock()
				vS, okS = nodeMap[keyS]
				muRW.RUnlock()
				if !okS {
					tnS = ChangeEdgeIDComing(tnS)
					anc <- tnS
				}

				//test code
				/*{
					var rb constructcf.KmerBnt
					rb.Seq = en.NodeS.Seq
					rb.Len = kmerlen - 1
					//revRb := constructcf.GetReadBntKmer(en.Edge.Utg.Ks[:kmerlen-1], 0, kmerlen)
					extNode := constructcf.ExtendKmerBnt2Byte(rb)
					fmt.Printf("[WriteEdgesToFn]base: %v, extNodeS: %v\n", en.Edge.Utg.Ks[kmerlen-1], extNode)
				}*/
			}

			if len(en.NodeE.Seq) > 0 {
				tnE = GetMinDBGNode(en.NodeE, kmerlen)
				copy(keyE[:], tnE.Seq)
				muRW.RLock()
				vE, okE = nodeMap[keyE]
				muRW.RUnlock()
				if !okE {
					tnE = ChangeEdgeIDComing(tnE)
					anc <- tnE
				}

				//test code
				/*{
					var rb constructcf.KmerBnt
					rb.Seq = en.NodeE.Seq
					rb.Len = kmerlen - 1
					extNode := constructcf.ExtendKmerBnt2Byte(rb)
					fmt.Printf("[WriteEdgesToFn]base: %v, extNodeE: %v\n", en.Edge.Utg.Ks[len(en.Edge.Utg.Ks)-kmerlen], extNode)
				}*/
			}
			if okS == false || okE == false {
				//fmt.Printf("[WriteEdgesToFn] okS: %v, okE: %v\n", okS, okE)
				continue
			}

			hasWrite := false
			{ // change nodeMap nodeS EdgeIDComing
				for j := 0; j < bnt.BaseTypeNum; j++ {
					if tnS.EdgeIDOutcoming[j] == math.MaxUint32 {
						if vS.EdgeIDOutcoming[j] == 1 {
							vS.EdgeIDOutcoming[j] = edgeID
						} else if vS.EdgeIDOutcoming[j] > 1 {
							hasWrite = true
						} else {
							log.Fatalf("[WriteEdgesToFn] the edge start node EdgeID not set 1, just == 0, vS: %v\n\ttnS: %v\n", vS, tnS)
						}
						break
					}

					if tnS.EdgeIDIncoming[j] == math.MaxUint32 {
						if vS.EdgeIDIncoming[j] == 1 {
							vS.EdgeIDIncoming[j] = edgeID
						} else if vS.EdgeIDIncoming[j] > 1 {
							hasWrite = true
						} else {
							log.Fatalf("[WriteEdgesToFn] the edge start node EdgeID not set 1, just == 0, vS: %v\n\ttnS: %v\n", vS, tnS)
						}
						break
					}
				}
			}

			if hasWrite {
				continue
			}

			{ // change nodeMap nodeE EdgeIDComing
				//hasWrite := false
				for j := 0; j < bnt.BaseTypeNum; j++ {
					if tnE.EdgeIDIncoming[j] == math.MaxUint32 {
						if vE.EdgeIDIncoming[j] == 1 {
							vE.EdgeIDIncoming[j] = edgeID
						} else if vE.EdgeIDIncoming[j] > 1 {
							hasWrite = true
						} else {
							log.Fatalf("[WriteEdgesToFn] the edge end node EdgeID not set 1, just == 0, vE: %v\n\ttnE: %v\n", vE, tnE)
						}
						break
					}

					if tnE.EdgeIDOutcoming[j] == math.MaxUint32 {
						if vE.EdgeIDOutcoming[j] == 1 {
							vE.EdgeIDOutcoming[j] = edgeID
						} else if vE.EdgeIDOutcoming[j] > 1 {
							hasWrite = true
						} else {
							log.Fatalf("[WriteEdgesToFn] the edge end node EdgeID not set 1, just == 0, vE: %v\n\ttnE: %v\n", vE, tnE)
						}
						break
					}
				}
			}

			if hasWrite {
				continue
			}

			// make sure same edge not write to file
			// if self cycle edge, vS == vE
			//fmt.Printf("[WriteEdgesToFn]changed vS: %v\n\t\tvE: %v\n", vS, vE)
			if vS.ID == vE.ID {
				// need merge EdgeIDcoming
				for j := 0; j < bnt.BaseTypeNum; j++ {
					if vS.EdgeIDIncoming[j] <= 1 && vE.EdgeIDIncoming[j] > 1 {
						vS.EdgeIDIncoming[j] = vE.EdgeIDIncoming[j]
					}
					if vS.EdgeIDOutcoming[j] <= 1 && vE.EdgeIDOutcoming[j] > 1 {
						vS.EdgeIDOutcoming[j] = vE.EdgeIDOutcoming[j]
					}
				}
				muRW.Lock()
				nodeMap[keyS] = vS
				muRW.Unlock()
			} else {
				muRW.Lock()
				nodeMap[keyS] = vS
				nodeMap[keyE] = vE
				muRW.Unlock()
			}
			ei.ID = edgeID
			ei.StartNID, ei.EndNID = vS.ID, vE.ID
			edgeID++
			WritefqRecord(edgesbuffp, ei)
			edgesNum++
			//fmt.Printf("[WriteEdgesToFn] vS:%v\n\tvE:%v\n", nodeMap[keyS], nodeMap[keyE])
			//fmt.Printf("[WriteEdgesToFn] the writed edgeID: %v, ei.StartNID: %v, ei.EndNID: %v\n\tei.Utg.Ks: %v\n", edgeID-1, ei.StartNID, ei.EndNID, ei.Utg.Ks)
		}
	}
	err = edgesbuffp.Flush()
	if err != nil {
		log.Fatalf("[WriteEdgesToFn] write edges to file : %v, err : %v\n", edgesbrfn, err)
	}
	err = brfp.Flush()
	if err != nil {
		log.Fatalf("[WriteEdgesToFn] write edges to file : %v, err : %v\n", edgesbrfn, err)
	}

	fmt.Printf("[WriteEdgesToFn] the writed file edges number is %d\n", edgesNum)
	//fmt.Printf("[WriteEdgesToFn] added nodes number is : %d\n", nodeID-oldNodeID)
	//nID = nodeID
	return
}

/*func ProcessAddedNode(cf cuckoofilter.CuckooFilter, nodeMap map[string]DBGNode, newNodeBntArr []constructcf.ReadBnt, wc chan DBGEdge, nodeID DBG_MAX_INT) (addedNodesNum, addedEdgesNum int) {

	InitialLen := len(newNodeBntArr)
	processAdded := 0
	totalProcessNUm := 0
	// for _, item := range newNodeBntArr {
	for i := 0; i < len(newNodeBntArr); i++ {
		totalProcessNUm++
		var node DBGNode
		node.ID = nodeID
		node.Seq = newNodeBntArr[i].Seq
		if _, ok := nodeMap[string(node.Seq)]; ok == false {
			nodeMap[string(node.Seq)] = node
			nodeID++
			addedNodesNum++
			// check if need added edges
			var rb constructcf.ReadBnt
			rb.Seq = node.Seq
			rb.Length = cf.Kmerlen - 1
			extRBnt := constructcf.ExtendReadBnt2Byte(rb)
			for i := 0; i < bnt.BaseTypeNum; i++ {
				bi := byte(i)
				if node.EdgeIDIncoming[i] == 0 {
					var nBnt constructcf.ReadBnt
					nBnt.Seq = append(nBnt.Seq, bi)
					nBnt.Seq = append(nBnt.Seq, extRBnt.Seq...)
					nBnt.Length = len(nBnt.Seq)
					ks := constructcf.GetReadBntKmer(nBnt, 0, cf.Kmerlen)
					rs := constructcf.ReverseComplet(ks)
					if ks.BiggerThan(rs) {
						ks, rs = rs, ks
					}
					count := cf.GetCountAllowZero(ks.Seq)
					if count >= MIN_KMER_COUNT {
						// get edge sequence
						edge, isNode := GetEdges(cf, nBnt, uint8(count), BACKWARD, MIN_KMER_COUNT)
						writedEdge := false
						if isNode == true {
							var tBnt constructcf.ReadBnt
							tBnt.Seq = edge.Utg.Ks[:cf.Kmerlen-1]
							tks := constructcf.GetReadBntKmer(tBnt, 0, cf.Kmerlen-1)
							trs := constructcf.ReverseComplet(tks)
							sks := tks
							if sks.BiggerThan((trs)) {
								sks = trs
							}
							if v, ok := nodeMap[string(sks.Seq)]; ok {
								c := edge.Utg.Ks[cf.Kmerlen-1]
								if reflect.DeepEqual(sks, tks) {
									// b := bnt.Base2Bnt[c]
									v.EdgeIDOutcoming[c] = 1
								} else {
									// b := bnt.Base2Bnt[bnt.BitNtRev[c]]
									v.EdgeIDIncoming[bnt.BntRev[c]] = 1
								}
								edge.StartNID = v.ID
								writedEdge = true
								nodeMap[string(sks.Seq)] = v
							} else { // is a new node, add to the newNodeBntArr
								newNodeBntArr = append(newNodeBntArr, sks)
								processAdded++
							}
						} else { // is a tip
							if len(edge.Utg.Ks) > 2*cf.Kmerlen {
								writedEdge = true
							}

						}
						if writedEdge == true {

							node.EdgeIDIncoming[i] = 1
							edge.EndNID = node.ID
							wc <- edge
							// fmt.Printf("[ProcessAddedNode] edge: %v\n", edge)
							addedEdgesNum++
						}
					}
				}

				if node.EdgeIDOutcoming[i] == 0 {
					var nBnt constructcf.ReadBnt
					nBnt.Seq = append(nBnt.Seq, extRBnt.Seq...)
					nBnt.Seq = append(nBnt.Seq, bi)
					nBnt.Length = len(nBnt.Seq)
					ks := constructcf.GetReadBntKmer(nBnt, 0, cf.Kmerlen)
					rs := constructcf.ReverseComplet(ks)
					if ks.BiggerThan(rs) {
						ks, rs = rs, ks
					}
					count := cf.GetCountAllowZero(ks.Seq)
					if count >= MIN_KMER_COUNT {
						edge, isNode := GetEdges(cf, nBnt, uint8(count), FORWARD, MIN_KMER_COUNT)
						writedEdge := false
						if isNode == true {
							var tBnt constructcf.ReadBnt
							tBnt.Seq = edge.Utg.Ks[len(edge.Utg.Ks)-cf.Kmerlen+1:]
							tks := constructcf.GetReadBntKmer(tBnt, 0, cf.Kmerlen-1)
							trs := constructcf.ReverseComplet(tks)
							sks := tks
							if sks.BiggerThan(trs) {
								sks = trs
							}
							if v, ok := nodeMap[string(sks.Seq)]; ok {
								c := edge.Utg.Ks[len(edge.Utg.Ks)-cf.Kmerlen]
								if reflect.DeepEqual(sks, tks) {
									// b := bnt.Base2Bnt[c]
									v.EdgeIDIncoming[c] = 1
								} else {
									// b := bnt.Base2Bnt[bnt.BitNtRev[c]]
									v.EdgeIDOutcoming[bnt.BntRev[c]] = 1
								}
								nodeMap[string(sks.Seq)] = v
								edge.EndNID = v.ID
								writedEdge = true
							} else {
								newNodeBntArr = append(newNodeBntArr, sks)
								processAdded++
							}
						} else { // is a tip
							if len(edge.Utg.Ks) > 2*cf.Kmerlen {
								writedEdge = true
							}
						}
						if writedEdge == true {
							node.EdgeIDOutcoming[i] = 1
							edge.StartNID = node.ID
							wc <- edge
							// fmt.Printf("[ProcessAddedNode] edge: %v\n", edge)
						}
					}
				}
			}
		}
	}

	// add a nil edge to wc, tell have not any more edge need to write
	var edge DBGEdge
	wc <- edge

	fmt.Printf("[ProcessAddedNode] initial newNodeBntArr len: %d, added node number: %d, at last newNodeBntArr len:%d, totalProcessNUm: %d\n", InitialLen, processAdded, len(newNodeBntArr), totalProcessNUm)

	return
} */

/*func cleanEdgeIDInNodeMap(nodeMap map[string]DBGNode) {
	for k, v := range nodeMap {
		for i := 0; i < bnt.BaseTypeNum; i++ {
			v.EdgeIDIncoming[i] = 0
			v.EdgeIDOutcoming[i] = 0
		}
		nodeMap[k] = v
	}
}*/

func GenerateDBGEdges(nodeMap map[[NODEMAP_KEY_LEN]uint64]DBGNode, cf cuckoofilter.CuckooFilter, edgesbrfn string, numCPU int, nodeID DBG_MAX_INT, MinKmerFreq uint16, MaxNGSReadLen int) (newNodeID DBG_MAX_INT, edgeID DBG_MAX_INT) {
	bufsize := 50000
	nc := make(chan DBGNode)
	wc := make(chan EdgeNode, bufsize)
	defer close(wc)
	readNodeMapFinishedC := make(chan int)
	nodeArr := make([]DBGNode, nodeID)
	idx := 0
	for _, value := range nodeMap {
		if len(value.Seq) > 0 && value.Flag == 0 {
			nodeArr[idx] = value
			idx++
		}
	}
	nodeArr = nodeArr[:idx]
	// Read DBGNode to the nc
	go ReadDBGNodeToChan(nodeArr, nodeMap, nc, readNodeMapFinishedC)
	// parallel construct edges from cuckoofilter
	for i := 0; i < numCPU; i++ {
		go paraGenerateDBGEdges(nc, cf, wc, MinKmerFreq, MaxNGSReadLen)
	}
	// collect added node and pass DBGNode to the nc
	anc := make(chan DBGNode, bufsize)
	//totalNodeNum := make(chan DBG_MAX_INT，1)
	go CollectAddedDBGNode(anc, nodeMap, nc, &nodeID, readNodeMapFinishedC)
	// write edges Seq to the file
	edgeID = WriteEdgesToFn(edgesbrfn, wc, numCPU, nodeMap, anc, cf.Kmerlen)
	newNodeID = nodeID
	// Change nodeMap monitor function
	//newNodeID, edgeID = ChangeNodeMap(nodeMap, anc, finishedC, nIEC, flagNIEC, cf.Kmerlen, nodeID)

	/* // cache the new added node Bnt info
	var newNodeBntArr []constructcf.ReadBnt
	finishedCount := 0
	for {
		nb := <-newNodeChan
		if nb.Length == 0 {
			finishedCount++
			if finishedCount == numCPU {
				break
			} else {
				continue
			}
		}

		newNodeBntArr = append(newNodeBntArr, nb)
	}

	// process added nodes' edges
	addedNodesNum, addedEdgesNum := ProcessAddedNode(cf, nodeMap, newNodeBntArr, wc, nodeID)
	newNodeID = nodeID + DBG_MAX_INT(addedNodesNum) */

	// clean set edgeID in the DBGNode
	//cleanEdgeIDInNodeMap(nodeMap)
	// fmt.Printf("[GenerateDBGEdges] added nodes number is : %d, added edges number is : %d\n", addedNodesNum, addedEdgesNum)
	return
}

func DBGStatWriter(DBGStatfn string, newNodeID, edgeID DBG_MAX_INT) {
	DBGStatfp, err := os.Create(DBGStatfn)
	if err != nil {
		log.Fatalf("[DBGStatWriter] file %s create error: %v\n", DBGStatfn, err)
	}
	defer DBGStatfp.Close()
	fmt.Fprintf(DBGStatfp, "nodes size:\t%v\n", newNodeID)
	fmt.Fprintf(DBGStatfp, "edges size:\t%v\n", edgeID)
}

func DBGStatReader(DBGStatfn string) (nodesSize, edgesSize DBG_MAX_INT) {
	DBGStatfp, err := os.Open(DBGStatfn)
	if err != nil {
		log.Fatalf("[DBGStatReader] file %s Open error: %v\n", DBGStatfn, err)
	}
	defer DBGStatfp.Close()
	if _, err = fmt.Fscanf(DBGStatfp, "nodes size:\t%v\n", &nodesSize); err != nil {
		log.Fatalf("[DBGStatReader] file: %v, nodes size parse error: %v\n", DBGStatfn, err)
	}
	if _, err = fmt.Fscanf(DBGStatfp, "edges size:\t%v\n", &edgesSize); err != nil {
		log.Fatalf("[DBGStatReader] file: %v, edges size parse error: %v\n", DBGStatfn, err)
	}

	return
}

func DBGInfoWriter(DBGInfofn string, edgesArrSize, nodesArrSize int) {
	DBGInfofp, err := os.Create(DBGInfofn)
	if err != nil {
		log.Fatalf("[DBGInfoWriter] file %s create error: %v\n", DBGInfofn, err)
	}
	defer DBGInfofp.Close()
	fmt.Fprintf(DBGInfofp, "edgesArr size:\t%v\n", edgesArrSize)
	fmt.Fprintf(DBGInfofp, "nodesArr size:\t%v\n", nodesArrSize)
}

func DBGInfoReader(DBGInfofn string) (edgesArrSize, nodesArrSize int) {
	DBGInfofp, err := os.Open(DBGInfofn)
	if err != nil {
		log.Fatalf("[DBGInfoReader] file %s Open error: %v\n", DBGInfofn, err)
	}
	defer DBGInfofp.Close()
	_, err = fmt.Fscanf(DBGInfofp, "edgesArr size:\t%v\n", &edgesArrSize)
	if err != nil {
		log.Fatalf("[edgesStatWriter] edgesArr size parse error: %v\n", err)
	}
	_, err = fmt.Fscanf(DBGInfofp, "nodesArr size:\t%v\n", &nodesArrSize)
	if err != nil {
		log.Fatalf("[edgesStatWriter] nodesArr size parse error: %v\n", err)
	}
	return edgesArrSize, nodesArrSize
}

func EdgesStatWriter(edgesStatfn string, edgesSize int) {
	edgesStatfp, err := os.Create(edgesStatfn)
	if err != nil {
		log.Fatalf("[EdgesStatWriter] file %s create error: %v\n", edgesStatfn, err)
	}
	defer edgesStatfp.Close()
	fmt.Fprintf(edgesStatfp, "edges size:\t%v\n", edgesSize)
}

func EdgesStatReader(edgesStatfn string) (edgesSize int) {
	edgesStatfp, err := os.Open(edgesStatfn)
	if err != nil {
		log.Fatalf("[NodesStatReader] file %s Open error: %v\n", edgesStatfn, err)
	}
	defer edgesStatfp.Close()
	_, err = fmt.Fscanf(edgesStatfp, "edges size:\t%v\n", &edgesSize)
	if err != nil {
		log.Fatalf("[edgesStatReader] edges size parse error: %v\n", err)
	}
	return
}

func NodeMapMmapWriter(nodeMap map[[NODEMAP_KEY_LEN]uint64]DBGNode, nodesfn string) {
	nodesfp, err := os.Create(nodesfn)
	if err != nil {
		log.Fatalf("[NodeMapMmapWriter] file %s create error, err: %v\n", nodesfn, err)
	}
	defer nodesfp.Close()
	enc := gob.NewEncoder(nodesfp)
	err = enc.Encode(nodeMap)
	if err != nil {
		log.Fatalf("[NodeMapMmapWriter] encode err: %v\n", err)
	}
}

func NodeMapMmapReader(nodesfn string) (nodeMap map[[NODEMAP_KEY_LEN]uint64]DBGNode) {
	nodesfp, err := os.Open(nodesfn)
	if err != nil {
		log.Fatalf("[NodeMapMmapReader] open file %s failed, err:%v\n", nodesfn, err)
	}
	defer nodesfp.Close()
	dec := gob.NewDecoder(nodesfp)
	err = dec.Decode(&nodeMap)
	if err != nil {
		log.Fatalf("[NodeMapMmapReader] decode failed, err: %v\n", err)
	}

	return
}

func NodesArrWriter(nodesArr []DBGNode, nodesfnbr string, fc chan<- int) error {
	nodesfp, err := os.Create(nodesfnbr)
	if err != nil {
		log.Fatalf("[NodesArrWriter] file %s create error, err: %v\n", nodesfnbr, err)
	}
	defer nodesfp.Close()
	cbrofp := cbrotli.NewWriter(nodesfp, cbrotli.WriterOptions{Quality: 1})
	defer cbrofp.Close()
	buffp := bufio.NewWriterSize(cbrofp, 1<<20) // 1<<25 == 2**25
	//var arr DBGNodeArr
	for i := 2; i < len(nodesArr); i++ {
		nd := nodesArr[i]
		if nd.GetDeleteFlag() > 0 || nd.ID < 2 {
			continue
		}
		id := uint32(nodesArr[i].ID)
		binary.Write(buffp, binary.LittleEndian, id)
		var IDcoming1, IDcoming2 [bnt.BaseTypeNum]uint32
		for j := 0; j < bnt.BaseTypeNum; j++ {
			IDcoming1[j] = uint32(nodesArr[i].EdgeIDIncoming[j])
			IDcoming2[j] = uint32(nodesArr[i].EdgeIDOutcoming[j])
		}
		//binary.Write(buffp, binary.LittleEndian, nd.ID)
		binary.Write(buffp, binary.LittleEndian, IDcoming1)
		binary.Write(buffp, binary.LittleEndian, IDcoming2)
		binary.Write(buffp, binary.LittleEndian, nd.Seq)
		//fmt.Printf("[NodesArrWriter]nID:%d Seq:%v\n", i, nd.Seq)
		//binary.Write(buffp, binary.LittleEndian, nd.ID)

		/*var ndi DBGNodeInfo
		ndi = DBGNodeInfo{
			ID:              uint32(nd.ID),
			EdgeIDIncoming:  []uint32{uint32(nd.EdgeIDIncoming[0]), uint32(nd.EdgeIDIncoming[1]), uint32(nd.EdgeIDIncoming[2]), uint32(nd.EdgeIDIncoming[3])},
			EdgeIDOutcoming: []uint32{uint32(nd.EdgeIDOutcoming[0]), uint32(nd.EdgeIDOutcoming[1]), uint32(nd.EdgeIDOutcoming[2]), uint32(nd.EdgeIDOutcoming[3])},
			Seq:             nd.Seq,
		}
		arr.Arr = append(arr.Arr, &ndi)*/
		/*ndi.ID =  nd.ID
		ndi.EdgeIDIncoming = nd.EdgeIDIncoming
		ndi.EdgeIDOutcoming = nd.EdgeIDOutcoming
		ndi.Seq = nd.Seq*/

	}
	/*data, err := proto.Marshal(&arr)
	if err != nil {
		log.Fatalf("[NodesArrWriter] Mashal data error: %v\n", err)
	}

	buffp.Write(data)*/
	/*js, err1 := json.Marshal(nodesArr)
	if err1 != nil {
		log.Fatalf("[NodesArrWriter] nodesArr failed to json.Marshal file: %s, err1: %v\n", nodesfnbr, err1)
	}
	err = binary.Write(buffp, binary.LittleEndian, js)
	if err != nil {
		log.Fatalf("[NodesArrWriter] nodesArr failed to Write file: %s, err: %v\n", nodesfnbr, err)
	}

	var err1, err2, err3, err4, err5, err6 error
	for i := 2; i < len(nodesArr); i++ {
		id := uint32(nodesArr[i].ID)
		err1 = binary.Write(buffp, binary.LittleEndian, id)
		var IDcoming1, IDcoming2 [bnt.BaseTypeNum]uint32
		for j := 0; j < bnt.BaseTypeNum; j++ {
			IDcoming1[j] = uint32(nodesArr[i].EdgeIDIncoming[j])
			IDcoming2[j] = uint32(nodesArr[i].EdgeIDOutcoming[j])
		}
		err2 = binary.Write(buffp, binary.LittleEndian, IDcoming1)
		err3 = binary.Write(buffp, binary.LittleEndian, IDcoming2)
		err4 = binary.Write(buffp, binary.LittleEndian, nodesArr[i].SubGID)
		err5 = binary.Write(buffp, binary.LittleEndian, nodesArr[i].Flag)
		err6 = binary.Write(buffp, binary.LittleEndian, nodesArr[i].Seq)
		if err1 != nil || err2 != nil || err3 != nil || err4 != nil || err5 != nil || err6 != nil {
			log.Fatalf("[NodesArrWriter] nodesArr failed to Write file: %s, err1: %v, err2: %v,err3: %v,err4:%v,err5:%v, err6: %v\n", nodesfnbr, err1, err2, err3, err4, err5, err6)
		}
	}*/

	/*err = binary.Write(buffp, binary.LittleEndian, nodesArr)
	if err != nil {
		log.Fatalf("[NodesArrWriter] nodesArr failed to Write file: %s, err: %v\n", nodesfnbr, err)
	}*/

	/*enc := gob.NewEncoder(nodesfp)
	err = enc.Encode(nodesArr)
	if err != nil {
		log.Fatalf("[NodeMapMmapWriter] encode err: %v\n", err)
	}*/
	if err1 := buffp.Flush(); err1 != nil {
		log.Fatalf("[NodesArrWriter] failed to flush file: %s, err: %v\n", nodesfnbr, err1)
	}

	if err2 := cbrofp.Flush(); err2 != nil {
		log.Fatalf("[NodesArrWriter] failed to flush file: %s, err: %v\n", nodesfnbr, err2)
	}

	fc <- 1
	close(fc)
	return err
}

func NodesArrReader(nodesfnbr string, nodesArr []DBGNode, kmerlen int, fc chan<- int) {
	nodesfp, err := os.Open(nodesfnbr)
	if err != nil {
		log.Fatalf("[NodesArrReader] open file %s failed, err:%v\n", nodesfnbr, err)
	}
	defer nodesfp.Close()
	brfp := cbrotli.NewReader(nodesfp)
	defer brfp.Close()
	buffp := bufio.NewReaderSize(brfp, 1<<20)

	NSeqLen := (kmerlen - 1 + (bnt.NumBaseInUint64 - 1)) / bnt.NumBaseInUint64
	count := 0
	for {
		var nd DBGNode
		var id uint32
		err1 := binary.Read(buffp, binary.LittleEndian, &id)
		nd.ID = DBG_MAX_INT(id)
		var IDcoming1, IDcoming2 [bnt.BaseTypeNum]uint32
		//err1 := binary.Read(buffp, binary.LittleEndian, &nd.ID)
		err2 := binary.Read(buffp, binary.LittleEndian, &IDcoming1)
		err3 := binary.Read(buffp, binary.LittleEndian, &IDcoming2)
		for j := 0; j < bnt.BaseTypeNum; j++ {
			nd.EdgeIDIncoming[j] = DBG_MAX_INT(IDcoming1[j])
			nd.EdgeIDOutcoming[j] = DBG_MAX_INT(IDcoming2[j])
		}
		seq := make([]uint64, NSeqLen)
		err4 := binary.Read(buffp, binary.LittleEndian, seq)
		if err1 == io.EOF {
			break
		} else if err1 != nil || err2 != nil || err3 != nil || err4 != nil {
			log.Fatalf("[NodesArrReader]nd: %v\n\t nodesArr failed to read file: %s, err1: %v, err2: %v,err3: %v,err4:%v\n", nd, nodesfnbr, err1, err2, err3, err4)
		}
		count++
		nd.Seq = seq
		//fmt.Printf("[NodesArrReader] nd: %v\n", nd)
		nodesArr[nd.ID] = nd
	}
	fmt.Printf("[NodesArrReader] recover nodes number : %v from file: %v\n", count, nodesfnbr)
	fc <- 1

	//buffp := bufio.NewReaderSize(brfp, 1<<21)
	/*buf, err := ioutil.ReadAll(brfp)
	if err != nil {
		log.Fatalf("[NodesArrReader] read file %s failed, err:%v\n", nodesfnbr, err)
	}
	var dbgNodeArr DBGNodeArr
	err = proto.Unmarshal(buf, &dbgNodeArr)
	if err != nil {
		log.Fatalf("[NodesArrReader] proto.Unmarshal() err:%v\n", err)
	}
	for _, ndi := range dbgNodeArr.Arr {
		var nd DBGNode
		id := ndi.ID
		if id < 2 || int(id) >= len(nodesArr) {
			log.Fatalf("[NodesArrReader] node ID: %v must between: [2, %v)\n", id, len(nodesArr))
		}
		nd.ID = DBG_MAX_INT(id)
		nd.EdgeIDIncoming[0], nd.EdgeIDIncoming[1], nd.EdgeIDIncoming[2], nd.EdgeIDIncoming[3] = DBG_MAX_INT(ndi.EdgeIDIncoming[0]), DBG_MAX_INT(ndi.EdgeIDIncoming[1]), DBG_MAX_INT(ndi.EdgeIDIncoming[2]), DBG_MAX_INT(ndi.EdgeIDIncoming[3])
		nd.EdgeIDOutcoming[0], nd.EdgeIDOutcoming[1], nd.EdgeIDOutcoming[2], nd.EdgeIDOutcoming[3] = DBG_MAX_INT(ndi.EdgeIDOutcoming[0]), DBG_MAX_INT(ndi.EdgeIDOutcoming[1]), DBG_MAX_INT(ndi.EdgeIDOutcoming[2]), DBG_MAX_INT(ndi.EdgeIDOutcoming[3])
		nd.Seq = ndi.Seq
		nodesArr[nd.ID] = nd
		//fmt.Printf("[NodesArrReader] nd: %v\n", nd)
	}*/

	/*json.Unmarshal(data, v)
	err4 = binary.Read(buffp, binary.LittleEndian)
	NSeqLen := (Kmerlen - 1 + (bnt.NumBaseInUint64 - 1)) / bnt.NumBaseInUint64
	count := 0
	var err1, err2, err3, err4, err5, err6 error
	for i := 2; i < len(nodesArr); i++ {
		var nd DBGNode
		nd.Seq = make([]uint64, NSeqLen)
		var id uint32
		err1 = binary.Read(buffp, binary.LittleEndian, id)
		nd.ID = DBG_MAX_INT(id)
		var IDcoming1, IDcoming2 [bnt.BaseTypeNum]uint32
		err2 = binary.Read(buffp, binary.LittleEndian, &IDcoming1)
		err3 = binary.Read(buffp, binary.LittleEndian, &IDcoming2)
		for j := 0; j < bnt.BaseTypeNum; j++ {
			nd.EdgeIDIncoming[j] = DBG_MAX_INT(IDcoming1[j])
			nd.EdgeIDOutcoming[j] = DBG_MAX_INT(IDcoming2[j])
		}
		err4 = binary.Read(buffp, binary.LittleEndian, nd.SubGID)
		err5 = binary.Read(buffp, binary.LittleEndian, nd.Flag)
		err6 = binary.Read(buffp, binary.LittleEndian, &nd.Seq)
		if err1 == io.EOF {
			break
		} else if err1 != nil || err2 != nil || err3 != nil || err4 != nil || err5 != nil || err6 != nil {
			log.Fatalf("[NodesArrReader]nd: %v\n\t nodesArr failed to read file: %s, err1: %v, err2: %v,err3: %v,err4:%v,err5:%v, err6: %v\n", nd, nodesfnbr, err1, err2, err3, err4, err5, err6)
		}
		fmt.Printf("[NodesArrReader] nd: %v\n", nd)
		nodesArr[nd.ID] = nd
		count++
	}
	fmt.Printf("[NodesArrReader] recover nodes number : %v from file: %v\n", count, nodesfnbr)
	dec := gob.NewDecoder(nodesfp)
	err = dec.Decode(&nodesArr)
	if err != nil {
		log.Fatalf("[NodeMapMmapReader] decode failed, err: %v\n", err)
	}*/
}

func NodeMap2NodeArr(nodeMap map[[NODEMAP_KEY_LEN]uint64]DBGNode, nodesArr []DBGNode) {
	naLen := DBG_MAX_INT(len(nodesArr))
	for _, v := range nodeMap {
		if v.ID >= naLen {
			log.Fatalf("[NodeMap2NodeArr] v.ID: %v >= nodesArr len: %v\n", v.ID, naLen)
		}
		if v.ID <= 1 || v.ID == math.MaxUint32 {
			continue
		}
		nodesArr[v.ID] = v
	}
}
func writeComplexNodesToFile(complexNodesbrFn string, wc chan DBGNode, numCPU int) (complexNodeNum int) {
	ckfp, err := os.Create(complexNodesbrFn)
	if err != nil {
		log.Fatal(err)
	}
	defer ckfp.Close()
	brfp := cbrotli.NewWriter(ckfp, cbrotli.WriterOptions{Quality: 1})
	defer brfp.Close()
	buffp := bufio.NewWriterSize(brfp, 1<<20)
	// if err != nil {
	// 	log.Fatal(err)
	// }
	//endFlagCount := 0

	// write complex node to the file
	// complexNodeNum := 0
	cpuCount := 0
	for {
		nd := <-wc
		if len(nd.Seq) == 0 {
			cpuCount++
			if cpuCount == numCPU {
				break
			} else {
				continue
			}
		}

		//fmt.Printf("[writeComplexNodesToFile]node: %v\n", nd)
		/*var ni DBGNodeInfo
		ni.Seq = nd.Seq
		ni.EdgeIDIncoming = make([]uint32, bnt.BaseTypeNum)
		ni.EdgeIDOutcoming = make([]uint32, bnt.BaseTypeNum)
		for j := 0; j < bnt.BaseTypeNum; j++ {
			ni.EdgeIDIncoming[j] = uint32(nd.EdgeIDIncoming[j]) + 1
			ni.EdgeIDOutcoming[j] = uint32(nd.EdgeIDOutcoming[j]) + 1
		}
		//var dnArr []DBGNodeInfo
		//dnArr = append(dnArr, ni)
		ob, err := proto.Marshal(&ni)
		if err != nil {
			log.Fatal("[writeComplexNodesToFile] failed to marshal: ", err)
		}
		p := &DBGNodeInfo{}
		err = proto.Unmarshal(ob, p)

		fmt.Printf("[writeComplexNodesToFile]*p: %v\n", *p)
		buffp.Write(ob)*/
		//buffp.Write(nd.EdgeIDIncoming)
		//buffp.Write(nd.EdgeIDOutcoming)
		if err := binary.Write(buffp, binary.LittleEndian, nd.Seq); err != nil {
			log.Fatalf("[writeComplexNodesToFile] write node seq to file err: %v\n", err)
		}
		if err := binary.Write(buffp, binary.LittleEndian, nd.EdgeIDIncoming); err != nil {
			log.Fatalf("[writeComplexNodesToFile] write node seq to file err: %v\n", err)
		}
		if err := binary.Write(buffp, binary.LittleEndian, nd.EdgeIDOutcoming); err != nil {
			log.Fatalf("[writeComplexNodesToFile] write node seq to file err: %v\n", err)
		}
		// *** test code ***
		/* extRB := constructcf.ExtendReadBnt2Byte(rb)
		fmt.Fprintf(os.Stderr, ">%v\n%v\n", complexNodeNum+1, Transform2Letters(extRB.Seq))
		*/
		// *** test code ***
		// ckgzfp.Write(rb.Seq)
		// ckgzfp.Write([]byte("\n"))
		complexNodeNum++
	}

	if err := buffp.Flush(); err != nil {
		log.Fatalf("[writeComplexNodesToFile] write to file: %s err: %v\n", complexNodesbrFn, err)
	}
	if err := brfp.Flush(); err != nil {
		log.Fatalf("[writeComplexNodesToFile] write to file: %s err: %v\n", complexNodesbrFn, err)
	}

	return
}
func CDBG(c cli.Command) {
	var cf cuckoofilter.CuckooFilter
	cf.NumItems = 0
	CDBGFunc(c, cf)
}

func CDBGFunc(c cli.Command, cf cuckoofilter.CuckooFilter) {
	fmt.Println(c.Flags(), c.Parent().Flags())

	// get set arguments
	// t0 := time.Now()
	numCPU, err := strconv.Atoi(c.Parent().Flag("t").String())
	if err != nil {
		log.Fatalf("[CDBG] the argument 't' set error, err: %v\n", err)
	}

	klen, err := strconv.Atoi(c.Parent().Flag("K").String())
	if err != nil {
		log.Fatalf("[CDBG] the argument 'K' set error, err: %v\n", err)
	}
	if klen >= NODEMAP_KEY_LEN*32 {
		log.Fatalf("[CDBG] the argument 'K' must small than [NODEMAP_KEY_LEN * 32]: %v\n", NODEMAP_KEY_LEN*32)
	}
	MinKmerFreq, err := strconv.Atoi(c.Flag("MinKmerFreq").String())
	if err != nil {
		log.Fatalf("[CDBG] the argument 'K' set error, err: %v\n", err)
	}
	MaxNGSReadLen, err := strconv.Atoi(c.Flag("MaxNGSReadLen").String())
	if err != nil {
		log.Fatalf("[CDBG] the argument 'K' set error, err: %v\n", err)
	}
	runtime.GOMAXPROCS(numCPU)
	prefix := c.Parent().Flag("p").String()
	// create cpu profile
	/*profileFn := prefix + ".CDBG.prof"
	cpuprofilefp, err := os.Create(profileFn)
	if err != nil {
		log.Fatalf("[CDBG] open cpuprofile file: %v failed\n", profileFn)
	}
	pprof.StartCPUProfile(cpuprofilefp)
	defer pprof.StopCPUProfile()*/
	// cfinfofn := prefix + ".cfInfo"
	// cf, err := cuckoofilter.RecoverCuckooFilterInfo(cfinfofn)
	// if err != nil {
	// 	log.Fatalf("[CDGB] cuckoofilter recover err: %v\n", err)
	// }

	// find complex Nodes
	t0 := time.Now()
	if cf.NumItems == 0 {
		cfInfofn := prefix + ".cf.Info"
		cf, err = cuckoofilter.RecoverCuckooFilterInfo(cfInfofn)
		if err != nil {
			log.Fatalf("[CDBG] Read CuckooFilter info file: %v err: %v\n", cfInfofn, err)
		}
		cf.Hash = make([]cuckoofilter.Bucket, cf.NumItems)
		fmt.Printf("[CDBG]cf.NumItems: %v, cf.Kmerlen: %v, len(cf.Hash): %v\n", cf.NumItems, cf.Kmerlen, len(cf.Hash))
		cffn := prefix + ".cf.Hash.br"
		err = cf.HashReader(cffn)
		if err != nil {
			log.Fatalf("[CDBG] Read CuckooFilter Hash file: %v err: %v\n", cffn, err)
		}
		fmt.Printf("[CDBG] Load CuckooFilter Struct used : %v\n", time.Now().Sub(t0))
	}

	t0 = time.Now()
	cf.GetStat()
	// fmt.Printf("[CDBG] cf.Hash[0]: %v\n", cf.Hash[0])
	//Kmerlen = cf.Kmerlen
	bufsize := 50000
	cs := make(chan constructcf.KmerBntBucket, numCPU*2)
	wc := make(chan DBGNode, bufsize)
	// read uniq kmers form file
	go ParaReadUniqKmer(prefix, cs, cf.Kmerlen)

	// identify complex Nodes
	for i := 0; i < numCPU; i++ {
		go paraLookupComplexNode(cs, wc, cf, uint16(MinKmerFreq))
	}

	// write complex Nodes to the file
	complexKmerbrfn := prefix + ".complexNode.br"
	complexNodeNum := writeComplexNodesToFile(complexKmerbrfn, wc, numCPU)
	fmt.Printf("[CDBG] found complex Node num is : %d\n", complexNodeNum)
	fmt.Printf("[CDBG] search for complex DBG Node used : %v\n", time.Now().Sub(t0))
	t0 = time.Now()

	// construct Node map
	NBntUint64Len := (cf.Kmerlen - 1 + bnt.NumBaseInUint64 - 1) / bnt.NumBaseInUint64
	//LenDBGNodeInfo := NBntUint64Len*8 + bnt.BaseTypeNum*2*4
	if NBntUint64Len > NODEMAP_KEY_LEN {
		log.Fatalf("[CDBG] nodeMap just allow max kmerlen is: %d, kmer set: %d\n", NODEMAP_KEY_LEN*32, cf.Kmerlen)
	}
	nodeMap := make(map[[NODEMAP_KEY_LEN]uint64]DBGNode)
	nodeID := constructNodeMap(complexKmerbrfn, nodeMap, NBntUint64Len)
	fmt.Printf("[CDBG] assgin nodeID to : %d\n", nodeID)
	// parallel generate edges and write to file
	edgebrfn := prefix + ".edges.fq.br"
	//numCPU = 1
	newNodeID, edgeID := GenerateDBGEdges(nodeMap, cf, edgebrfn, numCPU, nodeID, uint16(MinKmerFreq), MaxNGSReadLen)
	DBGStatfn := prefix + ".DBG.stat"
	DBGStatWriter(DBGStatfn, newNodeID, edgeID)
	// write DBG nodesArr to the file
	nodesArr := make([]DBGNode, newNodeID)
	NodeMap2NodeArr(nodeMap, nodesArr)
	nodesfn := prefix + ".nodes.Arr.br"
	fc := make(chan int, 1)
	NodesArrWriter(nodesArr, nodesfn, fc)
	fmt.Printf("[CDBG] Generate DBG edges used : %v\n", time.Now().Sub(t0))
}

func ParseEdge(ri constructcf.ReadInfo) (edge DBGEdge) {
	edge.ID = DBG_MAX_INT(ri.ID)
	flist := strings.Fields(string(ri.Anotition))
	if len(flist) > 0 {
		id, err2 := strconv.Atoi(flist[0])
		if err2 != nil {
			log.Fatalf("[ParseEdge]eID:%d edge StartID: %v not digits, please convert to digits...\n", edge.ID, flist[0])
		}
		edge.StartNID = DBG_MAX_INT(id)
	}
	if len(flist) > 1 {
		id, err2 := strconv.Atoi(flist[1])
		if err2 != nil {
			log.Fatalf("[ParseEdge]eID:%d edge EndID: %v not digits, please convert to digits...\n", edge.ID, flist[1])
		}
		edge.EndNID = DBG_MAX_INT(id)
	}

	if len(flist) > 2 && len(flist[2]) > 5 {
		s := flist[2][5:]
		list1 := strings.Split(s, ",")
		edge.NGSPathArr = make([][]EdgeFreq, len(list1)-1)
		for i, path := range list1[:len(list1)-1] {
			if len(path) <= 1 {
				continue
			}
			list2 := strings.Split(path[:len(path)-1], "-")
			ep := make([]EdgeFreq, len(list2))
			for j, id := range list2 {
				eID, err2 := strconv.Atoi(id)
				if err2 != nil {
					log.Fatalf("[ParseEdge]eID:%d path EID: %v not digits, string:%v len(list2):%d please convert to digits...\n", edge.ID, id, list2, len(list2))
				}
				ep[j].ID = DBG_MAX_INT(eID)
			}
			edge.NGSPathArr[i] = ep
		}
		//fmt.Printf("[ParseEdge]NGSPathArr: %v\n", edge.NGSPathArr)
	}
	if len(flist) > 3 {

	}

	//edge.Utg.Ks = make([]byte, len(ri.Seq))
	//copy(edge.Utg.Ks, ri.Seq)
	edge.Utg.Ks = ri.Seq
	edge.Utg.Kq = ri.Qual

	if len(edge.Utg.Ks) != len(edge.Utg.Kq) {
		log.Fatalf("[ParseEdge]eID:%d len(edge.Utg.Ks): %v != len(edge.Utg.Kq): %v\n", edge.ID, len(edge.Utg.Ks), len(edge.Utg.Ks))
	}

	return
}

func ReadEdgesFromFile(edgesbrfn string, edgesSize DBG_MAX_INT) (edgesArr []DBGEdge) {
	edgesArr = make([]DBGEdge, edgesSize)
	format := constructcf.GetReadsFileFormat(edgesbrfn)
	bufSize := (1 << 20)
	size := 1000
	fncs := make(chan []byte, 2)
	recordChan := make(chan constructcf.ReadInfo, size)

	go constructcf.ReadBrFile(edgesbrfn, fncs, bufSize)
	go constructcf.GetReadFileRecord(fncs, recordChan, format, bufSize, true)

	var edgesNum int
	for {
		ri, ok := <-recordChan
		if !ok {
			break
		}
		if len(ri.Seq) == 0 {
			continue
		}
		edge := ParseEdge(ri)
		//fmt.Printf("[ReadEdgesFromFile]edge : %v, ks len:%d\n", edge, len(edge.Utg.Ks))
		edgesArr[edge.ID] = edge
		edgesNum++
	}

	fmt.Printf("[ReadEdgesFromFile] found edge number is : %v\n", edgesNum)
	return
}

func GetRCUnitig(u Unitig) (ru Unitig) {
	ru.Ks = GetReverseCompByteArr(u.Ks)
	ru.Kq = make([]uint8, len(u.Kq))
	copy(ru.Kq, u.Kq)
	ReverseUint8Arr(ru.Kq)
	return
}

/*func RevNode(node DBGNode, kmerlen int) DBGNode {
	rnode := node
	var nBnt constructcf.KmerBnt
	nBnt.Seq = rnode.Seq
	nBnt.Len = kmerlen - 1
	rs := constructcf.ReverseComplet(nBnt)
	rnode.Seq = rs.Seq
	for i := 0; i < bnt.BaseTypeNum; i++ {
		rnode.EdgeIDIncoming[i] = node.EdgeIDOutcoming[bnt.BntRev[i]]
		rnode.EdgeIDOutcoming[i] = node.EdgeIDIncoming[bnt.BntRev[i]]
		// rnode.EdgeIDOutcoming[bnt.BntRev[i]] = node.EdgeIDOutcoming[bnt.BaseTypeNum-1-i], node.EdgeIDIncoming[i]
	}

	return rnode
}*/

/*func ReconstructConsistenceDBG(nodeMap map[string]DBGNode, edgesArr []DBGEdge) {
	for k, v := range nodeMap {
		stk := list.New()
		if v.GetProcessFlag() == 0 {
			v.SetProcessFlag()
			// fmt.Printf("[ReconstructConsistenceDBG] v: %v\n", v)
			nodeMap[k] = v
			stk.PushBack(v)
			for stk.Len() > 0 {
				// Pop a element from stack
				e := stk.Back()
				stk.Remove(e)
				node := e.Value.(DBGNode)
				// fmt.Printf("[ReconstructConsistenceDBG] Pop node: %v\n", node)
				// Processed flag
				if node.GetProcessFlag() != 1 {
					log.Fatalf("[ReconstructConsistenceDBG] node have not been set processed flag, node: %v\n", node)
				}
				for i := 0; i < bnt.BaseTypeNum; i++ {
					if node.EdgeIDOutcoming[i] > 0 {
						eid := node.EdgeIDOutcoming[i]
						if edgesArr[eid].GetProcessFlag() == 0 {
							if edgesArr[eid].StartNID != node.ID {
								// Debug code start
								// if edgesArr[eid].EndNID != node.ID {
								// 	log.Fatalf("[ReconstructConsistenceDBG] edgesArr[eid].EndNID != node.ID\n")
								// }
								// Debug code end
								// fmt.Printf("[ReconstructConsistenceDBG] before RCEdge edge: %v\n", edgesArr[eid])
								RCEdge(edgesArr, eid)
								// fmt.Printf("[ReconstructConsistenceDBG] after RCEdge edge: %v\n", edgesArr[eid])
							}
							if edgesArr[eid].EndNID > 0 {
								var tBnt constructcf.ReadBnt
								tBnt.Seq = edgesArr[eid].Utg.Ks[len(edgesArr[eid].Utg.Ks)-Kmerlen+1:]

								ks := constructcf.GetReadBntKmer(tBnt, 0, Kmerlen-1)
								rs := constructcf.ReverseComplet(ks)
								min := ks
								if min.BiggerThan(rs) {
									min = rs
								}
								if v2, ok := nodeMap[string(min.Seq)]; ok {
									var v2Bnt constructcf.ReadBnt
									v2Bnt.Seq = v2.Seq
									v2Bnt.Length = Kmerlen - 1
									if v2.GetProcessFlag() == 1 {
										if reflect.DeepEqual(v2Bnt, ks) == false {
											// log.Fatalf("[ReconstructConsistenceDBG] found not consistence node\n")
											fmt.Printf("[ReconstructConsistenceDBG] found not consistence node, edge: %v\nv1: %v\nv2: %v\n", edgesArr[eid], node, v2)
											if eid == 2870 {
												fmt.Printf("[ReconstructConsistenceDBG] edge: %v\n", edgesArr[8014])
											}
											fmt.Printf("[ReconstructConsistenceDBG] edge start: %v\nedge end: %v\n", edgesArr[eid].Utg.Ks[:Kmerlen], edgesArr[eid].Utg.Ks[len(edgesArr[eid].Utg.Ks)-Kmerlen:])
											var v1Bnt constructcf.ReadBnt
											v1Bnt.Seq = node.Seq
											v1Bnt.Length = Kmerlen - 1
											fmt.Printf("[ReconstructConsistenceDBG] v1.Seq: %v\n", constructcf.ExtendReadBnt2Byte(v1Bnt))
											fmt.Printf("[ReconstructConsistenceDBG] v2.Seq: %v\n", constructcf.ExtendReadBnt2Byte(v2Bnt))
										}
									} else {
										if reflect.DeepEqual(v2Bnt, ks) == false {
											v2 = RevNode(v2)
										}
										v2.SetProcessFlag()
										nodeMap[string(min.Seq)] = v2
										stk.PushBack(v2)

									}

								} else {
									log.Fatalf("[ReconstructConsistenceDBG] not found edge' end node, edge: %v\n", edgesArr[eid])
								}
							}
							edgesArr[eid].SetProcessFlag()
						}
					}

					if node.EdgeIDIncoming[i] > 0 {
						eid := node.EdgeIDIncoming[i]
						if edgesArr[eid].GetProcessFlag() == 0 {
							if edgesArr[eid].EndNID != node.ID {
								RCEdge(edgesArr, eid)
							}
							if edgesArr[eid].StartNID > 0 {
								var tBnt constructcf.ReadBnt
								tBnt.Seq = edgesArr[eid].Utg.Ks[:Kmerlen-1]
								ks := constructcf.GetReadBntKmer(tBnt, 0, Kmerlen-1)
								rs := constructcf.ReverseComplet(ks)
								min := ks
								if ks.BiggerThan(rs) {
									min = rs
								}

								if v2, ok := nodeMap[string(min.Seq)]; ok {
									var v2Bnt constructcf.ReadBnt
									v2Bnt.Seq = v2.Seq
									v2Bnt.Length = Kmerlen - 1
									if v2.GetProcessFlag() == 1 {
										if reflect.DeepEqual(v2Bnt, ks) == false {
											// log.Fatalf("[ReconstructConsistenceDBG] found not consistence node\n")
											fmt.Printf("[ReconstructConsistenceDBG] found not consistence node, edge: %v\nv1: %v\nv2: %v\n", edgesArr[eid], node, v2)
										}
									} else {
										if reflect.DeepEqual(v2Bnt, ks) == false {
											v2 = RevNode(v2)
										}
										v2.SetProcessFlag()
										nodeMap[string(min.Seq)] = v2
										stk.PushBack(v2)
									}
								} else {
									log.Fatalf("[ReconstructConsistenceDBG] not found edge' start node, edge: %v\n", edgesArr[eid])
								}
							}
							edgesArr[eid].SetProcessFlag()
						}
					}
				}
			}
		}
	}
}*/

/*func Coming2String(coming [bnt.BaseTypeNum]DBG_MAX_INT) (cs string) {
	for _, v := range coming {
		cs += " " + strconv.Itoa(int(v))
	}
	return
}*/

func GraphvizDBGArr(nodesArr []DBGNode, edgesArr []DBGEdge, graphfn string) {
	// create a new graph
	g := gographviz.NewGraph()
	g.SetName("G")
	g.SetDir(true)
	g.SetStrict(false)
	for _, v := range nodesArr {
		if v.ID < 2 || v.GetDeleteFlag() > 0 {
			continue
		}
		attr := make(map[string]string)
		attr["color"] = "Green"
		attr["shape"] = "record"
		var labels string
		//labels = "{<f0>" + strconv.Itoa(int(v.EdgeIDIncoming[0])) + "|<f1>" + strconv.Itoa(int(v.EdgeIDIncoming[1])) + "|<f2>" + strconv.Itoa(int(v.EdgeIDIncoming[2])) + "|<f3>" + strconv.Itoa(int(v.EdgeIDIncoming[3])) + "}|" + strconv.Itoa(int(v.ID)) + "|{<f0>" + strconv.Itoa(int(v.EdgeIDOutcoming[0])) + "|<f1>" + strconv.Itoa(int(v.EdgeIDOutcoming[1])) + "|<f2>" + strconv.Itoa(int(v.EdgeIDOutcoming[2])) + "|<f3>" + strconv.Itoa(int(v.EdgeIDOutcoming[3])) + "}"
		labels = "\"{" + strconv.Itoa(int(v.EdgeIDIncoming[0])) +
			"|" + strconv.Itoa(int(v.EdgeIDIncoming[1])) +
			"|" + strconv.Itoa(int(v.EdgeIDIncoming[2])) +
			"|" + strconv.Itoa(int(v.EdgeIDIncoming[3])) +
			"}|" + strconv.Itoa(int(v.ID)) +
			"| {" + strconv.Itoa(int(v.EdgeIDOutcoming[0])) +
			"|" + strconv.Itoa(int(v.EdgeIDOutcoming[1])) +
			"|" + strconv.Itoa(int(v.EdgeIDOutcoming[2])) +
			"|" + strconv.Itoa(int(v.EdgeIDOutcoming[3])) + "}\""
		attr["label"] = labels
		g.AddNode("G", strconv.Itoa(int(v.ID)), attr)
	}
	g.AddNode("G", "0", nil)
	//fmt.Printf("[GraphvizDBGArr] finished Add Nodes\n")

	for i := 1; i < len(edgesArr); i++ {
		e := edgesArr[i]
		if e.ID < 2 || e.GetDeleteFlag() > 0 {
			continue
		}
		attr := make(map[string]string)
		attr["color"] = "Blue"
		//labels := strconv.Itoa(int(e.ID)) + "len" + strconv.Itoa(len(e.Utg.Ks))
		//labels := strconv.Itoa(int(e.ID))
		labels := "\"ID:" + strconv.Itoa(int(e.ID)) + " len:" + strconv.Itoa(len(e.Utg.Ks)) + "\""
		attr["label"] = labels
		g.AddEdge(strconv.Itoa(int(e.StartNID)), strconv.Itoa(int(e.EndNID)), true, attr)
	}
	//fmt.Printf("[GraphvizDBGArr] finished Add edges\n")
	// output := graph.String()
	gfp, err := os.Create(graphfn)
	if err != nil {
		log.Fatalf("[GraphvizDBG] Create file: %s failed, err: %v\n", graphfn, err)
	}
	defer gfp.Close()
	gfp.WriteString(g.String())
}

func IsContainCycleEdge(nd DBGNode) bool {
	var arrIn, arrOut []DBG_MAX_INT
	for _, eID := range nd.EdgeIDIncoming {
		if eID > 1 {
			arrIn = append(arrIn, eID)
		}
	}
	for _, eID := range nd.EdgeIDOutcoming {
		if eID > 1 {
			arrOut = append(arrOut, eID)
		}
	}

	for _, eID1 := range arrIn {
		for _, eID2 := range arrOut {
			if eID1 == eID2 {
				return true
			}
		}
	}

	return false
}

func GetEdgeIDComing(coming [bnt.BaseTypeNum]DBG_MAX_INT) (num int, edgeID DBG_MAX_INT) {
	for _, v := range coming {
		if v > 1 {
			num++
			edgeID = v
		}
	}

	return
}

func ConcatEdges(u1, u2 Unitig, kmerlen int) (u Unitig) {
	u.Ks = make([]byte, len(u1.Ks)+len(u2.Ks)-(kmerlen-1))
	u.Kq = make([]uint8, len(u1.Kq)+len(u2.Kq)-(kmerlen-1))
	if !reflect.DeepEqual(u1.Ks[len(u1.Ks)-kmerlen+1:], u2.Ks[:kmerlen-1]) {
		fmt.Fprintf(os.Stderr, "[ConcatEdges] u1:%s, u2:%s can not concatenatable\n", Transform2Char(u1.Ks[len(u1.Ks)-kmerlen+1:]), Transform2Char(u2.Ks[:kmerlen-1]))
		u.Ks = u.Ks[:0]
		u.Kq = u.Kq[:0]
		return
	}
	if len(u1.Ks) != len(u1.Kq) {
		fmt.Fprintf(os.Stderr, "[ConcatEdges] len(u1.Ks):%d != len(u1.Kq):%d\n", len(u1.Ks), len(u1.Kq))
		u.Ks = u.Ks[:0]
		u.Kq = u.Kq[:0]
		return
	}
	if len(u2.Ks) != len(u2.Kq) {
		fmt.Fprintf(os.Stderr, "[ConcatEdges] len(u2.Ks):%d != len(u2.Kq):%d\n", len(u2.Ks), len(u2.Kq))
		u.Ks = u.Ks[:0]
		u.Kq = u.Kq[:0]
		return
	}

	copy(u.Ks[:len(u1.Ks)], u1.Ks)
	copy(u.Ks[len(u1.Ks):], u2.Ks[kmerlen-1:])

	copy(u.Kq[:len(u1.Kq)], u1.Kq)
	for i := 0; i < kmerlen-1; i++ {
		u.Kq[len(u1.Kq)-(kmerlen-1)+i] += u2.Kq[i]
	}
	copy(u.Kq[len(u1.Kq):], u2.Kq[kmerlen-1:])

	return
}

/*func ConcatEdges(edgesArr []DBGEdge, inID, outID, dstID DBG_MAX_INT) {
	// check if is connective
	var inBnt, outBnt constructcf.ReadBnt
	fmt.Printf("[ConcatEdges] inID: %d, outID: %d, dstID: %d\nedgesArr[inID]:%v\nedgesArr[outID]:%v\n", inID, outID, dstID, edgesArr[inID], edgesArr[outID])
	inBnt.Seq = edgesArr[inID].Utg.Ks[len(edgesArr[inID].Utg.Ks)-Kmerlen+1:]
	inBnt.Length = len(inBnt.Seq)
	outBnt.Seq = edgesArr[outID].Utg.Ks[:Kmerlen-1]
	outBnt.Length = len(outBnt.Seq)
	//fmt.Printf("[ConcatEdges] inID: %d, outID: %d, dstID: %d\nseq1:%v\nseq2:%v\n", inID, outID, dstID, inBnt, outBnt)
	if reflect.DeepEqual(inBnt.Seq, outBnt.Seq) == false {
		log.Fatalf("[ConcatEdges] two edges is not connective\n\tin: %v\n\tout: %v\n", edgesArr[inID], edgesArr[outID])
	}
	if dstID == inID {
		u1, u2 := edgesArr[inID].Utg, edgesArr[outID].Utg
		edgesArr[inID].Utg.Ks = append(u1.Ks, u2.Ks[Kmerlen-1:]...)
		for i := 0; i < Kmerlen-1; i++ {
			if u1.Kq[len(u1.Kq)-Kmerlen+1+i] < u2.Kq[i] {
				u1.Kq[len(u1.Kq)-Kmerlen+1+i] = u2.Kq[i]
			}
		}
		edgesArr[inID].Utg.Kq = append(u1.Kq, u2.Kq[Kmerlen-1:]...)
		// DeleteEdgeID(nodesArr, edgesArr[inID].EndNID, inID)
		edgesArr[inID].EndNID = edgesArr[outID].EndNID

	} else {
		u1, u2 := edgesArr[inID].Utg, edgesArr[outID].Utg
		seq := make([]byte, len(u1.Ks))
		copy(seq, u1.Ks)
		edgesArr[outID].Utg.Ks = append(seq, u2.Ks[Kmerlen-1:]...)
		qul := make([]uint8, len(u1.Kq))
		copy(qul, u1.Kq)
		for i := 0; i < Kmerlen-1; i++ {
			if u1.Kq[len(u1.Kq)-Kmerlen+1+i] < u2.Kq[i] {
				qul[len(qul)-Kmerlen+1+i] = u2.Kq[i]
			}
		}
		edgesArr[outID].Utg.Kq = append(qul, u2.Kq[Kmerlen-1:]...)
		// DeleteEdgeID(nodesArr, edgesArr[outID].StartNID, outID)
		edgesArr[outID].StartNID = edgesArr[inID].StartNID
	}
}*/

/*func substituteEdgeID(nodeMap map[[NODEMAP_KEY_LEN]uint64]DBGNode, nodekey []uint64, srcID, dstID DBG_MAX_INT, kmerlen int) bool {
	var nkB constructcf.KmerBnt
	nkB.Seq = nodekey
	ks := constructcf.GetReadBntKmer(nkB, 0, kmerlen-1)
	rs := constructcf.ReverseComplet(ks)
	min := ks
	if ks.BiggerThan(rs) {
		min = rs
	}
	suc := false
	if nv, ok := nodeMap[string(min.Seq)]; ok {
		for i := 0; i < bnt.BaseTypeNum; i++ {
			if nv.EdgeIDIncoming[i] == srcID {
				nv.EdgeIDIncoming[i] = dstID
				suc = true
				break
			}
			if nv.EdgeIDOutcoming[i] == srcID {
				nv.EdgeIDOutcoming[i] = dstID
				suc = true
				break
			}
		}
		if suc == true {
			nodeMap[string(min.Seq)] = nv
		}
	} else {
		log.Fatalf("[substituteEdgeID] not found correct node\n")
	}

	return suc

}*/

/*func ChangeIDMapPathAndJoinPathArr(IDMapPath map[DBG_MAX_INT]uint32, joinPathArr []Path, e1, e2 DBGEdge, v DBGNode) []Path {
	var p, p1, p2 Path
	if idx, ok := IDMapPath[e1.ID]; ok {
		p1.IDArr = append(p1.IDArr, joinPathArr[idx].IDArr...)
		joinPathArr[idx].IDArr = nil
	} else {
		p1.IDArr = append(p1.IDArr, e1.ID)
	}

	if idx, ok := IDMapPath[e2.ID]; ok {
		p2.IDArr = append(p2.IDArr, joinPathArr[idx].IDArr...)
		joinPathArr[idx].IDArr = nil
	} else {
		p2.IDArr = append(p2.IDArr, e2.ID)
	}

	if e1.EndNID == v.ID {
		p.IDArr = append(p.IDArr, p1.IDArr...)
		tp := p2
		if e2.EndNID == v.ID {
			tp.IDArr = GetReverseDBG_MAX_INTArr(p2.IDArr)
		}
		p.IDArr = append(p.IDArr, tp.IDArr...)
	} else {
		tp := p2
		if e2.StartNID == v.ID {
			tp.IDArr = GetReverseDBG_MAX_INTArr(p2.IDArr)
		}
		p.IDArr = append(p.IDArr, tp.IDArr...)
		p.IDArr = append(p.IDArr, p1.IDArr...)
	}
	joinPathArr = append(joinPathArr, p)
	ni := uint32(len(joinPathArr)-1)
	IDMapPath[e1.ID] = ni
	IDMapPath[p.IDArr[0]]

	return joinPathArr
} */

func CleanDBGNodeEdgeIDComing(nodesArr []DBGNode, nID DBG_MAX_INT) {
	for i := 0; i < bnt.BaseTypeNum; i++ {
		if nodesArr[nID].EdgeIDIncoming[i] == 1 {
			nodesArr[nID].EdgeIDIncoming[i] = 0
		}
		if nodesArr[nID].EdgeIDOutcoming[i] == 1 {
			nodesArr[nID].EdgeIDOutcoming[i] = 0
		}
	}
}

func CleanDBGEdgeIDComing(nodesArr []DBGNode) {
	for i, v := range nodesArr {
		if i < 2 || v.GetDeleteFlag() > 0 {
			continue
		}
		CleanDBGNodeEdgeIDComing(nodesArr, v.ID)
	}
}

/*func MakeSelfCycleEdgeOutcomingToIncoming(nodesArr []DBGNode, edgesArr []DBGEdge, opt Options) {
	for i, e := range edgesArr {
		if i < 2 || e.GetDeleteFlag() > 0 {
			continue
		}
		if e.StartNID == e.EndNID { // self cycle edge
			nd := nodesArr[e.StartNID]
			bn := constructcf.GetReadBntKmer(e.Utg.Ks, 0, opt.Kmer-1)
			if !reflect.DeepEqual(bn.Seq, nd.Seq) {
				e.Utg.Ks = GetReverseCompByteArr(e.Utg.Ks)
				ReverseByteArr(e.Utg.Kq)
			}
			if (nd.EdgeIDOutcoming[e.Utg.Ks[opt.Kmer-1]] != e.ID) || (nd.EdgeIDIncoming[e.Utg.Ks[len(e.Utg.Ks)-opt.Kmer]] != e.ID) {
				log.Fatalf("[MakeSelfCycleEdgeOutcomingToIncoming] error cycle edge set, e: %v\n\tv: %v\n", e, nd)
			}

		}
	}
}*/

func SetEdgeID(nodesArr []DBGNode, edgesArr []DBGEdge) []DBGEdge {
	nEA := make([]DBGEdge, 0, len(edgesArr)*2/3)
	edgeID := DBG_MAX_INT(0)
	for i, e := range edgesArr {
		if i < 2 {
			if len(nEA) != int(edgeID) {
				log.Fatalf("[SetEdgeID]len(nEA):%d != edgeID:%d\n", len(nEA), edgeID)
			}
			nEA = append(nEA, e)
			edgeID++
		}
		if e.GetDeleteFlag() > 0 {
			continue
		}
		if e.StartNID > 0 {
			if !SubstituteEdgeID(nodesArr, e.StartNID, e.ID, edgeID) {
				log.Fatalf("[SetEdgeID]v: %v\ne.ID: %v substitute by :%d failed\n", nodesArr[e.StartNID], e.ID, edgeID)
			}
		}
		if e.EndNID > 0 {
			if !SubstituteEdgeID(nodesArr, e.EndNID, e.ID, edgeID) {
				if e.StartNID != e.EndNID {
					log.Fatalf("[SetEdgeID]v: %v\ne.ID: %v substitute by :%d failed\n", nodesArr[e.EndNID], e.ID, edgeID)
				}
			}
		}
		e.ID = edgeID
		if len(nEA) != int(edgeID) {
			log.Fatalf("[SetEdgeID]len(nEA):%d != edgeID:%d\n", len(nEA), edgeID)
		}
		nEA = append(nEA, e)
		edgeID++
	}

	fmt.Printf("[SetEdgeID]len(edgeArr):%d edgeID:%v len(nEA):%d\n", len(edgesArr), edgeID, len(nEA))

	return nEA
}

func SetNodeID(nodesArr []DBGNode, edgesArr []DBGEdge) []DBGNode {
	nNA := make([]DBGNode, 0, len(nodesArr)*3/5)
	nodeID := DBG_MAX_INT(0)
	for i, n := range nodesArr {
		if i < 2 {
			if len(nNA) != int(nodeID) {
				log.Fatalf("[SetNodeID]len(nNA):%d != nodeID:%d\n", len(nNA), nodeID)
			}
			nNA = append(nNA, n)
			nodeID++
		}
		if n.GetDeleteFlag() > 0 {
			continue
		}
		if CountNumEdgeIDComing(n.EdgeIDIncoming) == 0 && CountNumEdgeIDComing(n.EdgeIDOutcoming) == 0 {
			continue
		}
		for j := 0; j < bnt.BaseTypeNum; j++ {
			eID1, eID2 := n.EdgeIDIncoming[j], n.EdgeIDOutcoming[j]
			if eID1 > 0 {
				if edgesArr[eID1].StartNID == n.ID {
					if edgesArr[eID1].StartNID == edgesArr[eID1].EndNID {
						edgesArr[eID1].EndNID = nodeID
					}
					edgesArr[eID1].StartNID = nodeID
				} else if edgesArr[eID1].EndNID == n.ID {
					if edgesArr[eID1].StartNID == edgesArr[eID1].EndNID {
						edgesArr[eID1].StartNID = nodeID
					}
					edgesArr[eID1].EndNID = nodeID
				} else {
					if edgesArr[eID1].StartNID != edgesArr[eID1].EndNID {
						log.Fatalf("[SetNodeID]e.ID:%v StartNID:%v,EndNID:%v substitute by :%v failed\n", eID1, edgesArr[eID1].StartNID, edgesArr[eID1].EndNID, nodeID)
					}
				}
			}
			if eID2 > 0 {
				if edgesArr[eID2].StartNID == n.ID {
					if edgesArr[eID2].StartNID == edgesArr[eID2].EndNID {
						edgesArr[eID2].EndNID = nodeID
					}
					edgesArr[eID2].StartNID = nodeID
				} else if edgesArr[eID2].EndNID == n.ID {
					if edgesArr[eID2].StartNID == edgesArr[eID2].EndNID {
						edgesArr[eID2].StartNID = nodeID
					}
					edgesArr[eID2].EndNID = nodeID
				} else {
					if edgesArr[eID2].StartNID != edgesArr[eID2].EndNID {
						log.Fatalf("[SetNodeID]e.ID:%v StartNID:%v,EndNID:%v substitute by :%v failed\n", eID2, edgesArr[eID2].StartNID, edgesArr[eID2].EndNID, nodeID)
					}
				}
			}
		}
		n.ID = nodeID
		if len(nNA) != int(nodeID) {
			log.Fatalf("[SetNodeID]len(nNA):%d != nodeID:%d\n", len(nNA), nodeID)
		}
		nNA = append(nNA, n)
		nodeID++
	}

	fmt.Printf("[SetNodeID]len(nodesArr):%d nodeID:%v len(nNA):%d\n", len(nodesArr), nodeID, len(nNA))

	return nNA
}

func CheckNGSPath(edgesArr []DBGEdge, nodesArr []DBGNode, arr []EdgeFreq) bool {
	idx := -1
	for i, ef := range arr {
		if edgesArr[ef.ID].GetDeleteFlag() > 0 {
			return false
		}
		if edgesArr[ef.ID].StartNID == edgesArr[ef.ID].EndNID {
			continue
		}
		if edgesArr[ef.ID].GetTwoEdgesCycleFlag() > 0 || IsTwoEdgesCyclePath(edgesArr, nodesArr, ef.ID) {
			continue
		}
		idx = i
		break
	}

	if idx < 0 {
		return false
	}
	e := &edgesArr[arr[idx].ID]
	if idx > 0 {
		e2 := &edgesArr[arr[idx-1].ID]
		nID := GetShareNodeID(e, e2)
		var strand bool
		if nID == e.StartNID {
			strand = PLUS
		} else {
			strand = MINUS
		}
		le := e
		for i := idx - 1; i >= 0; i-- {
			if edgesArr[arr[i].ID].GetDeleteFlag() > 0 {
				return false
			}
			eIDArr := GetNextEdgeIDArr(le, nodesArr, strand, BACKWARD)
			if !IsInDBG_MAX_INTArr(eIDArr, arr[i].ID) {
				return false
			}
			ne := &edgesArr[arr[i].ID]
			strand = GetNextEdgeStrand(le, ne, nodesArr, strand)
			le = ne
		}
	}

	if idx < len(arr)-1 {
		e2 := &edgesArr[arr[idx+1].ID]
		nID := GetShareNodeID(e, e2)
		var strand bool
		if nID == e.EndNID {
			strand = PLUS
		} else {
			strand = MINUS
		}
		le := e
		for i := idx + 1; i < len(arr); i++ {
			if edgesArr[arr[i].ID].GetDeleteFlag() > 0 {
				return false
			}
			eIDArr := GetNextEdgeIDArr(le, nodesArr, strand, FORWARD)
			if !IsInDBG_MAX_INTArr(eIDArr, arr[i].ID) {
				return false
			}
			ne := &edgesArr[arr[i].ID]
			strand = GetNextEdgeStrand(le, ne, nodesArr, strand)
			le = ne
		}
	}

	return true

}

func SmfyDBG(nodesArr []DBGNode, edgesArr []DBGEdge, opt Options, DeleteLowFreqEdgeFlag bool) {
	kmerlen := opt.Kmer
	deleteNodeNum, deleteEdgeNum := 0, 0
	longTipsEdgesNum := 0
	allKmerFreq := 0
	kmerNum := 0

	// clean DBG EdgeIDComing
	CleanDBGEdgeIDComing(nodesArr)

	// delete maybe short repeat edge than small than opt.MaxNGSReadLen
	for i := range edgesArr {
		if i < 2 {
			continue
		}
		e := &edgesArr[i]

		if e.GetDeleteFlag() > 0 {
			fmt.Printf("[SmfyDBG]has been deleted edge ID:%d el:%d\n", i, e.GetSeqLen())
			if e.StartNID > 0 {
				SubstituteEdgeID(nodesArr, e.StartNID, e.ID, 0)
			}
			if e.EndNID > 0 {
				SubstituteEdgeID(nodesArr, e.EndNID, e.ID, 0)
			}
			continue
		}

		var deleteFlag bool
		// remove short tips
		if e.StartNID == 0 && e.EndNID == 0 && len(e.Utg.Ks) < opt.MaxNGSReadLen {
			deleteFlag = true
		} else if e.StartNID == 0 && e.EndNID > 0 && len(e.Utg.Ks) < opt.MaxNGSReadLen {
			deleteFlag = true
			if !SubstituteEdgeID(nodesArr, e.EndNID, e.ID, 0) {
				fmt.Fprintf(os.Stderr, "[SmfyDBG]v: %v\ne.ID: %v substitute by 0 failed\n", nodesArr[e.EndNID], e.ID)
			}
			//fmt.Printf("[SmfyDBG]delete e.ID : %v, v: %v\n", e.ID, nodesArr[e.EndNID])
		} else if e.EndNID == 0 && e.StartNID > 0 && len(e.Utg.Ks) < opt.MaxNGSReadLen {
			deleteFlag = true
			if !SubstituteEdgeID(nodesArr, e.StartNID, e.ID, 0) {
				fmt.Fprintf(os.Stderr, "[SmfyDBG]v: %v\ne.ID: %v substitute by 0 failed\n", nodesArr[e.StartNID], e.ID)
			}
		}

		if deleteFlag {
			edgesArr[i].SetDeleteFlag()
			deleteEdgeNum++
			fmt.Printf("[SmfyDBG]delete tip edge ID:%d el:%d\n", i, e.GetSeqLen())
			continue
		} else {
			longTipsEdgesNum++
		}

		// remove low freq edge
		if DeleteLowFreqEdgeFlag {
			var totalFreq int
			for j := 0; j < e.GetSeqLen(); j++ {
				totalFreq += int(e.Utg.Kq[j])
			}
			allKmerFreq += totalFreq
			kmerNum += (e.GetSeqLen() - (kmerlen - 1))
			avgFreq := totalFreq / (e.GetSeqLen() - (kmerlen - 1))
			if avgFreq < opt.MinMapFreq {
				edgesArr[i].SetDeleteFlag()
				deleteEdgeNum++
				//fmt.Printf("[SmfyDBG]delete low freq:%d edge ID:%d el:%d\n",avgFreq, i, e.GetSeqLen())
				if e.StartNID > 0 {
					if !SubstituteEdgeID(nodesArr, e.StartNID, e.ID, 0) {
						fmt.Fprintf(os.Stderr, "[SmfyDBG]v: %v\ne.ID:%d substitute by 0 failed\n", nodesArr[e.StartNID], e.ID)
					}
				}
				if e.EndNID > 0 {
					if !SubstituteEdgeID(nodesArr, e.EndNID, e.ID, 0) {
						fmt.Fprintf(os.Stderr, "[SmfyDBG]v:%v\ne.ID:%d substitute by 0 failed\n", nodesArr[e.EndNID], e.ID)
					}
				}
				fmt.Fprintf(os.Stderr, "[SmfyDBG]delete low kmer freq eID:%d, avgFreq:%d\n", e.ID, avgFreq)
			}
		}
	}

	// remove samll cycle maybe repeat
	/*if e.StartNID > 0 && e.StartNID == e.EndNID && len(e.Utg.Ks) < opt.MaxNGSReadLen {
		v := nodesArr[e.StartNID]
		inNum, inID := GetEdgeIDComing(v.EdgeIDIncoming)
		outNum, outID := GetEdgeIDComing(v.EdgeIDOutcoming)
		if inNum == 1 && outNum == 1 && inID == outID {
			edgesArr[i].SetDeleteFlag()
			deleteEdgeNum++
			nodesArr[v.ID].SetDeleteFlag()
			deleteNodeNum++
		}
	}*/

	mergePathMap := make(map[DBG_MAX_INT][]EdgeFreq, 1000)

	for i, v := range nodesArr {
		if i < 2 {
			continue
		}
		if v.GetDeleteFlag() > 0 {
			fmt.Printf("[SmfyDBG]has deleted node ID:%d\n", i)
			continue
		}
		if IsContainCycleEdge(v) {
			//fmt.Printf("[SmfyDBG]node ID: %v is contain Cycle Edge\n", v.ID)
			continue
		}
		//fmt.Printf("[SmfyDBG] v: %v\n", v)
		inNum, inID := GetEdgeIDComing(v.EdgeIDIncoming)
		outNum, outID := GetEdgeIDComing(v.EdgeIDOutcoming)
		if inNum == 0 && outNum == 0 {
			nodesArr[i].SetDeleteFlag()
			deleteNodeNum++
		} else if inNum+outNum == 1 {
			/*id := inID
			if outNum == 1 {
				id = outID
			}
			//fmt.Printf("[SmfyDBG]v: %v,id: %v\n", v, id)
			//fmt.Printf("[SmfyDBG]edgesArr[%v]: %v\n",id, edgesArr[id])
			//if edgesArr[id].StartNID == v.ID {
			//	edgesArr[id].StartNID = 0
			//} else {
			//	edgesArr[id].EndNID = 0
			//}
			if len(edgesArr[id].Utg.Ks) < opt.MaxNGSReadLen {
				edgesArr[id].SetDeleteFlag()
				deleteEdgeNum++
				fmt.Printf("[SmfyDBG]delete tip edge ID:%d el:%d\n", i, edgesArr[id].GetSeqLen())
				if edgesArr[id].StartNID > 0 && int(edgesArr[id].StartNID) != i && !SubstituteEdgeID(nodesArr, edgesArr[id].StartNID, id, 0) {
					log.Fatalf("[SmfyDBG]v: %v\ne.ID: %v substitute by 0 failed\n", nodesArr[edgesArr[id].StartNID], id)
				} else if edgesArr[id].EndNID > 0 && int(edgesArr[id].EndNID) != i && !SubstituteEdgeID(nodesArr, edgesArr[id].EndNID, id, 0) {
					log.Fatalf("[SmfyDBG]v: %v\ne.ID: %v substitute by 0 failed\n", nodesArr[edgesArr[id].EndNID], id)
				}
			} else {
				if int(edgesArr[id].StartNID) == i {
					edgesArr[id].StartNID = DBG_MAX_INT(0)
				} else {
					edgesArr[id].EndNID = DBG_MAX_INT(0)
				}
				longTipsEdgesNum++
			}
			nodesArr[i].SetDeleteFlag()
			fmt.Printf("[SmfyDBG]delete tip node ID:%d eID:%d\n", i, id)
			deleteNodeNum++ */
		} else if inNum == 1 && outNum == 1 && inID != outID { // prevent cycle ring
			e1 := edgesArr[inID]
			e2 := edgesArr[outID]
			if e1.StartNID == e1.EndNID || e2.StartNID == e2.EndNID { // if  encounter cycle ring have same EdgeIDcoming
				continue
			}

			nodesArr[v.ID].SetDeleteFlag()
			deleteNodeNum++

			// set merge path
			{
				fmt.Printf("[SmfyDBG]eID1:%d eID2:%d\n", inID, outID)
				PrintEdgeInfo(&e1)
				PrintEdgeInfo(&e2)
				PrintNodeInfo(&v)
				mergePath := make([]EdgeFreq, 2)
				mergePath[0].ID, mergePath[1].ID = e1.ID, e2.ID
				mergePath[0].Freq, mergePath[1].Freq = 10, 10
				efa1, ok1 := mergePathMap[e1.ID]
				efa2, ok2 := mergePathMap[e2.ID]
				if ok1 && ok2 {
					if efa1[len(efa1)-1].ID == e1.ID {
						if efa2[0].ID == e2.ID {
							efa1 = append(efa1, efa2...)
						} else {
							efa1 = append(efa1, GetReverseEdgeFreqArr2(efa2)...)
						}
						delete(mergePathMap, e1.ID)
						delete(mergePathMap, e2.ID)
						mergePathMap[efa1[0].ID] = efa1
						mergePathMap[efa1[len(efa1)-1].ID] = efa1
					} else {
						if efa2[0].ID == e2.ID {
							efa1 = append(GetReverseEdgeFreqArr2(efa2), efa1...)
						} else {
							efa1 = append(efa2, efa1...)
						}
						delete(mergePathMap, e1.ID)
						delete(mergePathMap, e2.ID)
						mergePathMap[efa1[0].ID] = efa1
						mergePathMap[efa1[len(efa1)-1].ID] = efa1
					}
				} else if ok1 && !ok2 {
					if efa1[0].ID == e1.ID {
						efa1 = append(GetReverseEdgeFreqArr2(mergePath), efa1[1:]...)
					} else {
						efa1 = append(efa1, mergePath[1:]...)
					}
					mergePathMap[efa1[0].ID] = efa1
					mergePathMap[efa1[len(efa1)-1].ID] = efa1
					delete(mergePathMap, e1.ID)
				} else if !ok1 && ok2 {
					if efa2[0].ID == e2.ID {
						efa2 = append(mergePath, efa2[1:]...)
					} else {
						efa2 = append(efa2, mergePath[0])
					}
					mergePathMap[efa2[0].ID] = efa2
					mergePathMap[efa2[len(efa2)-1].ID] = efa2
					delete(mergePathMap, e2.ID)
				} else {
					mergePathMap[e1.ID] = mergePath
					mergePathMap[e2.ID] = mergePath
				}
			}

			//fmt.Printf("[SmfyDBG] e1: %v\n\te2: %v\n\tnd: %v\n", e1, e2, v)
			/*u1, u2 := e1.Utg, e2.Utg
			if e1.EndNID == v.ID {
				nID := e2.EndNID
				if e2.StartNID != v.ID {
					u2 = GetRCUnitig(u2)
					nID = e2.StartNID
				}
				//fmt.Printf("[SmfyDBG]v: %v\ne1.ID: %v, e1.StartNID: %v, e1.EndNID: %v, e2.ID:%v, e2.StartNID: %v, e2.EndNID: %v\n", v, e1.ID, e1.StartNID, e1.EndNID, e2.ID, e2.StartNID, e2.EndNID)
				edgesArr[inID].Utg = ConcatEdges(u1, u2, kmerlen)
				edgesArr[inID].EndNID = nID
				if nID > 0 && !SubstituteEdgeID(nodesArr, nID, e2.ID, e1.ID) {
					log.Fatalf("[SmfyDBG]v: %v\ne2.ID: %v substitute by e1.ID: %v failed, node: %v\n", v, e2.ID, e1.ID, nodesArr[nID])
				}
			} else {
				nID := e2.StartNID
				if e2.EndNID != v.ID {
					u2 = GetRCUnitig(u2)
					nID = e2.EndNID
				}
				//fmt.Printf("[SmfyDBG]v: %v\ne1.ID: %v, e1.StartNID: %v, e1.EndNID: %v, e2.ID:%v, e2.StartNID: %v, e2.EndNID: %v\n", v, e1.ID, e1.StartNID, e1.EndNID, e2.ID, e2.StartNID, e2.EndNID)
				edgesArr[inID].Utg = ConcatEdges(u2, u1, kmerlen)
				edgesArr[inID].StartNID = nID
				if nID > 0 && !SubstituteEdgeID(nodesArr, nID, e2.ID, e1.ID) {
					log.Fatalf("[SmfyDBG]v: %v\ne2.ID: %v substitute by e1.ID: %v failed, node: %v\n", v, e2.ID, e1.ID, nodesArr[nID])
				}
			}*/
		}
	}

	// reset DBG
	{
		for k, v := range mergePathMap {
			if k != v[0].ID {
				continue
			}
			fmt.Printf("[SmfyDBG]merged path:%v\n", v)
			var nID DBG_MAX_INT
			if edgesArr[v[0].ID].EndNID == edgesArr[v[1].ID].StartNID || edgesArr[v[0].ID].EndNID == edgesArr[v[1].ID].EndNID {
				nID = edgesArr[v[0].ID].EndNID
			} else {
				nID = edgesArr[v[0].ID].StartNID
			}
			e := ConcatEdgePathUtg(v, nID, edgesArr, nodesArr, kmerlen)
			if e.EndNID > 0 && !SubstituteEdgeID(nodesArr, e.EndNID, v[len(v)-1].ID, e.ID) {
				log.Fatalf("[SmfyDBG]v: %v\ne2.ID: %v substitute by e1.ID: %v failed, node: %v\n", nodesArr[e.EndNID], v[len(v)-1].ID, e.ID, nodesArr[e.EndNID])
			}
			edgesArr[e.ID] = e
			for _, ef := range v[1:] {
				edgesArr[ef.ID].SetDeleteFlag()
				deleteEdgeNum++
			}
		}
	}

	// reset NGSPath
	changedNGSPathCount := 0
	for i := range edgesArr {
		e := &edgesArr[i]
		if i < 2 || e.GetDeleteFlag() > 0 || len(e.NGSPathArr) == 0 {
			continue
		}
		for j := range e.NGSPathArr {
			for x := 0; x < len(e.NGSPathArr[j]); x++ {
				id := e.NGSPathArr[j][x].ID
				if mp, ok := mergePathMap[id]; ok {
					rmp := GetReverseEdgeFreqArr2(mp)
					al := len(e.NGSPathArr[j])
					e.NGSPathArr[j] = RePlaceNGSPath2(e.NGSPathArr[j], mp, rmp)
					if al != len(e.NGSPathArr[j]) {
						x = 0
					}
					changedNGSPathCount++
				}
			}
		}
	}

	// check code
	for i := range edgesArr {
		e := &edgesArr[i]
		if i < 2 || e.GetDeleteFlag() > 0 || len(e.NGSPathArr) == 0 {
			continue
		}
		idx := 0
		for _, arr := range e.NGSPathArr {
			var deleteFlag bool
			for _, ef := range arr {
				ne := &edgesArr[ef.ID]
				if ne.GetDeleteFlag() > 0 {
					fmt.Fprintf(os.Stderr, "[SmfyDBG]error!!!eID:%d id:%d delete in Arr:%v\n", e.ID, ef.ID, arr)
					deleteFlag = true
					break
				}
			}
			if !deleteFlag {
				if CheckNGSPath(edgesArr, nodesArr, arr) {
					e.NGSPathArr[idx] = arr
					idx++
				}
			}
		}
		e.NGSPathArr = e.NGSPathArr[:idx]
	}

	if DeleteLowFreqEdgeFlag {
		fmt.Printf("[SmfyDBG]kmer num : %d, edge kmer  average freq: %d\n", kmerNum, allKmerFreq/kmerNum)
	}

	fmt.Printf("[SmfyDBG]changed NGSPath number is : %d\n", changedNGSPathCount)
	fmt.Printf("[SmfyDBG]deleted nodes number is : %d\n", deleteNodeNum)
	fmt.Printf("[SmfyDBG]deleted edges number is : %d\n", deleteEdgeNum)
	fmt.Printf("[SmfyDBG]long tips number is : %d\n", longTipsEdgesNum)
}

func CheckDBGSelfCycle(nodesArr []DBGNode, edgesArr []DBGEdge, kmerlen int) {
	for _, e := range edgesArr {
		if e.ID < 2 || e.StartNID < 2 || e.EndNID < 2 || e.GetDeleteFlag() > 0 {
			continue
		}
		if e.StartNID == e.EndNID {
			//var nb constructcf.KmerBnt
			extNb := constructcf.ExtendKmerBnt2Byte(nodesArr[e.StartNID].Seq, kmerlen-1)
			if EqualByteArr(extNb, e.Utg.Ks[:kmerlen-1]) {

			} else {
				if CountEdgeIDComing(nodesArr[e.StartNID].EdgeIDIncoming, e.ID) == 2 {
					//fmt.Printf("[CheckDBGSelfCycle] self cycle edge two end in EdgeIDIncoming,len(e): %v, v: %v\n", len(e.Utg.Ks), nodesArr[e.StartNID])
					if EqualByteArr(extNb, e.Utg.Ks[len(e.Utg.Ks)-(kmerlen-1):]) {
						continue
					}
				}
				e.Utg.Ks = GetReverseCompByteArr(e.Utg.Ks)
				ReverseByteArr(e.Utg.Kq)
				edgesArr[e.ID] = e
				//fmt.Printf("[CheckDBGSelfCycle] ReverSeComp self cycle edge: %v\n", e.ID)
			}
		}
	}
}

func PrintTmpDBG(nodesArr []DBGNode, edgesArr []DBGEdge, prefix string) {
	nodesfn := prefix + ".tmp.nodes"
	edgesfn := prefix + ".tmp.edges"
	nodesfp, err := os.Create(nodesfn)
	if err != nil {
		log.Fatalf("[PrintTmpDBG] failed to create file: %s, err: %v\n", nodesfn, err)
	}
	defer nodesfp.Close()
	edgesfp, err := os.Create(edgesfn)
	if err != nil {
		log.Fatalf("[PrintTmpDBG] failed to create file: %s, err: %v\n", edgesfn, err)
	}
	defer edgesfp.Close()

	for _, v := range nodesArr {
		if v.ID < 2 || v.GetDeleteFlag() > 0 {
			continue
		}
		s := fmt.Sprintf("ID: %v, EdgeIncoming: %v, EdgeOutcoming: %v\n", v.ID, v.EdgeIDIncoming, v.EdgeIDOutcoming)
		nodesfp.WriteString(s)
	}

	for _, e := range edgesArr {
		if e.ID < 2 || e.GetDeleteFlag() > 0 {
			continue
		}
		s := fmt.Sprintf("ID: %v,len: %v, StartNID: %v, EndNID: %v\n", e.ID, len(e.Utg.Ks), e.StartNID, e.EndNID)
		edgesfp.WriteString(s)
	}

}

func Transform2QSeq(utg Unitig) alphabet.QLetters {
	if len(utg.Ks) != len(utg.Kq) {
		log.Fatalf("[transform2QSeq] len(ks):%d != len(kq):%d\n", len(utg.Ks), len(utg.Kq))
	}
	qs := make(alphabet.QLetters, len(utg.Ks))
	for i := 0; i < len(utg.Ks); i++ {
		var ql alphabet.QLetter
		ql.L = alphabet.Letter(bnt.BitNtCharUp[utg.Ks[i]])
		ql.Q = alphabet.Qphred(utg.Kq[i])
		qs[i] = ql
	}

	return qs
}

func Transform2Letters(ks []byte) alphabet.Letters {
	ls := make(alphabet.Letters, len(ks))
	for i, b := range ks {
		//ls = append(ls, alphabet.Letter(bnt.BitNtCharUp[b]))
		ls[i] = alphabet.Letter(bnt.BitNtCharUp[b])
	}

	return ls
}

func Transform2Char(ks []byte) []byte {
	cs := make([]byte, len(ks))
	for i, b := range ks {
		//fmt.Printf("[Transform2Char] ks[%v]: %c\n", i, b)
		//ls = append(ls, alphabet.Letter(bnt.BitNtCharUp[b]))
		cs[i] = bnt.BitNtCharUp[b]
	}
	return cs
}

func Transform2Char2(ks []byte) []byte {
	//cs := make([]byte, len(ks))
	for i, b := range ks {
		//fmt.Printf("[Transform2Char] ks[%v]: %c\n", i, b)
		//ls = append(ls, alphabet.Letter(bnt.BitNtCharUp[b]))
		/*if b > 3 {
			b = 0
		}*/
		ks[i] = bnt.BitNtCharUp[b]
	}
	return ks
}

func Transform2Qual(kq []byte) []byte {
	//cs := make([]byte, len(ks))
	for i, q := range kq {
		//fmt.Printf("[Transform2Char] ks[%v]: %c\n", i, b)
		//ls = append(ls, alphabet.Letter(bnt.BitNtCharUp[b]))
		/*if b > 3 {
			b = 0
		}*/
		kq[i] = q + 33
	}
	return kq
}

func Transform2BntByte(ks []byte) []byte {
	bs := make([]byte, len(ks))
	for i, c := range ks {
		bs[i] = bnt.Base2Bnt[c]
	}
	return bs
}

func Transform2Unitig(Seq alphabet.QLetters, qual bool) (utg Unitig) {

	utg.Ks = make([]byte, len(Seq))
	if qual {
		utg.Kq = make([]uint8, len(Seq))
	}
	for i, v := range Seq {
		utg.Ks[i] = bnt.Base2Bnt[v.L]
		if qual {
			utg.Kq[i] = uint8(v.Q)
		}
	}

	return utg
}

var AdpaterSeq = "cggccgcaaggggttcgcgtcagcgggtgttggcgggtgtcggggctggcttaactatgcggcatcagagcagattgtactgagagtgcaccatatgcggtgtgaaataccacacagatgcgtaaggagaaaataccgcatcaggcgccattcgccattcagctgcgcaactgttgggaagggcgatcggtgcgggcctc"

// Set default quality(default = 1)
func SetDefaultQual(seq Unitig) (new Unitig) {
	for i := 0; i < len(seq.Ks); i++ {
		seq.Ks[i] = bnt.Base2Bnt[seq.Ks[i]]
		seq.Kq = append(seq.Kq, uint8(1))
	}

	return seq
}

func GetEdgesByteArr(edgesArr []DBGEdge, wc chan<- []byte) {
	for _, e := range edgesArr {
		if e.ID > 1 && e.GetDeleteFlag() == 0 && e.GetSeqLen() > 0 {
			eb := make([]byte, 0, e.GetSeqLen()*2+30+len(e.NGSPathArr)*30)
			eb = append(eb, '@')
			eb = append(eb, strconv.Itoa(int(e.ID))...)
			eb = append(eb, '\t')
			eb = append(eb, strconv.Itoa(int(e.StartNID))...)
			eb = append(eb, '\t')
			eb = append(eb, strconv.Itoa(int(e.EndNID))...)
			eb = append(eb, "\tpath:"...)
			if len(e.NGSPathArr) > 0 {
				for _, ep := range e.NGSPathArr {
					if len(ep) < 3 {
						continue
					}
					for _, ef := range ep {
						eb = append(eb, strconv.Itoa(int(ef.ID))...)
						eb = append(eb, '-')
					}
					eb = append(eb, ',')
				}
			}

			eb = append(eb, "\tFreq:0\tlen:"...)
			eb = append(eb, strconv.Itoa(e.GetSeqLen())...)
			eb = append(eb, '\n')
			eb = append(eb, Transform2Char2(e.Utg.Ks)...)
			eb = append(eb, "\n+\n"...)
			eb = append(eb, Transform2Qual(e.Utg.Kq)...)
			eb = append(eb, '\n')
			wc <- eb
		}
	}

	close(wc)
}

func StoreEdgesToFn(edgesbrfn string, edgesArr []DBGEdge) {
	wc := make(chan []byte, 2000)
	go GetEdgesByteArr(edgesArr, wc)

	fp, err := os.Create(edgesbrfn)
	if err != nil {
		log.Fatalf("[StoreEdgesToFn] create file: %s failed, err: %v\n", edgesbrfn, err)
	}
	defer fp.Close()
	cbrofp := cbrotli.NewWriter(fp, cbrotli.WriterOptions{Quality: 1})
	defer cbrofp.Close()
	buffp := bufio.NewWriterSize(cbrofp, 1<<20)
	count := 0
	for {
		eb, ok := <-wc
		if !ok {
			break
		}
		buffp.Write(eb)
		count++
	}

	if err := buffp.Flush(); err != nil {
		log.Fatalf("[StoreEdgesToFn] failed to flush file: %s, err: %v\n", edgesbrfn, err)
	}

	if err := cbrofp.Flush(); err != nil {
		log.Fatalf("[StoreEdgesToFn] failed to flush file: %s, err: %v\n", edgesbrfn, err)
	}
	fmt.Printf("[StoreEdgesToFn] total write %d edges to file:%v\n", count, edgesbrfn)
}

func StoreMappingEdgesToFn(edgesfn string, edgesArr []DBGEdge, MaxMapEdgeLen int) {
	fp, err := os.Create(edgesfn)
	if err != nil {
		log.Fatalf("[StoreMappingEdgesToFn] create file: %s failed, err: %v\n", edgesfn, err)
	}
	defer fp.Close()

	fafp := fasta.NewWriter(fp, 80)
	for _, v := range edgesArr {
		if v.ID > 0 && v.GetDeleteFlag() == 0 {
			if len(v.Utg.Ks) > MaxMapEdgeLen {
				seq1 := linear.NewSeq("", nil, alphabet.DNA)
				seq2 := linear.NewSeq("", nil, alphabet.DNA)
				seq1.ID = strconv.Itoa(int(v.ID)) + "/1"
				seq2.ID = strconv.Itoa(int(v.ID)) + "/2"
				//fmt.Printf("[StoreMappingEdgesToFn] v.Utg.Ks[:MaxMapEdgeLen/2]: %v\n", v.Utg.Ks[:MaxMapEdgeLen/2])
				seq1.AppendLetters(Transform2Letters(v.Utg.Ks[:MaxMapEdgeLen/2])...)
				seq2.AppendLetters(Transform2Letters(v.Utg.Ks[len(v.Utg.Ks)-MaxMapEdgeLen/2:])...)
				ans := strconv.Itoa(int(v.StartNID)) + "\t" + strconv.Itoa(int(v.EndNID)) + "\tlen:" + strconv.Itoa(seq1.Len())
				seq1.Annotation.SetDescription(ans)
				seq2.Annotation.SetDescription(ans)
				_, err1 := fafp.Write(seq1)
				_, err2 := fafp.Write(seq2)
				if err1 != nil || err2 != nil {
					log.Fatalf("[StoreMappingEdgesToFn] write seq1: %v\n\tseq2: %v\n\terr1: %v\terr2: %v\n", seq1, seq2, err1, err2)
				}
			} else {
				//fmt.Printf("[StoreMappingEdgesToFn] v.Utg.Ks: %v\n", v.Utg.Ks)
				seq := linear.NewSeq("", nil, alphabet.DNA)
				seq.ID = strconv.Itoa(int(v.ID))
				la := Transform2Letters(v.Utg.Ks)
				seq.AppendLetters(la...)
				ans := strconv.Itoa(int(v.StartNID)) + "\t" + strconv.Itoa(int(v.EndNID)) + "\tlen:" + strconv.Itoa(seq.Len())
				seq.Annotation.SetDescription(ans)
				_, err := fafp.Write(seq)
				if err != nil {
					log.Fatalf("[StoreMappingEdgesToFn] write seq: %v; err: %v\n", seq, err)
				}
			}
		}
	}
}

func LoadEdgesfqFromFn(brfn string, edgesArr []DBGEdge, qual bool) {
	fp, err := os.Open(brfn)
	if err != nil {
		log.Fatalf("[LoadEdgesfaFromFn] open file: %s error: %v\n", brfn, err)
	}
	defer fp.Close()
	brfp := cbrotli.NewReader(fp)
	defer brfp.Close()
	buffp := bufio.NewReaderSize(brfp, 1<<20)
	fqfp := fastq.NewReader(buffp, linear.NewQSeq("", nil, alphabet.DNA, alphabet.Sanger))
	for {
		if s, err := fqfp.Read(); err != nil {
			if err == io.EOF {
				break
			} else {
				log.Fatalf("[LoadEdgesfqFromFn] read file: %s error: %v\n", brfn, err)
			}
		} else {
			l := s.(*linear.QSeq)
			var edge DBGEdge
			id, err := strconv.Atoi(l.ID)
			if err != nil {
				log.Fatalf("[LoadEdgesfqFromFn] parse Name:%s of fastq err: %v\n", l.ID, err)
			}
			edge.ID = DBG_MAX_INT(id)
			var ps string
			var freq int
			var lenKs int
			_, err = fmt.Sscanf(l.Description(), "%d\t%d\t%v\tFreq:%d\tlen:%d\n", &edge.StartNID, &edge.EndNID, &ps, &freq, &lenKs)
			if err != nil {
				log.Fatalf("[LoadEdgesfqFromFn] parse Description:%s of fastq err: %v\n", l.Description(), err)
			}
			if len(ps) > 5 {
				pa := strings.Split(ps[5:], ",")
				for _, p := range pa {
					var path []EdgeFreq
					for _, item := range strings.Split(p, "-") { // ps[:5] == "path:"
						if item == "" {
							continue
						}
						id, err := strconv.Atoi(item)
						if err != nil {
							log.Fatalf("[LoadEdgesFqFromFn] path: %v convert to int err: %v\n", p, err)
						}
						var ef EdgeFreq
						ef.ID = DBG_MAX_INT(id)
						path = append(path, ef)
					}
					//path.Freq = freq
					edge.NGSPathArr = append(edge.NGSPathArr, path)
				}
				//edge.PathMat = append(edge.PathMat, path)
			}
			edge.Utg = Transform2Unitig(l.Seq, qual)
			if edge.ID >= DBG_MAX_INT(len(edgesArr)) {
				log.Fatalf("[LoadEdgesfqFromFn] edge.ID:%v >= len(edgesArr):%d\n", edge.ID, len(edgesArr))
			} else if edgesArr[edge.ID].ID > 0 {
				log.Fatalf("[LoadEdgesfqFromFn] the position: %v in edgesArr has value:%v\n", edge.ID, edgesArr[edge.ID])
			}
			edgesArr[edge.ID] = edge
		}
	}

}

func Set(ea1, ea2 []DBG_MAX_INT) []DBG_MAX_INT {
	arr := make([]DBG_MAX_INT, len(ea1))
	copy(arr, ea1)
	for _, id := range ea2 {
		j := 0
		for ; j < len(arr); j++ {
			if arr[j] == id {
				break
			}
		}
		if j == len(arr) {
			arr = append(arr, id)
		}
	}
	return arr
}

func InterSecDBG_MAX_INTArr(ea1, ea2 []DBG_MAX_INT) (arr []DBG_MAX_INT) {
	for _, id := range ea2 {
		for j := 0; j < len(ea1); j++ {
			if ea1[j] == id {
				arr = append(arr, id)
				break
			}
		}
	}
	return arr
}

func IsEIDArrBubble(ea []DBG_MAX_INT, edgesArr []DBGEdge) (ok bool) {
	ok = true
	e1 := &edgesArr[ea[0]]
	if e1.StartNID == e1.EndNID {
		ok = false
		return
	}
	for _, eID := range ea[1:] {
		e2 := &edgesArr[eID]
		if (e1.StartNID == e2.StartNID && e1.EndNID == e2.EndNID) || (e1.StartNID == e2.EndNID && e1.EndNID == e2.StartNID) {

		} else {
			ok = false
			break
		}
	}
	return ok
}

type PathInfo struct {
	Path       []DBG_MAX_INT
	LastStrand bool
	Direction  uint8
	SeqLen     int
}

func GetNextPathArr(le *DBGEdge, nd DBGNode, MaxPathLen int, MaxBubbleSeqLen int, edgesArr []DBGEdge, nodesArr []DBGNode, kmerlen int) [][]DBG_MAX_INT {
	var pi PathInfo
	pi.Direction = FORWARD
	pi.Path = append(pi.Path, le.ID)
	if le.StartNID == nd.ID {
		pi.LastStrand = false
	} else {
		pi.LastStrand = true
	}
	var pathArr [][]DBG_MAX_INT
	pi.SeqLen = kmerlen - 1
	stk := make([]PathInfo, 0, 20)
	stk = append(stk, pi)

	for len(stk) > 0 {
		if len(stk) > 10 || len(pathArr) > 20 {
			pathArr = nil
			break
		}
		t := stk[len(stk)-1]
		stk = stk[:len(stk)-1]

		e := &edgesArr[t.Path[len(t.Path)-1]]
		eIDArr := GetNextEdgeIDArr(e, nodesArr, t.LastStrand, t.Direction)
		//if len(eIDArr) == 0 {
		//	pathArr = append(pathArr, t.Path[1:])
		//	continue
		//}
		for _, eID := range eIDArr {
			ne := &edgesArr[eID]
			var npi PathInfo
			npi.SeqLen = t.SeqLen + ne.GetSeqLen() - (kmerlen - 1)
			npi.Path = append(npi.Path, t.Path...)
			npi.Path = append(npi.Path, eID)
			npi.Direction = t.Direction
			npi.LastStrand = GetNextEdgeStrand(e, ne, nodesArr, t.LastStrand)
			if len(npi.Path) >= MaxPathLen || npi.SeqLen > MaxBubbleSeqLen {
				pathArr = append(pathArr, npi.Path[1:])
				//fmt.Printf("[GetNextPathArr]eID:%v SeqLen:%d Path:%v\n", e.ID, npi.SeqLen, npi.Path)
			} else {
				stk = append(stk, npi)
			}
		}
	}
	//fmt.Printf("[GetNextPathArr]pathArr:%v\n", pathArr)

	return pathArr
}

type BubblePath struct {
	Idx int
	//EPIPath []constructdbg.DBG_MAX_INT
	Path []DBG_MAX_INT
}

func IsInBubblePathArr(bubblePathArr []BubblePath, bp BubblePath) bool {
	for _, t := range bubblePathArr {
		if reflect.DeepEqual(t.Path, bp.Path) {
			return true
		}
	}
	return false
}

func SortPathArr(pathArr [][]DBG_MAX_INT, edgesArr []DBGEdge, kmerlen int) [][]DBG_MAX_INT {
	spa := make([][]DBG_MAX_INT, 0, len(pathArr))
	idxArr := make([]bool, len(pathArr))
	for i := 0; i < len(pathArr); i++ {
		max := math.MinInt32
		idx := -1
		for j := 0; j < len(pathArr); j++ {
			if idxArr[j] {
				continue
			}
			p := pathArr[j]
			if max < edgesArr[p[0]].GetSeqLen() {
				max = edgesArr[p[0]].GetSeqLen()
				idx = j
			} else if max == edgesArr[p[0]].GetSeqLen() {
				if len(pathArr[idx]) > 1 && len(p) > 1 {
					if edgesArr[pathArr[idx][1]].GetSeqLen() < edgesArr[p[1]].GetSeqLen() {
						idx = j
					}
				}
			}
		}
		if idx >= 0 {
			idxArr[idx] = true
			spa = append(spa, pathArr[idx])
		}
	}

	return spa
}

func IsContainBubPathArr(p []DBG_MAX_INT, bubPathArr []BubblePath) bool {
	ok := false
	for _, bp := range bubPathArr {
		if len(bp.Path) < len(p) && reflect.DeepEqual(bp.Path, p[:len(bp.Path)]) {
			ok = true
			break
		}
	}
	return ok
}

type PBInfo struct {
	Path []DBG_MAX_INT
	Idx  int
	Cum  int
}

type PBInfoArr []PBInfo

func (arr PBInfoArr) Len() int {
	return len(arr)
}

func (arr PBInfoArr) Less(i, j int) bool {
	return arr[i].Cum > arr[j].Cum
}

func (arr PBInfoArr) Swap(i, j int) {
	arr[i], arr[j] = arr[j], arr[i]
}

func GetBubblePathArr(pathArr [][]DBG_MAX_INT, edgesArr []DBGEdge, kmerlen int, MaxPathLen, MaxBubbleSeqLen int) (bubPathArr []BubblePath, num int) {
	if len(pathArr) == 0 {
		return
	}
	minBubEdgeLen := 2 * (kmerlen - 1) * 9 / 10
	pathArr = SortPathArr(pathArr, edgesArr, kmerlen)
	cumLen := make([]int, len(pathArr))
	idxArr := make([]int, len(pathArr))
	for j := 0; j < MaxPathLen; j++ {
		for i, p := range pathArr {
			if j >= len(p) {
				continue
			}
			if cumLen[i] >= minBubEdgeLen {
				continue
			}
			if j == 0 {
				cumLen[i] += edgesArr[p[j]].GetSeqLen()
			} else {
				cumLen[i] += edgesArr[p[j]].GetSeqLen() - (kmerlen - 1)
			}
			idxArr[i] = j
		}
	}

	var finished bool
	for x, p := range pathArr {
		if cumLen[x] < minBubEdgeLen || idxArr[x] >= len(p)-1 {
			continue
		}
		/*if cumLen[x] < 2*(kmerlen-1) {
			for w := idxArr[x]; w < len(p)-1; w++ {
				if cumLen[x] < 2*(kmerlen-1) {
					cumLen[x] += edgesArr[p[w]].GetSeqLen() - (kmerlen - 1)
					idxArr[x] = w
				} else {
					break
				}
			}
		} */
		for y := x + 1; y < len(pathArr); y++ {
			np := pathArr[y]
			if cumLen[y] < minBubEdgeLen || idxArr[y] >= len(np)-1 {
				continue
			}
			if p[0] == np[0] {
				continue
			}
			cy := cumLen[y]
			z := idxArr[y]
			for ; z < len(np)-1; z++ {
				if np[z+1] == p[idxArr[x]+1] {
					break
				}
				if cy < cumLen[x]*11/10 {
					cy += edgesArr[np[z+1]].GetSeqLen() - (kmerlen - 1)
				} else {
					break
				}
			}
			if z >= len(np)-1 || np[z+1] != p[idxArr[x]+1] {
				continue
			}
			if utils.AbsInt(cumLen[x]-cy) < Min(cumLen[x], cy)/10 {
				var bpx, bpy BubblePath
				bpx.Path = append(bpx.Path, pathArr[x][:idxArr[x]+1]...)
				bpy.Path = append(bpy.Path, pathArr[y][:z+1]...)
				bubPathArr = append(bubPathArr, bpx, bpy)
				for w := y + 1; w < len(pathArr); w++ {
					zp := pathArr[w]
					if cumLen[w] < minBubEdgeLen || idxArr[w] >= len(zp)-1 {
						continue
					}
					if p[0] == zp[0] {
						continue
					}
					cw := cumLen[w]
					v := idxArr[w]
					for ; v < len(zp)-1; v++ {
						if zp[v+1] == p[idxArr[x]+1] {
							break
						}
						if cw < cumLen[x]*11/10 {
							cw += edgesArr[zp[v+1]].GetSeqLen() - (kmerlen - 1)
						} else {
							break
						}
					}
					if v >= len(zp)-1 || zp[v+1] != p[idxArr[x]+1] {
						continue
					}
					if utils.AbsInt(cumLen[x]-cw) < Min(cumLen[x], cw)/10 {
						var bp BubblePath
						bp.Path = append(bp.Path, pathArr[w][:v+1]...)
						if !IsInBubblePathArr(bubPathArr, bp) {
							bubPathArr = append(bubPathArr, bp)
						}
					}
				}
				finished = true
				break
			}
		}
		if finished {
			break
		}
	}
	if len(bubPathArr) > 0 {
		for _, p := range pathArr {
			if IsContainBubPathArr(p, bubPathArr) {
				num++
			}
		}
		return
	}

	for j := idxArr[0]; j < len(pathArr[0])-2; j++ {
		if len(bubPathArr) == 0 && len(pathArr[0]) > idxArr[0]+2 {
			finished = false
			nparr := make([]PBInfo, len(pathArr))
			for i, p := range pathArr {
				if idxArr[i] >= len(p)-2 {
					continue
				}
				idxArr[i]++
				cumLen[i] += edgesArr[idxArr[i]].GetSeqLen() - (kmerlen - 1)
				nparr[i].Path = p
				nparr[i].Idx = idxArr[i]
				nparr[i].Cum = cumLen[i]
			}
			sort.Sort(PBInfoArr(nparr))
			for i, pb := range nparr {
				pathArr[i] = pb.Path
				idxArr[i] = pb.Idx
				cumLen[i] = pb.Cum
			}
			for x, p := range pathArr {
				if cumLen[x] < minBubEdgeLen || idxArr[x] >= len(p)-1 {
					continue
				}
				for y := x + 1; y < len(pathArr); y++ {
					np := pathArr[y]
					if cumLen[y] < minBubEdgeLen || idxArr[y] >= len(np)-1 {
						continue
					}
					if p[0] == np[0] {
						continue
					}
					cy := cumLen[y]
					z := idxArr[y]
					for ; z < len(np)-1; z++ {
						if np[z+1] == p[idxArr[x]+1] {
							break
						}
						if cy < cumLen[x]*11/10 {
							cy += edgesArr[np[z+1]].GetSeqLen() - (kmerlen - 1)
						} else {
							break
						}
					}
					if z >= len(np)-1 || np[z+1] != p[idxArr[x]+1] {
						continue
					}
					if utils.AbsInt(cumLen[x]-cy) < Min(cumLen[x], cy)/10 {
						var bpx, bpy BubblePath
						bpx.Path = append(bpx.Path, pathArr[x][:idxArr[x]+1]...)
						bpy.Path = append(bpy.Path, pathArr[y][:z+1]...)
						bubPathArr = append(bubPathArr, bpx, bpy)
						for w := y + 1; w < len(pathArr); w++ {
							zp := pathArr[w]
							if cumLen[w] < minBubEdgeLen || idxArr[w] >= len(zp)-1 {
								continue
							}
							if p[0] == zp[0] {
								continue
							}
							cw := cumLen[w]
							v := idxArr[w]
							for ; v < len(zp)-1; v++ {
								if zp[v+1] == p[idxArr[x]+1] {
									break
								}
								if cw < cumLen[x]*11/10 {
									cw += edgesArr[zp[v+1]].GetSeqLen() - (kmerlen - 1)
								} else {
									break
								}
							}
							if v >= len(zp)-1 || zp[v+1] != p[idxArr[x]+1] {
								continue
							}
							if utils.AbsInt(cumLen[x]-cw) < Min(cumLen[x], cw)/10 {
								var bp BubblePath
								bp.Path = append(bp.Path, pathArr[w][:v+1]...)
								if !IsInBubblePathArr(bubPathArr, bp) {
									bubPathArr = append(bubPathArr, bp)
								}
							}
						}
						finished = true
						break
					}
				}
				if finished {
					break
				}
			}
			if finished {
				break
			}
		} else {
			break
		}
	}

	if len(bubPathArr) > 0 {
		for _, p := range pathArr {
			if IsContainBubPathArr(p, bubPathArr) {
				num++
			}
		}
		return
	}

	return
}

// set the unique edge of edgesArr
func SetDBGEdgesUniqueFlag(edgesArr []DBGEdge, nodesArr []DBGNode, kmerlen int) (uniqueNum, bubbleRepeatNum, twoCycleNum, selfCycle, selfCycleSameComingNum, bubbleEdgeNum int) {
	MaxPathLen := 5
	MaxBubbleSeqLen := 3 * kmerlen
	TipRegionNum := 6
	tipNumArr := make([]int, TipRegionNum)
	tipLenSumArr := make([]int, TipRegionNum)
	LenBoundArr := [5]int{500, 1000, 2000, 5000, 10000}
	var bubbleNum, bubbleTotalLen, SD, longTipNum, tipTotalLen, selfCycleTotalLen, uniqueTotalLen, repeatNum, repeatTotalLen, bubbleRepeatLen int
	for i := range edgesArr {
		e := &edgesArr[i]
		if e.ID < 2 || e.GetDeleteFlag() > 0 {
			continue
		}
		//if e.GetTwoEdgesCycleFlag() > 0 {
		//	continue
		//}
		el := e.GetSeqLen()

		if e.StartNID == 0 || e.EndNID == 0 {
			longTipNum++
			tipTotalLen += el
			for j, lb := range LenBoundArr {
				if el < lb {
					tipNumArr[j]++
					tipLenSumArr[j] += el
					break
				} else if j == len(LenBoundArr)-1 {
					tipNumArr[j+1]++
					tipLenSumArr[j+1] += el
				}
			}
			//edgesArr[i].SetUniqueFlag()
			fmt.Fprintf(os.Stderr, "[SetDBGEdgesUniqueFlag]eID:%d len:%d tip\n", e.ID, el)
			//continue
		} else if e.StartNID == e.EndNID {
			f1 := IsInComing(nodesArr[e.StartNID].EdgeIDIncoming, e.ID)
			f2 := IsInComing(nodesArr[e.StartNID].EdgeIDOutcoming, e.ID)
			if !(f1 && f2) && (f1 || f2) {
				selfCycleSameComingNum++
				fmt.Fprintf(os.Stderr, "[SetDBGEdgesUniqueFlag]eID:%d len:%d selfCycleSameComing\n", e.ID, el)
			} else {
				fmt.Fprintf(os.Stderr, "[SetDBGEdgesUniqueFlag]eID:%d len:%d selfCycle\n", e.ID, el)
			}
			selfCycle++
			selfCycleTotalLen += el
			continue
		}

		var ea1, ea2 []DBG_MAX_INT
		ea1 = GetNextEdgeIDArr(e, nodesArr, PLUS, BACKWARD)
		ea2 = GetNextEdgeIDArr(e, nodesArr, PLUS, FORWARD)
		if len(ea1)+len(ea2) <= 1 {
			e.SetUniqueFlag()
			continue
		} else if len(ea1) == 1 && len(ea2) == 1 && ea1[0] == ea2[0] {
			e.SetTwoEdgesCycleFlag()
			twoCycleNum++
			fmt.Fprintf(os.Stderr, "[SetDBGEdgesUniqueFlag]eID:%d len:%d TwoEdgesCycle\n", e.ID, el)
			continue
		}

		if len(ea1) == 1 && len(ea2) == 1 {
			if ea, ok := IsBubbleEdge(e, nodesArr, edgesArr); ok {
				edgesArr[i].SetBubbleFlag()
				fmt.Fprintf(os.Stderr, "[SetDBGEdgesUniqueFlag]eID:%d len:%d Bubble\n", e.ID, el)
				bubbleEdgeNum++
				bubbleTotalLen += el
				//ea := GetOtherEArr(nodesArr[e.StartNID], e.ID)
				for _, id := range ea {
					if IsBubble(e.ID, id, edgesArr) {
						//fmt.Printf("[SetDBGEdgesUniqueFlag] Bubble edge: %v, length: %v\n", e.ID, len(e.Utg.Ks))
						if el > edgesArr[id].GetSeqLen() {
							bubbleNum++
							SD += el - edgesArr[id].GetSeqLen()
						}
					}
				}
			} else {
				e.SetUniqueFlag()
				uniqueNum++
				uniqueTotalLen += el
			}
			continue
		}

		if len(ea1) >= 1 && len(ea2) >= 1 && len(ea1)+len(ea2) > 2 {
			// test if bubble repeat
			oea1 := GetOtherEArr(nodesArr[e.StartNID], e.ID)
			oea2 := GetOtherEArr(nodesArr[e.EndNID], e.ID)
			ok := true
			if len(ea1) > 1 {
				if len(oea1) > 0 {
					ok = false
				} else {
					if !IsEIDArrBubble(ea1, edgesArr) {
						nd := nodesArr[e.StartNID]
						pathArr := GetNextPathArr(e, nd, MaxPathLen+1, MaxBubbleSeqLen, edgesArr, nodesArr, kmerlen)
						_, num := GetBubblePathArr(pathArr, edgesArr, kmerlen, MaxPathLen, MaxBubbleSeqLen)
						if num > 0 && num == len(pathArr) {

						} else {
							ok = false
						}
					}
				}
			}

			if len(ea2) > 1 {
				if len(oea2) > 0 {
					ok = false
				} else {
					if !IsEIDArrBubble(ea2, edgesArr) {
						nd := nodesArr[e.EndNID]
						pathArr := GetNextPathArr(e, nd, MaxPathLen+1, MaxBubbleSeqLen, edgesArr, nodesArr, kmerlen)
						_, num := GetBubblePathArr(pathArr, edgesArr, kmerlen, MaxPathLen, MaxBubbleSeqLen)
						if num > 0 && num == len(pathArr) {

						} else {
							ok = false
						}
					}
				}
			}
			if ok {
				e.SetBubbleRepeatFlag()
				fmt.Fprintf(os.Stderr, "[SetDBGEdgesUniqueFlag]eID:%d len:%d BubbleRepeat\n", e.ID, el)
				bubbleRepeatNum++
				bubbleRepeatLen += el
				continue
			}
		}
		repeatNum++
		repeatTotalLen += el
	}
	//bubbleNum /= 2
	if bubbleNum > 0 {
		fmt.Printf("[SetDBGEdgesUniqueFlag] Bubble Number: %v, average length:%v, SD: %v\n", bubbleEdgeNum, bubbleTotalLen/bubbleEdgeNum, SD/bubbleNum)
	}
	if uniqueNum > 0 {
		fmt.Printf("[SetDBGEdgesUniqueFlag] unique  Number: %v, average length:%v, \n", uniqueNum, uniqueTotalLen/uniqueNum)
	}
	if bubbleRepeatNum > 0 {
		fmt.Printf("[SetDBGEdgesUniqueFlag] bubble repeat  Number: %v, average length:%v, \n", bubbleRepeatNum, bubbleRepeatLen/bubbleRepeatNum)
	}
	if repeatNum > 0 {
		fmt.Printf("[SetDBGEdgesUniqueFlag] repeat  Number: %v, average length:%v, \n", repeatNum, repeatTotalLen/repeatNum)
	}
	if longTipNum > 0 {
		fmt.Printf("[SetDBGEdgesUniqueFlag] long tips Number: %v, average length: %v\n", longTipNum, tipTotalLen/longTipNum)
	}
	if selfCycle > 0 {
		fmt.Printf("[SetDBGEdgesUniqueFlag] self cycle edge average length: %v\n", selfCycleTotalLen/selfCycle)
	}
	// output tip info
	for i, lb := range LenBoundArr {
		fmt.Printf("[SetDBGEdgesUniqueFlag] Tip Length <%d: num is %d, length summary is %d\n", lb, tipNumArr[i], tipLenSumArr[i])
	}
	fmt.Printf("[SetDBGEdgesUniqueFlag] Tip Length >=%d: num is %d, length summary is %d\n", LenBoundArr[len(LenBoundArr)-1], tipNumArr[len(LenBoundArr)], tipLenSumArr[len(LenBoundArr)])
	return
}

func AddNodeInfo2DBGEdgeArr(edgesArr []DBGEdge, nodesArr []DBGNode) {
	for i := range edgesArr {
		e := &edgesArr[i]
		if e.ID < 2 || e.GetDeleteFlag() > 0 {
			continue
		}
		if e.EndNID > 0 {
			nd := nodesArr[e.EndNID]
			f1 := IsInComing(nd.EdgeIDIncoming, e.ID)
			f2 := IsInComing(nd.EdgeIDOutcoming, e.ID)
			if f1 || f2 {
				if (f1 && f2) || f1 {
					edgesArr[i].EdgeIDOutcoming = nd.EdgeIDOutcoming
				} else {
					edgesArr[i].EdgeIDOutcoming = nd.EdgeIDIncoming
				}
			}
		}

		if e.StartNID > 0 {
			nd := nodesArr[e.StartNID]
			f1 := IsInComing(nd.EdgeIDIncoming, e.ID)
			f2 := IsInComing(nd.EdgeIDOutcoming, e.ID)
			if f1 || f2 {
				if (f1 && f2) || f2 {
					edgesArr[i].EdgeIDIncoming = nd.EdgeIDIncoming
				} else {
					edgesArr[i].EdgeIDIncoming = nd.EdgeIDOutcoming
				}
			}
		}

		// check NGSPath
		if len(e.NGSPathArr) == 0 {
			continue
		}
		idx := 0
		for _, arr := range e.NGSPathArr {
			if CheckNGSPath(edgesArr, nodesArr, arr) {
				e.NGSPathArr[idx] = arr
				idx++
			}
		}
		if idx < len(e.NGSPathArr) {
			fmt.Printf("[AddNodeInfo2DBGEdgeArr]eID:%d len(NGSPathArr):%d delete NGSPathNum:%d\n", e.ID, len(e.NGSPathArr), len(e.NGSPathArr)-idx)
		}
		e.NGSPathArr = e.NGSPathArr[:idx]
	}

}

/*type EdgeMapInfo struct {
	ID DBG_MAX_INT
	Strand bool
}*/
type PathSeq struct {
	ID DBG_MAX_INT
	//NID    DBG_MAX_INT // the path end node ID
	Strand bool
	//Start, End int
}

type ReadMapInfo struct {
	ID           int64
	StartP, EndP int
	//NID          DBG_MAX_INT
	//Seq        []byte
	PathSeqArr []PathSeq
	//Strands    []bool
}

type AlignInfo struct {
	ID      int64
	EndPos  int
	Paths   []DBG_MAX_INT
	Strands []bool // strand for Paths
	Seq     []byte
}

func paraLoadNGSReads(brfn string, cs chan constructcf.ReadInfo, kmerLen int, we chan int) {
	var count int
	format := constructcf.GetReadsFileFormat(brfn)
	bufSize := (1 << 26)
	size := 1000
	fncs1 := make(chan []byte, 3)
	recordChan1 := make(chan constructcf.ReadInfo, size)

	go constructcf.ReadBrFile(brfn, fncs1, bufSize)
	go constructcf.GetReadFileRecord(fncs1, recordChan1, format, bufSize, true)

	for {
		ri, ok := <-recordChan1
		if !ok {
			break
		}
		cs <- ri
		count++
	}
	we <- count
}

func LoadNGSReads(cfgFn string, correct bool, cs chan constructcf.ReadInfo, numCPU, kmerLen int) {
	cfgInfo, err := constructcf.ParseCfg(cfgFn, false, correct)
	if err != nil {
		log.Fatal("[constructcf.ParseCfg] found err")
	}
	fmt.Println(cfgInfo)

	var totalNumReads int

	// iterate cfgInfo find fastq files
	we := make(chan int, numCPU)
	for i := 0; i < numCPU; i++ {
		we <- 0
	}
	//var numT int
	for _, lib := range cfgInfo.Libs {
		// seqProfile == 1 note Illumina
		if lib.AsmFlag != constructcf.AllState && lib.SeqProfile != 1 {
			continue
		}
		for _, fn := range lib.FnName {
			totalNumReads += <-we
			go paraLoadNGSReads(fn, cs, kmerLen, we)
		}
	}

	// check child routinue have finishedT
	for i := 0; i < numCPU; i++ {
		totalNumReads += <-we
	}
	fmt.Printf("[LoadNGSReads] total processed number reads : %v\n", totalNumReads)

	// send map goroutinues finish signal
	close(cs)
	//time.Sleep(time.Second * 10)
}

func writeAlignToFile(wrFn string, wc chan AlignInfo, numCPU int) {
	fp, err := os.Create(wrFn)
	if err != nil {
		log.Fatalf("[writeAlignToFile] failed to create file: %s, err:%v\n", wrFn, err)
	}
	defer fp.Close()
	buffp := bufio.NewWriterSize(fp, 1<<25)

	var finishedT int
	for {
		ai := <-wc
		if len(ai.Paths) == 0 {
			finishedT++
			if finishedT == numCPU {
				break
			} else {
				continue
			}
		}
		ps := make([]string, len(ai.Paths))
		for i := 0; i < len(ai.Paths); i++ {
			ps[i] = strconv.Itoa(int(ai.Paths[i]))
		}
		s := fmt.Sprintf("%d\t%s\n", ai.ID, strings.Join(ps, ":"))
		buffp.WriteString(s)
	}

	if err = buffp.Flush(); err != nil {
		log.Fatalf("[writeAlignToFile] failed to Flush buffer file: %s, err:%v\n", wrFn, err)
	}
}

func BiggerThan(kb, rb []byte) bool {
	if len(kb) != len(rb) {
		log.Fatalf("[BiggerThan] len(kb): %v != len(rb): %v\n", len(kb), len(rb))
	}
	for i, b := range kb {
		if b > rb[i] {
			return true
		} else if b < rb[i] {
			return false
		}
	}
	return false
}

func EqualByteArr(kb, rb []byte) bool {
	if len(kb) != len(rb) {
		log.Fatalf("[BiggerThan] len(kb): %v != len(rb): %v\n", len(kb), len(rb))
	}
	for i, b := range kb {
		if b != rb[i] {
			return false
		}
	}
	return true
}

func GetCuckoofilterDBGSampleSize(edgesArr []DBGEdge, winSize, kmerLen int64) (cfSize int64) {

	for _, e := range edgesArr {
		if e.GetDeleteFlag() > 0 || e.ID < 2 {
			continue
		}
		//fmt.Printf("[GetCuckoofilterDBGSampleSize] e: %v\n", e)
		el := int64(len(e.Utg.Ks))
		cfSize += ((el - kmerLen + winSize - 1) / winSize) + 1
		/*if el < 2*maxNGSReadLen-kmerLen { // get whole length sliding windowns
			cfSize += ((el - kmerLen + winSize - 1) / winSize) + 1
		} else { // just sample two ends
			cfSize += ((maxNGSReadLen-kmerLen+winSize-1)/winSize + 1) * 2
		}*/
	}
	return cfSize
}

func GetMinimizer(seq []byte, kmerlen int) (minSeq []byte, pos int, strand bool) {
	for i := 0; i < len(seq)-kmerlen+1; i++ {
		kb := seq[i : i+kmerlen]
		rb := GetReverseCompByteArr(kb)
		st := PLUS
		if BiggerThan(kb, rb) {
			kb, rb = rb, kb
			st = MINUS
		}
		if len(minSeq) == 0 || BiggerThan(minSeq, kb) {
			minSeq = kb
			pos = i
			strand = st
		}
	}

	return
}

type KmerInfo struct {
	Kmer   []byte
	Pos    int
	Strand bool
}

func FindMinKmerInfo(ka []KmerInfo) (min KmerInfo, idx int) {
	min = ka[0]
	idx = 0
	for i := 1; i < len(ka); i++ {
		ki := ka[i]
		if BiggerThan(min.Kmer, ki.Kmer) {
			min = ki
			idx = i
		}
	}
	return
}

func ConstructCFDBGMinimizers(cf CuckooFilter, edgesArr []DBGEdge, winSize int) (count int) {
	bufSize := 100
	buf := make([]KmerInfo, bufSize)
	for _, e := range edgesArr {
		if e.ID < 2 || e.GetDeleteFlag() > 0 {
			continue
		}
		ks := e.Utg.Ks
		sl := len(ks)
		rs := GetReverseCompByteArr(ks)
		idx := 0
		last := -1
		var ki, rki KmerInfo
		for j := 0; j < winSize-1 && j < sl-(cf.Kmerlen-1); j++ {
			ki.Kmer, ki.Pos, ki.Strand = ks[j:j+cf.Kmerlen], j, PLUS
			rki.Kmer, rki.Pos, rki.Strand = rs[sl-j-cf.Kmerlen:sl-j], j, MINUS
			if BiggerThan(ki.Kmer, rki.Kmer) {
				buf[idx] = rki
			} else {
				buf[idx] = ki
			}
			idx++
		}

		for j := winSize - 1; j < sl-(cf.Kmerlen-1); j++ {
			ki.Kmer, ki.Pos, ki.Strand = ks[j:j+cf.Kmerlen], j, PLUS
			rki.Kmer, rki.Pos, rki.Strand = rs[sl-j-cf.Kmerlen:sl-j], j, MINUS
			if BiggerThan(ki.Kmer, rki.Kmer) {
				buf[idx] = rki
			} else {
				buf[idx] = ki
			}
			idx++

			minKI, x := FindMinKmerInfo(buf[idx-winSize : idx])
			y := idx - winSize + x
			if last < 0 {
				if cf.Insert(minKI.Kmer, e.ID, uint32(minKI.Pos), minKI.Strand) == false {
					log.Fatalf("[constructCFDBGSample] Insert to the CuckooFilter of DBGSample false\n")
				}
				count++
				last = y
			} else {
				if y > last {
					if cf.Insert(minKI.Kmer, e.ID, uint32(minKI.Pos), minKI.Strand) == false {
						log.Fatalf("[constructCFDBGSample] Insert to the CuckooFilter of DBGSample false\n")
					}
					count++
					last = y
				}
			}

			// adjust buf
			if idx == bufSize {
				y := 0
				for x := idx - winSize; x < idx; x++ {
					buf[y] = buf[x]
					y++
				}
				last = last - (idx - winSize)
				idx = y
			}
		}
	}
	return count
}

// found kmer seed position in the DBG edges
func LocateSeedKmerCF(cf CuckooFilter, ri constructcf.ReadInfo, winSize int, edgesArr []DBGEdge, MaxStepNum int, ka []KmerInfo, rSeq []byte) (dbgK DBGKmer, pos int, strand bool) {
	//MaxStepNum := 10
	max := MaxStepNum * winSize
	if len(ri.Seq) < cf.Kmerlen+MaxStepNum*winSize {
		if len(ri.Seq) < cf.Kmerlen+winSize {
			return
		}
		max = (len(ri.Seq) - cf.Kmerlen)
		//fmt.Printf("[LocateSeedKmerCF] read: %v sequence length smaller than KmerLen(%v) + MaxStepNum[%v] * winSize(%v)\n", ri, cf.Kmerlen, MaxStepNum, winSize)
		//MaxStepNum = (len(ri.Seq) - cf.Kmerlen) / winSize
		//return
	}
	rl := len(ri.Seq)
	rSeq = GetReverseCompByteArr2(ri.Seq, rSeq)
	//ka := make([]KmerInfo, MaxStepNum*winSize)
	//fmt.Printf("[constructCFDBGSample] e: %v\n", e)
	kaIdx := 0
	last := -1
	for j := 0; j < max; j++ {
		var ki, rki KmerInfo
		ki.Kmer, ki.Pos, ki.Strand = ri.Seq[j:j+cf.Kmerlen], j, PLUS
		rki.Kmer, rki.Pos, rki.Strand = rSeq[rl-j-cf.Kmerlen:rl-j], j, MINUS
		//fmt.Printf("[LocateSeedKmerCF] j: %v, ki: %v\n", j, ki)
		//fmt.Printf("[LocateSeedKmerCF] j: %v, rki: %v\n", j, rki)
		if BiggerThan(ki.Kmer, rki.Kmer) {
			ka[kaIdx] = rki
		} else {
			ka[kaIdx] = ki
		}
		kaIdx++
		if kaIdx < winSize {
			continue
		}
		min, x := FindMinKmerInfo(ka[kaIdx-winSize : kaIdx])
		newOne := false
		if last < 0 {
			newOne = true
			last = kaIdx - winSize + x
		} else {
			if kaIdx-winSize+x > last {
				newOne = true
				last = kaIdx - winSize + x
			}
		}
		if newOne {
			pos = min.Pos
			strand = min.Strand
			//var count int
			//ki  := GetMinimizer(e.Utg.Ks[j:z], cf.Kmerlen)
			//fmt.Printf("e.ID: %v\tmaxLen: %v\tj: %v,pos: %v, strand: %v, kb: %v\n", e.ID, maxLen, j, pos, strand, kb)
			//kb := make([]byte, len(min.Kmer))
			//copy(kb, min.Kmer)
			dk := cf.Contain(min.Kmer)
			//fmt.Printf("[LocateSeedKmerCF] found seed kmers count: %v, dk: %v\n", dk.GetCount(), dk)
			if dk.GetCount() > 0 {
				dbgK = dk
				return
			}
		}
	}

	// search seed by read end partition
	/*if count > 1 {
		count = 0
		for i := 0; i < 2; i++ {
			startP := len(ri.Seq) - cf.Kmerlen - winSize + 1 - i*winSize
			kb, pos, strand = GetMinimizer(ri.Seq[startP:len(ri.Seq)-i*winSize], cf.Kmerlen)
			dbgK, count = cf.Lookup(kb, edgesArr)
			pos += uint32(startP)
			if count == 1 && dbgK.GetCount() > 0 {
				fmt.Printf("[LocateSeedKmerCF] found seed count: %v\n", count)
				return
			}
			if count > 1 {
				dbgK.setCFItem(0, 0)
				break
			}
		}
	} */

	/*// found the second window
	kb, pos, strand = GetMinimizer(ri.Seq[winSize:winSize+cf.Kmerlen+winSize-1], cf.Kmerlen)
	pos += uint32(winSize)
	dbgK = cf.Lookup(kb, edgesArr)*/

	return
}

const (
	IN  = true
	OUT = false
)

func IsInDBGNode(nd DBGNode, eID DBG_MAX_INT) bool {
	for i := 0; i < bnt.BaseTypeNum; i++ {
		if nd.EdgeIDIncoming[i] == eID || nd.EdgeIDOutcoming[i] == eID {
			return true
		}
	}

	return false
}

func GetNextMappingEID(nd DBGNode, e DBGEdge, base byte) DBG_MAX_INT {
	if e.StartNID == nd.ID {
		if IsInComing(nd.EdgeIDOutcoming, e.ID) {
			return nd.EdgeIDIncoming[base]
		} else {
			return nd.EdgeIDOutcoming[bnt.BntRev[base]]
		}
	} else {
		if IsInComing(nd.EdgeIDIncoming, e.ID) {
			return nd.EdgeIDOutcoming[base]
		} else {
			return nd.EdgeIDIncoming[bnt.BntRev[base]]
		}
	}
}

/*func GetNextMappingEdgeInfo(e DBGEdge, strand bool, direction uint8, ne DBGEdge, node DBGNode, kmerlen int) (bool, int) {
	var edgeNodeSeq [kmerlen-1]byte
	if direction == BACKWARD {
		if strand == PLUS {
			copy(edgeNodeSeq, e.Utg.Ks[:kmerlen-1])
		} else {

		}
	} else { // FORWARD

	}
}*/

func GetComingOtherEArr(coming [bnt.BaseTypeNum]DBG_MAX_INT, eID DBG_MAX_INT) (ea []DBG_MAX_INT) {
	for _, id := range coming {
		if id > 1 && id != eID {
			ea = append(ea, id)
		}
	}
	return
}

/*func GetSelfCycleNextMapEdgeInfo(eID DBG_MAX_INT, nd DBGNode, edgesArr []DBGEdge, nodeSeq []byte, kmerlen int, direction uint8, base byte, correct bool) (ID DBG_MAX_INT, pos int, strand bool) {
	var tmp constructcf.KmerBnt
	//kb.Seq, kb.Len = nodeSeq, len(nodeSeq)
	nodeBnt := constructcf.GetReadBntKmer(nodeSeq, 0, kmerlen-1)
	revNdSeq := constructcf.ReverseComplet(nodeBnt)
	tmp.Seq, tmp.Len = nd.Seq, kmerlen-1
	ndBntSeq := constructcf.ExtendKmerBnt2Byte(tmp)
	if reflect.DeepEqual(nodeBnt.Seq, nd.Seq) {
		if direction == BACKWARD {
			if nd.EdgeIDIncoming[base] > 1 {
				ID = nd.EdgeIDIncoming[base]
			} else {
				if correct {
					ea := GetComingOtherEArr(nd.EdgeIDIncoming, eID)
					if len(ea) == 1 {
						ID = ea[0]
					} else {
						return
					}
				} else {
					return
				}
			}
			ne := edgesArr[ID]
			if ne.EndNID == nd.ID && reflect.DeepEqual(nodeSeq, ne.Utg.Ks[len(ne.Utg.Ks)-(kmerlen-1):]) {
				pos = len(ne.Utg.Ks) - (kmerlen - 1)
				strand = PLUS
				return
			}
			if ne.StartNID == nd.ID && reflect.DeepEqual(nodeSeq, GetReverseCompByteArr(ne.Utg.Ks[:kmerlen-1])) {
				pos = kmerlen - 1
				strand = MINUS
			} else {
				log.Fatalf("[GetSelfCycleNextMapEdgeInfo] ne of node seq != nodeSeq, ne: %v\n\tnodeSeq: %v", ne, nodeSeq)
			}
		} else { // direction == FORWARD
			if nd.EdgeIDOutcoming[base] > 1 {
				ID = nd.EdgeIDOutcoming[base]
			} else {
				if correct {
					ea := GetComingOtherEArr(nd.EdgeIDOutcoming, eID)
					if len(ea) == 1 {
						ID = ea[0]
					} else {
						return
					}
				} else {
					return
				}
			}
			ne := edgesArr[ID]
			if ne.StartNID == nd.ID && reflect.DeepEqual(nodeSeq, ne.Utg.Ks[:kmerlen-1]) {
				pos = kmerlen - 1
				strand = PLUS
				return
			}
			if ne.EndNID == nd.ID && reflect.DeepEqual(nodeSeq, GetReverseCompByteArr(ne.Utg.Ks[len(ne.Utg.Ks)-(kmerlen-1):])) {
				pos = len(ne.Utg.Ks) - (kmerlen - 1)
				strand = MINUS
			} else {
				log.Fatalf("[GetSelfCycleNextMapEdgeInfo] ne of node seq != nodeSeq, ne: %v\n\tnodeSeq: %v", ne, nodeSeq)
			}
		}
	} else if reflect.DeepEqual(revNdSeq.Seq, nd.Seq) {
		base = bnt.BntRev[base]
		if direction == BACKWARD {
			if nd.EdgeIDOutcoming[base] > 1 {
				ID = nd.EdgeIDOutcoming[base]
			} else {
				if correct {
					ea := GetComingOtherEArr(nd.EdgeIDOutcoming, eID)
					if len(ea) == 1 {
						ID = ea[0]
					} else {
						return
					}
				} else {
					return
				}
			}
			ne := edgesArr[ID]
			if ne.StartNID == nd.ID && reflect.DeepEqual(nodeSeq, GetReverseCompByteArr(ne.Utg.Ks[:kmerlen-1])) {
				pos = kmerlen - 1
				strand = MINUS
				return
			}
			if ne.EndNID == nd.ID && reflect.DeepEqual(nodeSeq, ne.Utg.Ks[len(ne.Utg.Ks)-(kmerlen-1):]) {
				pos = len(ne.Utg.Ks) - (kmerlen - 1)
				strand = PLUS
			} else {
				log.Fatalf("[GetSelfCycleNextMapEdgeInfo] ne of node revNdSeq != nodeSeq, ne: %v\n\tnodeSeq: %v", ne, revNdSeq)
			}
		} else { // direction == FORWARD
			if nd.EdgeIDIncoming[base] > 1 {
				ID = nd.EdgeIDIncoming[base]
			} else {
				if correct {
					ea := GetComingOtherEArr(nd.EdgeIDIncoming, eID)
					if len(ea) == 1 {
						ID = ea[0]
					} else {
						return
					}
				} else {
					return
				}
			}
			ne := edgesArr[ID]
			if ne.EndNID == nd.ID && reflect.DeepEqual(nodeSeq, GetReverseCompByteArr(ne.Utg.Ks[len(ne.Utg.Ks)-(kmerlen-1):])) {
				pos = len(ne.Utg.Ks) - (kmerlen - 1)
				strand = MINUS
				return
			}
			if ne.StartNID == nd.ID && reflect.DeepEqual(nodeSeq, ne.Utg.Ks[:kmerlen-1]) {
				pos = kmerlen - 1
				strand = PLUS
			} else {
				log.Fatalf("[GetSelfCycleNextMapEdgeInfo] ne of node revNdSeq != nodeSeq, ne: %v\n\tnodeSeq: %v", ne, revNdSeq)
			}
		}
	} else {
		log.Fatalf("[GetSelfCycleNextMapEdgeInfo] nd.Seq != nodeSeq, nd.Seq: %v\n\tnodeSeq: %v", ndBntSeq, nodeSeq)
	}

	return
}

func MappingReadToEdgesBackWard(dk DBGKmer, ri constructcf.ReadInfo, rpos int, rstrand bool, edgesArr []DBGEdge, nodesArr []DBGNode, kmerlen int, correct bool) (errorNum int, ai AlignInfo) {
	var strand bool
	if dk.Strand != rstrand {
		dk.Pos += kmerlen
		strand = MINUS
	} else {
		strand = PLUS
	}

	ai.Seq = make([]byte, rpos+kmerlen)
	copy(ai.Seq[rpos:rpos+kmerlen], ri.Seq[rpos:rpos+kmerlen])
	ai.EndPos = math.MinInt32
	if rpos == 0 {
		var ps PathSeq
		ps.Start = int(dk.Pos) + kmerlen - 1
		ps.End = int(dk.Pos) - 1
		ai.Paths = append(ai.Paths, dk.ID)
		ai.Strands = append(ai.Strands, strand)
		if strand == PLUS {
			ai.EndPos = int(dk.Pos) - 1
		} else {
			ai.EndPos = int(dk.Pos)
		}
		return
	}
	for i := rpos - 1; i >= 0; {
		e := edgesArr[dk.ID]
		ai.Paths = append(ai.Paths, e.ID)
		ai.Strands = append(ai.Strands, strand)
		b := 0
		var j int
		if strand == PLUS {
			if int(dk.Pos) < rpos {
				b = rpos - int(dk.Pos)
			}
			j = int(dk.Pos) - 1
			for ; i >= b; i-- {
				//fmt.Printf("[MappingReadToEdgesBackWard] i: %v, j: %v, dk.Pos: %v, pos: %v, b: %v\n", i, j, dk.Pos, rpos, b)
				if ri.Seq[i] != e.Utg.Ks[j] {
					errorNum++
					if !correct {
						break
					}
				}
				ai.Seq[i] = e.Utg.Ks[j]
				j--
			}
		} else { // strand == MINUS
			if len(e.Utg.Ks)-int(dk.Pos) < rpos {
				b = rpos - (len(e.Utg.Ks) - int(dk.Pos))
			}
			j = int(dk.Pos)
			for ; i >= b; i-- {
				//fmt.Printf("[MappingReadToEdgesBackWard] i: %v, j: %v, dk.Pos: %v, pos: %v, b: %v\n", i, j, dk.Pos, rpos, b)
				if ri.Seq[i] != bnt.BntRev[e.Utg.Ks[j]] {
					errorNum++
					if !correct {
						break
					}
				}
				ai.Seq[i] = bnt.BntRev[e.Utg.Ks[j]]
				j++
			}
		}

		if !correct && errorNum > 0 {
			fmt.Printf("[MappingReadToEdgesBackWard]not perfect start i: %v,edge ID: %v,len(e.Utg.Ks):%v,  dk.Pos: %v, pos: %v, b: %v\n", i, dk.ID, len(e.Utg.Ks), dk.Pos, rpos, b)
			break
		}
		//fmt.Printf("[paraMapNGS2DBG] i: %v, j: %v, dk.Pos: %v, pos: %v, b: %v\n", i, j, dk.Pos, rpos, b)

		if i < 0 {
			ai.EndPos = j
			break
		}
		// find next edge
		rpos = i + 1
		var node DBGNode
		//var base byte
		// if is a self cycle edge
		if e.StartNID == e.EndNID {
			var pos int
			dk.ID, pos, strand = GetSelfCycleNextMapEdgeInfo(e.ID, nodesArr[e.StartNID], edgesArr, ai.Seq[rpos:rpos+(kmerlen-1)], kmerlen, BACKWARD, ri.Seq[i], correct)
			if dk.ID <= 1 {
				ai.EndPos = j
				ai.Seq = ai.Seq[rpos:]
				break
			}
			dk.Pos = pos
			continue
		}

		if strand == PLUS {
			if e.StartNID == 0 {
				break
			}
			node = nodesArr[e.StartNID]
			base := bnt.BntRev[ri.Seq[i]]
			if IsInComing(node.EdgeIDOutcoming, e.ID) && node.EdgeIDIncoming[ri.Seq[i]] > 1 {
				dk.ID = node.EdgeIDIncoming[ri.Seq[i]]
			} else if IsInComing(node.EdgeIDIncoming, e.ID) && node.EdgeIDOutcoming[base] > 1 {
				dk.ID = node.EdgeIDOutcoming[base]
			} else {
				if correct {
					ea := GetNearEdgeIDArr(node, e.ID)
					if len(ea) == 1 {
						dk.ID = ea[0]
					} else {
						ai.EndPos = j
						ai.Seq = ai.Seq[rpos:]
						fmt.Printf("[MappingReadToEdgesBackWard] not found next edge in node: %v\n", node)
						break
					}
				} else {
					fmt.Printf("[MappingReadToEdgesBackWard] not found next edge in node: %v\n", node)
					break
				}
			}
			ne := edgesArr[dk.ID]
			if ne.StartNID == ne.EndNID {
				nodeSeq := ai.Seq[rpos : rpos+(kmerlen-1)]
				if reflect.DeepEqual(nodeSeq, ne.Utg.Ks[len(ne.Utg.Ks)-(kmerlen-1):]) {
					dk.Pos = len(ne.Utg.Ks) - (kmerlen - 1)
					strand = PLUS
				} else if reflect.DeepEqual(nodeSeq, GetReverseCompByteArr(ne.Utg.Ks[:kmerlen-1])) {
					dk.Pos = kmerlen - 1
					strand = MINUS
				}
			} else {
				if ne.EndNID == node.ID {
					dk.Pos = len(ne.Utg.Ks) - (kmerlen - 1)
				} else {
					dk.Pos = kmerlen - 1
					strand = !strand
				}
			}
		} else { // strand == MINUS
			if e.EndNID == 0 {
				break
			}
			node = nodesArr[e.EndNID]
			base := bnt.BntRev[ri.Seq[i]]
			if IsInComing(node.EdgeIDIncoming, e.ID) && node.EdgeIDOutcoming[base] > 1 {
				dk.ID = node.EdgeIDOutcoming[base]
			} else if IsInComing(node.EdgeIDOutcoming, e.ID) && node.EdgeIDIncoming[ri.Seq[i]] > 1 {
				dk.ID = node.EdgeIDIncoming[ri.Seq[i]]
			} else {
				if correct {
					ea := GetNearEdgeIDArr(node, e.ID)
					if len(ea) == 1 {
						dk.ID = ea[0]
					} else {
						ai.EndPos = j
						ai.Seq = ai.Seq[rpos:]
						fmt.Printf("[MappingReadToEdgesBackWard] not found next edge in node: %v\n", node)
						break
					}
				} else {
					fmt.Printf("[MappingReadToEdgesBackWard] not found next edge in node: %v\n", node)
					break
				}
			}
			ne := edgesArr[dk.ID]
			if ne.StartNID == ne.EndNID {
				nodeSeq := ai.Seq[rpos : rpos+(kmerlen-1)]
				if reflect.DeepEqual(nodeSeq, ne.Utg.Ks[len(ne.Utg.Ks)-(kmerlen-1):]) {
					dk.Pos = len(ne.Utg.Ks) - (kmerlen - 1)
					strand = PLUS
				} else if reflect.DeepEqual(nodeSeq, GetReverseCompByteArr(ne.Utg.Ks[:kmerlen-1])) {
					dk.Pos = kmerlen - 1
					strand = MINUS
				}
			} else {
				if ne.StartNID == node.ID {
					dk.Pos = kmerlen - 1
				} else {
					dk.Pos = len(ne.Utg.Ks) - (kmerlen - 1)
					strand = !strand
				}
			}
		}
	}
	return
}*/

/*func MappingReadToEdgesForWard(dk DBGKmer, ri constructcf.ReadInfo, rpos int, rstrand bool, edgesArr []DBGEdge, nodesArr []DBGNode, Kmerlen int, correct bool) (errorNum int, ai AlignInfo) {
	var strand bool
	if dk.Strand == rstrand {
		dk.Pos += Kmerlen
		strand = PLUS
	} else {
		strand = MINUS
	}
	ai.Seq = make([]byte, len(ri.Seq))
	startPos := rpos - Kmerlen
	copy(ai.Seq[rpos-Kmerlen:rpos], ri.Seq[rpos-Kmerlen:rpos])
	ai.EndPos = math.MinInt32
	for i := rpos; i < len(ri.Seq); {
		e := edgesArr[dk.ID]
		ai.Paths = append(ai.Paths, e.ID)
		ai.Strands = append(ai.Strands, strand)
		//fmt.Printf("[MappingReadToEdgesForWard]ai: %v\n", ai)
		b := len(ri.Seq)
		var j int
		if strand == PLUS {
			if len(e.Utg.Ks)-int(dk.Pos) < len(ri.Seq)-rpos {
				b = rpos + (len(e.Utg.Ks) - int(dk.Pos))
			}
			j = int(dk.Pos)
			for ; i < b; i++ {
				//fmt.Printf("[MappingReadToEdgesForWard] i: %v, j: %v, dk.Pos: %v, pos: %v, b: %v\n", i, j, dk.Pos, rpos, b)
				if ri.Seq[i] != e.Utg.Ks[j] {
					errorNum++
					if !correct {
						break
					}
				}
				ai.Seq[i] = e.Utg.Ks[j]
				j++
			}
		} else { // strand == MINUS
			if len(ri.Seq)-rpos > int(dk.Pos) {
				b = rpos + int(dk.Pos)
			}
			j = int(dk.Pos) - 1
			for ; i < b; i++ {
				//fmt.Printf("[MappingReadToEdgesForWard] i: %v, j: %v, dk.Pos: %v, pos: %v, b: %v\n", i, j, dk.Pos, rpos, b)
				if ri.Seq[i] != bnt.BntRev[e.Utg.Ks[j]] {
					errorNum++
					if !correct {
						break
					}
				}
				ai.Seq[i] = bnt.BntRev[e.Utg.Ks[j]]
				j--
			}
		}

		if !correct && errorNum > 0 {
			fmt.Printf("[MappingReadToEdgesForWard]not perfect end i: %v,edge ID: %v,len(e.Utg.Ks): %v,  dk.Pos: %v, pos: %v, b: %v\n", i, dk.ID, len(e.Utg.Ks), dk.Pos, rpos, b)
			break
		}

		//fmt.Printf("[MappingReadToEdgesForWard]after alignment i: %v, j: %v, dk.Pos: %v, pos: %v, b: %v\n", i, j, dk.Pos, rpos, b)
		if i >= len(ri.Seq) {
			ai.EndPos = j
			ai.Seq = ai.Seq[startPos:]
			break
		}

		// find next edge
		rpos = i
		var node DBGNode
		var base byte
		// if is a self cycle edge
		if e.StartNID == e.EndNID {
			var pos int
			dk.ID, pos, strand = GetSelfCycleNextMapEdgeInfo(e.ID, nodesArr[e.StartNID], edgesArr, ai.Seq[rpos-(Kmerlen-1):rpos], Kmerlen, FORWARD, ri.Seq[i], correct)
			if dk.ID <= 1 {
				ai.EndPos = j
				break
			}
			dk.Pos = int(pos)
			continue
		}
		if strand == PLUS {
			if e.EndNID == 0 {
				break
			}
			node = nodesArr[e.EndNID]
			base = bnt.BntRev[ri.Seq[i]]
			if IsInComing(node.EdgeIDIncoming, e.ID) && node.EdgeIDOutcoming[ri.Seq[i]] > 1 {
				dk.ID = node.EdgeIDOutcoming[ri.Seq[i]]
			} else if IsInComing(node.EdgeIDOutcoming, e.ID) && node.EdgeIDIncoming[base] > 1 {
				dk.ID = node.EdgeIDIncoming[base]
			} else {
				if correct {
					ea := GetNearEdgeIDArr(node, e.ID)
					if len(ea) == 1 {
						dk.ID = ea[0]
					} else {
						ai.EndPos = j
						ai.Seq = ai.Seq[startPos:rpos]
						fmt.Printf("[MappingReadToEdgesForWard] not found next edge in node: %v\n", node)
						break
					}
				} else {
					fmt.Printf("[MappingReadToEdgesForWard] not found next edge in node: %v\n", node)
					break
				}
			}
			ne := edgesArr[dk.ID]
			if ne.StartNID == ne.EndNID {
				nodeSeq := ai.Seq[rpos-(Kmerlen-1) : rpos]
				if reflect.DeepEqual(nodeSeq, ne.Utg.Ks[:Kmerlen-1]) {
					dk.Pos = Kmerlen - 1
					strand = PLUS
				} else if reflect.DeepEqual(nodeSeq, GetReverseCompByteArr(ne.Utg.Ks[len(ne.Utg.Ks)-(Kmerlen-1):])) {
					dk.Pos = len(ne.Utg.Ks) - (Kmerlen - 1)
					strand = MINUS
				} else {
					log.Fatalf("[MappingReadToEdgesForWard] ne: %v not found proper start position\n", ne)
				}
			} else {
				if ne.StartNID == node.ID {
					dk.Pos = Kmerlen - 1
				} else {
					dk.Pos = len(ne.Utg.Ks) - (Kmerlen - 1)
					strand = !strand
				}
			}
		} else { // strand == MINUS
			if e.StartNID == 0 {
				break
			}
			node = nodesArr[e.StartNID]
			base = bnt.BntRev[ri.Seq[i]]
			if IsInComing(node.EdgeIDOutcoming, e.ID) && node.EdgeIDIncoming[base] > 1 {
				dk.ID = node.EdgeIDIncoming[base]
			} else if IsInComing(node.EdgeIDIncoming, e.ID) && node.EdgeIDOutcoming[ri.Seq[i]] > 1 {
				dk.ID = node.EdgeIDOutcoming[ri.Seq[i]]
			} else {
				if correct {
					ea := GetNearEdgeIDArr(node, e.ID)
					if len(ea) == 1 {
						dk.ID = ea[0]
					} else {
						ai.EndPos = j
						ai.Seq = ai.Seq[startPos:rpos]
						fmt.Printf("[MappingReadToEdgesForWard] not found next edge in node: %v\n", node)
						break
					}
				} else {
					fmt.Printf("[MappingReadToEdgesForWard] not found next edge in node: %v\n", node)
					break
				}
			}
			ne := edgesArr[dk.ID]
			if ne.StartNID == ne.EndNID {
				nodeSeq := ai.Seq[rpos-(Kmerlen-1) : rpos]
				if reflect.DeepEqual(nodeSeq, ne.Utg.Ks[:Kmerlen-1]) {
					dk.Pos = Kmerlen - 1
					strand = PLUS
				} else if reflect.DeepEqual(nodeSeq, GetReverseCompByteArr(ne.Utg.Ks[len(ne.Utg.Ks)-(Kmerlen-1):])) {
					dk.Pos = len(ne.Utg.Ks) - (Kmerlen - 1)
					strand = MINUS
				} else {
					log.Fatalf("[MappingReadToEdgesForWard] ne: %v not found proper start position\n", ne)
				}
			} else {
				if ne.EndNID == node.ID {
					dk.Pos = len(ne.Utg.Ks) - (Kmerlen - 1)
				} else {
					dk.Pos = Kmerlen - 1
					strand = !strand
				}
			}
		}
	}
	return
}*/

// parallel Map NGS reads to the DBG edges, then output alignment path for the DBG
/*func paraMapNGS2DBG(cs chan constructcf.ReadInfo, wc chan AlignInfo, nodesArr []DBGNode, edgesArr []DBGEdge, cf CuckooFilter, winSize int) {
	var notFoundSeedNum, mapOneEdgeNum, notPerfectNum int
	for {
		ri, ok := <-cs
		var ai AlignInfo
		if !ok {
			wc <- ai
			break
		}

		// found kmer seed position in the DBG edges
		dbgK, pos, strand := LocateSeedKmerCF(cf, ri, winSize, edgesArr)

		if dbgK.GetCount() == 0 { // not found in the cuckoofilter
			notFoundSeedNum++
			continue
		}

		// check can map two or more edges
		{
			el := len(edgesArr[dbgK.ID].Utg.Ks)
			if dbgK.Strand == strand {
				if dbgK.Pos > int(pos) && len(ri.Seq)-int(pos) < el-dbgK.Pos {
					mapOneEdgeNum++
					continue
				}
			} else {
				if int(pos) < el-(int(dbgK.Pos)+cf.Kmerlen) && len(ri.Seq)-(int(pos)+cf.Kmerlen) < int(dbgK.Pos) {
					mapOneEdgeNum++
					continue
				}
			}
		}

		// extend map to the DBG edges
		{
			// overAll := true // note if overall length read sequence match
			// map the start partition of read sequence
			errorNum, ar := MappingReadToEdgesBackWard(dbgK, ri, int(pos), strand, edgesArr, nodesArr, cf.Kmerlen, false)
			if errorNum > 0 {
				notPerfectNum++
				continue
			}
			ar.Paths = ReverseDBG_MAX_INTArr(ar.Paths)
			// map the end partition of read sequence
			errorNum, al := MappingReadToEdgesForWard(dbgK, ri, int(pos)+cf.Kmerlen, strand, edgesArr, nodesArr, cf.Kmerlen, false)

			if errorNum > 0 {
				//fmt.Printf("[paraMapNGS2DBG] ", a)
				notPerfectNum++
				continue
			}

			if len(ar.Paths) > 0 && len(al.Paths) > 0 && ar.Paths[len(ar.Paths)-1] == al.Paths[0] {
				ai.Paths = append(ar.Paths, al.Paths[1:]...)
				ai.ID = ri.ID
			} else {
				log.Fatalf("[paraMapNGS2DBG] ar: %v and al: %v not consis\n", ar, al)
			}

			// write to output
			if len(ai.Paths) > 1 {
				wc <- ai
			}
		}
	}
	fmt.Printf("[paraMapNGS2DBG] not found seed reads number is : %v\n", notFoundSeedNum)
	fmt.Printf("[paraMapNGS2DBG] map one edge reads number is : %v\n", mapOneEdgeNum)
	fmt.Printf("[paraMapNGS2DBG] not perfect mapping reads number is : %v\n", notPerfectNum)

	ai.ID = ri.ID
	e := edgesArr[dbgK.ID]
	var overAll bool // note if overall length match
	for i < len(ri.Seq)-k {
		var plus bool // if the read map to the edge plus strand
		var n DBGNode
		ek := e.Utg.Ks[dbgK.Pos : dbgK.Pos+uint32(k)]
		if reflect.DeepEqual(ri.Seq[i:i+k], ek) {
			plus = true
		} else if reflect.DeepEqual(ri.Seq[i:i+k], GetReverseCompByteArr(ek)) {
			plus = false
		} else { // not equal for the DBG edge
			log.Fatalf("[paraMapNGS2DBG] read kmer not equal to the DBG edge\nread info: %v\nedge info: %v\n", ri.Seq[i:i+k], ek)
		}
		//if output {
		//	fmt.Printf("[paraMapNGS2DBG] i : %v\tdbgK: %v\n\tread info: %v\n\tedge seq: %v\n\tedge info: ID: %v, len: %v, startNID: %v, endNID: %v\n", i, dbgK, ri.Seq[i:i+k], ek, e.ID, len(e.Utg.Ks), e.StartNID, e.EndNID)
		//}

		x := i + k
		if plus {
			y := int(dbgK.Pos) + k
			for ; y < len(e.Utg.Ks) && x < len(ri.Seq); y++ {
				if ri.Seq[x] != e.Utg.Ks[y] {
					break
				}
				x++
			}
			//if output {
			//	fmt.Printf("[paraMapNGS2DBG]plus i: %v,x: %v, dbgK: %v\n", i, x, dbgK)
			//}
			if x >= len(ri.Seq) {
				ai.Paths = append(ai.Paths, dbgK.ID)
				overAll = true
				break
			}
			if y < len(e.Utg.Ks) {
				break
			}
			ai.Paths = append(ai.Paths, dbgK.ID)
			// set next edge info
			if e.EndNID == 0 {
				break
			}
			n = nodesArr[e.EndNID]
			if IsInComing(n.EdgeIDIncoming, e.ID) {
				if n.EdgeIDOutcoming[ri.Seq[x]] > 0 {
					dbgK.ID = n.EdgeIDOutcoming[ri.Seq[x]]
				} else {
					break
				}
			} else {
				b := bnt.BntRev[ri.Seq[x]]
				//b := ri.Seq[x]
				if n.EdgeIDIncoming[b] > 0 {
					dbgK.ID = n.EdgeIDIncoming[b]
				} else {
					break
				}
			}

		} else { // if the strand as minus
			y := int(dbgK.Pos) - 1
			for ; y >= 0 && x < len(ri.Seq); y-- {
				b := bnt.BntRev[e.Utg.Ks[y]]
				if ri.Seq[x] != b {
					break
				}
				x++
			}
			//if output {
			//	fmt.Printf("[paraMapNGS2DBG]i: %v, x: %v, dbgK: %v\n", i, x, dbgK)
			//}
			if x >= len(ri.Seq) {
				ai.Paths = append(ai.Paths, dbgK.ID)
				break
			}
			if y >= 0 {
				break
			}
			ai.Paths = append(ai.Paths, dbgK.ID)
			// set next edge info
			if e.StartNID == 0 {
				break
			}
			n = nodesArr[e.StartNID]
			if IsInComing(n.EdgeIDOutcoming, e.ID) {
				//fmt.Printf("[paraMapNGS2DBG]x: %v, len(ri.Seq): %v\tri.Seq[x]: %v\n", x, len(ri.Seq), ri.Seq[x])
				b := bnt.BntRev[ri.Seq[x]]
				if n.EdgeIDIncoming[b] > 0 {
					dbgK.ID = n.EdgeIDIncoming[b]
				} else {
					break
				}
			} else {
				if n.EdgeIDOutcoming[ri.Seq[x]] > 0 {
					dbgK.ID = n.EdgeIDOutcoming[ri.Seq[x]]
				} else {
					break
				}
			}
		}

		e = edgesArr[dbgK.ID]
		if e.StartNID == n.ID {
			dbgK.Pos = 0
		} else {
			dbgK.Pos = uint32(len(e.Utg.Ks) - k)
		}
		i = x - k + 1
	}
}*/

/*func MapNGS2DBG(opt Options, nodesArr []DBGNode, edgesArr []DBGEdge, wrFn string) {
	// construct cuckoofilter of DBG sample
	cfSize := GetCuckoofilterDBGSampleSize(edgesArr, int64(opt.WinSize), int64(opt.MaxNGSReadLen), int64(opt.Kmer))
	fmt.Printf("[MapNGS2DBG] cfSize: %v\n", cfSize)
	cf := MakeCuckooFilter(uint64(cfSize*5), opt.Kmer)
	fmt.Printf("[MapNGS2DBG] cf.numItems: %v\n", cf.NumItems)
	count := ConstructCFDBGMinimizers(cf, edgesArr, opt.WinSize, int64(opt.MaxNGSReadLen))
	fmt.Printf("[MapNGS2DBG]construct Smaple of DBG edges cuckoofilter number is : %v\n", count)
	//if cfSize != int64(count) {
	//	log.Fatalf("[MapNGS2DBG]cfSize : %v != count : %v, please check\n", cfSize, count)
	//}

	numCPU := opt.NumCPU
	runtime.GOMAXPROCS(numCPU + 2)
	bufSize := 66000
	cs := make(chan constructcf.ReadInfo, bufSize)
	wc := make(chan AlignInfo, numCPU*20)
	defer close(wc)
	//defer close(cs)

	// Load NGS read from cfg
	fn := opt.CfgFn
	readT := (numCPU + 5 - 1) / 5
	go LoadNGSReads(fn, opt.Correct, cs, readT, opt.Kmer)
	for i := 0; i < numCPU; i++ {
		go paraMapNGS2DBG(cs, wc, nodesArr, edgesArr, cf, opt.WinSize)
	}

	// write function
	writeAlignToFile(wrFn, wc, numCPU)
}*/

func AtoiArr(sa []string) []DBG_MAX_INT {
	da := make([]DBG_MAX_INT, len(sa))
	for i, e := range sa {
		d, err := strconv.Atoi(e)
		if err != nil {
			log.Fatalf("[AtoiArr] string: %v convert to integer, err: %v\n", e, err)
		}
		da[i] = DBG_MAX_INT(d)
	}
	return da
}

// add to the DBGEdge pathMat
func AddPathToDBGEdge(edgesArr []DBGEdge, mapNGSFn string) {
	fp, err := os.Open(mapNGSFn)
	if err != nil {
		log.Fatalf("[AddPathToDBGEdge] open file: '%v' error, err : %v\n", mapNGSFn, err)
	}
	buffp := bufio.NewReaderSize(fp, 1<<25)
	defer fp.Close()
	m, err := buffp.ReadString('\n')
	for ; err == nil; m, err = buffp.ReadString('\n') {
		sa := strings.Split(m[:len(m)-1], "\t")
		pa := strings.Split(sa[1], ":")
		da := AtoiArr(pa)
		//fmt.Printf("sa: %v\npa: %v\n", sa, pa)
		for i, eID := range da {
			// statistics coverage depth for every reads pass the start and end position
			cd := uint16(1)
			if 0 < i && i < len(da)-1 {
				cd = 2
			}
			if edgesArr[eID].CovD < math.MaxUint16 {
				edgesArr[eID].CovD += cd
			}
			if edgesArr[eID].GetUniqueFlag() > 0 {
				edgesArr[eID].InsertPathToEdge(da, 1)
			}
		}
	}
	if err != io.EOF {
		log.Fatalf("[AddPathToDBGEdge] Failed to read file: %v, err : %v\n", mapNGSFn, err)
	}

	// at last edge average coverage = (start position coverage depth +  end position coverage depth)/2
	for i := 2; i < len(edgesArr); i++ {
		edgesArr[i].CovD /= 2
	}

}

func FindConsisPath(pID DBG_MAX_INT, e DBGEdge) (consisP Path) {
	var pm []Path
	for _, pa := range e.PathMat {
		p1 := IndexEID(pa.IDArr, pID)
		p2 := IndexEID(pa.IDArr, e.ID)
		if p1 >= 0 {
			var npa Path
			if p1 < p2 {
				npa.IDArr = GetReverseDBG_MAX_INTArr(pa.IDArr[:p2+1])
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
}

// coming denote node EdgeIDcoming, true is Outcoming, false is Incoming
func GetNearEdgeIDArr(nd DBGNode, eID DBG_MAX_INT, coming bool) (eArr []DBG_MAX_INT) {
	if eID <= 0 {
		log.Fatalf("[GetNearEdgeIDArr] eID must bigger than zero, eID: %v\n", eID)
	}
	if nd.ID < 2 {
		return
	}
	var ok bool
	if coming {
		for _, id := range nd.EdgeIDIncoming {
			if id == eID {
				ok = true
				break
			}
		}
	} else {
		for _, id := range nd.EdgeIDOutcoming {
			if id == eID {
				ok = true
				break
			}
		}
	}
	if !ok {
		log.Fatalf("[GetNearEdgeIDArr] coming: %v, not found eID: %v, in nd: %v\n", coming, eID, nd)
	}

	if coming {
		for _, id := range nd.EdgeIDOutcoming {
			if id > 1 {
				eArr = append(eArr, id)
			}
		}
	} else {
		for _, id := range nd.EdgeIDIncoming {
			if id > 1 {
				eArr = append(eArr, id)
			}
		}
	}

	return eArr
}

func GetNextEdgeStrand(e, ne *DBGEdge, nodesArr []DBGNode, strand bool) bool {
	var nd DBGNode
	if (e.StartNID == ne.StartNID) || (e.StartNID == ne.EndNID) {
		nd = nodesArr[e.StartNID]
	} else {
		nd = nodesArr[e.EndNID]
	}
	if e.ID == ne.ID {
		return strand
	}
	if e.StartNID == e.EndNID {
		if ne.StartNID != ne.EndNID {
			if CountEdgeIDComing(nd.EdgeIDOutcoming, e.ID) == 2 {
				if IsInComing(nd.EdgeIDIncoming, ne.ID) && (ne.EndNID == e.EndNID) {
					strand = !strand
				}
			} else if CountEdgeIDComing(nd.EdgeIDIncoming, e.ID) == 2 {
				if IsInComing(nd.EdgeIDOutcoming, ne.ID) && (ne.EndNID == e.EndNID) {
					strand = !strand
				}
			} else {
				if (IsInComing(nd.EdgeIDIncoming, ne.ID) && ne.StartNID == e.StartNID) || (IsInComing(nd.EdgeIDOutcoming, ne.ID) && ne.EndNID == e.EndNID) {
					strand = !strand
				}
			}
		}
	} else {
		if ne.StartNID == ne.EndNID {
			if CountEdgeIDComing(nd.EdgeIDOutcoming, ne.ID) == 2 {
				if IsInComing(nd.EdgeIDIncoming, e.ID) && ne.StartNID == e.StartNID {
					strand = !strand
				}
			} else if CountEdgeIDComing(nd.EdgeIDIncoming, ne.ID) == 2 {
				if IsInComing(nd.EdgeIDOutcoming, e.ID) && ne.StartNID == e.StartNID {
					strand = !strand
				}
			} else if (IsInComing(nd.EdgeIDIncoming, e.ID) && (e.StartNID == ne.StartNID)) || (IsInComing(nd.EdgeIDOutcoming, e.ID) && e.EndNID == ne.EndNID) {
				strand = !strand
			}
		} else {
			if (e.StartNID == ne.StartNID) || (e.EndNID == ne.EndNID) {
				strand = !strand
			}
		}
	}
	return strand
}

func GetNextEdgeStrand2(e, ne *DBGEdge, strand bool, direction uint8) bool {
	if e.ID == ne.ID {
		return strand
	}
	if e.StartNID == e.EndNID {
		if CountEdgeIDComing(ne.EdgeIDOutcoming, e.ID) == 2 {
			if e.StartNID == ne.StartNID {
				strand = !strand
			}
		} else if CountEdgeIDComing(ne.EdgeIDIncoming, e.ID) == 2 {
			if e.EndNID == ne.EndNID {
				strand = !strand
			}
		} else {
			if IsInComing(e.EdgeIDOutcoming, ne.ID) {
				if e.EndNID == ne.EndNID {
					strand = !strand
				}
			} else {
				if e.StartNID == ne.StartNID {
					strand = !strand
				}
			}
		}
	} else if ne.StartNID == ne.EndNID {
		if CountEdgeIDComing(e.EdgeIDOutcoming, ne.ID) == 2 {
			if e.StartNID == ne.StartNID {
				strand = !strand
			}
		} else if CountEdgeIDComing(e.EdgeIDIncoming, ne.ID) == 2 {
			if e.EndNID == ne.EndNID {
				strand = !strand
			}
		} else {
			if IsInComing(ne.EdgeIDOutcoming, e.ID) {
				if e.EndNID == ne.EndNID {
					strand = !strand
				}
			} else {
				if e.StartNID == ne.StartNID {
					strand = !strand
				}
			}
		}
	} else {
		if direction == FORWARD {
			if strand {
				if e.EndNID == ne.EndNID {
					strand = !strand
				}
			} else {
				if e.StartNID == ne.StartNID {
					strand = !strand
				}
			}
		} else {
			if strand {
				if e.StartNID == ne.StartNID {
					strand = !strand
				}
			} else {
				if e.EndNID == ne.EndNID {
					strand = !strand
				}
			}
		}
	}

	return strand
}

func FreqNumInDBG_MAX_INTArr(arr []DBG_MAX_INT, eID DBG_MAX_INT) (count int) {
	for _, id := range arr {
		if id == eID {
			count++
		}
	}

	return count
}

// merge DBGEdge's pathMat
/*func MergePathMat(edgesArr []DBGEdge, nodesArr []DBGNode, minMapFreq int) {
	for i, e := range edgesArr {
		//if e.GetUniqueFlag() > 0 {
		//	fmt.Printf("[MergePathMat]unique edge : %v\n", e)
		//}
		if i < 2 || e.GetDeleteFlag() > 0 || len(e.PathMat) == 0 || (e.GetUniqueFlag() == 0 && e.GetSemiUniqueFlag() == 0) {
			continue
		}

		// debug code
		if e.PathMat[0].Freq == 0 {
			log.Fatalf("[MergePathMat] e.ID: %v PathMat[0]: %v\n", e.ID, e.PathMat[0])
		}
		//fmt.Printf("[MergePathMat]e.PathMat : %v\n", e.PathMat)

		// check if is a node cycle edge
		if e.StartNID == e.EndNID {
			node := nodesArr[e.StartNID]
			n1, _ := GetEdgeIDComing(node.EdgeIDIncoming)
			n2, _ := GetEdgeIDComing(node.EdgeIDOutcoming)
			if n1 > 1 || n2 > 1 {
				edgesArr[e.ID].ResetUniqueFlag()
				edgesArr[e.ID].ResetSemiUniqueFlag()
			}
			continue
		}

		CheckPathDirection(edgesArr, nodesArr, e.ID)
		if edgesArr[e.ID].GetUniqueFlag() == 0 && edgesArr[e.ID].GetSemiUniqueFlag() == 0 {
			continue
		}
		if len(e.PathMat) == 1 {
			continue
		}

		// merge process
		var leftMax, rightMax int
		for _, p := range e.PathMat {
			idx := IndexEID(p.IDArr, e.ID)
			if idx > leftMax {
				leftMax = idx
			}
			if len(p.IDArr)-idx > rightMax {
				rightMax = len(p.IDArr) - idx
			}
		}
		al := leftMax + rightMax
		// copy e.PathMat to the pm
		pm := make([]Path, len(e.PathMat))
		for j, p := range e.PathMat {
			t := make([]DBG_MAX_INT, al)
			idx := IndexEID(p.IDArr, e.ID)
			copy(t[leftMax-idx:], p.IDArr)
			pm[j].Freq = p.Freq
			pm[j].IDArr = t
		}

		// alignment PathMat
		for j, p := range pm {
			fmt.Printf("[MergePathMat] pm[%v]: %v\n", j, p)
		}

		// find consis Path
		var path Path
		path.Freq = math.MaxInt
		// add left==0 and right==1 partition
		for z := 0; z < 2; z++ {
			var j, step int
			if z == 0 {
				j = leftMax
				step = -1
			} else {
				j = leftMax + 1
				step = 1
			}
			for ; ; j += step {
				if z == 0 {
					if j < 0 {
						break
					}
				} else {
					if j >= al {
						break
					}
				}
				suc := true // note found consis Path
				var freq int
				var id DBG_MAX_INT
				for k := 0; k < len(pm); k++ {
					if pm[k].IDArr[j] > 0 {
						if id == 0 {
							id = pm[k].IDArr[j]
							freq = pm[k].Freq
						} else {
							if id == pm[k].IDArr[j] {
								freq += pm[k].Freq
							} else {
								suc = false
								break
							}
						}
					}
				}
				if !suc {
					break
				}
				if freq >= minMapFreq {
					path.IDArr = append(path.IDArr, id)
					if path.Freq > freq {
						path.Freq = freq
					}
				} else {
					break
				}
			}

			if z == 0 {
				path.IDArr = GetReverseDBG_MAX_INTArr(path.IDArr)
			}
		}

		if path.Freq >= minMapFreq && len(path.IDArr) >= 2 {
			var pm []Path
			pm = append(pm, path)
			edgesArr[i].PathMat = pm
			fmt.Printf("[MergePathMat] edge ID: %v, merge path: %v\n", e.ID, edgesArr[i].PathMat)
		} else {
			edgesArr[i].PathMat = nil
		}
	}
}*/

func IsTwoEdgesCyclePath(edgesArr []DBGEdge, nodesArr []DBGNode, eID DBG_MAX_INT) bool {
	e := edgesArr[eID]
	if e.StartNID == 0 || e.EndNID == 0 {
		return false
	}
	var arr1, arr2 []DBG_MAX_INT
	if IsInComing(nodesArr[e.StartNID].EdgeIDIncoming, eID) {
		arr1 = GetNearEdgeIDArr(nodesArr[e.StartNID], eID, true)
	} else {
		arr1 = GetNearEdgeIDArr(nodesArr[e.StartNID], eID, false)
	}
	if IsInComing(nodesArr[e.EndNID].EdgeIDIncoming, eID) {
		arr2 = GetNearEdgeIDArr(nodesArr[e.EndNID], eID, true)
	} else {
		arr2 = GetNearEdgeIDArr(nodesArr[e.EndNID], eID, false)
	}
	if len(arr1) > 0 && len(arr2) > 0 {
		for _, id := range arr1 {
			for _, eID := range arr2 {
				if id == eID {
					return true
				}
			}
		}
	}

	return false
}
func GetOtherEdgeInTwoEdgesCycle(edgesArr []DBGEdge, nodesArr []DBGNode, eID DBG_MAX_INT) DBG_MAX_INT {
	e := edgesArr[eID]
	if e.StartNID == 0 || e.EndNID == 0 {
		return 0
	}
	var arr1, arr2 []DBG_MAX_INT
	if IsInComing(nodesArr[e.StartNID].EdgeIDIncoming, eID) {
		arr1 = GetNearEdgeIDArr(nodesArr[e.StartNID], eID, true)
	} else {
		arr1 = GetNearEdgeIDArr(nodesArr[e.StartNID], eID, false)
	}
	if IsInComing(nodesArr[e.EndNID].EdgeIDIncoming, eID) {
		arr2 = GetNearEdgeIDArr(nodesArr[e.EndNID], eID, true)
	} else {
		arr2 = GetNearEdgeIDArr(nodesArr[e.EndNID], eID, false)
	}
	//fmt.Printf("[GetOtherEdgeInTwoEdgesCycle]arr1: %v, arr2: %v\n", arr1, arr2)
	if len(arr1) > 0 && len(arr2) > 0 {
		for _, id := range arr1 {
			for _, eID := range arr2 {
				if id == eID {
					return id
				}
			}
		}
	}

	return 0
}

func ExtendPath(edgesArr []DBGEdge, nodesArr []DBGNode, e DBGEdge, minMappingFreq int, semi bool) (maxP Path) {
	var p Path
	p.IDArr = make([]DBG_MAX_INT, len(e.PathMat[0].IDArr))
	copy(p.IDArr, e.PathMat[0].IDArr)
	p.Freq = e.PathMat[0].Freq
	if p.Freq < minMappingFreq {
		return
	}
	idx := IndexEID(p.IDArr, e.ID)
	mutualArr := []DBG_MAX_INT{e.ID}
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
		if IsInDBG_MAX_INTArr(mutualArr, e2.ID) {
			fmt.Printf("[ExtendPath] e2.ID: %v in the mutualArr: %v\n", e2.ID, mutualArr)
			continue
		}
		p2 := e2.PathMat[0]
		eID1 := mutualArr[len(mutualArr)-1]
		e1 := edgesArr[eID1]
		j1 := IndexEID(p2.IDArr, e1.ID)
		j2 := IndexEID(p2.IDArr, e2.ID)
		fmt.Printf("[ExtendPath]BACKWARD e1 ID: %v,  e2 ID: %v, p2: %v\n", e1.ID, e2.ID, p2)
		if j1 < 0 || j2 < 0 {
			continue
		}
		if j1 < j2 {
			p2.IDArr = GetReverseDBG_MAX_INTArr(p2.IDArr)
			j1 = IndexEID(p2.IDArr, e1.ID)
			//j2 = IndexEID(p2.IDArr, e2.ID)
		}
		k1 := IndexEID(p.IDArr, e1.ID)
		//fmt.Printf("[ExtendPath]j1: %v, j2:%v,k1: %v, eID1: %v, eID2: %v, p: %v, p2: %v\n", j1, j2, k1, e1.ID, e2.ID, p, p2)
		if j1 >= k1 && len(p2.IDArr)-j1 <= len(p.IDArr)-k1 {
			if reflect.DeepEqual(p2.IDArr[j1-k1:], p.IDArr[:len(p2.IDArr)-(j1-k1)]) {
				mutualArr = append(mutualArr, e2.ID)
				na := make([]DBG_MAX_INT, len(p.IDArr)+j1-k1)
				copy(na[:j1], p2.IDArr[:j1])
				copy(na[j1:], p.IDArr[k1:])
				p.IDArr = na
				if p2.Freq < p.Freq {
					p.Freq = p2.Freq
				}
				idx = IndexEID(p.IDArr, e2.ID)
				i = -1
			}
		}
		fmt.Printf("[ExtendPath] p: %v\n", p)
	}

	ReverseDBG_MAX_INTArr(mutualArr)
	idx = IndexEID(p.IDArr, e.ID)
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
		if IsInDBG_MAX_INTArr(mutualArr, e2.ID) {
			fmt.Printf("[ExtendPath] e2.ID: %v in the mutualArr: %v\n", e2.ID, mutualArr)
			continue
		}
		p2 := e2.PathMat[0]
		eID1 := mutualArr[len(mutualArr)-1]
		e1 := edgesArr[eID1]
		j1 := IndexEID(p2.IDArr, e1.ID)
		j2 := IndexEID(p2.IDArr, e2.ID)
		if j1 < 0 || j2 < 0 {
			continue
		}
		if j2 < j1 {
			p2.IDArr = GetReverseDBG_MAX_INTArr(p2.IDArr)
			j1 = IndexEID(p2.IDArr, e1.ID)
			//j2 = IndexEID(p2.IDArr, e2.ID)
		}
		k1 := IndexEID(p.IDArr, e1.ID)
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
				idx = IndexEID(p.IDArr, e2.ID)
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
	i1 := IndexEID(p.IDArr, mutualArr[0])
	i2 := IndexEID(p.IDArr, mutualArr[len(mutualArr)-1])
	maxP.IDArr = p.IDArr[i1 : i2+1]
	maxP.Freq = p.Freq

	/*// set maxP paths subsequence by edge direction
	ce := edgesArr[maxP.IDArr[0]]
	if ce.StartNID > 0 {
		ea := GetNearEdgeIDArr(nodesArr[ce.StartNID], maxP.IDArr[0])
		if IsInDBG_MAX_INTArr(ea, maxP.IDArr[1]) {
			maxP.IDArr = GetReverseDBG_MAX_INTArr(maxP.IDArr)
		}
	}*/

	edgesArr[mutualArr[0]].SetProcessFlag()
	for _, id := range mutualArr[1:] {
		edgesArr[id].SetDeleteFlag()
		edgesArr[id].SetProcessFlag()
	}
	fmt.Printf("[ExtendPath] mutualArr: %v\n", mutualArr)
	fmt.Printf("[ExtendPath] maxP: %v\n", maxP)

	return
}

/* func ExtendPath(edgesArr []DBGEdge, nodesArr []DBGNode, e DBGEdge) (maxP Path) {
	maxP.IDArr = append(maxP.IDArr, e.PathMat[0].IDArr...)
	maxP.Freq = e.PathMat[0].Freq
	mutualArr := []DBG_MAX_INT{e.ID}
	fmt.Printf("[ExtendPath] edge info, eID: %v\t e.StartNID: %v\te.EndNID: %v\n", e.ID, e.StartNID, e.EndNID)
	if e.StartNID > 0 {
		nd := nodesArr[e.StartNID]
		ea := GetNearEdgeIDArr(nd, e.ID)
		fmt.Printf("[ExtendPath] ea: %v, maxP: %v, nd: %v\n", ea, maxP, nd)
		if len(ea) == 1 {
			maxP.IDArr = ReverseDBG_MAX_INTArr(maxP.IDArr)
			//p1 := IndexEID(maxP.IDArr, ea[0])
			p1 := IndexEID(maxP.IDArr, e.ID) + 1
			if p1 > 0 {
				for j := p1; j < len(maxP.IDArr); j++ {
					ne := edgesArr[maxP.IDArr[j]]
					var nextNID DBG_MAX_INT
					if ne.StartNID == nd.ID {
						nextNID = ne.EndNID
					} else {
						nextNID = ne.StartNID
					}
					if ne.GetDeleteFlag() > 0 || ne.GetUniqueFlag() == 0 || len(ne.PathMat) == 0 {
						nd = nodesArr[nextNID]
						continue
					}
					fmt.Printf("[ExtendPath] ne: %v\nne.PathMat[0]: %v\n", ne, ne.PathMat[0])
					var nP Path // next Path
					if ne.EndNID == nd.ID {
						nP.IDArr = GetReverseDBG_MAX_INTArr(ne.PathMat[0].IDArr)
					} else {
						nP.IDArr = make([]DBG_MAX_INT, len(ne.PathMat[0].IDArr))
						copy(nP.IDArr, ne.PathMat[0].IDArr)
					}
					nP.Freq = ne.PathMat[0].Freq
					fmt.Printf("[ExtendPath] eID1: %v, eID2: %v\n\tmaxP: %v\n\tnP: %v\n\tnd: %v\n", mutualArr[len(mutualArr)-1], ne.ID, maxP, nP, nd)
					if mutualReachable(maxP.IDArr, nP.IDArr, mutualArr[len(mutualArr)-1], ne.ID) {
						var suc bool
						maxP.IDArr, suc = FindConsistenceAndMergePath(maxP.IDArr, nP.IDArr, mutualArr[len(mutualArr)-1], ne.ID)
						if suc == false {
							log.Fatalf("[ExtendPath] Failed to merge two paths, maxP: %v\tnP: %v\n", maxP, nP)
						}
						mutualArr = append(mutualArr, ne.ID)
						if nP.Freq < maxP.Freq {
							maxP.Freq = nP.Freq
						}
						fmt.Printf("after merge, maxP: %v\n", maxP)
					}
					nd = nodesArr[nextNID]
				}
			}
			// reverse Array
			maxP.IDArr = ReverseDBG_MAX_INTArr(maxP.IDArr)
			mutualArr = ReverseDBG_MAX_INTArr(mutualArr)
		}
	}

	fmt.Printf("[ExtendPath] after extend previous edges, maxP: %v\n\te: %v\n", maxP, e)
	if e.EndNID > 0 {
		nd := nodesArr[e.EndNID]
		ea := GetNearEdgeIDArr(nd, e.ID)
		fmt.Printf("[ExtendPath] ea: %v, nd: %v\n", ea, nd)
		if len(ea) == 1 {
			//p1 := IndexEID(maxP.IDArr, ea[0])
			p1 := IndexEID(maxP.IDArr, e.ID) + 1
			if p1 > 0 {
				for j := p1; j < len(maxP.IDArr); j++ {
					ne := edgesArr[maxP.IDArr[j]]
					var nextNID DBG_MAX_INT
					if ne.StartNID == nd.ID {
						nextNID = ne.EndNID
					} else {
						nextNID = ne.StartNID
					}
					fmt.Printf("[ExtendPath] ne: %v\nne.PathMat: %v\n\tnd: %v\n", ne, ne.PathMat, nd)
					if ne.GetDeleteFlag() > 0 || ne.GetUniqueFlag() == 0 || len(ne.PathMat) == 0 {
						nd = nodesArr[nextNID]
						continue
					}
					var nP Path // next Path
					if ne.EndNID == nd.ID {
						nP.IDArr = GetReverseDBG_MAX_INTArr(ne.PathMat[0].IDArr)
					} else {
						nP.IDArr = make([]DBG_MAX_INT, len(ne.PathMat[0].IDArr))
						copy(nP.IDArr, ne.PathMat[0].IDArr)
					}
					nP.Freq = ne.PathMat[0].Freq
					fmt.Printf("[ExtendPath] eID1: %v, eID2: %v\n\tmaxP: %v\n\tnP: %v\n", mutualArr[len(mutualArr)-1], ne.ID, maxP, nP)
					if mutualReachable(maxP.IDArr, nP.IDArr, mutualArr[len(mutualArr)-1], ne.ID) {
						var suc bool
						maxP.IDArr, suc = FindConsistenceAndMergePath(maxP.IDArr, nP.IDArr, mutualArr[len(mutualArr)-1], ne.ID)
						if suc == false {
							log.Fatalf("[ExtendPath] Failed to merge two paths, maxP: %v\tne.PathMat[0]: %v\n", maxP, nP)
						}
						mutualArr = append(mutualArr, ne.ID)
						if nP.Freq < maxP.Freq {
							maxP.Freq = nP.Freq
						}
						fmt.Printf("after merge, maxP: %v\n", maxP)
					}
					nd = nodesArr[nextNID]
				}
			}
		}
	}

	// merge path and process DBG
	if len(mutualArr) == 1 {
		edgesArr[e.ID].SetProcessFlag()
		maxP.IDArr = nil
		maxP.Freq = 0
		//edgesArr[i].PathMat = nil
	} else {
		fmt.Printf("[ExtendPath]mutualArr: %v\t maxP: %v\n", mutualArr, maxP)
		eID1 := mutualArr[0]
		p1 := IndexEID(maxP.IDArr, eID1)
		eID2 := mutualArr[len(mutualArr)-1]
		p2 := IndexEID(maxP.IDArr, eID2)
		if IsTwoEdgesCyclePath(edgesArr, nodesArr, eID1) || IsTwoEdgesCyclePath(edgesArr, nodesArr, eID2) || len(maxP.IDArr[p1:p2+1]) <= 2 {
			maxP.IDArr = nil
			maxP.Freq = 0
			return maxP
		}
		//var np Path
		//np.IDArr = append(np.IDArr, maxP.IDArr[:p1+1]...)
		//np.IDArr = append(np.IDArr, maxP.IDArr[p2+1:]...)
		//np.Freq = maxP.Freq
		var np Path
		np.IDArr = append(np.IDArr, maxP.IDArr...)
		np.Freq = maxP.Freq
		edgesArr[eID1].PathMat = []Path{np}
		CheckPathDirection(edgesArr, nodesArr, eID1)
		maxP.IDArr = maxP.IDArr[p1 : p2+1]
		edgesArr[eID1].SetProcessFlag()
		edgesArr[eID2].SetProcessFlag()
		//edgesArr[eID].PathMat = []Path{maxP}
		//edgesArr[eID].SetProcessFlag()
		for _, id := range mutualArr[1 : len(mutualArr)-1] {
			//edgesArr[id].PathMat = nil
			edgesArr[id].SetDeleteFlag()
			edgesArr[id].SetProcessFlag()
		}
	}

	return maxP
} */

// find maximum path
func findMaxPath(edgesArr []DBGEdge, nodesArr []DBGNode, minMapFreq int, semi bool) (pathArr []Path) {
	for i, e := range edgesArr {
		if i < 2 || e.GetDeleteFlag() > 0 || e.GetProcessFlag() > 0 || len(e.PathMat) != 1 || e.GetTwoEdgesCycleFlag() > 0 {
			continue
		}

		if e.GetUniqueFlag() == 0 {
			continue
		}

		p := e.PathMat[0]
		if p.Freq < minMapFreq || len(p.IDArr) < 2 {
			//continue
			log.Fatalf("[findMaxPath] edgesArr[%v].PathMat: %v\t not contain useful info\n", i, p)
		}

		maxP := ExtendPath(edgesArr, nodesArr, e, minMapFreq, semi)
		if len(maxP.IDArr) > 1 {
			pathArr = append(pathArr, maxP)
		}

	}
	return pathArr
}

func GetLinkNodeID(p0, p1 DBGEdge) DBG_MAX_INT {
	if p0.StartNID == p1.StartNID || p0.StartNID == p1.EndNID {
		return p0.StartNID
	} else if p0.EndNID == p1.StartNID || p0.EndNID == p1.EndNID {
		return p0.EndNID
	} else {
		log.Fatalf("[GetLinkNodeID]not found link node ID p0: %v\n\tp1: %v\n", p0, p1)
	}
	return 0
}

func GetLinkPathArr(nodesArr []DBGNode, edgesArr []DBGEdge, nID DBG_MAX_INT) (p Path) {
	nd := nodesArr[nID]
	inNum, inID := GetEdgeIDComing(nd.EdgeIDIncoming)
	outNum, outID := GetEdgeIDComing(nd.EdgeIDOutcoming)
	if inNum != 1 || outNum != 1 {
		log.Fatalf("[GetLinkPathArr] inNum: %v or outNum: %v != 1\n", inNum, outNum)
	}
	nodesArr[nID].SetProcessFlag()
	arr := [2]DBG_MAX_INT{inID, outID}
	var dpArr [2][]DBG_MAX_INT
	for i, eID := range arr {
		dpArr[i] = append(dpArr[i], eID)
		id := edgesArr[eID].StartNID
		if id == nID {
			id = edgesArr[eID].EndNID
		}
		for id > 0 {
			nd := nodesArr[id]
			inNum, inID := GetEdgeIDComing(nd.EdgeIDIncoming)
			outNum, outID := GetEdgeIDComing(nd.EdgeIDOutcoming)
			if inNum == 1 && outNum == 1 {
				nodesArr[id].SetProcessFlag()
				if inID == eID {
					eID = outID
				} else {
					eID = inID
				}
				if IsTwoEdgesCyclePath(edgesArr, nodesArr, inID) || IsTwoEdgesCyclePath(edgesArr, nodesArr, outID) {
					break
				}
				dpArr[i] = append(dpArr[i], eID)
				if edgesArr[eID].StartNID == id {
					id = edgesArr[eID].EndNID
				} else {
					id = edgesArr[eID].StartNID
				}
			} else {
				break
			}
		}
	}
	p.IDArr = ReverseDBG_MAX_INTArr(dpArr[0])
	p.IDArr = append(p.IDArr, dpArr[1]...)

	return p
}

func CascadePath(p Path, edgesArr []DBGEdge, nodesArr []DBGNode, kmerlen int, changeDBG bool) DBGEdge {

	// keep self cycle path start node OUtcoming, end node Incoming
	e1, e2 := edgesArr[p.IDArr[0]], edgesArr[p.IDArr[1]]
	if e1.EndNID == e2.StartNID || e1.EndNID == e2.EndNID {
		if IsInComing(nodesArr[e1.StartNID].EdgeIDIncoming, p.IDArr[0]) {
			ReverseDBG_MAX_INTArr(p.IDArr)
		}
	} else {
		if IsInComing(nodesArr[e1.EndNID].EdgeIDIncoming, p.IDArr[0]) {
			ReverseDBG_MAX_INTArr(p.IDArr)
		}
	}

	p0 := edgesArr[p.IDArr[0]]
	p1 := edgesArr[p.IDArr[1]]
	//eStartNID, eEndNID := p0.StartNID, p0.EndNID
	var direction uint8
	nID := GetLinkNodeID(p0, p1)
	if p0.StartNID == nID {
		edgesArr[p.IDArr[0]].Utg = GetRCUnitig(edgesArr[p.IDArr[0]].Utg)
		//ReverseCompByteArr(edgesArr[p.IDArr[0]].Utg.Ks)
		//ReverseByteArr(edgesArr[p.IDArr[0]].Utg.Kq)
		edgesArr[p.IDArr[0]].StartNID, edgesArr[p.IDArr[0]].EndNID = edgesArr[p.IDArr[0]].EndNID, edgesArr[p.IDArr[0]].StartNID
		p0 = edgesArr[p.IDArr[0]]
	}
	if nID == p0.EndNID {
		direction = FORWARD
	} else {
		direction = BACKWARD
	}
	lastEID := p0.ID
	strand := true
	for _, eID := range p.IDArr[1:] {
		p0 = edgesArr[p0.ID]
		p1 = edgesArr[eID]
		//fmt.Printf("[CascadePath]p0.ID: %v\n\tp1.ID: %v, lastEID: %v, strand: %v, nID: %v\n", p0.ID, p1.ID, lastEID, strand,nID)
		//inID, outID := p0.ID, p1.ID
		//if direction == BACKWARD { inID, outID = outID, inID }
		//if IsInComing(nodesArr[nID].EdgeIDOutcoming, p1.ID) ==false &&  IsInComing(nodesArr[nID].EdgeIDIncoming, p1.ID) {
		//	inID, outID = outID, inID
		//}
		if direction == BACKWARD {
			if IsInComing(nodesArr[nID].EdgeIDOutcoming, lastEID) {
				if IsInComing(nodesArr[nID].EdgeIDIncoming, p1.ID) {
					if edgesArr[lastEID].StartNID == nID {
						if p1.EndNID != nID {
							strand = !strand
						}
					} else {
						if p1.StartNID != nID {
							strand = !strand
						}
					}
				} else {
					log.Fatalf("[CascadePath] nodesArr[%v]: %v\n\t, nID set error\n", nID, nodesArr[nID])
				}
			} else {
				if IsInComing(nodesArr[nID].EdgeIDOutcoming, p1.ID) {
					if edgesArr[lastEID].EndNID == nID {
						if p1.StartNID != nID {
							strand = !strand
						}
					} else {
						if p1.EndNID != nID {
							strand = !strand
						}
					}
				} else {
					log.Fatalf("[CascadePath]BACKWARD nodesArr[%v]: %v\n\t, nID set error\n", nID, nodesArr[nID])
				}
			}
			if !strand {
				p1.Utg = GetRCUnitig(p1.Utg)
			}
			fmt.Printf("[CascadePath]p0.ID: %v, p1.ID: %v, lastEID: %v, strand: %v, nID: Incoming: %v, Outcoming: %v\n", p0.ID, p1.ID, lastEID, strand, nodesArr[nID].EdgeIDIncoming, nodesArr[nID].EdgeIDOutcoming)
			edgesArr[p0.ID].Utg = ConcatEdges(p1.Utg, p0.Utg, kmerlen)
		} else {
			if IsInComing(nodesArr[nID].EdgeIDIncoming, lastEID) {
				if IsInComing(nodesArr[nID].EdgeIDOutcoming, p1.ID) {
					if edgesArr[lastEID].EndNID == nID {
						if p1.StartNID != nID {
							strand = !strand
						}
					} else {
						if p1.EndNID != nID {
							strand = !strand
						}
					}
				} else {
					log.Fatalf("[CascadePath] nodesArr[%v]: %v, nID set error\n", nID, nodesArr[nID])
				}
			} else {
				if IsInComing(nodesArr[nID].EdgeIDIncoming, p1.ID) {
					if edgesArr[lastEID].StartNID == nID {
						if p1.EndNID != nID {
							strand = !strand
						}
					} else {
						if p1.StartNID != nID {
							strand = !strand
						}
					}
				} else {
					log.Fatalf("[CascadePath] nodesArr[%v]: %v\n\t, nID set error\n", nID, nodesArr[nID])
				}
			}
			if !strand {
				p1.Utg = GetRCUnitig(p1.Utg)
			}
			fmt.Printf("[CascadePath]FORWARD p0.ID: %v, p1.ID: %v, lastEID: %v, strand: %v, nID: Incoming: %v, Outcoming: %v\n", p0.ID, p1.ID, lastEID, strand, nodesArr[nID].EdgeIDIncoming, nodesArr[nID].EdgeIDOutcoming)
			edgesArr[p0.ID].Utg = ConcatEdges(p0.Utg, p1.Utg, kmerlen)
		}

		if nID == edgesArr[p1.ID].StartNID {
			nID = edgesArr[p1.ID].EndNID
		} else {
			nID = edgesArr[p1.ID].StartNID
		}
		lastEID = p1.ID
	}

	if !changeDBG {
		return edgesArr[p0.ID]
	}

	if direction == FORWARD {
		if !SubstituteEdgeID(nodesArr, p0.EndNID, p0.ID, 0) {
			log.Fatalf("[ReconstructDBG] SubsitututeEdgeID for Merged edge error, node: %v\n\tedge:%v\n", nodesArr[p0.EndNID], p0)
		}
		if !SubstituteEdgeID(nodesArr, nID, p1.ID, p0.ID) {
			log.Fatalf("[ReconstructDBG] SubsitututeEdgeID for Merged edge error, node: %v\n\tedge:%v\n", nodesArr[nID], p1)
		}
		edgesArr[p0.ID].EndNID = nID
	} else {
		if !SubstituteEdgeID(nodesArr, p0.StartNID, p0.ID, 0) {
			log.Fatalf("[ReconstructDBG] SubsitututeEdgeID for Merged edge error, node: %v\n\tedge:%v\n", nodesArr[p0.StartNID], p0)
		}
		if !SubstituteEdgeID(nodesArr, nID, p1.ID, p0.ID) {
			log.Fatalf("[ReconstructDBG] SubsitututeEdgeID for Merged edge error, node: %v\n\tedge:%v\n", nodesArr[nID], p1)
		}
		edgesArr[p0.ID].StartNID = nID
	}
	return edgesArr[p0.ID]
}

// reset the nodesArr EdgeIDIncoming and EdgeIDOutcoming
func ResetDeleteEdgeIncoming(edgesArr []DBGEdge, nodesArr []DBGNode) {
	for _, e := range edgesArr {
		if e.ID < 2 || e.GetDeleteFlag() == 0 {
			continue
		}
		fmt.Printf("[ResetDeleteEdgeIncoming] deleted edge:%v\n", e)
		if e.StartNID > 0 {
			SubstituteEdgeID(nodesArr, e.StartNID, e.ID, 0)
			//log.Fatalf("[ResetDeleteEdgeIncoming] SubsitututeEdgeID for deleted edge error, nodesArr[%v]: %v\n\tedge:%v\n\tedge len: %v\n", e.StartNID, nodesArr[e.StartNID], e, e.GetSeqLen())
		}
		if e.EndNID > 0 {
			SubstituteEdgeID(nodesArr, e.EndNID, e.ID, 0)
			//log.Fatalf("[ResetDeleteEdgeIncoming] SubsitututeEdgeID for deleted edge error, nodesArr[%v]: %v\n\tedge:%v\n\tedge len: %v\n", e.EndNID, nodesArr[e.EndNID], e, e.GetSeqLen())
		}
	}
}

/*func ResetMergedEdgeNID(edgesArr []DBGEdge, nodesArr []DBGNode, joinPathArr []Path, IDMapPath map[DBG_MAX_INT]uint32, kmerlen int) {
	for _, e := range edgesArr {
		if e.ID == 0 {
			continue
		}
		if idx, ok := IDMapPath[e.ID]; ok {
			if e.ID == joinPathArr[idx].IDArr[0] {
				if e.StartNID > 0 {
					//SubstituteEdgeID(nodesArr, e.StartNID, e.ID, 0)
					if !IsInDBGNode(nodesArr[e.StartNID], e.ID) {
						//var sBnt constructcf.KmerBnt
						//sBnt.Seq = e.Utg.Ks[:kmerlen-1]
						ks := constructcf.GetReadBntKmer(e.Utg.Ks[:kmerlen-1], 0, kmerlen-1)
						ts := ks
						rs := constructcf.ReverseComplet(ks)
						if ks.BiggerThan(rs) {
							ks, rs = rs, ks
						}
						if reflect.DeepEqual(ks.Seq, nodesArr[e.StartNID].Seq) == false {
							log.Fatalf("[ResetMergedEdgeNID] ks != nodesArr[e.StartNID].Seq, ks: %v\n\tnodesArr[%v]: %v\n", ks, e.StartNID, nodesArr[e.StartNID])
						}
						c := e.Utg.Ks[kmerlen-1]
						if reflect.DeepEqual(ts, ks) {
							nodesArr[e.StartNID].EdgeIDOutcoming[c] = e.ID
						} else {
							nodesArr[e.StartNID].EdgeIDIncoming[bnt.BntRev[c]] = e.ID
						}
					}
				}
				if e.EndNID > 0 {
					//SubstituteEdgeID(nodesArr, e.EndNID, e.ID, 0)
					if !IsInDBGNode(nodesArr[e.EndNID], e.ID) {
						//var sBnt constructcf.ReadBnt
						//sBnt.Seq = e.Utg.Ks[e.GetSeqLen()-kmerlen+1:]
						ks := constructcf.GetReadBntKmer(e.Utg.Ks[e.GetSeqLen()-kmerlen+1:], 0, kmerlen-1)
						ts := ks
						rs := constructcf.ReverseComplet(ks)
						if ks.BiggerThan(rs) {
							ks, rs = rs, ks
						}
						if reflect.DeepEqual(ks.Seq, nodesArr[e.EndNID].Seq) == false {
							log.Fatalf("[ResetMergedEdgeNID] ks != nodesArr[e.StartNID].Seq, ks: %v\n\tnodesArr[%v]: %v\n", ks, e.StartNID, nodesArr[e.StartNID])
						}
						c := e.Utg.Ks[e.GetSeqLen()-kmerlen]
						if reflect.DeepEqual(ks, ts) {
							nodesArr[e.EndNID].EdgeIDIncoming[c] = e.ID
						} else {
							nodesArr[e.EndNID].EdgeIDOutcoming[bnt.BntRev[c]] = e.ID
						}
					}
				}
			}
		}
	}
}*/

func ResetNodeEdgecoming(edgesArr []DBGEdge, nodesArr []DBGNode) {
	for i, nd := range nodesArr {
		if nd.ID == 0 {
			continue
		}
		if nd.GetDeleteFlag() > 0 {
			nodesArr[i].ID = 0
			nodesArr[i].EdgeIDIncoming = [4]DBG_MAX_INT{0, 0, 0, 0}
			nodesArr[i].EdgeIDOutcoming = [4]DBG_MAX_INT{0, 0, 0, 0}
			nodesArr[i].Seq = nil
			continue
		}
		for j := 0; j < bnt.BaseTypeNum; j++ {
			if nd.EdgeIDIncoming[j] > 0 {
				eID := nd.EdgeIDIncoming[j]
				if edgesArr[eID].StartNID != nd.ID && edgesArr[eID].EndNID != nd.ID {
					nodesArr[i].EdgeIDIncoming[j] = 0
				}
			}
			if nd.EdgeIDOutcoming[j] > 0 {
				eID := nd.EdgeIDOutcoming[j]
				if edgesArr[eID].StartNID != nd.ID && edgesArr[eID].EndNID != nd.ID {
					nodesArr[i].EdgeIDOutcoming[j] = 0
				}
			}
		}
	}
}

func InitialDBGNodeProcessFlag(nodesArr []DBGNode) {
	for i, _ := range nodesArr {
		nodesArr[i].ResetProcessFlag()
	}
}

func ResetIDMapPath(IDMapPath map[DBG_MAX_INT]uint32, joinPathArr []Path, path Path) []Path {
	var np Path
	np.Freq = path.Freq
	joinIdx := -1
	for _, t := range path.IDArr {
		if idx, ok := IDMapPath[t]; ok {
			arr := joinPathArr[idx].IDArr
			if t == arr[len(arr)-1] {
				ReverseDBG_MAX_INTArr(arr)
			}
			np.IDArr = append(np.IDArr, arr...)
			delete(IDMapPath, t)
			if t == arr[0] {
				t = arr[len(arr)-1]
			} else {
				t = arr[0]
			}
			delete(IDMapPath, t)
			joinIdx = int(idx)
		} else {
			np.IDArr = append(np.IDArr, t)
		}
	}
	if joinIdx >= 0 {
		joinPathArr[joinIdx] = np
	} else {
		joinPathArr = append(joinPathArr, np)
		joinIdx = len(joinPathArr) - 1
	}
	IDMapPath[np.IDArr[0]] = uint32(joinIdx)
	IDMapPath[np.IDArr[len(np.IDArr)-1]] = uint32(joinIdx)

	return joinPathArr
}

func CleanEdgesArr(edgesArr []DBGEdge, nodesArr []DBGNode) (deleteNodeNum, deleteEdgeNum int) {
	for i, e := range edgesArr {
		if e.ID == 0 || e.GetDeleteFlag() > 0 || e.StartNID == 0 || e.EndNID == 0 {
			continue
		}
		var num int
		if IsInComing(nodesArr[e.StartNID].EdgeIDIncoming, e.ID) {
			arr := GetNearEdgeIDArr(nodesArr[e.StartNID], e.ID, true)
			num = len(arr)
		} else {
			arr := GetNearEdgeIDArr(nodesArr[e.StartNID], e.ID, false)
			num = len(arr)
		}
		if IsInComing(nodesArr[e.EndNID].EdgeIDIncoming, e.ID) {
			arr := GetNearEdgeIDArr(nodesArr[e.EndNID], e.ID, true)
			num += len(arr)
		} else {
			arr := GetNearEdgeIDArr(nodesArr[e.EndNID], e.ID, false)
			num += len(arr)
		}
		if num == 0 {
			deleteEdgeNum++
			edgesArr[i].SetDeleteFlag()
			fmt.Printf("[CleanEdgesArr]delete edge len: %v\t%v\n", len(edgesArr[i].Utg.Ks), edgesArr[i])
			if e.StartNID > 0 {
				nodesArr[e.StartNID].SetDeleteFlag()
				deleteNodeNum++
			}
			if e.EndNID > 0 {
				nodesArr[e.EndNID].SetDeleteFlag()
				deleteNodeNum++
			}
		}
	}

	return deleteNodeNum, deleteEdgeNum
}

func ReconstructDBG(edgesArr []DBGEdge, nodesArr []DBGNode, joinPathArr []Path, kmerlen int) {
	// join path that unique path
	for _, p := range joinPathArr {
		if len(p.IDArr) <= 1 {
			continue
		}
		fmt.Printf("[ReconstructDBG] p: %v\n", p)
		//if IsTwoEdgeCyclePath(path) { joinPathArr[i].IDArr = }

		CascadePath(p, edgesArr, nodesArr, kmerlen, true)
	}
	// Reset delete edges coming DBGNodes
	ResetDeleteEdgeIncoming(edgesArr, nodesArr)

	/*
		fmt.Printf("[ReconstructDBG] finished joinPathArr reconstruct\n")
		fmt.Printf("[ReconstructDBG] edgesArr[%v]: %v\n", 481, edgesArr[481])

		// reset the nodesArr EdgeIDIncoming and EdgeIDOutcoming
		ResetDeleteEdgeIncoming(edgesArr, nodesArr)
		fmt.Printf("[ReconstructDBG] after ResetDeleteEdgeIncoming edgesArr[%v]: %v\n", 481, edgesArr[481])
		ResetMergedEdgeNID(edgesArr, nodesArr, joinPathArr, IDMapPath)
		fmt.Printf("[ReconstructDBG] after ResetMergedEdgeNID edgesArr[%v]: %v\n", 481, edgesArr[481])

		// simplify DBG
		//var deleteNodeNum, deleteEdgeNum int
		// clean edgesArr
		deleteNodeNum, deleteEdgeNum := CleanEdgesArr(edgesArr, nodesArr)
		fmt.Printf("[ReconstructDBG] after CleanEdgesArr edgesArr[%v]: %v\n", 481, edgesArr[481])
		// clean nodesArr
		// initial nodesArr procesed Flag
		InitialDBGNodeProcessFlag(nodesArr)

		for i, nd := range nodesArr {
			if nd.ID == 0 || nd.GetDeleteFlag() > 0 || nd.GetProcessFlag() > 0 {
				continue
			}
			inNum, inID := GetEdgeIDComing(nd.EdgeIDIncoming)
			outNum, outID := GetEdgeIDComing(nd.EdgeIDOutcoming)
			if inNum == 0 && outNum == 0 {
				nodesArr[i].SetDeleteFlag()
				nodesArr[i].ID = 0
				deleteNodeNum++
			} else if inNum+outNum == 1 {
				id := inID
				if outNum == 1 {
					id = outID
				}
				if edgesArr[id].StartNID == nd.ID {
					edgesArr[id].StartNID = 0
				} else {
					edgesArr[id].EndNID = 0
				}
				deleteNodeNum++
			} else if inNum == 1 && outNum == 1 && inID != outID { // prevent cycle ring
				if edgesArr[inID].StartNID != nd.ID && edgesArr[inID].EndNID != nd.ID {
					edgesArr[inID].SetDeleteFlag()
					fmt.Printf("[ReconstructDBG] clean nodesArr deleted edge, edgesArr[%v]: %v\n\tnd: %v\n", inID, edgesArr[inID], nd)
					fmt.Printf("[ReconstructDBG] clean nodesArr deleted edge, edgesArr[%v]: %v\n", outID, edgesArr[outID])
					continue
				}
				if edgesArr[outID].StartNID != nd.ID && edgesArr[outID].EndNID != nd.ID {
					edgesArr[outID].SetDeleteFlag()
					fmt.Printf("[ReconstructDBG] clean nodesArr deleted edge, edgesArr[%v]: %v\n\tnd: %v\n", inID, edgesArr[inID], nd)
					fmt.Printf("[ReconstructDBG] clean nodesArr deleted edge, edgesArr[%v]: %v\n", outID, edgesArr[outID])
					continue
				}
				path := GetLinkPathArr(nodesArr, edgesArr, nd.ID)
				joinPathArr = ResetIDMapPath(IDMapPath, joinPathArr, path)
				fmt.Printf("[ReconstructDBG]clean nodesArr path: %v\n", path)

				for _, eID := range path.IDArr[1:] {
					edgesArr[eID].SetDeleteFlag()
				}
				// join path
				CascadePath(path, edgesArr, nodesArr)
			}
		}

		// reset the nodesArr EdgeIDIncoming and EdgeIDOutcoming
		ResetDeleteEdgeIncoming(edgesArr, nodesArr)
		ResetMergedEdgeNID(edgesArr, nodesArr, joinPathArr, IDMapPath)
		ResetNodeEdgecoming(edgesArr, nodesArr)

		fmt.Printf("[ReconstructDBG] delete edge number is : %v\n\tdelete node number is : %v\n\tmerge edge path number is : %v\n", deleteEdgeNum, deleteNodeNum, len(joinPathArr))
	*/
}

func findPathOverlap(jp Path, pathArr []DBG_MAX_INT, edgesArr []DBGEdge) (id DBG_MAX_INT, num int) {
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
		rp := GetReverseDBG_MAX_INTArr(jpArr)
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

/*func CheckPathDirection(edgesArr []DBGEdge, nodesArr []DBGNode, eID DBG_MAX_INT) {
	e := edgesArr[eID]
	step := 1
	var ea1, ea2 []DBG_MAX_INT
	if e.StartNID > 0 {
		if IsInComing(nodesArr[e.StartNID].EdgeIDIncoming, eID) {
			ea1 = GetNearEdgeIDArr(nodesArr[e.StartNID], eID, true)
		} else {
			ea1 = GetNearEdgeIDArr(nodesArr[e.StartNID], eID, false)
		}
	}
	if e.EndNID > 0 {
		if IsInComing(nodesArr[e.EndNID].EdgeIDIncoming, eID) {
			ea2 = GetNearEdgeIDArr(nodesArr[e.EndNID], eID, true)
		} else {
			ea2 = GetNearEdgeIDArr(nodesArr[e.EndNID], eID, false)
		}
	}

	// found two edge cycle
	if e.GetTwoEdgesCycleFlag() > 0 {
		ea1 = GetNearEdgeIDArr(nodesArr[e.EndNID], ea1[0])
		ea2 = GetNearEdgeIDArr(nodesArr[e.StartNID], ea2[0])
		if len(ea1) == 2 && len(ea2) == 2 {
			if ea1[0] == e.ID {
				ea1[0] = ea1[1]
			}
			ea1 = ea1[:1]
			if ea2[0] == e.ID {
				ea2[0] = ea2[1]
			}
			ea2 = ea2[:1]
			step = 2
		} else {
			edgesArr[eID].ResetUniqueFlag()
			edgesArr[eID].ResetSemiUniqueFlag()
			return
		}
	} else if IsIntersection(ea1, ea2) {
		edgesArr[eID].ResetUniqueFlag()
		edgesArr[eID].ResetSemiUniqueFlag()
		return
	}

	for i, p := range e.PathMat {
		idx := IndexEID(p.IDArr, eID)
		if idx < len(p.IDArr)-step {
			if IsInDBG_MAX_INTArr(ea1, p.IDArr[idx+step]) {
				ReverseDBG_MAX_INTArr(edgesArr[eID].PathMat[i].IDArr)
			}
		} else {
			if idx-step >= 0 && IsInDBG_MAX_INTArr(ea2, p.IDArr[idx-step]) {
				ReverseDBG_MAX_INTArr(edgesArr[eID].PathMat[i].IDArr)
			}
		}
	}

	//fmt.Printf("[CheckPathDirection] ea1: %v, ea2: %v, PathMat: %v\n", ea1, ea2, edgesArr[eID].PathMat)

}*/

func AdjustEIDPath(e DBGEdge, joinPathArr []Path, IDMapPath map[DBG_MAX_INT]uint32) (p Path) {
	ep := e.PathMat[0]
	idx := IndexEID(ep.IDArr, e.ID)
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
			a2 = GetReverseDBG_MAX_INTArr(a2)
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
			a2 := GetReverseDBG_MAX_INTArr(a2)
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
		/*else {
			log.Fatalf("[AdjustEIDPath] not found overlap Path, eID: %v, ep: %v\n\tjoinPathArr: %v\n", e.ID, a1, a2)
		}*/
	}
	if len(p.IDArr) == 0 {
		p.IDArr = ep.IDArr
	}
	fmt.Printf("[AdjustEIDPath] p: %v\n", p)
	return p
}

/*func AdjustPathMat(edgesArr []DBGEdge, nodesArr []DBGNode, joinPathArr []Path, IDMapPath map[DBG_MAX_INT]uint32) {

	// adjust pathMat
	for i, e := range edgesArr {
		if e.ID == 0 || e.GetDeleteFlag() > 0 || len(e.PathMat) == 0 {
			continue
		}
		if IsTwoEdgesCyclePath(edgesArr, nodesArr, e.ID) {
			continue
		}
		edgesArr[i].PathMat[0] = AdjustEIDPath(edgesArr[i], joinPathArr, IDMapPath)
		CheckPathDirection(edgesArr, nodesArr, e.ID)
		e = edgesArr[i]
		p := e.PathMat[0]
		fmt.Printf("[AdjustPathMat] e.PathMat: %v\n", e.PathMat)
		idx := IndexEID(p.IDArr, e.ID)
		arr := []DBG_MAX_INT{e.ID}
		if e.StartNID > 0 && idx > 0 {
			nd := nodesArr[e.StartNID]
			te := e
			rp := GetReverseDBG_MAX_INTArr(p.IDArr[:idx])
			for j := 0; j < len(rp); {
				eID := rp[j]
				ea := GetNearEdgeIDArr(nd, te.ID)
				if idx, ok := IDMapPath[eID]; ok {
					id, num := findPathOverlap(joinPathArr[idx], rp[j:], edgesArr)
					if id == 0 {
						log.Fatalf("[AdjustPathMat] not found overlap Path, eID: %v, p: %v\n\trp: %v\n", eID, joinPathArr[idx], rp)
					} else {
						arr = append(arr, id)
						j += num
					}
				} else {
					if IsInDBG_MAX_INTArr(ea, eID) {
						arr = append(arr, eID)
						j++
					} else {
						log.Fatalf("[AdjustPathMat]not found joinPathArr, edge ID: %v, PathMat: %v, eID: %v\n", e.ID, p, eID)
					}
				}
				te = edgesArr[arr[len(arr)-1]]
				if te.StartNID == nd.ID {
					nd = nodesArr[te.EndNID]
				} else {
					nd = nodesArr[te.StartNID]
				}
			}
			ReverseDBG_MAX_INTArr(arr)
		}

		if e.EndNID > 0 && idx < len(p.IDArr)-1 {
			nd := nodesArr[e.EndNID]
			te := e
			tp := p.IDArr[idx+1:]
			for j := 0; j < len(tp); {
				eID := tp[j]
				ea := GetNearEdgeIDArr(nd, te.ID)
				if idx, ok := IDMapPath[eID]; ok {
					id, num := findPathOverlap(joinPathArr[idx], tp[j:], edgesArr)
					if id == 0 {
						log.Fatalf("[AdjustPathMat] not found overlap Path, eID: %v, p: %v\n\ttp: %v\n", eID, joinPathArr[idx], tp)
					} else {
						arr = append(arr, id)
						j += num
					}
				} else {
					if IsInDBG_MAX_INTArr(ea, eID) {
						arr = append(arr, eID)
						j++
					} else {
						log.Fatalf("[AdjustPathMat] not found joinPathArr and nd: %v, edge ID: %v, PathMat: %v, eID: %v\n", nd, e.ID, p, eID)
					}
				}
				te = edgesArr[arr[len(arr)-1]]
				if te.StartNID == nd.ID {
					nd = nodesArr[te.EndNID]
				} else {
					nd = nodesArr[te.StartNID]
				}
			}
		}
		var np Path
		np.Freq = p.Freq
		np.IDArr = arr
		edgesArr[i].PathMat[0] = np
		fmt.Printf("[AdjustPathMat] before adjust path: %v\n\tafter path: %v\n", p, np)
	}
}*/

func IsContainBoundaryArr(a1, a2 []DBG_MAX_INT) bool {
	if reflect.DeepEqual(a1[:len(a2)], a2) {
		return true
	}
	if reflect.DeepEqual(a1[len(a1)-len(a2):], a2) {
		return true
	}
	ra2 := GetReverseDBG_MAX_INTArr(a2)
	if reflect.DeepEqual(a1[:len(ra2)], ra2) {
		return true
	}
	if reflect.DeepEqual(a1[len(a1)-len(ra2):], ra2) {
		return true
	}
	return false
}

// constuct map edge ID to the path
func ConstructIDMapPath(joinPathArr []Path) map[DBG_MAX_INT]uint32 {
	IDMapPath := make(map[DBG_MAX_INT]uint32)
	for i, p := range joinPathArr {
		sc := [2]DBG_MAX_INT{p.IDArr[0], p.IDArr[len(p.IDArr)-1]}
		for _, id := range sc {
			if idx, ok := IDMapPath[id]; ok {
				a := joinPathArr[idx]
				log.Fatalf("[ConstructIDMapPath] path: %v collison with : %v\n", a, p)
			} else {
				IDMapPath[id] = uint32(i)
			}
			/*fmt.Printf("[ConstructIDMapPath] path: %v collison with : %v\n", a, p)
			a1, a2 := p, a
			if len(a.IDArr) > len(p.IDArr) {
				a1, a2 = a, p
			}
			if IsContainBoundaryArr(a1.IDArr, a2.IDArr) {
				if len(a1.IDArr) > len(a.IDArr) {
					fmt.Printf("[ConstructIDMapPath] delete : %v\n", joinPathArr[idx])
					delete(IDMapPath, joinPathArr[idx].IDArr[0])
					delete(IDMapPath, joinPathArr[idx].IDArr[len(joinPathArr[idx].IDArr)-1])
					joinPathArr[idx].IDArr = nil
					joinPathArr[idx].Freq = 0
					IDMapPath[id] = uint32(i)
				} else {
					joinPathArr[i].IDArr = nil
					joinPathArr[i].Freq = 0
					if j == 0 {
						break
					} else {
						delete(IDMapPath, sc[0])
					}
				}
			} else {
				log.Fatalf("[ConstructIDMapPath] path: %v collison with : %v\n", a, p)
			} */

		}
	}

	return IDMapPath
}

func DeleteJoinPathArrEnd(edgesArr []DBGEdge, joinPathArr []Path) {
	for _, p := range joinPathArr {
		if p.Freq > 0 && len(p.IDArr) > 0 {
			eID := p.IDArr[len(p.IDArr)-1]
			edgesArr[eID].SetDeleteFlag()
		}
	}
}

/*func SimplifyByNGS(opt Options, nodesArr []DBGNode, edgesArr []DBGEdge, mapNGSFn string) {
	// add to the DBGEdge pathMat
	AddPathToDBGEdge(edgesArr, mapNGSFn)

	// merge pathMat
	MergePathMat(edgesArr, nodesArr, opt.MinMapFreq)

	// find maximum path
	joinPathArr := findMaxPath(edgesArr, nodesArr, opt.MinMapFreq, false) // just used for unique edge
	joinPathArr1 := findMaxPath(edgesArr, nodesArr, opt.MinMapFreq, true) // just used for semi-unique edge
	fmt.Printf("[SimplifyByNGS] findMaxPath number of the uinque edges : %v\n", len(joinPathArr))
	fmt.Printf("[SimplifyByNGS] findMaxPath number of  the semi-uinque edges : %v\n", len(joinPathArr1))
	// ReconstructDBG must first reconstruct uinque edges path , then process semi-unique edges,
	// because semi-unique path maybe been contained in the uinque path,
	// that will cause collison when  ReconstructDBG()
	joinPathArr = append(joinPathArr, joinPathArr1...)
	i := 125
	if edgesArr[i].GetDeleteFlag() == 0 {
		log.Fatalf("[SimplifyByNGS]edgesArr[%v]: %v\n", i, edgesArr[i])
	}
	// constuct map edge ID to the path
	//DeleteJoinPathArrEnd(edgesArr, joinPathArr)
	ReconstructDBG(edgesArr, nodesArr, joinPathArr, opt.Kmer)
	// debug code
	//graphfn := opt.Prefix + ".afterNGS.dot"
	//GraphvizDBGArr(nodesArr, edgesArr, graphfn)

	//AdjustPathMat(edgesArr, nodesArr, joinPathArr, IDMapPath)
}*/

func CheckInterConnectivity(edgesArr []DBGEdge, nodesArr []DBGNode) {
	for i := range edgesArr {
		e := &edgesArr[i]
		if e.ID < 2 || e.GetDeleteFlag() > 0 {
			continue
		}
		if e.StartNID > 0 && e.StartNID == e.EndNID && nodesArr[e.StartNID].GetDeleteFlag() > 0 {
			if e.GetSeqLen() < 1000 {
				e.SetDeleteFlag()
			} else {
				e.StartNID = DBG_MAX_INT(0)
				e.EndNID = DBG_MAX_INT(0)
			}
			continue
		}
		if e.StartNID > 0 && (nodesArr[e.StartNID].GetDeleteFlag() > 0 || !IsInDBGNode(nodesArr[e.StartNID], e.ID)) {
			log.Fatalf("[CheckInterConnectivity]edge check nd: %v\n\tedge: %v\n", nodesArr[e.StartNID], e)
		}
		if e.EndNID > 0 && (nodesArr[e.EndNID].GetDeleteFlag() > 0 || !IsInDBGNode(nodesArr[e.EndNID], e.ID)) {
			log.Fatalf("[CheckInterConnectivity]edge check nd: %v\n\tedge: %v\n", nodesArr[e.EndNID], e)
		}
	}

	for _, nd := range nodesArr {
		if nd.ID == 0 || nd.GetDeleteFlag() > 0 {
			continue
		}

		for j := 0; j < bnt.BaseTypeNum; j++ {
			if nd.EdgeIDIncoming[j] > 1 {
				eID := nd.EdgeIDIncoming[j]
				if edgesArr[eID].GetDeleteFlag() > 0 || (edgesArr[eID].StartNID != nd.ID && edgesArr[eID].EndNID != nd.ID) {
					log.Fatalf("[CheckInterConnectivity]node check nd: %v\n\tedge: %v\n", nd, edgesArr[eID])
				}
			}
			if nd.EdgeIDOutcoming[j] > 1 {
				eID := nd.EdgeIDOutcoming[j]
				if edgesArr[eID].GetDeleteFlag() > 0 || (edgesArr[eID].StartNID != nd.ID && edgesArr[eID].EndNID != nd.ID) {
					log.Fatalf("[CheckInterConnectivity]node check nd: %v\n\tedge: %v\n", nd, edgesArr[eID])
				}
			}
		}
	}
}

type Options struct {
	utils.ArgsOpt
	TipMaxLen     int
	WinSize       int
	MaxNGSReadLen int
	MinMapFreq    int
	Correct       bool
	Graph         bool
	//MaxMapEdgeLen int // max length of edge that don't need cut two flank sequence to map Long Reads
}

func checkArgs(c cli.Command) (opt Options, succ bool) {

	var ok bool
	opt.TipMaxLen, ok = c.Flag("TipMaxLen").Get().(int)
	if !ok {
		log.Fatalf("[checkArgs] argument 'TipMaxLen': %v set error\n ", c.Flag("TipMaxlen").String())
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
	if opt.MinMapFreq < 5 && opt.MinMapFreq >= 20 {
		log.Fatalf("[checkArgs] argument 'MinMapFreq': %v must 5 <= MinMapFreq < 20\n", c.Flag("MinMapFreq").String())
	}

	opt.Correct, ok = c.Flag("Correct").Get().(bool)
	if !ok {
		log.Fatalf("[checkArgs] argument 'Correct': %v set error\n ", c.Flag("Correct").String())
	}
	opt.Graph, ok = c.Flag("Graph").Get().(bool)
	if !ok {
		log.Fatalf("[checkArgs] argument 'Graph': %v set error\n ", c.Flag("Graph").String())
	}

	if opt.TipMaxLen == 0 {
		opt.TipMaxLen = opt.MaxNGSReadLen
	}

	/*opt.MaxMapEdgeLen, ok = c.Flag("MaxMapEdgeLen").Get().(int)
	if !ok {
		log.Fatalf("[checkArgs] argument 'MaxMapEdgeLen': %v set error\n ", c.Flag("MaxMapEdgeLen").String())
	}
	if opt.MaxMapEdgeLen < 2000 {
		log.Fatalf("[checkArgs] argument 'MaxMapEdgeLen': %v must bigger than 2000\n", c.Flag("MaxMapEdgeLen").String())
	}*/

	succ = true
	return opt, succ
}

func Smfy(c cli.Command) {

	t0 := time.Now()
	// check agruments
	gOpt, suc := utils.CheckGlobalArgs(c.Parent())
	if suc == false {
		log.Fatalf("[Smfy] check global Arguments error, opt: %v\n", gOpt)
	}
	opt := Options{gOpt, 0, 0, 0, 0, false, false}
	tmp, suc := checkArgs(c)
	if suc == false {
		log.Fatalf("[Smfy] check Arguments error, opt: %v\n", tmp)
	}
	opt.TipMaxLen = tmp.TipMaxLen
	opt.MaxNGSReadLen = tmp.MaxNGSReadLen
	opt.WinSize = tmp.WinSize
	opt.MinMapFreq = tmp.MinMapFreq
	opt.Correct = tmp.Correct
	opt.Graph = tmp.Graph
	//opt.MaxMapEdgeLen = tmp.MaxMapEdgeLen
	fmt.Printf("Arguments: %v\n", opt)

	// set package-level variable
	//Kmerlen = opt.Kmer

	DBGStatfn := opt.Prefix + ".DBG.stat"
	nodesSize, edgesSize := DBGStatReader(DBGStatfn)
	//nodesSize := len(nodeMap)
	fmt.Printf("[Smfy] len(nodesArr): %v, length of edge array: %v\n", nodesSize, edgesSize)
	nodesfn := opt.Prefix + ".nodes.Arr.br"
	nodesArr := make([]DBGNode, nodesSize)
	fc := make(chan int, 1)
	go NodesArrReader(nodesfn, nodesArr, opt.Kmer, fc)
	// read edges file
	edgesbrfn := opt.Prefix + ".edges.fq.br"
	edgesArr := ReadEdgesFromFile(edgesbrfn, edgesSize)
	<-fc
	//gfn1 := opt.Prefix + ".beforeSmfyDBG.dot"
	//GraphvizDBGArr(nodesArr, edgesArr, gfn1)

	//gfn := opt.Prefix + ".smfyDBG.dot"
	//GraphvizDBG(nodeMap, edgesArr, gfn)
	// reconstruct consistence De Bruijn Graph
	// ReconstructConsistenceDBG(nodeMap, edgesArr)
	for i := 0; i < 5; i++ {
		fmt.Printf("[Smfy] %v cycle smfy DBG....\n", i)
		SmfyDBG(nodesArr, edgesArr, opt, false)
	}
	edgesArr = SetEdgeID(nodesArr, edgesArr)
	nodesArr = SetNodeID(nodesArr, edgesArr)
	CheckDBGSelfCycle(nodesArr, edgesArr, opt.Kmer)
	CheckInterConnectivity(edgesArr, nodesArr)
	//MakeSelfCycleEdgeOutcomingToIncoming(nodesArr, edgesArr, opt)
	// set the unique edge of edgesArr
	uniqueNum, bubbleRepeatNum, twoEdgeCycleNum, selfCycleNum, selfCycleSameComingNum, bubbleEdgeNum := SetDBGEdgesUniqueFlag(edgesArr, nodesArr, opt.Kmer)
	fmt.Printf("[Smfy] the number of DBG Unique  Edges is : %d\n", uniqueNum)
	fmt.Printf("[Smfy] the number of DBG bubble repeat  Edges is : %d\n", bubbleRepeatNum)
	fmt.Printf("[Smfy] the number of DBG twoEdgeCycleNum  Edges is : %d\n", twoEdgeCycleNum)
	fmt.Printf("[Smfy] the number of DBG selfCycleNum  Edges is : %d\n", selfCycleNum)
	fmt.Printf("[Smfy] the number of DBG selfCycle and have same EdgeIDComing  Edges is : %d\n", selfCycleSameComingNum)
	fmt.Printf("[Smfy] the number of DBG Bubble  Edges is : %d\n", bubbleEdgeNum)

	// map Illumina reads to the DBG and find reads map path for simplify DBG
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
	fmt.Printf("[Smfy]smfy total used : %v\n", t2.Sub(t0))

	// Debug code
	// for i := 1; i < len(edgesArr); i++ {
	// 	fmt.Printf("[Smfy]edgesArr[%d]: %v\n", i, edgesArr[i])
	// 	if len(edgesArr[i].Utg.Ks) == len(edgesArr[i-1].Utg.Ks) {
	// 		fmt.Printf("[Smfy] outcoming: %v, %v, len: %d\n", edgesArr[i-1].Utg.Ks[Kmerlen-1], edgesArr[i].Utg.Ks[Kmerlen-1], len(edgesArr[i].Utg.Ks))
	// 	}
	// }
	// Debug code
	smfyNodesfn := opt.Prefix + ".nodes.smfy.Arr.br"
	fc1 := make(chan int, 1)
	go NodesArrWriter(nodesArr, smfyNodesfn, fc1)

	smfyEdgesfn := opt.Prefix + ".edges.smfy.fq.br"
	StoreEdgesToFn(smfyEdgesfn, edgesArr)
	<-fc1
	//mappingEdgefn := opt.Prefix + ".edges.mapping.fa"
	// StoreMappingEdgesToFn(mappingEdgefn, edgesArr, opt.MaxMapEdgeLen)
	//	adpaterEdgesfn := prefix + ".edges.adapter.fq"
	//	StoreEdgesToFn(adpaterEdgesfn, edgesArr, true)

	DBGInfofn := opt.Prefix + ".smfy.DBGInfo"
	DBGInfoWriter(DBGInfofn, len(edgesArr), len(nodesArr))

	//EdgesStatWriter(edgesStatfn, len(edgesArr))

	// output graphviz graph
	if opt.Graph {
		graphfn := opt.Prefix + ".afterSmfy.dot"
		GraphvizDBGArr(nodesArr, edgesArr, graphfn)
	}
}

func CountByte(ba []byte, s byte) (count int) {
	for _, c := range ba {
		if c == s {
			count++
		}
	}
	return count
}

func SplitEIDArr(ba []byte, s byte) [][]byte {
	n := CountByte(ba, s)
	a := make([][]byte, n+1)
	na := 0
	//bytes.Split(s, sep)
	start := 0
	for i, c := range ba {
		if c == s {
			if i > start {
				a[na] = ba[start:i]
				na++
			}
			start = i + 1
		}
	}
	if start < len(ba) {
		a[na] = ba[start:]
		na++
	}
	return a[:na]
}

func WriteNGSPath(pathfn string, wc chan ReadPath, numThreads int, writeNumC chan<- int) int {
	pathfp, err := os.Create(pathfn)
	if err != nil {
		log.Fatalf("[WriteNGSPath] file %s create error, err: %v\n", pathfn, err)
	}
	//sbuf := make([]byte, 0, 100)
	defer pathfp.Close()
	cbrofp := cbrotli.NewWriter(pathfp, cbrotli.WriterOptions{Quality: 1})
	defer cbrofp.Close()
	buffp := bufio.NewWriterSize(cbrofp, 1<<18) // 1<<25 == 2**25
	var finishNum int
	readNum := 0
	//var pathArr PathArr
	for {
		rp := <-wc
		path := rp.GetPath0()
		//extArr := rp.GetPath0()
		//fmt.Fprintf(os.Stderr, "[WriteNGSPath] path: %v\n", extArr)
		if len(path) == 0 {
			finishNum++
			if finishNum == numThreads {
				break
			}
		}
		if len(path) < 3 {
			continue
		}
		var s string
		for _, eID := range path {
			s += strconv.Itoa(int(eID)) + " "
		}
		s += "\n"
		buffp.WriteString(s)
		//pathArr.Arr = append(pathArr.Arr, &rp)
		readNum++
	}

	/*data, err := proto.Marshal(&pathArr)
	if err != nil {
		log.Fatalf("[WriteNGSPath] Mashal data error: %v\n", err)
	}

	buffp.Write(data)*/

	if err := buffp.Flush(); err != nil {
		log.Fatalf("[WriteNGSPath] failed to flush file: %s, err: %v\n", pathfn, err)
	}

	if err := cbrofp.Flush(); err != nil {
		log.Fatalf("[WriteNGSPath] failed to flush file: %s, err: %v\n", pathfn, err)
	}
	close(wc)
	//time.Sleep(time.Second * 10)
	fmt.Printf("[WriteNGSPath] write NGSPath num: %d to file: %s\n", readNum, pathfn)
	writeNumC <- readNum
	return readNum
}

func LoadNGSPathFromFile(pathfn string, cs chan<- []DBG_MAX_INT) {
	bufSize := (1 << 15)
	fncs := make(chan []byte, 2)
	go constructcf.ReadBrFile(pathfn, fncs, bufSize)
	buf := make([]byte, 0, 1<<16)
	//lb := make([]byte, 100)
	var lb []byte
	pos := 0
	var suc bool
	for {
		lb, pos, suc = constructcf.GetNextLine(buf, pos)
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
		idList := SplitEIDArr(lb, ' ')
		eIDArr := make([]DBG_MAX_INT, len(idList))
		for i, bs := range idList {
			//eID, err := utils.ByteArrInt(id)
			eID, err := strconv.Atoi(string(bs))
			if err != nil {
				log.Fatalf("[LoadPathFromFile]%s transform eID error\n", bs)
			}
			eIDArr[i] = DBG_MAX_INT(eID)
		}
		//fmt.Printf("[LoadPathFromFile]eIDArr:%v\n", eIDArr)
		cs <- eIDArr
	}
	close(cs)
	return
}

type MapingInfo struct {
	ID        uint32
	PathMat   [][][]uint32
	Anotition string
}

func LoadPathFromFile(pathfn string, cs chan<- MapingInfo) {
	bufSize := (1 << 15)
	fncs := make(chan []byte, 2)
	go constructcf.ReadBrFile(pathfn, fncs, bufSize)
	buf := make([]byte, 0, 1<<16)
	//lb := make([]byte, 100)
	var lb []byte
	pos := 0
	var suc bool
	for {
		lb, pos, suc = constructcf.GetNextLine(buf, pos)
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
		//eIDArr := make([]DBG_MAX_INT, len(colList))
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
						//eID, err := utils.ByteArrInt(id)
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
	*/

	return
}

func LoadNGSFile(fafn string, cs chan<- constructcf.ReadInfo, kmerlen int) {
	bufSize := (1 << 18)
	fncs := make(chan []byte, 3)
	format := constructcf.GetReadsFileFormat(fafn)

	go constructcf.ReadBrFile(fafn, fncs, bufSize)
	num := constructcf.GetReadFileRecord(fncs, cs, format, bufSize, true)
	fmt.Printf("[LoadNGSFile] processed %d reads fromfile: %s\n", num, fafn)
}

func MappingNGSRead(dk DBGKmer, ri constructcf.ReadInfo, rpos int, rstrand bool, edgesArr []DBGEdge, nodesArr []DBGNode, kmerlen int) (errorNum, mappingNum int, path []DBG_MAX_INT) {
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
			if len(e.Utg.Ks)-int(dk.Pos) < len(ri.Seq)-i {
				b = i + (len(e.Utg.Ks) - int(dk.Pos))
			}
			j = int(dk.Pos)
			for ; i < b && j < e.GetSeqLen(); i++ {
				if ri.Seq[i] != e.Utg.Ks[j] {
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
				if ri.Seq[i] != bnt.BntRev[e.Utg.Ks[j]] {
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
					eb = ne.Utg.Ks[dk.Pos]
				} else {
					dk.Pos = int32(ne.GetSeqLen() - kmerlen)
					eb = bnt.BntRev[ne.Utg.Ks[dk.Pos]]
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
}

func MapNGSReadFindPath(cs <-chan constructcf.ReadInfo, wc chan<- ReadPath, nodesArr []DBGNode, edgesArr []DBGEdge, cf CuckooFilter, WinSize int, kmerlen int) {
	var notFoundSeedNum, highErrorNum, shortMapingNum, allNum int
	var totalMapLen int
	MaxStepNum := 10
	ka := make([]KmerInfo, MaxStepNum*WinSize)
	rSeq := make([]byte, 300)
	for {
		ri, ok := <-cs
		if !ok {
			var t ReadPath
			wc <- t
			break
		}
		allNum++

		dbgK, pos, strand := LocateSeedKmerCF(cf, ri, WinSize, edgesArr, MaxStepNum, ka, rSeq)
		if dbgK.GetCount() == 0 { // not found in the cuckoofilter
			//fmt.Printf("[paraMapNGSAndMerge] read ID: %v not found seed!!!\n", pairRI[j].ID)
			notFoundSeedNum++
			continue
		}
		var rp ReadPath
		errorNum, mappingNum, path := MappingNGSRead(dbgK, ri, int(pos), strand, edgesArr, nodesArr, cf.Kmerlen)
		totalMapLen += mappingNum
		if errorNum > len(ri.Seq)/80 {
			highErrorNum++
			continue
		}
		if mappingNum < len(ri.Seq)*8/10 {
			shortMapingNum++
			continue
		}
		if len(path) < 2 {
			continue
		}
		rp.Path0 = make([]uint32, len(path))
		for j, id := range path {
			rp.Path0[j] = uint32(id)
		}
		wc <- rp
	}
	fmt.Printf("[MapNGSReadFindPath] not found seed read pair number is : %v,allNum: %v,  percent: %v\n", notFoundSeedNum, allNum, float32(notFoundSeedNum)/float32(allNum))
	fmt.Printf("[MapNGSReadFindPath] high error mapping read pair number is : %v, short mapping num: %v,avg mapping len: %v\n", highErrorNum, shortMapingNum, totalMapLen/(allNum-notFoundSeedNum))
}

func ParaMapReadsFile(fafn string, wc chan ReadPath, concurrentNum int, nodesArr []DBGNode, edgesArr []DBGEdge, cf CuckooFilter, opt Options, processT chan int) {
	bufSize := 1000
	cs := make(chan constructcf.ReadInfo, bufSize)
	for i := 0; i < concurrentNum; i++ {
		go MapNGSReadFindPath(cs, wc, nodesArr, edgesArr, cf, opt.WinSize, opt.Kmer)
	}
	LoadNGSFile(fafn, cs, opt.Kmer)
	processT <- 1
}

func MapingNGSFindDBGPath(opt Options, nodesArr []DBGNode, edgesArr []DBGEdge) {
	// construct cuckoofilter of DBG sample
	cfSize := GetCuckoofilterDBGSampleSize(edgesArr, int64(opt.WinSize), int64(opt.Kmer))
	fmt.Printf("[MapingNGSFindDBGPath] cfSize: %v\n", cfSize)
	cf := MakeCuckooFilter(uint64(cfSize*2), opt.Kmer)
	fmt.Printf("[MapingNGSFindDBGPath] cf.NumItems: %v\n", cf.NumItems)
	count := ConstructCFDBGMinimizers(cf, edgesArr, opt.WinSize)
	fmt.Printf("[MapingNGSFindDBGPath]construct Smaple of DBG edges cuckoofilter number is : %v\n", count)

	cfgInfo, err := constructcf.ParseCfg(opt.CfgFn, false, false)
	if err != nil {
		log.Fatal("[constructcf.ParseCfg] found err")
	}
	fmt.Printf("[MapingNGSFindDBGPath] cfgInfo: %v\n", cfgInfo)

	runtime.GOMAXPROCS(opt.NumCPU + 10)

	concurrentNum := opt.Kmer / 20
	totalNumT := (opt.NumCPU + concurrentNum - 1) / concurrentNum
	if concurrentNum < 2 || totalNumT < 2 {
		totalNumT = 1
		concurrentNum = 1
	}
	processT := make(chan int, totalNumT)
	writeNumC := make(chan int, 1)
	wc := make(chan ReadPath, 1000)
	for i := 0; i < totalNumT; i++ {
		processT <- 1
	}
	var fnNum int
	for _, lib := range cfgInfo.Libs {
		if lib.AsmFlag != constructcf.AllState && lib.SeqProfile != 1 {
			continue
		}
		fnNum += len(lib.FnName) + len(lib.SingleFn)
	}

	pathfn := opt.Prefix + ".NGSPath.br"
	go WriteNGSPath(pathfn, wc, concurrentNum*fnNum, writeNumC)
	for _, lib := range cfgInfo.Libs {
		if lib.AsmFlag != constructcf.AllState && lib.SeqProfile != 1 {
			continue
		}
		fnNum += len(lib.FnName) + len(lib.SingleFn)
		for i := 0; i < len(lib.FnName); i++ {
			<-processT
			go ParaMapReadsFile(lib.FnName[i], wc, concurrentNum, nodesArr, edgesArr, cf, opt, processT)
		}
		for i := 0; i < len(lib.SingleFn); i++ {
			<-processT
			go ParaMapReadsFile(lib.SingleFn[i], wc, concurrentNum, nodesArr, edgesArr, cf, opt, processT)
		}
	}

	for i := 0; i < totalNumT; i++ {
		<-processT
	}
	<-writeNumC
	//time.Sleep(time.Second * 10)
}

func Transform2Path(p32 []uint32) (path []DBG_MAX_INT) {
	path = make([]DBG_MAX_INT, len(p32))
	for i, id := range p32 {
		path[i] = DBG_MAX_INT(id)
	}
	return
}
func Transform2Path2(p32 []uint32, path []DBG_MAX_INT) []DBG_MAX_INT {
	la := len(p32)
	if len(path) < la {
		path = make([]DBG_MAX_INT, la)
	} else {
		path = path[:la]
	}
	for i, id := range p32 {
		path[i] = DBG_MAX_INT(id)
	}
	return path
}

func Transform2EdgeFreqArr(p32 []DBG_MAX_INT, path []EdgeFreq) []EdgeFreq {
	la := len(p32)
	if cap(path) < la {
		path = make([]EdgeFreq, la)
	} else {
		path = path[:la]
	}
	for i, id := range p32 {
		path[i].ID = id
		path[i].Freq = 1
	}
	return path
}

func EqualDBG_MAX_INTArr(arr1, arr2 []DBG_MAX_INT) bool {
	if len(arr1) != len(arr2) || len(arr1) == 0 || len(arr2) == 0 {
		return false
	}

	ok := true
	for i, c := range arr1 {
		if c != arr2[i] {
			ok = false
			break
		}
	}
	return ok
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

func IndexEdgeFreqArr(path []EdgeFreq, eID DBG_MAX_INT) (idx int) {
	idx = -1
	for i, ef := range path {
		if ef.ID == eID {
			idx = i
			break
		}
	}

	return
}

func EdgeFreqArr2DBG_MAX_INTArr(ep []EdgeFreq) []DBG_MAX_INT {
	arr := make([]DBG_MAX_INT, len(ep))
	for i, ef := range ep {
		arr[i] = ef.ID
	}

	return arr
}

func AddPathToEdge2(e *DBGEdge, path, rpath []EdgeFreq) *DBGEdge {
	ok := false
	eID := e.ID
	idx := IndexEdgeFreqArr(path, eID)
	ridx := len(path) - 1 - idx
	//ridx := IndexEdgeFreqArr(rpath, eID)
	//fmt.Printf("[AddPathToEdge2]eID:%d, NGSPathArr:%v, path:%v\n", e.ID, e.NGSPathArr, path)
	for i, np := range e.NGSPathArr {
		nidx := IndexEdgeFreqArr(np, eID)
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
}

/*
func AddPathToEdge(edgesArr []DBGEdge, nodesArr []DBGNode, e DBGEdge, path, rpath []DBG_MAX_INT, MaxPathSize int) {
	ok := false
	eID := e.ID
	bArr := GetNextEdgeIDArr(e, nodesArr, PLUS, BACKWARD)
	fArr := GetNextEdgeIDArr(e, nodesArr, PLUS, FORWARD)
	idx := IndexEID(path, eID)
	var p []DBG_MAX_INT
	if idx > 0 {
		id := path[idx-1]
		if IsInDBG_MAX_INTArr(bArr, id) {
			p = path
		} else {
			p = rpath
		}
	} else {
		id := path[idx+1]
		if IsInDBG_MAX_INTArr(fArr, id) {
			p = path
		} else {
			p = rpath
		}
	}
	idx = IndexEID(p, eID)
	pl := len(p)
	for i, np := range e.NGSPathArr {
		nidx := IndexEID(np.IDArr, eID)
		eo := false
		if pl-idx <= len(np.IDArr)-nidx {
			if pl-idx == 1 && len(np.IDArr)-nidx == 1 && EqualDBG_MAX_INTArr(p[idx:], np.IDArr[nidx:nidx+pl-idx]) {
				eo = true
			} else if pl-idx > 1 && EqualDBG_MAX_INTArr(p[idx:], np.IDArr[nidx:nidx+pl-idx]) {
				eo = true
			}
		} else {
			if len(np.IDArr)-nidx > 1 && EqualDBG_MAX_INTArr(p[idx:idx+len(np.IDArr)-nidx], np.IDArr[nidx:]) {
				eo = true
			}
		}
		if eo {
			eo = false
			if idx <= nidx {
				if EqualDBG_MAX_INTArr(p[:idx+1], np.IDArr[nidx-idx:nidx+1]) {
					eo = true
				}
			} else {
				if EqualDBG_MAX_INTArr(p[idx-nidx:idx+1], np.IDArr[:nidx+1]) {
					eo = true
				}
			}
			if eo {
				if idx > nidx || pl-idx > len(np.IDArr)-nidx {
					var al int
					if idx > nidx {
						al = idx
					} else {
						al = nidx
					}
					if pl-idx > len(np.IDArr)-nidx {
						al += pl - idx
					} else {
						al += len(np.IDArr) - nidx
					}
					arr := make([]DBG_MAX_INT, al)
					fa := make([]uint8, al)
					tidx := nidx
					if idx > nidx {
						copy(arr, p[:idx+1])
						copy(fa[idx-nidx:], edgesArr[eID].NGSPathArr[i].FreqArr)
						tidx = idx
					} else {
						copy(arr, edgesArr[eID].NGSPathArr[i].IDArr[:nidx+1])
						copy(fa, edgesArr[eID].NGSPathArr[i].FreqArr)
					}
					if pl-idx > len(np.IDArr)-nidx {
						copy(arr[tidx+1:], p[idx+1:])
					} else {
						copy(arr[tidx+1:], edgesArr[eID].NGSPathArr[i].IDArr[nidx+1:])
					}
					nidx = tidx
					edgesArr[eID].NGSPathArr[i].IDArr = arr
					edgesArr[eID].NGSPathArr[i].FreqArr = fa
				}
				// add freq
				np = edgesArr[eID].NGSPathArr[i]
				for j := nidx - idx; j < nidx+(pl-idx); j++ {
					f := np.FreqArr[j]
					if f < math.MaxUint8 {
						np.FreqArr[j]++
					}
				}
				edgesArr[eID].NGSPathArr[i] = np
				ok = true
				break
			}
		}
	}

	if !ok {
		var np PathFreq
		np.IDArr = make([]DBG_MAX_INT, pl)
		np.FreqArr = make([]uint8, pl)
		copy(np.IDArr, p)
		for j := range np.FreqArr {
			np.FreqArr[j]++
		}

		if len(e.NGSPathArr) == MaxPathSize {
			pa := e.NGSPathArr
			for x, pf := range pa {
				freq := GetPathFreq(pf)
				if freq == 1 {
					edgesArr[eID].NGSPathArr[x] = np
					ok = true
					break
				}
			}
		}

		if !ok {
			edgesArr[eID].NGSPathArr = append(edgesArr[eID].NGSPathArr, np)
		}
	}
	return
} */

func DeleteEID(ea []DBG_MAX_INT, eID DBG_MAX_INT) []DBG_MAX_INT {
	if len(ea) == 0 {
		return ea
	}
	na := make([]DBG_MAX_INT, 0, len(ea))
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

type IDFreq struct {
	ID   DBG_MAX_INT
	Freq int
}

func AddIDFreqArr(idFreqArr []IDFreq, eID DBG_MAX_INT, freq int) []IDFreq {
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
		if idFreq.Freq >= MinMapFreq {
			idFreqArr[count] = idFreq
			count++
		}
	}
	idFreqArr = idFreqArr[:count]
	return idFreqArr
}

func GetNumFreqInPathFreqArr(pa []PathFreq, eID DBG_MAX_INT) (cn int) {
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

func GetNumInPathFreqArr(pa []PathFreq, eID DBG_MAX_INT) (cn int, idx int) {
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

/*func GetUniquePathIdxArr(pa []PathFreq, pf PathFreq, eID DBG_MAX_INT, edgesArr []DBGEdge, EqualFlag bool) (idxArr, ix2Arr []int) {
	l1 := len(pf.IDArr)
	idxArr = make([]int, 0, l1/3+1)
	ix2Arr = make([]int, 0, l1/3+1)
	idx := IndexEID(pf.IDArr, eID)
	for x, id := range pf.IDArr {
		if !EqualFlag && x == idx {
			continue
		}

		e2 := edgesArr[id]
		if len(e2.NGSPathArr) == 0 {
			continue
		}
		c1, _ := GetNumInPathFreqArr(pa, id)
		if c1 == 1 {
			c2, ix2 := GetNumInPathFreqArr(e2.NGSPathArr, eID)
			if c2 == 1 {
				pf2 := e2.NGSPathArr[ix2]
				l2 := len(pf2.IDArr)
				if !EqualFlag {
					if l1 == l2 && ((pf.IDArr[0] == pf2.IDArr[0] && pf.IDArr[l1-1] == pf2.IDArr[l2-1]) || (pf.IDArr[0] == pf2.IDArr[l2-1] && pf.IDArr[l1-1] == pf2.IDArr[0])) {
						continue
					}
				}

				idx2 := IndexEID(pf2.IDArr, eID)
				y := IndexEID(pf2.IDArr, id)
				if x < idx {
					if y > idx2 {
						pf2.IDArr = GetReverseDBG_MAX_INTArr(pf2.IDArr)
						y = l2 - 1 - y
					}
					if y < x || l2-y > l1-x {
						continue
					}
				} else {
					if y < idx2 {
						pf2.IDArr = GetReverseDBG_MAX_INTArr(pf2.IDArr)
						y = l2 - 1 - y
					}
					if l1-x > l2-y || y > x {
						continue
					}
				}
				lmin := utils.MinInt(x, y)
				rmin := utils.MinInt(l1-x, l2-y)
				if EqualDBG_MAX_INTArr(pf.IDArr[x-lmin:x+rmin], pf2.IDArr[y-lmin:y+rmin]) {
					idxArr = append(idxArr, x)
					ix2Arr = append(ix2Arr, ix2)
				}
			}
		}
	}
	return
}*/

func GetShareNodeID(e, e1 *DBGEdge) DBG_MAX_INT {
	if e.StartNID == e1.StartNID || e.StartNID == e1.EndNID {
		return e.StartNID
	} else if e.EndNID == e1.StartNID || e.EndNID == e1.EndNID {
		return e.EndNID
	}
	var nID DBG_MAX_INT
	return nID
}

func GetNumInDBG_MAX_INTArr(arr []DBG_MAX_INT, eID DBG_MAX_INT) (count int) {
	for _, id := range arr {
		if id == eID {
			count++
		}
	}
	return
}

/*
func GetMaxUniqueRegion(idxArr []int, pf PathFreq, edgesArr []DBGEdge) (idx1, idx2 int) {
	idx1, idx2 = -1, -1
	for i, idx := range idxArr {
		eID := pf.IDArr[idx]
		e := edgesArr[eID]
		if e.StartNID == e.EndNID {
			continue
		}
		if len(e.NGSPathArr) == 0 {
			continue
		}

		if len(e.NGSPathArr) == 1 {
			idx1 = i
			break
		}
		if i >= len(idxArr)-1 {
			break
		}
		if GetNumInDBG_MAX_INTArr(pf.IDArr, eID) != 1 {
			continue
		}
		e1 := edgesArr[pf.IDArr[idx+1]]
		c0, _ := GetNumInPathFreqArr(e.NGSPathArr, e1.ID)
		if c0 == 1 {
			nID := GetShareNodeID(e, e1)
			c1 := 0
			if nID == e.EndNID {
				for _, pf := range e.NGSPathArr {
					if pf.IDArr[len(pf.IDArr)-1] == eID {
						continue
					} else {
						c1++
					}
				}
			} else {
				for _, pf := range e.NGSPathArr {
					if pf.IDArr[0] == eID {
						continue
					} else {
						c1++
					}
				}
			}
			if c1 == 1 {
				idx1 = i
				break
			}
		}
	}

	if idx1 < 0 || idx1 >= len(idxArr)-1 {
		return
	}

	for i := len(idxArr) - 1; i > idx1; i-- {
		idx := idxArr[i]
		eID := pf.IDArr[idx]
		e := edgesArr[eID]
		if e.StartNID == e.EndNID {
			continue
		}
		if len(e.NGSPathArr) == 0 {
			continue
		}

		if len(e.NGSPathArr) == 1 {
			idx2 = i
			break
		}
		if GetNumInDBG_MAX_INTArr(pf.IDArr, eID) != 1 {
			continue
		}
		e1 := edgesArr[pf.IDArr[idx-1]]
		c0, _ := GetNumInPathFreqArr(e.NGSPathArr, e1.ID)
		if c0 == 1 {
			nID := GetShareNodeID(e, e1)
			c1 := 0
			if nID == e.StartNID {
				for _, pf := range e.NGSPathArr {
					if pf.IDArr[0] == eID {
						continue
					} else {
						c1++
					}
				}
			} else {
				for _, pf := range e.NGSPathArr {
					if pf.IDArr[len(pf.IDArr)-1] == eID {
						continue
					} else {
						c1++
					}
				}
			}
			if c1 == 1 {
				idx2 = i
				break
			}
		}
	}

	return
}*/

func AddToIdxArr(idxArr []int, ix2Arr []int, idx, ix int) ([]int, []int) {
	idxArr = append(idxArr, idx)
	ix2Arr = append(ix2Arr, ix)
	for i := len(idxArr) - 2; i >= 0; i-- {
		if idxArr[i] > idxArr[i+1] {
			idxArr[i], idxArr[i+1] = idxArr[i+1], idxArr[i]
			ix2Arr[i], ix2Arr[i+1] = ix2Arr[i+1], ix2Arr[i]
		} else {
			break
		}
	}

	return idxArr, ix2Arr
}

func GetPathSeqLen(path []DBG_MAX_INT, edgesArr []DBGEdge, kmerlen int) (sl int) {
	sl += kmerlen - 1
	for _, eID := range path {
		e := edgesArr[eID]
		sl += e.GetSeqLen() - (kmerlen - 1)
	}
	return sl
}

func GetEdgeFreqPathSeqLen(path []EdgeFreq, edgesArr []DBGEdge, kmerlen int) (sl int) {
	sl += kmerlen - 1
	for _, ef := range path {
		e := edgesArr[ef.ID]
		sl += e.GetSeqLen() - (kmerlen - 1)
	}
	return sl
}

func ConcatEdgePath(path []DBG_MAX_INT, edgesArr []DBGEdge, nodesArr []DBGNode, kmerlen int) (u Unitig) {
	sl := GetPathSeqLen(path, edgesArr, kmerlen)
	fmt.Printf("[ConcatEdgePath]path:%v seqLen:%d\n", path, sl)
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
		copy(u.Ks[pos:], e.Utg.Ks)
		copy(u.Kq[pos:], e.Utg.Kq)
	} else {
		direction = FORWARD
		copy(u.Ks, e.Utg.Ks)
		copy(u.Kq, e.Utg.Kq)
		pos = e.GetSeqLen()
	}
	fmt.Printf("[ConcatEdgePath]direction:%v pos:%d\n\teID:%d %s\n", direction, pos, e.ID, Transform2Char(e.Utg.Ks))
	for i := 1; i < len(path); i++ {
		eID := path[i]
		ne := &edgesArr[eID]
		nl := ne.GetSeqLen()
		strand = GetNextEdgeStrand(e, ne, nodesArr, strand)
		if direction == FORWARD {
			if strand == PLUS {
				copy(u.Ks[pos:], ne.Utg.Ks[kn:])
				copy(u.Kq[pos:], ne.Utg.Kq[kn:])
			} else {
				copy(u.Ks[pos:], ne.Utg.Ks[:nl-kn])
				copy(u.Kq[pos:], ne.Utg.Kq[:nl-kn])
				ReverseCompByteArr(u.Ks[pos : pos+nl-kn])
				ReverseUint8Arr(u.Kq[pos : pos+nl-kn])
			}
			pos += nl - kn
		} else {
			tl := nl - kn
			if strand == PLUS {
				copy(u.Ks[pos-tl:pos], ne.Utg.Ks[:tl])
				copy(u.Kq[pos-tl:pos], ne.Utg.Kq[:tl])
			} else {
				copy(u.Ks[pos-tl:pos], ne.Utg.Ks[nl-tl:])
				copy(u.Kq[pos-tl:pos], ne.Utg.Kq[nl-tl:])
				ReverseCompByteArr(u.Ks[pos-tl : pos])
				ReverseUint8Arr(u.Kq[pos-tl : pos])
			}
			pos -= tl
		}
		if strand == PLUS {
			fmt.Printf("[ConcatEdgePath2]strand:%v pos:%d ne.StartNID:%d ne.EndNID:%d eID:%d %s\n", strand, pos, ne.StartNID, ne.EndNID, ne.ID, Transform2Char(ne.Utg.Ks))
		} else {
			fmt.Printf("[ConcatEdgePath2]strand:%v pos:%d ne.StartNID:%d ne.EndNID:%d eID:%d %s\n", strand, pos, ne.StartNID, ne.EndNID, ne.ID, Transform2Char(GetReverseCompByteArr(ne.Utg.Ks)))
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

func ConcatEdgePath2(path []DBG_MAX_INT, edgesArr []DBGEdge, nodesArr []DBGNode, kmerlen int) (u Unitig) {
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
		copy(u.Ks[pos:], e.Utg.Ks)
		copy(u.Kq[pos:], e.Utg.Kq)
	} else {
		direction = FORWARD
		copy(u.Ks, e.Utg.Ks)
		copy(u.Kq, e.Utg.Kq)
		pos = e.GetSeqLen()
	}
	rSeq := make([]byte, 0, 500)
	rQual := make([]byte, 0, 500)
	fmt.Printf("[ConcatEdgePath2]direction:%v pos:%d\n\teID:%d %s\n", direction, pos, e.ID, Transform2Char(e.Utg.Ks))
	for i := 1; i < len(path); i++ {
		eID := path[i]
		ne := &edgesArr[eID]
		nl := ne.GetSeqLen()
		strand = GetNextEdgeStrand(e, ne, nodesArr, strand)
		ks := ne.Utg.Ks
		kq := ne.Utg.Kq
		if strand == MINUS {
			rSeq = GetReverseCompByteArr2(ks, rSeq)
			rQual = GetReverseByteArr2(kq, rQual)
			ks = rSeq
			kq = rQual
		}
		if direction == FORWARD {
			copy(u.Ks[pos:], ks[kn:])
			copy(u.Kq[pos:], kq[kn:])
			pos += nl - kn
		} else {
			tl := nl - kn
			copy(u.Ks[pos-tl:pos], ks[:tl])
			copy(u.Kq[pos-tl:pos], kq[:tl])
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

func PrintEdgeInfo(e *DBGEdge) string {
	s := fmt.Sprintf("eID:%d StartNID:%d EndNID:%d Incoming:%v Outcoming:%v el:%d", e.ID, e.StartNID, e.EndNID, e.EdgeIDIncoming, e.EdgeIDOutcoming, e.GetSeqLen())
	return s
}

func PrintNodeInfo(nd *DBGNode) string {
	s := fmt.Sprintf("nID:%d Incoming:%v Outcoming:%v", nd.ID, nd.EdgeIDIncoming, nd.EdgeIDOutcoming)
	return s
}

func IsConnectNode(v DBGNode, eID1, eID2 DBG_MAX_INT) bool {
	var ok bool
	if IsInComing(v.EdgeIDIncoming, eID1) {
		if IsInComing(v.EdgeIDOutcoming, eID2) {
			ok = true
		}
	} else if IsInComing(v.EdgeIDOutcoming, eID1) {
		if IsInComing(v.EdgeIDIncoming, eID2) {
			ok = true
		}
	}
	return ok
}

func ConcatEdgePathUtg(path []EdgeFreq, nID DBG_MAX_INT, edgesArr []DBGEdge, nodesArr []DBGNode, kmerlen int) (e DBGEdge) {
	//sl := GetEdgeFreqPathSeqLen(path, edgesArr, kmerlen)
	//fmt.Printf("[ConcatEdgePath2]path:%v seqLen:%d\n", path, sl)
	//kn := kmerlen - 1
	e = edgesArr[path[0].ID]
	//e.Ks = make([]byte, sl)
	//e.Kq = make([]uint8, sl)
	e1, e2 := &edgesArr[path[0].ID], &edgesArr[path[1].ID]
	if e1.StartNID == e1.EndNID {
		if IsInComing(e1.EdgeIDIncoming, e2.ID) {
			if e2.EndNID == e1.EndNID {
				e.Utg = GetRCUnitig(e1.Utg)
				e.StartNID, e.EndNID = e.EndNID, e.StartNID
				e.EdgeIDIncoming, e.EdgeIDOutcoming = e.EdgeIDOutcoming, e.EdgeIDIncoming
			}
		} else {
			if e2.EndNID == e1.EndNID {
				e.Utg = GetRCUnitig(e1.Utg)
				e.StartNID, e.EndNID = e.EndNID, e.StartNID
				e.EdgeIDIncoming, e.EdgeIDOutcoming = e.EdgeIDOutcoming, e.EdgeIDIncoming
			}
		}
	} else if e1.EndNID != nID {
		e.Utg = GetRCUnitig(e1.Utg)
		e.StartNID, e.EndNID = e.EndNID, e.StartNID
		e.EdgeIDIncoming, e.EdgeIDOutcoming = e.EdgeIDOutcoming, e.EdgeIDIncoming
	}
	//fmt.Printf("[ConcatEdgePathUtg]nID:%d %s\n", nID, PrintEdgeInfo(&e))
	//strand := PLUS
	var ne DBGEdge
	lastE := &e
	strand := PLUS
	for i := 1; i < len(path); i++ {
		ne = edgesArr[path[i].ID]
		strand = GetNextEdgeStrand2(lastE, &ne, strand, FORWARD)
		if strand == MINUS {
			ne.Utg = GetRCUnitig(ne.Utg)
			ne.StartNID, ne.EndNID = ne.EndNID, ne.StartNID
			ne.EdgeIDIncoming, ne.EdgeIDOutcoming = ne.EdgeIDOutcoming, ne.EdgeIDIncoming
		}
		//fmt.Printf("[ConcatEdgePathUtg]strand:%t %s\n", strand, PrintEdgeInfo(&ne))
		e.Utg = ConcatEdges(e.Utg, ne.Utg, kmerlen)
		if len(e.Utg.Ks) == 0 {
			return
		}
		e.EndNID = ne.EndNID
		lastE = &edgesArr[path[i].ID]
	}

	return
}

func RePlaceNGSPath2(path, mergePath, revMP []EdgeFreq) []EdgeFreq {
	eID1, eID2 := mergePath[0].ID, mergePath[len(mergePath)-1].ID
	idx1 := IndexEdgeFreqArr(path, eID1)
	idx2 := IndexEdgeFreqArr(path, eID2)
	if idx1 >= 0 || idx2 >= 0 {
		fmt.Printf("[RePlaceNGSPath2]idx1:%d idx2:%d path:%v mergePath:%v\n", idx1, idx2, path, mergePath)
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
		fmt.Printf("\tchanged path:%v\n", path)
	}
	return path
}

/*
func RePlaceNGSPath(pf PathFreq, path []DBG_MAX_INT, eID DBG_MAX_INT) PathFreq {
	idx := IndexEID(pf.IDArr, path[0])
	idx2 := IndexEID(pf.IDArr, path[len(path)-1])
	ok := false
	if idx >= 0 && idx2 >= 0 && idx+idx2 > 0 {
		if idx > idx2 {
			path = GetReverseDBG_MAX_INTArr(path)
		}
		idx, idx2 = idx2, idx
		min := utils.MinInt(len(pf.IDArr)-idx, len(path))
		if EqualDBG_MAX_INTArr(pf.IDArr[idx:idx+min], path[:min]) {
			ok = true
			var npf PathFreq
			pl := idx + 1
			if len(pf.IDArr)-idx > len(path) {
				pl += len(pf.IDArr) - idx - len(path)
			}
			npf.IDArr = make([]DBG_MAX_INT, pl)
			copy(npf.IDArr[:idx], pf.IDArr[:idx])
			npf.IDArr[idx] = eID
			if len(pf.IDArr)-idx > len(path) {
				copy(npf.IDArr[idx+1:], pf.IDArr[idx+len(path):])
			}
			pf = npf
		}
	} else if idx >= 0 && idx2 < 0 {
		var pl int
		if len(pf.IDArr)-idx < len(path) {
			pl = len(pf.IDArr) - idx
		} else {
			pl = len(path)
		}
		if EqualDBG_MAX_INTArr(pf.IDArr[idx:idx+pl], path[:pl]) {
			ok = true
			var npf PathFreq
			npf.IDArr = make([]DBG_MAX_INT, idx+1)
			copy(npf.IDArr[:idx], pf.IDArr[:idx])
			npf.IDArr[idx] = eID
			pf = npf
		}
	} else if idx < 0 && idx2 >= 0 {
		var pl int
		if idx2 < len(path) {
			pl = idx2
		} else {
			pl = len(path)
		}
		if EqualDBG_MAX_INTArr(pf.IDArr[idx2-pl:idx2+1], path[len(path)-pl:]) {
			ok = true
			var npf PathFreq
			pl := 1
			if idx2+1 > len(pf.IDArr) {
				pl += len(pf.IDArr) - (idx2 + 1)
			}
			npf.IDArr = make([]DBG_MAX_INT, pl)
			npf.IDArr[0] = eID
			if idx2+1 < len(pf.IDArr) {
				copy(npf.IDArr[1:], pf.IDArr[idx2+1:])
			}
			pf = npf
		}
	}

	if !ok {
		fmt.Fprintf(os.Stderr, "[RePlaceNGSPath]pf:%v not found path:%v\n", pf, path)
	}
	return pf
}*/

func IsIntArr(id int, arr []int) bool {
	ok := false
	for _, t := range arr {
		if t == id {
			ok = true
			break
		}
	}
	return ok
}

func GetNumInPath0(pa []PathFreq, eID DBG_MAX_INT) (count int) {
	for _, pf := range pa {
		if pf.IDArr[0] == eID {
			count++
		}
	}
	return
}

func FindEdgeFreqArrBoundary(catPath []EdgeFreq) (idx1, idx2 int) {
	idx1, idx2 = -1, -1
	for i := 0; i <= len(catPath)/2; i++ {
		if catPath[i].ID != 0 {
			idx1 = i
			break
		}
	}
	for i := len(catPath) - 1; i >= len(catPath)/2; i-- {
		if catPath[i].ID != 0 {
			idx2 = i
			break
		}
	}
	return
}

func GetExtendUniquePath(e *DBGEdge, edgesArr []DBGEdge, minAvgEdgeFreq int) (efArr []EdgeFreq, eIDArr []DBG_MAX_INT) {
	max := 1000
	// extend left path
	ep := e.NGSPathArr[0]
	rArr := make([]EdgeFreq, 0, len(ep)*3)
	shareArr := make([]EdgeFreq, 0, 20)
	eIdx := IndexEdgeFreqArr(ep, e.ID)
	efArr = GetReverseEdgeFreqArr(ep[:eIdx+1], rArr)
	eIDArr = make([]DBG_MAX_INT, 0, 5)
	eIDArr = append(eIDArr, e.ID)
	fmt.Printf("[GetExtendUniquePath]eID:%v, ep:%v\n", e.ID, ep)
	//idx := IndexEdgeFreqArr(efArr, eIDArr[len(eIDArr)-1])
	for i := 1; i < len(efArr) && i < max; i++ {
		if int(efArr[i].Freq) < minAvgEdgeFreq*2 {
			break
		}
		eID := efArr[i].ID
		te := &edgesArr[eID]
		if !IsUniqueEdge2(te) || len(te.NGSPathArr) != 1 {
			continue
		}

		if CountEIDArr(eIDArr, te.ID) > 0 {
			break
		}
		if te.GetProcessFlag() > 0 {
			//idx := IndexEdgeFreqArr(efArr[:i], te.ID)
			fmt.Printf("[GetExtendUniquePath]error! eID:%v has processed,eIDArr:%v NGSPath:%v\n", te.ID, eIDArr, te.NGSPathArr)
			break
			//log.Fatalf("[GetExtendUniquePath]error! eID:%v has processed\n", te.ID)
		}
		avgFreq := GetAvgEdgeFreq(te.NGSPathArr[0])
		if avgFreq < minAvgEdgeFreq*2 {
			break
		}
		lastEID := eIDArr[len(eIDArr)-1]
		ep1 := te.NGSPathArr[0]
		idx1 := IndexEdgeFreqArr(ep1, lastEID)
		idx2 := IndexEdgeFreqArr(ep1, eID)
		if idx1 >= 0 && idx2 >= 0 {
			idx3 := IndexEdgeFreqArr(efArr, lastEID)
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
			fmt.Printf("[GetExtendUniquePath]eID:%v, ep1:%v\n\tefArr:%v\n", eID, ep1, efArr)
			eIDArr = append(eIDArr, te.ID)
		}
	}
	efArr = ReverseEdgeFreqArr(efArr)
	eIDArr = GetReverseDBG_MAX_INTArr(eIDArr)
	if eIdx < len(ep)-1 {
		efArr = append(efArr, ep[eIdx+1:]...)
	}

	// extend right path
	leftIdx := IndexEdgeFreqArr(efArr, eIDArr[len(eIDArr)-1])
	for i := leftIdx + 1; i < len(efArr) && i < max; i++ {
		if int(efArr[i].Freq) < minAvgEdgeFreq*2 {
			break
		}
		eID := efArr[i].ID
		te := &edgesArr[eID]
		if !IsUniqueEdge2(te) || len(te.NGSPathArr) != 1 {
			continue
		}
		if CountEIDArr(eIDArr, te.ID) > 0 {
			break
		}
		if te.GetProcessFlag() > 0 {
			//idx := IndexEdgeFreqArr(efArr[:i], te.ID)
			fmt.Printf("[GetExtendUniquePath]error! eID:%v has processed,eIDArr:%v NGSPath:%v \n", te.ID, eIDArr, te.NGSPathArr)
			break
			//log.Fatalf("[GetExtendUniquePath]error! eID:%v has processed\n", te.ID)
		}
		avgFreq := GetAvgEdgeFreq(te.NGSPathArr[0])
		if avgFreq < minAvgEdgeFreq*2 {
			break
		}
		lastEID := eIDArr[len(eIDArr)-1]
		ep1 := te.NGSPathArr[0]
		idx1 := IndexEdgeFreqArr(ep1, lastEID)
		idx2 := IndexEdgeFreqArr(ep1, eID)
		if idx1 >= 0 && idx2 >= 0 {
			idx3 := IndexEdgeFreqArr(efArr, lastEID)
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
			fmt.Printf("[GetExtendUniquePath]eID:%v, ep1:%v efArr:%v\n", te.ID, ep1, efArr)
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

func GetAvgEdgeFreq(ep []EdgeFreq) (avgFreq int) {
	avgFreq = 0
	for _, ef := range ep {
		avgFreq += int(ef.Freq)
	}
	avgFreq /= len(ep)

	return avgFreq
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

func AddNGSPathToDBG(NGSPathfn string, nodesArr []DBGNode, edgesArr []DBGEdge, opt Options) {
	bufSize := 10000
	minAvgEdgeFreq := opt.MinMapFreq
	cs := make(chan []DBG_MAX_INT, bufSize)
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
			/*if e.ID == 223169 {
				fmt.Printf("[AddNGSPathToDBG]pArr:%v rPath:%v\n", pArr, rPath)
			}*/
			AddPathToEdge2(e, pArr, rPath)
			/*if e.ID == 223169 {
				fmt.Printf("[AddNGSPathToDBG]NGSPath:%v\n", edgesArr[ef.ID].NGSPathArr)
			}*/
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
					x1, x2 := IndexEdgeFreqArr(p0, p1[0].ID), IndexEdgeFreqArr(p0, p1[len(p1)-1].ID)
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
			ix := IndexEdgeFreqArr(px, e.ID)
			rPath := GetReverseEdgeFreqArr(px, rPath)
			ixR := len(px) - 1 - ix
			for y, py := range e.NGSPathArr {
				if x == y || deleteFlagArr[y] {
					continue
				}
				iy := IndexEdgeFreqArr(py, e.ID)
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

	/*// check code
	{
		for i, e := range edgesArr {
			if e.GetUniqueFlag() > 0 {
				fmt.Printf("[AddNGSPathToDBG]eID:%d NGSPath:%v\n", i, e.NGSPathArr)
			}
		}
	}*/

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
		idx1, idx2 := IndexEdgeFreqArr(efArr, eIDArr[0]), IndexEdgeFreqArr(efArr, eIDArr[len(eIDArr)-1])

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
				idx := IndexEID(eIDArr, ef.ID)
				if idx >= 0 && idx+1 < len(eIDArr) {
					eID1 = eID2
					eID2 = eIDArr[idx+1]
					continue
				}
				pa := edgesArr[ef.ID].NGSPathArr
				nArr := pa[:0]
				for _, ep := range pa {
					idx5 := IndexEdgeFreqArr(ep, eID1)
					idx6 := IndexEdgeFreqArr(ep, eID2)
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
			if IsInDBG_MAX_INTArr(eIDArr[1:], id) {
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
				idx := IndexEID(pf.IDArr, eID)
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
				pfe.PF.IDArr = make([]DBG_MAX_INT, al)
				pfe.PF.FreqArr = make([]uint8, al)
				idx := IndexEID(pf.IDArr, eID)
				copy(pfe.PF.IDArr[leftMax-idx:], pf.IDArr)
				copy(pfe.PF.FreqArr[leftMax-idx:], pf.FreqArr)
				pfe.Start = leftMax - idx
				pfe.End = pfe.Start + len(pf.IDArr)
				pArr[j] = pfe
			}

			stk := make([]PathFreqExt, 0, 10)
			var pfe PathFreqExt
			pfe.PF.IDArr = make([]DBG_MAX_INT, al)
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
						if EqualDBG_MAX_INTArr(p1.PF.IDArr[j+1:p1.End], pfe.PF.IDArr[j+1:p1.End]) {
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
								npf.PF.IDArr = make([]DBG_MAX_INT, al)
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
			rpfe.PF.IDArr = make([]DBG_MAX_INT, al)
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
							start = utils.MaxInt(p1.Start, pfe.Start)
						} else {
							if p1.Start == leftMax-1 {
								start = p1.Start
							} else {
								start = utils.MinInt(p1.Start, pfe.Start)
							}
						}
						if EqualDBG_MAX_INTArr(p1.PF.IDArr[start+1:p1.End], pfe.PF.IDArr[start+1:p1.End]) {
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
								npf.PF.IDArr = make([]DBG_MAX_INT, al)
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
						pf.IDArr = make([]DBG_MAX_INT, l0+l1-1)
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
							idx2 := IndexEID(pf2.IDArr, eID)
							if len(pf2.IDArr)-idx2 >= len(pf.IDArr) && EqualDBG_MAX_INTArr(pf.IDArr, pf2.IDArr[idx2:idx2+len(pf.IDArr)]) {
								ok = false
								break
							} else if c0 == 1 && len(pf2.IDArr)-idx2 < len(pf.IDArr) && EqualDBG_MAX_INTArr(pf.IDArr[:len(pf2.IDArr)-idx2], pf2.IDArr[idx2:]) {
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
					idx := IndexEID(pf.IDArr, eID)
					fmt.Printf("[AddNGSPathToDBG]eID:%v pf:%v idxArr:%v\n", eID, pf, idxArr)
					var npf PathFreq
					if idxArr[0] < idx {
						x := idxArr[0]
						e2 := edgesArr[pf.IDArr[x]]
						pf2 := e2.NGSPathArr[ix2Arr[0]]
						idx2 := IndexEID(pf2.IDArr, eID)
						l2 := len(pf2.IDArr)
						y := IndexEID(pf2.IDArr, e2.ID)
						if y > idx2 {
							pf2.IDArr = GetReverseDBG_MAX_INTArr(pf2.IDArr)
							y = l2 - 1 - y
						}
						lmin := utils.MinInt(x, y)
						rmin := utils.MinInt(l1-x, l2-y)
						if !EqualDBG_MAX_INTArr(pf.IDArr[x-lmin:x+rmin], pf2.IDArr[y-lmin:y+rmin]) {
							log.Fatalf("[AddNGSPathToDBG]pf:%v pf2:%v not found consis path", pf.IDArr, pf2.IDArr)
						}
						npf.IDArr = make([]DBG_MAX_INT, y+l1-x)
						copy(npf.IDArr[:y], pf2.IDArr[:y])
						copy(npf.IDArr[y:], pf.IDArr[x:])
					}

					if idxArr[len(idxArr)-1] > idx {
						x := idxArr[len(idxArr)-1]
						e2 := edgesArr[pf.IDArr[x]]
						pf2 := e2.NGSPathArr[ix2Arr[len(ix2Arr)-1]]
						idx2 := IndexEID(pf2.IDArr, eID)
						l2 := len(pf2.IDArr)
						y := IndexEID(pf2.IDArr, e2.ID)
						if y < idx2 {
							pf2.IDArr = GetReverseDBG_MAX_INTArr(pf2.IDArr)
							y = l2 - 1 - y
						}
						//lmin := utils.MinInt(x, y)
						//rmin := utils.MinInt(l1-x, l2-y)
						if !EqualDBG_MAX_INTArr(pf.IDArr[x-y:], pf2.IDArr[:y+l1-x]) {
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
					rnpf.IDArr = GetReverseDBG_MAX_INTArr(npf.IDArr)
					for j, t := range idxArr {
						p2 := ix2Arr[j]
						id := pf.IDArr[t]
						e2 := edgesArr[id]
						pf2 := e2.NGSPathArr[p2]
						idx2 := IndexEID(pf2.IDArr, eID)
						fmt.Printf("[AddNGSPathToDBG]e2ID:%v e2 pf:%v\n", id, pf2)
						y := IndexEID(pf2.IDArr, id)
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
		pathMap := make(map[DBG_MAX_INT][]DBG_MAX_INT)
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
				npf.IDArr = make([]DBG_MAX_INT, len(pf.IDArr))
				copy(npf.IDArr, pf.IDArr)
				pf = npf
				idxArr, _ := GetUniquePathIdxArr(pa, pf, eID, edgesArr, EqualFlag)
				if len(idxArr) < 1 {
					continue
				}
				fmt.Printf("[AddNGSPathToDBG]j:%d pf:%v idxArr:%v\n", j, pf, idxArr)
				//idx := IndexEID(pf.IDArr, eID)
				//if len(edgesArr[eID].NGSPathArr) == 1 {
				//	idxArr, ix2Arr = AddToIdxArr(idxArr, ix2Arr, idx, j)
				//}
				idx1, idx2 := GetMaxUniqueRegion(idxArr, pf, edgesArr)
				fmt.Printf("[AddNGSPathToDBG]idxArr:%v, idx1:%d idx2:%d\n", idxArr, idx1, idx2)
				if idx1 >= 0 && idx2 >= 0 && idx1 < idx2 {
					start, end := idxArr[idx1], idxArr[idx2]
					path := make([]DBG_MAX_INT, end+1-start)
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
							if v[0] == path[0] && EqualDBG_MAX_INTArr(v, path) {

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
					e.Utg = utg
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
				npf.IDArr = make([]DBG_MAX_INT, len(pf.IDArr))
				copy(npf.IDArr, pf.IDArr)
				pf = npf
				idxArr, _ := GetUniquePathIdxArr(pa, pf, eID, edgesArr, EqualFlag)
				if len(idxArr) < 1 {
					continue
				}
				fmt.Printf("[AddNGSPathToDBG]j:%d pf:%v idxArr:%v\n", j, pf, idxArr)
				//idx := IndexEID(pf.IDArr, eID)
				//if len(edgesArr[eID].NGSPathArr) == 1 {
				//	idxArr, ix2Arr = AddToIdxArr(idxArr, ix2Arr, idx, j)
				//}
				idx1, idx2 := GetMaxUniqueRegion(idxArr, pf, edgesArr)
				fmt.Printf("[AddNGSPathToDBG]idxArr:%v, idx1:%d idx2:%d\n", idxArr, idx1, idx2)
				if idx1 >= 0 && idx2 >= 0 && idx1 < idx2 {
					start, end := idxArr[idx1], idxArr[idx2]
					path := make([]DBG_MAX_INT, end+1-start)
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
							if v[0] == path[0] && EqualDBG_MAX_INTArr(v, path) {

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
					e.Utg = utg
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
			edgesArr[eID].Utg = e0.Utg
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
							v = GetReverseDBG_MAX_INTArr(v)
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
	*/

	return
}

func MapNGS(c cli.Command) {

	t0 := time.Now()
	// check agruments
	gOpt, suc := utils.CheckGlobalArgs(c.Parent())
	if suc == false {
		log.Fatalf("[MapNGS] check global Arguments error, opt: %v\n", gOpt)
	}
	opt := Options{gOpt, 0, 0, 0, 0, false, false}
	tmp, suc := checkArgs(c)
	if suc == false {
		log.Fatalf("[MapNGS] check Arguments error, opt: %v\n", tmp)
	}
	DebugModel := true
	if DebugModel {
		profileFn := opt.Prefix + ".MapNGS.prof"
		cpuprofilefp, err := os.Create(profileFn)
		if err != nil {
			log.Fatalf("[MapNGS] open cpuprofile file: %v failed\n", profileFn)
		}
		pprof.StartCPUProfile(cpuprofilefp)
		defer pprof.StopCPUProfile()
	}

	opt.MaxNGSReadLen = tmp.MaxNGSReadLen
	opt.WinSize = tmp.WinSize
	opt.MinMapFreq = tmp.MinMapFreq
	fmt.Printf("Arguments: %v\n", opt)

	// read files and construt DBG
	DBGInfofn := opt.Prefix + ".smfy.DBGInfo"
	eSize, nSize := DBGInfoReader(DBGInfofn)
	nodesfn := opt.Prefix + ".nodes.smfy.Arr.br"
	nodesArr := make([]DBGNode, nSize)
	fc := make(chan int, 1)
	go NodesArrReader(nodesfn, nodesArr, opt.Kmer, fc)
	edgesfn := opt.Prefix + ".edges.smfy.fq.br"
	//LoadEdgesfqFromFn(edgesfn, edgesArr, true)
	edgesArr := ReadEdgesFromFile(edgesfn, DBG_MAX_INT(eSize))
	// check NodesArrReader() finished
	<-fc

	CheckInterConnectivity(edgesArr, nodesArr)

	uniqueNum, bubbleRepeatNum, twoCycleNum, selfCycle, selfCycleSameComingNum, bubbleEdgeNum := SetDBGEdgesUniqueFlag(edgesArr, nodesArr, opt.Kmer)
	fmt.Printf("[MapNGS] unique edge number is : %v, bubbleRepeatNum: %v, twoCycleNum:%v, selfCycle:%v, selfCycleSameComingNum: %v, bubbleEdgeNum: %v\n", uniqueNum, bubbleRepeatNum, twoCycleNum, selfCycle, selfCycleSameComingNum, bubbleEdgeNum)
	AddNodeInfo2DBGEdgeArr(edgesArr, nodesArr)
	//PrintTmpDBG(nodesArr, edgesArr, opt.Prefix)

	// map Illumina reads to the DBG and find reads map path
	NGSPathfn := opt.Prefix + ".NGSPath.br"
	if _, err := os.Stat(NGSPathfn); err != nil {
		MapingNGSFindDBGPath(opt, nodesArr, edgesArr)
		t2 := time.Now()
		fmt.Printf("[MapNGS]MapingNGSFindDBGPath: %v\n", t2.Sub(t0))
	}
	AddNGSPathToDBG(NGSPathfn, nodesArr, edgesArr, opt)
	SmfyDBG(nodesArr, edgesArr, opt, false)
	//edgesArr = SetEdgeID(nodesArr, edgesArr)
	//nodesArr = SetNodeID(nodesArr, edgesArr)
	CheckDBGSelfCycle(nodesArr, edgesArr, opt.Kmer)
	CheckInterConnectivity(edgesArr, nodesArr)

	//store DBG
	mapDBGNodesfn := opt.Prefix + ".nodes.MapDBG.Arr.br"
	fc2 := make(chan int, 1)
	go NodesArrWriter(nodesArr, mapDBGNodesfn, fc2)

	mapDBGEdgesfn := opt.Prefix + ".edges.MapDBG.fq.br"
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
	fmt.Printf("[MapNGS]MapNGS total used : %v\n", t2.Sub(t0))

}
