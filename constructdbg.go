package main

import (
	"bufio"
	"encoding/binary"
	"encoding/gob"
	"fmt"
	"io"
	"log"
	"math"
	"os"
	"reflect"
	"runtime"
	"sort"
	"strconv"
	"strings"
	"sync"
	"time"

	// "time"
	"github.com/awalterschulze/gographviz"
	"github.com/jwaldrip/odin/cli"
	"github.com/klauspost/compress/zstd"
)

type DBGNode struct {
	ID              uint32
	EdgeIDIncoming  [BaseTypeNum]uint32 // the ID of EdgeID inclming
	EdgeIDOutcoming [BaseTypeNum]uint32
	Seq             []uint64 // kmer node seq, use 2bit schema
	//SubGID          uint8    // Sub graph ID, most allow 2**8 -1 subgraphs
	Flag uint8 // from low~high, 1:Process, 2:
}

func (n *DBGNode) String() string {
	return fmt.Sprintf("ID:%d EdgeIncoming:%v EdgeOutcoming:%v\n", n.ID, n.EdgeIDIncoming, n.EdgeIDOutcoming)
}

type NodeInfo struct {
	ID              uint32
	EdgeIDIncoming  [BaseTypeNum]uint32 // the ID of EdgeID inclming
	EdgeIDOutcoming [BaseTypeNum]uint32
}

const (
	DIN             uint8 = 1 // note DBGNode and DBGEdge relationship
	DOUT            uint8 = 2 // note DBGNode and DBGEdge relationship
	PLUS, MINUS           = true, false
	NODEMAP_KEY_LEN       = 5
)

type NodeInfoByEdge struct {
	n1, n2 DBGNode
	c1, c2 uint8
	i1, i2 uint
	Flag   bool
	edgeID uint32
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
	IDArr  []uint32
	AltArr []uint32
	Freq   int
}

type EdgeFreq struct {
	ID   uint32
	Freq uint8
}

type PathFreq struct {
	IDArr   []uint32
	FreqArr []uint8
}

type Path1 struct {
	IDArr []uint32
	//AltArr []uint32
	Freq int
}

const seedLen = 21
const seqPositionBitNum = 64 - 1 - seedLen*NumBitsInBase
const saxAllowSeqLen = (1 << seqPositionBitNum) - 1
const kmerMask = (1 << (seqPositionBitNum + 1)) - 1
const positionMask = (1 << seqPositionBitNum) - 1
const kmerReSetPosition = (((1 << (seedLen * NumBitsInBase)) - 1) << (seqPositionBitNum + 1)) | 0x1
const kmerResetStrand = ((1 << (seedLen*NumBitsInBase + seqPositionBitNum)) - 1) << 1

type SeedInfo uint64

type KI struct {
	SeedInfo        // kmerBase|Position|strand [SeedLen*NumBitsInBase:seqPositionBitNum:1] just store 2-bits base, allow max 32 kmer base
	ID       uint32 // edgeID or reads ID
	//Info uint32 // the first high 31-bit for Position,and last bit for strand info
}

func (k SeedInfo) GetKmer() uint64 {
	return uint64(k >> (seqPositionBitNum + 1))
}
func (k *SeedInfo) SetKmer(kmer uint64) {
	//fmt.Printf("[SetKmer]kmer:%b\n", kmer)
	*k = SeedInfo(kmer<<(seqPositionBitNum+1)) | (*k & kmerMask)
	//fmt.Printf("[SetKmer]k.Kmer:%b\n", k.GetKmer())
}

func (k SeedInfo) GetPosition() uint64 {
	return uint64((k >> 1) & positionMask)
}
func (k *SeedInfo) SetPosition(p uint64) {
	if p >= (1 << seqPositionBitNum) {
		log.Fatalf("[SetPosition]Pos:%d must < %d\n", p, 1<<seqPositionBitNum)
	}
	//fmt.Printf("[SetPosition]Pos:%b\n", p)
	*k = SeedInfo((p << 1) | (uint64(*k) & kmerReSetPosition))
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
	ID              uint32
	StartNID        uint32              // start node ID
	EndNID          uint32              // end node ID
	EdgeIDIncoming  [BaseTypeNum]uint32 // incoming edgeIDArr
	EdgeIDOutcoming [BaseTypeNum]uint32 // outcoming edgeIDArr
	CovD            uint16              // Coverage Depth or number of  low freq minimers
	Flag            uint8               //
	Ks              []byte
	//Utg              Unitig
	NGSPathArr        [][]EdgeFreq
	PathMat, PathMat1 []Path // read Path matrix
	//GFAPath           []*[][][]uint32
	SeedInfoArr []SeedInfo // store seq seed array
}

func (e *DBGEdge) String() string {
	return fmt.Sprintf("eID:%d StartNID:%d EndNID:%d el:%d EdgeIDIncoming:%v EdgeIDOutcoming:%d\n", e.ID, e.StartNID, e.EndNID, len(e.Ks), e.EdgeIDIncoming, e.EdgeIDOutcoming)
}

/*
func Equaluint32(path1, path2 []uint32) bool {
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

/*func InsertPathToEdge(pathMat []Path, path []uint32, freq int) []Path {

	// check have been added
	added := false
	for i, v := range pathMat {
		if Equaluint32(v.IDArr, path) {
			pathMat[i].Freq += freq
			added = true
			break
		}
	}
	if added == false {
		var np Path
		np.IDArr = make([]uint32, len(path))
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
	return len(e.Ks)
}

//var MIN_KMER_COUNT uint16 = 3
var BACKWARD uint8 = 0
var FORWARD uint8 = 1

//var Kmerlen int = -1

func Reverseuint32Arr(path []uint32) []uint32 {
	al := len(path)
	dl := al / 2
	for i := 0; i < dl; i++ {
		path[i], path[al-1-i] = path[al-1-i], path[i]
	}

	return path
}

func GetKmerRecord(fncs <-chan []byte, cs chan<- KmerBntBucket, bufSize, kmerlen int, kmerNumC chan<- int) {
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
	var buck KmerBntBucket
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
			copy(buf[0:], buf[pos:])
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
		var rb KmerBnt
		rb.Seq = b
		rb.Len = kmerlen
		if buck.Count >= ReadSeqSize {
			cs <- buck
			var nb KmerBntBucket
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

type KmerPool struct {
	KmerArr [][]byte
	Cs      []byte
}

func GetUniqKmer(fn string, rbytesPool *sync.Pool, kmerArrPool *sync.Pool, kmerPoolChan chan<- KmerPool, kmerlen int) {
	//bufSize := WindowSize
	cs := make(chan []byte, 4)
	go ReadZstdFile(fn, rbytesPool, cs)
	var notFoundNewKmer []byte
	var num int
	fmt.Printf("[GetUniqKmer] begin processe file:%s...\n", fn)

	for {
		bs, ok := <-cs
		if !ok {
			if len(notFoundNewKmer) != 0 {
				log.Fatalf("[GetUniqKmer] to the file end, but len(notFoundNewRecord) = %d\n", len(notFoundNewKmer))
			}
			break
		}

		var kp KmerPool
		kp.KmerArr = kmerArrPool.Get().([][]byte)
		kp.KmerArr = kp.KmerArr[:0]
		kp.Cs = bs
		idx := 0
		if len(notFoundNewKmer) > 0 {
			idx += kmerlen - len(notFoundNewKmer)
			notFoundNewKmer = append(notFoundNewKmer, kp.Cs[0:idx]...)
			kp.KmerArr = append(kp.KmerArr, notFoundNewKmer)
		}
		for ; idx <= len(kp.Cs)-kmerlen; idx += kmerlen {
			kp.KmerArr = append(kp.KmerArr, kp.Cs[idx:idx+kmerlen])
		}
		num += len(kp.KmerArr)
		var nk []byte
		if idx < len(kp.Cs) {
			nk = append(nk, kp.Cs[idx:]...)
		}
		notFoundNewKmer = nk
		kmerPoolChan <- kp
	}

	// send read finish signal
	close(kmerPoolChan)
	fmt.Printf("[GetUniqKmer] processed kmer number is:%d,finished processed file:%s\n", num, fn)
	return
}

/*func UniDirectExtend(kb, rkb KmerBnt, cf CuckooFilter, min_kmer_count uint16, direction uint8) (ec int) {
	kmerlen := cf.Kmerlen
	//var nBnt ReadBnt
	//nSeq = make([]byte, kmerlen)
	if direction == FORWARD {
		//copy(nSeq[:kmerlen-1], nb.Seq)
		for i := 0; i < BaseTypeNum; i++ {
			b := uint64(i)
			//nSeq[kmerlen-1] = b
			ks := GetNextKmer(kb, b, kmerlen)
			rs := GetPreviousKmer(rkb, uint64(BntRev[b]), kmerlen)
			if ks.BiggerThan(rs) {
				ks = rs
			}
			if count := cf.GetCountAllowZero(ks.Seq); count >= min_kmer_count {
				ec++
			}
		}
	} else { // direction == BACKWARD
		//copy(nSeq[1:], nb.Seq)
		for i := 0; i < BaseTypeNum; i++ {
			b := uint64(i)
			//nSeq[0] = b
			ks := GetPreviousKmer(kb, b, kmerlen)
			rs := GetNextKmer(rkb, uint64(BntRev[b]), kmerlen)
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

func ExtendNodeKmer(nkb, rkb []byte, cf CuckooFilter, min_kmer_count uint16) (nd DBGNode) {
	kmerlen := cf.Kmerlen
	//rb2 = make([]byte, kmerlen)
	//var nBnt ReadBnt
	//nSeq = make([]byte, kmerlen+1)
	//nd.Seq = GetReadBntKmer(nkb.Seq, 0, kmerlen-1)
	//copy(nSeq[1:], nodeSeq)
	//rkb := ReverseComplet(nkb)
	//nd.Seq = nkb.Seq
	//fmt.Printf("[ExtendNodeKmer]nkb: %v\n\trkb: %v\n", nkb, rkb)
	var min []byte
	for i := 0; i < BaseTypeNum; i++ {
		bi := byte(i)
		//nSeq[0] = bi
		{
			//ks := GetPreviousKmer(nkb, bi, kmerlen)
			//rs := GetNextKmer(rkb, uint64(BntRev[bi]), kmerlen)
			//kb2 = NoAllocGetPreviousKmer(nkb, kb2, bi, kmerlen)
			//rb2 = NoAllocGetNextKmer(rkb, rb2, uint64(BntRev[bi]), kmerlen)
			nkb[0] = bi
			rkb[kmerlen] = BntRev[bi]
			min = nkb[:kmerlen]
			if BiggerThan(min, rkb[1:]) {
				min = rkb[1:]
			}
			count, _ := cf.Lookup(min)
			//fmt.Printf("[ExtendNodeKmer]min: %v, count: %v\n", min, count)
			if count >= min_kmer_count {
				//var nb ReadBnt
				//nb.Seq = append(nb.Seq, nSeq[:kmerlen-1]...)
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
			//nSeq[cf.Kmerlen] = bi
			//kb2 = NoAllocGetNextKmer(nkb, kb2, bi, kmerlen)
			//rb2 = NoAllocGetPreviousKmer(rkb, rb2, uint64(BntRev[bi]), kmerlen)
			nkb[kmerlen] = bi
			rkb[0] = BntRev[bi]
			min = nkb[1:]
			if BiggerThan(min, rkb[:kmerlen]) {
				min = rkb[:kmerlen]
			}
			count, _ := cf.Lookup(min)
			//fmt.Printf("[ExtendNodeKmer]min: %v, count: %v\n", min, count)
			if count >= min_kmer_count {
				//fmt.Printf("[ExtendNodeKmer]count: %v\n", count)
				//var nb ReadBnt
				//nb.Seq = append(nb.Seq, nSeq[2:]...)
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
	for i := 0; i < BaseTypeNum; i++ {
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
	var nkb KmerBnt
	nkb.Len = kmerlen - 1
	tmp := make([]uint64, len(nd.Seq))
	copy(tmp, nd.Seq)
	rs := ReverseCompletBnt(nd.Seq, kmerlen-1)
	//fmt.Printf("[GetMinDBGNode] node: %v\nRC node: %v\n", nodeBnt, rnode)
	if BiggerThanBnt(nd.Seq, rs) {
		minN.Seq = rs
		for i := 0; i < BaseTypeNum; i++ {
			minN.EdgeIDIncoming[i] = nd.EdgeIDOutcoming[BaseTypeNum-1-i]
			minN.EdgeIDOutcoming[i] = nd.EdgeIDIncoming[BaseTypeNum-1-i]
		}
	} else {
		minN.Seq = nd.Seq
		minN.EdgeIDIncoming = nd.EdgeIDIncoming
		minN.EdgeIDOutcoming = nd.EdgeIDOutcoming
	}
	return minN
}

func paraLookupComplexNode(cs chan KmerPool, rbytesPool, kmerArrPool *sync.Pool, wc chan []DBGNode, nodeArrPool *sync.Pool, cf CuckooFilter) {
	//fmt.Printf("[paraLookupComplexNode]goroutine\n")
	rk := make([]byte, cf.Kmerlen)
	kb := make([]byte, cf.Kmerlen+1)
	rkb := make([]byte, cf.Kmerlen+1)
	nodeArr := nodeArrPool.Get().([]DBGNode)
	nodeArr = nodeArr[:cap(nodeArr)]
	count := 0
	var min []byte
	for {
		kp, ok := <-cs
		if !ok {
			break
		}
		// if rsb.Count < ReadSeqSize {
		// 	fmt.Printf("rsb.ReadBuf length is : %d\n", len(rsb.ReadBuf))
		// }
		// if found kmer count is 1 , this kmer will be ignore, and skip this branch
		for _, k := range kp.KmerArr {
			rk = GetReverseCompletBytes(k, rk)
			//fmt.Printf("[paraLookupComplexNode] kb : %v\n", kb)
			for i := 0; i < 2; i++ { // check fisrt node of kmer
				copy(kb[1:], k[i:i+cf.Kmerlen-1])
				copy(rkb[1:], rk[1-i:cf.Kmerlen-i])
				//nkb, _ := DeleteLastBaseKmer(kb)
				//fmt.Printf("[paraLookupComplexNode] nkb : %v\n", nkb)
				//rkb := ReverseComplet(nkb)
				//fmt.Printf("[paraLookupComplexNode] rkb : %v\n", rkb)
				//extkb := ExtendKmerBnt2Byte(kb)
				//extnkb := ExtendKmerBnt2Byte(nkb)
				//extrkb := ExtendKmerBnt2Byte(rkb)
				//fmt.Printf("[paraLookupComplexNode]first kb: %v\n\tbase: %v,nkb: %v\n\trkb: %v\n", ExtendKmerBnt2Byte(kb), base, ExtendKmerBnt2Byte(nkb), ExtendKmerBnt2Byte(rkb))
				//nodeLength = len(nodeSeq)
				nd := ExtendNodeKmer(kb, rkb, cf, 1)
				var leftcount, rightcount int
				for i := 0; i < BaseTypeNum; i++ {
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
					min = k[i : i+cf.Kmerlen-1]
					if BiggerThan(min, rk[1-i:cf.Kmerlen-i]) {
						min = rk[1-i : cf.Kmerlen-i]
						for i := 0; i < BaseTypeNum; i++ {
							wd.EdgeIDIncoming[i] = nd.EdgeIDOutcoming[BaseTypeNum-1-i]
							wd.EdgeIDOutcoming[i] = nd.EdgeIDIncoming[BaseTypeNum-1-i]
						}
					} else {
						wd.EdgeIDIncoming = nd.EdgeIDIncoming
						wd.EdgeIDOutcoming = nd.EdgeIDOutcoming
					}
					wd.Seq = GetReadBntKmer(min, 0, cf.Kmerlen-1)
					//copy(wd.Seq, nkb.Seq)
					//tn := GetMinDBGNode(nd, cf.Kmerlen)
					//fmt.Printf("[paraLookupComplexNode] nd: %v\n", wd)
					//fmt.Printf("[paraLookupComplexNode] node: %v\n", wd)
					nodeArr[count] = wd
					count++
					if count == len(nodeArr) {
						nodeArr = nodeArr[:count]
						wc <- nodeArr
						nodeArr = nodeArrPool.Get().([]DBGNode)
						nodeArr = nodeArr[:cap(nodeArr)]
						count = 0
					}
				}
			}
		}
		kmerArrPool.Put(kp.KmerArr)
		rbytesPool.Put(kp.Cs)
	}
	if count > 0 {
		nodeArr = nodeArr[:count]
		wc <- nodeArr
	}
	var na []DBGNode
	na = nil
	wc <- na
	//fmt.Printf("[paraLookupComplexNode]finished goroutine\n")
}

func GetNodeRecord(complexKmerfn string, nodeArrChan chan<- []DBGNode, nodeArrPool *sync.Pool, NBntUint64Len int) {
	ckfp, err := os.Open(complexKmerfn)
	if err != nil {
		log.Fatal(err)
	}
	defer ckfp.Close()
	buffp := bufio.NewReaderSize(ckfp, 1<<16)
	nodeArr := nodeArrPool.Get().([]DBGNode)
	nodeArr = nodeArr[:cap(nodeArr)]
	//size := cap(nodeArr)
	count := 0
	for {
		if err = binary.Read(buffp, binary.LittleEndian, &nodeArr[count].Seq); err != nil {
			if err == io.EOF {
				break
			}
			log.Fatalf("[constructNodeMap] err: %v\n", err)
		}
		if err = binary.Read(buffp, binary.LittleEndian, &nodeArr[count].EdgeIDIncoming); err != nil {
			log.Fatalf("[constructNodeMap] err: %v\n", err)
		}
		if err = binary.Read(buffp, binary.LittleEndian, &nodeArr[count].EdgeIDOutcoming); err != nil {
			log.Fatalf("[constructNodeMap] err: %v\n", err)
		}
		count++
		if count == cap(nodeArr) {
			nodeArrChan <- nodeArr
			nodeArr = nodeArrPool.Get().([]DBGNode)
			nodeArr = nodeArr[:cap(nodeArr)]
			count = 0
		}
	}
	if count > 0 {
		nodeArr = nodeArr[:count]
		nodeArrChan <- nodeArr
	}
	close(nodeArrChan)
}

func constructNodeMap(complexKmerfn string, nodeMap *sync.Map, NBntUint64Len int) uint32 {
	nodeID := uint32(2)
	nodeArrChan := make(chan []DBGNode, 1)
	nodeArrPool := sync.Pool{New: func() interface{} {
		arr := make([]DBGNode, WindowSize/(NBntUint64Len*8+8*4))
		for i := 0; i < len(arr); i++ {
			arr[i].Seq = make([]uint64, NBntUint64Len)
		}
		return arr
	}}

	go GetNodeRecord(complexKmerfn, nodeArrChan, &nodeArrPool, NBntUint64Len)
	readnodeNum := 0
	var key [NODEMAP_KEY_LEN]uint64
	var value NodeInfo
	for {
		nodeArr, ok := <-nodeArrChan
		if !ok {
			break
		}
		// fmt.Fprintf(os.Stderr, "[constructNodeMap] node: %v\n", node)
		// n, err := ckfp.Read(b)
		// if n != NBntByteLen {
		// 	log.Fatalf("[constructNodeMap] read node seq err: n(%d) != NBntByteLen(%d\n)", n, NBntByteLen)
		// }
		for i := 0; i < len(nodeArr); i++ {
			nd := &nodeArr[i]
			copy(key[:], nd.Seq)
			value.EdgeIDIncoming = nd.EdgeIDIncoming
			value.EdgeIDOutcoming = nd.EdgeIDOutcoming
			value.ID = nodeID
			//fmt.Printf("[constructNodeMap] key: %v\n\tnode: %v\n", key, node)
			if _, ok := nodeMap.LoadOrStore(key, value); ok == false {
				nodeID++
				//fmt.Printf("[constructNodeMap] node: %v\n", node)
			} else {
				//fmt.Printf("[constructNodeMap] repeat node: %v\n", node)
			}
		}
		readnodeNum += len(nodeArr)
		nodeArrPool.Put(nodeArr)

	}

	fmt.Printf("[constructNodeMap] read node number is : %d\n", readnodeNum)

	return nodeID
}

func AddNodeToNodeMap(node DBGNode, nodeMap map[[NODEMAP_KEY_LEN]uint64]DBGNode, nodeID uint32) uint32 {
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

func CollectAddedDBGNode(anc <-chan DBGNode, nodeMap *sync.Map, nc chan<- DBGNode, readNodeMapFinishedC <-chan int) {
	narr := make([]DBGNode, 0, 1000)
loop:
	for {
		select {
		case nd := <-anc:
			narr = append(narr, nd)
		case <-readNodeMapFinishedC:
			break loop
		}
	}

	addedNum := 0
	narrFlag := make(chan int, 1)
	finished := false
	var key [NODEMAP_KEY_LEN]uint64
	var ni NodeInfo
	for {
		if addedNum < len(narr) && len(narrFlag) == 0 {
			narrFlag <- 1
		}
		select {
		case nd := <-anc:
			narr = append(narr, nd)
		case <-narrFlag:
			{
				copy(key[:], narr[addedNum].Seq)
				v, _ := nodeMap.Load(key)
				ni = v.(NodeInfo)
				added := false
				for i := 0; i < BaseTypeNum; i++ {
					if ni.EdgeIDIncoming[i] == 1 {
						added = true
						break
					}
					if ni.EdgeIDOutcoming[i] == 1 {
						added = true
						break
					}
				}
				if added {
					fmt.Printf("[CollectAddedDBGNode]narr[%d]:%v key:%v\n", addedNum, ni, key)
					narr[addedNum].EdgeIDIncoming = ni.EdgeIDIncoming
					narr[addedNum].EdgeIDOutcoming = ni.EdgeIDOutcoming
					nc <- narr[addedNum]
				}
				addedNum++
			}
		default:
			time.Sleep(time.Second)
			if len(anc) == 0 {
				finished = true
			}
		}

		if finished {
			break
		}
	}
	close(nc)
	//time.Sleep(time.Second * 5)
	fmt.Printf("[CollectAddedDBGNode] added node number is: %d\n", addedNum)
}

// ChangeNodeMap add new Node to the nodeMap and check node edge has been output
/*func ChangeNodeMap(nodeMap map[[NODEMAP_KEY_LEN]uint64]DBGNode, anc chan<- DBGNode, finishedC <-chan int, nIEC <-chan NodeInfoByEdge, flagNIEC chan<- NodeInfoByEdge, Kmerlen int, nodeID uint32) (nID uint32, edgeID uint32) {
	oldNodeID := nodeID
	edgeID = uint32(2)
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
					b := BntRev[nie.i2]
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

func GetEdgeIDComingCount(nd *DBGNode) (count int) {
	for i := 0; i < BaseTypeNum; i++ {
		if nd.EdgeIDIncoming[i] > 0 {
			count++
		}
		if nd.EdgeIDOutcoming[i] > 0 {
			count++
		}
	}
	return
}

var muRW sync.RWMutex

// ReadDBGNodeToChan  read DBG nodeMap and simultaneously add new node to the nodeMap
func ReadDBGNodeToChan(nodeMap *sync.Map, nc chan<- DBGNode, readNodeMapFinished chan<- int, NBntUint64Len int) {
	var nd DBGNode
	var ni NodeInfo
	var k [NODEMAP_KEY_LEN]uint64
	var edgeCount int
	nodeMap.Range(func(key, value interface{}) bool {
		nd.Seq = make([]uint64, NBntUint64Len)
		k = key.([NODEMAP_KEY_LEN]uint64)
		for i := 0; i < NBntUint64Len; i++ {
			nd.Seq[i] = k[i]
		}
		ni = value.(NodeInfo)
		nd.ID = ni.ID
		nd.EdgeIDIncoming = ni.EdgeIDIncoming
		nd.EdgeIDOutcoming = ni.EdgeIDOutcoming
		nc <- nd
		edgeCount += GetEdgeIDComingCount(&nd)
		return true
	})

	/*for _, v := range nodeMap {
		//fmt.Printf("[ReadDBGNodeToChan]v:%v\n", v)
		nc <- v
		//value.Flag = uint8(1)
		//nodeMap[string(value.Seq)] = value
		//}
	}*/

	// notice Function ChangeNodeMap() has finished read nodeMap
	//time.Sleep(time.Second)
	readNodeMapFinished <- 1
	close(readNodeMapFinished)
	fmt.Printf("[ReadDBGNodeToChan]edgeCount:%d\n", edgeCount)
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
		cseq[i] = BntRev[b]
	}
	return
}

func ReverseCompByteArr(seq []byte) {
	ls := len(seq)
	for i := 0; i < ls/2; i++ {
		seq[i], seq[ls-1-i] = BntRev[seq[ls-1-i]], BntRev[seq[i]]
	}
	if ls%2 != 0 {
		seq[ls/2] = BntRev[seq[ls/2]]
	}
}

func GetReverseCompByteArr(seq []byte) []byte {
	sl := len(seq)
	rv := make([]byte, sl)
	for i := 0; i < len(rv); i++ {
		rv[i] = BntRev[seq[sl-1-i]]
	}

	return rv
}

func GetReverseCompNtByteArr(seq []byte) []byte {
	sl := len(seq)
	rv := make([]byte, sl)
	for i := 0; i < len(rv); i++ {
		c := seq[sl-1-i]
		if c == 'A' {
			rv[i] = 'T'
		} else if c == 'C' {
			rv[i] = 'G'
		} else if c == 'G' {
			rv[i] = 'C'
		} else {
			rv[i] = 'A'
		}
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
	for i := 0; i < sl; i++ {
		rv[i] = BntRev[seq[sl-1-i]]
	}

	return rv
}

func ReverseUint8Arr(ua []uint8) {
	lua := len(ua)
	for i := 0; i < lua/2; i++ {
		ua[i], ua[lua-1-i] = ua[lua-1-i], ua[i]
	}
}

func GetEdges(cf CuckooFilter, nseq, nrseq []byte, bi byte, direction uint8, MIN_KMER_COUNT uint16, e *DBGEdge) (nd DBGNode) {
	nk := make([]byte, len(nseq)+2)
	rnk := make([]byte, len(nrseq)+2)
	ndkLen := len(nseq)
	//tb.Seq = make([]uint64, (kmerlen+NumBaseInUint64-1)/NumBaseInUint64)
	//seq := ExtendKmerBnt2Byte(kb)
	e.Ks = e.Ks[:0]
	//var node DBGNode
	if direction == FORWARD {
		e.Ks = append(e.Ks, nseq...)
		copy(nk[1:], nseq[1:])
		nk[1+ndkLen-1] = bi
		rnk[1] = BntRev[bi]
		copy(rnk[2:], nrseq[:ndkLen-1])
		for {
			e.Ks = append(e.Ks, bi)
			nd = ExtendNodeKmer(nk, rnk, cf, MIN_KMER_COUNT)
			var leftcount, rightcount int
			for i := 0; i < BaseTypeNum; i++ {
				if nd.EdgeIDIncoming[i] == 1 {
					leftcount++
				}
				if nd.EdgeIDOutcoming[i] == 1 {
					bi = uint8(i)
					rightcount++
				}
			}
			if rightcount == 0 {
				break
			} else if leftcount <= 1 && rightcount == 1 {
				copy(nk[1:], nk[2:1+ndkLen])
				nk[1+ndkLen-1] = bi
				copy(rnk[2:], rnk[1:1+ndkLen-1])
				rnk[1] = BntRev[bi]
				//nkb = GetNextKmer(nkb, uint64(baseBnt), cf.Kmerlen-1)
				//kb2 = NoAllocGetNextKmer(nkb, kb2, uint64(baseBnt), cf.Kmerlen-1)
				//rb2 = NoAllocGetPreviousKmer(nr, rb2, uint64(BntRev[baseBnt]), cf.Kmerlen-1)
			} else {
				if leftcount > 1 || rightcount > 1 {
					nd.Seq = GetReadBntKmer(nk[1:1+ndkLen], 0, ndkLen)
					//nd.Flag = leftcount + rightcount
					nd.EdgeIDIncoming[e.Ks[len(e.Ks)-ndkLen-1]] = math.MaxUint32
				}
				break
			}
		}
	} else {
		e.Ks = append(e.Ks, nseq...)
		ReverseByteArr(e.Ks[:ndkLen])
		nk[1] = bi
		copy(nk[2:], nseq[:ndkLen-1])
		copy(rnk[1:], nrseq[1:])
		rnk[1+ndkLen-1] = BntRev[bi]
		for {
			e.Ks = append(e.Ks, bi)
			nd = ExtendNodeKmer(nk, rnk, cf, MIN_KMER_COUNT)
			var leftcount, rightcount int
			for i := 0; i < BaseTypeNum; i++ {
				if nd.EdgeIDIncoming[i] == 1 {
					bi = uint8(i)
					leftcount++
				}
				if nd.EdgeIDOutcoming[i] == 1 {
					rightcount++
				}
			}
			if leftcount == 0 {
				ReverseByteArr(e.Ks)
				break
			} else if leftcount == 1 && rightcount <= 1 {
				copy(nk[2:], nk[1:ndkLen])
				nk[1] = bi
				copy(rnk[1:], rnk[2:1+ndkLen])
				rnk[1+ndkLen-1] = BntRev[bi]
				//nkb = GetPreviousKmer(nkb, uint64(baseBnt), cf.Kmerlen-1)
				//nr = GetNextKmer(nr, uint64(BntRev[baseBnt]), cf.Kmerlen-1)
				//kb2 = NoAllocGetPreviousKmer(nkb, kb2, uint64(baseBnt), cf.Kmerlen-1)
				//rb2 = NoAllocGetNextKmer(nr, rb2, uint64(BntRev[baseBnt]), cf.Kmerlen-1)
			} else {
				if leftcount > 1 || rightcount > 1 {
					nd.Seq = GetReadBntKmer(nk[1:1+ndkLen], 0, ndkLen)
					nd.EdgeIDOutcoming[e.Ks[len(e.Ks)-ndkLen-1]] = math.MaxUint32
				}
				// Reverse the seq and quality count
				ReverseByteArr(e.Ks)
				break
			}
		}
	}
	return
}

var mapSM sync.Mutex

//var mapRWMu sync.RWMutex

//type EdgeNode struct {
//	Edge         DBGEdge
//	NodeS, NodeE DBGNode // NodeS note Start Node, NodeE note End Node
//}

func paraGenerateDBGEdges(nc <-chan DBGNode, cf CuckooFilter, nodeMap *sync.Map, wc chan<- []DBGEdge, edgeArrPool *sync.Pool, anc chan<- DBGNode, nodeIDChan chan uint32, TipMaxLen int) {
	seq := make([]byte, cf.Kmerlen-1)
	rseq := make([]byte, cf.Kmerlen-1)
	NumNodeSeqLen := (cf.Kmerlen - 1 - 1 + NumBaseInUint64) / NumBaseInUint64
	bnt := make([]uint64, NumNodeSeqLen)
	ea := edgeArrPool.Get().([]DBGEdge)
	ea = ea[:cap(ea)]
	ec := 0 // count of ea
	var key [NODEMAP_KEY_LEN]uint64
	var ni NodeInfo
	var edgeCount, tipNum, tipTotalLen int
	var node, nd DBGNode
	var ok bool
	for {
		node, ok = <-nc
		if !ok {
			if ec > 0 {
				wc <- ea[:ec]
			}
			var ea []DBGEdge
			wc <- ea
			break
		}
		copy(bnt, node.Seq)
		// read edge seq from cuckoofilter
		//fmt.Printf("[paraGenerateDBGEdges]node:%v\n\tbnt:%v\n", node, bnt)
		ExtendKmerBnt2Byte2(bnt, seq)
		GetReverseCompletBytes(seq, rseq)
		//extRSeq := ExtendKmerBnt2Byte(kb)
		// leftcount, rightcount, _, _ := ExtendNodeKmer(extRBnt, cf, MIN_KMER_COUNT, FORWARD)
		for i := uint(0); i < BaseTypeNum; i++ {
			bi := byte(i)
			if node.EdgeIDIncoming[i] == 1 {
				edgeCount++
				// get edge sequence
				// init edge
				ea[ec].StartNID, ea[ec].EndNID = 0, 0
				nd = GetEdges(cf, seq, rseq, bi, BACKWARD, 1, &ea[ec])
				if len(nd.Seq) == 0 {
					tipTotalLen += len(ea[ec].Ks)
					tipNum++
				}
				//fmt.Printf("[paraGenerateDBGEdges]BACKWARD edge:%v\n\tnd:%v\n", edge, nd)
				//writedEdge := false
				if len(nd.Seq) > 0 || len(ea[ec].Ks) > TipMaxLen {
					hasW := false
					mapSM.Lock()
					if len(nd.Seq) > 0 {
						rs := ReverseCompletBnt(nd.Seq, cf.Kmerlen-1)
						if BiggerThanBnt(nd.Seq, rs) {
							copy(nd.Seq, rs)
							ReverseCompletIncoming(&nd)
						}
						copy(key[:], nd.Seq)
						ni.ID = nd.ID
						ni.EdgeIDIncoming = nd.EdgeIDIncoming
						ni.EdgeIDOutcoming = nd.EdgeIDOutcoming
						if v, ok := nodeMap.LoadOrStore(key, ni); ok {
							ni = v.(NodeInfo)
							for j := 0; j < BaseTypeNum; j++ {
								if nd.EdgeIDIncoming[j] == math.MaxUint32 {
									if ni.EdgeIDIncoming[j] == math.MaxUint32 {
										hasW = true
									} else {
										ni.EdgeIDIncoming[j] = math.MaxUint32
										nodeMap.Store(key, ni)
										ea[ec].StartNID = ni.ID
									}
									break
								}
								if nd.EdgeIDOutcoming[j] == math.MaxUint32 {
									if ni.EdgeIDOutcoming[j] == math.MaxUint32 {
										hasW = true
									} else {
										ni.EdgeIDOutcoming[j] = math.MaxUint32
										nodeMap.Store(key, ni)
										ea[ec].StartNID = ni.ID
									}
									break
								}
							}
						} else {
							id := <-nodeIDChan
							ni.ID = id
							id++
							nodeIDChan <- id
							ea[ec].StartNID = ni.ID
							nodeMap.Store(key, ni)
							nd.ID = ni.ID
							anc <- nd
						}
					}
					if !hasW {
						copy(key[:], node.Seq)
						v, _ := nodeMap.Load(key)
						ni = v.(NodeInfo)
						ni.EdgeIDIncoming[i] = math.MaxUint32
						nodeMap.Store(key, ni)
						ea[ec].EndNID = node.ID
						ec++
						if ec == cap(ea) {
							wc <- ea
							ea = edgeArrPool.Get().([]DBGEdge)
							ea = ea[:cap(ea)]
							ec = 0
						}
					}
					mapSM.Unlock()
				}
			}

			if node.EdgeIDOutcoming[i] == 1 {
				edgeCount++
				// get edge sequence
				ea[ec].StartNID, ea[ec].EndNID = 0, 0
				nd = GetEdges(cf, seq, rseq, bi, FORWARD, 1, &ea[ec])
				if len(nd.Seq) == 0 {
					tipTotalLen += len(ea[ec].Ks)
					tipNum++
				}
				//fmt.Printf("[paraGenerateDBGEdges]FORWARD edge:%v\n\tnd:%v\n", edge, nd)
				//writedEdge := false
				if len(nd.Seq) > 0 || len(ea[ec].Ks) > TipMaxLen {
					hasW := false
					mapSM.Lock()
					if len(nd.Seq) > 0 {
						rs := ReverseCompletBnt(nd.Seq, cf.Kmerlen-1)
						if BiggerThanBnt(nd.Seq, rs) {
							copy(nd.Seq, rs)
							ReverseCompletIncoming(&nd)
						}
						copy(key[:], nd.Seq)
						ni.ID = nd.ID
						ni.EdgeIDIncoming = nd.EdgeIDIncoming
						ni.EdgeIDOutcoming = nd.EdgeIDOutcoming
						if v, ok := nodeMap.LoadOrStore(key, ni); ok {
							ni = v.(NodeInfo)
							for j := 0; j < BaseTypeNum; j++ {
								if nd.EdgeIDIncoming[j] == math.MaxUint32 {
									if ni.EdgeIDIncoming[j] == math.MaxUint32 {
										hasW = true
									} else {
										ni.EdgeIDIncoming[j] = math.MaxUint32
										nodeMap.Store(key, ni)
										ea[ec].EndNID = ni.ID
									}
									break
								}
								if nd.EdgeIDOutcoming[j] == math.MaxUint32 {
									if ni.EdgeIDOutcoming[j] == math.MaxUint32 {
										hasW = true
									} else {
										ni.EdgeIDOutcoming[j] = math.MaxUint32
										nodeMap.Store(key, ni)
										ea[ec].EndNID = ni.ID
									}
									break
								}
							}
						} else {
							id := <-nodeIDChan
							ni.ID = id
							id++
							nodeIDChan <- id
							ea[ec].EndNID = ni.ID
							nodeMap.Store(key, ni)
							nd.ID = ni.ID
							anc <- nd
						}
					}
					if !hasW {
						copy(key[:], node.Seq)
						v, _ := nodeMap.Load(key)
						ni = v.(NodeInfo)
						ni.EdgeIDOutcoming[i] = math.MaxUint32
						nodeMap.Store(key, ni)
						ea[ec].StartNID = node.ID
						ec++
						if ec == cap(ea) {
							wc <- ea
							ea = edgeArrPool.Get().([]DBGEdge)
							ea = ea[:cap(ea)]
							ec = 0
						}
					}
					mapSM.Unlock()
				}
			}
		}
	}
	fmt.Printf("[paraGenerateDBGEdges]edgeCount:%d tipNum:%d avgLen: %d\n", edgeCount, tipNum, tipTotalLen/tipNum)
}

// WritefqRecord write one record of fastq to the file
func WritefaRecord(fp io.Writer, ei *DBGEdge) {
	fmt.Fprintf(fp, ">%d\t%d\t%d\n", ei.ID, ei.StartNID, ei.EndNID)
	// fmt.Fprintf(edgesgzfp, "%s\n", ei.Ks)
	for i, b := range ei.Ks {
		if ei.Ks[i] > 3 || ei.Ks[i] < 0 {
			log.Fatalf("[WriteEdgesToFn] not correct base: %d:%b\n", i, b)
		}
		ei.Ks[i] = BitNtCharUp[b]
	}
	fmt.Fprintf(fp, "%s\n", ei.Ks)
	// fmt.Printf("[WriteEdgesToFn] edge: %v\n", ei)
}

// WriteEdgesToFn write edges seq to the file
func WriteEdgesToFn(edgesfn string, wc <-chan []DBGEdge, edgeArrPool *sync.Pool, numCPU int, kmerlen int) (edgeID uint32) {
	//oldNodeID := nodeID
	edgeID = uint32(2)
	edgesNum := 0
	edgesfp, err := os.Create(edgesfn)
	if err != nil {
		log.Fatalf("[WriteEdgesToFn] Open file %s failed: %v\n", edgesfn, err)
	}
	defer edgesfp.Close()
	fp, err1 := zstd.NewWriter(edgesfp, zstd.WithEncoderCRC(false), zstd.WithEncoderConcurrency(1), zstd.WithEncoderLevel(1))
	defer fp.Close()
	if err1 != nil {
		log.Fatalf("[WriteEdgesToFn]open write file:%s err: %v\n", edgesfn, err1)
	}

	finishNum := 0
	for {
		ea := <-wc
		if len(ea) == 0 {
			//fmt.Printf("[WriteEdgesToFn]edge:%v\n", ei)
			finishNum++
			if finishNum == numCPU {
				break
			}
			continue
		}

		for i := 0; i < len(ea); i++ {
			e := &ea[i]
			e.ID = edgeID
			edgeID++
			WritefaRecord(fp, e)
		}
		edgesNum += len(ea)
		edgeArrPool.Put(ea)
		//fmt.Printf("[WriteEdgesToFn] the writed edgeID: %v, ei.StartNID: %v, ei.EndNID: %v\n\tei.Ks: %v\n", edgeID-1, ei.StartNID, ei.EndNID, ei.Ks)
	}

	fmt.Printf("[WriteEdgesToFn] the writed file edges number is %d\n", edgesNum)
	//fmt.Printf("[WriteEdgesToFn] added nodes number is : %d\n", nodeID-oldNodeID)
	//nID = nodeID
	return
}

/*func ProcessAddedNode(cf CuckooFilter, nodeMap map[string]DBGNode, newNodeBntArr []ReadBnt, wc chan DBGEdge, nodeID uint32) (addedNodesNum, addedEdgesNum int) {

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
			var rb ReadBnt
			rb.Seq = node.Seq
			rb.Length = cf.Kmerlen - 1
			extRBnt := ExtendReadBnt2Byte(rb)
			for i := 0; i < BaseTypeNum; i++ {
				bi := byte(i)
				if node.EdgeIDIncoming[i] == 0 {
					var nBnt ReadBnt
					nSeq = append(nSeq, bi)
					nSeq = append(nSeq, extRSeq...)
					nLength = len(nSeq)
					ks := GetReadBntKmer(nBnt, 0, cf.Kmerlen)
					rs := ReverseComplet(ks)
					if ks.BiggerThan(rs) {
						ks, rs = rs, ks
					}
					count := cf.GetCountAllowZero(ks.Seq)
					if count >= MIN_KMER_COUNT {
						// get edge sequence
						edge, isNode := GetEdges(cf, nBnt, uint8(count), BACKWARD, MIN_KMER_COUNT)
						writedEdge := false
						if isNode == true {
							var tBnt ReadBnt
							tSeq = edge.Ks[:cf.Kmerlen-1]
							tks := GetReadBntKmer(tBnt, 0, cf.Kmerlen-1)
							trs := ReverseComplet(tks)
							sks := tks
							if sks.BiggerThan((trs)) {
								sks = trs
							}
							if v, ok := nodeMap[string(sks.Seq)]; ok {
								c := edge.Ks[cf.Kmerlen-1]
								if reflect.DeepEqual(sks, tks) {
									// b := Base2Bnt[c]
									v.EdgeIDOutcoming[c] = 1
								} else {
									// b := Base2Bnt[BitNtRev[c]]
									v.EdgeIDIncoming[BntRev[c]] = 1
								}
								edge.StartNID = v.ID
								writedEdge = true
								nodeMap[string(sks.Seq)] = v
							} else { // is a new node, add to the newNodeBntArr
								newNodeBntArr = append(newNodeBntArr, sks)
								processAdded++
							}
						} else { // is a tip
							if len(edge.Ks) > 2*cf.Kmerlen {
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
					var nBnt ReadBnt
					nSeq = append(nSeq, extRSeq...)
					nSeq = append(nSeq, bi)
					nLength = len(nSeq)
					ks := GetReadBntKmer(nBnt, 0, cf.Kmerlen)
					rs := ReverseComplet(ks)
					if ks.BiggerThan(rs) {
						ks, rs = rs, ks
					}
					count := cf.GetCountAllowZero(ks.Seq)
					if count >= MIN_KMER_COUNT {
						edge, isNode := GetEdges(cf, nBnt, uint8(count), FORWARD, MIN_KMER_COUNT)
						writedEdge := false
						if isNode == true {
							var tBnt ReadBnt
							tSeq = edge.Ks[len(edge.Ks)-cf.Kmerlen+1:]
							tks := GetReadBntKmer(tBnt, 0, cf.Kmerlen-1)
							trs := ReverseComplet(tks)
							sks := tks
							if sks.BiggerThan(trs) {
								sks = trs
							}
							if v, ok := nodeMap[string(sks.Seq)]; ok {
								c := edge.Ks[len(edge.Ks)-cf.Kmerlen]
								if reflect.DeepEqual(sks, tks) {
									// b := Base2Bnt[c]
									v.EdgeIDIncoming[c] = 1
								} else {
									// b := Base2Bnt[BitNtRev[c]]
									v.EdgeIDOutcoming[BntRev[c]] = 1
								}
								nodeMap[string(sks.Seq)] = v
								edge.EndNID = v.ID
								writedEdge = true
							} else {
								newNodeBntArr = append(newNodeBntArr, sks)
								processAdded++
							}
						} else { // is a tip
							if len(edge.Ks) > 2*cf.Kmerlen {
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
		for i := 0; i < BaseTypeNum; i++ {
			v.EdgeIDIncoming[i] = 0
			v.EdgeIDOutcoming[i] = 0
		}
		nodeMap[k] = v
	}
}*/

func GenerateDBGEdges(nodeMap *sync.Map, cf CuckooFilter, edgesfn string, numCPU int, nodeID uint32, opt optionsCDBG, NBntUint64Len int) (newNodeID uint32, edgeID uint32) {
	//MinKmerFreq := opt.MinKmerFreq
	TipMaxLen := opt.TipMaxLen
	bufsize := 100
	nc := make(chan DBGNode, numCPU)
	wc := make(chan []DBGEdge, numCPU)
	defer close(wc)
	edgeArrPool := sync.Pool{New: func() interface{} {
		arr := make([]DBGEdge, bufsize)
		for i := 0; i < len(arr); i++ {
			arr[i].Ks = make([]byte, opt.Kmer*5)
		}
		return arr
	}}
	anc := make(chan DBGNode, numCPU)
	defer close(anc)
	readNodeMapFinishedC := make(chan int, 1)
	nodeIDChan := make(chan uint32, 1)
	nodeIDChan <- nodeID
	// Read DBGNode to the nc
	go ReadDBGNodeToChan(nodeMap, nc, readNodeMapFinishedC, NBntUint64Len)
	// parallel construct edges from cuckoofilter
	for i := 0; i < numCPU; i++ {
		go paraGenerateDBGEdges(nc, cf, nodeMap, wc, &edgeArrPool, anc, nodeIDChan, TipMaxLen)
	}
	// collect added node and pass DBGNode to the nc
	//totalNodeNum := make(chan uint32，1)
	go CollectAddedDBGNode(anc, nodeMap, nc, readNodeMapFinishedC)
	// write edges Seq to the file
	edgeID = WriteEdgesToFn(edgesfn, wc, &edgeArrPool, numCPU, cf.Kmerlen)
	newNodeID = <-nodeIDChan
	// Change nodeMap monitor function
	//newNodeID, edgeID = ChangeNodeMap(nodeMap, anc, finishedC, nIEC, flagNIEC, cf.Kmerlen, nodeID)

	/* // cache the new added node Bnt info
	var newNodeBntArr []ReadBnt
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
	newNodeID = nodeID + uint32(addedNodesNum) */

	// clean set edgeID in the DBGNode
	//cleanEdgeIDInNodeMap(nodeMap)
	// fmt.Printf("[GenerateDBGEdges] added nodes number is : %d, added edges number is : %d\n", addedNodesNum, addedEdgesNum)
	return
}

func DBGStatWriter(DBGStatfn string, newNodeID, edgeID uint32) {
	DBGStatfp, err := os.Create(DBGStatfn)
	if err != nil {
		log.Fatalf("[DBGStatWriter] file %s create error: %v\n", DBGStatfn, err)
	}
	defer DBGStatfp.Close()
	fmt.Fprintf(DBGStatfp, "nodes size:\t%v\n", newNodeID)
	fmt.Fprintf(DBGStatfp, "edges size:\t%v\n", edgeID)
}

func DBGStatReader(DBGStatfn string) (nodesSize, edgesSize uint32) {
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

func NodesArrWriter(nodesArr []DBGNode, nodesfn string, fc chan<- int) error {
	nodesfp, err := os.Create(nodesfn)
	if err != nil {
		log.Fatalf("[NodesArrWriter] file %s create error, err: %v\n", nodesfn, err)
	}
	defer nodesfp.Close()
	//cbrofp := cbrotli.NewWriter(nodesfp, cbrotli.WriterOptions{Quality: 1})
	//defer cbrofp.Close()
	buffp := bufio.NewWriterSize(nodesfp, 1<<16) // 1<<25 == 2**25
	//var arr DBGNodeArr
	for i := 2; i < len(nodesArr); i++ {
		nd := nodesArr[i]
		if nd.GetDeleteFlag() > 0 || nd.ID < 2 {
			continue
		}
		//binary.Write(buffp, binary.LittleEndian, id)
		binary.Write(buffp, binary.LittleEndian, nd.ID)
		binary.Write(buffp, binary.LittleEndian, nd.EdgeIDIncoming)
		binary.Write(buffp, binary.LittleEndian, nd.EdgeIDOutcoming)
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
		var IDcoming1, IDcoming2 [BaseTypeNum]uint32
		for j := 0; j < BaseTypeNum; j++ {
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
		log.Fatalf("[NodesArrWriter] failed to flush file: %s, err: %v\n", nodesfn, err1)
	}
	fc <- 1
	return err
}

func NodesArrReader(nodesfn string, nodesArr []DBGNode, kmerlen int, fc chan<- int) {
	nodesfp, err := os.Open(nodesfn)
	if err != nil {
		log.Fatalf("[NodesArrReader] open file %s failed, err:%v\n", nodesfn, err)
	}
	defer nodesfp.Close()
	buffp := bufio.NewReaderSize(nodesfp, 1<<16)

	NSeqLen := (kmerlen - 1 + (NumBaseInUint64 - 1)) / NumBaseInUint64
	count := 0
	var nd DBGNode
	for {
		err = binary.Read(buffp, binary.LittleEndian, &nd.ID)
		if err == io.EOF {
			break
		} else if err != nil {
			log.Fatalf("[NodesArrReader]nd:%v\n\tnodesArr failed to read file:%s, err:%v\n", nd, nodesfn, err)
		}
		err = binary.Read(buffp, binary.LittleEndian, &nd.EdgeIDIncoming)
		if err != nil {
			log.Fatalf("[NodesArrReader]nd:%v\n\tnodesArr failed to read file:%s, err:%v\n", nd, nodesfn, err)
		}
		err = binary.Read(buffp, binary.LittleEndian, &nd.EdgeIDOutcoming)
		if err != nil {
			log.Fatalf("[NodesArrReader]nd:%v\n\tnodesArr failed to read file:%s, err:%v\n", nd, nodesfn, err)
		}
		nd.Seq = make([]uint64, NSeqLen)
		err = binary.Read(buffp, binary.LittleEndian, &nd.Seq)
		if err != nil {
			log.Fatalf("[NodesArrReader]nd:%v\n\tnodesArr failed to read file:%s, err:%v\n", nd, nodesfn, err)
		}
		count++
		//fmt.Printf("[NodesArrReader]nd:%v\n", nd)
		nodesArr[nd.ID] = nd
	}
	fmt.Printf("[NodesArrReader] recover nodes number : %v from file: %v\n", count, nodesfn)
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
		nd.ID = uint32(id)
		nd.EdgeIDIncoming[0], nd.EdgeIDIncoming[1], nd.EdgeIDIncoming[2], nd.EdgeIDIncoming[3] = uint32(ndi.EdgeIDIncoming[0]), uint32(ndi.EdgeIDIncoming[1]), uint32(ndi.EdgeIDIncoming[2]), uint32(ndi.EdgeIDIncoming[3])
		nd.EdgeIDOutcoming[0], nd.EdgeIDOutcoming[1], nd.EdgeIDOutcoming[2], nd.EdgeIDOutcoming[3] = uint32(ndi.EdgeIDOutcoming[0]), uint32(ndi.EdgeIDOutcoming[1]), uint32(ndi.EdgeIDOutcoming[2]), uint32(ndi.EdgeIDOutcoming[3])
		nd.Seq = ndi.Seq
		nodesArr[nd.ID] = nd
		//fmt.Printf("[NodesArrReader] nd: %v\n", nd)
	}*/

	/*json.Unmarshal(data, v)
	err4 = binary.Read(buffp, binary.LittleEndian)
	NSeqLen := (Kmerlen - 1 + (NumBaseInUint64 - 1)) / NumBaseInUint64
	count := 0
	var err1, err2, err3, err4, err5, err6 error
	for i := 2; i < len(nodesArr); i++ {
		var nd DBGNode
		nd.Seq = make([]uint64, NSeqLen)
		var id uint32
		err1 = binary.Read(buffp, binary.LittleEndian, id)
		nd.ID = uint32(id)
		var IDcoming1, IDcoming2 [BaseTypeNum]uint32
		err2 = binary.Read(buffp, binary.LittleEndian, &IDcoming1)
		err3 = binary.Read(buffp, binary.LittleEndian, &IDcoming2)
		for j := 0; j < BaseTypeNum; j++ {
			nd.EdgeIDIncoming[j] = uint32(IDcoming1[j])
			nd.EdgeIDOutcoming[j] = uint32(IDcoming2[j])
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

func NodeMap2NodeArr(nodeMap *sync.Map, nodesArr []DBGNode, NBntUint64Len int) {
	naLen := uint32(len(nodesArr))
	var nd DBGNode
	var ni NodeInfo
	var k [NODEMAP_KEY_LEN]uint64
	nodeMap.Range(func(key, value interface{}) bool {
		k = key.([NODEMAP_KEY_LEN]uint64)
		ni = value.(NodeInfo)
		if ni.ID >= naLen {
			log.Fatalf("[NodeMap2NodeArr] ni.ID:%d >= nodesArr len:%d\n", ni.ID, naLen)
		}
		if ni.ID > 1 && ni.ID != math.MaxUint32 {
			nd.Seq = make([]uint64, NBntUint64Len)
			for i := 0; i < NBntUint64Len; i++ {
				nd.Seq[i] = k[i]
			}
			nd.ID = ni.ID
			nd.EdgeIDIncoming = ni.EdgeIDIncoming
			nd.EdgeIDOutcoming = ni.EdgeIDOutcoming
			nodesArr[ni.ID] = nd
		}
		return true
	})
}
func writeComplexNodesToFile(complexNodesFn string, wc <-chan []DBGNode, nodeArrPool *sync.Pool, numCPU int) (complexNodeNum int) {
	ckfp, err := os.Create(complexNodesFn)
	if err != nil {
		log.Fatal(err)
	}
	defer ckfp.Close()
	buffp := bufio.NewWriterSize(ckfp, 1<<18)
	// if err != nil {
	// 	log.Fatal(err)
	// }
	//endFlagCount := 0

	// write complex node to the file
	// complexNodeNum := 0
	cpuCount := 0
	for {
		ndArr := <-wc
		if ndArr == nil {
			//fmt.Printf("[writeComplexNodesToFile]ndArr == nil\n")
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
		ni.EdgeIDIncoming = make([]uint32, BaseTypeNum)
		ni.EdgeIDOutcoming = make([]uint32, BaseTypeNum)
		for j := 0; j < BaseTypeNum; j++ {
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
		for _, nd := range ndArr {
			if err := binary.Write(buffp, binary.LittleEndian, nd.Seq); err != nil {
				log.Fatalf("[writeComplexNodesToFile] write node seq to file err: %v\n", err)
			}
			if err := binary.Write(buffp, binary.LittleEndian, nd.EdgeIDIncoming); err != nil {
				log.Fatalf("[writeComplexNodesToFile] write node seq to file err: %v\n", err)
			}
			if err := binary.Write(buffp, binary.LittleEndian, nd.EdgeIDOutcoming); err != nil {
				log.Fatalf("[writeComplexNodesToFile] write node seq to file err: %v\n", err)
			}
		}
		nodeArrPool.Put(ndArr)
		complexNodeNum += len(ndArr)
		// *** test code ***
		/* extRB := ExtendReadBnt2Byte(rb)
		fmt.Fprintf(os.Stderr, ">%v\n%v\n", complexNodeNum+1, Transform2Letters(extRB.Seq))
		*/
		// *** test code ***
		// ckgzfp.Write(rb.Seq)
		// ckgzfp.Write([]byte("\n"))
	}

	if err := buffp.Flush(); err != nil {
		log.Fatalf("[writeComplexNodesToFile] write to file: %s err: %v\n", complexNodesFn, err)
	}
	fmt.Printf("[writeComplexNodesToFile] found complex Node num is:%d\n", complexNodeNum)
	return
}

type optionsCDBG struct {
	ArgsOpt
	CFSize        int
	TipMaxLen     int
	MinKmerFreq   int
	MaxNGSReadLen int
}

func checkArgsCDBG(c cli.Command) (opt optionsCDBG, suc bool) {
	gOpt, suc := CheckGlobalArgs(c.Parent())
	if suc == false {
		log.Fatalf("[checkArgsCDBG] check global Arguments error, opt:%v\n", gOpt)
	}
	opt.ArgsOpt = gOpt
	opt.CFSize = c.Flag("S").Get().(int)
	if opt.CFSize < 1024*1024 {
		log.Fatalf("[checkArgsCCF]the argument'S':%d must bigger than 1024 * 1024\n", opt.CFSize)
	}
	opt.TipMaxLen = c.Flag("TipMaxLen").Get().(int)
	opt.MinKmerFreq = c.Flag("MinKmerFreq").Get().(int)
	opt.MaxNGSReadLen = c.Flag("MaxNGSReadLen").Get().(int)
	suc = true
	return
}

func ReconstructCuckooFilter(rbytesPool *sync.Pool, kmerArrPool *sync.Pool, kmerPoolChan <-chan KmerPool, cf *CuckooFilter) {
	for {
		kp, ok := <-kmerPoolChan
		if !ok {
			break
		}
		for i := range kp.KmerArr {
			_, suc := cf.Insert(kp.KmerArr[i])
			if suc == false {
				log.Fatal("[ReconstructCuckooFilter] Insert to the CuckooFilter false")
			}
		}
		rbytesPool.Put(kp.Cs)
		kmerArrPool.Put(kp.KmerArr)
	}
}

func CleanNodeArr(nodesArr []DBGNode) {
	for i, nd := range nodesArr {
		for j := 0; j < BaseTypeNum; j++ {
			if nd.EdgeIDIncoming[j] == 1 {
				nodesArr[i].EdgeIDIncoming[j] = 0
			}
			if nd.EdgeIDOutcoming[j] == 1 {
				nodesArr[i].EdgeIDOutcoming[j] = 0
			}
		}
	}
}

func CDBG(c cli.Command) {
	opt, suc := checkArgsCDBG(c)
	if suc == false {
		log.Fatalf("[CDBG] check Arguments error, opt:%v\n", opt)
	}
	fmt.Printf("[CDBG]opt:%v\n", opt)

	runtime.GOMAXPROCS(opt.NumCPU)
	//prefix := c.Parent().Flag("p").String()
	// create cpu profile
	/*profileFn := prefix + ".CDBG.prof"
	cpuprofilefp, err := os.Create(profileFn)
	if err != nil {
		log.Fatalf("[CDBG] open cpuprofile file: %v failed\n", profileFn)
	}
	pprof.StartCPUProfile(cpuprofilefp)
	defer pprof.StopCPUProfile()*/
	// cfinfofn := prefix + ".cfInfo"
	// cf, err := RecoverCuckooFilterInfo(cfinfofn)
	// if err != nil {
	// 	log.Fatalf("[CDGB] cuckoofilter recover err: %v\n", err)
	// }
	t0 := time.Now()
	//Reconstruct cuckoofilter
	//kifn := opt.Prefix + ".uniqkmerseq.info"
	//cfSize, err := ReadUniqKmerseqInfo(kifn)
	//if err != nil {
	//	log.Fatalf("[CDBG] Read CuckooFilter info file: %v err: %v\n", kifn, err)
	//}
	cf := MakeCuckooFilter(uint64(opt.CFSize), opt.Kmer)
	// read uniq kmers form file
	kmerfn := opt.Prefix + ".uniqkmerseq.zst"
	bufSize := WindowSize
	rbytesPool := sync.Pool{New: func() interface{} {
		bytes := make([]byte, bufSize)
		return bytes
	}}

	kmerArrPool := sync.Pool{New: func() interface{} {
		arr := make([][]byte, bufSize/opt.Kmer+2)
		return arr
	}}

	{
		kmerPoolChan := make(chan KmerPool, opt.NumCPU)
		nc := opt.NumCPU - 2
		if nc < 1 {
			nc = 1
		}
		for i := 0; i < nc; i++ {
			go ReconstructCuckooFilter(&rbytesPool, &kmerArrPool, kmerPoolChan, &cf)
		}
		GetUniqKmer(kmerfn, &rbytesPool, &kmerArrPool, kmerPoolChan, cf.Kmerlen)
		time.Sleep(time.Second)
	}
	fmt.Printf("[CDBG] Reconstruct CuckooFilter Struct used : %v\n", time.Now().Sub(t0))

	t0 = time.Now()
	cf.GetStat()
	// fmt.Printf("[CDBG] cf.Hash[0]: %v\n", cf.Hash[0])
	//Kmerlen = cf.Kmerlen

	//find complex Nodes
	// write complex Nodes to the file
	NBntUint64Len := (cf.Kmerlen - 1 + NumBaseInUint64 - 1) / NumBaseInUint64
	if NBntUint64Len > NODEMAP_KEY_LEN {
		log.Fatalf("[CDBG] nodeMap just allow max kmerlen is: %d, kmer set: %d\n", NODEMAP_KEY_LEN*32, cf.Kmerlen)
	}

	var complexNodeNum int
	complexKmerfn := opt.Prefix + ".complexNode"
	if fi, err := os.Stat(complexKmerfn); err != nil {
		wc := make(chan []DBGNode, opt.NumCPU)
		defer close(wc)
		kmerPoolChan := make(chan KmerPool, opt.NumCPU)
		nodeArrPool := sync.Pool{New: func() interface{} {
			arr := make([]DBGNode, 1000)
			return arr
		}}

		go GetUniqKmer(kmerfn, &rbytesPool, &kmerArrPool, kmerPoolChan, cf.Kmerlen)
		// identify complex Nodes
		for i := 0; i < opt.NumCPU-1; i++ {
			go paraLookupComplexNode(kmerPoolChan, &rbytesPool, &kmerArrPool, wc, &nodeArrPool, cf)
		}

		//complexNodeNum := writeComplexNodesToFile(complexKmerfn, wc, &nodeArrPool, opt.NumCPU-1)
		complexNodeNum = writeComplexNodesToFile(complexKmerfn, wc, &nodeArrPool, opt.NumCPU-1)
		fmt.Printf("[CDBG] search for complex DBG Node used : %v\n", time.Now().Sub(t0))
		t0 = time.Now()
	} else {
		complexNodeNum = int(fi.Size()) / (NBntUint64Len*8 + 8*4)
	}

	// construct Node map
	edgefn := opt.Prefix + ".edges.fa.zst"
	var nodesArr []DBGNode
	{
		//LenDBGNodeInfo := NBntUint64Len*8 + BaseTypeNum*2*4
		//nodeMap := make(map[[NODEMAP_KEY_LEN]uint64]DBGNode, complexNodeNum*2/5)
		var nodeMap sync.Map
		nodeID := constructNodeMap(complexKmerfn, &nodeMap, NBntUint64Len)
		runtime.GC()
		fmt.Printf("[CDBG]complexNodeNum:%d  assgin nodeID:%d\n", complexNodeNum, nodeID)
		// parallel generate edges and write to file
		//numCPU = 1
		newNodeID, edgeID := GenerateDBGEdges(&nodeMap, cf, edgefn, opt.NumCPU-1, nodeID, opt, NBntUint64Len)
		cf.Hash = nil
		runtime.GC()
		DBGStatfn := opt.Prefix + ".DBG.stat"
		DBGStatWriter(DBGStatfn, newNodeID, edgeID)
		// write DBG nodesArr to the file
		nodesArr = make([]DBGNode, newNodeID)
		NodeMap2NodeArr(&nodeMap, nodesArr, NBntUint64Len)
	}
	runtime.GC()
	// finish edgeIDComing set of DBGNode
	SetDBGNodeComing(edgefn, nodesArr, opt.Kmer)
	nodesfn := opt.Prefix + ".nodes.Arr"
	fc := make(chan int, 1)
	defer close(fc)
	CleanNodeArr(nodesArr)
	NodesArrWriter(nodesArr, nodesfn, fc)
	fmt.Printf("[CDBG] Generate DBG edges used : %v\n", time.Now().Sub(t0))
}

func ParseEdge(ri ReadInfo) (edge DBGEdge) {
	edge.ID = uint32(ri.ID)
	flist := strings.Fields(string(ri.Anotition))
	if len(flist) > 0 {
		id, err2 := strconv.Atoi(flist[0])
		if err2 != nil {
			log.Fatalf("[ParseEdge]eID:%d edge StartID: %v not digits, please convert to digits...\n", edge.ID, flist[0])
		}
		edge.StartNID = uint32(id)
	}
	if len(flist) > 1 {
		id, err2 := strconv.Atoi(flist[1])
		if err2 != nil {
			log.Fatalf("[ParseEdge]eID:%d edge EndID: %v not digits, please convert to digits...\n", edge.ID, flist[1])
		}
		edge.EndNID = uint32(id)
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
				ep[j].ID = uint32(eID)
			}
			edge.NGSPathArr[i] = ep
		}
		//fmt.Printf("[ParseEdge]NGSPathArr: %v\n", edge.NGSPathArr)
	}
	if len(flist) > 3 {

	}
	edge.Ks = make([]byte, len(ri.Seq))
	copy(edge.Ks, ri.Seq)
	return
}

func ReadEdgesFromFile(edgesfn string, edgesSize uint32, kmerlen int) (edgesArr []DBGEdge) {
	edgesArr = make([]DBGEdge, edgesSize)
	//format := GetReadsFileFormat(edgesbrfn)
	bufSize := WindowSize
	cs := make(chan []byte, 2)

	rbytesPool := sync.Pool{New: func() interface{} {
		bytes := make([]byte, bufSize)
		return bytes
	}}

	size := bufSize / (kmerlen*2 + 20)
	RIArrPool := sync.Pool{New: func() interface{} {
		arr := make([]ReadInfo, size)
		return arr
	}}

	RIPoolChan := make(chan RIPool, 2)

	go ReadZstdFile(edgesfn, &rbytesPool, cs)
	go GetReadInfoBucket(edgesfn, cs, &RIArrPool, RIPoolChan, false)

	var edgesNum int
	var e DBGEdge
	for {
		riPool, ok := <-RIPoolChan
		if !ok {
			break
		}
		for i := 0; i < len(riPool.RIArr); i++ {
			e = ParseEdge(riPool.RIArr[i])
			edgesArr[e.ID] = e
		}
		edgesNum += len(riPool.RIArr)
		RIArrPool.Put(riPool.RIArr)
		rbytesPool.Put(riPool.Cs)
		//fmt.Printf("[ReadEdgesFromFile]edge : %v, ks len:%d\n", edge, len(edge.Ks))
	}
	fmt.Printf("[ReadEdgesFromFile] found edge number is : %v\n", edgesNum)
	return
}

func SetDBGNodeComing(edgesfn string, nodesArr []DBGNode, kmerlen int) {
	//format := GetReadsFileFormat(edgesbrfn)
	bufSize := WindowSize
	cs := make(chan []byte, 2)

	rbytesPool := sync.Pool{New: func() interface{} {
		bytes := make([]byte, bufSize)
		return bytes
	}}

	size := bufSize / (kmerlen*2 + 20)
	RIArrPool := sync.Pool{New: func() interface{} {
		arr := make([]ReadInfo, size)
		return arr
	}}

	RIPoolChan := make(chan RIPool, 2)

	go ReadZstdFile(edgesfn, &rbytesPool, cs)
	go GetReadInfoBucket(edgesfn, cs, &RIArrPool, RIPoolChan, false)

	var edgesNum int
	var e DBGEdge
	rseq := make([]byte, kmerlen-1)
	for {
		riPool, ok := <-RIPoolChan
		if !ok {
			break
		}
		for i := 0; i < len(riPool.RIArr); i++ {
			e = ParseEdge(riPool.RIArr[i])
			//fmt.Printf("[SetDBGNodeComing]e.ID:%d\n", e.ID)
			// selfCycle edge
			//if e.StartNID > 1 && e.StartNID == e.EndNID {

			//}
			if e.StartNID > 1 {
				GetReverseCompletBytes(e.Ks[:kmerlen-1], rseq)
				if BiggerThan(e.Ks[:kmerlen-1], rseq) {
					if nodesArr[e.StartNID].EdgeIDIncoming[BntRev[e.Ks[kmerlen-1]]] != math.MaxUint32 {
						if e.StartNID == e.EndNID && nodesArr[e.StartNID].EdgeIDIncoming[BntRev[e.Ks[kmerlen-1]]] == e.ID {

						} else {
							log.Fatalf("[SetDBGNodeComing]nd:%v\n\te:%v\n\tnd.Seq:%v", nodesArr[e.StartNID], e, ExtendKmerBnt2Byte(nodesArr[e.StartNID].Seq, kmerlen-1))
						}
					}
					nodesArr[e.StartNID].EdgeIDIncoming[BntRev[e.Ks[kmerlen-1]]] = e.ID
				} else {
					if nodesArr[e.StartNID].EdgeIDOutcoming[e.Ks[kmerlen-1]] != math.MaxUint32 {
						if e.StartNID == e.EndNID && nodesArr[e.StartNID].EdgeIDOutcoming[e.Ks[kmerlen-1]] == e.ID {

						} else {
							log.Fatalf("[SetDBGNodeComing]nd:%v\n\te:%v\n\tnd.Seq:%v", nodesArr[e.StartNID], e, ExtendKmerBnt2Byte(nodesArr[e.StartNID].Seq, kmerlen-1))
						}
					}
					nodesArr[e.StartNID].EdgeIDOutcoming[e.Ks[kmerlen-1]] = e.ID
				}
			}

			if e.EndNID > 1 {
				nk := e.Ks[len(e.Ks)-(kmerlen-1):]
				GetReverseCompletBytes(nk, rseq)
				if BiggerThan(nk, rseq) {
					if nodesArr[e.EndNID].EdgeIDOutcoming[BntRev[e.Ks[len(e.Ks)-kmerlen]]] != math.MaxUint32 {
						if e.StartNID == e.EndNID && nodesArr[e.EndNID].EdgeIDOutcoming[BntRev[e.Ks[len(e.Ks)-kmerlen]]] == e.ID {

						} else {
							log.Fatalf("[SetDBGNodeComing]nd:%v\n\te:%v\n\tnd.Seq:%v", nodesArr[e.EndNID], e, ExtendKmerBnt2Byte(nodesArr[e.EndNID].Seq, kmerlen-1))
						}
					}
					nodesArr[e.EndNID].EdgeIDOutcoming[BntRev[e.Ks[len(e.Ks)-kmerlen]]] = e.ID
				} else {
					if nodesArr[e.EndNID].EdgeIDIncoming[e.Ks[len(e.Ks)-kmerlen]] != math.MaxUint32 {
						if e.StartNID == e.EndNID && nodesArr[e.EndNID].EdgeIDIncoming[e.Ks[len(e.Ks)-kmerlen]] == e.ID {

						} else {
							log.Fatalf("[SetDBGNodeComing]nd:%v\n\te:%v\n\tnd.Seq:%v", nodesArr[e.EndNID], e, ExtendKmerBnt2Byte(nodesArr[e.EndNID].Seq, kmerlen-1))
						}
					}
					nodesArr[e.EndNID].EdgeIDIncoming[e.Ks[len(e.Ks)-kmerlen]] = e.ID
				}
			}

		}
		rbytesPool.Put(riPool.Cs)
		RIArrPool.Put(riPool.RIArr)
		edgesNum += len(riPool.RIArr)
		//fmt.Printf("[ReadEdgesFromFile]edge : %v, ks len:%d\n", edge, len(edge.Ks))
	}

	fmt.Printf("[SetDBGNodeComing] found edge number is : %v\n", edgesNum)
}

func GetRCUnitig(u Unitig) (ru Unitig) {
	ru.Ks = GetReverseCompByteArr(u.Ks)
	ru.Kq = make([]uint8, len(u.Kq))
	copy(ru.Kq, u.Kq)
	ReverseUint8Arr(ru.Kq)
	return
}

func RCUnitig(u *Unitig) {
	ReverseCompByteArr(u.Ks)
	ReverseUint8Arr(u.Kq)
	return
}

/*func RevNode(node DBGNode, kmerlen int) DBGNode {
	rnode := node
	var nBnt KmerBnt
	nSeq = rnode.Seq
	nLen = kmerlen - 1
	rs := ReverseComplet(nBnt)
	rnode.Seq = rs.Seq
	for i := 0; i < BaseTypeNum; i++ {
		rnode.EdgeIDIncoming[i] = node.EdgeIDOutcoming[BntRev[i]]
		rnode.EdgeIDOutcoming[i] = node.EdgeIDIncoming[BntRev[i]]
		// rnode.EdgeIDOutcoming[BntRev[i]] = node.EdgeIDOutcoming[BaseTypeNum-1-i], node.EdgeIDIncoming[i]
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
				for i := 0; i < BaseTypeNum; i++ {
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
								var tBnt ReadBnt
								tSeq = edgesArr[eid].Ks[len(edgesArr[eid].Ks)-Kmerlen+1:]

								ks := GetReadBntKmer(tBnt, 0, Kmerlen-1)
								rs := ReverseComplet(ks)
								min := ks
								if min.BiggerThan(rs) {
									min = rs
								}
								if v2, ok := nodeMap[string(min.Seq)]; ok {
									var v2Bnt ReadBnt
									v2Seq = v2.Seq
									v2Length = Kmerlen - 1
									if v2.GetProcessFlag() == 1 {
										if reflect.DeepEqual(v2Bnt, ks) == false {
											// log.Fatalf("[ReconstructConsistenceDBG] found not consistence node\n")
											fmt.Printf("[ReconstructConsistenceDBG] found not consistence node, edge: %v\nv1: %v\nv2: %v\n", edgesArr[eid], node, v2)
											if eid == 2870 {
												fmt.Printf("[ReconstructConsistenceDBG] edge: %v\n", edgesArr[8014])
											}
											fmt.Printf("[ReconstructConsistenceDBG] edge start: %v\nedge end: %v\n", edgesArr[eid].Ks[:Kmerlen], edgesArr[eid].Ks[len(edgesArr[eid].Ks)-Kmerlen:])
											var v1Bnt ReadBnt
											v1Seq = node.Seq
											v1Length = Kmerlen - 1
											fmt.Printf("[ReconstructConsistenceDBG] v1.Seq: %v\n", ExtendReadBnt2Byte(v1Bnt))
											fmt.Printf("[ReconstructConsistenceDBG] v2.Seq: %v\n", ExtendReadBnt2Byte(v2Bnt))
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
								var tBnt ReadBnt
								tSeq = edgesArr[eid].Ks[:Kmerlen-1]
								ks := GetReadBntKmer(tBnt, 0, Kmerlen-1)
								rs := ReverseComplet(ks)
								min := ks
								if ks.BiggerThan(rs) {
									min = rs
								}

								if v2, ok := nodeMap[string(min.Seq)]; ok {
									var v2Bnt ReadBnt
									v2Seq = v2.Seq
									v2Length = Kmerlen - 1
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

/*func Coming2String(coming [BaseTypeNum]uint32) (cs string) {
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
		//labels := strconv.Itoa(int(e.ID)) + "len" + strconv.Itoa(len(e.Ks))
		//labels := strconv.Itoa(int(e.ID))
		labels := "\"ID:" + strconv.Itoa(int(e.ID)) + " len:" + strconv.Itoa(len(e.Ks)) + "\""
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

func IsContainCycleEdge(nd *DBGNode, arr []uint32) bool {
	arr = arr[:0]
	for i := 0; i < BaseTypeNum; i++ {
		if nd.EdgeIDIncoming[i] > 1 {
			arr = append(arr, nd.EdgeIDIncoming[i])
		}
		if nd.EdgeIDOutcoming[i] > 1 {
			arr = append(arr, nd.EdgeIDOutcoming[i])
		}
	}
	for i := 0; i < len(arr)-1; i++ {
		for j := i + 1; j < len(arr); j++ {
			if arr[i] == arr[j] {
				return true
			}
		}
	}

	return false
}

func GetEdgeIDComing(coming [BaseTypeNum]uint32) (num int, edgeID uint32) {
	for _, v := range coming {
		if v > 1 {
			num++
			edgeID = v
		}
	}

	return
}

func ConcatEdges(u1, u2 []byte, kmerlen int) (u []byte) {
	u = make([]byte, len(u1)+len(u2)-(kmerlen-1))
	if !reflect.DeepEqual(u1[len(u1)-kmerlen+1:], u2[:kmerlen-1]) {
		fmt.Fprintf(os.Stderr, "[ConcatEdges] u1:%s, u2:%s can not concatenatable\n", u1[len(u1)-kmerlen+1:], u2[:kmerlen-1])
		u = u[:0]
		return
	}

	copy(u[:len(u1)], u1)
	copy(u[len(u1):], u2[kmerlen-1:])
	return
}

/*func ConcatEdges(edgesArr []DBGEdge, inID, outID, dstID uint32) {
	// check if is connective
	var inBnt, outBnt ReadBnt
	fmt.Printf("[ConcatEdges] inID: %d, outID: %d, dstID: %d\nedgesArr[inID]:%v\nedgesArr[outID]:%v\n", inID, outID, dstID, edgesArr[inID], edgesArr[outID])
	inSeq = edgesArr[inID].Ks[len(edgesArr[inID].Ks)-Kmerlen+1:]
	inLength = len(inSeq)
	outSeq = edgesArr[outID].Ks[:Kmerlen-1]
	outLength = len(outSeq)
	//fmt.Printf("[ConcatEdges] inID: %d, outID: %d, dstID: %d\nseq1:%v\nseq2:%v\n", inID, outID, dstID, inBnt, outBnt)
	if reflect.DeepEqual(inSeq, outSeq) == false {
		log.Fatalf("[ConcatEdges] two edges is not connective\n\tin: %v\n\tout: %v\n", edgesArr[inID], edgesArr[outID])
	}
	if dstID == inID {
		u1, u2 := edgesArr[inID]. edgesArr[outID].Utg
		edgesArr[inID].Ks = append(u1.Ks, u2.Ks[Kmerlen-1:]...)
		for i := 0; i < Kmerlen-1; i++ {
			if u1.Kq[len(u1.Kq)-Kmerlen+1+i] < u2.Kq[i] {
				u1.Kq[len(u1.Kq)-Kmerlen+1+i] = u2.Kq[i]
			}
		}
		edgesArr[inID].Kq = append(u1.Kq, u2.Kq[Kmerlen-1:]...)
		// DeleteEdgeID(nodesArr, edgesArr[inID].EndNID, inID)
		edgesArr[inID].EndNID = edgesArr[outID].EndNID

	} else {
		u1, u2 := edgesArr[inID]. edgesArr[outID].Utg
		seq := make([]byte, len(u1.Ks))
		copy(seq, u1.Ks)
		edgesArr[outID].Ks = append(seq, u2.Ks[Kmerlen-1:]...)
		qul := make([]uint8, len(u1.Kq))
		copy(qul, u1.Kq)
		for i := 0; i < Kmerlen-1; i++ {
			if u1.Kq[len(u1.Kq)-Kmerlen+1+i] < u2.Kq[i] {
				qul[len(qul)-Kmerlen+1+i] = u2.Kq[i]
			}
		}
		edgesArr[outID].Kq = append(qul, u2.Kq[Kmerlen-1:]...)
		// DeleteEdgeID(nodesArr, edgesArr[outID].StartNID, outID)
		edgesArr[outID].StartNID = edgesArr[inID].StartNID
	}
}*/

/*func substituteEdgeID(nodeMap map[[NODEMAP_KEY_LEN]uint64]DBGNode, nodekey []uint64, srcID, dstID uint32, kmerlen int) bool {
	var nkB KmerBnt
	nkB.Seq = nodekey
	ks := GetReadBntKmer(nkB, 0, kmerlen-1)
	rs := ReverseComplet(ks)
	min := ks
	if ks.BiggerThan(rs) {
		min = rs
	}
	suc := false
	if nv, ok := nodeMap[string(min.Seq)]; ok {
		for i := 0; i < BaseTypeNum; i++ {
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

/*func ChangeIDMapPathAndJoinPathArr(IDMapPath map[uint32]uint32, joinPathArr []Path, e1, e2 DBGEdge, v DBGNode) []Path {
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
			tp.IDArr = GetReverseUint32Arr(p2.IDArr)
		}
		p.IDArr = append(p.IDArr, tp.IDArr...)
	} else {
		tp := p2
		if e2.StartNID == v.ID {
			tp.IDArr = GetReverseUint32Arr(p2.IDArr)
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

func CleanDBGNodeEdgeIDComing(nd *DBGNode) {
	for i := 0; i < BaseTypeNum; i++ {
		if nd.EdgeIDIncoming[i] == 1 {
			nd.EdgeIDIncoming[i] = 0
		}
		if nd.EdgeIDOutcoming[i] == 1 {
			nd.EdgeIDOutcoming[i] = 0
		}
		if nd.EdgeIDIncoming[i] == math.MaxUint32 && nd.EdgeIDOutcoming[i] == math.MaxUint32 {
			log.Fatalf("[CleanDBGNodeEdgeIDComing]nd:%v\n", nd)
		}
	}
}

func CleanDBGEdgeIDComing(nodesArr []DBGNode) {
	for i := 0; i < len(nodesArr); i++ {
		if i < 2 || nodesArr[i].GetDeleteFlag() > 0 {
			continue
		}
		CleanDBGNodeEdgeIDComing(&nodesArr[i])
	}
}

/*func MakeSelfCycleEdgeOutcomingToIncoming(nodesArr []DBGNode, edgesArr []DBGEdge, opt Options) {
	for i, e := range edgesArr {
		if i < 2 || e.GetDeleteFlag() > 0 {
			continue
		}
		if e.StartNID == e.EndNID { // self cycle edge
			nd := nodesArr[e.StartNID]
			bn := GetReadBntKmer(e.Ks, 0, opt.Kmer-1)
			if !reflect.DeepEqual(bn.Seq, nd.Seq) {
				e.Ks = GetReverseCompByteArr(e.Ks)
				ReverseByteArr(e.Kq)
			}
			if (nd.EdgeIDOutcoming[e.Ks[opt.Kmer-1]] != e.ID) || (nd.EdgeIDIncoming[e.Ks[len(e.Ks)-opt.Kmer]] != e.ID) {
				log.Fatalf("[MakeSelfCycleEdgeOutcomingToIncoming] error cycle edge set, e: %v\n\tv: %v\n", e, nd)
			}

		}
	}
}*/

func SubstituteEdgeID(nodesArr []DBGNode, nID, srcEID, dstEID uint32) bool {
	var success bool
	for i := 0; i < BaseTypeNum; i++ {
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

func SetEdgeID(nodesArr []DBGNode, edgesArr []DBGEdge) []DBGEdge {
	//nEA := make([]DBGEdge, 0, len(edgesArr)*2/3)
	edgeID := uint32(0)
	for i, e := range edgesArr {
		if i < 2 {
			edgeID++
			continue
		}
		if e.ID < 2 || e.GetDeleteFlag() > 0 {
			continue
		}
		if e.StartNID > 0 {
			if !SubstituteEdgeID(nodesArr, e.StartNID, e.ID, edgeID) {
				log.Fatalf("[SetEdgeID]v:%v\ne.ID:%d substitute by:%d failed\n", nodesArr[e.StartNID], e.ID, edgeID)
			}
		}
		if e.EndNID > 0 {
			if !SubstituteEdgeID(nodesArr, e.EndNID, e.ID, edgeID) {
				if e.StartNID != e.EndNID {
					log.Fatalf("[SetEdgeID]v:%v\ne.ID:%d substitute by:%d failed\n", nodesArr[e.EndNID], e.ID, edgeID)
				}
			}
		}
		e.ID = edgeID
		edgesArr[edgeID] = e
		edgeID++
	}
	fmt.Printf("[SetEdgeID]len(edgeArr):%d edgeID:%d\n", len(edgesArr), edgeID)
	edgesArr = edgesArr[:edgeID]
	return edgesArr
}

func CountNumEdgeIDComing(eIDComing [4]uint32) (count int) {
	for _, id := range eIDComing {
		if id > 0 {
			count++
		}
	}
	return
}

func SetNodeID(nodesArr []DBGNode, edgesArr []DBGEdge) []DBGNode {
	//nNA := make([]DBGNode, 0, len(nodesArr)*3/5)
	nodeID := uint32(0)
	for i, n := range nodesArr {
		if i < 2 {
			nodeID++
			continue
		}
		if n.ID < 2 || n.GetDeleteFlag() > 0 {
			continue
		}
		if CountNumEdgeIDComing(n.EdgeIDIncoming) == 0 && CountNumEdgeIDComing(n.EdgeIDOutcoming) == 0 {
			continue
		}
		for j := 0; j < BaseTypeNum; j++ {
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
						log.Fatalf("[SetNodeID]e.ID:%v StartNID:%v,EndNID:%v substitute by:%d failed\n", eID1, edgesArr[eID1].StartNID, edgesArr[eID1].EndNID, nodeID)
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
						log.Fatalf("[SetNodeID]e.ID:%v StartNID:%v,EndNID:%v substitute by:%d failed\n", eID2, edgesArr[eID2].StartNID, edgesArr[eID2].EndNID, nodeID)
					}
				}
			}
		}
		n.ID = nodeID
		nodesArr[nodeID] = n
		nodeID++
	}
	fmt.Printf("[SetNodeID]len(nodesArr):%d nodeID:%d\n", len(nodesArr), nodeID)
	nodesArr = nodesArr[:nodeID]
	return nodesArr
}

func GetEIDArr(coming [BaseTypeNum]uint32, arr []uint32) []uint32 {
	for _, id := range coming {
		if id > 1 {
			arr = append(arr, id)
		}
	}
	return arr
}

func GetNextEdgeIDArr(e *DBGEdge, strand bool, direction uint8, eIDArr []uint32) []uint32 {
	eIDArr = eIDArr[:0]
	if direction == FORWARD {
		if strand == PLUS {
			eIDArr = GetEIDArr(e.EdgeIDOutcoming, eIDArr)
		} else {
			eIDArr = GetEIDArr(e.EdgeIDIncoming, eIDArr)
		}
	} else {
		if strand == PLUS {
			eIDArr = GetEIDArr(e.EdgeIDIncoming, eIDArr)
		} else {
			eIDArr = GetEIDArr(e.EdgeIDOutcoming, eIDArr)
		}
	}
	return eIDArr
}

func IsInuint32Arr(arr []uint32, eID uint32) bool {
	for _, id := range arr {
		if id == eID {
			return true
		}
	}
	return false
}

func IsInuint32Slice(arr [BaseTypeNum]uint32, eID uint32) bool {
	for _, id := range arr {
		if id == eID {
			return true
		}
	}
	return false
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
	eIDArr := make([]uint32, 0, BaseTypeNum)
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
			eIDArr = GetNextEdgeIDArr(le, strand, BACKWARD, eIDArr)
			if !IsInuint32Arr(eIDArr, arr[i].ID) {
				return false
			}
			ne := &edgesArr[arr[i].ID]
			strand = GetNextEdgeStrand2(le, ne, strand, BACKWARD)
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
			eIDArr = GetNextEdgeIDArr(le, strand, FORWARD, eIDArr)
			if !IsInuint32Arr(eIDArr, arr[i].ID) {
				return false
			}
			ne := &edgesArr[arr[i].ID]
			strand = GetNextEdgeStrand2(le, ne, strand, FORWARD)
			le = ne
		}
	}

	return true
}

func GetReverseEdgeFreqArr2(arr []EdgeFreq) []EdgeFreq {
	la := len(arr)
	rArr := make([]EdgeFreq, la)
	for i, ef := range arr {
		rArr[la-1-i] = ef
	}
	return rArr
}

func IsMergedNode(nd *DBGNode) bool {
	var inNum, outNum int
	for i := 0; i < BaseTypeNum; i++ {
		if nd.EdgeIDIncoming[i] > 1 {
			inNum++
		}
		if nd.EdgeIDOutcoming[i] > 1 {
			outNum++
		}
	}
	if inNum == 1 && outNum == 1 {
		return true
	} else {
		return false
	}
}

func SmfyDBG(nodesArr []DBGNode, edgesArr []DBGEdge, opt optionsSF) {
	kmerlen := opt.Kmer
	deleteNodeNum, deleteEdgeNum := 0, 0
	longTipsEdgesNum := 0
	// clean DBG EdgeIDComing
	CleanDBGEdgeIDComing(nodesArr)

	// delete maybe short repeat edge small than opt.MaxNGSReadLen
	var e *DBGEdge
	for i := range edgesArr {
		e = &edgesArr[i]
		if e.ID < 2 {
			continue
		}

		if e.GetDeleteFlag() > 0 {
			//deleteEdgeNum++
			fmt.Printf("[SmfyDBG]has been deleted edge:%s", e.String())
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
		if e.StartNID == 0 && e.EndNID == 0 && len(e.Ks) < opt.TipMaxLen {
			deleteFlag = true
		} else if e.StartNID == 0 && e.EndNID > 0 && len(e.Ks) < opt.TipMaxLen {
			deleteFlag = true
			if !SubstituteEdgeID(nodesArr, e.EndNID, e.ID, 0) {
				fmt.Fprintf(os.Stderr, "[SmfyDBG]v: %v\ne.ID: %v substitute by 0 failed\n", nodesArr[e.EndNID], e.ID)
			}
			//fmt.Printf("[SmfyDBG]delete e.ID : %v, v: %v\n", e.ID, nodesArr[e.EndNID])
		} else if e.EndNID == 0 && e.StartNID > 0 && len(e.Ks) < opt.TipMaxLen {
			deleteFlag = true
			if !SubstituteEdgeID(nodesArr, e.StartNID, e.ID, 0) {
				fmt.Fprintf(os.Stderr, "[SmfyDBG]v: %v\ne.ID: %v substitute by 0 failed\n", nodesArr[e.StartNID], e.ID)
			}
		}

		if deleteFlag {
			edgesArr[i].SetDeleteFlag()
			deleteEdgeNum++
			fmt.Printf("[SmfyDBG]delete tip edge:%s", e.String())
			edgesArr[i].Ks = nil
			continue
		}
		if (e.StartNID == 0 || e.EndNID == 0) && len(e.Ks) >= opt.TipMaxLen {
			longTipsEdgesNum++
		}
	}

	// remove samll cycle maybe repeat
	/*if e.StartNID > 0 && e.StartNID == e.EndNID && len(e.Ks) < opt.MaxNGSReadLen {
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

	var e1, e2 *DBGEdge
	comingArr := make([]uint32, BaseTypeNum*2)
	for i, v := range nodesArr {
		if v.ID < 2 || v.GetDeleteFlag() > 0 || IsContainCycleEdge(&v, comingArr) {
			continue
		}
		inNum, inID := GetEdgeIDComing(v.EdgeIDIncoming)
		outNum, outID := GetEdgeIDComing(v.EdgeIDOutcoming)
		if inNum == 0 && outNum == 0 {
			nodesArr[i].SetDeleteFlag()
			deleteNodeNum++
		} else if inNum+outNum == 1 {
			id := inID
			if outNum == 1 {
				id = outID
			}
			if edgesArr[id].GetDeleteFlag() > 0 {
				nodesArr[i].SetDeleteFlag()
				fmt.Printf("[SmfyDBG]delete tip node ID:%d eID:%d\n", i, id)
				deleteNodeNum++
				continue
			}
			if edgesArr[id].StartNID == edgesArr[id].EndNID {
				continue
			}
			//fmt.Printf("[SmfyDBG]v: %v,id: %v\n", v, id)
			//fmt.Printf("[SmfyDBG]edgesArr[%v]: %v\n",id, edgesArr[id])
			//if edgesArr[id].StartNID == v.ID {
			//	edgesArr[id].StartNID = 0
			//} else {
			//	edgesArr[id].EndNID = 0
			//}
			if len(edgesArr[id].Ks) < opt.TipMaxLen {
				edgesArr[id].SetDeleteFlag()
				deleteEdgeNum++
				fmt.Printf("[SmfyDBG]delete tip edge ID:%d el:%d\n", id, edgesArr[id].GetSeqLen())
				if edgesArr[id].StartNID > 0 && int(edgesArr[id].StartNID) != i && !SubstituteEdgeID(nodesArr, edgesArr[id].StartNID, id, 0) {
					log.Fatalf("[SmfyDBG]v: %v\ne.ID: %v substitute by 0 failed\n", nodesArr[edgesArr[id].StartNID], id)
				} else if edgesArr[id].EndNID > 0 && int(edgesArr[id].EndNID) != i && !SubstituteEdgeID(nodesArr, edgesArr[id].EndNID, id, 0) {
					log.Fatalf("[SmfyDBG]v: %v\ne.ID: %v substitute by 0 failed\n", nodesArr[edgesArr[id].EndNID], id)
				}
				edgesArr[id].Ks = nil
			} else {
				if edgesArr[id].StartNID == v.ID {
					edgesArr[id].StartNID = uint32(0)
				} else {
					edgesArr[id].EndNID = uint32(0)
				}
				longTipsEdgesNum++
			}
			nodesArr[i].SetDeleteFlag()
			fmt.Printf("[SmfyDBG]delete tip node ID:%d eID:%d\n", i, id)
			deleteNodeNum++
		}
		//fmt.Printf("[SmfyDBG] e1: %v\n\te2: %v\n\tnd: %v\n", e1, e2, v)
		/*u1, u2 := e1. e2.Utg
		if e1.EndNID == v.ID {
			nID := e2.EndNID
			if e2.StartNID != v.ID {
				u2 = GetRCUnitig(u2)
				nID = e2.StartNID
			}
			//fmt.Printf("[SmfyDBG]v: %v\ne1.ID: %v, e1.StartNID: %v, e1.EndNID: %v, e2.ID:%v, e2.StartNID: %v, e2.EndNID: %v\n", v, e1.ID, e1.StartNID, e1.EndNID, e2.ID, e2.StartNID, e2.EndNID)
			edgesArr[inID].= ConcatEdges(u1, u2, kmerlen)
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
			edgesArr[inID].= ConcatEdges(u2, u1, kmerlen)
			edgesArr[inID].StartNID = nID
			if nID > 0 && !SubstituteEdgeID(nodesArr, nID, e2.ID, e1.ID) {
				log.Fatalf("[SmfyDBG]v: %v\ne2.ID: %v substitute by e1.ID: %v failed, node: %v\n", v, e2.ID, e1.ID, nodesArr[nID])
			}
		}*/
	}

	// get merge path
	mergePathMap := make(map[uint32][]EdgeFreq, 1000)
	for _, v := range nodesArr {
		if v.ID < 2 || v.GetDeleteFlag() > 0 || IsContainCycleEdge(&v, comingArr) {
			continue
		}
		fmt.Printf("[SmfyDBG]v:%s\n", v.String())
		inNum, inID := GetEdgeIDComing(v.EdgeIDIncoming)
		outNum, outID := GetEdgeIDComing(v.EdgeIDOutcoming)
		if inNum == 1 && outNum == 1 && inID != outID { // prevent cycle ring
			fmt.Printf("[SmfyDBG]inNum:%d outNum:%d inID:%d outID:%d\n", inNum, outNum, inID, outID)
			e1 = &edgesArr[inID]
			e2 = &edgesArr[outID]
			if e1.StartNID == e1.EndNID || e2.StartNID == e2.EndNID { // if  encounter cycle ring have same EdgeIDcoming
				continue
			}
			nodesArr[v.ID].SetDeleteFlag()
			deleteNodeNum++

			// set merge path
			{
				//fmt.Printf("[SmfyDBG]merge EdgePair e1:%s\n\te2:%s\n\tnd:%s\n", PrintEdgeInfo(e1), PrintEdgeInfo(e2), PrintNodeInfo(&v))
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
		}
	}

	// reset DBG
	var mergeNum int
	{
		var em DBGEdge
		for k, v := range mergePathMap {
			if k != v[0].ID {
				continue
			}
			var d bool
			for _, ef := range v {
				if edgesArr[ef.ID].GetDeleteFlag() > 0 {
					d = true
					break
				}
			}
			if d {
				continue
			}
			mergeNum++
			fmt.Printf("[SmfyDBG]merged path:%v\n", v)
			em = ConcatEdgePathUtg(v, edgesArr, nodesArr, kmerlen)
			if em.EndNID > 0 && !SubstituteEdgeID(nodesArr, em.EndNID, v[len(v)-1].ID, em.ID) {
				log.Fatalf("[SmfyDBG]v: %v\ne2.ID: %v substitute by em.ID: %v failed, node: %v\n", nodesArr[em.EndNID], v[len(v)-1].ID, em.ID, nodesArr[em.EndNID])
			}
			edgesArr[em.ID] = em
			for _, ef := range v[1:] {
				edgesArr[ef.ID].SetDeleteFlag()
				edgesArr[ef.ID].Ks = nil
				deleteEdgeNum++
			}
		}
	}

	/*// reset NGSPath
	changedNGSPathCount := 0
	rmp := make([]EdgeFreq, 40)
	idArr := make([]uint32, 5)
	for i := range edgesArr {
		e = &edgesArr[i]
		if i < 2 || e.GetDeleteFlag() > 0 || len(e.NGSPathArr) == 0 {
			continue
		}
		for j := range e.NGSPathArr {
			idArr = idArr[:0]
			for x := 0; x < len(e.NGSPathArr[j]); x++ {
				id := e.NGSPathArr[j][x].ID
				if IsInUint32Arr(idArr, id) {
					continue
				}
				if mp, ok := mergePathMap[id]; ok {
					rmp = GetReverseEdgeFreqArr(mp, rmp)
					idArr = append(idArr, mp[0].ID, rmp[0].ID)
					al := len(e.NGSPathArr[j])
					e.NGSPathArr[j] = RePlaceNGSPath2(e.NGSPathArr[j], mp, rmp)
					if al != len(e.NGSPathArr[j]) {
						x = 1
					}
					changedNGSPathCount++
				}
			}
		}
	}*/

	//fmt.Printf("[SmfyDBG]changed NGSPath number is : %d\n", changedNGSPathCount)
	fmt.Printf("[SmfyDBG]merged path number is : %d\n", mergeNum)
	fmt.Printf("[SmfyDBG]deleted nodes number is : %d\n", deleteNodeNum)
	fmt.Printf("[SmfyDBG]deleted edges number is : %d\n", deleteEdgeNum)
	fmt.Printf("[SmfyDBG]long tips number is : %d\n", longTipsEdgesNum)
}

// check code
func CheckNGSPathArr(edgesArr []DBGEdge, nodesArr []DBGNode) {
	for i := range edgesArr {
		e := &edgesArr[i]
		if e.ID < 2 || e.GetDeleteFlag() > 0 || len(e.NGSPathArr) == 0 {
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
}

func CountEdgeIDComing(eIDComing [4]uint32, eID uint32) (count int) {
	for _, id := range eIDComing {
		if id == eID {
			count++
		}
	}
	return
}

func CheckDBGSelfCycle(nodesArr []DBGNode, edgesArr []DBGEdge, kmerlen int) {
	var extNb []byte
	for _, e := range edgesArr {
		if e.ID < 2 || e.StartNID < 2 || e.EndNID < 2 || e.GetDeleteFlag() > 0 {
			continue
		}
		if e.StartNID == e.EndNID {
			//var nb KmerBnt
			extNb = ExtendKmerBnt2Byte(nodesArr[e.StartNID].Seq, kmerlen-1)
			if EqualByteArr(extNb, e.Ks[:kmerlen-1]) {

			} else {
				if CountEdgeIDComing(nodesArr[e.StartNID].EdgeIDIncoming, e.ID) == 2 {
					//fmt.Printf("[CheckDBGSelfCycle] self cycle edge two end in EdgeIDIncoming,len(e): %v, v: %v\n", len(e.Ks), nodesArr[e.StartNID])
					if EqualByteArr(extNb, e.Ks[len(e.Ks)-(kmerlen-1):]) {
						continue
					}
				}
				ReverseCompByteArr(e.Ks)
				edgesArr[e.ID] = e
				//fmt.Printf("[CheckDBGSelfCycle] ReverSeComp self cycle edge: %v\n", e.ID)
			}
		}
	}
}

func PrintTmpNodesArr(nodesArr []DBGNode, prefix string) {
	nodesfn := prefix + ".tmp.nodes"
	nodesfp, err := os.Create(nodesfn)
	if err != nil {
		log.Fatalf("[PrintTmpNodesArr] failed to create file: %s, err: %v\n", nodesfn, err)
	}
	defer nodesfp.Close()

	for _, v := range nodesArr {
		if v.ID < 2 || v.GetDeleteFlag() > 0 {
			continue
		}
		//s := fmt.Sprintf("ID:%v EdgeIncoming:%v EdgeOutcoming:%v\n", v.ID, v.EdgeIDIncoming, v.EdgeIDOutcoming)
		nodesfp.WriteString(v.String())
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
		s := fmt.Sprintf("ID:%v EdgeIncoming:%v EdgeOutcoming:%v\n", v.ID, v.EdgeIDIncoming, v.EdgeIDOutcoming)
		nodesfp.WriteString(s)
	}

	for _, e := range edgesArr {
		if e.ID < 2 || e.GetDeleteFlag() > 0 {
			continue
		}
		s := fmt.Sprintf("ID:%d len:%d StartNID:%d EndNID:%d EdgeIDIncoming:%v EdgeIDIncoming:%v seq:%v\n", e.ID, len(e.Ks), e.StartNID, e.EndNID, e.EdgeIDIncoming, e.EdgeIDOutcoming, e.Ks)
		edgesfp.WriteString(s)
	}

}

func Transform2Char(ks []byte) []byte {
	cs := make([]byte, len(ks))
	for i, b := range ks {
		//fmt.Printf("[Transform2Char] ks[%v]: %c\n", i, b)
		//ls = append(ls, alphabet.Letter(BitNtCharUp[b]))
		cs[i] = BitNtCharUp[b]
	}
	return cs
}

func Transform2Char2(ks []byte) []byte {
	//cs := make([]byte, len(ks))
	for i, b := range ks {
		//fmt.Printf("[Transform2Char] ks[%v]: %c\n", i, b)
		//ls = append(ls, alphabet.Letter(BitNtCharUp[b]))
		/*if b > 3 {
			b = 0
		}*/
		ks[i] = BitNtCharUp[b]
	}
	return ks
}

func Transform2Qual(kq []byte) []byte {
	//cs := make([]byte, len(ks))
	for i, q := range kq {
		//fmt.Printf("[Transform2Char] ks[%v]: %c\n", i, b)
		//ls = append(ls, alphabet.Letter(BitNtCharUp[b]))
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
		bs[i] = Base2Bnt[c]
	}
	return bs
}

var AdpaterSeq = "cggccgcaaggggttcgcgtcagcgggtgttggcgggtgtcggggctggcttaactatgcggcatcagagcagattgtactgagagtgcaccatatgcggtgtgaaataccacacagatgcgtaaggagaaaataccgcatcaggcgccattcgccattcagctgcgcaactgttgggaagggcgatcggtgcgggcctc"

// Set default quality(default = 1)
func SetDefaultQual(seq Unitig) (new Unitig) {
	for i := 0; i < len(seq.Ks); i++ {
		seq.Ks[i] = Base2Bnt[seq.Ks[i]]
		seq.Kq = append(seq.Kq, uint8(1))
	}

	return seq
}

/*func GetEdgesByteArr(edgesArr []DBGEdge, wc chan<- []byte) {
	for _, e := range edgesArr {
		if e.ID > 1 && e.GetDeleteFlag() == 0 && e.GetSeqLen() > 0 {
			eb := make([]byte, 0, e.GetSeqLen()*2+30+len(e.NGSPathArr)*30)
			eb = append(eb, '>')
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
			eb = append(eb, Transform2Char2(e.Ks)...)
			eb = append(eb, "\n"...)
			wc <- eb
		}
	}

	close(wc)
}*/

func StoreEdgesToFn(edgesfn string, edgesArr []DBGEdge, TipMaxLen int) {
	//wc := make(chan []byte, 2000)
	//go GetEdgesByteArr(edgesArr, wc)

	edgesfp, err := os.Create(edgesfn)
	if err != nil {
		log.Fatalf("[StoreEdgesToFn] create file:%s failed, err:%v\n", edgesfn, err)
	}
	defer edgesfp.Close()
	fp, err1 := zstd.NewWriter(edgesfp, zstd.WithEncoderCRC(false), zstd.WithEncoderConcurrency(1), zstd.WithEncoderLevel(1))
	defer fp.Close()
	if err1 != nil {
		log.Fatalf("[WriteEdgesToFn]open write file:%s err: %v\n", edgesfn, err1)
	}
	count := 0
	for _, e := range edgesArr {
		if e.ID < 2 || e.GetDeleteFlag() > 0 || (e.StartNID == 0 && e.EndNID == 0 && len(e.Ks) < TipMaxLen) {
			continue
		}
		WritefaRecord(fp, &e)
		count++
	}
	fmt.Printf("[StoreEdgesToFn] total write %d edges to file:%s\n", count, edgesfn)
}

func Set(ea1, ea2 []uint32) []uint32 {
	arr := make([]uint32, len(ea1))
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

func InterSecuint32Arr(ea1, ea2 []uint32) (arr []uint32) {
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

func IsEIDArrBubble(ea []uint32, edgesArr []DBGEdge) (ok bool) {
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
	Path       []uint32
	LastStrand bool
	Direction  uint8
	SeqLen     int
}

func GetNextPathArr(le *DBGEdge, nd DBGNode, MaxPathLen int, MaxBubbleSeqLen int, edgesArr []DBGEdge, nodesArr []DBGNode, kmerlen int) [][]uint32 {
	var pi PathInfo
	pi.Direction = FORWARD
	pi.Path = append(pi.Path, le.ID)
	if le.StartNID == nd.ID {
		pi.LastStrand = false
	} else {
		pi.LastStrand = true
	}
	var pathArr [][]uint32
	pi.SeqLen = kmerlen - 1
	stk := make([]PathInfo, 0, 20)
	stk = append(stk, pi)

	eIDArr := make([]uint32, 0, BaseTypeNum)
	for len(stk) > 0 {
		if len(stk) > 10 || len(pathArr) > 20 {
			pathArr = nil
			break
		}
		t := stk[len(stk)-1]
		stk = stk[:len(stk)-1]

		e := &edgesArr[t.Path[len(t.Path)-1]]
		eIDArr = GetNextEdgeIDArr(e, t.LastStrand, t.Direction, eIDArr)
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
			npi.LastStrand = GetNextEdgeStrand2(e, ne, t.LastStrand, t.Direction)
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
	//EPIPath []constructdbg.uint32
	Path []uint32
}

func IsInBubblePathArr(bubblePathArr []BubblePath, bp BubblePath) bool {
	for _, t := range bubblePathArr {
		if reflect.DeepEqual(t.Path, bp.Path) {
			return true
		}
	}
	return false
}

func SortPathArr(pathArr [][]uint32, edgesArr []DBGEdge, kmerlen int) [][]uint32 {
	spa := make([][]uint32, 0, len(pathArr))
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

func IsContainBubPathArr(p []uint32, bubPathArr []BubblePath) bool {
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
	Path []uint32
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

func GetBubblePathArr(pathArr [][]uint32, edgesArr []DBGEdge, kmerlen int, MaxPathLen, MaxBubbleSeqLen int) (bubPathArr []BubblePath, num int) {
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
			if AbsInt(cumLen[x]-cy) < Min(cumLen[x], cy)/10 {
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
					if AbsInt(cumLen[x]-cw) < Min(cumLen[x], cw)/10 {
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
					if AbsInt(cumLen[x]-cy) < Min(cumLen[x], cy)/10 {
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
							if AbsInt(cumLen[x]-cw) < Min(cumLen[x], cw)/10 {
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

func IsInComing(eIDcoming [4]uint32, eID uint32) bool {
	for _, id := range eIDcoming {
		if id == eID {
			return true
		}
	}
	return false
}

func IsBubbleEdge(e *DBGEdge, nodesArr []DBGNode) (ea []uint32, ok bool) {
	if (e.StartNID < 2) || (e.EndNID < 2) || e.StartNID == e.EndNID {
		return
	}
	var ea1, ea2 []uint32

	if IsInComing(nodesArr[e.StartNID].EdgeIDIncoming, e.ID) {
		ea1 = GetComingOtherEArr(nodesArr[e.StartNID].EdgeIDIncoming, e.ID)
	} else {
		ea1 = GetComingOtherEArr(nodesArr[e.StartNID].EdgeIDOutcoming, e.ID)
	}
	if IsInComing(nodesArr[e.EndNID].EdgeIDIncoming, e.ID) {
		ea2 = GetComingOtherEArr(nodesArr[e.EndNID].EdgeIDIncoming, e.ID)
	} else {
		ea2 = GetComingOtherEArr(nodesArr[e.EndNID].EdgeIDOutcoming, e.ID)
	}

	ea = InterSecuint32Arr(ea1, ea2)
	if len(ea) > 0 {
		ok = true
	}
	return
}

func GetOtherEArr(node DBGNode, eID uint32) (eIDArr []uint32) {
	var direction uint8
	if eID <= 0 {
		log.Fatalf("[GetOtherEArr] eID must bigger than zero, eID: %v\n", eID)
	}
	for i := 0; i < BaseTypeNum; i++ {
		if node.EdgeIDIncoming[i] == eID {
			direction = BACKWARD
			break
		}
		if node.EdgeIDOutcoming[i] == eID {
			direction = FORWARD
			break
		}
	}
	if direction != BACKWARD && direction != FORWARD {
		log.Fatalf("[GetOtherEArrEID] direction not set\n")
	}
	for i := 0; i < BaseTypeNum; i++ {
		if direction == FORWARD {
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

// set the unique edge of edgesArr
func SetDBGEdgesUniqueFlag(edgesArr []DBGEdge, nodesArr []DBGNode, kmerlen int) (uniqueNum, bubbleRepeatNum, twoCycleNum, selfCycle, selfCycleSameComingNum, bubbleEdgeNum int) {
	MaxPathLen := 5
	MaxBubbleSeqLen := 3 * kmerlen
	TipRegionNum := 6
	tipNumArr := make([]int, TipRegionNum)
	tipLenSumArr := make([]int, TipRegionNum)
	LenBoundArr := [5]int{500, 1000, 2000, 5000, 10000}
	var bubbleNum, bubbleTotalLen, SD, longTipNum, tipTotalLen, selfCycleTotalLen, uniqueTotalLen, repeatNum, longRepeatNum, repeatTotalLen, bubbleRepeatLen int
	ea1 := make([]uint32, 0, BaseTypeNum)
	ea2 := make([]uint32, 0, BaseTypeNum)
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

		ea1 = GetNextEdgeIDArr(e, PLUS, BACKWARD, ea1)
		ea2 = GetNextEdgeIDArr(e, PLUS, FORWARD, ea2)
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
			if ea, ok := IsBubbleEdge(e, nodesArr); ok {
				edgesArr[i].SetBubbleFlag()
				fmt.Fprintf(os.Stderr, "[SetDBGEdgesUniqueFlag]eID:%d len:%d Bubble\n", e.ID, el)
				bubbleEdgeNum++
				bubbleTotalLen += el
				//ea := GetOtherEArr(nodesArr[e.StartNID], e.ID)
				for _, id := range ea {
					if IsBubble(e, &edgesArr[id], nodesArr) {
						//fmt.Printf("[SetDBGEdgesUniqueFlag] Bubble edge: %v, length: %v\n", e.ID, len(e.Ks))
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
		if el > kmerlen*2 {
			longRepeatNum++
		} else {
			repeatNum++
		}
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
		fmt.Printf("[SetDBGEdgesUniqueFlag]repeat Number:%d longRepeatNum(>%d):%d average length:%d\n", repeatNum, 2*kmerlen, longRepeatNum, repeatTotalLen/(repeatNum+longRepeatNum))
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

func ResetComing(coming *[4]uint32) {
	coming[0], coming[1], coming[2], coming[3] = 0, 0, 0, 0
}

func AddNodeInfo2DBGEdgeArr(edgesArr []DBGEdge, nodesArr []DBGNode) {
	for i := range edgesArr {
		e := &edgesArr[i]
		ResetComing(&e.EdgeIDIncoming)
		ResetComing(&e.EdgeIDOutcoming)
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
	ID uint32
	Strand bool
}*/
type PathSeq struct {
	ID uint32
	//NID    uint32 // the path end node ID
	Strand bool
	//Start, End int
}

type ReadMapInfo struct {
	ID           int
	StartP, EndP int
	//NID          uint32
	//Seq        []byte
	PathSeqArr []PathSeq
	//Strands    []bool
}

type AlignInfo struct {
	ID      int64
	EndPos  int
	Paths   []uint32
	Strands []bool // strand for Paths
	Seq     []byte
}

/*func GetCycleEdgePath(pa []PathSeq) (ca []PathSeq) {
	for i := len(pa); i > 0; i-- {
		id := pa[i].ID
		for j := i - 1; j >= 0; j-- {
			if pa[j].ID == id {
				ca = pa[j : i+1]
				return
			}
		}
	}
	return
}

func IsInPathSeqArr(pa []PathSeq, id uint32) (ok bool) {
	for _, ps := range pa {
		if ps.ID == id {
			ok = true
			break
		}
	}
	return ok
}

func SumPathSeqLen(pa []PathSeq, edgesArr []DBGEdge, kmerlen int) (sl int) {
	for _, ps := range pa {
		sl += len(edgesArr[ps.ID].Ks) - (kmerlen - 1)
	}
	return sl
}*/

/*func paraLoadNGSReads(brfn string, cs chan ReadInfo, kmerLen int, we chan int) {
	var count int
	format := GetReadsFileFormat(brfn)
	bufSize := (1 << 26)
	size := 1000
	fncs1 := make(chan []byte, 3)
	recordChan1 := make(chan ReadInfo, size)

	go ReadBrFile(brfn, fncs1, bufSize)
	go GetReadFileRecord(fncs1, recordChan1, format, bufSize, true)

	for {
		ri, ok := <-recordChan1
		if !ok {
			break
		}
		cs <- ri
		count++
	}
	we <- count
}*/

/*func LoadNGSReadFiles(cfgFn string, correct bool, cs chan ReadInfo, numCPU, kmerLen int) {
	cfgInfo, err := ParseCfg(cfgFn, false, correct)
	if err != nil {
		log.Fatal("[ParseCfg] found err")
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
		if lib.AsmFlag != AllState && lib.SeqProfile != 1 {
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
}*/

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

func EqualByteArr(kb, rb []byte) bool {
	if len(kb) != len(rb) {
		log.Fatalf("[EqualByteArr] len(kb): %v != len(rb): %v\n", len(kb), len(rb))
	}
	for i, b := range kb {
		if b != rb[i] {
			return false
		}
	}
	return true
}

func GetCuckoofilterDBGSampleSize(edgesArr []DBGEdge, winSize, kmerlen int) (cfSize int) {
	var e *DBGEdge
	for i := range edgesArr {
		e = &edgesArr[i]
		if e.GetDeleteFlag() > 0 || e.ID < 2 {
			continue
		}
		//fmt.Printf("[GetCuckoofilterDBGSampleSize] e: %v\n", e)
		cfSize += ((len(e.Ks) - kmerlen + winSize - 1) / winSize) + 1
		/*if el < 2*maxNGSReadLen-kmerLen { // get whole length sliding windowns
			cfSize += ((el - kmerLen + winSize - 1) / winSize) + 1
		} else { // just sample two ends
			cfSize += ((maxNGSReadLen-kmerLen+winSize-1)/winSize + 1) * 2
		}*/
	}
	return cfSize
}

func GetMinimizer(seq []byte, kmerlen int) (minSeq []byte, pos int, strand bool) {
	var kb, rb []byte
	for i := 0; i < len(seq)-kmerlen+1; i++ {
		kb = seq[i : i+kmerlen]
		rb = GetReverseCompByteArr(kb)
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

func FindMinDBGKmerSeq(da []DBGKmerSeq) (min DBGKmerSeq, idx int) {
	min = da[0]
	idx = 0
	for i := 1; i < len(da); i++ {
		if BiggerThan(min.ks, da[i].ks) {
			min = da[i]
			idx = i
		}
	}
	return
}

type DBGKmerSeq struct {
	dk DBGKmer
	ks []byte
}

func ConstructCFDBGMinimizers(cf *CuckooFilterDBGKmer, edgesArr []DBGEdge, kmerlen, winSize int) {
	bufSize := winSize * 10
	buf := make([]DBGKmerSeq, bufSize)
	rs := make([]byte, 5000)
	var ds, drs DBGKmerSeq
	for _, e := range edgesArr {
		if e.ID < 2 || e.GetDeleteFlag() > 0 {
			continue
		}
		ks := e.Ks
		sl := len(ks)
		rs = GetReverseCompByteArr2(ks, rs)
		idx := 0
		last := -1

		for j := 0; j <= sl-kmerlen; j++ {
			ds.ks = ks[j : j+kmerlen]
			ds.dk.ID, ds.dk.Pos, ds.dk.Strand = e.ID, uint32(j), PLUS
			drs.ks = rs[sl-j-kmerlen : sl-j]
			drs.dk.ID, drs.dk.Pos, drs.dk.Strand = e.ID, uint32(j), MINUS
			if BiggerThan(ds.ks, drs.ks) {
				buf[idx] = drs
			} else {
				buf[idx] = ds
			}
			idx++

			if j < winSize-1 {
				continue
			}

			if last < 0 || idx-last > winSize {
				minDS, x := FindMinDBGKmerSeq(buf[idx-winSize : idx])
				y := idx - winSize + x
				//fmt.Printf("[constructCFDBGSample]dk:%v\n", minDS)
				if cf.Insert(minDS.ks, minDS.dk) == false {
					log.Fatalf("[ConstructCFDBGMinimizers] Insert to the CuckooFilter of DBGSample false, cf.Count:%d\n", cf.Count)
				}
				cf.Count++
				last = y
			} else {
				if BiggerThan(buf[last].ks, buf[idx-1].ks) {
					//fmt.Printf("[constructCFDBGSample]dk:%v\n", buf[idx-1])
					if cf.Insert(buf[idx-1].ks, buf[idx-1].dk) == false {
						log.Fatalf("[ConstructCFDBGMinimizers] Insert to the CuckooFilter of DBGSample false, cf.Count:%d\n", cf.Count)
					}
					cf.Count++
					last = idx - 1
				}
			}

			// adjust buf
			if idx == bufSize {
				copy(buf[:bufSize-last], buf[last:])
				idx = bufSize - last
				last = 0
			}
		}
	}
}

// found kmer seed position in the DBG edges
func LocateSeedKmerCF(cf CuckooFilterDBGKmer, ri ReadInfo, winSize, kmerlen int, edgesArr []DBGEdge, ka []DBGKmerSeq, rSeq []byte, dbgKArr []DBGKmer) ([]DBGKmer, DBGKmer) {
	//MaxStepNum := 10
	sl := len(ri.Seq)
	rSeq = GetReverseCompByteArr2(ri.Seq, rSeq)
	//ka := make([]KmerInfo, MaxStepNum*winSize)
	//fmt.Printf("[constructCFDBGSample] e: %v\n", e)
	idx := 0
	var ds, drs DBGKmerSeq
	for j := 0; j < winSize; j++ {
		ds.ks = ri.Seq[j : j+kmerlen]
		ds.dk.Pos, ds.dk.Strand = uint32(j), PLUS
		drs.ks = rSeq[sl-j-kmerlen : sl-j]
		drs.dk.Pos, drs.dk.Strand = uint32(j), MINUS
		if BiggerThan(drs.ks, ds.ks) {
			ka[idx] = ds
		} else {
			ka[idx] = drs
		}
		idx++
	}
	min, x := FindMinDBGKmerSeq(ka[:idx])
	dbgKArr = cf.Lookup(min.ks, dbgKArr)
	return dbgKArr, ka[x].dk
}

const (
	IN  = true
	OUT = false
)

func IsInDBGNode(nd DBGNode, eID uint32) bool {
	for i := 0; i < BaseTypeNum; i++ {
		if nd.EdgeIDIncoming[i] == eID || nd.EdgeIDOutcoming[i] == eID {
			return true
		}
	}

	return false
}

func GetNextMappingEID(nd DBGNode, e DBGEdge, base byte) uint32 {
	if e.StartNID == nd.ID {
		if IsInComing(nd.EdgeIDOutcoming, e.ID) {
			return nd.EdgeIDIncoming[base]
		} else {
			return nd.EdgeIDOutcoming[BntRev[base]]
		}
	} else {
		if IsInComing(nd.EdgeIDIncoming, e.ID) {
			return nd.EdgeIDOutcoming[base]
		} else {
			return nd.EdgeIDIncoming[BntRev[base]]
		}
	}
}

/*func GetNextMappingEdgeInfo(e DBGEdge, strand bool, direction uint8, ne DBGEdge, node DBGNode, kmerlen int) (bool, int) {
	var edgeNodeSeq [kmerlen-1]byte
	if direction == BACKWARD {
		if strand == PLUS {
			copy(edgeNodeSeq, e.Ks[:kmerlen-1])
		} else {

		}
	} else { // FORWARD

	}
}*/

func GetComingOtherEArr(coming [BaseTypeNum]uint32, eID uint32) (ea []uint32) {
	for _, id := range coming {
		if id > 1 && id != eID {
			ea = append(ea, id)
		}
	}
	return
}

/*func GetSelfCycleNextMapEdgeInfo(eID uint32, nd DBGNode, edgesArr []DBGEdge, nodeSeq []byte, kmerlen int, direction uint8, base byte, correct bool) (ID uint32, pos int, strand bool) {
	var tmp KmerBnt
	//kb.Seq, kb.Len = nodeSeq, len(nodeSeq)
	nodeBnt := GetReadBntKmer(nodeSeq, 0, kmerlen-1)
	revNdSeq := ReverseComplet(nodeBnt)
	tmp.Seq, tmp.Len = nd.Seq, kmerlen-1
	ndBntSeq := ExtendKmerBnt2Byte(tmp)
	if reflect.DeepEqual(nodeSeq, nd.Seq) {
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
			if ne.EndNID == nd.ID && reflect.DeepEqual(nodeSeq, ne.Ks[len(ne.Ks)-(kmerlen-1):]) {
				pos = len(ne.Ks) - (kmerlen - 1)
				strand = PLUS
				return
			}
			if ne.StartNID == nd.ID && reflect.DeepEqual(nodeSeq, GetReverseCompByteArr(ne.Ks[:kmerlen-1])) {
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
			if ne.StartNID == nd.ID && reflect.DeepEqual(nodeSeq, ne.Ks[:kmerlen-1]) {
				pos = kmerlen - 1
				strand = PLUS
				return
			}
			if ne.EndNID == nd.ID && reflect.DeepEqual(nodeSeq, GetReverseCompByteArr(ne.Ks[len(ne.Ks)-(kmerlen-1):])) {
				pos = len(ne.Ks) - (kmerlen - 1)
				strand = MINUS
			} else {
				log.Fatalf("[GetSelfCycleNextMapEdgeInfo] ne of node seq != nodeSeq, ne: %v\n\tnodeSeq: %v", ne, nodeSeq)
			}
		}
	} else if reflect.DeepEqual(revNdSeq.Seq, nd.Seq) {
		base = BntRev[base]
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
			if ne.StartNID == nd.ID && reflect.DeepEqual(nodeSeq, GetReverseCompByteArr(ne.Ks[:kmerlen-1])) {
				pos = kmerlen - 1
				strand = MINUS
				return
			}
			if ne.EndNID == nd.ID && reflect.DeepEqual(nodeSeq, ne.Ks[len(ne.Ks)-(kmerlen-1):]) {
				pos = len(ne.Ks) - (kmerlen - 1)
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
			if ne.EndNID == nd.ID && reflect.DeepEqual(nodeSeq, GetReverseCompByteArr(ne.Ks[len(ne.Ks)-(kmerlen-1):])) {
				pos = len(ne.Ks) - (kmerlen - 1)
				strand = MINUS
				return
			}
			if ne.StartNID == nd.ID && reflect.DeepEqual(nodeSeq, ne.Ks[:kmerlen-1]) {
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

func MappingReadToEdgesBackWard(dk DBGKmer, ri ReadInfo, rpos int, rstrand bool, edgesArr []DBGEdge, nodesArr []DBGNode, kmerlen int, correct bool) (errorNum int, ai AlignInfo) {
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
				if ri.Seq[i] != e.Ks[j] {
					errorNum++
					if !correct {
						break
					}
				}
				ai.Seq[i] = e.Ks[j]
				j--
			}
		} else { // strand == MINUS
			if len(e.Ks)-int(dk.Pos) < rpos {
				b = rpos - (len(e.Ks) - int(dk.Pos))
			}
			j = int(dk.Pos)
			for ; i >= b; i-- {
				//fmt.Printf("[MappingReadToEdgesBackWard] i: %v, j: %v, dk.Pos: %v, pos: %v, b: %v\n", i, j, dk.Pos, rpos, b)
				if ri.Seq[i] != BntRev[e.Ks[j]] {
					errorNum++
					if !correct {
						break
					}
				}
				ai.Seq[i] = BntRev[e.Ks[j]]
				j++
			}
		}

		if !correct && errorNum > 0 {
			fmt.Printf("[MappingReadToEdgesBackWard]not perfect start i: %v,edge ID: %v,len(e.Ks):%v,  dk.Pos: %v, pos: %v, b: %v\n", i, dk.ID, len(e.Ks), dk.Pos, rpos, b)
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
			base := BntRev[ri.Seq[i]]
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
				if reflect.DeepEqual(nodeSeq, ne.Ks[len(ne.Ks)-(kmerlen-1):]) {
					dk.Pos = len(ne.Ks) - (kmerlen - 1)
					strand = PLUS
				} else if reflect.DeepEqual(nodeSeq, GetReverseCompByteArr(ne.Ks[:kmerlen-1])) {
					dk.Pos = kmerlen - 1
					strand = MINUS
				}
			} else {
				if ne.EndNID == node.ID {
					dk.Pos = len(ne.Ks) - (kmerlen - 1)
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
			base := BntRev[ri.Seq[i]]
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
				if reflect.DeepEqual(nodeSeq, ne.Ks[len(ne.Ks)-(kmerlen-1):]) {
					dk.Pos = len(ne.Ks) - (kmerlen - 1)
					strand = PLUS
				} else if reflect.DeepEqual(nodeSeq, GetReverseCompByteArr(ne.Ks[:kmerlen-1])) {
					dk.Pos = kmerlen - 1
					strand = MINUS
				}
			} else {
				if ne.StartNID == node.ID {
					dk.Pos = kmerlen - 1
				} else {
					dk.Pos = len(ne.Ks) - (kmerlen - 1)
					strand = !strand
				}
			}
		}
	}
	return
}*/

/*func MappingReadToEdgesForWard(dk DBGKmer, ri ReadInfo, rpos int, rstrand bool, edgesArr []DBGEdge, nodesArr []DBGNode, Kmerlen int, correct bool) (errorNum int, ai AlignInfo) {
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
			if len(e.Ks)-int(dk.Pos) < len(ri.Seq)-rpos {
				b = rpos + (len(e.Ks) - int(dk.Pos))
			}
			j = int(dk.Pos)
			for ; i < b; i++ {
				//fmt.Printf("[MappingReadToEdgesForWard] i: %v, j: %v, dk.Pos: %v, pos: %v, b: %v\n", i, j, dk.Pos, rpos, b)
				if ri.Seq[i] != e.Ks[j] {
					errorNum++
					if !correct {
						break
					}
				}
				ai.Seq[i] = e.Ks[j]
				j++
			}
		} else { // strand == MINUS
			if len(ri.Seq)-rpos > int(dk.Pos) {
				b = rpos + int(dk.Pos)
			}
			j = int(dk.Pos) - 1
			for ; i < b; i++ {
				//fmt.Printf("[MappingReadToEdgesForWard] i: %v, j: %v, dk.Pos: %v, pos: %v, b: %v\n", i, j, dk.Pos, rpos, b)
				if ri.Seq[i] != BntRev[e.Ks[j]] {
					errorNum++
					if !correct {
						break
					}
				}
				ai.Seq[i] = BntRev[e.Ks[j]]
				j--
			}
		}

		if !correct && errorNum > 0 {
			fmt.Printf("[MappingReadToEdgesForWard]not perfect end i: %v,edge ID: %v,len(e.Ks): %v,  dk.Pos: %v, pos: %v, b: %v\n", i, dk.ID, len(e.Ks), dk.Pos, rpos, b)
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
			base = BntRev[ri.Seq[i]]
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
				if reflect.DeepEqual(nodeSeq, ne.Ks[:Kmerlen-1]) {
					dk.Pos = Kmerlen - 1
					strand = PLUS
				} else if reflect.DeepEqual(nodeSeq, GetReverseCompByteArr(ne.Ks[len(ne.Ks)-(Kmerlen-1):])) {
					dk.Pos = len(ne.Ks) - (Kmerlen - 1)
					strand = MINUS
				} else {
					log.Fatalf("[MappingReadToEdgesForWard] ne: %v not found proper start position\n", ne)
				}
			} else {
				if ne.StartNID == node.ID {
					dk.Pos = Kmerlen - 1
				} else {
					dk.Pos = len(ne.Ks) - (Kmerlen - 1)
					strand = !strand
				}
			}
		} else { // strand == MINUS
			if e.StartNID == 0 {
				break
			}
			node = nodesArr[e.StartNID]
			base = BntRev[ri.Seq[i]]
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
				if reflect.DeepEqual(nodeSeq, ne.Ks[:Kmerlen-1]) {
					dk.Pos = Kmerlen - 1
					strand = PLUS
				} else if reflect.DeepEqual(nodeSeq, GetReverseCompByteArr(ne.Ks[len(ne.Ks)-(Kmerlen-1):])) {
					dk.Pos = len(ne.Ks) - (Kmerlen - 1)
					strand = MINUS
				} else {
					log.Fatalf("[MappingReadToEdgesForWard] ne: %v not found proper start position\n", ne)
				}
			} else {
				if ne.EndNID == node.ID {
					dk.Pos = len(ne.Ks) - (Kmerlen - 1)
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
/*func paraMapNGS2DBG(cs chan ReadInfo, wc chan AlignInfo, nodesArr []DBGNode, edgesArr []DBGEdge, cf CuckooFilter, winSize int) {
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
			el := len(edgesArr[dbgK.ID].Ks)
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
			ar.Paths = Reverseuint32Arr(ar.Paths)
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
		ek := e.Ks[dbgK.Pos : dbgK.Pos+uint32(k)]
		if reflect.DeepEqual(ri.Seq[i:i+k], ek) {
			plus = true
		} else if reflect.DeepEqual(ri.Seq[i:i+k], GetReverseCompByteArr(ek)) {
			plus = false
		} else { // not equal for the DBG edge
			log.Fatalf("[paraMapNGS2DBG] read kmer not equal to the DBG edge\nread info: %v\nedge info: %v\n", ri.Seq[i:i+k], ek)
		}
		//if output {
		//	fmt.Printf("[paraMapNGS2DBG] i : %v\tdbgK: %v\n\tread info: %v\n\tedge seq: %v\n\tedge info: ID: %v, len: %v, startNID: %v, endNID: %v\n", i, dbgK, ri.Seq[i:i+k], ek, e.ID, len(e.Ks), e.StartNID, e.EndNID)
		//}

		x := i + k
		if plus {
			y := int(dbgK.Pos) + k
			for ; y < len(e.Ks) && x < len(ri.Seq); y++ {
				if ri.Seq[x] != e.Ks[y] {
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
			if y < len(e.Ks) {
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
				b := BntRev[ri.Seq[x]]
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
				b := BntRev[e.Ks[y]]
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
				b := BntRev[ri.Seq[x]]
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
			dbgK.Pos = uint32(len(e.Ks) - k)
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
	cs := make(chan ReadInfo, bufSize)
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

func AtoiArr(sa []string) []uint32 {
	da := make([]uint32, len(sa))
	for i, e := range sa {
		d, err := strconv.Atoi(e)
		if err != nil {
			log.Fatalf("[AtoiArr] string: %v convert to integer, err: %v\n", e, err)
		}
		da[i] = uint32(d)
	}
	return da
}

// coming denote node EdgeIDcoming, true is Outcoming, false is Incoming
func GetNearEdgeIDArr(nd DBGNode, eID uint32, coming bool) (eArr []uint32) {
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

/*func GetNextEdgeStrand(e, ne *DBGEdge, nodesArr []DBGNode, strand bool) bool {
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
}*/

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

func GetComingEdgeArr(edgeIDComing [BaseTypeNum]uint32) (ea []uint32) {
	for _, id := range edgeIDComing {
		if id > 1 {
			ea = append(ea, id)
		}
	}
	return
}

func FreqNumInuint32Arr(arr []uint32, eID uint32) (count int) {
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
			idx := IndexUint32(p.IDArr, e.ID)
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
			t := make([]uint32, al)
			idx := IndexUint32(p.IDArr, e.ID)
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
				var id uint32
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
				path.IDArr = GetReverseUint32Arr(path.IDArr)
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

func IsTwoEdgesCycle(e1, e2 *DBGEdge) bool {
	if e1.StartNID == e1.EndNID || e2.StartNID == e2.EndNID {
		return false
	}
	if e1.StartNID == e2.EndNID && e2.StartNID == e1.EndNID {
		return true
	}
	return false
}

func IsTwoEdgesCyclePath(edgesArr []DBGEdge, nodesArr []DBGNode, eID uint32) bool {
	e := edgesArr[eID]
	if e.StartNID == 0 || e.EndNID == 0 {
		return false
	}
	var arr1, arr2 []uint32
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
func GetOtherEdgeInTwoEdgesCycle(edgesArr []DBGEdge, nodesArr []DBGNode, eID uint32) uint32 {
	e := edgesArr[eID]
	if e.StartNID == 0 || e.EndNID == 0 {
		return 0
	}
	var arr1, arr2 []uint32
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

/* func ExtendPath(edgesArr []DBGEdge, nodesArr []DBGNode, e DBGEdge) (maxP Path) {
	maxP.IDArr = append(maxP.IDArr, e.PathMat[0].IDArr...)
	maxP.Freq = e.PathMat[0].Freq
	mutualArr := []uint32{e.ID}
	fmt.Printf("[ExtendPath] edge info, eID: %v\t e.StartNID: %v\te.EndNID: %v\n", e.ID, e.StartNID, e.EndNID)
	if e.StartNID > 0 {
		nd := nodesArr[e.StartNID]
		ea := GetNearEdgeIDArr(nd, e.ID)
		fmt.Printf("[ExtendPath] ea: %v, maxP: %v, nd: %v\n", ea, maxP, nd)
		if len(ea) == 1 {
			maxP.IDArr = Reverseuint32Arr(maxP.IDArr)
			//p1 := IndexUint32(maxP.IDArr, ea[0])
			p1 := IndexUint32(maxP.IDArr, e.ID) + 1
			if p1 > 0 {
				for j := p1; j < len(maxP.IDArr); j++ {
					ne := edgesArr[maxP.IDArr[j]]
					var nextNID uint32
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
						nP.IDArr = GetReverseUint32Arr(ne.PathMat[0].IDArr)
					} else {
						nP.IDArr = make([]uint32, len(ne.PathMat[0].IDArr))
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
			maxP.IDArr = Reverseuint32Arr(maxP.IDArr)
			mutualArr = Reverseuint32Arr(mutualArr)
		}
	}

	fmt.Printf("[ExtendPath] after extend previous edges, maxP: %v\n\te: %v\n", maxP, e)
	if e.EndNID > 0 {
		nd := nodesArr[e.EndNID]
		ea := GetNearEdgeIDArr(nd, e.ID)
		fmt.Printf("[ExtendPath] ea: %v, nd: %v\n", ea, nd)
		if len(ea) == 1 {
			//p1 := IndexUint32(maxP.IDArr, ea[0])
			p1 := IndexUint32(maxP.IDArr, e.ID) + 1
			if p1 > 0 {
				for j := p1; j < len(maxP.IDArr); j++ {
					ne := edgesArr[maxP.IDArr[j]]
					var nextNID uint32
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
						nP.IDArr = GetReverseUint32Arr(ne.PathMat[0].IDArr)
					} else {
						nP.IDArr = make([]uint32, len(ne.PathMat[0].IDArr))
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
		p1 := IndexUint32(maxP.IDArr, eID1)
		eID2 := mutualArr[len(mutualArr)-1]
		p2 := IndexUint32(maxP.IDArr, eID2)
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
/*func findMaxPath(edgesArr []DBGEdge, nodesArr []DBGNode, minMapFreq int, semi bool) (pathArr []Path) {
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
}*/

func GetLinkNodeID(p0, p1 DBGEdge) uint32 {
	if p0.StartNID == p1.StartNID || p0.StartNID == p1.EndNID {
		return p0.StartNID
	} else if p0.EndNID == p1.StartNID || p0.EndNID == p1.EndNID {
		return p0.EndNID
	} else {
		log.Fatalf("[GetLinkNodeID]not found link node ID p0: %v\n\tp1: %v\n", p0, p1)
	}
	return 0
}

func GetLinkPathArr(nodesArr []DBGNode, edgesArr []DBGEdge, nID uint32) (p Path) {
	nd := nodesArr[nID]
	inNum, inID := GetEdgeIDComing(nd.EdgeIDIncoming)
	outNum, outID := GetEdgeIDComing(nd.EdgeIDOutcoming)
	if inNum != 1 || outNum != 1 {
		log.Fatalf("[GetLinkPathArr] inNum: %v or outNum: %v != 1\n", inNum, outNum)
	}
	nodesArr[nID].SetProcessFlag()
	arr := [2]uint32{inID, outID}
	var dpArr [2][]uint32
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
	p.IDArr = Reverseuint32Arr(dpArr[0])
	p.IDArr = append(p.IDArr, dpArr[1]...)

	return p
}

func CascadePath(p Path, edgesArr []DBGEdge, nodesArr []DBGNode, kmerlen int, changeDBG bool) DBGEdge {

	// keep self cycle path start node OUtcoming, end node Incoming
	e1, e2 := edgesArr[p.IDArr[0]], edgesArr[p.IDArr[1]]
	if e1.EndNID == e2.StartNID || e1.EndNID == e2.EndNID {
		if IsInComing(nodesArr[e1.StartNID].EdgeIDIncoming, p.IDArr[0]) {
			Reverseuint32Arr(p.IDArr)
		}
	} else {
		if IsInComing(nodesArr[e1.EndNID].EdgeIDIncoming, p.IDArr[0]) {
			Reverseuint32Arr(p.IDArr)
		}
	}

	p0 := edgesArr[p.IDArr[0]]
	p1 := edgesArr[p.IDArr[1]]
	//eStartNID, eEndNID := p0.StartNID, p0.EndNID
	var direction uint8
	nID := GetLinkNodeID(p0, p1)
	if p0.StartNID == nID {
		ReverseCompByteArr(edgesArr[p.IDArr[0]].Ks)
		//ReverseByteArr(edgesArr[p.IDArr[0]].Kq)
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
				ReverseCompByteArr(p1.Ks)
			}
			fmt.Printf("[CascadePath]p0.ID: %v, p1.ID: %v, lastEID: %v, strand: %v, nID: Incoming: %v, Outcoming: %v\n", p0.ID, p1.ID, lastEID, strand, nodesArr[nID].EdgeIDIncoming, nodesArr[nID].EdgeIDOutcoming)
			edgesArr[p0.ID].Ks = ConcatEdges(p1.Ks, p0.Ks, kmerlen)
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
				ReverseCompByteArr(p1.Ks)
			}
			fmt.Printf("[CascadePath]FORWARD p0.ID: %v, p1.ID: %v, lastEID: %v, strand: %v, nID: Incoming: %v, Outcoming: %v\n", p0.ID, p1.ID, lastEID, strand, nodesArr[nID].EdgeIDIncoming, nodesArr[nID].EdgeIDOutcoming)
			edgesArr[p0.ID].Ks = ConcatEdges(p0.Ks, p1.Ks, kmerlen)
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

/*func ResetMergedEdgeNID(edgesArr []DBGEdge, nodesArr []DBGNode, joinPathArr []Path, IDMapPath map[uint32]uint32, kmerlen int) {
	for _, e := range edgesArr {
		if e.ID == 0 {
			continue
		}
		if idx, ok := IDMapPath[e.ID]; ok {
			if e.ID == joinPathArr[idx].IDArr[0] {
				if e.StartNID > 0 {
					//SubstituteEdgeID(nodesArr, e.StartNID, e.ID, 0)
					if !IsInDBGNode(nodesArr[e.StartNID], e.ID) {
						//var sBnt KmerBnt
						//sSeq = e.Ks[:kmerlen-1]
						ks := GetReadBntKmer(e.Ks[:kmerlen-1], 0, kmerlen-1)
						ts := ks
						rs := ReverseComplet(ks)
						if ks.BiggerThan(rs) {
							ks, rs = rs, ks
						}
						if reflect.DeepEqual(ks.Seq, nodesArr[e.StartNID].Seq) == false {
							log.Fatalf("[ResetMergedEdgeNID] ks != nodesArr[e.StartNID].Seq, ks: %v\n\tnodesArr[%v]: %v\n", ks, e.StartNID, nodesArr[e.StartNID])
						}
						c := e.Ks[kmerlen-1]
						if reflect.DeepEqual(ts, ks) {
							nodesArr[e.StartNID].EdgeIDOutcoming[c] = e.ID
						} else {
							nodesArr[e.StartNID].EdgeIDIncoming[BntRev[c]] = e.ID
						}
					}
				}
				if e.EndNID > 0 {
					//SubstituteEdgeID(nodesArr, e.EndNID, e.ID, 0)
					if !IsInDBGNode(nodesArr[e.EndNID], e.ID) {
						//var sBnt ReadBnt
						//sSeq = e.Ks[e.GetSeqLen()-kmerlen+1:]
						ks := GetReadBntKmer(e.Ks[e.GetSeqLen()-kmerlen+1:], 0, kmerlen-1)
						ts := ks
						rs := ReverseComplet(ks)
						if ks.BiggerThan(rs) {
							ks, rs = rs, ks
						}
						if reflect.DeepEqual(ks.Seq, nodesArr[e.EndNID].Seq) == false {
							log.Fatalf("[ResetMergedEdgeNID] ks != nodesArr[e.StartNID].Seq, ks: %v\n\tnodesArr[%v]: %v\n", ks, e.StartNID, nodesArr[e.StartNID])
						}
						c := e.Ks[e.GetSeqLen()-kmerlen]
						if reflect.DeepEqual(ks, ts) {
							nodesArr[e.EndNID].EdgeIDIncoming[c] = e.ID
						} else {
							nodesArr[e.EndNID].EdgeIDOutcoming[BntRev[c]] = e.ID
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
			nodesArr[i].EdgeIDIncoming = [4]uint32{0, 0, 0, 0}
			nodesArr[i].EdgeIDOutcoming = [4]uint32{0, 0, 0, 0}
			nodesArr[i].Seq = nil
			continue
		}
		for j := 0; j < BaseTypeNum; j++ {
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
	for i := range nodesArr {
		nodesArr[i].ResetProcessFlag()
	}
}

func ResetIDMapPath(IDMapPath map[uint32]uint32, joinPathArr []Path, path Path) []Path {
	var np Path
	np.Freq = path.Freq
	joinIdx := -1
	for _, t := range path.IDArr {
		if idx, ok := IDMapPath[t]; ok {
			arr := joinPathArr[idx].IDArr
			if t == arr[len(arr)-1] {
				Reverseuint32Arr(arr)
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
			fmt.Printf("[CleanEdgesArr]delete edge len: %v\t%v\n", len(edgesArr[i].Ks), edgesArr[i])
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

/*func CheckPathDirection(edgesArr []DBGEdge, nodesArr []DBGNode, eID uint32) {
	e := edgesArr[eID]
	step := 1
	var ea1, ea2 []uint32
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
		idx := IndexUint32(p.IDArr, eID)
		if idx < len(p.IDArr)-step {
			if IsInuint32Arr(ea1, p.IDArr[idx+step]) {
				Reverseuint32Arr(edgesArr[eID].PathMat[i].IDArr)
			}
		} else {
			if idx-step >= 0 && IsInuint32Arr(ea2, p.IDArr[idx-step]) {
				Reverseuint32Arr(edgesArr[eID].PathMat[i].IDArr)
			}
		}
	}

	//fmt.Printf("[CheckPathDirection] ea1: %v, ea2: %v, PathMat: %v\n", ea1, ea2, edgesArr[eID].PathMat)

}*/

/*func AdjustPathMat(edgesArr []DBGEdge, nodesArr []DBGNode, joinPathArr []Path, IDMapPath map[uint32]uint32) {

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
		idx := IndexUint32(p.IDArr, e.ID)
		arr := []uint32{e.ID}
		if e.StartNID > 0 && idx > 0 {
			nd := nodesArr[e.StartNID]
			te := e
			rp := GetReverseUint32Arr(p.IDArr[:idx])
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
					if IsInuint32Arr(ea, eID) {
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
			Reverseuint32Arr(arr)
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
					if IsInuint32Arr(ea, eID) {
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

// constuct map edge ID to the path
func ConstructIDMapPath(joinPathArr []Path) map[uint32]uint32 {
	IDMapPath := make(map[uint32]uint32)
	for i, p := range joinPathArr {
		sc := [2]uint32{p.IDArr[0], p.IDArr[len(p.IDArr)-1]}
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
	var e *DBGEdge
	for i := range edgesArr {
		e = &edgesArr[i]
		if e.ID < 2 || e.GetDeleteFlag() > 0 {
			continue
		}
		if e.StartNID > 0 && e.StartNID == e.EndNID && nodesArr[e.StartNID].GetDeleteFlag() > 0 {
			if e.GetSeqLen() < 1000 {
				e.SetDeleteFlag()
			} else {
				e.StartNID = uint32(0)
				e.EndNID = uint32(0)
			}
			continue
		}
		if e.StartNID > 0 {
			if nodesArr[e.StartNID].GetDeleteFlag() > 0 {
				e.StartNID = 0
			} else if !IsInDBGNode(nodesArr[e.StartNID], e.ID) {
				log.Fatalf("[CheckInterConnectivity]edge check nd: %v\n\tedge: %v\n", nodesArr[e.StartNID], e)
			}
		}
		if e.EndNID > 0 {
			if nodesArr[e.EndNID].GetDeleteFlag() > 0 {
				e.EndNID = 0
			} else if !IsInDBGNode(nodesArr[e.EndNID], e.ID) {
				log.Fatalf("[CheckInterConnectivity]edge check nd: %v\n\tedge: %v\n", nodesArr[e.EndNID], e)
			}
		}
	}

	for _, nd := range nodesArr {
		if nd.ID == 0 || nd.GetDeleteFlag() > 0 {
			continue
		}

		for j := 0; j < BaseTypeNum; j++ {
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

type optionsSF struct {
	ArgsOpt
	TipMaxLen     int
	WinSize       int
	MaxNGSReadLen int
	MinMapFreq    int
	Correct       bool
	ReNameID      bool
	Graph         bool
	//MaxMapEdgeLen int // max length of edge that don't need cut two flank sequence to map Long Reads
}

func checkArgsSF(c cli.Command) (opt optionsSF, succ bool) {
	gOpt, suc := CheckGlobalArgs(c.Parent())
	if suc == false {
		log.Fatalf("[checkArgsSF] check global Arguments error, opt: %v\n", gOpt)
	}
	opt.ArgsOpt = gOpt
	var ok bool
	opt.TipMaxLen, ok = c.Flag("TipMaxLen").Get().(int)
	if !ok {
		log.Fatalf("[checkArgsSF] argument 'TipMaxLen': %v set error\n ", c.Flag("TipMaxlen").String())
	}
	//opt.TipMaxLen = tmp

	//tmp, err = strconv.Atoi(c.Flag("WinSize").String())
	//tmp = c.Flag("WinSize")
	opt.WinSize, ok = c.Flag("WinSize").Get().(int)
	if !ok {
		log.Fatalf("[checkArgsSF] argument 'WinSize': %v set error\n ", c.Flag("WinSize").String())
	}
	if opt.WinSize < 1 || opt.WinSize > 100 {
		log.Fatalf("[checkArgsSF] argument 'WinSize': %v must between 1~100\n", c.Flag("WinSize"))
	}
	opt.MaxNGSReadLen, ok = c.Flag("MaxNGSReadLen").Get().(int)
	if !ok {
		log.Fatalf("[checkArgsSF] argument 'MaxNGSReadLen': %v set error\n ", c.Flag("MaxNGSReadLen").String())
	}
	if opt.MaxNGSReadLen < opt.Kmer+50 {
		log.Fatalf("[checkArgsSF] argument 'MaxNGSReadLen': %v must bigger than K+50\n", c.Flag("MaxNGSReadLen").String())
	}

	opt.MinMapFreq, ok = c.Flag("MinMapFreq").Get().(int)
	if !ok {
		log.Fatalf("[checkArgsSF] argument 'MinMapFreq': %v set error\n ", c.Flag("MinMapFreq").String())
	}
	if opt.MinMapFreq < 5 && opt.MinMapFreq >= 20 {
		log.Fatalf("[checkArgsSF] argument 'MinMapFreq': %v must 5 <= MinMapFreq < 20\n", c.Flag("MinMapFreq").String())
	}

	opt.Correct, ok = c.Flag("Correct").Get().(bool)
	if !ok {
		log.Fatalf("[checkArgsSF] argument 'Correct': %v set error\n ", c.Flag("Correct").String())
	}
	//opt.ReNameID, ok = c.Flag("ReNameID").Get().(bool)
	//if !ok {
	//	log.Fatalf("[checkArgsSF] argument 'Correct': %v set error\n ", c.Flag("ReNameID").String())
	//}
	opt.Graph, ok = c.Flag("Graph").Get().(bool)
	if !ok {
		log.Fatalf("[checkArgsSF] argument 'Graph': %v set error\n ", c.Flag("Graph").String())
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
	opt, suc := checkArgsSF(c)
	if suc == false {
		log.Fatalf("[Smfy] check Arguments error, opt: %v\n", opt)
	}
	//opt.MaxMapEdgeLen = tmp.MaxMapEdgeLen
	fmt.Printf("Arguments: %v\n", opt)

	// set package-level variable
	//Kmerlen = opt.Kmer

	DBGStatfn := opt.Prefix + ".DBG.stat"
	nodesSize, edgesSize := DBGStatReader(DBGStatfn)
	//nodesSize := len(nodeMap)
	fmt.Printf("[Smfy] len(nodesArr): %v, length of edge array: %v\n", nodesSize, edgesSize)
	nodesfn := opt.Prefix + ".nodes.Arr"
	nodesArr := make([]DBGNode, nodesSize)
	fc := make(chan int, 1)
	go NodesArrReader(nodesfn, nodesArr, opt.Kmer, fc)
	// read edges file
	edgesfn := opt.Prefix + ".edges.fa.zst"
	edgesArr := ReadEdgesFromFile(edgesfn, edgesSize, opt.Kmer)
	<-fc

	//PrintTmpNodesArr(nodesArr, opt.Prefix+".beforeSmfy")
	//gfn1 := opt.Prefix + ".beforeSmfyDBG.dot"
	//GraphvizDBGArr(nodesArr, edgesArr, gfn1)

	//gfn := opt.Prefix + ".smfyDBG.dot"
	//GraphvizDBG(nodeMap, edgesArr, gfn)
	// reconstruct consistence De Bruijn Graph
	// ReconstructConsistenceDBG(nodeMap, edgesArr)
	for i := 0; i < 3; i++ {
		fmt.Printf("[Smfy] %v cycle smfy DBG....\n", i)
		SmfyDBG(nodesArr, edgesArr, opt)
	}
	/*if opt.ReNameID {
	edgesArr = SetEdgeID(nodesArr, edgesArr)
	nodesArr = SetNodeID(nodesArr, edgesArr)
	}*/
	CheckDBGSelfCycle(nodesArr, edgesArr, opt.Kmer)
	CheckInterConnectivity(edgesArr, nodesArr)
	AddNodeInfo2DBGEdgeArr(edgesArr, nodesArr)
	//MakeSelfCycleEdgeOutcomingToIncoming(nodesArr, edgesArr, opt)
	// set the unique edge of edgesArr
	//PrintTmpDBG(nodesArr, edgesArr, opt.Prefix+".afterSmfy")
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
	// 	if len(edgesArr[i].Ks) == len(edgesArr[i-1].Ks) {
	// 		fmt.Printf("[Smfy] outcoming: %v, %v, len: %d\n", edgesArr[i-1].Ks[Kmerlen-1], edgesArr[i].Ks[Kmerlen-1], len(edgesArr[i].Ks))
	// 	}
	// }
	// Debug code
	smfyNodesfn := opt.Prefix + ".nodes.smfy.Arr"
	fc1 := make(chan int, 1)
	defer close(fc1)
	go NodesArrWriter(nodesArr, smfyNodesfn, fc1)

	smfyEdgesfn := opt.Prefix + ".edges.smfy.fa.zst"
	StoreEdgesToFn(smfyEdgesfn, edgesArr, opt.TipMaxLen)
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

func Transform2Path(p32 []uint32) (path []uint32) {
	path = make([]uint32, len(p32))
	for i, id := range p32 {
		path[i] = uint32(id)
	}
	return
}
func Transform2Path2(p32 []uint32, path []uint32) []uint32 {
	la := len(p32)
	if len(path) < la {
		path = make([]uint32, la)
	} else {
		path = path[:la]
	}
	for i, id := range p32 {
		path[i] = uint32(id)
	}
	return path
}

func Transform2EdgeFreqArr(p32 []uint32, path []EdgeFreq) []EdgeFreq {
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

func Equaluint32Arr(arr1, arr2 []uint32) bool {
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

/*func GetUniquePathIdxArr(pa []PathFreq, pf PathFreq, eID uint32, edgesArr []DBGEdge, EqualFlag bool) (idxArr, ix2Arr []int) {
	l1 := len(pf.IDArr)
	idxArr = make([]int, 0, l1/3+1)
	ix2Arr = make([]int, 0, l1/3+1)
	idx := IndexUint32(pf.IDArr, eID)
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

				idx2 := IndexUint32(pf2.IDArr, eID)
				y := IndexUint32(pf2.IDArr, id)
				if x < idx {
					if y > idx2 {
						pf2.IDArr = GetReverseUint32Arr(pf2.IDArr)
						y = l2 - 1 - y
					}
					if y < x || l2-y > l1-x {
						continue
					}
				} else {
					if y < idx2 {
						pf2.IDArr = GetReverseUint32Arr(pf2.IDArr)
						y = l2 - 1 - y
					}
					if l1-x > l2-y || y > x {
						continue
					}
				}
				lmin := MinInt(x, y)
				rmin := MinInt(l1-x, l2-y)
				if Equaluint32Arr(pf.IDArr[x-lmin:x+rmin], pf2.IDArr[y-lmin:y+rmin]) {
					idxArr = append(idxArr, x)
					ix2Arr = append(ix2Arr, ix2)
				}
			}
		}
	}
	return
}*/

func GetShareNodeID(e, e1 *DBGEdge) uint32 {
	if e.StartNID == e1.StartNID || e.StartNID == e1.EndNID {
		return e.StartNID
	} else if e.EndNID == e1.StartNID || e.EndNID == e1.EndNID {
		return e.EndNID
	}
	var nID uint32
	return nID
}

func GetShareNodeIDArr(e, e1 *DBGEdge) (arr []uint32) {
	if e.StartNID == e1.StartNID || e.StartNID == e1.EndNID {
		arr = append(arr, e.StartNID)
	}
	if e.EndNID == e1.StartNID || e.EndNID == e1.EndNID {
		arr = append(arr, e.EndNID)
	}
	return arr
}

func GetNumInuint32Arr(arr []uint32, eID uint32) (count int) {
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
		if GetNumInuint32Arr(pf.IDArr, eID) != 1 {
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
		if GetNumInuint32Arr(pf.IDArr, eID) != 1 {
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

func GetPathSeqLen(path []uint32, edgesArr []DBGEdge, kmerlen int) (sl int) {
	for _, eID := range path {
		sl += edgesArr[eID].GetSeqLen() - (kmerlen - 1)
	}
	if sl > 0 {
		sl += kmerlen - 1
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

func PrintEdgeInfo(e *DBGEdge) string {
	s := fmt.Sprintf("eID:%d StartNID:%d EndNID:%d Incoming:%v Outcoming:%v el:%d", e.ID, e.StartNID, e.EndNID, e.EdgeIDIncoming, e.EdgeIDOutcoming, e.GetSeqLen())
	return s
}

func PrintNodeInfo(nd *DBGNode) string {
	s := fmt.Sprintf("nID:%d Incoming:%v Outcoming:%v", nd.ID, nd.EdgeIDIncoming, nd.EdgeIDOutcoming)
	return s
}

func IsConnectNode(v DBGNode, eID1, eID2 uint32) bool {
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

func IsLinearNode(nd *DBGNode, eID1, eID2 uint32) bool {
	inNum, inID := GetEdgeIDComing(nd.EdgeIDIncoming)
	outNum, outID := GetEdgeIDComing(nd.EdgeIDOutcoming)
	if inNum == 1 && outNum == 1 {
		if inID == eID1 && outID == eID2 {
			return true
		} else if inID == eID2 && outID == eID1 {
			return true
		}
	}
	return false
}

// this version just used for smfy, incoming == outcoming == 1
func ConcatEdgePathUtg(path []EdgeFreq, edgesArr []DBGEdge, nodesArr []DBGNode, kmerlen int) (e DBGEdge) {
	sl := GetEdgeFreqPathSeqLen(path, edgesArr, kmerlen)
	e.Ks = make([]byte, 0, sl)
	e.ID = path[0].ID
	e0, e1 := &edgesArr[e.ID], &edgesArr[path[1].ID]
	narr := GetShareNodeIDArr(e0, e1)
	if len(narr) > 1 {
		if IsLinearNode(&nodesArr[narr[1]], e0.ID, e1.ID) {
			narr[0] = narr[1]
		}
	}
	var strand bool
	if edgesArr[e0.ID].EndNID == narr[0] {
		strand = PLUS
		e.StartNID = e0.StartNID
		e.EdgeIDIncoming = e0.EdgeIDIncoming
	} else {
		strand = MINUS
		e.StartNID = e0.EndNID
		e.EdgeIDIncoming = e0.EdgeIDOutcoming
	}
	e.Ks = append(e.Ks, e0.Ks...)
	if !strand {
		ReverseCompByteArr(e.Ks)
	}
	//fmt.Printf("[ConcatEdgePathUtg]eID:%d strand:%v\n", e0.ID, strand)
	var e2 *DBGEdge
	for i := 1; i < len(path); i++ {
		e2 = &edgesArr[path[i].ID]
		strand = GetNextEdgeStrand2(e0, e2, strand, FORWARD)
		//fmt.Printf("[ConcatEdgePathUtg]eID:%d strand:%v\n", e2.ID, strand)
		if strand {
			e.Ks = append(e.Ks, e2.Ks[kmerlen-1:]...)
			e.EndNID = e2.EndNID
			e.EdgeIDOutcoming = e2.EdgeIDOutcoming
		} else {
			e.Ks = append(e.Ks, e2.Ks[:len(e2.Ks)-(kmerlen-1)]...)
			ReverseCompByteArr(e.Ks[len(e.Ks)-(len(e2.Ks)-(kmerlen-1)):])
			e.EndNID = e2.StartNID
			e.EdgeIDOutcoming = e2.EdgeIDIncoming
		}
		e0 = e2
	}
	//fmt.Printf("[ConcatEdgePathUtg]seqLen:%d len(e.Ks):%d e:%v\n", sl, len(e.Ks), e)
	return
}

func ConcatEdgePathUtg2(path []EdgeFreq, edgesArr []DBGEdge, nodesArr []DBGNode, kmerlen int) (e DBGEdge) {
	/*e1, e2 := &edgesArr[path[0].ID], &edgesArr[path[1].ID]
	if e1.StartNID == e1.EndNID {
		if IsInComing(e1.EdgeIDIncoming, e2.ID) {
			if e2.EndNID == e1.EndNID {
				e.= GetRCUnitig(e1.

			}
		} else {
			if e2.EndNID == e1.EndNID {
				e.= GetRCUnitig(e1.
				e.StartNID, e.EndNID = e.EndNID, e.StartNID
				e.EdgeIDIncoming, e.EdgeIDOutcoming = e.EdgeIDOutcoming, e.EdgeIDIncoming
			}
		}
	} else if e1.EndNID != nID {
		e.= GetRCUnitig(e1.
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
			ne.= GetRCUnitig(ne.
			ne.StartNID, ne.EndNID = ne.EndNID, ne.StartNID
			ne.EdgeIDIncoming, ne.EdgeIDOutcoming = ne.EdgeIDOutcoming, ne.EdgeIDIncoming
		}
		//fmt.Printf("[ConcatEdgePathUtg]strand:%t %s\n", strand, PrintEdgeInfo(&ne))
		e.= ConcatEdges(e. ne. kmerlen)
		if len(e.Ks) == 0 {
			return
		}
		e.EndNID = ne.EndNID
		lastE = &edgesArr[path[i].ID]
	}*/
	return
}

/*
func RePlaceNGSPath(pf PathFreq, path []uint32, eID uint32) PathFreq {
	idx := IndexUint32(pf.IDArr, path[0])
	idx2 := IndexUint32(pf.IDArr, path[len(path)-1])
	ok := false
	if idx >= 0 && idx2 >= 0 && idx+idx2 > 0 {
		if idx > idx2 {
			path = GetReverseUint32Arr(path)
		}
		idx, idx2 = idx2, idx
		min := MinInt(len(pf.IDArr)-idx, len(path))
		if Equaluint32Arr(pf.IDArr[idx:idx+min], path[:min]) {
			ok = true
			var npf PathFreq
			pl := idx + 1
			if len(pf.IDArr)-idx > len(path) {
				pl += len(pf.IDArr) - idx - len(path)
			}
			npf.IDArr = make([]uint32, pl)
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
		if Equaluint32Arr(pf.IDArr[idx:idx+pl], path[:pl]) {
			ok = true
			var npf PathFreq
			npf.IDArr = make([]uint32, idx+1)
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
		if Equaluint32Arr(pf.IDArr[idx2-pl:idx2+1], path[len(path)-pl:]) {
			ok = true
			var npf PathFreq
			pl := 1
			if idx2+1 > len(pf.IDArr) {
				pl += len(pf.IDArr) - (idx2 + 1)
			}
			npf.IDArr = make([]uint32, pl)
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

func GetNumInPath0(pa []PathFreq, eID uint32) (count int) {
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

func GetAvgEdgeFreq(ep []EdgeFreq) (avgFreq int) {
	avgFreq = 0
	for _, ef := range ep {
		avgFreq += int(ef.Freq)
	}
	avgFreq /= len(ep)

	return avgFreq
}
