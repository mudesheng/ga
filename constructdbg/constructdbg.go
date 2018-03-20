package constructdbg

import (
	"bufio"
	"compress/gzip"
	"encoding/binary"
	"encoding/gob"
	"fmt"
	"io"
	"log"
	"math"
	"os"
	"reflect"
	"runtime"
	"strconv"
	"strings"
	"sync"
	"time"

	"github.com/mudesheng/ga/bnt"
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

type DBG_MAX_INT uint32

type DBGNode struct {
	ID              DBG_MAX_INT
	EdgeIDIncoming  [bnt.BaseTypeNum]DBG_MAX_INT // the ID of EdgeID inclming
	EdgeIDOutcoming [bnt.BaseTypeNum]DBG_MAX_INT
	Seq             []byte // kmer node seq, use 2bit schema
	SubGID          uint8  // Sub graph ID, most allow 2**8 -1 subgraphs
	Flag            uint8  // from low~high, 1:Process, 2:
}

const (
	DIN         uint8 = 1 // note DBGNode and DBGEdge relationship
	DOUT        uint8 = 2 // note DBGNode and DBGEdge relationship
	PLUS, MINUS       = true, false
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
	IDArr []DBG_MAX_INT
	Freq  int
}

type DBGEdge struct {
	ID       DBG_MAX_INT
	StartNID DBG_MAX_INT // start node ID
	EndNID   DBG_MAX_INT // end node ID
	CovD     uint16      // Coverage Depth
	Flag     uint8       //
	Utg      Unitig
	PathMat  []Path // read Path matrix
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

func (e *DBGEdge) SetUniqueFlag() {
	e.Flag = e.Flag | 0x4
}

func (e *DBGEdge) ResetUniqueFlag() {
	e.Flag = e.Flag & 0xFB
}

func (e *DBGEdge) SetProcessFlag() {
	e.Flag = e.Flag | 0x1
}

func (e *DBGEdge) ResetProcessFlag() {
	e.Flag = e.Flag & 0xFE
}

func (e *DBGEdge) SetDeleteFlag() {
	e.Flag = e.Flag | 0x2
}

func (e *DBGEdge) GetSeqLen() int {
	return len(e.Utg.Ks)
}

var MIN_KMER_COUNT uint16 = 2
var BACKWARD uint8 = 1
var FORWARD uint8 = 2
var Kmerlen int = -1

func ReverseDBG_MAX_INTArr(path []DBG_MAX_INT) []DBG_MAX_INT {
	for i := 0; i < len(path)/2; i++ {
		path[i], path[len(path)-1-i] = path[len(path)-1-i], path[i]
	}

	return path
}

func readUniqKmer(uniqkmergzfn string, cs chan constructcf.ReadBnt, kmerlen, numCPU int) {
	uniqkmergzfp, err := os.Open(uniqkmergzfn)
	if err != nil {
		log.Fatal(err)
	}
	ukgzfp, err := gzip.NewReader(uniqkmergzfp)
	if err != nil {
		log.Fatal(err)
	}
	ukbuffp := bufio.NewReader(ukgzfp)
	defer ukgzfp.Close()
	defer uniqkmergzfp.Close()
	var processNumKmer int
	// var rsb constructcf.ReadSeqBucket
	eof := false
	KBntByteNum := (kmerlen + bnt.NumBaseInByte - 1) / bnt.NumBaseInByte
	for !eof {
		// var kmer []byte
		b := make([]byte, KBntByteNum)
		// n, err := ukbuffp.Read(b)
		err := binary.Read(ukbuffp, binary.LittleEndian, b)
		/*if n < KBntByteNum && err == nil {
			n1, err1 := ukbuffp.Read(b[n:])
			if err1 != nil {
				log.Fatalf("[readUniqKmer] read kmer seq err1: %v\n", err1)
			}
			n += n1
		}*/
		// fmt.Printf("[readUniqKmer]len(b): %v,  b: %v\n", len(b), b)
		if len(b) != KBntByteNum || err != nil {
			if err == io.EOF {
				eof = true
				continue
			} else if len(b) != KBntByteNum {
				log.Fatalf("[readUniqKmer] read kmer seq length: %v != KBntByteNum[%v]\n", len(b), KBntByteNum)
			} else {
				log.Fatalf("[readUniqKmer] read kmer seq err: %v\n", err)
			}
			// fmt.Printf("[readUniqKmer]: read %d bytes\n", n)
		}
		// if n != KBntByteNum {
		// 	// log.Fatalf("[readUniqKmer] read kmer seq err: n(%d) != KBntByteNum(%d)\n", n, KBntByteNum)
		// }
		var rb constructcf.ReadBnt
		rb.Seq = b
		// if len(rb.Seq) != kmerlen {
		// 	log.Fatalf("[readUniqKmer] len(rb.Seq) != kmerlen\n")
		// } else {
		rb.Length = kmerlen
		// }
		cs <- rb
		processNumKmer++
	}

	fmt.Printf("[readUniqKmer] total read kmer number is : %d\n", processNumKmer)
	// send read finish signal
	close(cs)
}

func UniDirectExtend(nb constructcf.ReadBnt, cf cuckoofilter.CuckooFilter, min_kmer_count uint16, direction uint8) (ec int) {
	var nBnt constructcf.ReadBnt
	nBnt.Seq = make([]byte, Kmerlen)
	if direction == FORWARD {
		copy(nBnt.Seq[:Kmerlen-1], nb.Seq)
		for i := 0; i < bnt.BaseTypeNum; i++ {
			b := byte(i)
			nBnt.Seq[Kmerlen-1] = b
			ks := constructcf.GetReadBntKmer(nBnt, 0, Kmerlen)
			rs := constructcf.ReverseComplet(ks)
			if ks.BiggerThan(rs) {
				ks, rs = rs, ks
			}
			if count := cf.GetCountAllowZero(ks.Seq); count >= min_kmer_count {
				ec++
			}
		}
	} else { // direction == BACKWARD
		copy(nBnt.Seq[1:], nb.Seq)
		for i := 0; i < bnt.BaseTypeNum; i++ {
			b := byte(i)
			nBnt.Seq[0] = b
			ks := constructcf.GetReadBntKmer(nBnt, 0, Kmerlen)
			rs := constructcf.ReverseComplet(ks)
			if ks.BiggerThan(rs) {
				ks, rs = rs, ks
			}
			if count := cf.GetCountAllowZero(ks.Seq); count >= min_kmer_count {
				ec++
			}
		}
	}

	return
}

func ExtendNodeKmer(nodeBnt constructcf.ReadBnt, cf cuckoofilter.CuckooFilter, min_kmer_count uint16, direction uint8, base byte) (nd DBGNode, baseBnt byte, baseCount uint16) {
	var nBnt constructcf.ReadBnt
	nBnt.Seq = make([]byte, Kmerlen+1)
	nd.Seq = constructcf.GetReadBntKmer(nodeBnt, 0, Kmerlen-1).Seq
	copy(nBnt.Seq[1:], nodeBnt.Seq)
	for i := 0; i < bnt.BaseTypeNum; i++ {
		bi := byte(i)
		if direction == BACKWARD && bi == base {
			nd.EdgeIDIncoming[bi] = 1
		} else {
			nBnt.Seq[0] = bi
			ks := constructcf.GetReadBntKmer(nBnt, 0, Kmerlen)
			rs := constructcf.ReverseComplet(ks)
			if ks.BiggerThan(rs) {
				ks, rs = rs, ks
			}
			// fmt.Printf("[ExtendNodeKmer]ks: %v, rs: %v\n", ks, rs)
			count := cf.GetCountAllowZero(ks.Seq)
			if count >= min_kmer_count {
				var nb constructcf.ReadBnt
				nb.Seq = append(nb.Seq, nBnt.Seq[:Kmerlen-1]...)
				nb.Length = len(nb.Seq)
				if ec := UniDirectExtend(nb, cf, min_kmer_count, BACKWARD); ec > 0 {
					if direction == FORWARD {
						baseBnt = uint8(i)
						baseCount = count
					}
					//fmt.Printf("[ExtendNodeKmer] BACKWARD bi: %v, ks.Seq: %v\n", bi, ks.Seq)
					nd.EdgeIDIncoming[bi] = 1
				}
			}
		}

		if direction == FORWARD && bi == base {
			nd.EdgeIDOutcoming[bi] = 1
		} else {
			nBnt.Seq[cf.Kmerlen] = bi
			ks := constructcf.GetReadBntKmer(nBnt, 1, cf.Kmerlen)
			rs := constructcf.ReverseComplet(ks)
			if ks.BiggerThan(rs) {
				ks, rs = rs, ks
			}
			count := cf.GetCountAllowZero(ks.Seq)
			if count >= min_kmer_count {
				var nb constructcf.ReadBnt
				nb.Seq = append(nb.Seq, nBnt.Seq[2:]...)
				nb.Length = len(nb.Seq)
				if ec := UniDirectExtend(nb, cf, min_kmer_count, FORWARD); ec > 0 {
					if direction == BACKWARD {
						baseBnt = uint8(i)
						baseCount = count
					}
					//fmt.Printf("[ExtendNodeKmer] FORWARD bi: %v, ks.Seq: %v\n", bi, ks.Seq)
					nd.EdgeIDOutcoming[bi] = 1
				}
			}
		}
	}

	return
}

func GetMinDBGNode(nd DBGNode, kmerlen int) (minN DBGNode) {
	var nodeBnt constructcf.ReadBnt
	nodeBnt.Seq = nd.Seq
	nodeBnt.Length = kmerlen - 1
	rnode := constructcf.ReverseComplet(nodeBnt)
	//fmt.Printf("[GetMinDBGNode] node: %v\nRC node: %v\n", nodeBnt, rnode)
	if nodeBnt.BiggerThan(rnode) {
		minN.Seq = rnode.Seq
		for i := 0; i < bnt.BaseTypeNum; i++ {
			minN.EdgeIDIncoming[i] = nd.EdgeIDOutcoming[bnt.BaseTypeNum-1-i]
			minN.EdgeIDOutcoming[i] = nd.EdgeIDIncoming[bnt.BaseTypeNum-1-i]
		}
	} else {
		minN = nd
	}
	return
}

func paraLookupComplexNode(cs chan constructcf.ReadBnt, wc chan DBGNode, cf cuckoofilter.CuckooFilter) {
	// var wrsb constructcf.ReadSeqBucket
	for {
		rb, ok := <-cs
		if !ok {
			var nd DBGNode
			wc <- nd
			break
		}
		// if rsb.Count < constructcf.ReadSeqSize {
		// 	fmt.Printf("rsb.ReadBuf length is : %d\n", len(rsb.ReadBuf))
		// }
		// if found kmer count is 1 , this kmer will be ignore, and skip this branch
		extRBnt := constructcf.ExtendReadBnt2Byte(rb)
		{ // check fisrt node of kmer
			var nodeBnt constructcf.ReadBnt
			nodeBnt.Seq = make([]byte, cf.Kmerlen-1)
			copy(nodeBnt.Seq, extRBnt.Seq[:cf.Kmerlen-1])
			nodeBnt.Length = len(nodeBnt.Seq)
			// fmt.Printf("[paraLookupComplexNode] nodeBnt : %v\n", nodeBnt)
			nd, _, _ := ExtendNodeKmer(nodeBnt, cf, MIN_KMER_COUNT, FORWARD, extRBnt.Seq[cf.Kmerlen-1])
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
				tn := GetMinDBGNode(nd, cf.Kmerlen)
				//fmt.Printf("[paraLookupComplexNode] node: %v\nmin of RC node: %v\n", nd, tn)
				wc <- tn
			}
		}

		{ // check second node of kmer
			var nodeBnt constructcf.ReadBnt
			nodeBnt.Seq = make([]byte, cf.Kmerlen-1)
			copy(nodeBnt.Seq, extRBnt.Seq[1:])
			nodeBnt.Length = len(nodeBnt.Seq)
			nd, _, _ := ExtendNodeKmer(nodeBnt, cf, MIN_KMER_COUNT, BACKWARD, extRBnt.Seq[0])
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
				tn := GetMinDBGNode(nd, cf.Kmerlen)
				wc <- tn
			}
		}
	}
}

func constructNodeMap(complexKmergzfn string) (map[string]DBGNode, DBG_MAX_INT) {
	nodeID := DBG_MAX_INT(1)
	nodeMap := make(map[string]DBGNode)
	ckfp, err := os.Open(complexKmergzfn)
	if err != nil {
		log.Fatalf("[constructNodeMap] open %s file failed: %v\n", complexKmergzfn, err)
	}
	defer ckfp.Close()
	ckgzfp, err := gzip.NewReader(ckfp)
	if err != nil {
		log.Fatalf("[constructNodeMap] gzip.NewReader err: %v\n", err)
	}
	defer ckgzfp.Close()
	//ckbuffp := bufio.NewReader(ckgzfp)
	NBntByteLen := (Kmerlen - 1 + bnt.BaseTypeNum - 1) / bnt.BaseTypeNum
	eof := false
	readnodeNum := 0
	for !eof {
		var node DBGNode
		node.Seq = make([]byte, NBntByteLen)
		err := binary.Read(ckgzfp, binary.LittleEndian, node.Seq)
		if err != nil {
			if err == io.EOF {
				eof = true
				continue
			} else {
				log.Fatalf("[constructNodeMap] err: %v\n", err)
			}
		}
		if err := binary.Read(ckgzfp, binary.LittleEndian, &node.EdgeIDIncoming); err != nil {
			log.Fatalf("[constructNodeMap] read file: %v err\n", complexKmergzfn)
		}
		if err := binary.Read(ckgzfp, binary.LittleEndian, &node.EdgeIDOutcoming); err != nil {
			log.Fatalf("[constructNodeMap] read file: %v err\n", complexKmergzfn)
		}
		// fmt.Fprintf(os.Stderr, "[constructNodeMap] node: %v\n", node)
		// n, err := ckfp.Read(b)
		// if n != NBntByteLen {
		// 	log.Fatalf("[constructNodeMap] read node seq err: n(%d) != NBntByteLen(%d\n)", n, NBntByteLen)
		// }
		readnodeNum++
		node.ID = nodeID
		if _, ok := nodeMap[string(node.Seq)]; ok == false {
			nodeMap[string(node.Seq)] = node
			nodeID++
		} else {
			fmt.Printf("[constructNodeMap] repeat node: %v\n", node)
		}
	}

	fmt.Printf("[constructNodeMap] read node number is : %d\n", readnodeNum)

	return nodeMap, nodeID
}

func AddNodeToNodeMap(node DBGNode, nodeMap map[string]DBGNode, nodeID DBG_MAX_INT) DBG_MAX_INT {
	if node.Flag != 1 {
		log.Fatalf("[AddNodeToNodeMap] found node.Flag: %v != 1\n", node.Flag)
	}
	if _, ok := nodeMap[string(node.Seq)]; ok == false {
		node.ID = nodeID
		nodeMap[string(node.Seq)] = node
		nodeID++
	} else {
		log.Fatalf("[AddNodeToNodeMap] node: %v has been exist in the nodeMap\n", node)
	}

	return nodeID
}

func CollectAddedDBGNode(anc <-chan DBGNode, nc chan<- DBGNode, readNodeMapFinishedC <-chan int) {
	var narr []DBGNode
loop:
	for {
		select {
		case nd := <-anc:
			narr = append(narr, nd)
		case <-readNodeMapFinishedC:
			break loop
		}
	}

	for j := 0; j < 2; j++ {
		for i := 0; i < len(narr); i++ {
			select {
			case nd := <-anc:
				narr = append(narr, nd)
			}
			nc <- narr[i]
		}
		if j == 0 {
			time.Sleep(2 * time.Second)
		}
	}
	fmt.Printf("[CollectAddedDBGNode] added node number is: %v\n", len(narr))
	close(nc)
}

// add new Node to the nodeMap and check node edge has been output
func ChangeNodeMap(nodeMap map[string]DBGNode, anc chan<- DBGNode, finishedC <-chan int, nIEC <-chan NodeInfoByEdge, flagNIEC chan<- NodeInfoByEdge, Kmerlen int, nodeID DBG_MAX_INT) (nID DBG_MAX_INT, edgeID DBG_MAX_INT) {
	oldNodeID := nodeID
	edgeID = DBG_MAX_INT(2)
loop:
	for {
		select {
		case <-finishedC:
			break loop
		case nie := <-nIEC:
			v1 := nodeMap[string(nie.n1.Seq)]
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
			nodeMap[string(v1.Seq)] = v1
			nd := nie.n2
			if len(nd.Seq) > 0 {
				tn := GetMinDBGNode(nd, Kmerlen)
				v2, ok := nodeMap[string(tn.Seq)]
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
				if ok {
					nodeMap[string(v2.Seq)] = v2
					nie.n2 = v2
				} else {
					v2.Flag = 1
					nodeID = AddNodeToNodeMap(v2, nodeMap, nodeID)
					//fmt.Printf("[ChangeNodeMap] v2: %v\nnodeID: %v\n", v2, nodeID-1)
					t := nodeMap[string(v2.Seq)]
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
}

// ReadDBGNodeToChan, read DBG nodeMap and simultaneously add new node to the nodeMap
func ReadDBGNodeToChan(nodeMap map[string]DBGNode, nc chan<- DBGNode, readNodeMapFinished chan<- int) {

	for _, value := range nodeMap {
		if value.Flag == 0 {
			nc <- value
		}
		//value.Flag = uint8(1)
		//nodeMap[string(value.Seq)] = value
		//}
	}

	// notice Function ChangeNodeMap() has finished read nodeMap
	time.Sleep(time.Second * 5)
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

func ReverseUint8Arr(ua []uint8) {
	lua := len(ua)
	for i := 0; i < lua/2; i++ {
		ua[i], ua[lua-1-i] = ua[lua-1-i], ua[i]
	}
}

func GetEdges(cf cuckoofilter.CuckooFilter, nBnt constructcf.ReadBnt, count uint8, direction uint8, MIN_KMER_COUNT uint16) (edge DBGEdge, nd DBGNode) {
	if direction == FORWARD {
		edge.Utg.Ks = append(edge.Utg.Ks, nBnt.Seq...)
		edge.Utg.Kq = make([]uint8, cf.Kmerlen-1)
		edge.Utg.Kq = append(edge.Utg.Kq, count)
		var tbnt constructcf.ReadBnt
		tbnt.Seq = nBnt.Seq[1:]
		tbnt.Length = len(tbnt.Seq)
		bi := nBnt.Seq[0]
		for {
			node, baseBnt, baseCount := ExtendNodeKmer(tbnt, cf, MIN_KMER_COUNT, BACKWARD, bi)
			var leftcount, rightcount int
			for i := 0; i < bnt.BaseTypeNum; i++ {
				if node.EdgeIDIncoming[i] == 1 {
					leftcount++
				}
				if node.EdgeIDOutcoming[i] == 1 {
					rightcount++
				}
			}
			if baseCount >= MIN_KMER_COUNT && leftcount == 1 && rightcount == 1 {
				edge.Utg.Ks = append(edge.Utg.Ks, baseBnt)
				edge.Utg.Kq = append(edge.Utg.Kq, uint8(baseCount))
				bi = tbnt.Seq[0]
				tbnt.Seq = tbnt.Seq[1:]
				tbnt.Seq = append(tbnt.Seq, baseBnt)
			} else {
				if leftcount > 1 || rightcount > 1 {
					nd = node
				}
				break
			}
		}
	} else {
		edge.Utg.Ks = append(edge.Utg.Ks, nBnt.Seq[0])
		edge.Utg.Kq = append(edge.Utg.Kq, count)
		var tbnt constructcf.ReadBnt
		bi := nBnt.Seq[cf.Kmerlen-1]
		tbnt.Seq = nBnt.Seq[:cf.Kmerlen-1]
		tbnt.Length = len(tbnt.Seq)
		for {
			node, baseBnt, baseCount := ExtendNodeKmer(tbnt, cf, MIN_KMER_COUNT, FORWARD, bi)
			var leftcount, rightcount int
			for i := 0; i < bnt.BaseTypeNum; i++ {
				if node.EdgeIDIncoming[i] == 1 {
					leftcount++
				}
				if node.EdgeIDOutcoming[i] == 1 {
					rightcount++
				}
			}
			if baseCount >= MIN_KMER_COUNT && leftcount == 1 && rightcount == 1 {
				edge.Utg.Ks = append(edge.Utg.Ks, baseBnt)
				edge.Utg.Kq = append(edge.Utg.Kq, uint8(baseCount))
				bi = tbnt.Seq[cf.Kmerlen-2]
				var seq []byte
				seq = append(seq, baseBnt)
				seq = append(seq, tbnt.Seq[:cf.Kmerlen-2]...)
				tbnt.Seq = seq
			} else {
				if leftcount > 1 || rightcount > 1 {
					nd = node
				}
				// Reverse the seq and quality count
				ReverseByteArr(edge.Utg.Ks)
				ReverseUint8Arr(edge.Utg.Kq)

				edge.Utg.Ks = append(edge.Utg.Ks, nBnt.Seq[1:]...)
				tmp := make([]uint8, cf.Kmerlen-1)
				edge.Utg.Kq = append(edge.Utg.Kq, tmp...)
				break
			}
		}
	}

	return
}

var mu sync.Mutex

func paraGenerateDBGEdges(nc <-chan DBGNode, nIEC chan<- NodeInfoByEdge, flagNIEC <-chan NodeInfoByEdge, cf cuckoofilter.CuckooFilter, wc chan DBGEdge) {
	for {
		node, ok := <-nc
		if !ok {
			var e DBGEdge
			wc <- e
			break
		}
		// read edge seq from cuckoofilter
		var rb constructcf.ReadBnt
		rb.Seq = node.Seq
		rb.Length = cf.Kmerlen - 1
		extRBnt := constructcf.ExtendReadBnt2Byte(rb)
		// leftcount, rightcount, _, _ := ExtendNodeKmer(extRBnt, cf, MIN_KMER_COUNT, FORWARD)
		// fmt.Printf("[paraGenerateDBGEdges] leftcount: %d, rightcount: %d\n", leftcount, rightcount)
		for i := uint(0); i < bnt.BaseTypeNum; i++ {
			bi := byte(i)
			if node.EdgeIDIncoming[i] == 1 {
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
				if count < MIN_KMER_COUNT {
					log.Fatalf("[paraGenerateDBGEdges] found count[%v]: < [%v], node: %v", count, MIN_KMER_COUNT, node)
				}
				// get edge sequence
				edge, nd := GetEdges(cf, nBnt, uint8(count), BACKWARD, MIN_KMER_COUNT)
				//writedEdge := false
				if len(nd.Seq) > 0 || len(edge.Utg.Ks) > 2*cf.Kmerlen {
					var nIE NodeInfoByEdge
					nIE.n1 = node
					nIE.c1 = DIN
					nIE.i1 = i
					nIE.n2 = nd
					nIE.c2 = DOUT
					c := edge.Utg.Ks[cf.Kmerlen-1]
					nIE.i2 = uint(c)
					mu.Lock()
					//fmt.Printf("[paraGenerateDBGEdges]nodeID: %v, Incoming nIE: %v\n", node.ID, nIE)
					nIEC <- nIE
					fNIE := <-flagNIEC
					//fmt.Printf("[paraGenerateDBGEdges]nodeID: %v, Incoming fNIE: %v\n", node.ID, fNIE)
					mu.Unlock()

					if fNIE.Flag == true {
						edge.StartNID = fNIE.n2.ID
						edge.EndNID = fNIE.n1.ID
						edge.ID = fNIE.edgeID
						wc <- edge
					}
				}
			}

			if node.EdgeIDOutcoming[i] == 1 {
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
				if count < MIN_KMER_COUNT {
					log.Fatalf("[paraGenerateDBGEdges] found count[%v]: < [%v], node: %v", count, MIN_KMER_COUNT, node)
				}
				edge, nd := GetEdges(cf, nBnt, uint8(count), FORWARD, MIN_KMER_COUNT)
				//fmt.Printf("[paraGenerateDBGEdges]nodeID: %v, edge: %v\nnode: %v\n", node.ID, edge, nd)
				// writedEdge := false
				if len(nd.Seq) == 0 && len(edge.Utg.Ks) <= 2*cf.Kmerlen {
					continue
				}
				var nIE NodeInfoByEdge
				nIE.n1 = node
				nIE.c1 = DOUT
				nIE.i1 = i
				nIE.n2 = nd
				nIE.c2 = DIN
				c := edge.Utg.Ks[len(edge.Utg.Ks)-cf.Kmerlen]
				nIE.i2 = uint(c)
				mu.Lock()
				//fmt.Printf("[paraGenerateDBGEdges]nodeID: %v, Outcoming nIE: %v\n", node.ID, nIE)
				nIEC <- nIE
				fNIE := <-flagNIEC
				//fmt.Printf("[paraGenerateDBGEdges]nodeID: %v, Outcoming fNIE: %v\n", node.ID, fNIE)
				mu.Unlock()

				if fNIE.Flag == true {
					edge.StartNID = fNIE.n1.ID
					edge.EndNID = fNIE.n2.ID
					edge.ID = fNIE.edgeID
					wc <- edge
				}
			}
		}
	}
}

// write edges seq to the file
func WriteEdgesToFn(edgesfn string, wc <-chan DBGEdge, numCPU int, finishC chan<- int) {
	edgesNum := 0
	edgesfp, err := os.Create(edgesfn)
	if err != nil {
		log.Fatalf("[WriteEdgesToFn] Open file %s failed: %v\n", edgesfn, err)
	}
	edgesbuffp := bufio.NewWriter(edgesfp)
	defer edgesfp.Close()

	// edgesgzfp := gzip.NewWriter(edgesbuffp)
	// defer edgesgzfp.Close()
	finishNum := 0
	for {
		ei := <-wc
		if len(ei.Utg.Ks) == 0 {
			if ei.StartNID != 0 || ei.EndNID != 0 {
				log.Fatalf("[WriteEdgesToFn] err edge: %v\n", ei)
			}
			fmt.Printf("[WriteEdgesToFn] end edge: %v\n", ei)
			finishNum++
			if finishNum == numCPU {
				finishC <- 1
				break
			}
			continue
		}
		edgesNum++
		fmt.Fprintf(edgesbuffp, ">%d\t%d\t%d\n", ei.ID, ei.StartNID, ei.EndNID)
		// fmt.Fprintf(edgesgzfp, "%s\n", ei.Utg.Ks)
		for i := 0; i < len(ei.Utg.Ks); i++ {
			if ei.Utg.Ks[i] > 3 || ei.Utg.Ks[i] < 0 {
				log.Fatalf("[WriteEdgesToFn] not correct base: %d:%d\n", i, ei.Utg.Ks[i])
			}
			fmt.Fprintf(edgesbuffp, "%c", bnt.BitNtCharUp[ei.Utg.Ks[i]])
		}
		fmt.Fprintf(edgesbuffp, "\n")
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
	edgesbuffp.Flush()

	fmt.Printf("[WriteEdgesToFn] the writed file edges number is %d\n", edgesNum)
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

func GenerateDBGEdges(nodeMap map[string]DBGNode, cf cuckoofilter.CuckooFilter, edgesfn string, numCPU int, nodeID DBG_MAX_INT) (newNodeID DBG_MAX_INT, edgeID DBG_MAX_INT) {
	bufsize := 50
	nc := make(chan DBGNode)
	wc := make(chan DBGEdge, bufsize)
	readNodeMapFinishedC := make(chan int)
	finishedC := make(chan int)
	nIEC := make(chan NodeInfoByEdge)
	flagNIEC := make(chan NodeInfoByEdge)
	// Read DBGNode to the nc
	go ReadDBGNodeToChan(nodeMap, nc, readNodeMapFinishedC)
	// parallel construct edges from cuckoofilter
	for i := 0; i < numCPU; i++ {
		go paraGenerateDBGEdges(nc, nIEC, flagNIEC, cf, wc)
	}
	// write edges Seq to the file
	go WriteEdgesToFn(edgesfn, wc, numCPU, finishedC)

	// collect added node and pass DBGNode to the nc
	anc := make(chan DBGNode)
	go CollectAddedDBGNode(anc, nc, readNodeMapFinishedC)

	// Change nodeMap monitor function
	newNodeID, edgeID = ChangeNodeMap(nodeMap, anc, finishedC, nIEC, flagNIEC, cf.Kmerlen, nodeID)

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

func NodeMapMmapWriter(nodeMap map[string]DBGNode, nodesfn string) {
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

func NodeMapMmapReader(nodesfn string) (nodeMap map[string]DBGNode) {
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

func NodesArrWriter(nodesArr []DBGNode, nodesfn string) {
	nodesfp, err := os.Create(nodesfn)
	if err != nil {
		log.Fatalf("[NodesArrWriter] file %s create error, err: %v\n", nodesfn, err)
	}
	defer nodesfp.Close()
	enc := gob.NewEncoder(nodesfp)
	err = enc.Encode(nodesArr)
	if err != nil {
		log.Fatalf("[NodeMapMmapWriter] encode err: %v\n", err)
	}
}

func NodesArrReader(nodesfn string) (nodesArr []DBGNode) {
	nodesfp, err := os.Open(nodesfn)
	if err != nil {
		log.Fatalf("[NodeMapMmapReader] open file %s failed, err:%v\n", nodesfn, err)
	}
	defer nodesfp.Close()
	dec := gob.NewDecoder(nodesfp)
	err = dec.Decode(&nodesArr)
	if err != nil {
		log.Fatalf("[NodeMapMmapReader] decode failed, err: %v\n", err)
	}

	return
}

func NodeMap2NodeArr(nodeMap map[string]DBGNode, nodesArr []DBGNode) {
	naLen := DBG_MAX_INT(len(nodesArr))
	for _, v := range nodeMap {
		if v.ID >= naLen {
			log.Fatalf("[NodeMap2NodeArr] v.ID: %v >= nodesArr len: %v\n", v.ID, naLen)
		}
		nodesArr[v.ID] = v
	}
}
func writeComplexNodesToFile(complexNodesgzFn string, wc chan DBGNode, numCPU int) (complexNodeNum int) {
	ckfp, err := os.Create(complexNodesgzFn)
	if err != nil {
		log.Fatal(err)
	}
	defer ckfp.Close()
	ckgzfp := gzip.NewWriter(ckfp)
	// if err != nil {
	// 	log.Fatal(err)
	// }
	endFlagCount := 0

	// write complex node to the file
	// complexNodeNum := 0
	for {
		nd := <-wc
		if len(nd.Seq) == 0 {
			endFlagCount++
			if endFlagCount == numCPU {
				break
			} else {
				continue
			}
		}

		fmt.Printf("node: %v\n", nd)
		if err := binary.Write(ckgzfp, binary.LittleEndian, nd.Seq); err != nil {
			log.Fatalf("[CDBG] write node seq to file err: %v\n", err)
		}
		if err := binary.Write(ckgzfp, binary.LittleEndian, nd.EdgeIDIncoming); err != nil {
			log.Fatalf("[CDBG] write node seq to file err: %v\n", err)
		}
		if err := binary.Write(ckgzfp, binary.LittleEndian, nd.EdgeIDOutcoming); err != nil {
			log.Fatalf("[CDBG] write node seq to file err: %v\n", err)
		}
		// *** test code ***
		/* extRB := constructcf.ExtendReadBnt2Byte(rb)
		fmt.Fprintf(os.Stderr, ">%v\n%v\n", complexNodeNum+1, transform2Letters(extRB.Seq))
		*/
		// *** test code ***
		// ckgzfp.Write(rb.Seq)
		// ckgzfp.Write([]byte("\n"))
		complexNodeNum += 1
	}
	ckgzfp.Close()

	return
}

func CDBG(c cli.Command) {
	fmt.Println(c.Flags(), c.Parent().Flags())

	// get set arguments
	// t0 := time.Now()
	numCPU, err := strconv.Atoi(c.Parent().Flag("t").String())
	if err != nil {
		log.Fatalf("[CDBG] the argument 't' set error, err: %v\n", err)
	}
	runtime.GOMAXPROCS(numCPU)
	prefix := c.Parent().Flag("p").String()
	// cfinfofn := prefix + ".cfInfo"
	// cf, err := cuckoofilter.RecoverCuckooFilterInfo(cfinfofn)
	// if err != nil {
	// 	log.Fatalf("[CDGB] cuckoofilter recover err: %v\n", err)
	// }

	// find complex Nodes
	cfmmapfn := prefix + ".cfmmap"
	cf, err := cuckoofilter.MmapReader(cfmmapfn)
	if err != nil {
		log.Fatalf("[CDBG] Read CuckooFilter mmap file: %v err: %v\n", cfmmapfn, err)
	}
	// fmt.Printf("[CDBG] cf: %v\n", cf)
	// fmt.Printf("[CDBG] cf.Hash[0]: %v\n", cf.Hash[0])
	Kmerlen = cf.Kmerlen
	bufsize := 100
	cs := make(chan constructcf.ReadBnt, bufsize)
	wc := make(chan DBGNode, bufsize)
	uniqkmergzfn := prefix + ".uniqkmerseq.gz"
	// read uniq kmers form file
	go readUniqKmer(uniqkmergzfn, cs, cf.Kmerlen, numCPU)

	// identify complex Nodes
	for i := 0; i < numCPU; i++ {
		go paraLookupComplexNode(cs, wc, cf)
	}

	// write complex Nodes to the file
	complexKmergzfn := prefix + ".complexNode.gz"
	complexNodeNum := writeComplexNodesToFile(complexKmergzfn, wc, numCPU)
	fmt.Printf("[CDBG] found complex Node num is : %d\n", complexNodeNum)

	// construct Node map
	nodeMap, nodeID := constructNodeMap(complexKmergzfn)
	fmt.Printf("[CDBG] assgin nodeID to : %d\n", nodeID)
	// parallel generate edges and write to file
	edgefn := prefix + ".edges"
	newNodeID, edgeID := GenerateDBGEdges(nodeMap, cf, edgefn, numCPU, nodeID)
	DBGStatfn := prefix + ".DBG.stat"
	DBGStatWriter(DBGStatfn, newNodeID, edgeID)
	// write DBG node map to the file
	nodesfn := prefix + ".nodes.mmap"
	NodeMapMmapWriter(nodeMap, nodesfn)
}

func ParseEdge(edgesbuffp *bufio.Reader) (edge DBGEdge, err error) {
	line1, err1 := edgesbuffp.ReadString('\n')
	line2, err2 := edgesbuffp.ReadString('\n')
	line3, err3 := edgesbuffp.ReadString('\n')
	if err1 != nil || err2 != nil || err3 != nil {
		if err1 == io.EOF {
			//err1 = nil
			err = io.EOF
			return
		} else {
			log.Fatalf("[ParseEdge] Read edge found err1,err2,err3: %v,%v,%v\n\tline1: %v\n\tline2: %v\n\tline3: %v\n", err1, err2, err3, line1, line2, line3)
		}
	}
	if _, err4 := fmt.Sscanf(string(line1), ">%d\t%d\t%d\n", &edge.ID, &edge.StartNID, &edge.EndNID); err4 != nil {
		log.Fatalf("[ParseEdge] fmt.Sscaf line1:%s err:%v\n", line1, err4)
	}
	// _, err4 = fmt.Sscanf(string(line2), "%s\n", &edge.Utg.Ks)
	// if err4 != nil {
	// 	log.Fatalf("[ParseEdge] Sscaf line2 err:%v\n", err4)
	// }
	// change char base to Bnt
	sline2 := string(line2[:len(line2)-1])
	edge.Utg.Ks = make([]byte, len(sline2))
	for i, v := range sline2 {
		edge.Utg.Ks[i] = bnt.Base2Bnt[v]
	}
	sline3 := string(line3[:len(line3)-1])
	edge.Utg.Kq = make([]uint8, len(sline3))
	for i, v := range sline3 {
		q := v - 33
		edge.Utg.Kq[i] = uint8(q)
	}
	if len(edge.Utg.Ks) != len(edge.Utg.Kq) {
		log.Fatalf("[ParseEdge] len(edge.Utg.Ks): %v != len(edge.Utg.Kq): %v\n", len(edge.Utg.Ks), len(edge.Utg.Ks))
	}

	return
}

func ReadEdgesFromFile(edgesfn string, edgesSize DBG_MAX_INT) (edgesArr []DBGEdge) {
	edgesArr = make([]DBGEdge, edgesSize)
	edgesfp, err := os.Open(edgesfn)
	if err != nil {
		log.Fatalf("[ReadEdgesFromFile] open file %s failed, err: %v\n", edgesfn, err)
	}
	defer edgesfp.Close()
	// edgesgzfp, err := gzip.NewReader(edgesfp)
	// if err != nil {
	// 	log.Fatalf("[ReadEdgesFromFile] read gz file: %s failed, err: %v\n", edgesfn, err)
	// }
	// defer edgesgzfp.Close()
	edgesbuffp := bufio.NewReader(edgesfp)

	var edgesNum int
	for {
		edge, err := ParseEdge(edgesbuffp)
		// num++
		//fmt.Printf("[ParseEdge] edge.ID: %v\n", edge.ID)
		if err == io.EOF {
			break
		} else if err != nil {
			log.Fatalf("[ParseEdge] file: %v, parse edge err: %v\n", edgesfn, err)
		}
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

func RevNode(node DBGNode) DBGNode {
	rnode := node
	var rBnt constructcf.ReadBnt
	rBnt.Seq = rnode.Seq
	rBnt.Length = Kmerlen - 1
	rs := constructcf.ReverseComplet(rBnt)
	rnode.Seq = rs.Seq
	for i := 0; i < bnt.BaseTypeNum; i++ {
		rnode.EdgeIDIncoming[i] = node.EdgeIDOutcoming[bnt.BntRev[i]]
		rnode.EdgeIDOutcoming[i] = node.EdgeIDIncoming[bnt.BntRev[i]]
		// rnode.EdgeIDOutcoming[bnt.BntRev[i]] = node.EdgeIDOutcoming[bnt.BaseTypeNum-1-i], node.EdgeIDIncoming[i]
	}

	return rnode
}

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
		if v.GetDeleteFlag() > 0 {
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

	for i := 1; i < len(edgesArr); i++ {
		e := edgesArr[i]
		if e.ID == 0 || e.GetDeleteFlag() > 0 {
			continue
		}
		attr := make(map[string]string)
		attr["color"] = "Blue"
		labels := "\"ID:" + strconv.Itoa(int(e.ID)) + " len:" + strconv.Itoa(len(e.Utg.Ks)) + "\""
		//labels := strconv.Itoa(int(e.ID)) + "len" + strconv.Itoa(len(e.Utg.Ks))
		//labels := strconv.Itoa(int(e.ID))
		attr["label"] = labels
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

func GetEdgeIDComing(coming [bnt.BaseTypeNum]DBG_MAX_INT) (num int, edgeID DBG_MAX_INT) {
	for _, v := range coming {
		if v > 0 {
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
		log.Fatalf("[ConcatEdges] u1: %v, u2: %v can not concatenatable\n", u1.Ks[len(u1.Ks)-kmerlen+1:], u2.Ks[:kmerlen-1])
	}
	if len(u1.Ks) != len(u1.Kq) {
		log.Fatalf("[ConcatEdges] len(u1.Ks): %v != len(u1.Kq): %v\n", len(u1.Ks), len(u1.Kq))
	}
	if len(u2.Ks) != len(u2.Kq) {
		log.Fatalf("[ConcatEdges] len(u2.Ks): %v != len(u2.Kq): %v\n", len(u2.Ks), len(u2.Kq))
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

func substituteEdgeID(nodeMap map[string]DBGNode, nodekey []byte, srcID, dstID DBG_MAX_INT) bool {
	var nkB constructcf.ReadBnt
	nkB.Seq = nodekey
	ks := constructcf.GetReadBntKmer(nkB, 0, Kmerlen-1)
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

}

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

func SmfyDBG(nodesArr []DBGNode, edgesArr []DBGEdge, opt Options) {
	kmerlen := opt.Kmer
	deleteNodeNum, deleteEdgeNum := 0, 0
	longTipsEdgesNum := 0
	for i, v := range nodesArr {
		if i < 1 || v.GetDeleteFlag() > 0 {
			continue
		}
		if IsInComing(v.EdgeIDIncoming, 1) || IsInComing(v.EdgeIDOutcoming, 1) {
			SubstituteEdgeID(nodesArr, v.ID, 1, 0)
			v = nodesArr[v.ID]
		}
		inNum, inID := GetEdgeIDComing(v.EdgeIDIncoming)
		outNum, outID := GetEdgeIDComing(v.EdgeIDOutcoming)
		if inNum == 0 && outNum == 0 {
			nodesArr[i].SetDeleteFlag()
			deleteNodeNum++
		} else if inNum+outNum == 1 {
			longTipsEdgesNum++
			id := inID
			if outNum == 1 {
				id = outID
			}
			if edgesArr[id].StartNID == v.ID {
				edgesArr[id].StartNID = 0
			} else {
				edgesArr[id].EndNID = 0
			}
			nodesArr[i].SetDeleteFlag()
			deleteNodeNum++
		} else if inNum == 1 && outNum == 1 && inID != outID { // prevent cycle ring
			e1 := edgesArr[inID]
			e2 := edgesArr[outID]
			//fmt.Printf("[SmfyDBG] e1: %v\n\te2: %v\n\tnd: %v\n", e1, e2, v)
			u1, u2 := e1.Utg, e2.Utg
			if e1.EndNID == v.ID {
				nID := e2.EndNID
				if e2.StartNID != v.ID {
					u2 = GetRCUnitig(u2)
					nID = e2.StartNID
				}
				edgesArr[inID].Utg = ConcatEdges(u1, u2, kmerlen)
				edgesArr[inID].EndNID = nID
				if nID > 0 && !SubstituteEdgeID(nodesArr, nID, e2.ID, e1.ID) {
					log.Fatalf("[SmfyDBG] e2.ID: %v substitute by e1.ID: %v failed, node: %v\n", e2.ID, e1.ID, nodesArr[nID])
				}
			} else {
				nID := e2.StartNID
				if e2.EndNID != v.ID {
					u2 = GetRCUnitig(u2)
					nID = e2.EndNID
				}
				edgesArr[inID].Utg = ConcatEdges(u2, u1, kmerlen)
				edgesArr[inID].StartNID = nID
				if nID > 0 && !SubstituteEdgeID(nodesArr, nID, e2.ID, e1.ID) {
					log.Fatalf("[SmfyDBG] e2.ID: %v substitute by e1.ID: %v failed, node: %v\n", e2.ID, e1.ID, nodesArr[nID])
				}
			}

			edgesArr[outID].SetDeleteFlag()
			deleteEdgeNum++
			nodesArr[v.ID].SetDeleteFlag()
			deleteNodeNum++
		}
	}

	// delete maybe short repeat edge than small than opt.MaxNGSReadLen
	for i, e := range edgesArr {
		if i < 2 || e.GetDeleteFlag() > 0 {
			continue
		}

		if e.StartNID == 0 && e.EndNID == 0 && len(e.Utg.Ks) < opt.MaxNGSReadLen {
			edgesArr[i].SetDeleteFlag()
			deleteEdgeNum++
		}
		// remove samll cycle maybe repeat
		if e.StartNID > 0 && e.StartNID == e.EndNID && len(e.Utg.Ks) < opt.MaxNGSReadLen {
			v := nodesArr[e.StartNID]
			inNum, inID := GetEdgeIDComing(v.EdgeIDIncoming)
			outNum, outID := GetEdgeIDComing(v.EdgeIDOutcoming)
			if inNum == outNum && inID == outID {
				edgesArr[i].SetDeleteFlag()
				deleteEdgeNum++
				nodesArr[v.ID].SetDeleteFlag()
				deleteNodeNum++
			}
		}
	}

	fmt.Printf("[SmfyDBG]deleted nodes number is : %d\n", deleteNodeNum)
	fmt.Printf("[SmfyDBG]deleted edges number is : %d\n", deleteEdgeNum)
	fmt.Printf("[SmfyDBG]long tips number is : %d\n", longTipsEdgesNum)
}

func transform2QSeq(utg Unitig) alphabet.QLetters {
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

func transform2Letters(ks []byte) alphabet.Letters {
	ls := make(alphabet.Letters, len(ks))
	for i, b := range ks {
		//ls = append(ls, alphabet.Letter(bnt.BitNtCharUp[b]))
		ls[i] = alphabet.Letter(bnt.BitNtCharUp[b])
	}

	return ls
}

func transform2Unitig(Seq alphabet.QLetters, qual bool) (utg Unitig) {

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
func StoreEdgesToFn(edgesfn string, edgesArr []DBGEdge) {
	fp, err := os.Create(edgesfn)
	if err != nil {
		log.Fatalf("[StoreEdgesToFn] create file: %s failed, err: %v\n", edgesfn, err)
	}
	defer fp.Close()

	fqfp := fastq.NewWriter(fp)
	for _, v := range edgesArr {
		if v.ID > 0 && v.GetDeleteFlag() == 0 {
			seq := linear.NewQSeq("", nil, alphabet.DNA, alphabet.Sanger)
			seq.ID = strconv.Itoa(int(v.ID))
			// Add start and end adapter seq
			qs := transform2QSeq(v.Utg)
			seq.AppendQLetters(qs...)
			var path string
			/*if len(v.PathMat) > 0 && len(v.PathMat[0].IDArr) > 1 {
				for _, id := range v.PathMat[0].IDArr {
					path += strconv.Itoa(int(id)) + "-"
				}
				path = path[:len(path)-1]
			}*/
			ans := strconv.Itoa(int(v.StartNID)) + "\t" + strconv.Itoa(int(v.EndNID)) + "\tpath:" + path + "\tlen:" + strconv.Itoa(seq.Len())
			seq.Annotation.SetDescription(ans)
			_, err := fqfp.Write(seq)
			if err != nil {
				log.Fatalf("[StoreEdgesToFn] write seq: %v; err: %v\n", seq, err)
			}
		}
	}
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
				seq1.AppendLetters(transform2Letters(v.Utg.Ks[:MaxMapEdgeLen/2])...)
				seq2.AppendLetters(transform2Letters(v.Utg.Ks[len(v.Utg.Ks)-MaxMapEdgeLen/2:])...)
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
				la := transform2Letters(v.Utg.Ks)
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

func LoadEdgesfqFromFn(fn string, edgesArr []DBGEdge, qual bool) {
	fp, err := os.Open(fn)
	if err != nil {
		log.Fatalf("[LoadEdgesfaFromFn] open file: %s error: %v\n", fn, err)
	}
	defer fp.Close()
	fqfp := fastq.NewReader(fp, linear.NewQSeq("", nil, alphabet.DNA, alphabet.Sanger))
	for {
		if s, err := fqfp.Read(); err != nil {
			if err == io.EOF {
				break
			} else {
				log.Fatalf("[LoadEdgesfqFromFn] read file: %s error: %v\n", fn, err)
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
			var lenKs int
			_, err = fmt.Sscanf(l.Description(), "%v\t%v\t%v\tlen:%d\n", &edge.StartNID, &edge.EndNID, &ps, &lenKs)
			if err != nil {
				log.Fatalf("[LoadEdgesfqFromFn] parse Description:%s of fastq err: %v\n", l.Description(), err)
			}
			if len(ps) > 5 {
				var path Path
				for _, item := range strings.Split(ps[5:], "-") { // ps[:5] == "path:"
					id, err := strconv.Atoi(item)
					if err != nil {
						log.Fatalf("[LoadEdgesFqFromFn] path: %v convert to int err: %v\n", ps, err)
					}
					path.IDArr = append(path.IDArr, DBG_MAX_INT(id))
				}
				edge.PathMat = append(edge.PathMat, path)
			}
			edge.Utg = transform2Unitig(l.Seq, qual)
			if edge.ID >= DBG_MAX_INT(len(edgesArr)) {
				log.Fatalf("[LoadEdgesfqFromFn] edge.ID:%v >= len(edgesArr):%d\n", edge.ID, len(edgesArr))
			} else if edgesArr[edge.ID].ID > 0 {
				log.Fatalf("[LoadEdgesfqFromFn] the position: %v in edgesArr has value:%v\n", edge.ID, edgesArr[edge.ID])
			}
			edgesArr[edge.ID] = edge
		}

	}
}

// set the unique edge of edgesArr
func SetDBGEdgesUniqueFlag(edgesArr []DBGEdge, nodesArr []DBGNode) (uniqueNum int) {
	for i, e := range edgesArr {
		if e.ID < 2 || e.GetDeleteFlag() > 0 {
			continue
		}

		var ea1, ea2 []DBG_MAX_INT
		if e.StartNID > 0 {
			nd := nodesArr[e.StartNID]
			ea1 = GetNearEdgeIDArr(nd, e.ID)
		}
		if e.EndNID > 0 {
			nd := nodesArr[e.EndNID]
			ea2 = GetNearEdgeIDArr(nd, e.ID)
		}

		if len(ea1) <= 1 && len(ea2) <= 1 {
			edgesArr[i].SetUniqueFlag()
			uniqueNum++
		}
	}
	return uniqueNum
}

type ReadInfo struct {
	ID  int64
	Seq []byte
}
type AlignInfo struct {
	ID    int64
	Paths []DBG_MAX_INT
}

func paraLoadNGSReads(fn string, cs chan ReadInfo, kmerLen int, we chan int) {
	fp, err := os.Open(fn)
	if err != nil {
		log.Fatalf("[paraLoadNGSReads] %v\n", err)
	}

	// brfp, err := brotli.NewReader(fp, nil)
	gzfp, err := gzip.NewReader(fp)
	if err != nil {
		log.Fatalf("[paraLoadNGSReads] %v\n", err)
	}
	defer gzfp.Close()
	// defer fp.Close()
	// buffp := bufio.NewReader(gzfp)
	fqfp := fastq.NewReader(gzfp, linear.NewQSeq("", nil, alphabet.DNA, alphabet.Sanger))
	s, err := fqfp.Read()
	var count int
	for ; err == nil; s, err = fqfp.Read() {
		var ri ReadInfo
		fq := s.(*linear.QSeq)
		if len(fq.Seq) < kmerLen+10 { // read length is short for found DBG paths
			continue
		}
		id, err := strconv.Atoi(fq.Name())
		if err != nil {
			log.Fatalf("[paraLoadNGSReads] load fn: '%v' file, read ID: %v not digits, please convert to digits...\n", fn, fq.Name())
		}
		ri.ID = int64(id)
		ri.Seq = transform2Unitig(fq.Seq, false).Ks
		cs <- ri
		count++
	}
	if err != io.EOF {
		log.Fatalf("[LoadNGSReads] Failed to read file %v, err: %v\n", fn, err)
	}
	we <- count

}

func LoadNGSReads(cfgFn string, cs chan ReadInfo, numCPU, kmerLen int) {
	cfgInfo, err := constructcf.ParseCfg(cfgFn)
	if err != nil {
		log.Fatal("[constructcf.ParseCfg] found err")
	}
	fmt.Println(cfgInfo)

	// iterate cfgInfo find fastq files
	we := make(chan int)
	var numT int
	for _, lib := range cfgInfo.Libs {
		// seqProfile == 1 note Illumina
		if lib.AsmFlag != constructcf.AllState && lib.SeqProfile != 1 {
			continue
		}
		for _, fn := range lib.FnName {
			go paraLoadNGSReads(fn, cs, kmerLen, we)
			numT++
		}
	}

	// check child routinue have finishedT
	var totalNumReads int
	for i := 0; i < numT; i++ {
		totalNumReads += <-we
	}
	fmt.Printf("[LoadNGSReads] total processed number reads : %v\n", totalNumReads)

	// send map goroutinues finish signal
	close(cs)
}

func writeAlignToFile(wrFn string, wc chan AlignInfo, numCPU int) {
	fp, err := os.Create(wrFn)
	if err != nil {
		log.Fatalf("[writeAlignToFile] failed to create file: %s, err:%v\n", wrFn, err)
	}
	defer fp.Close()

	var finishedT int
	for {
		ai := <-wc
		if ai.ID < 0 {
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
		fp.WriteString(s)
	}
}

func BiggerThan(kb, rb []byte) bool {
	for i, b := range kb {
		if b > rb[i] {
			return true
		} else if b < rb[i] {
			return false
		}
	}
	return false
}

func getCuckoofilterDBGSampleSize(edgesArr []DBGEdge, winSize, maxNGSReadLen, kmerLen int) (cfSize int64) {

	for _, e := range edgesArr[2:] {
		if e.GetDeleteFlag() > 0 {
			continue
		}
		el := len(e.Utg.Ks)
		if el < 2*maxNGSReadLen-kmerLen { // get whole length sliding windowns
			cfSize += int64((el-kmerLen+winSize-1)/winSize) + 1
		} else { // just sample two ends
			cfSize += int64((maxNGSReadLen-kmerLen+winSize-1)/winSize+1) * 2
		}
	}
	return cfSize
}

func GetMinimizer(seq []byte, kmerlen int) (minSeq []byte, pos uint32, strand bool) {
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
			pos = uint32(i)
			strand = st
		}
	}

	return
}

func constructCFDBGMinimizers(cf CuckooFilter, edgesArr []DBGEdge, winSize, maxNGSReadLen int) (count int) {
	for _, e := range edgesArr {
		if e.ID < 2 || e.GetDeleteFlag() > 0 {
			continue
		}
		el := len(e.Utg.Ks)
		maxLen := el
		if el > 2*maxNGSReadLen-cf.Kmerlen {
			maxLen = maxNGSReadLen
		}
		//fmt.Printf("[constructCFDBGSample] edge length: %v\n", len(e.Utg.Ks))
		for j := 0; j < maxLen-cf.Kmerlen+1; j += winSize {
			z := j + cf.Kmerlen + winSize - 1
			if z > maxLen { // make the boundary kmer in the Samples
				z = maxLen
			}
			//fmt.Printf("el: %v\tmaxGNSReadLen: %v\tcf.Kmerlen: %v\tmaxLen: %v\tj: %v\t", el, maxNGSReadLen, cf.Kmerlen, maxLen, j)
			kb, pos, strand := GetMinimizer(e.Utg.Ks[j:z], cf.Kmerlen)

			suc := cf.Insert(kb, e.ID, uint32(j)+pos, strand)
			if suc == false {
				log.Fatalf("[constructCFDBGSample] Insert to the CuckooFilter of DBGSample false\n")
			}
			count++
		}

		if el > 2*maxNGSReadLen-cf.Kmerlen {
			for j := el - maxNGSReadLen; j < el-cf.Kmerlen+1; j += winSize {
				z := j + cf.Kmerlen + winSize - 1
				if z > el { // make the boundary kmer in the Samples
					z = el
				}
				kb, pos, strand := GetMinimizer(e.Utg.Ks[j:z], cf.Kmerlen)

				suc := cf.Insert(kb, e.ID, uint32(j)+pos, strand)
				if suc == false {
					log.Fatalf("[constructCFDBGSample] Insert to the CuckooFilter of DBGSample false\n")
				}
				count++
			}
		}
		//fmt.Printf("[constructCFDBGSample] count : %v\n", count)
	}

	return count
}

// found kmer seed position in the DBG edges
func LocateSeedKmerCF(cf CuckooFilter, ri ReadInfo, winSize int, edgesArr []DBGEdge) (dbgK DBGKmer, pos uint32, strand bool) {
	if len(ri.Seq) < cf.Kmerlen+2*winSize {
		fmt.Printf("[LocateSeedKmerCF] read: %v sequence length smaller than KmerLen(%v) + 2 * winSize(%v)\n", ri, cf.Kmerlen, winSize)
	}
	var kb []byte
	count := 0
	for i := 0; i < 2; i++ {
		kb, pos, strand = GetMinimizer(ri.Seq[i*winSize:i*winSize+cf.Kmerlen+winSize-1], cf.Kmerlen)
		dbgK, count = cf.Lookup(kb, edgesArr)
		pos += uint32(i * winSize)
		if count == 1 && dbgK.GetCount() > 0 {
			return
		}
		if count > 1 {
			//fmt.Printf("[LocateSeedKmerCF] found seed count: %v\n", count)
			dbgK.setCFItem(0, 0)
			break
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

// parallel Map NGS reads to the DBG edges, then output alignment path for the DBG
func paraMapNGS2DBG(cs chan ReadInfo, wc chan AlignInfo, nodesArr []DBGNode, edgesArr []DBGEdge, cf CuckooFilter, winSize int) {
	var notFoundSeedNum, mapOneEdgeNum, notPerfectNum int
	for {
		ri, ok := <-cs
		var ai AlignInfo
		if !ok {
			ai.ID = -1
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
				if dbgK.Pos > pos && len(ri.Seq)-int(pos) < el-int(dbgK.Pos) {
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
			ai.ID = ri.ID
			overAll := true // note if overall length read sequence match
			// map the start partition of read sequence
			{
				dk := dbgK
				rpos := int(pos)
				if dk.Strand != strand {
					dk.Pos += uint32(cf.Kmerlen)
				}
				for i := rpos - 1; i >= 0; {
					e := edgesArr[dk.ID]
					b := 0
					if dk.Strand == strand {
						if int(dk.Pos) < rpos {
							b = rpos - int(dk.Pos)
						}
						j := int(dk.Pos) - 1
						for ; i >= b; i-- {
							//fmt.Printf("[paraMapNGS2DBG] i: %v, j: %v, dk.Pos: %v, pos: %v, b: %v\n", i, j, dk.Pos, rpos, b)
							if ri.Seq[i] != e.Utg.Ks[j] {
								overAll = false
								break
							}
							j--
						}
					} else { // dk.Strand != strand
						if len(e.Utg.Ks)-int(dk.Pos) < rpos {
							b = rpos - (len(e.Utg.Ks) - int(dk.Pos))
						}
						j := int(dk.Pos)
						for ; i >= b; i-- {
							//fmt.Printf("[paraMapNGS2DBG] i: %v, j: %v, dk.Pos: %v, pos: %v, b: %v\n", i, j, dk.Pos, rpos, b)
							if ri.Seq[i] != bnt.BntRev[e.Utg.Ks[j]] {
								overAll = false
								break
							}
							j++
						}
					}

					if !overAll {
						fmt.Printf("[paraMapNGS2DBG]not perfect start i: %v,edge ID: %v,len(e.Utg.Ks):%v,  dk.Pos: %v, pos: %v, b: %v\n", i, dk.ID, len(e.Utg.Ks), dk.Pos, rpos, b)
						break
					}
					ai.Paths = append(ai.Paths, e.ID)

					if i < 0 {
						break
					}
					// find next edge
					rpos = i + 1
					if dk.Strand == strand {
						if e.StartNID == 0 {
							break
						}
						n := nodesArr[e.StartNID]
						if IsInComing(n.EdgeIDOutcoming, e.ID) {
							if n.EdgeIDIncoming[ri.Seq[i]] > 1 {
								dk.ID = n.EdgeIDIncoming[ri.Seq[i]]
								e = edgesArr[dk.ID]
								if e.EndNID == n.ID {
									dk.Pos = uint32(len(e.Utg.Ks) - (cf.Kmerlen - 1))
								} else {
									dk.Strand = !dk.Strand
									dk.Pos = uint32(cf.Kmerlen) - 1
								}
							} else {
								fmt.Printf("[paraMapNGS2DBG] not found next edge in node: %v\n", n)
								break
							}
						} else {
							b := bnt.BntRev[ri.Seq[i]]
							if n.EdgeIDOutcoming[b] > 1 {
								dk.ID = n.EdgeIDOutcoming[b]
								e = edgesArr[dk.ID]
								if e.StartNID == n.ID {
									dk.Strand = !dk.Strand
									dk.Pos = uint32(cf.Kmerlen) - 1
								} else {
									dk.Pos = uint32(len(e.Utg.Ks) - (cf.Kmerlen - 1))
								}
							} else {
								fmt.Printf("[paraMapNGS2DBG] not found next edge in node: %v\n", n)
								break
							}
						}
					} else {
						if e.EndNID == 0 {
							break
						}
						n := nodesArr[e.EndNID]
						if IsInComing(n.EdgeIDIncoming, e.ID) {
							b := bnt.BntRev[ri.Seq[i]]
							if n.EdgeIDOutcoming[b] > 1 {
								dk.ID = n.EdgeIDOutcoming[b]
								e = edgesArr[dk.ID]
								if e.StartNID == n.ID {
									dk.Pos = uint32(cf.Kmerlen) - 1
								} else {
									dk.Strand = !dk.Strand
									dk.Pos = uint32(len(e.Utg.Ks) - (cf.Kmerlen - 1))
								}
							} else {
								fmt.Printf("[paraMapNGS2DBG] not found next edge in node: %v\n", n)
								break
							}
						} else { // e.ID in the n.EdgeIDOutcoming
							if n.EdgeIDIncoming[ri.Seq[i]] > 1 {
								dk.ID = n.EdgeIDIncoming[ri.Seq[i]]
								e = edgesArr[dk.ID]
								if e.EndNID == n.ID {
									dk.Strand = !dk.Strand
									dk.Pos = uint32(len(e.Utg.Ks) - (cf.Kmerlen - 1))
								} else {
									dk.Pos = uint32(cf.Kmerlen) - 1
								}
							} else {
								fmt.Printf("[paraMapNGS2DBG] not found next edge in node: %v\n", n)
								break
							}
						}
					}
				}
			}
			if !overAll {
				notPerfectNum++
				continue
			}
			ai.Paths = ReverseDBG_MAX_INTArr(ai.Paths)
			// map the end partition of read sequence
			var path []DBG_MAX_INT
			{
				dk := dbgK
				rpos := int(pos) + cf.Kmerlen
				if dk.Strand == strand {
					dk.Pos += uint32(cf.Kmerlen)
				}
				for i := rpos; i < len(ri.Seq); {
					e := edgesArr[dk.ID]
					b := len(ri.Seq)
					if dk.Strand == strand {
						if len(e.Utg.Ks)-int(dk.Pos) < len(ri.Seq)-rpos {
							b = rpos + (len(e.Utg.Ks) - int(dk.Pos))
						}
						j := int(dk.Pos)
						for ; i < b; i++ {
							//fmt.Printf("[paraMapNGS2DBG] i: %v, j: %v, dk.Pos: %v, pos: %v, b: %v\n", i, j, dk.Pos, rpos, b)
							if ri.Seq[i] != e.Utg.Ks[j] {
								overAll = false
								break
							}
							j++
						}
					} else { // dbgK.Strand != strand
						if len(ri.Seq)-rpos > int(dk.Pos) {
							b = rpos + int(dk.Pos)
						}
						j := int(dk.Pos) - 1
						for ; i < b; i++ {
							//fmt.Printf("[paraMapNGS2DBG] i: %v, j: %v, dk.Pos: %v, pos: %v, b: %v\n", i, j, dk.Pos, rpos, b)
							if ri.Seq[i] != bnt.BntRev[e.Utg.Ks[j]] {
								overAll = false
								break
							}
							j--
						}
					}

					if overAll == false {
						fmt.Printf("[paraMapNGS2DBG]not perfect end i: %v,edge ID: %v,len(e.Utg.Ks): %v,  dk.Pos: %v, pos: %v, b: %v\n", i, dk.ID, len(e.Utg.Ks), dk.Pos, rpos, b)
						break
					}
					path = append(path, e.ID)

					if i == len(ri.Seq) {
						break
					}

					// find next edge
					rpos = i
					if dk.Strand == strand {
						if e.EndNID == 0 {
							break
						}
						n := nodesArr[e.EndNID]
						if IsInComing(n.EdgeIDIncoming, e.ID) {
							if n.EdgeIDOutcoming[ri.Seq[i]] > 1 {
								dk.ID = n.EdgeIDOutcoming[ri.Seq[i]]
								e = edgesArr[dk.ID]
								if e.StartNID == n.ID {
									dk.Pos = uint32(cf.Kmerlen) - 1
								} else {
									dk.Strand = !dk.Strand
									dk.Pos = uint32(len(e.Utg.Ks) - (cf.Kmerlen - 1))
								}
							} else {
								fmt.Printf("[paraMapNGS2DBG] not found next edge in node: %v\n", n)
								break
							}
						} else { // e.ID in the n.EdgeIDOutcoming
							base := bnt.BntRev[ri.Seq[i]]
							if n.EdgeIDIncoming[base] > 1 {
								dk.ID = n.EdgeIDIncoming[base]
								e = edgesArr[dk.ID]
								if e.EndNID == n.ID {
									dk.Strand = !dk.Strand
									dk.Pos = uint32(len(e.Utg.Ks) - (cf.Kmerlen - 1))
								} else {
									dk.Pos = uint32(cf.Kmerlen) - 1
								}
							} else {
								fmt.Printf("[paraMapNGS2DBG] not found next edge in node: %v\n", n)
								break
							}
						}
					} else { // dk.Strand != strand
						if e.StartNID == 0 {
							break
						}
						n := nodesArr[e.StartNID]
						if IsInComing(n.EdgeIDOutcoming, e.ID) {
							base := bnt.BntRev[ri.Seq[i]]
							if n.EdgeIDIncoming[base] > 1 {
								dk.ID = n.EdgeIDIncoming[base]
								e = edgesArr[dk.ID]
								if e.EndNID == n.ID {
									dk.Pos = uint32(len(e.Utg.Ks) - (cf.Kmerlen - 1))
								} else {
									dk.Strand = !dk.Strand
									dk.Pos = uint32(cf.Kmerlen - 1)
								}
							} else {
								fmt.Printf("[paraMapNGS2DBG] not found next edge in node: %v\n", n)
								break
							}
						} else { // e.ID in the n.EdgeIDIncoming
							if n.EdgeIDOutcoming[ri.Seq[i]] > 1 {
								dk.ID = n.EdgeIDOutcoming[ri.Seq[i]]
								e = edgesArr[dk.ID]
								if e.StartNID == n.ID {
									dk.Strand = !dk.Strand
									dk.Pos = uint32(cf.Kmerlen - 1)
								} else {
									dk.Pos = uint32(len(e.Utg.Ks) - (cf.Kmerlen - 1))
								}
							} else {
								fmt.Printf("[paraMapNGS2DBG] not found next edge in node: %v\n", n)
								break
							}
						}
					}
				}
			}

			if !overAll {
				//fmt.Printf("[paraMapNGS2DBG] ", a)
				notPerfectNum++
				continue
			}

			if len(ai.Paths) > 0 && len(path) > 0 && ai.Paths[len(ai.Paths)-1] == path[0] {
				ai.Paths = append(ai.Paths, path[1:]...)
			} else {
				ai.Paths = append(ai.Paths, path...)
			}

			// write to output
			if overAll && len(ai.Paths) > 1 {
				wc <- ai
			}
		}
	}
	fmt.Printf("[paraMapNGS2DBG] not found seed reads number is : %v\n", notFoundSeedNum)
	fmt.Printf("[paraMapNGS2DBG] map one edge reads number is : %v\n", mapOneEdgeNum)
	fmt.Printf("[paraMapNGS2DBG] not perfect mapping reads number is : %v\n", notPerfectNum)

	/*ai.ID = ri.ID
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
	} */
}

func MapNGS2DBG(opt Options, nodesArr []DBGNode, edgesArr []DBGEdge, wrFn string) {
	// construct cuckoofilter of DBG sample
	cfSize := getCuckoofilterDBGSampleSize(edgesArr, opt.WinSize, opt.MaxNGSReadLen, opt.Kmer)
	fmt.Printf("[MapNGS2DBG] cfSize: %v\n", cfSize)
	cf := MakeCuckooFilter(uint64(cfSize*3), opt.Kmer)
	fmt.Printf("[MapNGS2DBG] cf.numItems: %v\n", cf.numItems)
	count := constructCFDBGMinimizers(cf, edgesArr, opt.WinSize, opt.MaxNGSReadLen)
	fmt.Printf("[MapNGS2DBG]construct Smaple of DBG edges cuckoofilter number is : %v\n", count)
	//if cfSize != int64(count) {
	//	log.Fatalf("[MapNGS2DBG]cfSize : %v != count : %v, please check\n", cfSize, count)
	//}

	numCPU := opt.NumCPU
	runtime.GOMAXPROCS(numCPU + 2)
	bufSize := numCPU * 20
	cs := make(chan ReadInfo, bufSize)
	wc := make(chan AlignInfo, bufSize)
	defer close(wc)
	//defer close(cs)

	// Load NGS read from cfg
	fn := opt.CfgFn
	go LoadNGSReads(fn, cs, numCPU, opt.Kmer)
	for i := 0; i < numCPU; i++ {
		go paraMapNGS2DBG(cs, wc, nodesArr, edgesArr, cf, opt.WinSize)
	}

	// write function
	writeAlignToFile(wrFn, wc, numCPU)
}

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
	buffp := bufio.NewReader(fp)
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

func GetNearEdgeIDArr(nd DBGNode, eID DBG_MAX_INT) (eArr []DBG_MAX_INT) {
	if eID <= 0 {
		log.Fatalf("[GetNearEdgeIDArr] eID must bigger than zero, eID: %v\n", eID)
	}
	in := false
	for _, id := range nd.EdgeIDIncoming {
		if id == eID {
			in = true
			break
		}
	}
	if in {
		for _, id := range nd.EdgeIDOutcoming {
			if id > 0 {
				eArr = append(eArr, id)
			}
		}
	} else {
		for _, id := range nd.EdgeIDIncoming {
			if id > 0 {
				eArr = append(eArr, id)
			}
		}
	}

	return eArr
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
func MergePathMat(edgesArr []DBGEdge, nodesArr []DBGNode, minMapFreq int) {
	for i, e := range edgesArr {
		//if e.GetUniqueFlag() > 0 {
		//	fmt.Printf("[MergePathMat]unique edge : %v\n", e)
		//}
		if i < 2 || e.GetDeleteFlag() > 0 || len(e.PathMat) == 0 || e.GetUniqueFlag() == 0 {
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
			}
			continue
		}

		CheckPathDirection(edgesArr, nodesArr, e.ID)
		if edgesArr[e.ID].GetUniqueFlag() == 0 {
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
		/*for j, p := range pm {
			fmt.Printf("[MergePathMat] pm[%v]: %v\n", j, p)
		} */

		// find consis Path
		var path Path
		path.Freq = math.MaxInt64
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
}

func IsTwoEdgesCyclePath(edgesArr []DBGEdge, nodesArr []DBGNode, eID DBG_MAX_INT) bool {
	e := edgesArr[eID]
	if e.StartNID == 0 || e.EndNID == 0 {
		return false
	}
	arr1 := GetNearEdgeIDArr(nodesArr[e.StartNID], eID)
	arr2 := GetNearEdgeIDArr(nodesArr[e.EndNID], eID)
	if len(arr1) == 1 && len(arr2) == 1 && arr1[0] == arr2[0] {
		return true
	}

	return false
}

func ExtendPath(edgesArr []DBGEdge, nodesArr []DBGNode, e DBGEdge, minMappingFreq int) (maxP Path) {
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
		if e2.GetUniqueFlag() == 0 || len(e2.PathMat) != 1 || e2.PathMat[0].Freq < minMappingFreq {
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
		if j1 < j2 {
			p2.IDArr = GetReverseDBG_MAX_INTArr(p2.IDArr)
			j1 = IndexEID(p2.IDArr, e1.ID)
			j2 = IndexEID(p2.IDArr, e2.ID)
		}
		k1 := IndexEID(p.IDArr, e1.ID)
		fmt.Printf("[ExtendPath]j1: %v, j2:%v,k1: %v, eID1: %v, eID2: %v, p: %v, p2: %v\n", j1, j2, k1, e1.ID, e2.ID, p, p2)
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
				i = 0
			}
		}
	}

	ReverseDBG_MAX_INTArr(mutualArr)
	idx = IndexEID(p.IDArr, e.ID)
	// find right path
	for i := len(p.IDArr) - 1; i > idx; i-- {
		e2 := edgesArr[p.IDArr[i]]
		if e2.GetUniqueFlag() == 0 || len(e2.PathMat) != 1 || e2.PathMat[0].Freq < minMappingFreq {
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
			j2 = IndexEID(p2.IDArr, e2.ID)
		}
		k1 := IndexEID(p.IDArr, e1.ID)
		fmt.Printf("[ExtendPath] j1: %v, j2: %v, k1: %v, eID1: %v, eID2: %v, p: %v, p2: %v\n", j1, j2, k1, e1.ID, e2.ID, p, p2)
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
				i = len(p.IDArr) - 1
			}
		}
	}

	if len(mutualArr) <= 1 {
		edgesArr[e.ID].SetProcessFlag()
		return
	}

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
	//fmt.Printf("[ExtendPath] mutualArr: %v\n", mutualArr)
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
func findMaxPath(edgesArr []DBGEdge, nodesArr []DBGNode, minMapFreq int) (pathArr []Path) {
	for i, e := range edgesArr {
		if i < 2 || e.GetDeleteFlag() > 0 || e.GetProcessFlag() > 0 || e.GetUniqueFlag() == 0 || len(e.PathMat) != 1 {
			continue
		}

		p := e.PathMat[0]
		if p.Freq < minMapFreq || len(p.IDArr) < 2 {
			continue
			//log.Fatalf("[findMaxPath] edgesArr[%v].PathMat: %v\t not contain useful info\n", i, p)
		}

		maxP := ExtendPath(edgesArr, nodesArr, e, minMapFreq)
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

func CascadePath(p Path, edgesArr []DBGEdge, nodesArr []DBGNode, kmerlen int) {
	p0 := edgesArr[p.IDArr[0]]
	p1 := edgesArr[p.IDArr[1]]
	//eStartNID, eEndNID := p0.StartNID, p0.EndNID
	var direction uint8
	nID := GetLinkNodeID(p0, p1)
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

func ResetMergedEdgeNID(edgesArr []DBGEdge, nodesArr []DBGNode, joinPathArr []Path, IDMapPath map[DBG_MAX_INT]uint32) {
	for _, e := range edgesArr {
		if e.ID == 0 {
			continue
		}
		if idx, ok := IDMapPath[e.ID]; ok {
			if e.ID == joinPathArr[idx].IDArr[0] {
				if e.StartNID > 0 {
					//SubstituteEdgeID(nodesArr, e.StartNID, e.ID, 0)
					if !IsInDBGNode(nodesArr[e.StartNID], e.ID) {
						var sBnt constructcf.ReadBnt
						sBnt.Seq = e.Utg.Ks[:Kmerlen-1]
						ks := constructcf.GetReadBntKmer(sBnt, 0, Kmerlen-1)
						ts := ks
						rs := constructcf.ReverseComplet(ks)
						if ks.BiggerThan(rs) {
							ks, rs = rs, ks
						}
						if reflect.DeepEqual(ks.Seq, nodesArr[e.StartNID].Seq) == false {
							log.Fatalf("[ResetMergedEdgeNID] ks != nodesArr[e.StartNID].Seq, ks: %v\n\tnodesArr[%v]: %v\n", ks, e.StartNID, nodesArr[e.StartNID])
						}
						c := e.Utg.Ks[Kmerlen-1]
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
						var sBnt constructcf.ReadBnt
						sBnt.Seq = e.Utg.Ks[e.GetSeqLen()-Kmerlen+1:]
						ks := constructcf.GetReadBntKmer(sBnt, 0, Kmerlen-1)
						ts := ks
						rs := constructcf.ReverseComplet(ks)
						if ks.BiggerThan(rs) {
							ks, rs = rs, ks
						}
						if reflect.DeepEqual(ks.Seq, nodesArr[e.EndNID].Seq) == false {
							log.Fatalf("[ResetMergedEdgeNID] ks != nodesArr[e.StartNID].Seq, ks: %v\n\tnodesArr[%v]: %v\n", ks, e.StartNID, nodesArr[e.StartNID])
						}
						c := e.Utg.Ks[e.GetSeqLen()-Kmerlen]
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
}

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
		arr := GetNearEdgeIDArr(nodesArr[e.StartNID], e.ID)
		num := len(arr)
		arr = GetNearEdgeIDArr(nodesArr[e.EndNID], e.ID)
		num += len(arr)
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

		CascadePath(p, edgesArr, nodesArr, kmerlen)
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

func CheckPathDirection(edgesArr []DBGEdge, nodesArr []DBGNode, eID DBG_MAX_INT) {
	e := edgesArr[eID]
	var ea1, ea2 []DBG_MAX_INT
	if e.StartNID > 0 {
		ea1 = GetNearEdgeIDArr(nodesArr[e.StartNID], eID)
	}
	if e.EndNID > 0 {
		ea2 = GetNearEdgeIDArr(nodesArr[e.EndNID], eID)
	}

	// found two edge cycle
	if len(ea1) == 1 && len(ea2) == 1 && ea1[0] == ea2[0] {
		edgesArr[eID].ResetUniqueFlag()
		return
	}

	for i, p := range e.PathMat {
		idx := IndexEID(p.IDArr, eID)
		if idx < len(p.IDArr)-1 {
			if IsInDBG_MAX_INTArr(ea1, p.IDArr[idx+1]) {
				ReverseDBG_MAX_INTArr(edgesArr[eID].PathMat[i].IDArr)
			}
		} else {
			if IsInDBG_MAX_INTArr(ea2, p.IDArr[idx-1]) {
				ReverseDBG_MAX_INTArr(edgesArr[eID].PathMat[i].IDArr)
			}
		}
	}

	//fmt.Printf("[CheckPathDirection] ea1: %v, ea2: %v, PathMat: %v\n", ea1, ea2, edgesArr[eID].PathMat)

}

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

func AdjustPathMat(edgesArr []DBGEdge, nodesArr []DBGNode, joinPathArr []Path, IDMapPath map[DBG_MAX_INT]uint32) {

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
}

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

func SimplifyByNGS(opt Options, nodesArr []DBGNode, edgesArr []DBGEdge, mapNGSFn string) {
	// add to the DBGEdge pathMat
	AddPathToDBGEdge(edgesArr, mapNGSFn)

	// merge pathMat
	MergePathMat(edgesArr, nodesArr, opt.MinMapFreq)

	// find maximum path
	joinPathArr := findMaxPath(edgesArr, nodesArr, opt.MinMapFreq)
	/*i := 125
	if edgesArr[i].GetDeleteFlag() == 0 {
		log.Fatalf("[SimplifyByNGS]edgesArr[%v]: %v\n", i, edgesArr[i])
	}*/
	// constuct map edge ID to the path
	//DeleteJoinPathArrEnd(edgesArr, joinPathArr)
	ReconstructDBG(edgesArr, nodesArr, joinPathArr, opt.Kmer)
	// debug code
	//graphfn := opt.Prefix + ".afterNGS.dot"
	//GraphvizDBGArr(nodesArr, edgesArr, graphfn)

	//AdjustPathMat(edgesArr, nodesArr, joinPathArr, IDMapPath)
}

func CheckInterConnectivity(edgesArr []DBGEdge, nodesArr []DBGNode) {
	for _, e := range edgesArr {
		if e.ID < 2 || e.GetDeleteFlag() > 0 {
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
	//MaxMapEdgeLen int // max length of edge that don't need cut two flank sequence to map Long Reads
}

func checkArgs(c cli.Command) (opt Options, succ bool) {

	var ok bool
	opt.TipMaxLen, ok = c.Flag("tipMaxLen").Get().(int)
	if !ok {
		log.Fatalf("[checkArgs] argument 'tipMaxLen': %v set error\n ", c.Flag("tipMaxlen").String())
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
	opt := Options{gOpt, 0, 0, 0, 0}
	tmp, suc := checkArgs(c)
	if suc == false {
		log.Fatalf("[Smfy] check Arguments error, opt: %v\n", tmp)
	}
	opt.TipMaxLen = tmp.TipMaxLen
	opt.MaxNGSReadLen = tmp.MaxNGSReadLen
	opt.WinSize = tmp.WinSize
	opt.MinMapFreq = tmp.MinMapFreq
	//opt.MaxMapEdgeLen = tmp.MaxMapEdgeLen
	fmt.Printf("Arguments: %v\n", opt)

	// set package-level variable
	Kmerlen = opt.Kmer

	// read nodes file and transform to array mode for more quckly access
	nodesfn := opt.Prefix + ".nodes.mmap"
	nodeMap := NodeMapMmapReader(nodesfn)
	DBGStatfn := opt.Prefix + ".DBG.stat"
	nodesSize, edgesSize := DBGStatReader(DBGStatfn)
	//nodesSize := len(nodeMap)
	fmt.Printf("[Smfy] len(nodeMap): %v, length of edge array: %v\n", nodesSize, edgesSize)
	// read edges file
	edgesfn := opt.Prefix + ".edges"
	edgesArr := ReadEdgesFromFile(edgesfn, edgesSize)
	gfn1 := opt.Prefix + ".beforeSmfyDBG.dot"
	GraphvizDBG(nodeMap, edgesArr, gfn1)

	//gfn := opt.Prefix + ".smfyDBG.dot"
	//GraphvizDBG(nodeMap, edgesArr, gfn)
	// reconstruct consistence De Bruijn Graph
	// ReconstructConsistenceDBG(nodeMap, edgesArr)

	nodesArr := make([]DBGNode, nodesSize)
	NodeMap2NodeArr(nodeMap, nodesArr)
	nodeMap = nil // nodeMap any more used

	// set the unique edge of edgesArr
	uniqueNum := SetDBGEdgesUniqueFlag(edgesArr, nodesArr)
	//CheckInterConnectivity(edgesArr, nodesArr)
	fmt.Printf("[Smfy] the number of DBG Unique  Edges is : %d\n", uniqueNum)

	t1 := time.Now()
	// map Illumina reads to the DBG and find reads map path for simplify DBG
	wrFn := opt.Prefix + ".smfy.NGSAlignment"
	MapNGS2DBG(opt, nodesArr, edgesArr, wrFn)
	//CheckInterConnectivity(edgesArr, nodesArr)
	SimplifyByNGS(opt, nodesArr, edgesArr, wrFn)

	CheckInterConnectivity(edgesArr, nodesArr)

	// simplify DBG
	//IDMapPath := ConstructIDMapPath(joinPathArr)
	SmfyDBG(nodesArr, edgesArr, opt)
	//AdjustPathMat(edgesArr, nodesArr, joinPathArr, IDMapPath)

	t2 := time.Now()
	fmt.Printf("[Smfy] total used : %v, MapNGS2DBG used : %v\n", t2.Sub(t0), t2.Sub(t1))
	graphfn := opt.Prefix + ".afterNGS.dot"
	GraphvizDBGArr(nodesArr, edgesArr, graphfn)

	// Debug code
	// for i := 1; i < len(edgesArr); i++ {
	// 	fmt.Printf("[Smfy]edgesArr[%d]: %v\n", i, edgesArr[i])
	// 	if len(edgesArr[i].Utg.Ks) == len(edgesArr[i-1].Utg.Ks) {
	// 		fmt.Printf("[Smfy] outcoming: %v, %v, len: %d\n", edgesArr[i-1].Utg.Ks[Kmerlen-1], edgesArr[i].Utg.Ks[Kmerlen-1], len(edgesArr[i].Utg.Ks))
	// 	}
	// }
	// Debug code
	// output graphviz graph

	smfyEdgesfn := opt.Prefix + ".edges.smfy.fq"
	StoreEdgesToFn(smfyEdgesfn, edgesArr)
	//mappingEdgefn := opt.Prefix + ".edges.mapping.fa"
	// StoreMappingEdgesToFn(mappingEdgefn, edgesArr, opt.MaxMapEdgeLen)
	//	adpaterEdgesfn := prefix + ".edges.adapter.fq"
	//	StoreEdgesToFn(adpaterEdgesfn, edgesArr, true)
	smfyNodesfn := opt.Prefix + ".nodes.smfy.Arr"
	NodesArrWriter(nodesArr, smfyNodesfn)
	DBGInfofn := opt.Prefix + ".smfy.DBGInfo"
	DBGInfoWriter(DBGInfofn, len(edgesArr), len(nodesArr))
	//EdgesStatWriter(edgesStatfn, len(edgesArr))
}
