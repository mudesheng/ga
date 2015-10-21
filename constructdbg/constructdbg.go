package constructdbg

import (
	"bufio"
	"compress/gzip"
	"container/list"
	"encoding/binary"
	"encoding/gob"
	"fmt"
	"ga/bnt"
	"ga/constructcf"
	"ga/cuckoofilter"
	"io"
	"log"
	"os"
	"runtime"
	"strconv"
	"strings"
	// "time"

	"github.com/awalterschulze/gographviz"
	"github.com/biogo/biogo/alphabet"
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

func (n *DBGNode) GetProcessFlag() uint8 {
	return n.Flag & 0x1
}

func (n *DBGNode) SetProcessFlag(f uint8) {
	// n.Flag = n.Flag & 0xFE
	if f == 1 {
		n.Flag |= 0x1
	} else if f != 0 {
		log.Fatalf("[SetProcessFlag] input flag error: %d\n", f)
	}
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
	Flag     uint8
	Utg      Unitig
	PathMat  []Path // read Path matrix
}

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
}

func (e *DBGEdge) InsertPathToEdge(path []DBG_MAX_INT, freq int) {

	// check have been added
	added := false
	for i, v := range e.PathMat {
		if EqualDBG_MAX_INT(v.IDArr, path) {
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

func InsertPathToEdge(pathMat []Path, path []DBG_MAX_INT, freq int) []Path {

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
}

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

func (e *DBGEdge) SetProcessFlag() {
	e.Flag = e.Flag | 0x1
}

func (e *DBGEdge) ResetProcessFlag() {
	e.Flag = e.Flag & 0xFE
}

func (e *DBGEdge) SetDeleteFlag() {
	e.Flag = e.Flag | 0x2
}

var MIN_KMER_COUNT uint16 = 2
var BACKWARD uint8 = 1
var FORWARD uint8 = 2
var Kmerlen int

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
	defer uniqkmergzfp.Close()
	ukgzfp, err := gzip.NewReader(uniqkmergzfp)
	if err != nil {
		log.Fatal(err)
	}
	defer ukgzfp.Close()
	ukbuffp := bufio.NewReader(ukgzfp)
	var processNumKmer int
	// var rsb constructcf.ReadSeqBucket
	eof := false
	KBntByteNum := (kmerlen + bnt.NumBaseInByte - 1) / bnt.NumBaseInByte
	for !eof {
		//var kmer []byte
		b := make([]byte, KBntByteNum)
		err := binary.Read(ukbuffp, binary.LittleEndian, b)
		if err != nil {
			if err == io.EOF {
				err = nil
				eof = true
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
		// fmt.Printf("[readUniqKmer] rb: %v\n", rb)
		// }
		cs <- rb
		processNumKmer += 1
	}

	fmt.Printf("[readUniqKmer] total read kmer number is : %d\n", processNumKmer)
	// send read finish signal
	for i := 0; i < numCPU; i++ {
		var rb constructcf.ReadBnt
		rb.Length = 0
		cs <- rb
	}
}

func ExtendNodeKmer(nodeBnt constructcf.ReadBnt, cf cuckoofilter.CuckooFilter, min_kmer_count uint16, direction uint8) (leftcount, rightcount int, baseBnt byte, baseCount uint16) {
	var nBnt constructcf.ReadBnt
	nBnt.Seq = make([]byte, Kmerlen+1)
	copy(nBnt.Seq[1:], nodeBnt.Seq)
	for i := 0; i < bnt.BaseTypeNum; i++ {
		bi := byte(i)
		nBnt.Seq[0] = bi
		ks := constructcf.GetReadBntKmer(nBnt, 0, Kmerlen)
		rs := constructcf.ReverseComplet(ks)
		if ks.BiggerThan(rs) {
			ks, rs = rs, ks
		}
		// fmt.Printf("[ExtendNodeKmer]ks: %v, rs: %v\n", ks, rs)
		count := cf.GetCountAllowZero(ks.Seq)
		if count >= min_kmer_count {
			if direction == BACKWARD {
				baseBnt = uint8(i)
				baseCount = count
			}
			leftcount++
		}
		nBnt.Seq[cf.Kmerlen] = bi
		ks = constructcf.GetReadBntKmer(nBnt, 1, cf.Kmerlen)
		rs = constructcf.ReverseComplet(ks)
		if ks.BiggerThan(rs) {
			ks, rs = rs, ks
		}
		count = cf.GetCountAllowZero(ks.Seq)
		if count >= min_kmer_count {
			if direction == FORWARD {
				baseBnt = uint8(i)
				baseCount = count
			}
			rightcount++
		}
	}

	return
}

func paraLookupComplexNode(cs chan constructcf.ReadBnt, wc chan constructcf.ReadBnt, cf cuckoofilter.CuckooFilter) {
	// var wrsb constructcf.ReadSeqBucket
	for {
		rb := <-cs
		if rb.Length == 0 {
			wc <- rb
			break
		} else {
			// if rsb.Count < constructcf.ReadSeqSize {
			// 	fmt.Printf("rsb.ReadBuf length is : %d\n", len(rsb.ReadBuf))
			// }
			// if found kmer count is 1 , this kmer will be ignore, and skip this branch
			extRBnt := constructcf.ExtendReadBnt2Byte(rb)
			var nodeBnt constructcf.ReadBnt
			nodeBnt.Seq = extRBnt.Seq[:cf.Kmerlen-1]
			nodeBnt.Length = len(nodeBnt.Seq)
			// fmt.Printf("[paraLookupComplexNode] nodeBnt : %v\n", nodeBnt)
			leftcount, rightcount, _, _ := ExtendNodeKmer(nodeBnt, cf, MIN_KMER_COUNT, FORWARD)
			// if leftcount > 0 || rightcount > 1 || (leftcount == 0 && rightcount == 0) {
			if leftcount > 1 || rightcount > 1 || (leftcount == 0 && rightcount == 1) ||
				(leftcount == 1 && rightcount == 0) {
				node := constructcf.GetReadBntKmer(nodeBnt, 0, cf.Kmerlen-1)
				rnode := constructcf.ReverseComplet(node)
				if node.BiggerThan(rnode) {
					node, rnode = rnode, node
				}
				// fmt.Printf("[paraLookupComplexNode] leftcount: %d, rightcount : %d\n", leftcount, rightcount)
				wc <- node
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
	// ckbuffp := bufio.NewReader(ckgzfp)
	NBntByteLen := (Kmerlen - 1 + bnt.BaseTypeNum - 1) / bnt.BaseTypeNum
	eof := false
	readnodeNum := 0
	for !eof {
		b := make([]byte, NBntByteLen)
		err := binary.Read(ckgzfp, binary.LittleEndian, b)
		// n, err := ckfp.Read(b)
		// if n != NBntByteLen {
		// 	log.Fatalf("[constructNodeMap] read node seq err: n(%d) != NBntByteLen(%d\n)", n, NBntByteLen)
		// }
		if err != nil {
			if err == io.EOF {
				err = nil
				eof = true
				break
			} else {
				log.Fatalf("[constructNodeMap] err: %v\n", err)
			}
		}
		readnodeNum++
		var node DBGNode
		node.ID = nodeID
		node.Seq = b
		if _, ok := nodeMap[string(node.Seq)]; ok == false {
			nodeMap[string(node.Seq)] = node
			nodeID++
		}
	}

	fmt.Printf("[constructNodeMap] read node number is : %d\n", readnodeNum)

	return nodeMap, nodeID
}

func ReadDBGNodeToChan(nodeMap map[string]DBGNode, nc chan DBGNode, numCPU int) {

	for _, value := range nodeMap {
		nc <- value
	}

	for i := 0; i < numCPU; i++ {
		var n DBGNode
		n.ID = 0
		nc <- n
	}
}

func ReverseByteArr(ba []byte) {
	lba := len(ba)
	for i := 0; i < lba/2; i++ {
		ba[i], ba[lba-1-i] = ba[lba-1-i], ba[i]
	}
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
func ReverseUint8Arr(ua []uint8) {
	lua := len(ua)
	for i := 0; i < lua/2; i++ {
		ua[i], ua[lua-1-i] = ua[lua-1-i], ua[i]
	}
}

func GetEdges(cf cuckoofilter.CuckooFilter, nBnt constructcf.ReadBnt, count uint8, direction uint8, MIN_KMER_COUNT uint16) (edge DBGEdge, isNode bool) {
	if direction == FORWARD {
		edge.Utg.Ks = append(edge.Utg.Ks, nBnt.Seq...)
		edge.Utg.Kq = make([]uint8, cf.Kmerlen-1)
		edge.Utg.Kq = append(edge.Utg.Kq, count)
		var tbnt constructcf.ReadBnt
		tbnt.Seq = nBnt.Seq[1:]
		tbnt.Length = len(tbnt.Seq)
		for {
			leftcount, rightcount, baseBnt, baseCount := ExtendNodeKmer(tbnt, cf, MIN_KMER_COUNT, FORWARD)
			if leftcount == 1 && rightcount == 1 {
				edge.Utg.Ks = append(edge.Utg.Ks, baseBnt)
				edge.Utg.Kq = append(edge.Utg.Kq, uint8(baseCount))
				tbnt.Seq = tbnt.Seq[1:]
				tbnt.Seq = append(tbnt.Seq, baseBnt)
			} else {
				if leftcount > 1 || rightcount > 1 {
					isNode = true
				}
				break
			}
		}
	} else {
		edge.Utg.Ks = append(edge.Utg.Ks, nBnt.Seq[0])
		edge.Utg.Kq = append(edge.Utg.Kq, count)
		var tbnt constructcf.ReadBnt
		tbnt.Seq = nBnt.Seq[:cf.Kmerlen-1]
		tbnt.Length = len(tbnt.Seq)
		for {
			leftcount, rightcount, baseBnt, baseCount := ExtendNodeKmer(tbnt, cf, MIN_KMER_COUNT, BACKWARD)
			if leftcount == 1 && rightcount == 1 {
				edge.Utg.Ks = append(edge.Utg.Ks, baseBnt)
				edge.Utg.Kq = append(edge.Utg.Kq, uint8(baseCount))
				var seq []byte
				seq = append(seq, baseBnt)
				seq = append(seq, tbnt.Seq[:cf.Kmerlen-2]...)
				tbnt.Seq = seq
			} else {
				if leftcount > 1 || rightcount > 1 {
					isNode = true
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
func paraGenerateDBGEdges(nc chan DBGNode, newNodeChan chan constructcf.ReadBnt, nodeMap map[string]DBGNode, cf cuckoofilter.CuckooFilter, wc chan DBGEdge) {
	for {
		node := <-nc
		if node.ID == 0 {
			var rb constructcf.ReadBnt
			rb.Length = 0
			newNodeChan <- rb
			break
		}

		// read edge seq from cuckoofilter
		var rb constructcf.ReadBnt
		rb.Seq = node.Seq
		rb.Length = cf.Kmerlen - 1
		extRBnt := constructcf.ExtendReadBnt2Byte(rb)
		// leftcount, rightcount, _, _ := ExtendNodeKmer(extRBnt, cf, MIN_KMER_COUNT, FORWARD)
		// fmt.Printf("[paraGenerateDBGEdges] leftcount: %d, rightcount: %d\n", leftcount, rightcount)
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
				// fmt.Printf("[paraGenerateDBGEdges] in: %d, count : %v\n", i, count)
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
						if sks.BiggerThan(trs) {
							sks = trs
						}
						if v, ok := nodeMap[string(sks.Seq)]; ok {
							c := edge.Utg.Ks[cf.Kmerlen-1]
							if sks.Equal(tks) {
								// b := bnt.Base2Bnt[c]
								// if b > bnt.BaseTypeNum-1 {
								// 	fmt.Printf("[paraGenerateDBGEdges]edge: %v, b:%v, c:%v\n", edge, b, c)
								// }
								v.EdgeIDOutcoming[c] = 1
							} else {
								// b := bnt.Base2Bnt[bnt.BitNtRev[c]]
								v.EdgeIDIncoming[bnt.BntRev[c]] = 1
							}
							edge.StartNID = v.ID
							writedEdge = true
							nodeMap[string(sks.Seq)] = v
						} else { // this is new node, write to chan
							newNodeChan <- sks
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
				// fmt.Printf("[paraGenerateDBGEdges] out: %d, count : %v\n", i, count)
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
							if sks.Equal(tks) {
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
							newNodeChan <- sks
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
					}
				}
			}
		}
	}
}

// write edges seq to the file
func WriteEdgesToFn(edgesfn string, wc chan DBGEdge) {
	edgesNum := 0
	edgesfp, err := os.Create(edgesfn)
	if err != nil {
		log.Fatalf("[WriteEdgesToFn] Open file %s failed: %v\n", edgesfn, err)
	}
	defer edgesfp.Close()
	edgesbuffp := bufio.NewWriter(edgesfp)
	defer edgesbuffp.Flush()
	// edgesgzfp := gzip.NewWriter(edgesbuffp)
	// defer edgesgzfp.Close()
	for {
		ei := <-wc
		if len(ei.Utg.Ks) == 0 {
			if ei.StartNID != 0 || ei.EndNID != 0 {
				log.Fatalf("[WriteEdgesToFn] err edge: %v\n", ei)
			}
			fmt.Printf("[WriteEdgesToFn] end edge: %v\n", ei)
			break
		}
		edgesNum++
		fmt.Fprintf(edgesbuffp, ">%d\t%d\t%d\n", edgesNum, ei.StartNID, ei.EndNID)
		// fmt.Fprintf(edgesgzfp, "%s\n", ei.Utg.Ks)
		for i := 0; i < len(ei.Utg.Ks); i++ {
			if ei.Utg.Ks[i] > 3 || ei.Utg.Ks[i] < 0 {
				log.Fatalf("[WriteEdgesToFn] not correct base: %d:%d\n", i, ei.Utg.Ks[i])
			}
			fmt.Fprintf(edgesbuffp, "%c", bnt.BitNtCharUp[ei.Utg.Ks[i]])
		}
		fmt.Fprintf(edgesbuffp, "%s", "\n")
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

	fmt.Printf("[WriteEdgesToFn] the writed file edges number is %d\n", edgesNum)
}

func ProcessAddedNode(cf cuckoofilter.CuckooFilter, nodeMap map[string]DBGNode, newNodeBntArr []constructcf.ReadBnt, wc chan DBGEdge, nodeID DBG_MAX_INT) (addedNodesNum, addedEdgesNum int) {

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
								if sks.Equal(tks) {
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
								if sks.Equal(tks) {
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
}

func cleanEdgeIDInNodeMap(nodeMap map[string]DBGNode) {
	for k, v := range nodeMap {
		for i := 0; i < bnt.BaseTypeNum; i++ {
			v.EdgeIDIncoming[i] = 0
			v.EdgeIDOutcoming[i] = 0
		}
		nodeMap[k] = v
	}
}

func GenerateDBGEdges(nodeMap map[string]DBGNode, cf cuckoofilter.CuckooFilter, edgesfn string, numCPU int, nodeID DBG_MAX_INT) (newNodeID DBG_MAX_INT) {
	bufsize := 50
	nc := make(chan DBGNode, bufsize)
	wc := make(chan DBGEdge, bufsize)
	newNodeChan := make(chan constructcf.ReadBnt, bufsize)
	// Read DBGNode to thw nc
	go ReadDBGNodeToChan(nodeMap, nc, numCPU)
	// parallel construct edges from cuckoofilter
	for i := 0; i < numCPU; i++ {
		go paraGenerateDBGEdges(nc, newNodeChan, nodeMap, cf, wc)
	}
	// write edges Seq to the file
	go WriteEdgesToFn(edgesfn, wc)

	// cache the new added node Bnt info
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
	newNodeID = nodeID + DBG_MAX_INT(addedNodesNum)

	// clean set edgeID in the DBGNode
	cleanEdgeIDInNodeMap(nodeMap)
	fmt.Printf("[GenerateDBGEdges] added nodes number is : %d, added edges number is : %d\n", addedNodesNum, addedEdgesNum)

	return

}

func NodesStatWriter(nodesStatfn string, newNodeID DBG_MAX_INT) {
	nodesStatfp, err := os.Create(nodesStatfn)
	if err != nil {
		log.Fatalf("[NodesStatWriter] file %s create error: %v\n", nodesStatfn, err)
	}
	defer nodesStatfp.Close()
	fmt.Fprintf(nodesStatfp, "nodes size:\t%v\n", newNodeID)
}

func NodesStatReader(nodesStatfn string) (nodesSize DBG_MAX_INT) {
	nodesStatfp, err := os.Open(nodesStatfn)
	if err != nil {
		log.Fatalf("[NodesStatReader] file %s Open error: %v\n", nodesStatfn, err)
	}
	defer nodesStatfp.Close()
	_, err = fmt.Fscanf(nodesStatfp, "nodes size:\t%v\n", &nodesSize)
	if err != nil {
		log.Fatalf("[NodesStatReader] nodes size parse error: %v\n", err)
	}
	return
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
			log.Fatalf("[NodeMap2NodeArr] v.ID:%u >= nodesArr len: %u\n", v.ID, naLen)
		}
		nodesArr[v.ID] = v
	}
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
	wc := make(chan constructcf.ReadBnt, bufsize)
	uniqkmergzfn := prefix + ".uniqkmerseq.gz"
	// read uniq kmers form file
	go readUniqKmer(uniqkmergzfn, cs, cf.Kmerlen, numCPU)

	// identify complex Nodes
	for i := 0; i < numCPU; i++ {
		go paraLookupComplexNode(cs, wc, cf)
	}

	complexKmergzfn := prefix + ".complexNode.gz"
	ckfp, err := os.Create(complexKmergzfn)
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
	complexNodeNum := 0
	for {
		rb := <-wc
		if rb.Length == 0 {
			endFlagCount++
			if endFlagCount == numCPU {
				break
			} else {
				continue
			}
		}

		err := binary.Write(ckgzfp, binary.LittleEndian, rb.Seq)
		if err != nil {
			log.Fatalf("[CDBG] write node seq to file err: %v\n", err)
		}
		// ckgzfp.Write(rb.Seq)
		// ckgzfp.Write([]byte("\n"))
		complexNodeNum += 1
	}
	ckgzfp.Close()
	fmt.Printf("[CDBG] found complex Node num is : %d\n", complexNodeNum)

	// construct Node map
	nodeMap, nodeID := constructNodeMap(complexKmergzfn)
	fmt.Printf("[CDBG] assgin nodeID to  : %d\n", nodeID)
	// parallel generate edges and write to file
	edgefn := prefix + ".edges"
	newNodeID := GenerateDBGEdges(nodeMap, cf, edgefn, numCPU, nodeID)
	nodesStatfn := prefix + ".nodes.stat"
	NodesStatWriter(nodesStatfn, newNodeID)
	// write DBG node map to the file
	nodesfn := prefix + ".nodes.mmap"
	NodeMapMmapWriter(nodeMap, nodesfn)
}

func ParseEdge(edgesbuffp *bufio.Reader) (edge DBGEdge, err error) {

	err = nil
	line1, err1 := edgesbuffp.ReadBytes('\n')
	line2, err2 := edgesbuffp.ReadBytes('\n')
	line3, err3 := edgesbuffp.ReadBytes('\n')
	if err1 != nil || err2 != nil || err3 != nil {
		if err1 == io.EOF {
			err3 = nil
			err = io.EOF
			return
		} else {
			log.Fatal("[ParseEdge] Read edge found err\n")
		}
	}
	var tmp int
	_, err4 := fmt.Sscanf(string(line1), ">%d\t%d\t%d\n", &tmp, &edge.StartNID, &edge.EndNID)
	if err4 != nil {
		log.Fatalf("[ParseEdge] Sscaf line1:%s err:%v\n", line1, err4)
	}
	// _, err4 = fmt.Sscanf(string(line2), "%s\n", &edge.Utg.Ks)
	// if err4 != nil {
	// 	log.Fatalf("[ParseEdge] Sscaf line2 err:%v\n", err4)
	// }
	// change char base to Bnt
	sline2 := string(line2[:len(line2)-1])
	for _, v := range sline2 {
		edge.Utg.Ks = append(edge.Utg.Ks, bnt.Base2Bnt[v])
	}
	sline3 := string(line3[:len(line3)-1])
	for _, v := range sline3 {
		q := v - 33
		edge.Utg.Kq = append(edge.Utg.Kq, uint8(q))
	}

	return
}

func ReadEdgesFromFile(nodeMap map[string]DBGNode, edgesfn string) (edgesArr []DBGEdge) {
	var niledge DBGEdge
	edgesArr = append(edgesArr, niledge)
	edgeID := DBG_MAX_INT(len(edgesArr))
	collisionNum := 0
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

	for {
		edge, err := ParseEdge(edgesbuffp)
		if err == io.EOF {
			break
		} else if err != nil {
			log.Fatalf("[ParseEdge] parse edge err: %v\n", err)
		}
		// if len(edgesArr) > 0 && int(edgesArr[len(edgesArr)-1].ID) != len(edgesArr) {
		// 	log.Fatalf("[ReadEdgesFromFile] edgesArr[len(edgesArr)-1](%v) != len(edgesArr):%d\n", edgesArr[len(edgesArr)-1], len(edgesArr))
		// }
		// fmt.Printf("[ParseEdge] edge: %v\n", edge)
		collisiontag := false
		var collisionEdgeID DBG_MAX_INT
		// check start and end node, 1 note start, 2 note end
		if edge.StartNID > 0 {
			var sBnt constructcf.ReadBnt
			sBnt.Seq = edge.Utg.Ks[:Kmerlen-1]
			ks1 := constructcf.GetReadBntKmer(sBnt, 0, Kmerlen-1)
			ts1 := ks1
			rs1 := constructcf.ReverseComplet(ks1)
			if ks1.BiggerThan(rs1) {
				ks1, rs1 = rs1, ks1
			}
			v1, ok := nodeMap[string(ks1.Seq)]
			// fmt.Printf("[ParseEdge] v1: %v\n", v1)
			if ok && v1.ID == edge.StartNID {
				c := edge.Utg.Ks[Kmerlen-1]
				if ks1.Equal(ts1) {
					// b := bnt.Base2Bnt[c]
					if v1.EdgeIDOutcoming[c] > 0 {
						collisiontag = true
						collisionEdgeID = v1.EdgeIDOutcoming[c]
					} else {
						v1.EdgeIDOutcoming[c] = edgeID
					}
				} else {
					// b := bnt.Base2Bnt[bnt.BitNtRev[c]]
					if v1.EdgeIDIncoming[bnt.BntRev[c]] > 0 {
						collisiontag = true
						collisionEdgeID = v1.EdgeIDIncoming[bnt.BntRev[c]]
					} else {
						v1.EdgeIDIncoming[bnt.BntRev[c]] = edgeID
					}
				}
				if collisiontag == false {
					nodeMap[string(ks1.Seq)] = v1
				}
			} else {
				log.Fatalf("[ReadEdgesFromFile] edge.StartNID not consistence with nodeMap\nedge: %v\nok: %v\n", edge, ok)
			}

		}

		if collisiontag == false {
			if edge.EndNID > 0 {
				var eBnt constructcf.ReadBnt
				eBnt.Seq = edge.Utg.Ks[len(edge.Utg.Ks)-Kmerlen+1:]
				ks2 := constructcf.GetReadBntKmer(eBnt, 0, Kmerlen-1)
				ts2 := ks2
				rs2 := constructcf.ReverseComplet(ks2)
				if ks2.BiggerThan(rs2) {
					ks2, rs2 = rs2, ks2
				}
				v2, ok := nodeMap[string(ks2.Seq)]
				// fmt.Printf("[ParseEdge] v2: %v\n", v2)
				if ok && v2.ID == edge.EndNID {
					c := edge.Utg.Ks[len(edge.Utg.Ks)-Kmerlen]
					// var collisionEdgeID DBG_MAX_INT
					if ks2.Equal(ts2) {
						// b := bnt.Base2Bnt[c]
						if v2.EdgeIDIncoming[c] > 0 {
							collisiontag = true
							collisionEdgeID = v2.EdgeIDIncoming[c]
						} else {
							v2.EdgeIDIncoming[c] = edgeID
						}
					} else {
						// b := bnt.Base2Bnt[bnt.BitNtRev[c]]
						if v2.EdgeIDOutcoming[bnt.BntRev[c]] > 0 {
							collisiontag = true
							collisionEdgeID = v2.EdgeIDOutcoming[bnt.BntRev[c]]
						} else {
							v2.EdgeIDOutcoming[bnt.BntRev[c]] = edgeID
						}
					}
					if collisiontag == false {
						nodeMap[string(ks2.Seq)] = v2
					} else {
						edge.EndNID = 0
						collisionNum++
						fmt.Printf("v2: %v\nlen(collsionedge.Utg.Ks) = %d, collsionedge = %v\nlen(edge.Utg.Ks) = %d, edge = %v\n\n", v2, len(edgesArr[collisionEdgeID].Utg.Ks), edgesArr[collisionEdgeID], len(edge.Utg.Ks), edge)
					}
				} else {
					log.Fatalf("[ReadEdgesFromFile] edge.EndNID not consistence with nodeMap\n")
				}

			}
			edge.ID = edgeID
			edgeID++
			edgesArr = append(edgesArr, edge)
			if edge.Utg.Ks[0] > 3 {
				log.Fatalf("[ReadEdgesFromFile] edge seq need transform bnt: %v\n", edge)
			}

		} else {
			if (edgesArr[collisionEdgeID].StartNID == edge.StartNID || edgesArr[collisionEdgeID].StartNID == edge.EndNID) && (edgesArr[collisionEdgeID].EndNID == edge.StartNID || edgesArr[collisionEdgeID].EndNID == edge.EndNID) {

			} else {
				fmt.Printf("[ReadEdgesFromFile] found collision edges, len(collsionedge.Utg.Ks) = %d, collsionedge = %v\nlen(edge.Utg.Ks) = %d, edge = %v\n\n", len(edgesArr[collisionEdgeID].Utg.Ks), edgesArr[collisionEdgeID], len(edge.Utg.Ks), edge)
				collisionNum++
			}
		}
	}
	fmt.Printf("[ReadEdgesFromFile] found collision edge number is : %d, percent: %f\n", collisionNum, float64(collisionNum)/float64((len(edgesArr)-1)))
	return
}

func RCEdge(edgesArr []DBGEdge, eid DBG_MAX_INT) {
	edgesArr[eid].StartNID, edgesArr[eid].EndNID = edgesArr[eid].EndNID, edgesArr[eid].StartNID
	ReverseCompByteArr(edgesArr[eid].Utg.Ks)
	ReverseUint8Arr(edgesArr[eid].Utg.Kq)
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

func ReconstructConsistenceDBG(nodeMap map[string]DBGNode, edgesArr []DBGEdge) {
	for k, v := range nodeMap {
		stk := list.New()
		if v.GetProcessFlag() == 0 {
			v.SetProcessFlag(uint8(1))
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
										if v2Bnt.Equal(ks) == false {
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
										if v2Bnt.Equal(ks) == false {
											v2 = RevNode(v2)
										}
										v2.SetProcessFlag(uint8(1))
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
										if v2Bnt.Equal(ks) == false {
											// log.Fatalf("[ReconstructConsistenceDBG] found not consistence node\n")
											fmt.Printf("[ReconstructConsistenceDBG] found not consistence node, edge: %v\nv1: %v\nv2: %v\n", edgesArr[eid], node, v2)
										}
									} else {
										if v2Bnt.Equal(ks) == false {
											v2 = RevNode(v2)
										}
										v2.SetProcessFlag(uint8(1))
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
}

/*func Coming2String(coming [bnt.BaseTypeNum]DBG_MAX_INT) (cs string) {
	for _, v := range coming {
		cs += " " + strconv.Itoa(int(v))
	}
	return
}*/

func GraphvizDBG(nodeMap map[string]DBGNode, edgesArr []DBGEdge, graphfn string) {
	// create a new graph
	g := gographviz.NewGraph()
	g.SetName("G")
	g.SetDir(true)
	g.SetStrict(false)
	for _, v := range nodeMap {
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

func GetEdgeIDComing(coming [bnt.BaseTypeNum]DBG_MAX_INT) (num int, edgeID DBG_MAX_INT) {
	for _, v := range coming {
		if v > 0 {
			num++
			edgeID = v
		}
	}

	return
}

func ConcatEdges(edgesArr []DBGEdge, inID, outID, dstID DBG_MAX_INT) {
	// check if is connective
	var inBnt, outBnt constructcf.ReadBnt
	inBnt.Seq = edgesArr[inID].Utg.Ks[len(edgesArr[inID].Utg.Ks)-Kmerlen+1:]
	inBnt.Length = len(inBnt.Seq)
	outBnt.Seq = edgesArr[outID].Utg.Ks[:Kmerlen-1]
	outBnt.Length = len(outBnt.Seq)
	fmt.Printf("[ConcatEdges] inID: %d, outID: %d, dstID: %d\nseq1:%v\nseq2:%v\n", inID, outID, dstID, inBnt, outBnt)
	if inBnt.Equal(outBnt) == false {
		log.Fatalf("[ConcatEdges] two edges is not connective\n")
	}
	if dstID == inID {
		edgesArr[inID].Utg.Ks = append(edgesArr[inID].Utg.Ks, edgesArr[outID].Utg.Ks[Kmerlen-1:]...)
		for i := 0; i < Kmerlen-1; i++ {
			if edgesArr[inID].Utg.Kq[len(edgesArr[inID].Utg.Kq)-Kmerlen+1+i] < edgesArr[outID].Utg.Kq[i] {
				edgesArr[inID].Utg.Kq[len(edgesArr[inID].Utg.Kq)-Kmerlen+1+i] = edgesArr[outID].Utg.Kq[i]
			}
		}
		edgesArr[inID].Utg.Kq = append(edgesArr[inID].Utg.Kq, edgesArr[outID].Utg.Kq[Kmerlen-1:]...)
		// DeleteEdgeID(nodesArr, edgesArr[inID].EndNID, inID)
		edgesArr[inID].EndNID = edgesArr[outID].EndNID

	} else {
		seq := make([]byte, len(edgesArr[inID].Utg.Ks))
		copy(seq, edgesArr[inID].Utg.Ks)
		edgesArr[outID].Utg.Ks = append(seq, edgesArr[outID].Utg.Ks[Kmerlen-1:]...)
		qul := make([]uint8, len(edgesArr[inID].Utg.Kq))
		copy(qul, edgesArr[inID].Utg.Kq)
		for i := 0; i < Kmerlen-1; i++ {
			if edgesArr[inID].Utg.Kq[len(edgesArr[inID].Utg.Kq)-Kmerlen+1+i] < edgesArr[outID].Utg.Kq[i] {
				qul[len(qul)-Kmerlen+1+i] = edgesArr[outID].Utg.Kq[i]
			}
		}
		edgesArr[outID].Utg.Kq = append(qul, edgesArr[outID].Utg.Kq[Kmerlen-1:]...)
		// DeleteEdgeID(nodesArr, edgesArr[outID].StartNID, outID)
		edgesArr[outID].StartNID = edgesArr[inID].StartNID

	}
}

func SubstituteEdgeID(nodeMap map[string]DBGNode, nodekey []byte, srcID, dstID DBG_MAX_INT) bool {
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
		log.Fatalf("[SubstituteEdgeID] not found correct node\n")
	}

	return suc

}

func SmfyDBG(nodeMap map[string]DBGNode, edgesArr []DBGEdge) {

	deleteNodeNum, deleteEdgeNum := 0, 0
	// longTipsEdgesNum := 0
	for k, v := range nodeMap {
		inNum, inID := GetEdgeIDComing(v.EdgeIDIncoming)
		outNum, outID := GetEdgeIDComing(v.EdgeIDOutcoming)
		if inNum == 0 && outNum == 0 {
			delete(nodeMap, k)
			deleteNodeNum++
		} else if inNum == 1 && outNum == 1 && inID != outID { // prevent cycle ring
			e1 := edgesArr[inID]
			e2 := edgesArr[outID]
			if e1.StartNID == v.ID {
				RCEdge(edgesArr, inID)
			}
			if e2.EndNID == v.ID {
				RCEdge(edgesArr, outID)
			}
			ConcatEdges(edgesArr, inID, outID, inID)
			// edgesArr[inID].EndNID = edgesArr[outID].EndNID
			edgesArr[outID].SetDeleteFlag()
			deleteEdgeNum++
			if edgesArr[outID].EndNID > 0 {
				if SubstituteEdgeID(nodeMap, edgesArr[inID].Utg.Ks[len(edgesArr[inID].Utg.Ks)-Kmerlen+1:], outID, inID) == false {
					log.Fatalf("[SmfyDBG] SubstituteEdgeID failed\n")
				}
			}
			delete(nodeMap, k)
			deleteNodeNum++
		}
	}

	fmt.Printf("[SmfyDBG]deleted nodes number is : %d\n", deleteNodeNum)
	fmt.Printf("[SmfyDBG]deleted edges number is : %d\n", deleteEdgeNum)
}

func transform2QSeq(utg Unitig) (qs alphabet.QLetters) {
	if len(utg.Ks) != len(utg.Kq) {
		log.Fatalf("[transform2QSeq] len(ks):%d != len(kq):%d\n", len(utg.Ks), len(utg.Kq))
	}
	for i := 0; i < len(utg.Ks); i++ {
		var ql alphabet.QLetter
		ql.L = alphabet.Letter(bnt.BitNtCharUp[utg.Ks[i]])
		ql.Q = alphabet.Qphred(utg.Kq[i])
		qs = append(qs, ql)
	}

	return
}

func transform2Unitig(Seq alphabet.QLetters, utg *Unitig, lenKs int) {

	utg.Ks = make([]byte, lenKs)
	utg.Kq = make([]uint8, lenKs)
	for i, v := range Seq {
		utg.Ks[i] = bnt.Base2Bnt[v.L]
		utg.Kq[i] = uint8(v.Q)
	}
}

var AdpaterSeq = "cggccgcaaggggttcgcgtcagcgggtgttggcgggtgtcggggctggcttaactatgcggcatcagagcagattgtactgagagtgcaccatatgcggtgtgaaataccacacagatgcgtaaggagaaaataccgcatcaggcgccattcgccattcagctgcgcaactgttgggaagggcgatcggtgcgggcctc"

func SetDefaultQual(Adapter Unitig) (new Unitig) {
	for i := 0; i < len(Adapter.Ks); i++ {
		Adapter.Ks[i] = bnt.Base2Bnt[Adapter.Ks[i]]
		Adapter.Kq = append(Adapter.Kq, uint8(1))
	}

	return Adapter
}
func StoreEdgesToFn(edgesfn string, edgesArr []DBGEdge, addAdapter bool) {
	fp, err := os.Create(edgesfn)
	if err != nil {
		log.Fatalf("[StoreEdgesToFn] create file: %s failed, err: %v\n", edgesfn, err)
	}
	defer fp.Close()
	var Adapter Unitig
	if addAdapter == true {
		UpAdpaterSeq := strings.ToUpper(AdpaterSeq)
		Adapter.Ks = []byte(UpAdpaterSeq)
		Adapter = SetDefaultQual(Adapter)
		fmt.Printf("[StoreEdgesToFn] Unitig of Adapter:%v\n", Adapter)
	}

	fqfp := fastq.NewWriter(fp)
	for _, v := range edgesArr {
		if v.ID > 0 && v.GetDeleteFlag() == 0 {
			seq := linear.NewQSeq("", nil, alphabet.DNA, alphabet.Sanger)
			seq.ID = strconv.Itoa(int(v.ID))
			// Add start and end adapter seq
			qs := transform2QSeq(v.Utg)
			if addAdapter == true {
				as := transform2QSeq(Adapter)
				seq.AppendQLetters(as...)
				seq.AppendQLetters(qs...)
				seq.AppendQLetters(as...)
			} else {
				seq.AppendQLetters(qs...)
			}
			ans := strconv.Itoa(int(v.StartNID)) + "\t" + strconv.Itoa(int(v.EndNID)) + "\tlen:" + strconv.Itoa(seq.Len())
			seq.Annotation.SetDescription(ans)
			_, err := fqfp.Write(seq)
			if err != nil {
				log.Fatalf("[StoreEdgesToFn] write seq: %v; err: %v\n", seq, err)
			}
		}
	}
}

func LoadEdgesfqFromFn(fn string, edgesArr []DBGEdge) {
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
			id, err := strconv.Atoi(l.Name())
			if err != nil {
				log.Fatalf("[LoadEdgesfqFromFn] parse Name:%s of fastq err: %v\n", l.Name(), err)
			}
			edge.ID = DBG_MAX_INT(id)
			var lenKs int
			_, err = fmt.Sscanf(l.Description(), "%v\t%v\tlen:%d\n", &edge.StartNID, &edge.EndNID, &lenKs)
			if err != nil {
				log.Fatalf("[LoadEdgesfqFromFn] parse Description:%s of fastq err: %v\n", l.Description(), err)
			}
			transform2Unitig(l.Seq, &edge.Utg, lenKs)
			if edge.ID >= DBG_MAX_INT(len(edgesArr)) {
				log.Fatalf("[LoadEdgesfqFromFn] edge.ID:%v >= len(edgesArr):%d\n", edge.ID, len(edgesArr))
			} else if edgesArr[edge.ID].ID > 0 {
				log.Fatalf("[LoadEdgesfqFromFn] the position: %u in edgesArr has value:%v\n", edge.ID, edgesArr[edge.ID])
			}
			edgesArr[edge.ID] = edge
		}

	}
}

func Smfy(c cli.Command) {
	k := c.Parent().Flag("K").String()
	var err error = nil
	Kmerlen, err = strconv.Atoi(k)
	if err != nil {
		log.Fatalf("[Smfy] argument 'K' set error\n")
	}
	prefix := c.Parent().Flag("p").String()
	nodesfn := prefix + ".nodes.mmap"
	nodeMap := NodeMapMmapReader(nodesfn)
	fmt.Printf("[Smfy] len(nodeMap): %v\n", len(nodeMap))
	// read edges file
	edgesfn := prefix + ".edges"
	edgesArr := ReadEdgesFromFile(nodeMap, edgesfn)
	// Debug code
	// for i := 1; i < len(edgesArr); i++ {
	// 	fmt.Printf("[Smfy]edgesArr[%d]: %v\n", i, edgesArr[i])
	// 	if len(edgesArr[i].Utg.Ks) == len(edgesArr[i-1].Utg.Ks) {
	// 		fmt.Printf("[Smfy] outcoming: %v, %v, len: %d\n", edgesArr[i-1].Utg.Ks[Kmerlen-1], edgesArr[i].Utg.Ks[Kmerlen-1], len(edgesArr[i].Utg.Ks))
	// 	}
	// }
	// Debug code
	// output graphviz graph
	graphfn := prefix + ".dot"
	// reconstruct consistence De Bruijn Graph
	// ReconstructConsistenceDBG(nodeMap, edgesArr)
	SmfyDBG(nodeMap, edgesArr)
	GraphvizDBG(nodeMap, edgesArr, graphfn)
	smfyEdgesfn := prefix + ".edges.smfy.fq"
	StoreEdgesToFn(smfyEdgesfn, edgesArr, false)
	adpaterEdgesfn := prefix + ".edges.adapter.fq"
	StoreEdgesToFn(adpaterEdgesfn, edgesArr, true)
	edgesStatfn := prefix + ".edges.stat"
	EdgesStatWriter(edgesStatfn, len(edgesArr))
	smfyNodesfn := prefix + ".nodes.smfy.mmap"
	NodeMapMmapWriter(nodeMap, smfyNodesfn)
}
