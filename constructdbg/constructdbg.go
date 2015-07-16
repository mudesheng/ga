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
	// "time"

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

func (n DBGNode) GetProcessFlag() uint8 {
	return n.Flag & 0x1
}

func (n DBGNode) SetProcessFlag(f uint8) {
	n.Flag = n.Flag & 0xFE
	if f == 1 {
		n.Flag = n.Flag | 0x1
	} else if f != 0 {
		log.Fatalf("[SetProcessFlag] input flag error: %d\n", f)
	}
}

type Unitig struct {
	Ks []byte
	Kq []uint8
}

type DBGEdge struct {
	ID       DBG_MAX_INT
	StartNID DBG_MAX_INT // start node ID
	EndNID   DBG_MAX_INT // end node ID
	Flag     uint8
	Utg      Unitig
}

func (e DBGEdge) GetProcessFlag() uint8 {
	return e.Flag & 0x1
}

func (e DBGEdge) SetProcessFlag(f uint8) {
	e.Flag = e.Flag & 0xFE
	if f == 1 {
		e.Flag = e.Flag | 0x1
	} else if f != 0 {
		log.Fatalf("[SetProcessFlag] input flag error: %d\n", f)
	}
}

var MIN_KMER_COUNT uint16 = 2
var BACKWARD uint8 = 1
var FORWARD uint8 = 2
var Kmerlen int

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
				break
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
		seq[i], seq[ls-1-i] = bnt.BitNtRev[seq[ls-1-i]], bnt.BitNtRev[seq[i]]
	}
	if ls%2 != 0 {
		seq[ls/2] = bnt.BitNtRev[seq[ls/2]]
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
							newNodeChan <- tks
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
							newNodeChan <- tks
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
	edgesgzfp := gzip.NewWriter(edgesfp)
	defer edgesgzfp.Close()
	for {
		ei := <-wc
		if len(ei.Utg.Ks) == 0 {
			break
		}
		edgesNum++
		fmt.Fprintf(edgesgzfp, ">%d\t%d\t%d\n", edgesNum, ei.StartNID, ei.EndNID)
		// fmt.Fprintf(edgesgzfp, "%s\n", ei.Utg.Ks)
		for i := 0; i < len(ei.Utg.Ks); i++ {
			if ei.Utg.Ks[i] > 3 || ei.Utg.Ks[i] < 0 {
				log.Fatalf("[WriteEdgesToFn] not correct base: %d:%d\n", i, ei.Utg.Ks[i])
			}
			fmt.Fprintf(edgesgzfp, "%c", bnt.BitNtCharUp[ei.Utg.Ks[i]])
		}
		fmt.Fprintf(edgesgzfp, "%s", "\n")
		// write quality to the file
		for i := 0; i < len(ei.Utg.Kq); i++ {
			fmt.Fprintf(edgesgzfp, "%c", ei.Utg.Kq[i]+33)
		}
		fmt.Fprintf(edgesgzfp, "\n")
		// fmt.Printf("[WriteEdgesToFn] edge: %v\n", ei)
	}

	fmt.Printf("[WriteEdgesToFn] the writed file edges number is %d\n", edgesNum)
}

func ProcessAddedNode(cf cuckoofilter.CuckooFilter, nodeMap map[string]DBGNode, newNodeBntArr []constructcf.ReadBnt, wc chan DBGEdge, nodeID DBG_MAX_INT) (addedNodesNum, addedEdgesNum int) {

	InitialLen := len(newNodeBntArr)
	processAdded := 0
	for _, item := range newNodeBntArr {
		var node DBGNode
		node.ID = nodeID
		node.Seq = item.Seq
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
								newNodeBntArr = append(newNodeBntArr, tks)
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
								newNodeBntArr = append(newNodeBntArr, tks)
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
						}
					}
				}
			}
		}
	}

	// add a nil edge to wc, tell have not any more edge need to write
	var edge DBGEdge
	wc <- edge

	fmt.Printf("[ProcessAddedNode] initial newNodeBntArr len: %d, added node number: %d, at last newNodeBntArr len:%d\n", InitialLen, processAdded, len(newNodeBntArr))

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

func GenerateDBGEdges(nodeMap map[string]DBGNode, cf cuckoofilter.CuckooFilter, edgesfn string, numCPU int, nodeID DBG_MAX_INT) {
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

	// clean set edgeID in the DBGNode
	cleanEdgeIDInNodeMap(nodeMap)
	fmt.Printf("[GenerateDBGEdges] added nodes number is : %d, added edges number is : %d\n", addedNodesNum, addedEdgesNum)

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
	err = dec.Decode(nodeMap)
	if err != nil {
		log.Fatalf("[NodeMapMmapReader] decode failed, err: %v\n", err)
	}

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
	cfmmapfn := prefix + ".cfmmap"
	cf, err := cuckoofilter.MmapReader(cfmmapfn)
	if err != nil {
		log.Fatalf("[CDBG] CuckooFilter mmap reader err: %v\n", err)
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
	edgegzfn := prefix + ".edges.gz"
	GenerateDBGEdges(nodeMap, cf, edgegzfn, numCPU, nodeID)
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
		if err3 == io.EOF {
			err3 = nil
			err = io.EOF
		} else {
			log.Fatal("[ReadEdgesFromFile] Read edge found err\n")
		}
	}
	var tmp int
	_, err4 := fmt.Sscanf(string(line1), ">%d\t%d\t%d\n", tmp, edge.StartNID, edge.EndNID)
	if err4 != nil {
		log.Fatalf("[ReadEdgesFromFile] Sscaf line1 err:%v\n", err4)
	}
	_, err4 = fmt.Sscanf(string(line2), "%s\n", edge.Utg.Ks)
	if err4 != nil {
		log.Fatalf("[ReadEdgesFromFile] Sscaf line2 err:%v\n", err4)
	}
	sline3 := string(line3[:len(line3)-1])
	for _, v := range sline3 {
		q := v - '0'
		edge.Utg.Kq = append(edge.Utg.Kq, uint8(q))
	}

	return
}

func ReadEdgesFromFile(nodeMap map[string]DBGNode, edgesfn string) (edgesArr []DBGEdge) {
	edgeID := DBG_MAX_INT(1)
	var niledge DBGEdge
	edgesArr = append(edgesArr, niledge)
	collisionNum := 0
	edgesfp, err := os.Open(edgesfn)
	if err != nil {
		log.Fatalf("[ReadEdgesFromFile] open file %s failed, err: %v\n", edgesfn, err)
	}
	defer edgesfp.Close()
	edgesgzfp, err := gzip.NewReader(edgesfp)
	if err != nil {
		log.Fatalf("[ReadEdgesFromFile] read gz file: %s failed, err: %v\n", edgesfn, err)
	}
	defer edgesgzfp.Close()
	edgesbuffp := bufio.NewReader(edgesgzfp)

	for {
		edge, err := ParseEdge(edgesbuffp)
		if err == io.EOF {
			break
		}
		collisiontag := false
		// check start and end node, 1 note start, 2 note end
		var sBnt constructcf.ReadBnt
		sBnt.Seq = edge.Utg.Ks[:Kmerlen-1]
		ks1 := constructcf.GetReadBntKmer(sBnt, 0, Kmerlen-1)
		ts1 := ks1
		rs1 := constructcf.ReverseComplet(ks1)
		if ks1.BiggerThan(rs1) {
			ks1, rs1 = rs1, ks1
		}
		v1, ok := nodeMap[string(ks1.Seq)]
		if ok && v1.ID == edge.StartNID {
			c := edge.Utg.Ks[Kmerlen-1]
			if ks1.Equal(ts1) {
				b := bnt.Base2Bnt[c]
				if v1.EdgeIDOutcoming[b] > 0 {
					collisiontag = true
				} else {
					v1.EdgeIDOutcoming[b] = edgeID
				}
			} else {
				b := bnt.Base2Bnt[bnt.BitNtRev[c]]
				if v1.EdgeIDIncoming[b] > 0 {
					collisiontag = true
				} else {
					v1.EdgeIDIncoming[b] = edgeID
				}
			}
		} else {
			log.Fatalf("[ReadEdgesFromFile] edge.StartNID not consistence with nodeMap\n")
		}

		if collisiontag == false {
			var eBnt constructcf.ReadBnt
			eBnt.Seq = edge.Utg.Ks[len(edge.Utg.Ks)-Kmerlen+1:]
			ks2 := constructcf.GetReadBntKmer(eBnt, 0, Kmerlen-1)
			ts2 := ks2
			rs2 := constructcf.ReverseComplet(ks2)
			if ks2.BiggerThan(rs2) {
				ks2, rs2 = rs2, ks2
			}
			if v2, ok := nodeMap[string(ks2.Seq)]; ok && v2.ID == edge.EndNID {
				c := edge.Utg.Ks[len(edge.Utg.Ks)-Kmerlen]
				if ks2.Equal(ts2) {
					b := bnt.Base2Bnt[c]
					if v2.EdgeIDIncoming[b] > 0 {
						collisiontag = true
					} else {
						v2.EdgeIDIncoming[b] = edgeID
					}
				} else {
					b := bnt.Base2Bnt[bnt.BitNtRev[c]]
					if v2.EdgeIDOutcoming[b] > 0 {
						collisiontag = true
					} else {
						v2.EdgeIDOutcoming[b] = edgeID
					}
				}
				if collisiontag == false {
					nodeMap[string(ks1.Seq)] = v1
					nodeMap[string(ks2.Seq)] = v2
					edge.ID = edgeID
					edgesArr = append(edgesArr, edge)
					edgeID++
				}
			} else {
				log.Fatalf("[ReadEdgesFromFile] edge.EndNID not consistence with nodeMap\n")
			}
		}

		if collisiontag == true {
			collisionNum++
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

func RevNode(nodeMap map[string]DBGNode, rBnt constructcf.ReadBnt) {
	v := nodeMap[string(rBnt.Seq)]
	rs := constructcf.ReverseComplet(rBnt)
	v.Seq = rs.Seq
	for i := 0; i < bnt.BaseTypeNum; i++ {
		v.EdgeIDIncoming[i], v.EdgeIDOutcoming[bnt.BaseTypeNum-1-i] = v.EdgeIDOutcoming[bnt.BaseTypeNum-1-i], v.EdgeIDIncoming[i]
	}
}

func ReconstructConsistenceDBG(nodeMap map[string]DBGNode, edgesArr []DBGEdge) {

	for k, v := range nodeMap {
		stk := list.New()
		if v.GetProcessFlag() == 0 {
			v.SetProcessFlag(uint8(1))
			nodeMap[k] = v
			stk.PushBack(v)
			for stk.Len() > 0 {
				// Pop a element from stack
				e := stk.Back()
				stk.Remove(e)
				node := e.Value.(DBGNode)
				// Processed flag
				if node.GetProcessFlag() != 1 {
					log.Fatalf("[ReconstructConsistenceDBG] node have not been set processed flag\n")
				}
				for i := 0; i < bnt.BaseTypeNum; i++ {
					if node.EdgeIDOutcoming[i] > 0 {
						eid := node.EdgeIDOutcoming[i]
						if edgesArr[eid].GetProcessFlag() == 0 {
							if edgesArr[eid].StartNID != node.ID {
								RCEdge(edgesArr, eid)
							}
							if edgesArr[eid].EndNID > 0 {
								var tBnt constructcf.ReadBnt
								tBnt.Seq = edgesArr[eid].Utg.Ks[len(edgesArr[eid].Utg.Ks)-Kmerlen+1:]

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
											log.Fatalf("[ReconstructConsistenceDBG] found not consistence node\n")
										}
									} else {
										if v2Bnt.Equal(ks) == false {
											RevNode(nodeMap, min)
										}
										nodeMap[string(min.Seq)].SetProcessFlag(uint8(1))
										stk.PushBack(nodeMap[string(min.Seq)])

									}

								} else {
									log.Fatalf("[ReconstructConsistenceDBG] not found node\n")
								}
							}
							edgesArr[eid].SetProcessFlag(uint8(1))
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
											log.Fatalf("[ReconstructConsistenceDBG] found not consistence node\n")
										}
									} else {
										if v2Bnt.Equal(ks) == false {
											RevNode(nodeMap, min)
										}
										nodeMap[string(min.Seq)].SetProcessFlag(uint8(1))
										stk.PushBack(nodeMap[string(min.Seq)])
									}
								} else {
									log.Fatalf("[ReconstructConsistenceDBG] not found node\n")
								}
							}
							edgesArr[eid].SetProcessFlag(uint8(1))
						}
					}
				}
			}
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
	nodesfn := prefix + "nodes.mmap"
	nodeMap := NodeMapMmapReader(nodesfn)
	// read edges file
	edgesfn := prefix + "edges.gz"
	edgesArr := ReadEdgesFromFile(nodeMap, edgesfn)

	// reconstruct consistence De Bruijn Graph
	ReconstructConsistenceDBG(nodeMap, edgesArr)
	// SmfyDBG(nodemap, edgesArr)
}
