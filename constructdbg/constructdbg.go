package constructdbg

import (
	"bufio"
	"compress/gzip"
	"fmt"
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
	for !eof {
		//var kmer []byte
		line, err := ukbuffp.ReadBytes('\n')
		if err != nil {
			if err == io.EOF {
				err = nil
				eof = true
			} else {
				log.Fatal(err)
			}
		}
		var rb constructcf.ReadBnt
		rb.Seq = line[:len(line)-1]
		if len(rb.Seq) != kmerlen {
			log.Fatalf("[readUniqKmer] len(rb.Seq) != kmerlen\n")
		} else {
			rb.Length = kmerlen
		}
		cs <- rb
	}
	// send read finish signal
	for i := 0; i < numCPU; i++ {
		var rb constructcf.ReadBnt
		rb.Length = 0
		cs <- rb
	}
}

func paraLookupComplexNode(cs chan constructcf.ReadBnt, wc chan constructcf.ReadBnt, cf cuckoofilter.CuckooFilter) {
	var wrsb constructcf.ReadSeqBucket
	for {
		rsb := <-cs

		if rsb.End == true {
			wc <- wrsb
			wc <- rsb
			break
		} else {
			if rsb.Count < constructcf.ReadSeqSize {
				fmt.Printf("rsb.ReadBuf length is : %d\n", len(rsb.ReadBuf))
			}
			leftcount := 0
			rightcount := 0
			for i := 0; i < rsb.Count; i++ {
				extRBnt := constructcf.ExtendReadBnt2Byte(rsb.ReadBuf[i])
				var nBnt constructcf.ReadBnt
				nBnt.Seq = make([]byte, cf.Kmerlen+1)
				copy(nBnt.Seq, extRBnt.Seq)
				for j := 0; j < 4; j++ {
					bj := byte(j)
					if bj != extRBnt.Seq[0] {
						nBnt.Seq[0] = bj
						ks := constructcf.GetReadBntKmer(nBnt, 0, cf.Kmerlen)
						rs := constructcf.ReverseComplet(ks)
						if ks.BiggerThan(rs) {
							ks, rs = rs, ks
						}
						if cf.Lookup(ks.Seq) {
							leftcount++
						}
					}
					nBnt.Seq[cf.Kmerlen] = bj
					ks := constructcf.GetReadBntKmer(nBnt, 1, cf.Kmerlen)
					rs := constructcf.ReverseComplet(ks)
					if ks.BiggerThan(rs) {
						ks, rs = rs, ks
					}
					if cf.Lookup(ks.Seq) {
						rightcount++
					}
				}
				if leftcount > 0 || rightcount > 1 || (leftcount == 0 && rightcount == 0) {
					node := constructcf.GetReadBntKmer(extRBnt, 1, cf.Kmerlen-1)
					if wrsb.Count >= constructcf.ReadSeqSize {
						wc <- wrsb
						wrsb.Count = 0
					}
					wrsb.ReadBuf[wrsb.Count] = node
				}
			}
		}
	}
}

func CDBG(c cli.Command) {
	fmt.Println(c.Flags(), c.Parent().Flags())

	// get set arguments
	// t0 := time.Now()
	numCPU, err := strconv.Atoi(c.Parent().Flag("t").String())
	if err != nil {
		log.Fatal(err)
	}
	runtime.GOMAXPROCS(numCPU)
	prefix := c.Parent().Flag("p").String()
	cfinfofn := prefix + ".cfInfo"
	cf, err := cuckoofilter.RecoverCuckooFilterInfo(cfinfofn)
	if err != nil {
		log.Fatal(err)
	}
	cfmmapfn := prefix + ".cfmmap"
	err = cf.MmapReader(cfmmapfn)
	if err != nil {
		log.Fatal(err)
	}
	bufsize := 100
	cs := make(chan constructcf.ReadBnt, bufsize)
	defer cs.Close()
	wc := make(chan constructcf.ReadBnt, bufsize)
	defer wc.Close()
	uniqkmergzfn := prefix + ".uniqkmerseq.gz"
	// read uniq kmers form file
	go readUniqKmer(uniqkmergzfn, cs, cf.Kmerlen, numCPU-1)

	// identify complex Nodes
	for i := 0; i < numCPU-1; i++ {
		go paraLookupComplexNode(cs, wc, cf)
	}

	complexKmerfn := prefix + ".complexNode.gz"
	ckfp, err := os.Create(complexKmerfn)
	if err != nil {
		log.Fatal(err)
	}
	defer ckfp.Close()
	ckgzfp := gzip.NewWriter(ckfp)
	if err != nil {
		log.Fatal(err)
	}
	defer ckgzfp.Close()
	endFlagCount := 0

	// write complex node to the file
	for {
		rb := <-wc
		if rb.Length == 0 {
			endFlagCount++
			if endFlagCount == numCPU-1 {
				break
			} else {
				continue
			}
		}

		ckgzfp.Write(rb.Seq)
	}

}
