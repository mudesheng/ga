package constructcf

import (
	"GA/bnt"
	"GA/cuckoofilter"
	"bufio"
	"compress/gzip"
	"fmt"
	"github.com/jwaldrip/odin/cli"
	"io"
	"log"
	"os"
	"runtime"
	"strconv"
	"strings"
	"time"
)

const (
	AllState    = 1
	ScaffState  = 2
	GapState    = 3
	EndString   = "end"
	ReadSeqSize = 100
)

type LibInfo struct {
	name          string // name of library
	numberReads   int64  // the number of reads
	diverse       uint8
	paired        uint8
	asmFlag       uint8 // denote which assembly phase used, note 1 used for all step of assembly pipeline, note 2 used for scaffold phase only, 3 used for filling gap only
	seqProfile    uint8 // denote the data origin
	qualBenchmark uint8 // the benchmark of quality score, strength encourage use phred+33
	totalBasesNum int64
	readLen       int
	insertSize    int // paired read insert size
	insertSD      int // Standard Deviation
	//	fnNum         int      // the number of files
	fnName []string // the files name slice
}

type CfgInfo struct {
	maxRdLen int // maximum read length
	minRdLen int // minimum read length
	//	n        int // the number of library
	libs []LibInfo
}

type ReadBnt struct {
	seq    []byte
	length int // length of sequence
}

type ReadSeqBucket struct {
	end     bool
	readBuf [ReadSeqSize]ReadBnt
	count   int
}

/*type ChanStruct struct {
	readseq chan []string
	state   chan int // note if goroutinue can return
} */

func (s1 ReadBnt) BiggerThan(s2 ReadBnt) bool {
	if s1.length < s2.length {
		return false
	} else if s1.length > s2.length {
		return true
	} else {
		for i := 0; i < len(s1.seq); i++ {
			if s1.seq[i] < s2.seq[i] {
				return false
			} else if s1.seq[i] > s2.seq[i] {
				return true
			}
		}
	}

	return false
}

func parseCfg(fn string) (cfgInfo CfgInfo, e error) {
	var inFile *os.File
	var err error
	if inFile, err = os.Open(fn); err != nil {
		log.Fatal(err)
	}
	var libInfo LibInfo
	defer inFile.Close()
	reader := bufio.NewReader(inFile)
	eof := false
	for !eof {
		var line string
		line, err = reader.ReadString('\n')
		if err == io.EOF {
			err = nil
			eof = true
		} else if err != nil {
			log.Fatal(err)
		}
		fields := strings.Fields(line)
		if len(fields) == 0 {
			continue
		}
		var v int
		switch fields[0] {
		case "[global_setting]":
		case "[LIB]":
			if libInfo.name != "" {
				cfgInfo.libs = append(cfgInfo.libs, libInfo)
				var nli LibInfo
				libInfo = nli
			}
		case "max_rd_len":
			v, err = strconv.Atoi(fields[2])
			cfgInfo.maxRdLen = v
		case "min_rd_len":
			v, err = strconv.Atoi(fields[2])
			cfgInfo.minRdLen = v
		case "name":
			libInfo.name = fields[2]
		case "avg_insert_len":
			v, err = strconv.Atoi(fields[2])
			libInfo.insertSize = v
		case "insert_SD":
			v, err = strconv.Atoi(fields[2])
			libInfo.insertSD = v
		case "diverse_rd_len":
			v, err = strconv.Atoi(fields[2])
			libInfo.diverse = uint8(v)
		case "asm_flag":
			v, err = strconv.Atoi(fields[2])
			libInfo.asmFlag = uint8(v)
		case "seq_profile":
			v, err = strconv.Atoi(fields[2])
			libInfo.seqProfile = uint8(v)
		case "qual_benchmark":
			v, err = strconv.Atoi(fields[2])
			libInfo.qualBenchmark = uint8(v)
		case "f1", "f2":
			libInfo.fnName = append(libInfo.fnName, fields[2])
		default:
			if fields[0][0] != '#' && fields[0][0] != ';' {
				log.Fatalf("noknown line: %s", line)
			}
		}
		if err != nil {
			e = err
			return
		}
	}
	if libInfo.name != "" {
		cfgInfo.libs = append(cfgInfo.libs, libInfo)
	}

	return
}

func ExtendReadBnt2Byte(rb ReadBnt) (extRB ReadBnt) {
	extRB.length = rb.length
	extRB.seq = make([]byte, extRB.length)
	var crb ReadBnt
	crb.length = rb.length
	crb.seq = make([]byte, len(rb.seq))
	copy(crb.seq, rb.seq)

	for i := crb.length - 1; i >= 0; i-- {
		base := crb.seq[i/bnt.NumBaseInByte] & bnt.BaseMask
		extRB.seq[i] = base
		crb.seq[i/bnt.NumBaseInByte] >>= bnt.NumBitsInBase
	}

	return extRB
}

func GetReadBntKmer(extRBnt ReadBnt, startPos int, kmerlen int) (rbk ReadBnt) {
	rbk.length = kmerlen
	rbk.seq = make([]byte, (rbk.length+bnt.NumBaseInByte-1)/bnt.NumBaseInByte)

	for i := 0; i < kmerlen; i++ {
		rbk.seq[i/bnt.NumBaseInByte] <<= bnt.NumBitsInBase
		rbk.seq[i/bnt.NumBaseInByte] |= extRBnt.seq[i+startPos]
	}

	return rbk
}

func ReverseComplet(ks ReadBnt) (rs ReadBnt) {
	tmp := ks
	tmp.seq = make([]byte, len(ks.seq))
	copy(tmp.seq, ks.seq)
	rs.length = tmp.length
	rs.seq = make([]byte, len(tmp.seq))
	for i := tmp.length - 1; i >= 0; i-- {
		base := tmp.seq[i/bnt.NumBaseInByte] & bnt.BaseMask
		rs.seq[(rs.length-i-1)/bnt.NumBaseInByte] <<= bnt.NumBitsInBase
		rs.seq[(rs.length-i-1)/bnt.NumBaseInByte] |= (^base & bnt.BaseMask)
		tmp.seq[i/bnt.NumBaseInByte] >>= bnt.NumBitsInBase
	}
	return rs
}

func paraConstructCF(cf cuckoofilter.CuckooFilter, cs chan ReadSeqBucket, wc chan ReadSeqBucket) {

	var wrsb ReadSeqBucket
	for {
		rsb := <-cs
		if rsb.end == true {
			wc <- wrsb
			wc <- rsb
			break
		} else {
			if rsb.count < ReadSeqSize {
				fmt.Printf("rsb.readBuf length is :%d\n", len(rsb.readBuf))
			}
			for i := 0; i < rsb.count; i++ {
				extRBnt := ExtendReadBnt2Byte(rsb.readBuf[i])
				for j := 0; j < extRBnt.length-cf.Kmerlen+1; j++ {
					ks := GetReadBntKmer(extRBnt, j, cf.Kmerlen)
					rs := ReverseComplet(ks)
					if ks.BiggerThan(rs) {
						ks, rs = rs, ks
					}
					count, suc := cf.Insert(ks.seq)
					if suc == false {
						log.Fatal("Insert to the CuckooFilter false")
					}
					if count == 1 {
						if wrsb.count >= ReadSeqSize {
							wc <- wrsb
							var nrsb ReadSeqBucket
							wrsb = nrsb
						}
						wrsb.readBuf[wrsb.count] = ks
						wrsb.count++
					}
				}
			}
		}
	}
}

func writeKmer(wrfn string, we chan int, wc chan ReadSeqBucket, numCPU int) {
	outfp, err := os.Create(wrfn)
	if err != nil {
		log.Fatal(err)
	}
	defer outfp.Close()
	gzwriter := gzip.NewWriter(outfp)
	if err != nil {
		log.Fatal(err)
	}
	defer gzwriter.Close()

	endFlagCount := 0
	for {
		rsb := <-wc
		if rsb.end == true {
			endFlagCount++
			if endFlagCount == numCPU {
				we <- 1
				break
			} else {
				continue
			}
		}

		for i := 1; i < rsb.count; i++ {
			gzwriter.Write(rsb.readBuf[i].seq)
		}
	}
}

func Trans2Byte(s string) (rb ReadBnt) {
	rb.length = len(s)
	rb.seq = make([]byte, (rb.length+bnt.NumBaseInByte-1)/bnt.NumBaseInByte)
	for i := 0; i < rb.length; i++ {
		b := bnt.Base2Bnt[s[i]]
		if b > 3 {
			fmt.Printf("found input sequence base '%c' not belong 'ACTG/actg', please check\n", s[i])
			log.Fatal("error found")
		}
		rb.seq[i/bnt.NumBaseInByte] <<= bnt.NumBitsInBase
		rb.seq[i/bnt.NumBaseInByte] |= b
	}

	return rb
}

func CCF(c cli.Command) {
	fmt.Println(c.Flags(), c.Parent().Flags())

	//argsCheck(c)
	fnName := c.Parent().Flag("C").String()
	cfgInfo, err := parseCfg(string(fnName))
	if err != nil {
		log.Fatal("[parseCfg] found err")
	}
	fmt.Println(cfgInfo)

	// make CuckooFilter
	var size, kmerlen int
	if size, err = strconv.Atoi(c.Flag("S").String()); err != nil {
		log.Fatal("flag 'S' set error")
	}
	if kmerlen, err = strconv.Atoi(c.Parent().Flag("K").String()); err != nil {
		log.Fatal("flag 'K' set error")
	}
	t0 := time.Now()
	cf := cuckoofilter.MakeCuckooFilter(uint64(size), kmerlen)
	numCPU, err := strconv.Atoi(c.Parent().Flag("t").String())
	runtime.GOMAXPROCS(numCPU + 2)
	bufsize := 10
	cs := make(chan ReadSeqBucket, bufsize)
	wc := make(chan ReadSeqBucket, bufsize)
	we := make(chan int)
	defer close(we)
	defer close(wc)
	defer close(cs)
	for i := 0; i < numCPU; i++ {
		go paraConstructCF(cf, cs, wc)
	}
	// write goroutinue
	wrfn := c.Parent().Flag("p").String() + ".kmerseq.gz"
	go writeKmer(wrfn, we, wc, numCPU)

	var processNumReads int
	var rsb ReadSeqBucket
	//var bucketCount int
	// iteration read cfgInfo files
	for _, lib := range cfgInfo.libs {
		if lib.asmFlag == AllState {
			for _, fn := range lib.fnName {
				infile, err := os.Open(fn)
				if err != nil {
					log.Fatal(err)
				}
				defer infile.Close()
				gzreader, err := gzip.NewReader(infile)
				bufin := bufio.NewReader(gzreader)
				defer func() {
					if err != nil {
						log.Fatal(err)
					}
					gzreader.Close()
				}()

				eof := false
				count := 0
				for !eof {
					var line string
					line, err = bufin.ReadString('\n')
					if err != nil {
						if err == io.EOF {
							err = nil
							eof = true
						} else {
							log.Fatal(err)
						}
					}

					if count%4 == 1 {
						s := strings.TrimSpace(line)
						bs := Trans2Byte(s)
						if rsb.count >= ReadSeqSize {
							cs <- rsb
							var nsb ReadSeqBucket
							rsb = nsb
						}
						rsb.readBuf[rsb.count] = bs
						rsb.count++
						//fmt.Printf("%s\n", s)
						processNumReads++
						if processNumReads%(1024*1024) == 0 {
							fmt.Printf("processed reads number: %d\n", processNumReads)
						}
					}
					count++
				}
			}
		}
	}
	cs <- rsb
	// send read finish signal
	for i := 0; i < numCPU; i++ {
		var nrsb ReadSeqBucket
		nrsb.end = true
		cs <- nrsb
	}

	<-we // end signal from wirte goroutinue

	// output stat
	cf.GetStat()
	t1 := time.Now()
	fmt.Printf("construct CuckooFilter took %v to run\n", t1.Sub(t0))

	if err != nil {
		log.Fatal("error found in the CCF")
	}
}
