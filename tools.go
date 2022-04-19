package main

import (
	"bytes"
	"fmt"
	"log"
	"math/rand"
	"os"
	"runtime"
	"strconv"
	"strings"
	"sync"

	"github.com/jwaldrip/odin/cli"
	"github.com/klauspost/compress/zstd"
)

type optionsAONTQC struct {
	Input  string
	Output string
	Paf    string
}

type optionsSimulateNGS struct {
	Input  string
	Output string
	Step   int
	ID     int
}

type MQ struct {
	ID int
	QC int
}

func checkArgsAddONTQC(c cli.Command) (opt optionsAONTQC, suc bool) {
	suc = true
	opt.Input = c.Flag("input").Get().(string)
	opt.Output = c.Flag("output").Get().(string)
	opt.Paf = c.Flag("paf").Get().(string)
	return opt, suc
}

func checkArgsSimulateNGS(c cli.Command) (opt optionsSimulateNGS, suc bool) {
	suc = true
	opt.Input = c.Flag("input").Get().(string)
	opt.Output = c.Flag("output").Get().(string)
	opt.Step = c.Flag("Step").Get().(int)
	opt.ID = c.Flag("ID").Get().(int)
	return opt, suc
}

func GetMQ(s string) (mq MQ) {
	sa := strings.Split(s, "\t")
	//fmt.Printf("[GetMQ]s:%v\n", s)
	mq.ID, _ = strconv.Atoi(sa[0])
	ml, _ := strconv.Atoi(sa[9])
	tl, _ := strconv.Atoi(sa[10])
	mq.QC = ml * 100 / tl
	return
}

func GetMQArr(cs <-chan []byte, rbytesPool *sync.Pool) []MQ {
	arr := make([]MQ, 0, 1000)
	notFoundNewRecord := make([]byte, 0, 500)
	for {
		bs, ok := <-cs
		if !ok {
			if len(notFoundNewRecord) != 0 {
				log.Fatalf("[GetMQArr] to the file end, but len(notFoundNewRecord) = %d\n", len(notFoundNewRecord))
			}
			break
		}
		sz := bytes.IndexByte(bs, '\n')
		if sz < 0 {
			notFoundNewRecord = append(notFoundNewRecord, bs...)
			continue
		} else {
			notFoundNewRecord = append(notFoundNewRecord, bs[:sz]...)
		}
		mq := GetMQ(string(notFoundNewRecord))
		if len(arr) > 0 {
			if mq.ID != arr[len(arr)-1].ID {
				arr = append(arr, mq)
			}
		} else {
			arr = append(arr, mq)
		}

		notFoundNewRecord = notFoundNewRecord[:0]
		sz++
		for sz >= 0 && sz < len(bs) {
			idx := bytes.IndexByte(bs[sz:], '\n')
			if idx < 0 {
				notFoundNewRecord = append(notFoundNewRecord, bs[sz+1:]...)
				break
			}
			mq := GetMQ(string(bs[sz : sz+idx]))
			if mq.ID != arr[len(arr)-1].ID {
				arr = append(arr, mq)
			}
			sz += idx + 1
		}
	}
	fmt.Printf("[GetMQArr]MQArr:%v\n", arr)
	return arr
}

func AddMQ(input, output string, cs <-chan []byte, rbytesPool *sync.Pool, MQArr []MQ) (readNum int) {
	format := GetReadsFileFormat(input)
	var notFoundNewRecord []byte
	fmt.Printf("[AddMQ] begin processe file: %v...\n", input)
	riArrPool := sync.Pool{New: func() interface{} {
		arr := make([]ReadInfo, WindowSize/10000)
		return arr
	}}
	outfp, err := os.Create(output)
	if err != nil {
		log.Fatal(err)
	}
	defer outfp.Close()

	zw, err1 := zstd.NewWriter(outfp, zstd.WithEncoderCRC(false), zstd.WithEncoderConcurrency(1), zstd.WithEncoderLevel(1))
	if err1 != nil {
		log.Fatalf("[AddMQ]open write file:%s err: %v\n", output, err)
	}
	defer zw.Close()
	writeCount := 0
	mqArrIdx := 0
	for {
		bs, ok := <-cs
		if !ok {
			if len(notFoundNewRecord) != 0 {
				log.Fatalf("[AddMQ] to the file end, but len(notFoundNewRecord) = %d\n", len(notFoundNewRecord))
			}
			break
		}

		var rp RIPool
		rp.RIArr = riArrPool.Get().([]ReadInfo)
		rp.RIArr = rp.RIArr[:0]
		rp.Cs = bs
		rp.NotFoundRecordBytes = notFoundNewRecord
		idx := GetReadFileReadInfo(&rp, format, false, false)
		//fmt.Printf("[GetMQArr]idx:%d rp:%v\n", idx, rp)
		if idx < 0 {
			notFoundNewRecord = append(notFoundNewRecord, rp.Cs...)
			rbytesPool.Put(rp.Cs)
			riArrPool.Put(rp.RIArr)
			continue
		}
		readNum += len(rp.RIArr)
		var nf []byte
		if idx < len(rp.Cs) {
			nf = append(nf, rp.Cs[idx:]...)
			//sp.Cs = sp.Cs[:idx]
		}
		notFoundNewRecord = nf
		//fmt.Printf("[GetReadSeqBucket] processed reads number is: %d,finished processed file: %v\n", len(sp.Cs),len(notFoundNewRecord), len(sp.SeqArr))
		for i := 0; i < len(rp.RIArr); i++ {
			ri := rp.RIArr[i]
			if ri.ID != MQArr[mqArrIdx].ID {
				log.Fatalf("[AddMQ]ri.ID:%d MQArr[%d].ID:%d\n", ri.ID, mqArrIdx, MQArr[mqArrIdx].ID)
			}
			fmt.Fprintf(zw, ">%d\tMQ=%d\n%s\n", ri.ID, MQArr[mqArrIdx].QC, ri.Seq)
			mqArrIdx++
			writeCount++
		}
		rbytesPool.Put(rp.Cs)
		riArrPool.Put(rp.RIArr)
	}
	// send read finish signal
	fmt.Printf("[AddMQ] processed reads number is: %d,finished processed file: %v\n", readNum, input)
	return
}

func SimulateReads(input, output string, cs <-chan []byte, rbytesPool *sync.Pool, step, id int) (writeCount int) {
	format := GetReadsFileFormat(input)
	var notFoundNewRecord []byte
	fmt.Printf("[SimulateReads] begin processe file: %v...\n", input)
	riArrPool := sync.Pool{New: func() interface{} {
		arr := make([]ReadInfo, WindowSize/10000)
		return arr
	}}
	outfp, err := os.Create(output)
	if err != nil {
		log.Fatal(err)
	}
	defer outfp.Close()

	zw, err1 := zstd.NewWriter(outfp, zstd.WithEncoderCRC(false), zstd.WithEncoderConcurrency(1), zstd.WithEncoderLevel(1))
	if err1 != nil {
		log.Fatalf("[SimulateReads]open write file:%s err: %v\n", output, err)
	}
	defer zw.Close()
	readNum := 0
	for {
		bs, ok := <-cs
		if !ok {
			if len(notFoundNewRecord) != 0 {
				log.Fatalf("[SimulateReads] to the file end, but len(notFoundNewRecord) = %d\n", len(notFoundNewRecord))
			}
			break
		}

		var rp RIPool
		rp.RIArr = riArrPool.Get().([]ReadInfo)
		rp.RIArr = rp.RIArr[:0]
		rp.Cs = bs
		rp.NotFoundRecordBytes = notFoundNewRecord
		idx := GetReadFileReadInfo(&rp, format, false, false)
		fmt.Printf("[SimulateReads]idx:%d len(rp.NotFoundRecordBytes):%d\n", idx, len(rp.NotFoundRecordBytes))
		if idx <= 0 {
			notFoundNewRecord = append(notFoundNewRecord, rp.Cs...)
			rbytesPool.Put(rp.Cs)
			riArrPool.Put(rp.RIArr)
			continue
		}
		readNum += len(rp.RIArr)
		var nf []byte
		if idx < len(rp.Cs) {
			nf = append(nf, rp.Cs[idx:]...)
			//sp.Cs = sp.Cs[:idx]
		}
		notFoundNewRecord = nf
		fmt.Printf("[SimulateReads]ID:%d len:%d\n", rp.RIArr[0].ID, len(rp.RIArr[0].Seq))
		fmt.Printf("[SimulateReads]size:%d\n", len(rp.RIArr))
		for i := 0; i < len(rp.RIArr); i++ {
			ri := rp.RIArr[i]
			rSeq := GetReverseCompNtByteArr(ri.Seq)
			librarySize := 350 + rand.Intn(100)
			for j := 0; j+librarySize < len(ri.Seq); j += step {
				//fmt.Printf("[SimulateReads]librarySize:%d\n", librarySize)
				fmt.Fprintf(zw, ">%d\n%s\n>%d\n%s\n", id, ri.Seq[j:j+250], id, rSeq[len(rSeq)-(j+librarySize):len(rSeq)-(j+librarySize)+250])
				writeCount++
				id++
				librarySize = 350 + rand.Intn(100)
			}
		}
		rbytesPool.Put(rp.Cs)
		riArrPool.Put(rp.RIArr)
	}
	// send read finish signal
	fmt.Printf("[SimulateReads] processed reads number is: %d,finished processed file: %v\n", readNum, input)
	return
}

func SplitSeq(input, output string, cs <-chan []byte, rbytesPool *sync.Pool) (readNum int) {
	format := GetReadsFileFormat(input)
	var notFoundNewRecord []byte
	fmt.Printf("[SplitSeq] begin processe file: %v...\n", input)
	riArrPool := sync.Pool{New: func() interface{} {
		arr := make([]ReadInfo, WindowSize/10000)
		return arr
	}}
	outfp, err := os.Create(output)
	if err != nil {
		log.Fatal(err)
	}
	defer outfp.Close()

	zw, err1 := zstd.NewWriter(outfp, zstd.WithEncoderCRC(false), zstd.WithEncoderConcurrency(1), zstd.WithEncoderLevel(1))
	if err1 != nil {
		log.Fatalf("[SplitSeq]open write file:%s err: %v\n", output, err)
	}
	defer zw.Close()
	totalLen := 0
	id := 1
	var finished bool
	for {
		bs, ok := <-cs
		if !ok {
			if len(notFoundNewRecord) != 0 {
				log.Fatalf("[SplitSeq] to the file end, but len(notFoundNewRecord) = %d\n", len(notFoundNewRecord))
			}
			break
		}

		var rp RIPool
		rp.RIArr = riArrPool.Get().([]ReadInfo)
		rp.RIArr = rp.RIArr[:0]
		rp.Cs = bs
		rp.NotFoundRecordBytes = notFoundNewRecord
		idx := GetReadFileReadInfo(&rp, format, false, false)
		//fmt.Printf("[GetMQArr]idx:%d rp:%v\n", idx, rp)
		if idx < 0 {
			notFoundNewRecord = append(notFoundNewRecord, rp.Cs...)
			rbytesPool.Put(rp.Cs)
			riArrPool.Put(rp.RIArr)
			continue
		}
		readNum += len(rp.RIArr)
		if idx < len(rp.Cs) {
			notFoundNewRecord = append(notFoundNewRecord[:0], rp.Cs[idx:]...)
			//sp.Cs = sp.Cs[:idx]
		}
		//fmt.Printf("[GetReadSeqBucket] processed reads number is: %d,finished processed file: %v\n", len(sp.Cs),len(notFoundNewRecord), len(sp.SeqArr))
		for i := 0; i < len(rp.RIArr); i++ {
			ri := rp.RIArr[i]
			sIdx := 0
			for sIdx < len(ri.Seq) {
				sl := rand.Intn(550000) + 10000
				end := sIdx + sl
				if end > len(ri.Seq) {
					break
				}
				fmt.Fprintf(zw, ">%d\n%s\n", id, ri.Seq[sIdx:end])
				totalLen += sl
				if totalLen > 5500000000 {
					finished = true
					break
				}
				sIdx += sl / 2
				id++
			}
			if finished {
				break
			}
		}
		rbytesPool.Put(rp.Cs)
		riArrPool.Put(rp.RIArr)
		if finished {
			break
		}
	}
	// send read finish signal
	fmt.Printf("[SplitSeq] processed reads number is: %d,finished processed file: %v\n", readNum, input)
	return
}

func AddONTQC(c cli.Command) {
	//fmt.Println(c.Flags(), c.Parent().Flags())
	opt, suc := checkArgsAddONTQC(c)
	if suc == false {
		log.Fatalf("[AddONTQC] check Arguments error, opt:%v\n", opt)
	}
	fmt.Printf("[AddONTQC]opt:%v\n", opt)

	//t0 := time.Now()
	runtime.GOMAXPROCS(2)
	bufSize := WindowSize

	rbytesPool := sync.Pool{New: func() interface{} {
		bytes := make([]byte, bufSize)
		return bytes
	}}

	cs := make(chan []byte, 4)
	if opt.Paf == "" {
		go ReadZstdFile(opt.Input, &rbytesPool, cs)
		SplitSeq(opt.Input, opt.Output, cs, &rbytesPool)
		return
	}
	// parse paf file
	go ReadZstdFile(opt.Paf, &rbytesPool, cs)
	MQArr := GetMQArr(cs, &rbytesPool)

	// process seq file
	cs1 := make(chan []byte, 4)
	go ReadZstdFile(opt.Input, &rbytesPool, cs1)
	AddMQ(opt.Input, opt.Output, cs1, &rbytesPool, MQArr)
	return
}

func SimulateNGS(c cli.Command) {
	//fmt.Println(c.Flags(), c.Parent().Flags())
	opt, suc := checkArgsSimulateNGS(c)
	if suc == false {
		log.Fatalf("[SimulateNGS] check Arguments error, opt:%v\n", opt)
	}
	fmt.Printf("[SimulateNGS]opt:%v\n", opt)

	//t0 := time.Now()
	runtime.GOMAXPROCS(2)
	bufSize := 2 << 22

	rbytesPool := sync.Pool{New: func() interface{} {
		bytes := make([]byte, bufSize)
		return bytes
	}}

	// simulate NGS reads
	cs := make(chan []byte, 4)
	go ReadZstdFile(opt.Input, &rbytesPool, cs)
	num := SimulateReads(opt.Input, opt.Output, cs, &rbytesPool, opt.Step, opt.ID)
	fmt.Printf("[SimulateNGS]write output file number:%d\n", num)
	return
}
