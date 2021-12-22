package main

import (
	"errors"
	"log"
	"math"
	"unsafe"

	"github.com/jwaldrip/odin/cli"
	//"fmt"
	//"os"
)

//DbgMaxInt use marker DBG node ID and DBG edge
//type DbgMaxInt uint32
const MaxKmerAllow = NODEMAP_KEY_LEN * 32

type ArgsOpt struct {
	Prefix     string
	Kmer       int
	NumCPU     int
	CfgFn      string
	Cpuprofile string
}

// return global arguments and check if successed
func CheckGlobalArgs(c cli.Command) (opt ArgsOpt, succ bool) {
	opt.Prefix = c.Flag("p").String()
	if opt.Prefix == "" {
		log.Fatalf("[CheckGlobalArgs] args 'p' not set\n")
	}
	opt.CfgFn = c.Flag("C").String()
	if opt.CfgFn == "" {
		log.Fatalf("[CheckGlobalArgs] args 'C' not set\n")
	}
	opt.Cpuprofile = c.Flag("cpuprofile").String()
	if opt.Cpuprofile == "" {
		log.Fatalf("[CheckGlobalArgs] args 'cpuprofile' not set\n")
	}

	var ok bool
	opt.Kmer, ok = c.Flag("K").Get().(int)
	if !ok {
		log.Fatalf("[CheckGlobalArgs] args 'K' : %v set error\n", c.Flag("K").String())
	}
	if opt.Kmer >= NODEMAP_KEY_LEN*32 {
		log.Fatalf("[CheckGlobalArgs] the argument 'K':%d must small than [NODEMAP_KEY_LEN * 32]:%d\n", opt.Kmer, NODEMAP_KEY_LEN*32)
	} else if opt.Kmer%2 != 1 {
		log.Fatalf("[CheckGlobalArgs] the argument 'K':%d must odd Number\n", opt.Kmer, NODEMAP_KEY_LEN*32)
	}
	opt.NumCPU, ok = c.Flag("t").Get().(int)
	if !ok {
		log.Fatalf("[CheckGlobalArgs] args 't': %v set error\n", c.Flag("t").String())
	}
	return opt, true
}

func AbsInt(a int) int {
	if a < 0 {
		return -a
	} else {
		return a
	}
}

func MaxInt(a, b int) int {
	if a > b {
		return a
	} else {
		return b
	}
}

func MaxUint8(a, b uint8) uint8 {
	if a > b {
		return a
	} else {
		return b
	}
}

func MinInt(a, b int) int {
	if a > b {
		return b
	} else {
		return a
	}
}

func MaxUint32(a, b uint32) uint32 {
	if a > b {
		return a
	} else {
		return b
	}
}

func MinUint32(a, b uint32) uint32 {
	if a > b {
		return b
	} else {
		return a
	}
}

func Min(x, y int) int {
	if x < y {
		return x
	}
	return y
}

func ByteArrInt(id []byte) (d int, err error) {
	//err := nil
	for _, c := range id {
		if c < '0' || c > '9' {
			err = errors.New("can't convert to digit...")
			return d, err
		}
		d = d*10 + int(c-'0')
	}
	return d, nil
}

func Bytes2String(b []byte) string {
	return *(*string)(unsafe.Pointer(&b))
}

func BytesEqual(a, b []byte) bool {
	if len(a) != len(b) {
		return false
	}
	return Bytes2String(a) == Bytes2String(b)
}

func BytesEqual2(a, b []byte) bool {
	if len(a) != len(b) {
		return false
	}
	for i := 0; i < len(a); i++ {
		if a[i] != b[i] {
			return false
		}
	}
	return true
}

func GetMinLenByteArr(ba [][]byte) (min int) {
	if len(ba) == 0 {
		min = 0
		return
	}

	min = math.MaxInt32
	for _, a := range ba {
		if len(a) < min {
			min = len(a)
		}
	}
	return
}

func CountEqualByte(a, b []byte, strand bool) (count int) {
	if strand {
		for i, c := range b {
			if a[i] == c {
				count++
			}
		}
	} else {
		al := len(a)
		for i, c := range b {
			if BntRev[a[al-1-i]] == c {
				count++
			}
		}
	}

	return
}

func GetReverseUint32Arr(arr []uint32) []uint32 {
	la := len(arr)
	rarr := make([]uint32, la)
	for i := la - 1; i >= 0; i-- {
		rarr[la-i-1] = arr[i]
	}

	return rarr

}

func GetReverseUint32Arr2(arr, rArr []uint32) []uint32 {
	la := len(arr)
	if len(rArr) < la {
		rArr = make([]uint32, la)
	} else {
		rArr = rArr[:la]
	}
	for i, id := range arr {
		rArr[la-1-i] = id
	}
	return rArr
}

func IndexUint32(arr []uint32, eID uint32) int {
	for i, id := range arr {
		if id == eID {
			return i
		}
	}
	return -1
}
