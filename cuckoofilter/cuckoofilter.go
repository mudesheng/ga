package cuckoofilter

// #include<stdint.h>
/*
uint16_t CompareAndSwapUint16(uint16_t *addr, uint16_t old, uint16_t new)
{
    return __sync_val_compare_and_swap(addr, old, new);
}*/
import "C"

import (
	"bufio"
	"crypto/sha1"
	"fmt"
	"io"
	"log"
	"os"
	//"log"
	//"sync/atomic"
	//"strings"
	"encoding/gob"
	// "syscall"
	"unsafe"
)

const (
	NUM_FP_BITS = 14     // number of fingerprint  bits occpied
	NUM_C_BITS  = 2      // number of count equal sizeof(uint16)*8 - NUM_FP_BITS
	FPMASK      = 0x3FFF // mask other info only set fingerprint = (1<<NUM_FP_BITS) -1
	CMASK       = 0x3    // set count bits field = (1<<NUM_C_BITS) -1
	MAX_C       = (1 << NUM_C_BITS) - 1
)

const BucketSize = 4
const MaxLoad = 0.95
const KMaxCount = 512

//func CompareAndSwapUint16(addr *uint16, old uint16, new uint16) (swapped bool)

type CFItem uint16

/* type cfitem struct {
  fingerprint uint16:NUM_FP_BITS
  count uint16:NUM_C_BITS
} */

type Bucket struct {
	bucket [BucketSize]CFItem
}

type CuckooFilter struct {
	hash     []Bucket
	numItems uint64
	Kmerlen  int
}

var countItems uint64

func upperpower2(x uint64) uint64 {
	x--
	x |= x >> 1
	x |= x >> 2
	x |= x >> 4
	x |= x >> 8
	x |= x >> 16
	x |= x >> 32
	x++

	return x
}

func MakeCuckooFilter(maxNumKeys uint64, kmerLen int) (cf CuckooFilter) {
	numBuckets := upperpower2(maxNumKeys) / BucketSize
	frac := float64(maxNumKeys / numBuckets / BucketSize)
	if frac > MaxLoad {
		numBuckets <<= 1
	}

	cf.hash = make([]Bucket, numBuckets)
	cf.numItems = numBuckets * BucketSize
	cf.Kmerlen = kmerLen

	return cf

}

func (cf CuckooFilter) IndexHash(v uint64) uint64 {
	v >>= NUM_FP_BITS

	return v % uint64(len(cf.hash))
}

func FingerHash(v uint64) uint64 {
	v &= FPMASK
	return v
}

func (cf CuckooFilter) AltIndex(index uint64, finger uint64) uint64 {
	index ^= finger

	return index % uint64(len(cf.hash))
}

func combineCFItem(fp uint64, count uint64) (cfi CFItem) {
	if count > MAX_C {
		panic("count bigger than CFItem allowed")
	}
	cfi = CFItem(fp)
	cfi <<= NUM_C_BITS
	cfi |= CFItem(count)
	//fmt.Printf("fp: %d, count: %d, cfi: %d\n", fp, count, uint16(cfi))
	return cfi
}

func (cfi CFItem) GetCount() uint16 {
	return uint16(cfi) & CMASK
}

func (cfi *CFItem) setCount(count uint16) {
	nc := uint16(*cfi) >> NUM_C_BITS
	nc <<= NUM_C_BITS
	nc |= count

	*cfi = CFItem(nc)
}

func (cfi CFItem) GetFinger() uint16 {
	return uint16(cfi >> NUM_C_BITS)
}

func (cfi CFItem) EqualFP(rcfi CFItem) bool {
	if (uint16(cfi) >> NUM_C_BITS) == (uint16(rcfi) >> NUM_C_BITS) {
		return true
	} else {
		return false
	}
}

// return oldcount, oldfinger, added
func (cfi *CFItem) AddCount() (int, uint64, bool) {
	for {
		oc := *cfi
		count := oc.GetCount()
		if count < MAX_C {
			nc := oc
			nc.setCount(count + 1)
			if C.CompareAndSwapUint16((*C.uint16_t)(cfi), C.uint16_t(oc), C.uint16_t(nc)) == C.uint16_t(oc) {
				return int(count), 0, true
			}
		} else {
			return MAX_C, 0, true
		}
	}
}

func (b Bucket) Contain(fingerprint uint16) bool {
	for _, item := range b.bucket {
		fp := item.GetFinger()
		if fp == fingerprint {
			return true
		}
	}

	return false
}

func (b *Bucket) AddBucket(fp uint64, kickout bool) (int, uint64, bool) {
	cfi := combineCFItem(fp, 1)
	for i := 0; i < BucketSize; i++ {
		//fmt.Printf("i: %d\n", i)
		for {
			oi := b.bucket[i]
			if oi.GetCount() == 0 {
				addr := (*C.uint16_t)(&b.bucket[i])
				if C.CompareAndSwapUint16(addr, C.uint16_t(oi), C.uint16_t(cfi)) == C.uint16_t(oi) {
					//fmt.Printf("[AddBucket]cfi.finger: %v, i: %d, finger: %v\n", cfi.GetFinger(), i, b.bucket[i].GetFinger())
					countItems++
					return 0, 0, true
				}
			} else {
				if b.bucket[i].EqualFP(cfi) {
					return b.bucket[i].AddCount()
				} else {
					break
				}
			}
		}
		//fmt.Printf("\tafter: %d\n", b.bucket[i])
		//if added == true {
		//	break
		//}
	}

	//fmt.Printf("kikcout: %t", kickout)
	if kickout {
		ci := fp & CMASK
		oldfinger := b.bucket[ci].GetFinger()
		b.bucket[ci] = cfi
		return 0, uint64(oldfinger), false
	} else {
		return 0, 0, false
	}
}

// return last count of kmer and if successed added
func (cf CuckooFilter) Add(index uint64, fingerprint uint64) (int, bool) {
	ci := index
	oldcount := -1
	cfinger := fingerprint
	for count := 0; count < KMaxCount; count++ {
		kickout := count > 0
		b := &cf.hash[ci]
		tc, oldfinger, added := b.AddBucket(cfinger, kickout)
		if added == true {
			if oldcount == -1 {
				oldcount = tc
			}
			return oldcount, true
		}
		if kickout {
			if oldcount == -1 {
				oldcount = tc
			}
			cfinger = oldfinger
		}
		//fmt.Printf("cycle : %d\n", count)

		ci = cf.AltIndex(ci, cfinger)
	}
	return -1, false
}

func hk2uint64(hk [sha1.Size]byte) (v uint64) {
	for i := 0; i <= len(hk)-8; i += 8 {
		hkp := unsafe.Pointer(&hk[i])
		t := (*uint64)(hkp)
		v ^= *t
	}
	if len(hk)%8 > 0 {
		hkp := unsafe.Pointer(&hk[len(hk)-8])
		t := (*uint64)(hkp)
		v ^= *t
	}

	return v
}

// return last count of kmer fingerprint and have been successed inserted
func (cf CuckooFilter) Insert(kb []byte) (int, bool) {
	hk := sha1.Sum(kb)
	v := hk2uint64(hk)
	//fmt.Printf("%v\t", v)
	index := cf.IndexHash(v)
	fingerprint := FingerHash(v)
	//fmt.Printf("[cf.Insert]%v\t%v\n", index, fingerprint)
	//fmt.Printf(" sizeof cuckoofilter.hash[0] : %d\n", unsafe.Sizeof(cf.hash[0]))

	return cf.Add(index, fingerprint)
}

func (cf CuckooFilter) Lookup(kb []byte) bool {
	hk := sha1.Sum(kb)
	v := hk2uint64(hk)
	index := cf.IndexHash(v)
	fingerprint := FingerHash(v)

	if cf.hash[index].Contain(uint16(fingerprint)) {
		return true
	} else {
		index = cf.AltIndex(index, fingerprint)
		return cf.hash[index].Contain(uint16(fingerprint))
	}
}

func (cf CuckooFilter) GetCount(kb []byte) uint16 {
	hk := sha1.Sum(kb)
	v := hk2uint64(hk)
	index := cf.IndexHash(v)
	fingerprint := FingerHash(v)
	for _, item := range cf.hash[index].bucket {
		fmt.Printf("index: %v, finger: %v\n", index, item.GetFinger())
		if item.GetFinger() == uint16(fingerprint) {
			return item.GetCount()
		}
	}
	// if not return , find another position
	index = cf.AltIndex(index, fingerprint)
	for _, item := range cf.hash[index].bucket {
		fmt.Printf("index: %v, finger: %v\n", index, item.GetFinger())
		if item.GetFinger() == uint16(fingerprint) {
			return item.GetCount()
		}
	}

	// not found in the CuckooFilter
	panic("not found in the CuckooFilter")
}

func (cf CuckooFilter) GetStat() {
	var ca [4]int
	for _, b := range cf.hash {
		for _, e := range b.bucket {
			count := e.GetCount()
			ca[count]++
		}
	}
	fmt.Printf("count statisticas : %v\n", ca)
	fmt.Printf("cuckoofilter numItems : %d, countItems: %d, load: %f\n", cf.numItems, countItems, float64(countItems)/float64(cf.numItems))
}

func (cf CuckooFilter) MmapWriter(cfmmapfn string) error {
	cfmmapfp, err := os.Create(cfmmapfn)
	if err != nil {
		return err
	}
	defer cfmmapfp.Close()

	// write cuckoofilter to the memory map file
	enc := gob.NewEncoder(cfmmapfp)
	err = enc.Encode(cf)
	if err != nil {
		return err
	}
	/*n, err := cfmmapfp.Write([]byte(cf.hash[0].bucket[0]))
	if err != nil {
		return err
	}
	if n != cf.numItems*unsafe.Sizeof(Bucket) {
		log.Fatal("write byte number is not equal cf size")
	}*/

	return nil
}

func (cf CuckooFilter) WriteCuckooFilterInfo(cfinfofn string) error {
	cfinfofp, err := os.Create(cfinfofn)
	if err != nil {
		return err
	}
	defer cfinfofp.Close()

	_, err = cfinfofp.WriteString(fmt.Sprintf("numItems\t%d\n", cf.numItems))
	if err != nil {
		return err
	}
	_, err = cfinfofp.WriteString(fmt.Sprintf("Kmerlen\t%d\n", cf.Kmerlen))
	if err != nil {
		return err
	}

	return nil
}

func RecoverCuckooFilterInfo(cfinfofn string) (cf CuckooFilter, err error) {
	var cfinfofp *os.File
	cfinfofp, err = os.Open(cfinfofn)
	if err != nil {
		return cf, err
	}
	defer cfinfofp.Close()
	cfinfobuf := bufio.NewReader(cfinfofp)
	var line string
	// eof := false
	line, err = cfinfobuf.ReadString('\n')
	if err != nil {
		if err == io.EOF {
			// eof = true
			err = nil
		} else {
			return cf, err
		}
	}
	// var num int
	_, err = fmt.Sscanf(line, "numItems\t%d\n", cf.numItems)
	if err != nil {
		return cf, err
	}
	line, err = cfinfobuf.ReadString('\n')
	_, err = fmt.Sscanf(line, "Kmerlen\t%d\n", cf.Kmerlen)
	if err != io.EOF {
		return cf, err
	} else {
		err = nil
	}

	return
}

func (cf CuckooFilter) MmapReader(cfmmapfn string) error {
	if cf.numItems <= 0 {
		log.Fatal("CuckooFilter number is <=0, please check")
	}
	cfmmapfp, err := os.Open(cfmmapfn)
	if err != nil {
		return err
	}
	defer cfmmapfp.Close()
	fmt.Println(cf)
	dec := gob.NewDecoder(cfmmapfp)
	err = dec.Decode(&cf)
	fmt.Println(cf)
	if err != nil {
		log.Fatal(err)
	}

	// mmap, err := syscall.Mmap(cfmmapfp.Fd(), 0, cf.numItems*unsafe.Sizeof(cf.hash[0].bucket[0]), syscall.PROT_READ, syscall.MAP_SHARED)
	// if err != nil {
	// 	log.Fatal(err)
	// }
	// cf.hash = []Bucket(mmap)

	return nil
}
