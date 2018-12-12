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
	"fmt"
	"io"
	"log"
	"math"
	"os"
	//"log"

	//"strings"
	"encoding/gob"

	"github.com/mudesheng/highwayhash"
	// "syscall"
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

// KEY used for highwayhash function, must len(KEY) == 32
//var KEY = []byte{35, 158, 189, 243, 123, 39, 95, 219, 58, 253, 127, 163, 91, 235, 248, 177, 139, 67, 229, 171, 195, 81, 95, 149, 191, 249, 148, 45, 155, 235}
var KEY = []uint64{0xBD4CCC325BEFCA6F, 0xA89A58CE65E641FF, 0xAE093FEF1F84E3E7, 0xFB4297E8C586EE2D}

func CompareAndSwapUint16(addr *uint16, old uint16, new uint16) (swapped bool) {
	a := (*C.uint16_t)(addr)
	return C.CompareAndSwapUint16(a, C.uint16_t(old), C.uint16_t(new)) == C.uint16_t(old)
}

type CFItem uint16

/* type cfitem struct {
  fingerprint uint16:NUM_FP_BITS
  count uint16:NUM_C_BITS
} */

type Bucket struct {
	Bkt [BucketSize]CFItem
}

type CuckooFilter struct {
	Hash     []Bucket
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

// "MakeCuckooFilter is for construct Cuckoo Filter"
func MakeCuckooFilter(maxNumKeys uint64, kmerLen int) (cf CuckooFilter) {
	numBuckets := upperpower2(maxNumKeys) / BucketSize
	/*frac := float64(maxNumKeys) / numBuckets / BucketSize
	if frac > MaxLoad {
		numBuckets <<= 1
	} */

	cf.Hash = make([]Bucket, numBuckets)
	cf.numItems = numBuckets * BucketSize
	cf.Kmerlen = kmerLen
	fmt.Printf("[MakeCuckooFilter]cf items number is: %d\n", cf.numItems)

	return cf

}

func (cf CuckooFilter) IndexHash(v uint64) uint64 {
	v = (v >> 37) ^ (v >> 27) ^ v

	return v % uint64(len(cf.Hash))
}

func FingerHash(v uint64) uint16 {
	v = (v >> 47) ^ (v >> 33) ^ (v >> 19) ^ (v >> 13) ^ v
	v &= FPMASK
	return uint16(v)
}

func (cf CuckooFilter) AltIndex(index uint64, finger uint16) uint64 {
	index ^= uint64(finger)

	return index % uint64(len(cf.Hash))
}

func combineCFItem(fp uint16, count uint16) (cfi CFItem) {
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
			a := (*uint16)(cfi)
			if CompareAndSwapUint16(a, uint16(oc), uint16(nc)) {
				return int(count), 0, true
			}
		} else {
			return MAX_C, 0, true
		}
	}
}

func (b Bucket) Contain(fingerprint uint16) bool {
	for _, item := range b.Bkt {
		if item.GetCount() > 0 && item.GetFinger() == fingerprint {
			return true
		}
	}

	return false
}

// return kickout CFItem and bool if successed added and new added fingerprint count before added
func (b *Bucket) AddBucket(cfi CFItem, kickout bool) (CFItem, bool, int) {
	for i := 0; i < BucketSize; i++ {
		//fmt.Printf("i: %d\n", i)
		for {
			oi := b.Bkt[i]
			if oi.GetCount() == 0 {
				a := (*uint16)(&b.Bkt[i])
				if CompareAndSwapUint16(a, uint16(oi), uint16(cfi)) {
					//addr := (*C.uint16_t)(&b.Bkt[i])
					//C.CompareAndSwapUint16(addr, C.uint16_t(oi), C.uint16_t(cfi)) == C.uint16_t(oi) {
					//fmt.Printf("[AddBucket]cfi.finger: %v, i: %d, finger: %v\n", cfi.GetFinger(), i, b.Bkt[i].GetFinger())
					countItems++
					return CFItem(0), true, 0
				}
			} else {
				if oi.GetCount() > 0 && b.Bkt[i].EqualFP(cfi) {
					oc, _, _ := b.Bkt[i].AddCount()
					return CFItem(0), true, oc
				} else {
					break
				}
			}
		}
		//fmt.Printf("\tafter: %d\n", b.Bkt[i])
		//if added == true {
		//	break
		//}
	}

	//fmt.Printf("kikcout: %t", kickout)
	if kickout {
		min := uint16(math.MaxUint16)
		idx := -1
		for j := BucketSize - 1; j >= 0; j-- {
			c := b.Bkt[j].GetCount()
			if c < min {
				min = c
				idx = j
			}
		}
		var oi CFItem
		for {
			oi = b.Bkt[idx]
			a := (*uint16)(&b.Bkt[idx])
			if CompareAndSwapUint16(a, uint16(oi), uint16(cfi)) {
				break
			}
		}
		return oi, true, 0
	} else {
		return CFItem(0), false, 0
	}
}

// return if successed added
func (cf CuckooFilter) Add(index uint64, fingerprint uint16) (oldcount int, succ bool) {
	ci := index
	cfi := combineCFItem(fingerprint, 1)
	for count := 0; count < KMaxCount; count++ {
		kickout := count > 0
		b := &cf.Hash[ci]
		old, added, oc := b.AddBucket(cfi, kickout)
		if count <= 1 { // add the new fingerprint, set oldcount
			oldcount = oc
		}
		if added == true && old == 0 {
			succ = true
			return oldcount, succ
		}
		if old.GetCount() > 0 {
			cfi = old
		}
		//fmt.Printf("cycle : %d\n", count)

		ci = cf.AltIndex(ci, cfi.GetFinger())
	}
	return oldcount, succ
}

/*func hk2uint64(hk [sha1.Size]byte) (v uint64) {
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
}*/

// return last count of kmer fingerprint and have been successed inserted
func (cf CuckooFilter) Insert(kb []uint64) (int, bool) {
	//hk := sha1.Sum(kb)
	//v := hk2uint64(hk)
	v := highwayhash.SumInput64Arr64(kb, KEY)

	//fmt.Printf("%v\t", v)
	index := cf.IndexHash(v)
	fingerprint := FingerHash(v)
	//fmt.Printf("[cf.Insert]%v\t%v\n", index, fingerprint)
	//fmt.Printf(" sizeof cuckoofilter.Hash[0] : %d\n", unsafe.Sizeof(cf.Hash[0]))

	return cf.Add(index, fingerprint)
}

func (cf CuckooFilter) Lookup(kb []uint64) bool {
	//hk := sha1.Sum(kb)
	//v := hk2uint64(hk)
	v := highwayhash.SumInput64Arr64(kb, KEY)
	index := cf.IndexHash(v)
	fingerprint := FingerHash(v)

	if cf.Hash[index].Contain(fingerprint) {
		return true
	} else {
		index = cf.AltIndex(index, fingerprint)
		return cf.Hash[index].Contain(fingerprint)
	}
}

func (cf CuckooFilter) GetCount(kb []uint64) uint16 {
	//hk := sha1.Sum(kb)
	//v := hk2uint64(hk)
	v := highwayhash.SumInput64Arr64(kb, KEY)
	index := cf.IndexHash(v)
	fingerprint := FingerHash(v)
	for _, item := range cf.Hash[index].Bkt {
		// fmt.Printf("index: %v, finger: %v\n", index, item.GetFinger())
		if item.GetFinger() == fingerprint {
			return item.GetCount()
		}
	}
	// if not return , find another position
	index = cf.AltIndex(index, fingerprint)
	for _, item := range cf.Hash[index].Bkt {
		// fmt.Printf("index: %v, finger: %v\n", index, item.GetFinger())
		if item.GetFinger() == fingerprint {
			return item.GetCount()
		}
	}

	// not found in the CuckooFilter
	panic("not found in the CuckooFilter")
}

/* allow function return zero version if kb not found in the CuckooFilter */
func (cf CuckooFilter) GetCountAllowZero(kb []uint64) uint16 {
	//hk := sha1.Sum(kb)
	//v := hk2uint64(hk)
	v := highwayhash.SumInput64Arr64(kb, KEY)
	// fmt.Printf("[GetCountAllowZero] len(cf.Hash): %d\n", len(cf.Hash))
	// fmt.Printf("[GetCountAllowZero] cf.Hash[0]: %v\n", cf.Hash[0])
	index := cf.IndexHash(v)
	fingerprint := FingerHash(v)
	for _, item := range cf.Hash[index].Bkt {
		// fmt.Printf("index: %v, finger: %v\n", index, item.GetFinger())
		if item.GetFinger() == fingerprint {
			return item.GetCount()
		}
	}
	// if not return , find another position
	index = cf.AltIndex(index, fingerprint)
	for _, item := range cf.Hash[index].Bkt {
		// fmt.Printf("index: %v, finger: %v\n", index, item.GetFinger())
		if item.GetFinger() == fingerprint {
			return item.GetCount()
		}
	}

	// not found in the CuckooFilter, return zero
	// panic("not found in the CuckooFilter")
	return 0
}

func (cf CuckooFilter) GetStat() {
	var ca [4]int
	for _, b := range cf.Hash {
		for _, e := range b.Bkt {
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
	/*n, err := cfmmapfp.Write([]byte(cf.Hash[0].Bkt[0]))
	if err != nil {
		return err
	}
	if n != cf.numItems*unsafe.Sizeof(bucket) {
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
	// fmt.Printf("[RecoverCuckooFilterInfo] line: %s\n", line)
	_, err = fmt.Sscanf(line, "numItems\t%d\n", &cf.numItems)
	if err != nil {
		return cf, err
	}
	line, err = cfinfobuf.ReadString('\n')
	// fmt.Printf("[RecoverCuckooFilterInfo] line: %s\n", line)
	_, err = fmt.Sscanf(line, "Kmerlen\t%d\n", &cf.Kmerlen)
	if err != io.EOF {
		return cf, err
	} else {
		err = nil
	}

	return
}

func MmapReader(cfmmapfn string) (cf CuckooFilter, err error) {
	// if cf.numItems <= 0 {
	// 	log.Fatal("CuckooFilter number is <=0, please check")
	// }
	cfmmapfp, err := os.Open(cfmmapfn)
	if err != nil {
		return cf, err
	}
	defer cfmmapfp.Close()
	dec := gob.NewDecoder(cfmmapfp)
	err = dec.Decode(&cf)
	if err != nil {
		log.Fatalf("[MmapReader] err : %v\n", err)
	}
	// fmt.Printf("[MmapReader] cf.Kmerlen: %d\n", cf.Kmerlen)
	// fmt.Printf("[MmapReader] len(cf.Hash): %d\n", len(cf.Hash))
	// fmt.Printf("[MmapReader] cf.Hash[0]: %d\n", cf.Hash[0])

	// mmap, err := syscall.Mmap(cfmmapfp.Fd(), 0, cf.numItems*unsafe.Sizeof(cf.Hash[0].Bkt[0]), syscall.PROT_READ, syscall.MAP_SHARED)
	// if err != nil {
	// 	log.Fatal(err)
	// }
	// cf.Hash = []Bucket(mmap)

	return cf, nil
}
