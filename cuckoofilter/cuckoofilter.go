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
	"math/rand"
	"os"
	"sync/atomic"
	"time"
	"unsafe"
	//"log"

	//"strings"
	"encoding/binary"
	"encoding/gob"

	"github.com/dgryski/go-metro"
	//"github.com/mudesheng/ga/cbrotli"
	"github.com/google/brotli/go/cbrotli"
	//"github.com/mudesheng/highwayhash"
	// "syscall"
)

const (
	NUM_FP_BITS = 13     // number of fingerprint  bits occpied
	NUM_C_BITS  = 3      // number of count equal sizeof(uint16)*8 - NUM_FP_BITS
	FPMASK      = 0x1FFF // mask other info only set fingerprint = (1<<NUM_FP_BITS) -1
	CMASK       = 0x7    // set count bits field = (1<<NUM_C_BITS) -1
	MAX_C       = (1 << NUM_C_BITS) - 1
)

const BucketSize = 4
const MaxLoad = 0.95
const KMaxCount = 10000

// KEY used for highwayhash function, must len(KEY) == 32
//var KEY = []byte{35, 158, 189, 243, 123, 39, 95, 219, 58, 253, 127, 163, 91, 235, 248, 177, 139, 67, 229, 171, 195, 81, 95, 149, 191, 249, 148, 45, 155, 235}
var KEY = []uint64{0xBD4CCC325BEFCA6F, 0xA89A58CE65E641FF, 0xAE093FEF1F84E3E7, 0xFB4297E8C586EE2D}

func CompareAndSwapUint16(addr *uint16, old uint16, new uint16) (swapped bool) {
	//a := (*C.uint16_t)(addr)
	//return C.CompareAndSwapUint16(a, C.uint16_t(old), C.uint16_t(new)) == C.uint16_t(old)
	return atomic.CompareAndSwapPointer((*unsafe.Pointer)(unsafe.Pointer(&addr)), unsafe.Pointer(&old), unsafe.Pointer(&new))
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
	NumItems uint64
	Kmerlen  int
}

var CountItems uint64

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
	cf.NumItems = numBuckets
	cf.Kmerlen = kmerLen
	fmt.Printf("[MakeCuckooFilter]cf items number is: %d\n", cf.NumItems)

	return cf

}

func (cf CuckooFilter) IndexHash(v uint64) uint64 {
	//v = (v >> 37) ^ (v >> 27) ^ v

	return v % cf.NumItems
}

func FingerHash(x []uint64) uint16 {
	//v = (v >> 47) ^ (v >> 33) ^ (v >> 19) ^ (v >> 13) ^ v
	m := uint64(0xc6a4a7935bd1e995)
	var hash uint64
	BITNUM := uint32(13)
	MASK := (uint64(1<<BITNUM) - 1)
	hash = 0x5bd1e995
	for i := len(x) - 1; i >= 0; i-- {
		a := x[i]
		hash += (a & MASK)
		hash *= m
		a >>= BITNUM // NUM_FP_BITS
		hash += (a & MASK)
		hash *= m
		a >>= BITNUM // NUM_FP_BITS
		hash += (a & MASK)
		hash *= m
		a >>= BITNUM // NUM_FP_BITS
		hash += (a & MASK)
		hash *= m
		a >>= BITNUM
		hash += (a & MASK)
		hash *= m
	}

	return uint16(hash & FPMASK)
}

func FingerPrint(data []byte) uint16 {
	hash := metro.Hash64(data, 1335)

	return uint16(hash%FPMASK + 1)
}

func (cf CuckooFilter) AltIndex(index uint64, finger uint16) uint64 {
	fp := make([]byte, 2)
	fp[0] = byte(finger >> 8)
	fp[1] = byte(finger & 255)
	hash := uint64(metro.Hash64(fp, 1337))
	return (index ^ hash) % cf.NumItems
	//index ^= uint64(finger)

	//return index
	//return index % uint64(cf.NumItems)
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

// return oldcount,  added
func (cfi *CFItem) NewAdd(ncfi CFItem) (oldCount int, succ bool) {
	oc := uint16(0)
	a := (*uint16)(cfi)
	if CompareAndSwapUint16(a, oc, uint16(ncfi)) {
		oldCount = 0
		succ = true
	} else {
		oldCount = int(cfi.GetCount())
		succ = false
	}
	return
}

// return oldcount,  added
func (cfi *CFItem) AddCount() (int, bool) {
	for {
		oc := *cfi
		count := oc.GetCount()
		if count < MAX_C {
			nc := oc
			nc.setCount(count + 1)
			a := (*uint16)(cfi)
			if CompareAndSwapUint16(a, uint16(oc), uint16(nc)) {
				return int(count), true
			}
		} else {
			return MAX_C, true
		}
	}
}

func (cf CuckooFilter) Switch(hashIdx uint64, bIdx int, nc CFItem) (CFItem, uint64) {
	a := (*uint16)(&cf.Hash[hashIdx].Bkt[bIdx])
	for {
		oc := *a
		if CompareAndSwapUint16(a, oc, uint16(nc)) {
			fp := CFItem(oc).GetFinger()
			return CFItem(oc), cf.AltIndex(hashIdx, fp)
		}
	}
	//return CFItem(0), uint64(0)
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
					CountItems++
					return CFItem(0), true, 0
				}
			} else {
				if oi.GetCount() > 0 && b.Bkt[i].EqualFP(cfi) {
					oc, _ := b.Bkt[i].AddCount()
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

/*// return if successed added
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
}*/

// return if successed added
func (cf CuckooFilter) Add1(index uint64, fingerprint uint16) (oldcount int, succ bool) {
	i1 := index
	i2 := cf.AltIndex(i1, fingerprint)
	// found same fingerprint and just need add fingerprint count
	for i := 0; i < BucketSize; i++ {
		c1, c2 := cf.Hash[i1].Bkt[i], cf.Hash[i2].Bkt[i]
		if c1 > 0 && c1.GetFinger() == fingerprint {
			return cf.Hash[i1].Bkt[i].AddCount()
		}
		if c2 > 0 && c2.GetFinger() == fingerprint {
			return cf.Hash[i2].Bkt[i].AddCount()
		}
	}
	// if not found same fingerprint, need add bucket
	newCFI := combineCFItem(fingerprint, 1)
	for i := 0; i < BucketSize; i++ {
		if cf.Hash[i1].Bkt[i] == 0 {
			_, ok := cf.Hash[i1].Bkt[i].NewAdd(newCFI)
			if ok {
				//fmt.Printf("[Add1]hashPos: %v, newCFI: %v, Bkt: %v\n", i1, newCFI, cf.Hash[i1])
				return 0, true
			}
		}
		if cf.Hash[i2].Bkt[i] == 0 {
			_, ok := cf.Hash[i2].Bkt[i].NewAdd(newCFI)
			if ok {
				//fmt.Printf("[Add1]hashPos: %v, newCFI: %v, Bkt: %v\n", i2, newCFI, cf.Hash[i2])
				return 0, true
			}
		}
	}

	// need kickout
	var kickCFI CFItem
	var hashPos uint64
	rand.Seed(time.Now().UnixNano())
	idx := rand.Intn(2 * BucketSize)
	if idx < BucketSize {
		//fmt.Printf("[Add1]hashPos: %v, need kickout CFI: %v, Bkt: %v\n", i1, newCFI, cf.Hash[i1])
		kickCFI, hashPos = cf.Switch(i1, idx, newCFI)
		//fmt.Printf("[Add1]hashPos: %v, need kickout CFI: %v, Bkt: %v\n", i1, kickCFI, cf.Hash[i1])
	} else {
		//fmt.Printf("[Add1]hashPos: %v, need kickout CFI: %v, Bkt: %v\n", i2, newCFI, cf.Hash[i2])
		kickCFI, hashPos = cf.Switch(i2, idx-BucketSize, newCFI)
		//fmt.Printf("[Add1]hashPos: %v, need kickout CFI: %v, Bkt: %v\n", i2, kickCFI, cf.Hash[i2])
	}
	oldcount = 0
	for i := 0; i < KMaxCount; i++ {
		//fmt.Printf("[Add1]hashPos: %v, kickCFI: %v, Bkt: %v\n", hashPos, kickCFI, cf.Hash[hashPos])
		for j := 0; j < BucketSize; j++ {
			if cf.Hash[hashPos].Bkt[j] == 0 {
				_, ok := cf.Hash[hashPos].Bkt[j].NewAdd(kickCFI)
				if ok {
					kickCFI = CFItem(0)
					break
				}
			}
		}
		if kickCFI == 0 {
			succ = true
			return
		}
		idx = rand.Intn(BucketSize)
		kickCFI, hashPos = cf.Switch(hashPos, idx, kickCFI)
	}
	succ = false
	return
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
const FnvPrime = uint64(1099511628211)
const FnvOffsetBias = uint64(14695981039346656037)
const MASK31 = (1 << 31) - 1

func HashUint64Arr(arr []uint64, len int) uint64 {
	hash := FnvOffsetBias
	for i := 0; i < len; i++ {
		hash *= ((arr[i] & MASK31) + FnvOffsetBias)
		hash *= FnvPrime
		hash *= ((arr[i] >> 31) + FnvOffsetBias)
		hash *= FnvPrime
	}
	return hash
}

var key = []uint64{uint64(0x32b87ef98934abcf), uint64(0x4f4f28a7931afc1b), uint64(0xf56a8aa814946e08), uint64(0xefebe623f41cd45d)}

func Transform(kb []byte) []byte {
	tkb := make([]byte, len(kb))
	for i := 0; i < len(kb); i++ {
		tkb[i] = (kb[i] + 7) * 23
	}
	return tkb
}

// getIndicesAndFingerprint returns the 2 bucket indices and fingerprint to be used
func (cf CuckooFilter) GetIndicesAndFingerprint(data []byte) (uint64, uint64, uint16) {
	hash := metro.Hash64(data, 1337)
	f := FingerPrint(data)
	i1 := uint64(hash) % cf.NumItems
	i2 := cf.AltIndex(i1, f)
	return i1, i2, f
}

func (cf CuckooFilter) insert(cfi CFItem, i uint64) (int, bool) {
	_, ok, count := cf.Hash[i].AddBucket(cfi, false)
	return count, ok
}

func randi(i1, i2 uint64) uint64 {
	if rand.Intn(2) == 0 {
		return i1
	}
	return i2
}

func (cf CuckooFilter) reinsert(cfi CFItem, i uint64) (int, bool) {
	for k := 0; k < KMaxCount; k++ {
		j := rand.Intn(BucketSize)
		cfi, i = cf.Switch(i, j, cfi)

		// look in the alternate location for that random element
		count, ok := cf.insert(cfi, i)
		if ok {
			return count, ok
		}
	}
	return 0, false
}

func (cf CuckooFilter) Insert(kb []byte) (int, bool) {
	i1, i2, fp := cf.GetIndicesAndFingerprint(kb)
	cfi := combineCFItem(fp, 1)
	count, ok := cf.insert(cfi, i1)
	if ok {
		return count, ok
	}
	count, ok = cf.insert(cfi, i2)
	if ok {
		return count, ok
	}

	return cf.reinsert(cfi, randi(i1, i2))
}

// return last count of kmer fingerprint and have been successed inserted
/*func (cf CuckooFilter) Insert(kb []byte) (int, bool) {
	//tkb := Transform(kb)
	hash := metro.Hash64(kb, 1337)
	fingerprint := FingerPrint(kb)
	//hash := highwayhash.SumInput64Arr64(kb, key)

	//hk := sha1.Sum(kb)
	//v := hk2uint64(hk)
	//hash := HashUint64Arr(kb, len(kb))
	//fmt.Printf("%v\t", v)
	index := cf.IndexHash(hash)
	//fmt.Printf("[cf.Insert]index: %v\tfinger: %v\n", index, fingerprint)
	//fmt.Printf(" sizeof cuckoofilter.Hash[0] : %d\n", unsafe.Sizeof(cf.Hash[0]))

	return cf.Add1(index, fingerprint)
}*/

func (cf CuckooFilter) Lookup(kb []byte) bool {
	//tkb := Transform(kb)
	//fingerprint := FingerHash(kb)
	hash := metro.Hash64(kb, 1337)
	fingerprint := FingerPrint(kb)
	//hk := sha1.Sum(kb)
	//v := hk2uint64(hk)
	//hash := HashUint64Arr(kb, len(kb))
	//hash := highwayhash.SumInput64Arr64(kb, key)
	index := cf.IndexHash(hash)

	if cf.Hash[index].Contain(fingerprint) {
		return true
	} else {
		index = cf.AltIndex(index, fingerprint)
		return cf.Hash[index].Contain(fingerprint)
	}
}

func (cf CuckooFilter) GetCount(kb []byte) uint16 {
	//fingerprint := FingerHash(kb)
	//tkb := Transform(kb)
	hash := metro.Hash64(kb, 1337)
	fingerprint := FingerPrint(kb)
	//hk := sha1.Sum(kb)
	//v := hk2uint64(hk)
	//v := highwayhash.SumInput64Arr64(kb, KEY)
	//hash := HashUint64Arr(kb, len(kb))
	//hash := highwayhash.SumInput64Arr64(kb, key)
	index := cf.IndexHash(hash)
	for _, item := range cf.Hash[index].Bkt {
		// fmt.Printf("index: %v, finger: %v\n", index, item.GetFinger())
		if item > 0 && item.GetFinger() == fingerprint {
			return item.GetCount()
		}
	}
	// if not return , find another position
	index = cf.AltIndex(index, fingerprint)
	for _, item := range cf.Hash[index].Bkt {
		// fmt.Printf("index: %v, finger: %v\n", index, item.GetFinger())
		if item > 0 && item.GetFinger() == fingerprint {
			return item.GetCount()
		}
	}

	// not found in the CuckooFilter
	panic("not found in the CuckooFilter")
}

/* allow function return zero version if kb not found in the CuckooFilter */
func (cf CuckooFilter) GetCountAllowZero(kb []byte) uint16 {
	//hk := sha1.Sum(kb)
	//v := hk2uint64(hk)
	//v := highwayhash.SumInput64Arr64(kb, KEY)
	//hash := HashUint64Arr(kb, len(kb))
	// fmt.Printf("[GetCountAllowZero] len(cf.Hash): %d\n", len(cf.Hash))
	// fmt.Printf("[GetCountAllowZero] cf.Hash[0]: %v\n", cf.Hash[0])
	//hash := highwayhash.SumInput64Arr64(kb, key)
	//fingerprint := FingerHash1(hash)
	//tkb := Transform(kb)
	hash := metro.Hash64(kb, 1337)
	fingerprint := FingerPrint(kb)
	index := cf.IndexHash(hash)
	//fmt.Printf("[cf.GetCountAllowZero]kb: %v\n\tindex: %v, finger: %v\n", kb, index, fingerprint)
	for _, item := range cf.Hash[index].Bkt {
		if item > 0 && item.GetFinger() == fingerprint {
			//fmt.Printf("[cf.GetCountAllowZero]i: %v, finger: %v\n", index, item.GetFinger())
			return item.GetCount()
		}
	}
	// if not return , find another position
	index = cf.AltIndex(index, fingerprint)
	for _, item := range cf.Hash[index].Bkt {
		//fmt.Printf("[cf.GetCountAllowZero]i: %v, finger: %v\n", index, item.GetFinger())
		if item > 0 && item.GetFinger() == fingerprint {
			return item.GetCount()
		}
	}

	// not found in the CuckooFilter, return zero
	// panic("not found in the CuckooFilter")
	return 0
}

func (cf CuckooFilter) GetStat() {
	var ca [MAX_C + 1]int
	for _, b := range cf.Hash {
		for _, e := range b.Bkt {
			count := e.GetCount()
			ca[count]++
		}
	}
	var total int
	for i := 1; i < MAX_C+1; i++ {
		total += ca[i]
	}
	fmt.Printf("count statisticas : %v\n", ca)
	fmt.Printf("cuckoofilter numItems : %d, CountItems: %d, load: %f\n", cf.NumItems, total, float64(total)/float64(cf.NumItems*BucketSize))
}

func (cf CuckooFilter) MmapWriter(cfmmapfn string) error {
	cfmmapfp, err := os.Create(cfmmapfn)
	if err != nil {
		return err
	}
	defer cfmmapfp.Close()

	cbrofp := cbrotli.NewWriter(cfmmapfp, cbrotli.WriterOptions{Quality: 1})
	defer cbrofp.Close()
	buffp := bufio.NewWriterSize(cbrofp, 1<<25) // 1<<24 == 2**24

	// write cuckoofilter to the memory map file
	enc := gob.NewEncoder(buffp)
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
	if err := buffp.Flush(); err != nil {
		log.Fatalf("[MmapWriter] failed to flush file: %s, err: %v\n", cfmmapfn, err)
	}

	if err := cbrofp.Flush(); err != nil {
		log.Fatalf("[MmapWriter] failed to flush file: %s, err: %v\n", cfmmapfn, err)
	}

	return nil
}

func (cf CuckooFilter) HashWriter(cffn string) (err error) {
	cffp, err := os.Create(cffn)
	if err != nil {
		return err
	}
	defer cffp.Close()

	cbrofp := cbrotli.NewWriter(cffp, cbrotli.WriterOptions{Quality: 1})
	defer cbrofp.Close()
	buffp := bufio.NewWriterSize(cbrofp, 1<<25) // 1<<25 == 2**25
	// write cuckoofilter to the memory map file
	err = binary.Write(buffp, binary.LittleEndian, cf.Hash)
	if err != nil {
		return err
	}
	/*for _, b := range cf.Hash {
		err = binary.Write(buffp, binary.LittleEndian, b)
		if err != nil {
			return err
		}
	}*/

	/*if n != cf.numItems*unsafe.Sizeof(bucket) {
		log.Fatal("write byte number is not equal cf size")
	}*/
	if err := buffp.Flush(); err != nil {
		log.Fatalf("[MmapWriter] failed to flush file: %s, err: %v\n", cffn, err)
	}

	if err := cbrofp.Flush(); err != nil {
		log.Fatalf("[MmapWriter] failed to flush file: %s, err: %v\n", cffn, err)
	}

	return nil
}

func (cf CuckooFilter) WriteCuckooFilterInfo(cfinfofn string) error {
	cfinfofp, err := os.Create(cfinfofn)
	if err != nil {
		return err
	}
	defer cfinfofp.Close()

	_, err = cfinfofp.WriteString(fmt.Sprintf("NumItems\t%d\n", cf.NumItems))
	if err != nil {
		return err
	}
	_, err = cfinfofp.WriteString(fmt.Sprintf("Kmerlen\t%d\n", cf.Kmerlen))
	if err != nil {
		return err
	}

	return nil
}

func RecoverCuckooFilterInfo(cfinfofn string) (CuckooFilter, error) {
	var cfinfofp *os.File
	var err error
	var cf CuckooFilter
	cfinfofp, err = os.Open(cfinfofn)
	if err != nil {
		return cf, err
	}
	defer cfinfofp.Close()
	cfinfobuf := bufio.NewReaderSize(cfinfofp, 1<<25)
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
	_, err = fmt.Sscanf(line, "NumItems\t%d\n", &cf.NumItems)
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

	if cf.NumItems == 0 || cf.NumItems%2 != 0 {
		log.Fatalf("[RecoverCuckooFilterInfo] cf.NumItems: %v, error\n", cf.NumItems)
	}

	return cf, err
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
	brfp := cbrotli.NewReader(cfmmapfp)
	defer brfp.Close()
	buffp := bufio.NewReaderSize(brfp, 1<<25)

	dec := gob.NewDecoder(buffp)
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

func (cf CuckooFilter) HashReader(cffn string) error {
	// if cf.numItems <= 0 {
	// 	log.Fatal("CuckooFilter number is <=0, please check")
	// }
	cffp, err := os.Open(cffn)
	if err != nil {
		return err
	}
	defer cffp.Close()
	brfp := cbrotli.NewReader(cffp)
	defer brfp.Close()
	buffp := bufio.NewReaderSize(brfp, 1<<25)

	// read hash array
	fmt.Printf("[HashReader]cf.NumItems: %v, len(cf.Hash): %v, cf.Kmerlen: %v\n", cf.NumItems, len(cf.Hash), cf.Kmerlen)
	err = binary.Read(buffp, binary.LittleEndian, cf.Hash)
	if err != nil {
		return err
	}
	/*for i := 0; i < int(cf.NumItems); i++ {
		err := binary.Read(buffp, binary.LittleEndian, &cf.Hash[i])
		if err != nil {
			return err
		}
	}*/

	return nil
}
