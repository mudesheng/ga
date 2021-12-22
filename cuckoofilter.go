package main

// #include<stdint.h>
/*
uint16_t CompareAndSwapUint16(uint16_t *addr, uint16_t old, uint16_t new)
{
    return __sync_val_compare_and_swap(addr, old, new);
}

uint8_t CompareAndSwapUint8(uint8_t *addr, uint8_t old, uint8_t new)
{
    return __sync_val_compare_and_swap(addr, old, new);
}
*/
import "C"

import (
	"fmt"
	"math/bits"
	"math/rand"

	//"sync/atomic"

	//"log"

	//"strings"

	"github.com/cespare/xxhash"
	//"github.com/klauspost/compress/zstd/internal/xxhash"
	//"github.com/mudesheng/highwayhash"
	// "-syscall"
)

const (
	//NumFpBits number bits for Fingerprint
	NumFpBits = 13 // number of fingerprint  bits occpied
	//NumCBits number bits for freq Count
	NumCBits = 3 // number of count equal sizeof(uint16)*8 - NumFpBits
	//FpMask mask other info only set fingerprint = (1<<NumFpBits) -1
	MaxC          = (1 << NumCBits) - 1
	CMask         = MaxC
	FpMask        = (1 << NumFpBits) - 1
	MaxCFKmerFreq = 3
	FPDBGKmerLen  = 16
	//FpMask = 0xFFF8
	//CMask set count bits field = (1<<NumCBits) -1
	//MaxC max allowded for freq Count
)

const BucketSize = 4
const MaxLoad = 0.95
const KMaxCount = 10000

// KEY used for highwayhash function, must len(KEY) == 32
//var KEY = []byte{35, 158, 189, 243, 123, 39, 95, 219, 58, 253, 127, 163, 91, 235, 248, 177, 139, 67, 229, 171, 195, 81, 95, 149, 191, 249, 148, 45, 155, 235}
//var KEY = []uint64{0xBD4CCC325BEFCA6F, 0xA89A58CE65E641FF, 0xAE093FEF1F84E3E7, 0xFB4297E8C586EE2D}
func CompareAndSwapUint8(addr *uint8, old uint8, new uint8) (swapped bool) {
	//fmt.Printf("[CompareAndSwapUint16]*addr:%d old:%d new:%d\n", *addr, old, new)
	a := (*C.uint8_t)(addr)
	return C.CompareAndSwapUint8(a, C.uint8_t(old), C.uint8_t(new)) == C.uint8_t(old)
}

func CompareAndSwapUint16(addr *uint16, old uint16, new uint16) (swapped bool) {
	//fmt.Printf("[CompareAndSwapUint16]*addr:%d old:%d new:%d\n", *addr, old, new)
	a := (*C.uint16_t)(addr)
	return C.CompareAndSwapUint16(a, C.uint16_t(old), C.uint16_t(new)) == C.uint16_t(old)
	//av := (*atomic.Uint16)(unsafe.Pointer(addr))
	//return av.CompareAndSwap(old, new)
	//return atomic.CompareAndSwapPointer((*unsafe.Pointer)(unsafe.Pointer(&addr)), unsafe.Pointer(&old), unsafe.Pointer(&new))
}

var masks [65]uint64

func init() {
	for i := uint64(0); i <= 64; i++ {
		masks[i] = (1 << i) - 1
	}
}

//type CFItem uint16

/* type cfitem struct {
  fingerprint uint16:NumFpBits
  count uint16:NumCBits
} */

type Bucket [BucketSize]uint16

type DBGKmer struct {
	ID     uint32
	Pos    uint32
	FP     uint16
	Strand bool
	//Freq uint16
}

type CuckooFilter struct {
	Hash      []Bucket
	Count     uint
	BucketPow uint
	Kmerlen   int
}

type CuckooFilterDBGKmer struct {
	Hash      [][BucketSize]DBGKmer
	Count     uint
	BucketPow uint
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
	cf.Count = 0
	cf.Kmerlen = kmerLen
	cf.BucketPow = uint(bits.TrailingZeros(uint(numBuckets)))
	fmt.Printf("[MakeCuckooFilter]numBuckets:%d cf BucketPow:%d\n", len(cf.Hash), cf.BucketPow)
	return cf

}

func MakeCuckooFilterDBGKmer(maxNumKeys uint64) (cf CuckooFilterDBGKmer) {
	numBuckets := upperpower2(maxNumKeys) / BucketSize
	/*frac := float64(maxNumKeys) / numBuckets / BucketSize
	if frac > MaxLoad {
		numBuckets <<= 1
	} */
	cf.Hash = make([][BucketSize]DBGKmer, numBuckets)
	cf.Count = 0
	cf.BucketPow = uint(bits.TrailingZeros(uint(numBuckets)))
	fmt.Printf("[MakeCuckooFilter]numBuckets:%d cf BucketPow:%d\n", len(cf.Hash), cf.BucketPow)
	return cf
}

/*func FingerHash(x []uint64) uint16 {
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
		a >>= BITNUM // NumFpBits
		hash += (a & MASK)
		hash *= m
		a >>= BITNUM // NumFpBits
		hash += (a & MASK)
		hash *= m
		a >>= BITNUM // NumFpBits
		hash += (a & MASK)
		hash *= m
		a >>= BITNUM
		hash += (a & MASK)
		hash *= m
	}

	return uint16(hash & FpMask)
}*/

func combineFpC(fp uint16, count uint16) (fc uint16) {
	if count > MaxC {
		panic("count bigger than CFItem allowed")
	}
	//fmt.Printf("fp: %d, count: %d, cfi: %d\n", fp, count, uint16(cfi))
	return (fp << NumCBits) | count
}

func GetCount(fc uint16) uint16 {
	return fc & CMask
}

func GetFinger(fc uint16) uint16 {
	return fc >> NumCBits
}

func EqualFP(fc1, fc2 uint16) bool {
	if (fc1 >> NumCBits) == (fc2 >> NumCBits) {
		return true
	} else {
		return false
	}
}

func (cf CuckooFilter) Switch(hashIdx uint64, bIdx int, nc uint16) uint16 {
	var a *uint16
	a = &cf.Hash[hashIdx][bIdx]
	for {
		oc := *a
		if CompareAndSwapUint16(a, oc, nc) {
			return oc
		}
	}
	//return CFItem(0), uint64(0)
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
//const FnvPrime = uint64(1099511628211)
//const FnvOffsetBias = uint64(14695981039346656037)
//const MASK31 = (1 << 31) - 1

/*func HashUint64Arr(arr []uint64, len int) uint64 {
	hash := FnvOffsetBias
	for i := 0; i < len; i++ {
		hash *= ((arr[i] & MASK31) + FnvOffsetBias)
		hash *= FnvPrime
		hash *= ((arr[i] >> 31) + FnvOffsetBias)
		hash *= FnvPrime
	}
	return hash
}*/

//var key = []uint64{uint64(0x32b87ef98934abcf), uint64(0x4f4f28a7931afc1b), uint64(0xf56a8aa814946e08), uint64(0xefebe623f41cd45d)}

/*func Transform(kb []byte) []byte {
	tkb := make([]byte, len(kb))
	for i := 0; i < len(kb); i++ {
		tkb[i] = (kb[i] + 7) * 23
	}
	return tkb
}*/

func (b *Bucket) insert(fc uint16) bool {
	for i, tfc := range b {
		if tfc == 0 {
			var a *uint16
			a = &(b[i])
			for {
				oc := *a
				if oc != 0 {
					break
				}
				if CompareAndSwapUint16(a, oc, fc) {
					return true
				}
				//fmt.Printf("[insert] oc:%d fc:%d\n", oc, fc)
			}
		}
	}
	return false
}

func (cf *CuckooFilter) reinsert(fp uint16, i uint64) bool {
	//fp = (fp << NumCBits) | 1
	for k := 0; k < KMaxCount; k++ {
		j := rand.Intn(BucketSize)
		fp = cf.Switch(i, j, fp)
		// look in the alternate location for that random element
		i = getAltIndex(fp, i)
		if cf.Hash[i].insert(fp) {
			return true
		}
	}
	return false
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

func (b *Bucket) getFingerprintIndex(fp uint16) (int, uint16) {
	for i, tfp := range b {
		if (tfp >> NumCBits) == fp {
			return i, tfp & CMask
		}
	}
	return -1, 0
}

func getDBGKmer(bc [BucketSize]DBGKmer, fp uint16, da []DBGKmer) []DBGKmer {
	for _, d := range bc {
		if d.ID == 0 {
			break
		}
		if d.FP == fp {
			da = append(da, d)
		}
	}
	return da
}

func getAltIndex(fp uint16, i uint64) uint64 {
	return i ^ uint64(fp)
}

func getFingerprint(hash uint64) uint16 {
	// Use least significant bits for fingerprint.
	return uint16(hash & FpMask)
}

func getFingerprintDBGKmer(hash uint64) uint16 {
	// Use least significant bits for fingerprint.
	return uint16(hash & masks[FPDBGKmerLen])
}

func getIndexAndFingerprint(hash uint64, bucketPow uint) (uint64, uint16) {
	fp := getFingerprint(hash)
	// Use most significant bits for deriving index.
	//i1 := (hash >> NumFpBits) & masks[bucketPow]
	//i1 := (((hash >> NumFpBits) & masks[bucketPow]) ^ (hash >> (NumFpBits + bucketPow))) & masks[bucketPow]
	m := 64 - NumFpBits - bucketPow
	fmt.Printf("[getIndexAndFingerprint]bucketPow:%d m:%d\n", bucketPow, m)
	i1 := (hash >> (m + NumFpBits)) ^ (((hash >> NumFpBits) & masks[m]) << ((bucketPow - m) >> 1))
	return i1, fp
}

func getIndexAndFingerprintDBGKmer(hash uint64, bucketPow uint) (uint64, uint16) {
	fp := getFingerprintDBGKmer(hash)
	// Use most significant bits for deriving index.
	//i1 := uint64(hash>>FPDBGKmerLen) & masks[bucketPow]
	//i1 := (((hash >> FPDBGKmerLen) & masks[bucketPow]) ^ (hash >> (FPDBGKmerLen + bucketPow))) & masks[bucketPow]
	m := 64 - FPDBGKmerLen - bucketPow
	i1 := (hash >> (m + FPDBGKmerLen)) ^ (((hash >> FPDBGKmerLen) & masks[m]) << ((bucketPow - m) >> 1))
	return i1, fp
}

func (cf *CuckooFilter) Lookup(kb []byte) (uint16, bool) {
	hash := xxhash.Sum64(kb)
	i1, fp := getIndexAndFingerprint(hash, cf.BucketPow)
	j, c := cf.Hash[i1].getFingerprintIndex(fp)
	if j >= 0 {
		return c, true
	}
	i2 := getAltIndex(fp, i1)
	j, c = cf.Hash[i2].getFingerprintIndex(fp)
	if j >= 0 {
		return c, true
	}
	return 0, false
}

func (cf *CuckooFilterDBGKmer) Lookup(kb []byte, da []DBGKmer) []DBGKmer {
	hash := xxhash.Sum64(kb)
	da = da[:0]
	i1, fp := getIndexAndFingerprintDBGKmer(hash, cf.BucketPow)
	da = getDBGKmer(cf.Hash[i1], fp, da)
	i2 := getAltIndex(fp, i1)
	if i2 != i1 {
		da = getDBGKmer(cf.Hash[i2], fp, da)
	}
	//fmt.Printf("[cf.Lookup]i1:%d i2:%d da:%v\n", i1, i2, da)
	return da
}

var addFalseCount int

func (cf *CuckooFilter) Add1(i uint64, j int, fp uint16) bool {
	var a *uint16
	a = &cf.Hash[i][j]
	for {
		oc := *a
		if (oc >> NumCBits) == fp {
			if (oc & CMask) < MaxCFKmerFreq {
				fc := (fp << NumCBits) | ((oc & CMask) + 1)
				if CompareAndSwapUint16(a, oc, fc) {
					return true
				}
			} else {
				return true
			}
		} else {
			return false
		}
	}
	return false
}

func AddDBGKmer(b *[BucketSize]DBGKmer, d DBGKmer) (ok bool) {
	for i := 0; i < BucketSize; i++ {
		if b[i].ID == 0 {
			b[i] = d
			ok = true
			break
		}
	}
	return
}

func (cf *CuckooFilter) Insert(kb []byte) (uint16, bool) {
	hash := xxhash.Sum64(kb)
	fmt.Printf("[cf.Insert]hash:%X\n", hash)
	i1, fp := getIndexAndFingerprint(hash, cf.BucketPow)
	fmt.Printf("[cf.Insert]i1:%X fp:%X\n", i1, fp)
	// lookup
	j, c := cf.Hash[i1].getFingerprintIndex(fp)
	if j >= 0 {
		if c < MaxCFKmerFreq {
			if !cf.Add1(i1, j, fp) {
				addFalseCount++
			}
		}
		return c, true
	}
	i2 := getAltIndex(fp, i1)
	j, c = cf.Hash[i2].getFingerprintIndex(fp)
	if j > -1 {
		if c < MaxCFKmerFreq {
			if !cf.Add1(i2, j, fp) {
				addFalseCount++
			}
		}
		return c, true
	}

	fc := (fp << NumCBits) | 1
	if cf.Hash[i1].insert(fc) {
		cf.Count++
		return 0, true
	}
	//i2 := getAltIndex(fp, i1)
	if cf.Hash[i2].insert(fc) {
		cf.Count++
		return 0, true
	}

	// select index
	i := i1
	if ((hash >> (NumFpBits + cf.BucketPow)) & 1) == 1 {
		i = i2
	}
	if cf.reinsert(fc, i) {
		cf.Count++
		return 0, true
	}
	return 0, false
}

func (cf *CuckooFilterDBGKmer) reinsert(i uint64, dk DBGKmer) bool {
	//fp = (fp << NumCBits) | 1
	for k := 0; k < KMaxCount; k++ {
		j := rand.Intn(BucketSize)
		//fp = cf.Switch(i, j, fp)
		cf.Hash[i][j], dk = dk, cf.Hash[i][j]
		// look in the alternate location for that random element
		i = getAltIndex(dk.FP, i)
		for x := 0; x < BucketSize; x++ {
			if cf.Hash[i][x].ID == 0 {
				cf.Hash[i][x] = dk
				return true
			}
		}
	}
	return false
}

func (cf *CuckooFilterDBGKmer) Insert(kb []byte, dk DBGKmer) bool {
	//fmt.Printf("[cf.Insert]dk:%v\n", dk)
	hash := xxhash.Sum64(kb)
	i1, fp := getIndexAndFingerprintDBGKmer(hash, cf.BucketPow)
	dk.FP = fp
	if AddDBGKmer(&cf.Hash[i1], dk) {
		//cf.Count++
		return true
	}
	i2 := getAltIndex(fp, i1)
	if AddDBGKmer(&cf.Hash[i2], dk) {
		//cf.Count++
		return true
	}

	// select index
	i := i1
	if (hash>>(FPDBGKmerLen+cf.BucketPow))&1 == 1 {
		i = i2
	}

	if cf.reinsert(i, dk) {
		//cf.Count++
		return true
	}
	return false
}

// allow function return zero version if kb not found in the CuckooFilter
/*func (cf *CuckooFilter) GetCountAllowZero(kb []byte) uint16 {
	//hk := sha1.Sum(kb)
	//v := hk2uint64(hk)
	//v := highwayhash.SumInput64Arr64(kb, KEY)
	//hash := HashUint64Arr(kb, len(kb))
	// fmt.Printf("[GetCountAllowZero] len(cf.Hash): %d\n", len(cf.Hash))
	// fmt.Printf("[GetCountAllowZero] cf.Hash[0]: %v\n", cf.Hash[0])
	//hash := highwayhash.SumInput64Arr64(kb, key)
	//fingerprint := FingerHash1(hash)
	//tkb := Transform(kb)
	//hash := metro.Hash64(kb, 1337)
	hash := xxhash.Sum64(kb)
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
}*/

func (cf *CuckooFilter) GetStat() {
	var ca [MaxC + 1]int
	for _, b := range cf.Hash {
		for _, e := range b {
			count := e & CMask
			ca[count]++
		}
	}
	var total int
	for i := 1; i < MaxC+1; i++ {
		total += ca[i]
	}
	fmt.Printf("[GetStat]addFalseCount:%d cf contained Count:%d count statisticas:%v\n", addFalseCount, cf.Count, ca)
	fmt.Printf("[GetStat]cuckoofilter numItems:%d CountItems:%d load:%f\n", len(cf.Hash), total, float64(total)/float64(len(cf.Hash)*BucketSize))
}

/*func (cf CuckooFilter) MmapWriter(cfmmapfn string) error {
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

	if err := buffp.Flush(); err != nil {
		log.Fatalf("[MmapWriter] failed to flush file: %s, err: %v\n", cfmmapfn, err)
	}

	if err := cbrofp.Flush(); err != nil {
		log.Fatalf("[MmapWriter] failed to flush file: %s, err: %v\n", cfmmapfn, err)
	}

	return nil
}*/

/*func (cf CuckooFilter) HashWriter(cffn string) (err error) {
	cffp, err := os.Create(cffn)
	if err != nil {
		return err
	}
	defer cffp.Close()

	//cbrofp := cbrotli.NewWriter(cffp, cbrotli.WriterOptions{Quality: 1})
	//defer cbrofp.Close()
	buffp := bufio.NewWriterSize(cbrofp, 1<<25) // 1<<25 == 2**25
	// write cuckoofilter to the memory map file
	err = binary.Write(buffp, binary.LittleEndian, cf.Hash)
	if err != nil {
		return err
	}
	if err := buffp.Flush(); err != nil {
		log.Fatalf("[MmapWriter] failed to flush file: %s, err: %v\n", cffn, err)
	}

	if err := cbrofp.Flush(); err != nil {
		log.Fatalf("[MmapWriter] failed to flush file: %s, err: %v\n", cffn, err)
	}

	return nil
}*/

/*func (cf CuckooFilter) WriteCuckooFilterInfo(cfinfofn string) error {
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
}*/

/*func RecoverCuckooFilterInfo(cfinfofn string) (CuckooFilter, error) {
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
}*/

/*func MmapReader(cfmmapfn string) (cf CuckooFilter, err error) {
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
}*/

/*func (cf CuckooFilter) HashReader(cffn string) error {
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
	return nil
}*/
