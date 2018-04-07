package constructdbg

import (
	"crypto/sha1"
	"fmt"
	"math/rand"
	"reflect"

	//"log"

	//"strings"

	// "syscall"
	"unsafe"
)

const (
	NUM_FP_BITS = 15                 // number of fingerprint  bits occpied
	NUM_C_BITS  = 1                  // number of count equal sizeof(uint16)*8 - NUM_FP_BITS
	FPMASK      = 1<<NUM_FP_BITS - 1 // mask other info only set fingerprint = (1<<NUM_FP_BITS) -1
	CMASK       = 1<<NUM_C_BITS - 1  // set count bits field = (1<<NUM_C_BITS) -1
	MAX_C       = 1<<NUM_C_BITS - 1
)

const BucketSize = 4
const MaxLoad = 0.95
const KMaxCount = 512

//const MAXFREQ = math.MaxUint16

//func CompareAndSwapUint16(addr *uint16, old uint16, new uint16) (swapped bool)

type CFItem uint16

type DBGKmer struct {
	Item   CFItem
	Strand bool
	//Freq uint16
	ID  DBG_MAX_INT
	Pos uint32
}

/* type cfitem struct {
  fingerprint uint16:NUM_FP_BITS
  count uint16:NUM_C_BITS
} */

type Bucket struct {
	Bkt [BucketSize]DBGKmer
}

type CuckooFilter struct {
	Hash     []Bucket
	NumItems uint64
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
	}*/

	cf.Hash = make([]Bucket, numBuckets)
	cf.NumItems = numBuckets * BucketSize
	cf.Kmerlen = kmerLen

	return cf

}

func (cf CuckooFilter) IndexHash(v uint64) uint64 {
	v >>= NUM_FP_BITS

	return v % uint64(len(cf.Hash))
}

func FingerHash(v uint64) uint64 {
	v &= FPMASK
	return v
}

func (cf CuckooFilter) AltIndex(index uint64, finger uint64) uint64 {
	index ^= finger

	return index % uint64(len(cf.Hash))
}

/*func combineCFItem(fp uint64, count uint64) (cfi CFItem) {
	if count > MAX_C {
		panic("count bigger than CFItem allowed")
	}
	cfi = CFItem(fp)
	cfi <<= NUM_C_BITS
	cfi |= CFItem(count)
	//fmt.Printf("fp: %d, count: %d, cfi: %d\n", fp, count, uint16(cfi))
	return cfi
}*/

func (dbgK DBGKmer) GetCount() uint16 {
	return uint16(dbgK.Item) & CMASK
}

func (dbgK *DBGKmer) setCFItem(finger, count uint16) {
	if count > MAX_C {
		panic("[setCFItem] count bigger than CFItem allowed")
	}
	nc := finger << NUM_C_BITS
	nc |= count

	dbgK.Item = CFItem(nc)
}

func (dbgK DBGKmer) GetFinger() uint16 {
	return uint16(dbgK.Item) >> NUM_C_BITS
}

/*func (dbgK DBGKmer) EqualFP(sec DBGKmer) bool {
	if dbgK.GetFinger() == sec.GetFinger() {
		return true
	} else {
		return false
	}
} */

// return oldcount, oldfinger, added
/*func (dbgK *DBGKmer) AddFreq(f uint16) bool {
	for {
		of := dbgK.Freq
		if of < math.MaxUint16 {
			nc := of + f
			a := (*uint16)(&dbgK.Freq)
			if cuckoofilter.CompareAndSwapUint16(a, uint16(of), uint16(nc)) {
				return true
			}
		} else {
			return true
		}
	}

	return false
} */

/*func (b Bucket) AddFreq(dbgK DBGKmer, c int) bool {
	for i, d := range b.Bkt {
		if d.Item == dbgK.Item && d.ID == dbgK.ID && d.Pos == dbgK.Pos {
			b.Bkt[i].AddFreq(uint16(c))
			return true
		}
	}
	return false
}*/

func (cf CuckooFilter) Contain(v uint64) (arr []DBGKmer) {
	index := cf.IndexHash(v)
	fingerprint := FingerHash(v)
	for _, dk := range cf.Hash[index].Bkt {
		if dk.GetCount() > 0 && dk.GetFinger() == uint16(fingerprint) {
			arr = append(arr, dk)
		}
	}
	index2 := cf.AltIndex(index, fingerprint)
	for _, dk := range cf.Hash[index2].Bkt {
		if dk.GetCount() > 0 && dk.GetFinger() == uint16(fingerprint) {
			arr = append(arr, dk)
		}
	}
	return arr
}

func (b *Bucket) AddBucket(dbgK DBGKmer, kickout bool) (old DBGKmer, suc bool) {
	//fmt.Printf("[AddBucket]b.Bkt: %v\n\tdbgK: %v\n", b.Bkt, dbgK)
	for i := 0; i < BucketSize; i++ {
		//fmt.Printf("i: %d\n", i)
		if b.Bkt[i].GetCount() == 0 {
			b.Bkt[i] = dbgK
			suc = true
			return old, suc
		}
	}
	//else {
	//	if b.Bkt[i].EqualFP(dbgK) {
	//		log.Fatalf("[AddBucket] found conflicting DBG Kmer, please increase size cuckoofilter of DBG Sample\n")
	//	}
	//}

	//fmt.Printf("kikcout: %t", kickout)
	if kickout {
		ci := rand.Uint32() % BucketSize
		old = b.Bkt[ci]
		b.Bkt[ci] = dbgK
		suc = true
		return old, suc
	}
	return old, suc
}

// return last count of kmer and if successed added
func (cf CuckooFilter) Add(index uint64, dbgK DBGKmer) bool {
	ci := index
	gK := dbgK
	_, added := cf.Hash[ci].AddBucket(gK, false)
	if added {
		return true
	}
	// choose alter index
	ci = cf.AltIndex(ci, uint64(gK.GetFinger()))
	for count := 0; count < KMaxCount; count++ {
		//b := &cf.Hash[ci]
		old, _ := cf.Hash[ci].AddBucket(gK, true)
		if old.GetCount() == 0 {
			return true
		}
		gK = old
		//fmt.Printf("[Add]cycle : %d\n\tb: %v\n", count, b)

		ci = cf.AltIndex(ci, uint64(gK.GetFinger()))
	}
	return false
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
func (cf CuckooFilter) Insert(kb []byte, id DBG_MAX_INT, pos uint32, strand bool) bool {
	hk := sha1.Sum(kb)
	v := hk2uint64(hk)
	//fmt.Printf("%v\t", v)
	index := cf.IndexHash(v)
	fingerprint := FingerHash(v)
	var dbgK DBGKmer
	dbgK.setCFItem(uint16(fingerprint), 1)
	dbgK.ID = id
	dbgK.Pos = pos
	dbgK.Strand = strand
	//fmt.Printf("[cf.Insert]fingerprint: %v\tindex: %v\tdbgK.%v\n", fingerprint, index, dbgK)
	//fmt.Printf(" sizeof cuckoofilter.Hash[0] : %d\n", unsafe.Sizeof(cf.Hash[0]))

	return cf.Add(index, dbgK)
}

func (cf CuckooFilter) Lookup(kb []byte, edgesArr []DBGEdge) (dbgK DBGKmer, count int) {
	hk := sha1.Sum(kb)
	v := hk2uint64(hk)
	da := cf.Contain(v)
	//fmt.Printf("[cf.Lookup] da: %v\n", da)
	if len(da) == 0 {
		return
	} else if len(da) == 1 {
		d := da[0]
		eb := edgesArr[d.ID].Utg.Ks[d.Pos : d.Pos+uint32(cf.Kmerlen)]
		//fmt.Printf("[cf.Lookup]\n\tkb: %v\n\teb: %v\n", kb, eb)
		if d.Strand == MINUS {
			eb = GetReverseCompByteArr(eb)
		}
		if reflect.DeepEqual(kb, eb) {
			dbgK = d
			count++
		} else {
			fmt.Printf("[cf.Lookup] found cf Item, but seq(%v) not same as read(%v)\n", eb, kb)
		}
	} else { // len(da) > 1
		// count := 0
		//fmt.Printf("[cf.Lookup]\n\tkb: %v\n", kb)
		for _, d := range da {
			//fmt.Printf("[Lookup] fingerprint: %v\td: %v\n\tedgesArr[%v]: %v\n", d.GetFinger(), d, d.ID, edgesArr[d.ID])
			eb := edgesArr[d.ID].Utg.Ks[d.Pos : d.Pos+uint32(cf.Kmerlen)]
			if d.Strand == MINUS {
				eb = GetReverseCompByteArr(eb)
			}
			//fmt.Printf("\teb: %v\n", eb)
			if reflect.DeepEqual(kb, eb) {
				dbgK = d
				count++
			}
		}
	}

	// add the frequency of kmer
	/*if dbgK.ID > 0 {
		if cf.Hash[index].AddFreq(dbgK, 1) == false {
			cf.Hash[index2].AddFreq(dbgK, 1)
		}
	} */

	return
}

/*func (cf CuckooFilter) GetCount(kb []byte) uint16 {
	hk := sha1.Sum(kb)
	v := hk2uint64(hk)
	index := cf.IndexHash(v)
	fingerprint := FingerHash(v)
	for _, item := range cf.Hash[index].Bkt {
		// fmt.Printf("index: %v, finger: %v\n", index, item.GetFinger())
		if item.GetFinger() == uint16(fingerprint) {
			return item.GetCount()
		}
	}
	// if not return , find another position
	index = cf.AltIndex(index, fingerprint)
	for _, item := range cf.Hash[index].Bkt {
		// fmt.Printf("index: %v, finger: %v\n", index, item.GetFinger())
		if item.GetFinger() == uint16(fingerprint) {
			return item.GetCount()
		}
	}

	// not found in the CuckooFilter
	panic("not found in the CuckooFilter")
}

// allow function return zero version if kb not found in the CuckooFilter
func (cf CuckooFilter) GetCountAllowZero(kb []byte) uint16 {
	hk := sha1.Sum(kb)
	v := hk2uint64(hk)
	// fmt.Printf("[GetCountAllowZero] len(cf.Hash): %d\n", len(cf.Hash))
	// fmt.Printf("[GetCountAllowZero] cf.Hash[0]: %v\n", cf.Hash[0])
	index := cf.IndexHash(v)
	fingerprint := FingerHash(v)
	for _, item := range cf.Hash[index].Bkt {
		// fmt.Printf("index: %v, finger: %v\n", index, item.GetFinger())
		if item.GetFinger() == uint16(fingerprint) {
			return item.GetCount()
		}
	}
	// if not return , find another position
	index = cf.AltIndex(index, fingerprint)
	for _, item := range cf.Hash[index].Bkt {
		// fmt.Printf("index: %v, finger: %v\n", index, item.GetFinger())
		if item.GetFinger() == uint16(fingerprint) {
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
// n, err := cfmmapfp.Write([]byte(cf.Hash[0].Bkt[0]))
// if err != nil {
// 	return err
// }
// if n != cf.numItems*unsafe.Sizeof(bucket) {
// 	log.Fatal("write byte number is not equal cf size")
// }

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
}*/
