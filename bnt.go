package main

const (
	NumBaseInByte   = 4
	NumBaseInUint64 = NumBaseInByte * 8
	NumBitsInBase   = 2
	BaseTypeNum     = 4
	// BaseMask      = ((1 << NumBitsInBase) - 1) // 0x3
	BaseMask = 0x3
)

var Base2Bnt = [256]uint8{
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 1, 4, 0, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 1, 4, 0, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4}

var BitNtCharLow = [4]byte{'c', 'a', 't', 'g'}
var BitNtCharUp = [4]byte{'C', 'A', 'T', 'G'}
var BntRev = [4]uint8{3, 2, 1, 0}
var BitNtRev = [4]byte{'G', 'T', 'A', 'C'}

func TransformSeq(bntArr []byte) (seq string) {
	ta := make([]byte, len(bntArr))
	for i, b := range bntArr {
		ta[i] = BitNtCharUp[b]
	}
	return string(ta)
}
