package utils

import (
	"fmt"
	"strconv"
	"testing"
	"unsafe"
)

func test(t *testing.T) {

}

func TestBytesEqual(t *testing.T) {
	var a int
	var b ArgsOpt
	var c string
	fmt.Printf("sizeof(int):%d\n", unsafe.Sizeof(a))
	fmt.Printf("sizeof(ArgsOpt):%d\n", unsafe.Sizeof(b))
	fmt.Printf("sizeof(string):%d\n", unsafe.Sizeof(c))
}

func Benchmark_Byte2String(b *testing.B) {
	x := []byte("Hello Gopher! Hello Gopher! Hello Gopher!")
	for i := 0; i < b.N; i++ {
		_ = Bytes2String(x)
	}
}

func Benchmark_BytesEqual(t *testing.B) {
	a := []byte("Gopher!HelloGopher!HelloGopher!Gopher!HelloGopher!HelloGopher!")
	b := []byte("Gopher!HelloGopher!HelloGopher!Gopher!HelloGopher!HelloGopher!")
	for i := 0; i < t.N; i++ {
		BytesEqual(a, b)
	}
}

func Benchmark_BytesEqual2(t *testing.B) {
	a := []byte("Gopher!HelloGopher!HelloGopher!Gopher!HelloGopher!HelloGopher!")
	b := []byte("Gopher!HelloGopher!HelloGopher!Gopher!HelloGopher!HelloGopher!")
	for i := 0; i < t.N; i++ {
		BytesEqual2(a, b)
	}
}

func Benchmark_ByteArrInt(t *testing.B) {
	a := []byte("5432786379334")
	for i := 0; i < t.N; i++ {
		ByteArrInt(a)
	}
}

func Benchmark_strconv(t *testing.B) {
	a := 5432786379334
	for i := 0; i < t.N; i++ {
		strconv.Itoa(a)
	}
}
