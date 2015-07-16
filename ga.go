package main

import (
	//	"GA/cuckoofilter"
	"ga/constructcf"
	"ga/constructdbg"
	//	"fmt"
	"github.com/jwaldrip/odin/cli"
)

const Kmerdef = 57

type GAArgs struct {
	cfg    string
	kmer   int
	prefix string
	numCPU int
	cfSize int64
}

var app = cli.New("1.0.0", "Graph Assembler for complex genome", func(c cli.Command) {})

//var gaargs GAArgs

func init() {
	app.DefineStringFlag("C", "ga.cfg", "configure file")
	app.DefineIntFlag("K", Kmerdef, "kmer length")
	app.DefineStringFlag("p", "ga", "prefix of the output file")
	app.DefineIntFlag("t", 4, "number of CPU used")
	ccf := app.DefineSubCommand("ccf", "construct cukcoofilter", constructcf.CCF)
	{
		ccf.DefineInt64Flag("S", 0, "the Size number of items cuckoofilter set")
	}
	cdbg := app.DefineSubCommand("cdbg", "construct De bruijn Graph", constructdbg.CDBG)
	{
		cdbg.DefineIntFlag("tipMaxLen", Kmerdef*2, "Maximum tip length")
	}

	smfy := app.DefineSubCommand("cdbg", "construct De bruijn Graph", constructdbg.Smfy)
	{
		smfy.DefineIntFlag("tipMaxLen", Kmerdef*2, "Maximum tip length")
	}

}

func main() {
	app.Start()
}
