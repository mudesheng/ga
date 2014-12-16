package main

import (
	//	"GA/cuckoofilter"
	"GA/constructcf"
	//	"fmt"
	"github.com/jwaldrip/odin/cli"
)

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
	app.DefineIntFlag("K", 57, "kmer length")
	app.DefineStringFlag("p", "ga", "prefix of the output file")
	app.DefineIntFlag("t", 4, "number of CPU used")
	ccf := app.DefineSubCommand("ccf", "construct cukcoofilter", constructcf.CCF)
	{
		ccf.DefineInt64Flag("S", 0, "the number of item cuckoofilter set")

	}

}

func main() {
	app.Start()
}
