package main

import _ "net/http/pprof"
import (
	//	"GA/cuckoofilter"
	"ga/constructcf"
	"ga/constructdbg"
	"ga/findPath"
	"log"
	"net/http"
	// "strconv"
	//	"fmt"
	"github.com/jwaldrip/odin/cli"
)

const Kmerdef = 203

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
	go func() {
		log.Println(http.ListenAndServe("localhost:6060", nil))
	}()
	app.DefineStringFlag("C", "ga.cfg", "configure file")
	app.DefineIntFlag("K", Kmerdef, "kmer length")
	// app.DefineStringFlag("p", "K"+strconv.Itoa(Kmerdef), "prefix of the output file")
	app.DefineStringFlag("p", "./test/t20150708/K203", "prefix of the output file")
	app.DefineIntFlag("t", 1, "number of CPU used")
	ccf := app.DefineSubCommand("ccf", "construct cukcoofilter", constructcf.CCF)
	{
		ccf.DefineInt64Flag("S", 0, "the Size number of items cuckoofilter set")
	}
	cdbg := app.DefineSubCommand("cdbg", "construct De bruijn Graph", constructdbg.CDBG)
	{
		cdbg.DefineIntFlag("tipMaxLen", Kmerdef*2, "Maximum tip length")
	}

	smfy := app.DefineSubCommand("smfy", "construct De bruijn Graph", constructdbg.Smfy)
	{
		smfy.DefineIntFlag("tipMaxLen", Kmerdef*2, "Maximum tip length")
	}
	// find short read mapping
	fspath := app.DefineSubCommand("fspath", "Parse short read path", findPath.FSpath)
	{
		fspath.DefineIntFlag("tipMaxLen", Kmerdef*2, "Maximum tip length")
	}
	// find long read mapping
	flpath := app.DefineSubCommand("flpath", "Parse long read path", findPath.FLpath)
	{
		flpath.DefineIntFlag("tipMaxLen", Kmerdef*2, "Maximum tip length")
	}
	// merge find short and long read mapping path
	fpath := app.DefineSubCommand("fpath", "Merge Parse short and long read path", findPath.Fpath)
	{
		fpath.DefineIntFlag("tipMaxLen", Kmerdef*2, "Maximum tip length")
	}
}

func main() {
	app.Start()
}
