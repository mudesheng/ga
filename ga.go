package main

import (
	"log"
	"net/http"
	_ "net/http/pprof"

	"github.com/mudesheng/ga/constructcf"
	//"./constructdbg"
	//"./findPath"

	"github.com/jwaldrip/odin/cli"
	"github.com/mudesheng/ga/constructdbg"
	"github.com/mudesheng/ga/deconstructdbg"
	"github.com/mudesheng/ga/preprocess"
)

//"./mapDBG"
// "strconv"
//	"fmt"

const Kmerdef = 203

type GAArgs struct {
	cfg         string
	cpuproffile string
	kmer        int
	prefix      string
	numCPU      int
	cfSize      int64
}

var app = cli.New("1.0.0", "Graph Assembler for complex genome", func(c cli.Command) {})

//var gaargs GAArgs

func init() {
	go func() {
		log.Println(http.ListenAndServe("localhost:6060", nil))
	}()
	app.DefineStringFlag("C", "ga.cfg", "configure file")
	app.DefineStringFlag("cpuprofile", "cpu.prof", "write cpu profile to file")
	app.DefineIntFlag("K", Kmerdef, "kmer length")
	// app.DefineStringFlag("p", "K"+strconv.Itoa(Kmerdef), "prefix of the output file")
	app.DefineStringFlag("p", "./test/t20150708/K203", "prefix of the output file")
	app.DefineIntFlag("t", 1, "number of CPU used")
	pp := app.DefineSubCommand("pp", "correct Illumina sequence reads and link pair end reads to single merged read", preprocess.Correct)
	{
		//pp.DefineInt64Flag("k", 89, "correct cukcoofilter kmer used")
		pp.DefineInt64Flag("S", 0, "the Size number of items cuckoofilter set")
		pp.DefineIntFlag("tipMaxLen", 0, "Maximum tip length, default[0] for MaxNGSReadLen")
		pp.DefineIntFlag("WinSize", 5, "th size of sliding window for DBG edge Sample")
		pp.DefineIntFlag("MaxNGSReadLen", 250, "Max NGS Read Length")
		//pp.DefineIntFlag("tipMaxLen", Kmerdef*2, "Maximum tip length(-K * 2)")
	}
	ccf := app.DefineSubCommand("ccf", "construct cukcoofilter", constructcf.CCF)
	{
		ccf.DefineInt64Flag("S", 0, "the Size number of items cuckoofilter set")
	}
	cdbg := app.DefineSubCommand("cdbg", "construct De bruijn Graph", constructdbg.CDBG)
	{
		cdbg.DefineIntFlag("tipMaxLen", Kmerdef*2, "Maximum tip length(-K * 2)")
	}

	smfy := app.DefineSubCommand("smfy", "find Illumina reads path and simplify De bruijn Graph", constructdbg.Smfy)
	{
		smfy.DefineIntFlag("tipMaxLen", 0, "Maximum tip length, default[0] for MaxNGSReadLen")
		smfy.DefineIntFlag("WinSize", 10, "th size of sliding window for DBG edge Sample")
		smfy.DefineIntFlag("MaxNGSReadLen", 450, "Max NGS Read Length")
		smfy.DefineIntFlag("MinMapFreq", 5, "Minimum reads Mapping Frequent")
		//smfy.DefineIntFlag("MaxMapEdgeLen", 2000, "Max Edge length for mapping Long Reads")
	}
	decontdbg := app.DefineSubCommand("decdbg", "deconstruct DBG using Long Reads Mapping info", deconstructdbg.DeconstructDBG)
	{
		decontdbg.DefineIntFlag("MinCov", 2, "Mininum coverage by long reads")
		decontdbg.DefineIntFlag("WinSize", 10, "th size of sliding window for DBG edge Sample")
		decontdbg.DefineIntFlag("MaxNGSReadLen", 550, "Max NGS Read Length")
		decontdbg.DefineIntFlag("MinMapFreq", 5, "Minimum reads Mapping Frequent")
		decontdbg.DefineIntFlag("ExtLen", 1000, "Extend Path length for distingush most probable path")
		decontdbg.DefineStringFlag("LongReadFile", "ONT.fa", "Oxford Nanopore Technology long reads file")

	}
	// mapping long read to the DBG
	mapDBG := app.DefineSubCommand("mapDBG", "mapping long read to the DBG", constructdbg.MapDBG)
	{
		mapDBG.DefineIntFlag("Seed", 15, "the seek length(must <=16)")
		mapDBG.DefineIntFlag("Width", 5, "band width for found min kmer")
	}
	// find short read mapping
	fspath := app.DefineSubCommand("fspath", "Parse short read path", constructdbg.FSpath)
	{
		fspath.DefineIntFlag("tipMaxLen", Kmerdef*2, "Maximum tip length")
	}
	// find long read mapping
	flpath := app.DefineSubCommand("flpath", "Parse long read path", constructdbg.FLpath)
	{
		flpath.DefineIntFlag("tipMaxLen", Kmerdef*2, "Maximum tip length")
	}
	// merge find short and long read mapping path
	fpath := app.DefineSubCommand("fpath", "Merge Parse short and long read path", constructdbg.Fpath)
	{
		fpath.DefineIntFlag("tipMaxLen", Kmerdef*2, "Maximum tip length")
	}
}

func main() {
	app.Start()
}
