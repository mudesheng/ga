package main

import (
	"log"
	"math"
	"net/http"
	_ "net/http/pprof"

	"github.com/mudesheng/ga/constructcf"
	"github.com/mudesheng/ga/constructdbg"
	//"./constructdbg"
	//"./findPath"

	"github.com/jwaldrip/odin/cli"
	//"github.com/mudesheng/ga/constructdbg"
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
		log.Println(http.ListenAndServe("localhost:6090", nil))
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
		pp.DefineIntFlag("tipMaxLen", 250, "Maximum tip length, default for MaxNGSReadLen")
		pp.DefineIntFlag("WinSize", 5, "th size of sliding window for DBG edge Sample")
		pp.DefineIntFlag("MaxNGSReadLen", 250, "Max NGS Read Length")
		pp.DefineBoolFlag("Correct", false, "Correct NGS Read")
		pp.DefineBoolFlag("Merge", false, "Merge pair reads")
		pp.DefineIntFlag("MinKmerFreq", 3, "Min Kmer Freq allown store")
		//pp.DefineIntFlag("tipMaxLen", Kmerdef*2, "Maximum tip length(-K * 2)")
	}
	ccf := app.DefineSubCommand("ccf", "construct cukcoofilter", constructcf.CCF)
	{
		ccf.DefineInt64Flag("S", 0, "the Size number of items cuckoofilter set")
		ccf.DefineBoolFlag("Correct", false, "Correct NGS Read and merge pair reads")
		ccf.DefineBoolFlag("Merge", false, "merge pair reads")
		ccf.DefineIntFlag("MinKmerFreq", 2, "Min Kmer Freq allown store")
	}
	cdbg := app.DefineSubCommand("cdbg", "construct De bruijn Graph", constructdbg.CDBG)
	{
		cdbg.DefineIntFlag("tipMaxLen", 550, "Maximum tip length(-K * 2)")
		cdbg.DefineIntFlag("MinKmerFreq", 2, "Min Kmer Freq allown Extend")
		cdbg.DefineIntFlag("MaxNGSReadLen", 550, "Max NGS Read Length")
	}

	smfy := app.DefineSubCommand("smfy", "Simplify De bruijn Graph", constructdbg.Smfy)
	{
		smfy.DefineIntFlag("TipMaxLen", 550, "Maximum tip length, default[0] for MaxNGSReadLen")
		smfy.DefineIntFlag("WinSize", 5, "th size of sliding window for DBG edge Sample")
		smfy.DefineIntFlag("MaxNGSReadLen", 550, "Max NGS Read Length")
		smfy.DefineIntFlag("MinMapFreq", 3, "Minimum reads Mapping Frequent")
		smfy.DefineBoolFlag("Correct", false, "Correct NGS Read and merge pair reads")
		smfy.DefineBoolFlag("Graph", false, "output dot graph file")
		//smfy.DefineIntFlag("MaxMapEdgeLen", 2000, "Max Edge length for mapping Long Reads")
	}
	mapNGS := app.DefineSubCommand("MapNGS", "find Illumina reads path", constructdbg.MapNGS)
	{
		mapNGS.DefineIntFlag("TipMaxLen", 0, "Maximum tip length, default[0] for MaxNGSReadLen")
		mapNGS.DefineIntFlag("WinSize", 5, "th size of sliding window for DBG edge Sample")
		mapNGS.DefineIntFlag("MaxNGSReadLen", 550, "Max NGS Read Length")
		mapNGS.DefineIntFlag("MinMapFreq", 3, "Minimum reads Mapping Frequent")
		mapNGS.DefineBoolFlag("Correct", false, "Correct NGS Read and merge pair reads")
		mapNGS.DefineBoolFlag("Graph", false, "output dot graph file")
	}

	decontdbg := app.DefineSubCommand("decdbg", "deconstruct DBG using Long Reads Mapping info", deconstructdbg.DeconstructDBG)
	{
		decontdbg.DefineIntFlag("MinDepth", 2, "Mininum coverage by long reads")
		decontdbg.DefineIntFlag("AvgDepth", 40, "average coverage estimate by long reads")
		decontdbg.DefineIntFlag("MinUniqueEdgeLen", 200, "min allow unique edge merge clean bubble")
		decontdbg.DefineIntFlag("AvgReadLen", 12000, "average long read length")
		decontdbg.DefineIntFlag("WinSize", 14, "the size of sliding window for DBG edge Sample[SeedLen* 2/3]")
		decontdbg.DefineIntFlag("MaxNGSReadLen", 550, "Max NGS Read Length")
		decontdbg.DefineIntFlag("MinMapFreq", 2, "Minimum reads Mapping Frequent")
		decontdbg.DefineIntFlag("ExtLen", 1000, "Extend Path length for distingush most probable path")
		decontdbg.DefineFloat64Flag("MinChainScoreIdentityPercent", 0.2, "Min Chain Score Identity Percent[0~1]")
		decontdbg.DefineStringFlag("LongReadFile", "ONT.fa", "Oxford Nanopore Technology long reads file")
		decontdbg.DefineBoolFlag("Correct", false, "Correct NGS Read and merge pair reads[false]")
		decontdbg.DefineBoolFlag("Haplotype", true, "Enable Haplotype model[true]")
		decontdbg.DefineBoolFlag("Debug", false, "Enable Debug model[false]")
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
	// merge find short and long read mapping path
	extractpairend := app.DefineSubCommand("extractpairend", "Extract mixtrue pair end reads to two single file", constructcf.ExtractPairEnd)
	{
		extractpairend.DefineStringFlag("input", "/dev/stdin", "input *.fq.br file name")
		extractpairend.DefineStringFlag("IDFile", "", "read ID list file")
		extractpairend.DefineStringFlag("prefix", "paired_output", "prefix of output file")
		extractpairend.DefineStringFlag("format", "fa", "output format[fq|fa]")
		extractpairend.DefineIntFlag("startID", 1, "start read ID")
		extractpairend.DefineIntFlag("SplitRecordNum", math.MaxUint32, "split big file to mutli files by SplitRecordNum")
		//extractpairend.DefineIntFlag("splitNum", 1, "split output file number")
	}
	filterlong := app.DefineSubCommand("filterlong", "filter long reads", constructcf.FilterLong)
	{
		filterlong.DefineStringFlag("input", "", "input *.fq.gz file name")
		filterlong.DefineStringFlag("output", "output.fa.br", "output file name")
		filterlong.DefineIntFlag("startID", 0, "start read ID")
		filterlong.DefineIntFlag("minLen", 5000, "filter by read seq length")
		filterlong.DefineIntFlag("minMeanQuality", 7, "filter by mean quality")
	}
}

func main() {
	app.Start()
}
