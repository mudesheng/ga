package main

import (
	"math"
	"net/http"
	"net/http/pprof"

	"github.com/jwaldrip/odin/cli"
)

//"./mapDBG"
// "strconv"
//	"fmt"

//Kmerdef set default kmerlen
const Kmerdef int = 203

//GAArgs used collect global args
type gaArgs struct {
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
		r := http.NewServeMux()
		r.Handle("/debug/pprof/", http.HandlerFunc(pprof.Index))
		//log.Println(http.ListenAndServe("localhost:6090", nil))
		http.ListenAndServe("localhost:8090", nil)
	}()
	app.DefineStringFlag("C", "ga.cfg", "configure file")
	app.DefineStringFlag("cpuprofile", "cpu.prof", "write cpu profile to file")
	app.DefineIntFlag("K", Kmerdef, "kmer length")
	// app.DefineStringFlag("p", "K"+strconv.Itoa(Kmerdef), "prefix of the output file")
	app.DefineStringFlag("p", "./test/t20150708/K203", "prefix of the output file")
	app.DefineIntFlag("t", 2, "number of CPU used")
	pp := app.DefineSubCommand("pp", "correct Illumina sequence reads and link pair end reads to single merged read", Correct)
	{
		//pp.DefineInt64Flag("k", 89, "correct cukcoofilter kmer used")
		pp.DefineIntFlag("S", 0, "the Size number of items cuckoofilter set")
		pp.DefineIntFlag("TipMaxLen", 250, "Maximum tip length, default for MaxNGSReadLen")
		pp.DefineIntFlag("WinSize", 8, "th size of sliding window for DBG edge Sample")
		pp.DefineIntFlag("MaxNGSReadLen", 250, "Max NGS Read Length")
		pp.DefineIntFlag("MinNGSReadLen", 250, "Max NGS Read Length")
		pp.DefineBoolFlag("Correct", false, "Correct NGS Read")
		pp.DefineBoolFlag("Merge", false, "Merge pair reads")
		pp.DefineIntFlag("MinKmerFreq", 3, "Min Kmer Freq allown store")
		//pp.DefineIntFlag("tipMaxLen", Kmerdef*2, "Maximum tip length(-K * 2)")
	}
	ccf := app.DefineSubCommand("ccf", "construct cukcoofilter", CCF)
	{
		ccf.DefineIntFlag("S", 0, "the Size number of items cuckoofilter set")
		ccf.DefineBoolFlag("Correct", false, "Correct NGS Read and merge pair reads")
		ccf.DefineBoolFlag("Merge", false, "merge pair reads")
		ccf.DefineIntFlag("MinKmerFreq", 3, "Min Kmer Freq allown store")
	}
	cdbg := app.DefineSubCommand("cdbg", "construct De bruijn Graph", CDBG)
	{
		cdbg.DefineIntFlag("S", 0, "the Size number of items cuckoofilter set")
		cdbg.DefineIntFlag("TipMaxLen", 250, "Maximum tip length(-K * 2)")
		cdbg.DefineIntFlag("MinKmerFreq", 3, "Min Kmer Freq allown Extend")
		cdbg.DefineIntFlag("MaxNGSReadLen", 250, "Max NGS Read Length")
	}

	smfy := app.DefineSubCommand("smfy", "Simplify De bruijn Graph", Smfy)
	{
		smfy.DefineIntFlag("TipMaxLen", 250, "Maximum tip length, default[0] for MaxNGSReadLen")
		smfy.DefineIntFlag("WinSize", 8, "th size of sliding window for DBG edge Sample")
		smfy.DefineIntFlag("MaxNGSReadLen", 250, "Max NGS Read Length")
		smfy.DefineIntFlag("MinMapFreq", 2, "Minimum reads Mapping Frequent")
		smfy.DefineBoolFlag("Correct", false, "Correct NGS Read and merge pair reads")
		//smfy.DefineBoolFlag("ReNameID", false, "Rename edge and node ID, make array compact")
		smfy.DefineBoolFlag("Graph", false, "output dot graph file")
		//smfy.DefineIntFlag("MaxMapEdgeLen", 2000, "Max Edge length for mapping Long Reads")
	}
	mapNGS := app.DefineSubCommand("MapNGS", "find Illumina reads path", MapNGS)
	{
		mapNGS.DefineIntFlag("TipMaxLen", 0, "Maximum tip length, default[0] for MaxNGSReadLen")
		mapNGS.DefineIntFlag("WinSize", 8, "th size of sliding window for DBG edge Sample")
		mapNGS.DefineIntFlag("MaxNGSReadLen", 550, "Max NGS Read Length")
		mapNGS.DefineIntFlag("MinNGSReadLen", 250, "Min NGS Read Length")
		mapNGS.DefineIntFlag("MinMapFreq", 2, "Minimum reads Mapping Frequent")
	}

	decontdbg := app.DefineSubCommand("decdbg", "deconstruct DBG using Long Reads Mapping info", DeconstructDBG)
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
		decontdbg.DefineStringFlag("LongReadFile", "ONT.fa.zst", "Oxford Nanopore Technology long reads file")
		decontdbg.DefineBoolFlag("Correct", false, "Correct NGS Read and merge pair reads[false]")
		decontdbg.DefineBoolFlag("Haplotype", false, "Enable Haplotype model[false]")
		decontdbg.DefineBoolFlag("Debug", false, "Enable Debug model[false]")
	}
	// mapping long read to the DBG
	/*mapDBG := app.DefineSubCommand("mapDBG", "mapping long read to the DBG", MapDBG)
	{
		mapDBG.DefineIntFlag("Seed", 15, "the seek length(must <=16)")
		mapDBG.DefineIntFlag("Width", 5, "band width for found min kmer")
	}
	// find short read mapping
	fspath := app.DefineSubCommand("fspath", "Parse short read path", FSpath)
	{
		fspath.DefineIntFlag("tipMaxLen", Kmerdef*2, "Maximum tip length")
	}
	// find long read mapping
	flpath := app.DefineSubCommand("flpath", "Parse long read path", FLpath)
	{
		flpath.DefineIntFlag("tipMaxLen", Kmerdef*2, "Maximum tip length")
	}
	// merge find short and long read mapping path
	fpath := app.DefineSubCommand("fpath", "Merge Parse short and long read path", Fpath)
	{
		fpath.DefineIntFlag("tipMaxLen", Kmerdef*2, "Maximum tip length")
	}*/
	// merge find short and long read mapping path
	extractpairend := app.DefineSubCommand("extractpairend", "Extract mixtrue pair end reads to two single file", ExtractPairEnd)
	{
		extractpairend.DefineStringFlag("read1", "/dev/stdin", "input *1.fq.gz file name")
		extractpairend.DefineStringFlag("read2", "/dev/stdin", "input *2.fq.gz file name")
		extractpairend.DefineStringFlag("IDFile", "", "read ID list file")
		extractpairend.DefineStringFlag("prefix", "paired_output", "prefix of output file")
		extractpairend.DefineStringFlag("format", "fa", "output format[fq|fa]")
		extractpairend.DefineIntFlag("startID", 1, "start read ID")
		extractpairend.DefineIntFlag("SplitRecordNum", math.MaxUint32, "split big file to mutli files by SplitRecordNum")
		//extractpairend.DefineIntFlag("splitNum", 1, "split output file number")
	}
	filterlong := app.DefineSubCommand("filterlong", "filter long reads", FilterLong)
	{
		filterlong.DefineStringFlag("input", "", "input *.fq.gz file name")
		filterlong.DefineStringFlag("output", "output.fa.br", "output file name")
		filterlong.DefineIntFlag("startID", 0, "start read ID")
		filterlong.DefineIntFlag("minLen", 5000, "filter by read seq length")
		filterlong.DefineIntFlag("minMeanQuality", 7, "filter by mean quality")
	}
	addONTQC := app.DefineSubCommand("addONTQC", "add ONT reads Quality Score", AddONTQC)
	{
		addONTQC.DefineStringFlag("input", "", "input *.fa.zst file name")
		addONTQC.DefineStringFlag("output", "output.fa.zst", "output file name")
		addONTQC.DefineStringFlag("paf", "", "ONT paf file")

	}
	simulateNGS := app.DefineSubCommand("simulateNGS", "simulate NGS reads", SimulateNGS)
	{
		simulateNGS.DefineStringFlag("input", "", "input *.fa.zst file name")
		simulateNGS.DefineStringFlag("output", "output.fa.zst", "output file name")
		simulateNGS.DefineIntFlag("ID", 1, "read start ID")
		simulateNGS.DefineIntFlag("Step", 10, "Step move by reference seq")

	}
}

func main() {
	app.Start()
}
