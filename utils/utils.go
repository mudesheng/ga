package utils

import (
	"log"

	"github.com/jwaldrip/odin/cli"
	//"fmt"
	//"os"
)

type ArgsOpt struct {
	Prefix     string
	Kmer       int
	NumCPU     int
	CfgFn      string
	Cpuprofile string
}

// return global arguments and check if successed
func CheckGlobalArgs(c cli.Command) (opt ArgsOpt, succ bool) {
	opt.Prefix = c.Flag("p").String()
	if opt.Prefix == "" {
		log.Fatalf("[CheckGlobalArgs] args 'p' not set\n")
	}
	opt.CfgFn = c.Flag("C").String()
	if opt.CfgFn == "" {
		log.Fatalf("[CheckGlobalArgs] args 'C' not set\n")
	}
	opt.Cpuprofile = c.Flag("cpuprofile").String()
	if opt.Cpuprofile == "" {
		log.Fatalf("[CheckGlobalArgs] args 'cpuprofile' not set\n")
	}

	var ok bool
	opt.Kmer, ok = c.Flag("K").Get().(int)
	if !ok {
		log.Fatalf("[CheckGlobalArgs] args 'K' : %v set error\n", c.Flag("K").String())
	}
	opt.NumCPU, ok = c.Flag("t").Get().(int)
	if !ok {
		log.Fatalf("[CheckGlobalArgs] args 't': %v set error\n", c.Flag("t").String())
	}
	return opt, true
}

func AbsInt(a int) int {
	if a < 0 {
		return -a
	} else {
		return a
	}
}
