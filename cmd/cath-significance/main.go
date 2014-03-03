package main

import (
	"flag"
	"fmt"
	"math"
	"math/rand"
	"runtime/pprof"
	"strconv"
	"strings"
	"sync"
	"time"

	"github.com/BurntSushi/ty/fun"

	"github.com/TuftsBCB/fragbag/bow"

	"github.com/TuftsBCB/tools/util"
)

var flagRandoms = 1000

func init() {
	flag.IntVar(&flagRandoms, "n", flagRandoms,
		"The number of random permutations to use.")
	util.FlagUse("cpu", "cpuprof")
	util.FlagParse(
		"out-file vectors-csv 'comma-sep-cath-domains' ...",
		"Computes the probability that each cluster of cath domains is\n"+
			"significantly different with respect to their Fragbag vectors.")
	util.AssertLeastNArg(3)
}

func main() {
	if len(util.FlagCpuProf) > 0 {
		f := util.CreateFile(util.FlagCpuProf)
		pprof.StartCPUProfile(f)
		defer f.Close()
		defer pprof.StopCPUProfile()
	}
	vectors := readVectors(util.Arg(1))
	groups := readCathGroups(util.Args()[2:])
	out := util.CreateFile(util.Arg(0))
	defer out.Close()

	pf := func(format string, v ...interface{}) {
		// fmt.Printf(format, v...) 
		fmt.Fprintf(out, format, v...)
	}

	type labeledPval struct {
		Name1, Name2 string
		Pval         float32
	}
	b := stdb(vectors, groups)
	pairs := combinations(len(groups))
	dopairs, pvals := make(chan pair), make(chan labeledPval)
	wg := new(sync.WaitGroup)
	for i := 0; i < util.FlagCpu; i++ {
		wg.Add(1)
		go func() {
			for p := range dopairs {
				g1, g2 := groups[p.i], groups[p.j]
				b1, b2 := b[p.i], b[p.j]
				bm1, bm2 := bmean(b1, b2)
				bw := delta(bm1, bm2)

				randws := make([]float32, 1000)
				for i := range randws {
					randws[i] = delta(shuffle_mean_rows(b1, b2))
				}

				bigger := 0
				for _, rw := range randws {
					if rw >= bw {
						bigger++
					}
				}
				pval := float32(bigger) / float32(len(randws))
				pvals <- labeledPval{g1.Name, g2.Name, pval}
			}
			wg.Done()
		}()
	}

	done := make(chan struct{})
	go func() {
		sig, cutoff := 0, 0.05/float32(len(pairs))
		for pval := range pvals {
			if pval.Pval < cutoff {
				sig++
			}
			pf("%s\t%s\t%f\n", pval.Name1, pval.Name2, pval.Pval)
		}
		pf("significant\t%d/%d\t(cutoff: %f)\n", sig, len(pairs), cutoff)
		done <- struct{}{}
	}()
	for _, p := range pairs {
		dopairs <- p
	}
	close(dopairs)
	wg.Wait()
	close(pvals)
	<-done
}

func delta(b1, b2 bow.Bow) float32 {
	if b1.Len() != b2.Len() {
		panic("bow lengths don't match")
	}
	maxdiff := float32(0)
	for i := range b1.Freqs {
		d := float32(math.Abs(float64(b1.Freqs[i] - b2.Freqs[i])))
		if d > maxdiff {
			maxdiff = d
		}
	}
	return maxdiff
}

func bmean(row1, row2 []bow.Bow) (bow.Bow, bow.Bow) {
	m1, m2 := bow.NewBow(row1[0].Len()), bow.NewBow(row2[0].Len())
	if m1.Len() != m2.Len() {
		panic("bow lengths don't match")
	}

	for _, b1 := range row1 {
		for c := range b1.Freqs {
			m1.Freqs[c] += b1.Freqs[c]
		}
	}
	for c := range m1.Freqs {
		m1.Freqs[c] /= float32(len(row1))
	}

	for _, b2 := range row2 {
		for c := range b2.Freqs {
			m2.Freqs[c] += b2.Freqs[c]
		}
	}
	for c := range m2.Freqs {
		m2.Freqs[c] /= float32(len(row2))
	}

	return m1, m2
}

func shuffle_mean_rows(row1, row2 []bow.Bow) (bow.Bow, bow.Bow) {
	all := make([]bow.Bow, len(row1) + len(row2))
	copy(all, row1)
	copy(all[len(row1):], row2)

	rng := rand.New(rand.NewSource(time.Now().UnixNano()))
	fun.ShuffleGen(all, rng)

	shuffled1, shuffled2 := all[:len(row1)], all[len(row1):]
	return bmean(shuffled1, shuffled2)
}

func stdb(vectors map[string]bow.Bow, groups []group) [][]bow.Bow {
	var columns [][]float32
	for _, g := range groups {
		for _, domain := range g.Domains {
			if b, ok := vectors[domain]; ok {
				if columns == nil {
					columns = make([][]float32, b.Len())
				}
				for c, freq := range b.Freqs {
					columns[c] = append(columns[c], freq)
				}
			} else {
				fmt.Printf("Could not find %s\n", domain)
			}
		}
	}

	stddevs := make([]float32, len(columns))
	for i := range columns {
		stddevs[i] = stddev(columns[i])
	}

	var bowgroups [][]bow.Bow
	for _, g := range groups {
		var bows []bow.Bow
		for _, domain := range g.Domains {
			if b, ok := vectors[domain]; ok {
				for i := range b.Freqs {
					b.Freqs[i] /= stddevs[i]
				}
				bows = append(bows, b)
			}
		}
		bowgroups = append(bowgroups, bows)
	}
	return bowgroups
}

type group struct {
	Name    string
	Domains []string
}

func readCathGroups(args []string) []group {
	var groups []group
	for _, arg := range args {
		pieces := strings.Split(arg, ":")
		groups = append(groups, group{pieces[0], strings.Split(pieces[1], ",")})
	}
	return groups
}

func readVectors(fpath string) map[string]bow.Bow {
	f := util.OpenFile(fpath)
	defer f.Close()

	bows := make(map[string]bow.Bow, 5000)
	for _, line := range util.ReadLines(f) {
		fields := strings.Fields(line)
		b := bow.NewBow(len(fields[1:]))
		for _, sfreq := range fields[1:] {
			freq, err := strconv.ParseFloat(sfreq, 32)
			util.Assert(err)
			b.Freqs = append(b.Freqs, float32(freq))
		}
		bows[fields[0]] = b
	}
	return bows
}

func stddev(nums []float32) float32 {
	avg := mean(nums)
	diffsum := float32(0)
	for _, n := range nums {
		diffsum += (avg - n) * (avg - n)
	}
	return float32(math.Sqrt(float64(diffsum / float32(len(nums)))))
}

func mean(nums []float32) float32 {
	sum := float32(0)
	for _, n := range nums {
		sum += n
	}
	return sum / float32(len(nums))
}

type pair struct {
	i, j int
}

func combinations(length int) []pair {
	var ps []pair
	for i := 0; i < length; i++ {
		for j := i + 1; j < length; j++ {
			ps = append(ps, pair{i, j})
		}
	}
	return ps
}
