package main

import (
	"bufio"
	"fmt"
	"runtime/pprof"
	"strconv"
	"strings"

	"github.com/BurntSushi/intern"
	"github.com/TuftsBCB/tools/util"
)

type inDomains struct {
	in    *intern.Interner
	ids   []string
	atoms []intern.Atom
}

type aligner struct {
	matrix  *intern.Table
	outpath string
}

var domains *inDomains

func init() {
	util.FlagUse("cpu", "cpuprof")
	util.FlagParse(
		"cath-domain-labels best-of-all-matrix"+
			"matrix-file out-file [ matrix-file out-file ... ]",
		"Computes the AUC of each aligner matrix given with respect to the\n"+
			"'best-of-all' matrix given. Each AUC is written to a separate\n"+
			"out-file. The sizes of all matrices must be exactly equivalent.")
	util.AssertLeastNArg(4)
	if util.NArg()%2 != 0 {
		util.Fatalf("There must be an 'out-file' for each 'matrix-file.'")
	}
}

func main() {
	if len(util.FlagCpuProf) > 0 {
		f := util.CreateFile(util.FlagCpuProf)
		pprof.StartCPUProfile(f)
		defer f.Close()
		defer pprof.StopCPUProfile()
	}

	domains = readDomains(util.Arg(0))
	boa := readMatrix(util.Arg(1))
	aligners := make([]aligner, 0)
	for i := 2; i < util.NArg(); i += 2 {
		aligners = append(aligners, aligner{
			readMatrix(util.Arg(i)),
			util.Arg(i + 1),
		})
	}
	fmt.Println(dist(boa, "12asB0", "153l00"))
	fmt.Println(dist(boa, "153l00", "16pk01"))
	fmt.Println(dist(boa, "9gafC2", "9wgaB4"))
	fmt.Println(dist(aligners[0].matrix, "12asB0", "153l00"))
	fmt.Println(dist(aligners[0].matrix, "153l00", "16pk01"))
	fmt.Println(dist(aligners[0].matrix, "9gafC2", "9wgaB4"))
	fmt.Printf("Total CATH domains: %d\n", len(domains.ids))
	domains.removeOldDomains()
	fmt.Printf("Total valid CATH domains: %d\n", len(domains.ids))
}

func dist(matrix *intern.Table, d1, d2 string) float64 {
	return matrix.Get(matrix.Atom(d1), matrix.Atom(d2))
}

func readMatrix(fpath string) *intern.Table {
	var (
		err  error
		fval float64
		sval string
	)
	tab := intern.NewTableInterner(domains.in)
	scanner := bufio.NewScanner(util.OpenFile(fpath))
	for i := 0; scanner.Scan(); i++ {
		// It'd be much simpler to use Split here, but let's be quicker.
		// In particular, avoid allocating.
		// Also, we're dealing with the line as a string since it's quicker
		// than using bytes and converting each number to a string for
		// strconv.ParseFloat.
		line := scanner.Text()
		bstart, j := 0, -1
		for bend, b := range scanner.Text() {
			// This actually skips the very last element in the table, but
			// it's OK because the value at [k, k] is always 0.
			switch {
			case b == ' ' || b == '\n' || bend + 1 == len(line):
				sval = line[bstart:bend]
				bstart = bend + 1
				j++
				// falls down to process this value
			default:
				continue
			}
			if j > i && len(sval) > 0 { // upper triangular
				fval, err = strconv.ParseFloat(sval, 64)
				if err != nil {
					panic(err)
				}
				tab.Set(domains.atoms[i], domains.atoms[j], fval)
			}
		}
	}
	util.Assert(scanner.Err())
	return tab
}

func readDomains(fpath string) *inDomains {
	domains := &inDomains{
		intern.NewInterner(),
		make([]string, 0, 2000),
		make([]intern.Atom, 0, 2000),
	}

	scanner := bufio.NewScanner(util.OpenFile(fpath))
	for scanner.Scan() {
		d := strings.Fields(scanner.Text())[0]
		a := domains.in.Atom(d)
		domains.ids = append(domains.ids, d)
		domains.atoms = append(domains.atoms, a)
	}
	util.Assert(scanner.Err())
	return domains
}

func (d *inDomains) removeOldDomains() {
	ids := make([]string, 0, len(d.ids))
	atoms := make([]intern.Atom, 0, len(d.atoms))
	for i, domain := range d.ids {
		cathPath := util.CathPath(domain)
		if util.Exists(cathPath) {
			ids = append(ids, d.ids[i])
			atoms = append(atoms, d.atoms[i])
		}
	}
	d.ids, d.atoms = ids, atoms
}
