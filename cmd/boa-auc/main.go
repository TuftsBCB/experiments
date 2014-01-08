package main

import (
	"bufio"
	"flag"
	"fmt"
	"math"
	path "path/filepath"
	"runtime/pprof"
	"strconv"
	"strings"

	"github.com/BurntSushi/intern"
	"github.com/BurntSushi/ty/fun"

	"github.com/TuftsBCB/fragbag/bow"
	"github.com/TuftsBCB/fragbag/bowdb"
	"github.com/TuftsBCB/tools/util"
)

var flagThreshold = float64(0)

var bowdbSearch = bowdb.SearchOptions{
	Limit:  -1,
	Min:    0,
	Max:    math.MaxFloat64,
	SortBy: bowdb.SortByCosine,
	Order:  bowdb.OrderAsc,
}

type inDomains struct {
	in    *intern.Interner
	ids   []string
	atoms []intern.Atom
}

type aligner struct {
	matrix  *intern.Table
	outpath string
}

type flib struct {
	db      *bowdb.DB
	bowed   []bow.Bowed // indexed by atom
	outpath string
}

func init() {
	flag.Float64Var(&flagThreshold, "threshold", flagThreshold,
		"Set the distance threshold to use when computing AUC.")
	util.FlagUse("cpu", "cpuprof")
	util.FlagParse(
		"cath-domain-labels best-of-all-matrix"+
			"(bowdb | matrix-file) out-file "+
			"[ (bowdb | matrix-file) out-file ... ]",
		"Computes the AUC of each aligner matrix (or BOW database) given\n"+
			"with respect to the 'best-of-all' matrix given. Each AUC is\n"+
			"written to a separate out-file. The sizes of all matrices must\n"+
			"be exactly equivalent.\n"+
			"Files are interpreted as BOW databases if they have a '.bowdb'\n"+
			"file extension.")
	util.AssertLeastNArg(4)
	if util.NArg()%2 != 0 {
		util.Fatalf("There must be an out file for each matrix or bowdb file.")
	}
}

func main() {
	if len(util.FlagCpuProf) > 0 {
		f := util.CreateFile(util.FlagCpuProf)
		pprof.StartCPUProfile(f)
		defer f.Close()
		defer pprof.StopCPUProfile()
	}

	// Read all CATH domains, the best-of-all matrix, and the matrix for
	// each aligner.
	domains := readDomains(util.Arg(0))
	boa := readMatrix(domains, util.Arg(1))
	aligners := make([]aligner, 0)
	flibs := make([]flib, 0)
	for i := 2; i < util.NArg(); i += 2 {
		fpath := util.Arg(i)
		if path.Ext(fpath) == ".bowdb" {
			db := util.OpenBowDB(fpath)
			records, err := db.ReadAll()
			util.Assert(err)

			bowed := make([]bow.Bowed, domains.in.Len())
			for _, b := range records {
				if !domains.in.Exists(b.Id) {
					util.Fatalf("Found ID in bowdb that isn't in the list "+
						"of CATH domains provided: %s", b.Id)
				}
				bowed[domains.in.Atom(b.Id)] = b
			}
			flibs = append(flibs, flib{db, bowed, util.Arg(i + 1)})
		} else {
			aligners = append(aligners, aligner{
				readMatrix(domains, fpath),
				util.Arg(i + 1),
			})
		}
	}
	// Now remove CATH domains that don't have a corresponding structure file.
	// We don't do this initially since the matrix files are indexed with
	// respect to all CATH domains (includings ones without structure).
	// This is an artifact of the fact that the matrices were generated with
	// a very old version of CATH.
	domains.removeOldDomains()

	if a := matrixAuc(domains, boa, boa, flagThreshold); a != 1.0 {
		util.Fatalf("Something is wrong. The AUC of the best-of-all matrix "+
			"with respect to itself is %f, but it should be 1.0.", a)
	}

	if len(aligners) > 0 {
		fmt.Println("Computing AUC for aligners...")
		writeAuc := func(aligner aligner) struct{} {
			w := util.CreateFile(aligner.outpath)
			a := matrixAuc(domains, boa, aligner.matrix, flagThreshold)
			fmt.Fprintf(w, "%f\n", a)
			return struct{}{}
		}
		fun.ParMap(writeAuc, aligners)
	}
	if len(flibs) > 0 {
		fmt.Println("Computing AUC for bowdbs...")
		writeAuc := func(flib flib) struct{} {
			w := util.CreateFile(flib.outpath)
			a := flibAuc(domains, boa, flib, flagThreshold)
			fmt.Fprintf(w, "%f\n", a)
			return struct{}{}
		}
		fun.ParMap(writeAuc, flibs)
	}
}

func dist(matrix *intern.Table, d1, d2 string) float64 {
	return matrix.Get(matrix.Atom(d1), matrix.Atom(d2))
}

func readMatrix(domains *inDomains, fpath string) *intern.Table {
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
			case b == ' ' || b == '\n' || bend+1 == len(line):
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
		d = stripExt(path.Base(util.CathPath(d)))
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

func stripExt(s string) string {
	i := strings.Index(s, ".")
	if i == -1 {
		return s
	}
	return s[0:i]
}
