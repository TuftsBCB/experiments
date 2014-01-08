package main

import (
	"sort"

	"github.com/BurntSushi/intern"
)

func matrixAuc(
	domains *inDomains,
	boa, test *intern.Table,
	threshold float64,
) float64 {
	sumAuc := float64(0)
	for _, query := range domains.atoms {
		tps, tns := tpsAndTns(domains, boa, query, threshold)
		ranked := rankMatrix(domains, test, query)
		sumAuc += auc(rocCurve(ranked, tps, tns))
	}
	return sumAuc / float64(len(domains.atoms))
}

func flibAuc(
	domains *inDomains,
	boa *intern.Table,
	test flib,
	threshold float64,
) float64 {
	sumAuc := float64(0)
	for _, query := range domains.atoms {
		tps, tns := tpsAndTns(domains, boa, query, threshold)
		ranked := rankFlib(domains, test, query)
		sumAuc += auc(rocCurve(ranked, tps, tns))
	}
	return sumAuc / float64(len(domains.atoms))
}

type point struct {
	x, y float64
}

func auc(points []point) (a float64) {
	for i := 1; i < len(points); i++ {
		pt0, pt1 := points[i-1], points[i]
		a += (pt1.x - pt0.x) * (pt1.y + pt0.y)
	}
	return a / 2.0
}

func rocCurve(ranked []intern.Atom, tps, tns []bool) []point {
	tpsLen, tnsLen := atomSetLen(tps), atomSetLen(tns)
	points := make([]point, len(ranked)+1)
	tpIntLen, tnIntLen := 0, 0
	for i := 0; i < len(points); i++ {
		if i > 0 {
			h := ranked[i-1]
			if tps[h] {
				tpIntLen++
			} else if tns[h] {
				tnIntLen++
			}
		}
		tpRate := rate(tpsLen, tpIntLen)
		tnRate := rate(tnsLen, tnIntLen)
		points[i] = point{tnRate, tpRate}
	}
	return points
}

func tpsAndTns(
	domains *inDomains,
	boa *intern.Table,
	query intern.Atom,
	threshold float64,
) (tps []bool, tns []bool) {
	tps = make([]bool, domains.in.Len())
	tns = make([]bool, domains.in.Len())
	for _, query2 := range domains.atoms {
		if boa.Get(query, query2) <= threshold {
			tps[query2] = true
		} else {
			tns[query2] = true
		}
	}
	return
}

func atomSetLen(set []bool) (count int) {
	for _, b := range set {
		if b {
			count++
		}
	}
	return
}

func rate(goldLen, intersectionLen int) float64 {
	if goldLen == 0 || intersectionLen == 0 {
		return 0
	}
	return float64(intersectionLen) / float64(goldLen)
}

func intersectionLens(tps, tns []bool, hits []intern.Atom) (ctps, ctns int) {
	for _, hit := range hits {
		if tps[hit] {
			ctps++
		} else if tns[hit] {
			ctns++
		}
	}
	return
}

type rankedAtoms struct {
	atoms []intern.Atom
	dists []float64
}

func (r *rankedAtoms) Len() int {
	return len(r.atoms)
}

func (r *rankedAtoms) Less(i, j int) bool {
	return r.dists[i] < r.dists[j]
}

func (r *rankedAtoms) Swap(i, j int) {
	r.atoms[i], r.atoms[j] = r.atoms[j], r.atoms[i]
	r.dists[i], r.dists[j] = r.dists[j], r.dists[i]
}

func rankMatrix(
	domains *inDomains,
	matrix *intern.Table,
	query intern.Atom,
) []intern.Atom {
	r := &rankedAtoms{
		atoms: make([]intern.Atom, len(domains.atoms)),
		dists: make([]float64, len(domains.atoms)),
	}
	for i, query2 := range domains.atoms {
		r.atoms[i] = query2
		r.dists[i] = matrix.Get(query, query2)
	}
	sort.Sort(r)
	return r.atoms
}

func rankFlib(
	domains *inDomains,
	db flib,
	query intern.Atom,
) []intern.Atom {
	atoms := make([]intern.Atom, 0, len(domains.atoms))
	for _, sr := range db.db.Search(bowdbSearch, db.bowed[query]) {
		atoms = append(atoms, domains.in.Atom(sr.Id))
	}
	return atoms
}
