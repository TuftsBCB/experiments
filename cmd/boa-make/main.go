package main

import (
	"bufio"
	"bytes"
	"fmt"
	"io"
	"strconv"

	"github.com/TuftsBCB/tools/util"
)

func init() {
	util.FlagParse(
		"make-best-of-all out-file matrix-file [ matrix-file ... ]",
		"Combines the input matrices into a 'best-of-all' matrix by\n"+
			"using the best (lowest) SAS score for each pair. All given\n"+
			"matrices must be exactly the same size.")
	util.AssertLeastNArg(2)
}

func main() {
	saveto := util.CreateFile(util.Arg(0))
	defer saveto.Close()

	w := func(format string, v ...interface{}) {
		_, err := fmt.Fprintf(saveto, format, v...)
		util.Assert(err)
	}

	var fmats []*bufio.Reader
	for _, fmat := range util.Args()[1:] {
		fmats = append(fmats, bufio.NewReader(util.OpenFile(fmat)))
	}
LOOP:
	for {
		var columns int
		scores := make([][]float64, len(fmats)) // matrix -> fields -> sas score
		for i, fmat := range fmats {
			line, err := fmat.ReadBytes('\n')
			if len(line) == 0 && err == io.EOF {
				break LOOP
			} else if err != io.EOF {
				util.Assert(err)
			}

			fields := bytes.Fields(line)
			columns = len(fields)
			scores[i] = make([]float64, columns)
			for j, sas := range fields {
				scores[i][j], err = strconv.ParseFloat(string(sas), 64)
				util.Assert(err)
			}
		}

		before := ""
		for j := 0; j < columns; j++ {
			best := scores[0][j]
			for i := 1; i < len(scores); i++ {
				if scores[i][j] < best {
					best = scores[i][j]
				}
			}
			if best == 0 {
				w("%s0", before)
			} else {
				w("%s%f", before, best)
			}
			before = " "
		}
		w("\n")
	}
}
