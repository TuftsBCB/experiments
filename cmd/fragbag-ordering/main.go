// test fragbag-ordering does an all-against-all search of the specified BOW
// database, and outputs the ordering of each search.
package main

import (
	"fmt"

	"github.com/TuftsBCB/fragbag/bowdb"
	"github.com/TuftsBCB/tools/util"
)

func init() {
	util.FlagParse("frag-lib-path output-file", "")
	util.AssertNArg(2)
}

func main() {
	db := util.OpenBowDB(util.Arg(0))
	out := util.CreateFile(util.Arg(1))

	printf := func(format string, v ...interface{}) {
		fmt.Fprintf(out, format, v...)
	}

	// Set our search options.
	bowOpts := bowdb.SearchDefault
	bowOpts.Limit = -1

	printf("QueryID\tResultID\tCosine\tEuclid\n")
	entries, err := db.ReadAll()
	util.Assert(err, "Could not read BOW database entries")

	for _, entry := range entries {
		results := db.Search(bowOpts, entry)

		for _, result := range results {
			printf("%s\t%s\t%0.4f\t%0.4f\n",
				entry.Id, result.Bowed.Id, result.Cosine, result.Euclid)
		}
		printf("\n")
	}
	util.Assert(out.Close())
	util.Assert(db.Close())
}

func max(a, b int) int {
	if a > b {
		return a
	}
	return b
}

func min(a, b int) int {
	if a < b {
		return a
	}
	return b
}
