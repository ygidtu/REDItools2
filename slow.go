package main

import (
	"os"
	"strings"
	"sync"

	"github.com/biogo/hts/bam"
	"github.com/golang-collections/collections/set"
)

func worker(wg *sync.WaitGroup, refs chan *Region, w chan string, omopolymericPositions, splicePositions, targetPositions map[string]*set.Set) {
	defer wg.Done()

	for {
		ref, ok := <-refs
		if !ok {
			break
		}

		chrRef, err := fetchFasta(ref)
		if err != nil {
			sugar.Warnf("try to modify %s", ref.Chrom)
			if strings.HasPrefix(ref.Chrom, "chr") {
				chrRef, err = fetchFasta(&Region{Chrom: strings.ReplaceAll(ref.Chrom, "chr", "")})
			} else {
				chrRef, err = fetchFasta(&Region{Chrom: "chr" + ref.Chrom})
			}
		}
		if err != nil {
			sugar.Fatal(err)
		}

		chunks, err := fetchBam(ref)
		if err != nil {
			sugar.Fatal(err)
		}

		ifh, err := os.Open(conf.File)
		//Panic if something went wrong:
		if err != nil {
			sugar.Fatal(err)
		}
		defer ifh.Close()

		//Create a new BAM reader with maximum
		//concurrency:
		bamReader, err := bam.NewReader(ifh, 1)
		if err != nil {
			sugar.Fatal(err)
		}
		defer bamReader.Close()

		iter, err := bam.NewIterator(bamReader, chunks)
		if err != nil {
			sugar.Fatal(err)
		}

		defer iter.Close()

		lastEnd := 0
		total := 0
		edits := make(map[int]*EditsInfo)
		for iter.Next() {
			rec := iter.Record()
			record := NewRecord(rec)

			if err := iter.Error(); err != nil {
				sugar.Fatal(err)
			}

			if lastEnd == 0 {
				lastEnd = record.End
			}

			total++

			if !filterReads(record) {
				continue
			}

			edits = updateEdits(edits, record, chrRef, ref.Chrom)

			// if record.Start > lastEnd {
			// 	getColumn(edits, []map[string]*set.Set{omopolymericPositions, splicePositions}, targetPositions, w)

			// 	edits = make(map[int]*EditsInfo)
			// 	lastEnd = record.End
			// }
		}

		getColumn(edits, []map[string]*set.Set{omopolymericPositions, splicePositions}, targetPositions, w)

		sugar.Debugf("read %d reads from %s", total, ref)
	}

	w <- "done"
}

func slowMode(wg *sync.WaitGroup, w chan string,
	omopolymericPositions, spicePositions, targetPositions map[string]*set.Set) {

	wg.Add(1)
	go writer(w, wg)

	refs := make(chan *Region)

	for i := 0; i < conf.Process; i++ {
		go worker(wg, refs, w, omopolymericPositions, spicePositions, targetPositions)
		wg.Add(1)
	}

	if region.Empty() {
		if references, err := fetchBamRefs(); err == nil {
			for _, r := range references {
				refs <- &Region{Chrom: r}
			}
		}
	} else {
		refs <- region
	}

	close(refs)
}
