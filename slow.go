package main

import (
	"github.com/biogo/hts/sam"
	"github.com/golang-collections/collections/set"
	"strings"
	"sync"
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

		iter, err := fetchBam(ref)
		if err != nil {
			sugar.Fatal(err)
		}

		lastEnd := 0
		total := 0
		edits := make(map[int]*EditsInfo)
		for iter.Next() {
			record := NewRecord(iter.Record())

			if lastEnd == 0 {
				lastEnd = record.End
			}

			total++

			if !filterReads(record) {
				continue
			}

			start, index := 0, 0
			for _, i := range record.Cigar {
				if i.Type() == sam.CigarMatch {
					for j := 1; j <= i.Len(); j++ {
						at := index + j - 1

						if at >= record.Seq.Length {
							sugar.Errorf("%s - %d - %d - %d", record.Cigar, index, at, record.Seq.Length)
							sugar.Fatal(record.SeqString())
						}

						genomic := start + record.Start

						if _, ok := edits[genomic]; !ok {
							edits[genomic] = NewEditsInfo(ref.Chrom, chrRef[genomic], genomic + 1)
						}

						edits[genomic].AddReads(record, at)
						start++
					}
					index += i.Len()
				} else if i.Type() != sam.CigarDeletion &&
					i.Type() != sam.CigarHardClipped &&
					i.Type() != sam.CigarInsertion &&
					i.Type() != sam.CigarSoftClipped {
					start += i.Len()
				}
			}

			if record.Start > lastEnd {
				getColumn(edits, []map[string]*set.Set{omopolymericPositions, splicePositions}, targetPositions, w)

				edits = make(map[int]*EditsInfo)
				lastEnd = record.End
			}
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
