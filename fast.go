package main

import (
	"github.com/biogo/hts/sam"
	"github.com/golang-collections/collections/set"
	"strings"
	"sync"
)

func workerFast(
	wg *sync.WaitGroup, refs chan *ChanChunk, w chan string,
	omopolymericPositions, splicePositions, targetPositions map[string]*set.Set,
	chrRefs map[string][]byte) {
	defer wg.Done()

	for {
		ref, ok := <-refs
		if !ok {
			break
		}

		chrRef, ok := chrRefs[ref.Ref]
		if !ok {
			sugar.Errorf("failed to get %s from reference", ref.Ref)
			continue
		}

		iter, err := fetchBamFast(ref.Chunks)
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
							// GL000220.1
							if genomic >= len(chrRef) {
								sugar.Error(record.Record)
								sugar.Fatalf("%s: genomic - 1[%d] >= len(chrRef)[%d]", ref.Ref, genomic, len(chrRef))
							}
							edits[genomic] = NewEditsInfo(ref.Ref, chrRef[genomic], genomic + 1)
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
	}

	w <- "done"
}

func fastMode(
	wg *sync.WaitGroup, w chan string,
	omopolymericPositions, spicePositions, targetPositions map[string]*set.Set) {
	var lock sync.Mutex
	refs := make(chan *ChanChunk)

	references, err := fetchBamRefsFast()
	if err != nil {
		sugar.Fatal(err)
	}

	sugar.Infof("load reference from %s", conf.Reference)
	chrRefs := make(map[string][]byte)
	for ref, _ := range references {
		wg.Add(1)
		go func(ref string, wg *sync.WaitGroup, lock *sync.Mutex) {
			defer wg.Done()
			temp, err := fetchFasta(&Region{Chrom: ref})
			if err != nil {
				sugar.Warnf("try to modify %s", ref)
				if strings.HasPrefix(ref, "chr") {
					temp, err = fetchFasta(&Region{Chrom: strings.ReplaceAll(ref, "chr", "")})
				} else {
					temp, err = fetchFasta(&Region{Chrom: "chr" + ref})
				}
			}
			if err != nil {
				sugar.Fatal(err)
			}
			lock.Lock()
			chrRefs[ref] = temp
			lock.Unlock()
		}(ref, wg, &lock)
	}

	wg.Wait()

	wg.Add(1)
	go writer(w, wg)

	for i := 0; i < conf.Process; i++ {
		go workerFast(wg, refs, w, omopolymericPositions, spicePositions, targetPositions, chrRefs)
		wg.Add(1)
	}

	for ref, chunks := range references {
		sugar.Infof("read reads from %s", ref)
		for _, c := range chunks {
			sugar.Debug(c)
			refs <- c
		}
	}

	close(refs)
}
