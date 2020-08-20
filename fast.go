package main

import (
	"os"
	"strings"
	"sync"

	"github.com/biogo/hts/bam"

	"github.com/golang-collections/collections/set"
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

		ifh, err := os.Open(conf.File)
		//Panic if something went wrong:
		if err != nil {
			sugar.Fatal(err)
		}
		//defer ifh.Close()

		//Create a new BAM reader with maximum
		//concurrency:
		bamReader, err := bam.NewReader(ifh, 1)
		if err != nil {
			sugar.Fatal(err)
		}
		//defer bamReader.Close()

		iter, err := bam.NewIterator(bamReader, ref.Chunks)
		if err != nil {
			sugar.Fatal(err)
		}
		//defer iter.Close()

		lastEnd := 0
		total := 0
		edits := make(map[int]*EditsInfo)
		for iter.Next() {
			record := NewRecord(iter.Record())

			if record.Start < ref.Start && record.End > ref.End {
				continue
			}

			if lastEnd == 0 {
				lastEnd = record.End
			}

			total++

			if !filterReads(record) {
				continue
			}

			if record.Start > lastEnd {

				if _, ok := edits[151678717]; ok {
					sugar.Infof("lastEnd: %d; record: %v; total: %d; ref: %v", lastEnd, record, total, ref)
				}

				getColumn(edits, []map[string]*set.Set{omopolymericPositions, splicePositions}, targetPositions, w)

				edits = make(map[int]*EditsInfo)
			}

			edits = updateEdits(edits, record, chrRef, ref.Ref)

			if record.End > lastEnd {
				lastEnd = record.End
			}
		}

		iter.Close()
		bamReader.Close()
		ifh.Close()

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
	refChan := make(chan string)

	for i := 0; i < conf.Process; i++ {
		wg.Add(1)

		go func(refChan chan string, wg *sync.WaitGroup, lock *sync.Mutex) {
			defer wg.Done()

			for {
				ref, ok := <-refChan

				if !ok {
					break
				}

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
			}

		}(refChan, wg, &lock)
	}

	for ref := range references {
		refChan <- ref
	}
	close(refChan)
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
			//sugar.Debug(c)
			refs <- c
		}
	}

	close(refs)
}
