package main

import (
	"github.com/biogo/hts/bam"
	"github.com/golang-collections/collections/set"
	"os"
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

			//// INDEX keep the position of record base position
			//// START keep the genomic position
			//start, index := record.Start, 0
			//for _, i := range record.Cigar {
			//	if i.Type() != sam.CigarDeletion &&
			//		i.Type() != sam.CigarHardClipped &&
			//		i.Type() != sam.CigarInsertion {
			//
			//		if i.Type() == sam.CigarMatch {
			//			for j := 0; j < i.Len(); j++ {
			//
			//				if _, ok := edits[start]; !ok {
			//					edits[start] = NewEditsInfo(ref.Chrom, chrRef[start], start + 1)
			//				}
			//
			//				edits[start].AddReads(record, index)
			//
			//				index++
			//				start++
			//			}
			//		} else {
			//			if i.Type() != sam.CigarSoftClipped {
			//				start += i.Len()
			//			}
			//
			//			if i.Type() != sam.CigarSkipped {
			//				index += i.Len()
			//			}
			//		}
			//	}
			//}

			edits = updateEdits(edits, record, chrRef, ref.Chrom)

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
