package main

import (
	"bufio"
	"compress/gzip"
	"os"
	"strings"
	"sync"

	"github.com/biogo/hts/sam"
	"github.com/golang-collections/collections/set"
	"github.com/voxelbrain/goptions"
)

func writer(w chan string, wg *sync.WaitGroup) {
	sugar.Infof("write into %s", conf.Output)
	mode := os.O_CREATE | os.O_WRONLY
	if conf.Append {
		mode = mode | os.O_APPEND
	} else {
		mode = mode | os.O_TRUNC
	}

	f, err := os.OpenFile(conf.Output, mode, 0644)
	if err != nil {
		sugar.Fatalf("failed to open %s: %s", conf.Output, err.Error())
	}
	defer f.Close()

	var gwriter *gzip.Writer
	var writer *bufio.Writer

	if strings.HasSuffix(conf.Output, "gz") {
		gwriter = gzip.NewWriter(f)
		defer gwriter.Flush()
		defer gwriter.Close()
	} else {
		writer = bufio.NewWriter(f)
		defer writer.Flush()
	}

	if !conf.RemoveHeader {
		if gwriter != nil {
			_, err = gwriter.Write([]byte(strings.Join(getHeader(), "\t") + "\n"))
		} else {
			_, err = writer.WriteString(strings.Join(getHeader(), "\t") + "\n")
		}

		if err != nil {
			sugar.Fatal(err)
		} else {
			sugar.Debug("write title")
		}
	}

	done := 0

	for {
		line, ok := <-w
		if !ok {
			break
		}

		if line == "done" {
			done++

			if done >= conf.Process {
				break
			} else {
				continue
			}
		}

		if gwriter == nil {
			_, _ = writer.WriteString(line + "\n")
		} else {
			gwriter.Write([]byte(line + "\n"))
		}

	}

	_ = f.Sync()
	wg.Done()
}

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

						// record.QueryPosition = append(record.QueryPosition, at)

						genomic := start + record.Start

						if _, ok := edits[genomic]; !ok {
							edits[genomic] = NewEditsInfo(ref.Chrom, chrRef[genomic-1], genomic)
						}

						edits[genomic].AddReads(record, at)
						start++
					}
					index += i.Len()
				} else if i.Type() != sam.CigarDeletion || i.Type() != sam.CigarHardClipped || i.Type() != sam.CigarInsertion {
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

func main() {
	conf = defaultConfig()
	goptions.ParseAndFail(&conf)

	setLogger(conf.Debug, conf.Log)

	if conf.Version {
		sugar.Infof("current version: %v", VERSION)
		os.Exit(0)
	}

	if conf.File == "" {
		sugar.Fatal("bam file is mandatory. Please, provide one (-f|--file)")
	}

	if conf.Reference == "" {
		sugar.Fatal("An input reference file is mandatory. Please, provide one (-r|--reference)")
	}

	if conf.Dna && conf.BedFile == "" {
		sugar.Fatal("When analyzing DNA-Seq files it is mandatory to provide a BED file containing the positions of target regions (-B|--bed_file)")
	}

	timer := InitTimer()
	timer.Tic()
	region = decodeRegion(conf.Region)

	targetPositions, err := loadBedFile()
	if err != nil {
		sugar.Error(err)
	}

	if conf.CreateOmopolymericFile {
		if err := createOmopolymericPositions(); err != nil {
			sugar.Error(err)
		}
	}

	omopolymericPositions, err := loadOmopolymericPositions()
	if err != nil {
		sugar.Error(err)
	}

	spicePositions, err := loadSplicingPositions()
	if err != nil {
		sugar.Error(err)
	}

	sugar.Infof("Fetching data from bam %s", conf.File)
	sugar.Infof("Narrowing REDItools to region %s", region.String())

	var wg sync.WaitGroup
	w := make(chan string)
	wg.Add(1)
	go writer(w, &wg)
	refs := make(chan *Region)

	for i := 0; i < conf.Process; i++ {
		go worker(&wg, refs, w, omopolymericPositions, spicePositions, targetPositions)
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
	wg.Wait()

	timer.Toc()
	sugar.Info(timer.TicToc())
}
