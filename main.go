package main

import (
	"bufio"
	"compress/gzip"
	"os"
	"strings"
	"sync"

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
		os.Exit(1)
	}
	//defer f.Close()

	var writer *bufio.Writer
	var gwriter *gzip.Writer
	if strings.HasSuffix(conf.Output, "gz") {
		gwriter = gzip.NewWriter(f)
		writer = bufio.NewWriter(gwriter)
	} else {
		writer = bufio.NewWriter(f)
	}
	//defer writer.Flush()

	if !conf.RemoveHeader {

		_, err = writer.WriteString(strings.Join(getHeader(), "\t") + "\n")

		if err != nil {
			sugar.Fatal(err)
		} else {
			sugar.Debug("write title")
		}
	}

	done := 0

	for {
		line := <-w

		if line == "done" {
			done++

			if done >= conf.Process {
				break
			} else {
				continue
			}
		}

		_, _ = writer.WriteString(line + "\n")
	}

	writer.Flush()

	if gwriter != nil {
		gwriter.Flush()
		gwriter.Close()
	}

	f.Close()
	wg.Done()
}

func main() {
	conf = defaultConfig()
	goptions.ParseAndFail(&conf)

	setLogger(conf.Debug, conf.Log)

	sugar.Debugf("options: %v", conf)

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

	if conf.Mode {
		slowMode(&wg, w, omopolymericPositions, spicePositions, targetPositions)
	} else {
		slowMode(&wg, w, omopolymericPositions, spicePositions, targetPositions)
	}

	wg.Wait()

	timer.Toc()
	sugar.Info(timer.TicToc())
}
