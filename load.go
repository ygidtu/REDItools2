package main

import (
	"bufio"
	"compress/gzip"
	"github.com/biogo/biogo/alphabet"
	"github.com/biogo/biogo/io/seqio"
	"github.com/biogo/biogo/io/seqio/fasta"
	"github.com/biogo/biogo/seq/linear"
	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/bgzf"
	"github.com/biogo/hts/fai"
	"github.com/cheggaaa/pb/v3"
	"github.com/golang-collections/collections/set"
	"github.com/pkg/errors"
	"io"
	"math"
	"os"
	"regexp"
	"strconv"
	"strings"
)

// loadBedFile is function that load bed files
func loadBedFile() (map[string]*set.Set, error) {
	targetPosition := make(map[string]*set.Set)
	if conf.BedFile == "" {
		return targetPosition, nil
	}
	sugar.Infof("Loading target positions from file %s (region:%s)", conf.BedFile, region.String())

	stats, err := os.Stat(conf.BedFile)
	if os.IsNotExist(err) {
		return targetPosition, errors.New(conf.BedFile + " not exists")
	}

	f, err := os.Open(conf.BedFile)
	if err != nil {
		return targetPosition, errors.Wrapf(err, "failed to open %s", conf.BedFile)
	}
	defer f.Close()

	bar := pb.New64(stats.Size())
	r := bar.NewProxyReader(f)

	reader := bufio.NewReader(r)
	if strings.HasSuffix(conf.BedFile, ".gz") {
		gr, err := gzip.NewReader(r)
		if err != nil {
			return targetPosition, errors.Wrapf(err, "failed to open %s", conf.BedFile)
		}
		reader = bufio.NewReader(gr)
	}
	totalPositions := 0
	for {
		line, err := reader.ReadString('\n')
		if err != nil {
			if err == io.EOF {
				break
			} else {
				return targetPosition, errors.Wrapf(err, "failed to read from %s", conf.BedFile)
			}
		}

		fields := strings.Split(strings.TrimSpace(line), "\t")

		chr := fields[0]

		if !region.Empty() && strings.ReplaceAll(chr, "chr", "") != strings.ReplaceAll(region.Chrom, "chr", "") {
			continue
		}

		start, err := strconv.Atoi(fields[1])
		if err != nil {
			return targetPosition, errors.Wrapf(err, "failed to convert %s: %s", fields[1], line)
		}

		end, err := strconv.Atoi(fields[2])
		if err != nil {
			end = start
		}

		intersectionStart := 0
		if region.Start > intersectionStart {
			intersectionStart = region.Start
		}
		if start > intersectionStart {
			intersectionStart = start
		}

		intersectionEnd := math.MaxInt64
		if region.End > 0 {
			intersectionEnd = region.End
		}
		if end < intersectionEnd {
			intersectionEnd = end
		}

		if intersectionEnd < intersectionStart {
			continue
		}

		// Add target positions
		temp, ok := targetPosition[chr]
		if !ok {
			temp = &set.Set{}
		}

		j := intersectionStart
		for j <= intersectionEnd {
			temp.Insert(j)
			j++
		}

		targetPosition[chr] = temp
		totalPositions++
	}

	bar.Finish()
	sugar.Infof("### TARGET POSITIONS ###: %d", totalPositions)
	return targetPosition, nil
}

func loadChromosomesFromFai(path string) ([]string, error) {
	res := make([]string, 0, 0)
	f, err := os.Open(path)
	if err != nil {
		return res, errors.Wrapf(err, "failed to open %s", path)
	}

	if r, err := fai.ReadFrom(f); err != nil {
		return res, err
	} else {
		for _, i := range r {
			res = append(res, i.Name)
		}
	}
	return res, nil
}

func createOmopolymericPositions() error {
	if conf.OmopolymericFile == "" {
		return nil
	}
	sugar.Infof("Creating omopolymeric positions (span=%d}) from reference file %s", conf.OmopolymericSpan, conf.Reference)

	index := conf.Reference + ".fai"
	sugar.Infof("Loading chromosome names from index file %s", index)

	chromosomes, err := loadChromosomesFromFai(index)
	if err != nil {
		return err
	}
	sugar.Infof("%d chromosome names found", len(chromosomes))

	positions := make([]*Omopolymeric, 0, 0)

	f, err := os.Open(conf.Reference)
	if err != nil {
		return err
	}

	for _, chromosome := range chromosomes {
		// fasta.Reader requires a known type template to fill
		// with FASTA data. Here we use *linear.Seq.
		template := linear.NewSeq(chromosome, nil, alphabet.DNAredundant)
		r := fasta.NewReader(f, template)

		// Make a seqio.Scanner to simplify iterating over a
		// stream of data.
		sc := seqio.NewScanner(r)

		// Iterate through each sequence in a multifasta and examine the
		// ID, description and sequence data.
		for sc.Next() {
			// Get the current sequence and type assert to *linear.Seq.
			// While this is unnecessary here, it can be useful to have
			// the concrete type.
			s := sc.Seq().(*linear.Seq)

			equals := 0
			var last alphabet.Letter
			for i, b := range s.Seq {
				if b == last {
					equals += 1
				} else {
					if equals >= conf.OmopolymericSpan {
						positions = append(positions, &Omopolymeric{
							Chromosome: chromosome, NegEquals: i - equals,
							I: i, Equals: equals, Last: last,
						})
					}
					equals = 1
				}
				last = b
			}
		}
	}

	sugar.Infof("%s total omopolymeric positions found.", len(positions))

	sugar.Infof("Writing omopolymeric positions to file: %s.", conf.OmopolymericFile)

	f, err = os.OpenFile(conf.OmopolymericFile, os.O_CREATE|os.O_TRUNC|os.O_WRONLY, 0755)
	if err != nil {
		return err
	}
	defer f.Close()

	writer := bufio.NewWriter(f)

	if _, err := writer.WriteString("#" + strings.Join([]string{"Chromomosome", "Start", "End", "Length", "Symbol"}, "\t") + "\n"); err != nil {
		return err
	}

	for _, position := range positions {
		if _, err := writer.WriteString(position.String() + "\n"); err != nil {
			return err
		}
	}
	return nil
}

func loadOmopolymericPositions() (map[string]*set.Set, error) {
	positions := make(map[string]*set.Set)

	if conf.OmopolymericFile == "" {
		return positions, nil
	}

	sugar.Infof("Loading omopolymeric positions from file %s", conf.OmopolymericFile)

	f, err := os.Open(conf.OmopolymericFile)
	if err != nil {
		return positions, errors.Wrapf(err, "failed to open %s", conf.OmopolymericFile)
	}
	defer f.Close()

	r := bufio.NewReader(f)

	total := 0
	for {
		line, err := r.ReadString('\n')
		if err != nil {
			if err == io.EOF {
				break
			} else {
				return positions, err
			}
		}

		if strings.HasPrefix(line, "#") {
			continue
		}

		fields := strings.Split(strings.TrimSpace(line), "\t")
		if region.Chrom == "" || fields[0] == region.Chrom {
			chrom := fields[0]
			f, err := strconv.Atoi(fields[1])
			if err != nil {
				return positions, errors.Wrapf(err, "fields to convert %v", fields)
			}
			t, err := strconv.Atoi(fields[2])
			if err != nil {
				return positions, errors.Wrapf(err, "fields to convert %v", fields)
			}

			if region.Start > 0 && region.Start > f {
				f = region.Start
			}
			if region.End > 0 && region.End < t {
				t = region.End
			}

			temp, ok := positions[chrom]
			if !ok {
				temp = &set.Set{}
			}

			j := f
			for j < t {
				temp.Insert(j)
				total++
				j++
			}
		} else if len(positions) > 0 {
			break
		}
	}

	sugar.Infof("[%v] %d total omopolymeric positions found.", region, total)
	return positions, nil
}

func loadSplicingPositions() (map[string]*set.Set, error) {
	res := make(map[string]*set.Set)
	if conf.SplicingFile == "" {
		return res, nil
	}

	sugar.Infof("Loading known splice sites from file %s", conf.SplicingFile)

	stats, err := os.Stat(conf.SplicingFile)
	if os.IsNotExist(err) {
		return res, errors.New(conf.BedFile + " not exists")
	}

	f, err := os.Open(conf.SplicingFile)
	if err != nil {
		return res, errors.Wrapf(err, "failed to open %s", conf.SplicingFile)
	}
	defer f.Close()

	bar := pb.New64(stats.Size())
	r := bar.NewProxyReader(f)

	reader := bufio.NewReader(r)
	if strings.HasSuffix(conf.SplicingFile, ".gz") {
		gr, err := gzip.NewReader(r)
		if err != nil {
			return res, errors.Wrapf(err, "failed to open %s", conf.SplicingFile)
		}
		reader = bufio.NewReader(gr)
	}

	total := 0

	for {
		line, err := reader.ReadString('\n')
		if err != nil {
			if err == io.EOF {
				break
			} else {
				return res, err
			}
		}

		lines := regexp.MustCompile("\\s+").Split(strings.TrimSpace(line), -1)
		chrom := lines[0]

		temp, ok := res[chrom]
		if !ok {
			temp = &set.Set{}
		}

		st, tp := lines[4], lines[3]
		cc, err := strconv.Atoi(lines[1])
		if err != nil {
			return res, err
		}

		total += conf.SplicingSpan

		if (st == "+" && tp == "D") || (st == "-" && tp == "A") {
			j := 0
			for j < conf.SplicingSpan {
				temp.Insert(cc + (j + 1))
			}
		} else if (st == "+" && tp == "A") || (st == "-" && tp == "D") {
			j := 0
			for j < conf.SplicingSpan {
				temp.Insert(cc - (j + 1))
			}
		}
		res[chrom] = temp
	}

	return res, nil
}

func fetchBamRefs() ([]string, error) {
	res := make([]string, 0, 0)
	ifh, err := os.Open(conf.File)
	//Panic if something went wrong:
	if err != nil {
		return res, err
	}

	//Create a new BAM reader with maximum
	//concurrency:
	bamReader, err := bam.NewReader(ifh, 0)
	if err != nil {
		return res, err
	}

	for _, ref := range bamReader.Header().Refs() {
		res = append(res, ref.Name())
	}
	return res, nil
}

func fetchFasta(region *Region) ([]byte, error) {
	ifh, err := os.Open(conf.Reference)
	if err != nil {
		return nil, err
	}
	defer ifh.Close()

	idxF, err := os.Open(conf.Reference + ".fai")
	if err != nil {
		return nil, err
	}
	defer idxF.Close()

	idx, err := fai.ReadFrom(idxF)

	if err != nil {
		return nil, err
	}

	r := fai.NewFile(ifh, idx)

	sugar.Debugf("load %s from %s", region, conf.Reference)
	seq, err := r.Seq(region.Chrom)
	if err != nil {
		return nil, err
	}
	if region.Start > 0 && region.End > 0 {
		seq, err = r.SeqRange(region.Chrom, region.Start, region.End)
		if err != nil {
			return nil, err
		}
	}

	res := make([]byte, seq.Length, seq.Length)
	_, err = seq.Read(res)
	return res, err
}

func fetchBam(region *Region) (*bam.Iterator, error) {
	sugar.Infof("read reads from %s", region.String())
	ifh, err := os.Open(conf.File)
	//Panic if something went wrong:
	if err != nil {
		return nil, err
	}

	idxF, err := os.Open(conf.File + ".bai")
	if err != nil {
		return nil, err
	}
	defer idxF.Close()

	idx, err := bam.ReadIndex(idxF)
	if err != nil {
		return nil, err
	}

	//Create a new BAM reader with maximum
	//concurrency:
	bamReader, err := bam.NewReader(ifh, 1)
	if err != nil {
		return nil, err
	}

	chunks := make([]bgzf.Chunk, 0, 0)

	for _, ref := range bamReader.Header().Refs() {
		if region.Chrom != "" && ref.Name() == region.Chrom {
			if region.Start == 0 && region.End == 0 {
				if stats, ok := idx.ReferenceStats(ref.ID()); ok {
					chunks = append(chunks, stats.Chunk)
				}
			} else if region.Chrom != "" && region.Start > 0 && region.End > 0 {
				if tempChunks, err := idx.Chunks(ref, region.Start, region.End); err != nil {
					chunks = append(chunks, tempChunks...)
				}
			}
		} else if region.Chrom == "" {
			if stats, ok := idx.ReferenceStats(ref.ID()); ok {
				chunks = append(chunks, stats.Chunk)
			}
		}
	}

	return bam.NewIterator(bamReader, chunks)
}
