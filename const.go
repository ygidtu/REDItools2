package main

import (
	"bytes"
	"fmt"
	"regexp"
	"strconv"
	"strings"
	"time"

	"github.com/biogo/biogo/alphabet"
	"github.com/biogo/hts/bgzf"
	"github.com/biogo/hts/sam"
	"github.com/voxelbrain/goptions"
	"go.uber.org/zap"
)

var logger *zap.Logger
var sugar *zap.SugaredLogger
var conf config
var region *Region

const (
	// VERSION is just the version number
	VERSION = "0.1.3"

	// DefaultBaseQuality as name says
	DefaultBaseQuality = 30
)

type config struct {
	Config                 string  `goptions:"-c, --config, description='Config file'"`
	Version                bool    `goptions:"-v, --version, description='Show version'"`
	Mode                   bool    `goptions:"--fast, description='Usage fast mode with higher memory usage'"`
	Debug                  bool    `goptions:"--debug, description='Show debug info'" long:"debug" description:"Config file"`
	File                   string  `goptions:"-f, --file, description='The bam file to be analyzed'"`
	Output                 string  `goptions:"-o, --output-file, description='The output statistics file'"`
	Process                int     `goptions:"-p, --process, description='How many process to use'"`
	Strict                 bool    `goptions:"-S, --strict, description='Activate strict mode: only sites with edits will be included in the output'"`
	Strand                 int     `goptions:"-s, --strand, description='Strand: this can be 0 (unstranded), 1 (secondstrand oriented) or 2 (firststrand oriented)'"`
	Append                 bool    `goptions:"-a, --append, description='Appends results to file (and creates if not existing)'"`
	Reference              string  `goptions:"-r, --reference, description='The reference FASTA file'"`
	Region                 string  `goptions:"-g, --region, description='The region of the bam file to be analyzed'"`
	OmopolymericFile       string  `goptions:"-m, --omopolymeric-file, description='The file containing the omopolymeric positions'"`
	CreateOmopolymericFile bool    `goptions:"-c, --create-omopolymeric-file, description='Whether to create the omopolymeric span'"`
	OmopolymericSpan       int     `goptions:"--omopolymeric-span, description='The omopolymeric span'"`
	SplicingFile           string  `goptions:"--splicing-file, description='The file containing the splicing sites positions'"`
	SplicingSpan           int     `goptions:"--splicing-span, description='The splicing span'"`
	MinReadLength          int     `goptions:"--min-read-length, description='The minimum read length. Reads whose length is below this value will be discarded.'"`
	MinReadQuality         int     `goptions:"-q, --min-read-quality, description='The minimum read quality. Reads whose mapping quality is below this value will be discarded.'"`
	MinBaseQuality         int     `goptions:"--min-base-quality, description='The minimum base quality. Bases whose quality is below this value will not be included in the analysis.'"`
	MinBasePosition        int     `goptions:"--min-base-position, description='The minimum base position. Bases which reside in a previous position (in the read) will not be included in the analysis.'"`
	MaxBasePosition        int     `goptions:"--max-base-position, description='The maximum base position. Bases which reside in a further position (in the read) will not be included in the analysis.'"`
	MinColumnLength        int     `goptions:"-l, --min-column-length, description='The minimum length of editing column (per position). Positions whose columns have length below this value will not be included in the analysis.'"`
	MinEditsPerNucleotide  int     `goptions:"--min-edits-per-nucleotide, description='The minimum number of editing for events each nucleotide (per position). Positions whose columns have bases with less than min-edits-per-base edits will not be included in the analysis.'"`
	MinEdits               int     `goptions:"--min-edits, description='The minimum number of editing events (per position). Positions whose columns have bases with less than (min-edits-per-base edits) will not be included in the analysis.'"`
	MaxEditsPerNucleotide  int     `goptions:"--max-edits-per-nucleotides, description='The maximum number of editing nucleotides, from 0 to 4 (per position). Positions whose columns have more than (max-editing-nucleotides) will not be included in the analysis.'"`
	StrandConfidence       int     `goptions:"-T, --strand-confidence, description='Strand inference type 1:maxValue 2:useConfidence [1]; maxValue: the most prominent strand count will be used; useConfidence: strand is assigned if over a prefixed frequency confidence (-TV option)'"`
	StrandCorrection       bool    `goptions:"-C, --strand-corection, description='Strand correction. Once the strand has been inferred, only bases according to this strand will be selected.'"`
	StrandConfidenceValue  float64 `goptions:"--strand-confidence-value, description='Strand confidence [0.70]'"`
	RemoveHeader           bool    `goptions:"-H, --remove-header, description='Do not include header in output file'"`
	Dna                    bool    `goptions:"-N, --dna, description='Run REDItools 2.0 on DNA-Seq data'"`
	BedFile                string  `goptions:"-B, --bed_file, description='Path of BED file containing target regions'"`

	Log  string        `goptions:"--log, description='Save log to file'"`
	Help goptions.Help `goptions:"--help, description='Show this help'"`
}

func defaultConfig() config {
	return config{
		Strict: false, Strand: 0, Process: 1,
		OmopolymericSpan: 5, SplicingSpan: 4,
		MinReadLength: 30, MinReadQuality: 20,
		MinBaseQuality: 30, MinBasePosition: 0,
		MaxBasePosition: 0, MinColumnLength: 1,
		MinEditsPerNucleotide: 1, MinEdits: 1,
		MaxEditsPerNucleotide: 100, StrandConfidenceValue: 0.7,
		StrandConfidence: 1, Output: "reditools2_table.gz",
	}
}

// Region is data struct handle the input region,
// format: 1:100-200:+
type Region struct {
	Chrom  string
	Start  int
	End    int
	Strand string
}

// String is function convert region to string format
func (region *Region) String() string {
	if region.Start == 0 && region.End == 0 {
		return region.Chrom
	} else if region.Strand == "" {
		return fmt.Sprintf("%s:%d-%d", region.Chrom, region.Start, region.End)
	}
	return fmt.Sprintf("%s:%d-%d:%s", region.Chrom, region.Start, region.End, region.Strand)
}

// Empty is function that check whether region is empty
func (region *Region) Empty() bool {
	return region.Chrom == ""
}

func decodeRegion(region string) *Region {
	res := &Region{}
	if region == "" {
		return res
	}

	regions := regexp.MustCompile("[:-]").Split(region, -1)

	if regexp.MustCompile("\\w+:\\d+-\\d+(:[+-])?").MatchString(region) {
		res.Chrom = regions[0]
		i, _ := strconv.Atoi(regions[1])
		res.Start = i

		i, _ = strconv.Atoi(regions[2])
		res.End = i
	} else {
		res.Chrom = region
	}

	return res
}

// ChanChunk is struct that transfer chunks between goroutines
type ChanChunk struct {
	Ref    string
	Start  int
	End    int
	Chunks []bgzf.Chunk
}

// Omopolymeric is struct handle the Omopolymeric data
type Omopolymeric struct {
	Chromosome string
	NegEquals  int
	I          int
	Equals     int
	Last       alphabet.Letter
}

func (o *Omopolymeric) String() string {
	return fmt.Sprintf("%s\t%d\t%d\t%d\t%v", o.Chromosome, o.NegEquals, o.I, o.Equals, o.Last)
}

func getHeader() []string {
	return []string{
		"Region", "Position", "Reference",
		"Strand", "Coverage-q30", "MeanQ",
		"BaseCount[A,C,G,T]", "AllSubs",
		"Frequency", "gCoverage-q30", "gMeanQ",
		"gBaseCount[A,C,G,T]", "gAllSubs", "gFrequency",
	}
}

// Record is a wrap of
type Record struct {
	Name        string
	Chrom       string
	Start       int
	End         int
	Ref         string
	Alt         string
	MapQ        byte
	Cigar       sam.Cigar
	Flags       sam.Flags
	Seq         sam.Seq
	Qual        []byte
	Record      *sam.Record
	QueryLength int
}

// NewRecord is function that create a pointer to record
func NewRecord(record *sam.Record) *Record {

	queryLength := 0
	for _, c := range record.Cigar {
		if c.Type() == sam.CigarMatch {
			queryLength += c.Len()
		}
	}

	rec := &Record{
		Name: record.Name,
		MapQ: record.MapQ, Cigar: record.Cigar,
		Flags: record.Flags, Chrom: record.Ref.Name(),
		Seq: record.Seq, Record: record, Qual: record.Qual,
		Start: record.Start(), End: record.End(),
		QueryLength: queryLength,
	}

	return rec
}

// String as name says
func (r *Record) String() string {
	return fmt.Sprintf("%s:%d-%d:%s %s %s", r.Chrom, r.Start, r.End, r.Strand(), r.Ref, r.Alt)
}

// SeqString is function that convert sequence from []byte to string
func (r *Record) SeqString() string {
	return string(r.Seq.Expand())
}

// IsRead1 is true if this is read1
func (r *Record) IsRead1() bool {
	return r.Flags&sam.Read1 != 0
}

// IsRead2 true if this is read2
func (r *Record) IsRead2() bool {
	return r.Flags&sam.Read2 != 0
}

// IsReverse true if this is reversed
func (r *Record) IsReverse() bool {
	return r.Flags&sam.Reverse != 0
}

// Strand is function that calculate the strand based on read1/read2 and reverse
func (r *Record) Strand() string {

	if r.IsRead2() {
		if (conf.Strand == 2 && r.IsReverse()) || (conf.Strand != 2 && !r.IsReverse()) {
			return "-"
		}
		return "+"
	}
	// read1 and single ends
	if (conf.Strand == 1 && r.IsReverse()) || (conf.Strand != 1 && !r.IsReverse()) {
		return "-"
	}
	return "+"
}

// QualityAt is used to return the quality of specific base
func (r *Record) QualityAt(i int) byte {
	if i >= len(r.Qual) {
		sugar.Fatalf("%d -  %d - %v", i, len(r.Qual), r.Qual)
	}
	if i < len(r.Qual) {
		return r.Qual[i]
	}
	return DefaultBaseQuality
}

func complement(sequence byte) byte {
	data := map[byte]byte{'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
	return data[sequence]
}

// SeqAt is used to return base in this position
func (r *Record) SeqAt(i int) byte {
	return r.Seq.At(i)
}

// EditsInfo is struct that keep the information for all reads that aligned to this position
type EditsInfo struct {
	LastChr   string
	Ref       byte
	Pos       int
	Edits     []byte // slice of upper case base from reads in this position
	Counter   map[byte]int
	Variants  map[byte]int
	Strand    string
	Qualities []byte
	Total     int
}

// NewEditsInfo as names says
func NewEditsInfo(lastChr string, ref byte, pos int) *EditsInfo {
	return &EditsInfo{
		LastChr: lastChr, Ref: bytes.ToUpper([]byte{ref})[0], Edits: []byte{}, Pos: pos,
		Counter:  map[byte]int{},
		Variants: map[byte]int{},
		Strand:   "", Qualities: []byte{}, Total: 0,
	}
}

// AddReads as name says add record to this genomic position
func (e *EditsInfo) AddReads(record *Record, at int) {

	if at < conf.MinBasePosition || (conf.MaxBasePosition > 0 && record.QueryLength-at < conf.MaxBasePosition) {
		//sugar.Debugf("failed at base position: %v, min: %v; max: %v", at, conf.MinBasePosition, conf.MaxBasePosition)
		return
	}

	if at >= record.Seq.Length {
		sugar.Fatalf("%s - %d", record.SeqString(), at)
	}

	if record.QualityAt(at) < byte(conf.MinBaseQuality) {
		//sugar.Debugf("failed at base position: %v, min: %v;", record.QualityAt(at), conf.MinBaseQuality)
		return
	}

	seq := record.SeqAt(at)

	if seq != 'N' && e.Ref != 'N' {
		if _, ok := e.Counter[seq]; !ok {
			e.Counter[seq] = 0
		}

		e.Counter[seq]++
		e.Qualities = append(e.Qualities, record.QualityAt(at)) //
		e.Edits = append(e.Edits, bytes.ToUpper([]byte{seq})[0])
		e.Strand += record.Strand()

		if e.Ref != seq {
			if _, ok := e.Variants[seq]; !ok {
				e.Variants[seq] = 0
			}
			e.Variants[seq]++
		}
	}

	e.Total++
}

// MeanQ as name says, return mean quality of add reads
func (e *EditsInfo) MeanQ() float64 {
	res := 0.0
	for _, i := range e.Qualities {
		res += float64(i)
	}

	return res / float64(maxInt(len(e.Qualities), 1))
}

// Valid is function that check whether this edits info is valid
func (e *EditsInfo) Valid() bool {
	if e.Ref == 'N' {
		return false
	}

	if e.MeanQ()-float64(conf.MinReadQuality) < 0 {
		return false
	}

	if len(e.Edits) < conf.MinColumnLength || len(e.Edits) == 0 {
		return false
	}

	numOfEdits := 0
	for _, edit := range e.Variants {
		// sugar.Infof("%s - %s: %d", string(e.Ref), string(seq), edit)

		if edit < conf.MinEditsPerNucleotide || (conf.MaxEditsPerNucleotide > 0 && edit > conf.MaxEditsPerNucleotide) {
			return false
		}
		numOfEdits += edit
	}

	if numOfEdits < conf.MinEdits {
		return false
	}

	return true
}

func (e *EditsInfo) variantsStr() (string, string) {

	seqs := make([]string, 0, 0)
	frequencies := make([]string, 0, 0)

	for key, val := range e.Variants {
		if key != e.Ref && val > 0 {
			seqs = append(seqs, string(key))
			frequencies = append(frequencies, fmt.Sprintf("%.2f", float64(val)/float64(maxInt(e.Total, 1))))
		}
	}

	return strings.Join(seqs, ","), strings.Join(frequencies, ",")
}

// GetStrand as name says, calculate the strand based on all aligned reads
func (e *EditsInfo) GetStrand() string {

	if !conf.Dna {
		strand := vStrand(e.Strand)
		if strand == "-" {
			e.Edits = complementAll(e.Edits)
		}

		if (strand == "+" || strand == "-") && conf.StrandCorrection {
			e.Edits, strand, e.Qualities = normByStrand(e.Edits, e.Strand, e.Qualities, strand)
		}

		return strand

	}
	return "*"
}

func (e *EditsInfo) String() string {

	seqs, frequencies := e.variantsStr()

	return strings.Join([]string{
		e.LastChr, fmt.Sprintf("%d", e.Pos),
		string(e.Ref), e.GetStrand(), fmt.Sprintf("%d", len(e.Edits)),
		fmt.Sprintf("%.2f", e.MeanQ()),
		fmt.Sprintf("[%d, %d, %d, %d]", e.Counter['A'], e.Counter['C'], e.Counter['G'], e.Counter['T']),
		seqs, frequencies, "-", "-", "-", "-", "-",
	}, "\t")
}

//TicTocTimer is structure for timer
type TicTocTimer struct {
	duration time.Duration
	start    time.Time
	repeats  int64
}

//InitTimer is constructor with default values for timer
func InitTimer() *TicTocTimer {
	return &TicTocTimer{duration: 0, start: time.Now(), repeats: 0}
}

// Tic is start timer
func (timer *TicTocTimer) Tic() {
	timer.start = time.Now()
}

//Toc is pause timer
func (timer *TicTocTimer) Toc() {
	timer.duration += time.Now().Sub(timer.start)
	timer.repeats++
}

//TicToc is total time of timer
func (timer *TicTocTimer) TicToc() time.Duration {
	return timer.duration
}
