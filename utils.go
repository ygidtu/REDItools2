package main

import (
	"sort"
	"strings"

	"github.com/golang-collections/collections/set"
)

func maxInt(i, j int) int {
	if i > j {
		return i
	}
	return j
}

func withinInterval(i int) bool {
	if region.Empty() {
		return true
	}
	return region.Start <= i && i <= region.End
}

//func updateReads(reads map[int][]Record, i int) (map[int][]Record, error) {
//	//timer := InitTimer()
//	//timer.Tic()
//	res := make(map[int][]Record)
//
//	for _, values := range reads {
//		for _, read := range values {
//			if len(read.Cigar) == 0 {
//				continue
//			}
//
//			if read.Pos >= i {
//				continue
//			}
//
//			block := read.Cigar[0]
//			op := block.Type()
//			switch op {
//			case sam.CigarSoftClipped: // S
//				read.Cigar = read.Cigar[1:]
//				if len(read.Cigar) > 0 {
//					block = read.Cigar[0]
//				} else {
//					block = sam.NewCigarOp(sam.CigarMismatch, -1)
//				}
//			case sam.CigarSkipped: // N
//				read.Pos += block.Len()
//				read.Cigar = read.Cigar[1:]
//
//				read.Ref = ""
//				read.Alt = ""
//				continue
//			}
//
//			if block.Len() != -1 && op == sam.CigarInsertion { // I
//				read.QueryAlignmentIndex += block.Len()
//				read.Ref = ""
//				read.Alt = read.SeqAtString(read.QueryAlignmentIndex)
//				read.Cigar = read.Cigar[1:]
//				if len(read.Cigar) > 0 {
//					block = read.Cigar[0]
//				} else {
//					block = sam.NewCigarOp(sam.CigarMismatch, -1)
//				}
//			} else if block.Len() != -1 {
//				switch op {
//				case sam.CigarMatch: // M
//					read.Pos += 1
//					read.ReferenceIndex += 1
//					read.QueryAlignmentIndex += 1
//
//					read.Ref = read.RefAtString(read.ReferenceIndex)
//					read.Alt = read.SeqAtString(read.QueryAlignmentIndex)
//
//					if block.Len() == 1 {
//						read.Cigar = read.Cigar[1:]
//					}
//				case sam.CigarDeletion:  // D
//					read.Pos += block.Len()
//					read.Ref = ""
//					read.Alt = ""
//					read.Cigar = read.Cigar[1:]
//				}
//			}
//
//			if read.Alt != "" {
//				read.Alt = strings.ToUpper(read.Alt)
//			}
//			if read.Ref != "" {
//				read.Ref = strings.ToUpper(read.Ref)
//			}
//
//			temp, ok := res[read.Pos]
//			if !ok {
//				temp = make([]Record, 0, 0)
//			}
//			temp = append(temp, read)
//			res[read.Pos] = temp
//		}
//	}
//	//timer.Toc()
//	//sugar.Debugf("update reads spend: %v", timer.TicToc())
//	return res, nil
//}

func prop(tot, va int) float64 {
	if tot == 0 {
		return 0.0
	}

	return float64(va) / float64(tot)
}

func vStrand(strands string) string {

	data := map[rune]int{'+': 0, '-': 0, '*': 0}

	for _, j := range strands {
		data[j]++
	}

	if conf.StrandConfidence != 0 {
		if prop(data['+']+data['-'], data['+']) >= conf.StrandConfidenceValue {
			return "+"
		} else if prop(data['+']+data['-'], data['-']) >= conf.StrandConfidenceValue {
			return "-"
		} else {
			return "*"
		}
	} else {
		if data['+'] == data['-'] && data['*'] == 0 {
			return "+"
		}
	}

	res, num := "", 0
	for i, j := range data {
		if j > num {
			num = j
			res = string(i)
		}
	}

	return res
}

func complementAll(sequence []byte) []byte {
	data := map[byte]byte{'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

	res := make([]byte, 0, 0)

	for _, i := range sequence {
		res = append(res, data[i])
	}

	return res
}

func normByStrand(seq_ []byte, strand_ string, sequal_ []byte) ([]byte, string, []byte) {

	seq, strand, qual := make([]byte, 0, 0), make([]string, 0, 0), make([]byte, 0, 0)

	for i, j := range seq_ {
		if strand_[i] == j {
			seq = append(seq, j)
			strand = append(strand, string(strand_[i]))
			qual = append(qual, sequal_[i])
		}
	}

	return seq, strings.Join(strand, ""), qual
}

func frequency(data []string) map[string]int {
	res := make(map[string]int)
	for _, j := range data {
		if temp, ok := res[j]; ok {
			temp++
			res[j] = temp
		} else {
			res[j] = 1
		}
	}
	return res
}

func getColumn(edits map[int]*EditsInfo, positions []map[string]*set.Set, targetPositions map[string]*set.Set, w chan string) {

	timer := InitTimer()
	timer.Tic()

	keys := make([]int, 0, 0)
	for key := range edits {
		keys = append(keys, key)
	}
	sort.Ints(keys)

	for _, key := range keys {
		edit, ok := edits[key]
		if !ok {
			continue
		}

		// if key == 1115863 {
		// 	sugar.Info(edit)
		// }

		if edit.Ref == byte('N') {
			continue
		}

		for _, pos := range positions {
			if temp, ok := pos[edit.LastChr]; ok {
				if temp.Has(edit.Pos) {
					continue
				}
			}
		}

		if temp, ok := targetPositions[edit.LastChr]; ok {
			if !temp.Has(edit.Pos) {
				continue
			}
		}

		if edit.Valid() {
			w <- edit.String()
		}
	}
}
