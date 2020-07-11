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
