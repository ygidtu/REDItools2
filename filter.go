package main

import "github.com/biogo/hts/sam"

func filterReads(record *Record) bool {
	if record.Flags&sam.Unmapped != 0 {
		// sugar.Debugf("Unmapped: %d", record.Flags)
		return false
	}

	if record.Flags&sam.QCFail != 0 {
		// sugar.Debugf("QCFail: %d", record.Flags)
		return false
	}

	if int(record.MapQ) < conf.MinReadQuality {
		// sugar.Debugf("MapQ %d - %d", record.MapQ, conf.MinReadQuality)
		return false
	}

	if record.Flags&sam.Secondary != 0 || record.Flags&sam.Supplementary != 0 {
		return false
	}

	if record.Flags&sam.Duplicate != 0 {
		return false
	}

	if record.Flags&sam.Paired != 0 &&
		(record.Flags != sam.Paired|sam.ProperPair|sam.MateReverse|sam.Read1) && // 99
		(record.Flags != sam.Paired|sam.ProperPair|sam.MateReverse|sam.Read2) && // 163
		(record.Flags != sam.Paired|sam.ProperPair|sam.Reverse|sam.Read2) && // 83
		(record.Flags != sam.Paired|sam.ProperPair|sam.Reverse|sam.Read1) { // 147
		return true
	}

	if _, ok := record.Record.Tag([]byte("SA")); ok {
		return false
	}

	return true
}
