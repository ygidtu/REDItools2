#!/usr/bin/env python

'''
Created on 09 gen 2017

@author: flati
'''

import argparse
import datetime
import gzip
import logging
import os
import re

import sys
from collections import Counter
from collections import defaultdict
from multiprocessing import Pool, Lock

import pysam
from sortedcontainers import SortedSet
from tqdm import tqdm


DATE_FORMAT = "%Y-%m-%d %H:%M"
CIGAR_TAGS = ["M", "I", "D", "N", "S", "H", "P", "=", "X", "B"]

# M	BAM_CMATCH	0
# I	BAM_CINS	1
# D	BAM_CDEL	2
# N	BAM_CREF_SKIP	3
# S	BAM_CSOFT_CLIP	4
# H	BAM_CHARD_CLIP	5
# P	BAM_CPAD	6
# =	BAM_CEQUAL	7
# X	BAM_CDIFF	8
# B	BAM_CBACK	9

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(levelname)s %(message)s', datefmt=DATE_FORMAT)


def delta(t2, t1):
    delta = t2 - t1
    hours, remainder = divmod(delta.seconds, 3600)
    minutes, seconds = divmod(remainder, 60)

    return "%02d:%02d:%02d" % (hours, minutes, seconds)


def print_reads(reads, i):
    total = 0
    for key in reads:
        total += len(reads[key])
        logging.info("E[i=" + str(key) + "][" + str(len(reads[key])) + "] strand=" + str(strand))
        for read in reads[key]:
            index = read["alignment_index"]

            logging.info("\tR:" + str(read["reference"]) + " [r1=" + str(read["object"].is_read1) + ", r2=" + str(
                read["object"].is_read2) + ", reverse=" + str(read["object"].is_reverse) + ", pos=" + str(
                read["pos"]) + ", alignment_index=" + str(index) + ", reference_start=" + str(
                read["object"].reference_start) + " , align_start=" + str(
                read["object"].query_alignment_start) + ", cigar=" + str(read["cigar"]) + ", cigar_list=" + str(
                read["cigar_list"]) + ", " + str(len(read["query_qualities"])) + ", " + str(
                read["query_qualities"]) + "]")
            logging.info("\tQ:" + str(read["sequence"]))
    logging.info("READS[i=" + str(i) + "] = " + str(total))


def update_reads(reads, i):
    logging.debug("UPDATING READS IN POSITION " + str(i))

    pos_based_read_dictionary = {}

    total = 0

    for ending_position in reads:
        for read in reads[ending_position]:

            cigar_list = read["cigar_list"]
            if len(cigar_list) == 0:
                continue

            if read["pos"] >= i:
                continue

            total += 1

            block = cigar_list[0]
            op = block[1]

            if op == "S":

                del cigar_list[0]

                if not cigar_list:
                    block = None
                else:
                    block = cigar_list[0]
                    op = block[1]

            elif op == "N":
                read["pos"] += block[0]
                del cigar_list[0]

                read["ref"] = None
                read["alt"] = None
                read["qual"] = DEFAULT_BASE_QUALITY

                continue

            if block is not None and op == "I":
                n = block[0]

                read["alignment_index"] += n
                read["ref"] = None
                read["alt"] = read["sequence"][read["alignment_index"]]
                del cigar_list[0]

                if not cigar_list:
                    block = None
                else:
                    block = cigar_list[0]
                    op = block[1]

            if block is not None:
                n = block[0]

                # D I M N S
                if op == "M":

                    read["pos"] += 1

                    block[0] -= 1
                    read["reference_index"] += 1
                    read["alignment_index"] += 1

                    read["ref"] = read["reference"][read["reference_index"]]
                    read["alt"] = read["sequence"][read["alignment_index"]]

                    if block[0] == 0:
                        del cigar_list[0]

                elif op == "D":
                    read["pos"] += n
                    read["ref"] = None
                    read["alt"] = None
                    del cigar_list[0]

                if read["query_qualities"] is not None:
                    read["qual"] = read["query_qualities"][read["alignment_index"]]

            p = read["pos"]
            if p not in pos_based_read_dictionary: pos_based_read_dictionary[p] = []
            pos_based_read_dictionary[p].append(read)

    logging.debug("READS UPDATED IN POSITION " + str(i) + ":" + str(total))

    return pos_based_read_dictionary


def get_column(reads, splice_positions, last_chr, omopolymeric_positions, target_positions, i):
    if splice_positions:
        if i in splice_positions[last_chr]:
            logging.debug("[SPLICE_SITE] Discarding position ({}, {}) because in splice site".format(last_chr, i))
            return None

    if omopolymeric_positions:
        if i in omopolymeric_positions[last_chr]:
            logging.debug("[OMOPOLYMERIC] Discarding position ({}, {}) because omopolymeric".format(last_chr, i))
            return None

    if target_positions:
        if (last_chr in target_positions and i not in target_positions[last_chr]) or (
                "chr" + last_chr in target_positions and i not in target_positions["chr" + last_chr]):
            logging.debug(
                "[TARGET POSITIONS] Discarding position ({}, {}) because not in target positions".format(last_chr, i))
            return None

    edits_no = 0
    edits = []
    ref = None

    r1r2distribution = defaultdict(int)

    strand_column = []
    qualities = []
    for key, reads_val in reads.items():
        for read in reads_val:

            pos = read["alignment_index"]

            # Se il carattere e' nelle prime X posizioni della read
            if pos < MIN_BASE_POSITION:
                logging.debug("APPLIED BASE FILTER [MIN_BASE_POSITION]")
                continue

            # Se il carattere e' nelle ultime Y posizioni della read
            if read["length"] - pos < MAX_BASE_POSITION:
                logging.debug("APPLIED BASE FILTER [MAX_BASE_POSITION]")
                continue

            # Se la qualita' e' < Q
            # if read["query_qualities"][read["alignment_index"]] < MIN_BASE_QUALITY:
            if read["qual"] < MIN_BASE_QUALITY:
                logging.debug("APPLIED BASE FILTER [MIN_BASE_QUALITY] {} {} {} {} {}".format(
                    str(read["query_qualities"]), pos,
                    str(read["query_qualities"][pos]), MIN_BASE_QUALITY, read)
                )
                continue

            if read["pos"] != i:
                logging.debug("[OUT_OF_RANGE] SKIPPING READ i=" + str(i) + " but READ=" + str(read["pos"]))
                continue

            logging.debug("GET_COLUMN  Q_NAME=" + str(read["object"].query_name) +
                          " READ1=" + str(read["object"].is_read1) +
                          " REVERSE=" + str(read["object"].is_reverse) +
                          " i=" + str(i) + " READ=" + str(read)
                          )

            if read["ref"] is None:
                logging.debug("[INVALID] SKIPPING READ i=" + str(i) + " BECAUSE REF is None", read)
                continue
            if read["alt"] is None:
                logging.debug("[INVALID] SKIPPING READ i=" + str(i) + " BECAUSE ALT is None", read)
                continue

            ref = read["ref"].upper()
            alt = read["alt"].upper()

            logging.debug("\tBEF={} {}".format(ref, alt))

            if ref == "N" or alt == "N":
                continue

            logging.debug("\tLAT={} {}".format(ref, alt))

            edits.append(alt)

            q = read["qual"]
            qualities.append(q)

            strand_column.append(read["strand"])

            if alt != ref:
                edits_no += 1

            r1r2distribution[
                ("R1" if read["object"].is_read1 else "R2") + ("-REV" if read["object"].is_reverse else "")] += 1

    if not IS_DNA:
        vstrand = 2
        if strand != 0:
            vstrand = vstand(''.join(strand_column))
            if vstrand == "+":
                vstrand = 1
            elif vstrand == "-":
                vstrand = 0
            elif vstrand == "*":
                vstrand = 2

        if vstrand == 0:
            edits = complement_all(edits)
            ref = complement(ref)

        if vstrand in [0, 1] and strand_correction:
            edits, strand_column, qualities, qualities_positions = normByStrand(edits, strand_column, qualities,
                                                                                vstrand)

        logging.debug(vstrand, ''.join(strand_column))
    else:
        vstrand = "*"

    logging.debug(r1r2distribution)

    passed = len(edits)

    counter = defaultdict(int)
    for e in edits:
        counter[e] += 1

    mean_q = 0
    logging.debug("Qualities[i=" + str(i) + "]=" + str(qualities))

    if len(qualities) > 0:
        # mean_q = numpy.mean(qualities)
        mean_q = float(sum(qualities)) / max(len(qualities), 1)

    if len(counter) == 0:
        logging.debug(
            "[EMPTY] Discarding position ({}, {}) because the associated counter is empty".format(last_chr, i))
        return None

    # [A,C,G,T]
    distribution = [counter['A'] if 'A' in counter else 0,
                    counter['C'] if 'C' in counter else 0,
                    counter['G'] if 'G' in counter else 0,
                    counter['T'] if 'T' in counter else 0]
    ref_count = counter[ref] if ref in counter else 0

    non_zero = 0
    for el in counter:
        if el != ref and counter[el] > 0:
            non_zero += 1

    variants = []

    ratio = 0.0

    for el in sorted(counter.items(), key=lambda x: x[1], reverse=True):
        if el[0] == ref:
            continue
        else:
            variants.append(el[0])
            if ratio == 0.0:
                ratio = (float)(el[1]) / (el[1] + ref_count)

    edits_info = {
        "edits": edits,
        "distribution": distribution,
        "mean_quality": mean_q,
        "counter": counter,
        "non_zero": non_zero,
        "edits_no": edits_no,
        "ref": ref,
        "variants": variants,
        "frequency": ratio,
        "passed": passed,
        "strand": vstrand
    }

    # Check that the column passes the filters
    if not filter_column(edits_info, i): return None

    return edits_info;


def normByStrand(seq_, strand_, squal_, mystrand_):
    st = '+'
    if mystrand_ == 0:
        st = '-'
    seq, strand, qual, squal = [], [], [], ''
    for i in range(len(seq_)):
        if strand_[i] == st:
            seq.append(seq_[i])
            strand.append(strand_[i])
            qual.append(squal_[i])
            squal += chr(squal_[i])
    return seq, strand, qual, squal


def get_strand(read):
    global strand

    raw_read = read["object"]

    if (strand == 1 and (
            (raw_read.is_read1 and raw_read.is_reverse) or (raw_read.is_read2 and not raw_read.is_reverse))) or (
            strand == 2 and (
            (raw_read.is_read1 and not raw_read.is_reverse) or (raw_read.is_read2 and raw_read.is_reverse))):
        return "-"

    return "+"


def filter_read(read):
    # Get the flag of the read
    f = read.flag

    # Se la read non e' mappata (FLAG 77 o 141)
    if f == 77 or f == 141:
        logging.debug("APPLIED FILTER [NOT_MAPPED] f={}".format(str(f)))
        return False

    # Se la read non passa i quality controls (FLAG 512)
    if read.is_qcfail:
        logging.debug("APPLIED FILTER [QC_FAIL]")
        return False

    # Se la read ha un MAPQ < di 30
    if read.mapping_quality < MIN_QUALITY:
        logging.debug("APPLIED FILTER [MAPQ] {} MIN={}".format(read.mapping_quality, MIN_QUALITY))
        return False

    # Se la read ha una lunghezza < XX
    if read.query_length < MIN_READ_LENGTH:
        logging.debug("APPLIED FILTER [MIN_READ_LENGTH] {} MIN={}".format(read.query_length, MIN_READ_LENGTH))
        return False

    # Se la read non mappa in modo unico (FLAG 256 o 2048)
    if read.is_secondary or read.is_supplementary:
        logging.debug("APPLIED FILTER [IS_SECONDARY][IS_SUPPLEMENTARY]")
        return False

    # Se la read e' un duplicato di PCR (FLAG 1024)
    if read.is_duplicate:
        logging.debug("APPLIED FILTER [IS_DUPLICATE]")
        return False

    # Se la read e' paired-end ma non mappa in modo proprio (FLAG diversi da 99/147(+-) o 83/163(-+))
    # 99 = 1+2+32+64 = PAIRED+PROPER_PAIR+MREVERSE+READ1 (+-)
    if read.is_paired and not (f == 99 or f == 147 or f == 83 or f == 163):
        logging.debug("APPLIED FILTER [NOT_PROPER]")
        return False

    if read.has_tag('SA'):
        logging.debug("APPLIED FILTER [CHIMERIC_READ]")
        return False

    return True


def filter_base(read):
    pos = read["alignment_index"]

    # Se il carattere e' nelle prime X posizioni della read
    if pos < MIN_BASE_POSITION:
        logging.debug("APPLIED BASE FILTER [MIN_BASE_POSITION]")
        return False

    # Se il carattere e' nelle ultime Y posizioni della read
    if read["length"] - pos < MAX_BASE_POSITION:
        logging.debug("APPLIED BASE FILTER [MAX_BASE_POSITION]")
        return False

    # Se la qualita' e' < Q
    if "qual" not in read:
        logging.debug("APPLIED BASE FILTER [QUAL MISSING] {} {}".format(pos, read))
        return False

    if read["qual"] < MIN_BASE_QUALITY:
        logging.debug("APPLIED BASE FILTER [MIN_BASE_QUALITY] {} {} {} {} {}".format(str(read["query_qualities"]), pos,
                                                                                     str(read["query_qualities"][pos]),
                                                                                     MIN_BASE_QUALITY, read))
        return False

    return True


def filter_column(column, i):
    edits = column["edits"]

    if column["mean_quality"] < MIN_QUALITY:
        logging.debug("DISCARDING COLUMN i={} {} [MIN_MEAN_COLUMN_QUALITY]".format(i, column))
        return False

    # Se il numero di caratteri e' < X
    if len(edits) < MIN_COLUMN_LENGTH:
        logging.debug("DISCARDING COLUMN i={} {} [MIN_COLUMN_LENGTH]".format(i, len(edits)))
        return False

    counter = column["counter"]
    ref = column["ref"]

    # (per ogni variazione) se singolarmente il numero delle basi che supportano la variazione e' < X
    for edit in counter:
        if edit != ref and counter[edit] < MIN_EDITS_SINGLE:
            logging.debug(
                "DISCARDING COLUMN i={} c({})={} [MIN_EDITS_SINGLE] {}".format(i, edit, counter[edit], counter))
            return False

    # Se esistono  multipli cambi rispetto al reference
    if len(counter.keys()) > MAX_CHANGES:
        logging.debug("DISCARDING COLUMN i={} changes={} [MULTIPLE_CHANGES] {}".format(i, len(counter.keys()), column))
        return False

    # Se tutte le sostituzioni sono < Y
    if column["edits_no"] < MIN_EDITS_NO:
        logging.debug("DISCARDING COLUMN i={} {} [MIN_EDITS_NO]".format(i, column["edits_no"]))
        return False

    return True


def load_omopolymeric_positions(positions, input_file, region):
    if input_file is None: return

    logging.info("Loading omopolymeric positions from file {}".format(input_file))

    chromosome = None
    start = None
    end = None

    if region is not None:
        if len(region) >= 1:
            chromosome = region[0]
        if len(region) >= 2:
            start = region[1]
        if len(region) >= 3:
            end = region[2]

    lines_read = 0
    total = 0

    logging.info("Loading omopolymeric positions of {} between {} and {}".format(chromosome, start, end))

    try:
        reader = open(input_file, "r")

        for line in reader:
            if line.startswith("#"):
                continue

            lines_read += 1
            if lines_read % 500000 == 0:
                logging.info("{} lines read.".format(lines_read))

            fields = line.rstrip().split("\t")
            if chromosome is None or fields[0] == chromosome:
                chrom = fields[0]
                f = int(fields[1])
                t = int(fields[2])

                if start is not None: f = max(start, f)
                if end is not None: t = min(t, end)

                if chrom not in positions:
                    positions[chrom] = SortedSet()

                for i in range(f, t):
                    positions[chrom].add(i)
                    total += 1

            elif positions:
                break

        reader.close()
    except IOError as e:
        logging.error("[{}] Omopolymeric positions file not found at {}. Error: {}".format(region, input_file, e))

    logging.info("[{}] {} total omopolymeric positions found.".format(region, total))


def load_chromosome_names(index_file):
    names = []

    with open(index_file, "r") as lines:
        for line in lines:
            names.append(line.split("\t")[0])

    return names


def load_splicing_file(splicing_file):
    splice_positions = {}

    if splicing_file is None: return splice_positions

    logging.info('Loading known splice sites from file {}'.format(splicing_file))

    if splicing_file.endswith("gz"):
        f = gzip.open(splicing_file, "r")

    else:
        f = open(splicing_file, "r")

    total = 0
    total_array = {}

    for i in f:
        l = i.strip().split()
        chrom = l[0]

        if chrom not in splice_positions:
            splice_positions[chrom] = SortedSet()
            total_array[chrom] = 0

        st, tp, cc = l[4], l[3], int(l[1])

        total += SPLICING_SPAN
        total_array[chrom] += SPLICING_SPAN

        if st == '+' and tp == 'D':
            for j in range(SPLICING_SPAN):
                splice_positions[chrom].add(cc + (j + 1))
        if st == '+' and tp == 'A':
            for j in range(SPLICING_SPAN):
                splice_positions[chrom].add(cc - (j + 1))
        if st == '-' and tp == 'D':
            for j in range(SPLICING_SPAN):
                splice_positions[chrom].add(cc - (j + 1))
        if st == '-' and tp == 'A':
            for j in range(SPLICING_SPAN):
                splice_positions[chrom].add(cc + (j + 1))

    f.close()

    logging.info('Loaded {} positions from file {}'.format(total, splicing_file))
    logging.info('\tPartial:{}'.format(total_array))

    return splice_positions


def create_omopolymeric_positions(reference_file, omopolymeric_file):
    tic = datetime.datetime.now()

    logging.info(
        "Creating omopolymeric positions (span={}) from reference file {}".format(OMOPOLYMERIC_SPAN, reference_file))

    index_file = reference_file + ".fai"
    logging.info("Loading chromosome names from index file {}".format(index_file))
    chromosomes = load_chromosome_names(index_file)
    logging.info("{} chromosome names found".format(str(len(chromosomes))))

    positions = []

    try:
        # Opening reference fasta file
        logging.info("Opening reference file {}.".format(reference_file))
        fasta_reader = pysam.FastaFile(reference_file)
        logging.info("Reference file {} opened.".format(reference_file))

        for chromosome in chromosomes:
            logging.info("Loading reference sequence for chromosome {}".format(chromosome))
            sequence = fasta_reader.fetch(chromosome).lower()
            logging.info("Reference sequence for chromosome {} loaded (len: {})".format(chromosome, len(sequence)))

            equals = 0
            last = None
            for i, b in enumerate(sequence):

                if b == last:
                    equals += 1
                else:
                    if equals >= OMOPOLYMERIC_SPAN:
                        positions.append((chromosome, i - equals, i, equals, last))

                    equals = 1

                last = b

        fasta_reader.close()
        logging.info("Reference file {} closed.".format(reference_file))

    except ValueError as e:
        logging.error("Error in reading reference file {}: message={}".format(reference_file, e))
    except IOError:
        logging.error("The reference file {} could not be opened.".format(reference_file))

    logging.info("{} total omopolymeric positions found.".format(len(positions)))

    toc = datetime.datetime.now()
    logging.info("Time to produce all the omopolymeric positions: {}".format(toc - tic))

    logging.info("Writing omopolymeric positions to file: {}.".format(omopolymeric_file))
    writer = open(omopolymeric_file, "w")
    writer.write("#" + "\t".join(["Chromomosome", "Start", "End", "Length", "Symbol"]) + "")
    for position in positions:
        writer.write("\t".join([str(x) for x in position]) + "")
    writer.close()
    logging.info("Omopolymeric positions written into file: {}.".format(omopolymeric_file))


def init(samfile, region):
    logging.info("Opening bamfile within region=" + str(region))

    if not region:
        with pysam.AlignmentFile(samfile) as r:
            return r.references
    if len(region) == 0:
        return [region[0]]
    try:
        return [region]
    except ValueError:
        region[0] = region[0].replace("chr", "")
        return [region]


def within_interval(i, region):
    if region is None or len(region) <= 1:
        return True

    else:
        start = region[1]
        end = region[2]
        return start <= i <= end


def get_header():
    return [
        "Region", "Position", "Reference",
        "Strand", "Coverage-q30", "MeanQ",
        "BaseCount[A,C,G,T]", "AllSubs",
        "Frequency", "gCoverage-q30", "gMeanQ",
        "gBaseCount[A,C,G,T]", "gAllSubs", "gFrequency"
    ]


def load_target_positions(bed_file, region):
    logging.info("Loading target positions from file {} (region:{})".format(bed_file, region))

    target_positions = {}

    extension = os.path.splitext(bed_file)[1]
    if extension == ".gz":
        handler = gzip.open(bed_file, "r")
    else:
        handler = open(bed_file, "r")

    read = 0
    total_positions = 0
    total = Counter()
    with handler as file:
        for line in file:
            read += 1
            fields = line.strip().split("\t")
            chr = fields[0]
            if read % 10000000 == 0:
                logging.info("[{1}] {0} total lines read. Total positions: {2}".format(read, datetime.datetime.now(),
                                                                                       total_positions))

            if region is not None and chr.replace("chr", "") != region[0].replace("chr", ""): continue

            start = int(fields[1]) - 1

            try:
                end = int(fields[2]) - 1
            except:
                end = start  # In case the file has 2 columns only or the third column is not an integer

            intersection_start = max(region[1] if region is not None and len(region) > 1 else 0, start)
            intersection_end = min(region[2] if region is not None and len(region) > 2 else sys.maxint, end)

            # If the target region does not intersect the currently analyzed region
            if intersection_end < intersection_start: continue

            # Add target positions
            if chr not in target_positions: target_positions[chr] = SortedSet()
            for i in range(intersection_start, intersection_end + 1):
                target_positions[chr].add(i)
                total[chr] += 1
                total_positions += 1

    logging.info("### TARGET POSITIONS ###")
    logging.info(total)
    logging.info("TOTAL POSITIONS:", sum(total.values()))

    return target_positions


def is_upstream(first: pysam.AlignedSegment, second: pysam.AlignedSegment) -> bool:
    u"""
    this function is used to determine whether first reads is upstream of second reads.
    Means: two reads do not have any overlap 
    """
    if first is None:
        return True

    if first.reference_name != second.reference_name:
        return first.reference_name < second.reference_name

    return first.reference_end < second.reference_start


def single_process(args):
    ref, samfile, reference_file, strand, region, splice_positions, omopolymeric_positions, target_positions, outfile, = args
    results = []

    logging.info("reading reference: %s" % ref)
    with pysam.FastaFile(reference_file) as reader:
        try:
            chr_ref = reader.fetch(ref)
        except KeyError as err:
            logging.error(err)
            if ref.startswith("chr"):
                chr_ref = reader.fetch(ref.replace("chr", ""))
            else:
                chr_ref = reader.fetch("chr" + ref)

    reads = {}
    with pysam.AlignmentFile(samfile) as reader:
        last_read = []

        if isinstance(ref, str):
            it = reader.fetch(ref)
        else:
            it = reader.fetch(ref[0], ref[1], ref[2])
        for read in it:
            if len(last_read) == 0:
                last_read = [read.reference_start, read.reference_end]

            # Check that the read passes the filters
            if not filter_read(read):
                continue

            ref_pos = [x[1] for x in read.get_aligned_pairs() if x[0] is not None and x[1] is not None]
            ref_seq = ''.join([chr_ref[x] for x in ref_pos]).upper()

            t = "*"

            if not IS_DNA:
                if read.is_read1:
                    if strand == 1:
                        t = '-' if read.is_reverse else '+'
                    else:
                        t = '+' if read.is_reverse else '-'
                elif read.is_read2:
                    if strand == 2:
                        t = '-' if read.is_reverse else '+'
                    else:
                        t = '+' if read.is_reverse else '-'
                else:  # for single ends
                    if strand == 1:
                        t = '-' if read.is_reverse else '+'
                    else:
                        t = '+' if read.is_reverse else '-'

            qualities = read.query_qualities
            if qualities is None:
                qualities = [DEFAULT_BASE_QUALITY for x in range(0, len(ref_seq))]

            item = {
                "pos": read.reference_start - 1,
                "alignment_index": read.query_alignment_start - 1,
                "reference_index": -1,
                "query_alignment_start": read.query_alignment_start,
                "object": read,
                "reference": ref_seq,
                "reference_len": len(ref_seq),
                "sequence": read.query_sequence,
                "chromosome": read.reference_name,
                "query_qualities": qualities,
                "qualities_len": len(qualities),
                "length": read.query_length,
                "cigar": read.cigarstring,
                "strand": t
            }

            cigar_list = [[x[1], CIGAR_TAGS[x[0]]] for x in read.cigartuples]
            item["cigar_list"] = cigar_list

            end_position = read.reference_end  # pos[-1]
            if end_position not in reads:
                reads[end_position] = []

            logging.debug("Adding item=" + str(item))
            reads[end_position].append(item)

            if read.reference_start != last_read[0]:
                results.append([ref, last_read[0], last_read[1] + 1, reads])
                last_read = [read.reference_start, read.reference_end]
                reads = {}
            else:
                last_read[1] = max(last_read[1], read.reference_end)

        if len(reads) > 0:
            results.append([ref, last_read[0], last_read[1] + 1, reads])

    lock.acquire()
    with gzip.open(outfile, "a") if outfile.endswith("gz") else open(outfile, "a") as w:
        for line in results:
            ref, start, end, reads = line
            for j in range(start, end):
                # pos_based_read_dictionary,
                column = get_column(
                    update_reads(reads, j), splice_positions, ref,
                    omopolymeric_positions, target_positions, j
                )
                if column is not None and within_interval(j, region) and not (strict_mode and column["non_zero"] == 0):
                    w.write("\t".join([
                        ref,
                        str(j),
                        column["ref"],
                        str(column["strand"]),
                        str(column["passed"]),
                        "{0:.2f}".format(column["mean_quality"]),
                        str(column["distribution"]),
                        " ".join([column["ref"] + el for el in column["variants"]]) if column["non_zero"] >= 1 else "-",
                        "{0:.2f}".format(column["frequency"]),
                        "\t".join(['-', '-', '-', '-', '-'])
                    ]) + "\n")
    lock.release()


def init_pool(l):
    global lock
    lock = l


def analyze(options: dict, n_jobs: int = 10):
    global DEBUG
    global activate_debug

    logging.info("PYSAM VERSION: {}".format(pysam.__version__))

    bamfile = options["bamfile"]
    region = options["region"]
    reference_file = options["reference"]
    output = options["output"]
    append = options["append"]
    omopolymeric_file = options["omopolymeric_file"]
    splicing_file = options["splicing_file"]
    create_omopolymeric_file = options["create_omopolymeric_file"]
    bed_file = options["bed_file"] if "bed_file" in options else None

    LAUNCH_TIME = datetime.datetime.now()
    logging.info("[" + str(region) + "] START=" + str(LAUNCH_TIME))

    target_positions = {}
    if bed_file is not None:
        target_positions = load_target_positions(bed_file, region)

    omopolymeric_positions = {}
    if create_omopolymeric_file is True:
        if omopolymeric_file is not None:
            create_omopolymeric_positions(reference_file, omopolymeric_file)
        else:
            logging.error(
                "You asked to create the omopolymeric file, but you did not specify any output file. Exiting.")
            return

    load_omopolymeric_positions(omopolymeric_positions, omopolymeric_file, region)

    splice_positions = []

    if splicing_file:
        splice_positions = load_splicing_file(splicing_file)

    # Take the time
    tic = datetime.datetime.now()
    first_tic = tic

    total = 0

    # Open the iterator
    logging.info("Fetching data from bam {}".format(bamfile))
    logging.info("Narrowing REDItools to region {}".format(region))

    if output is not None:
        outputfile = output
    else:
        prefix = os.path.basename(bamfile)
        if region is not None:
            prefix += "_" + '_'.join([str(x) for x in region])
        outputfile = prefix + "_reditools2_table.gz"

    mode = "a" if append else "w"

    with gzip.open(outputfile, mode) if outputfile.endswith("gz") else open(outputfile, mode) as writer:
        if not options["remove_header"]:
            writer.write("\t".join(get_header()) + "")

    cmds = []
    for ref in init(bamfile, region):
        cmds.append([
            ref, bamfile, reference_file, strand, region,
            splice_positions, omopolymeric_positions, target_positions,
            outputfile,
        ])

    lock = Lock()
    with Pool(min(n_jobs, len(cmds)), initializer=init_pool, initargs=(lock,)) as p:
        p.map(single_process, cmds)

    tac = datetime.datetime.now()
    logging.info(str(total) + " total reads read")
    logging.info("END=" + str(tac) + "\t[" + delta(tac, tic) + "]")
    logging.info(
        "FINAL END=" + str(tac) +
        " START=" + str(first_tic) + "\t" + str(region) +
        "\t[TOTAL COMPUTATION=" + delta(tac, first_tic) +
        "] [LAUNCH TIME:" + str(LAUNCH_TIME) +
        "] [TOTAL RUN=" + delta(tac, LAUNCH_TIME) +
        "] [READS=" + str(total) + "]")


complement_map = {"A": "T", "T": "A", "C": "G", "G": "C"}


def complement(b):
    return complement_map[b]


def complement_all(sequence):
    return ''.join([complement_map[l] for l in sequence])


def prop(tot, va):
    try:
        av = float(va) / tot
    except (ZeroDivisionError, AttributeError):
        av = 0.0
    return av


def vstand(strand):  # strand='+-+-+-++++++-+++'

    vv = [(strand.count('+'), '+'), (strand.count('-'), '-'), (strand.count('*'), '*')]
    if vv[0][0] == 0 and vv[1][0] == 0: return '*'
    if use_strand_confidence:
        # flag che indica se usare il criterio 2, altrimenti usa il criterio 1
        totvv = sum([x[0] for x in vv[:2]])
        if prop(totvv, vv[0][0]) >= strand_confidence_value:
            return '+'
            # strand_confidence_value e' il valore soglia, compreso tra 0 e 1, default 0.7
        if prop(totvv, vv[1][0]) >= strand_confidence_value:
            return '-'
        return '*'
    else:
        if vv[0][0] == vv[1][0] and vv[2][0] == 0:
            return '+'
        return max(vv)[1]


def parse_options():
    # Options parsing
    parser = argparse.ArgumentParser(description='REDItools 2.0')
    parser.add_argument('-f', '--file', help='The bam file to be analyzed')
    parser.add_argument('-o', '--output-file', help='The output statistics file')
    parser.add_argument('-S', '--strict', default=False, action='store_true',
                        help='Activate strict mode: only sites with edits will be included in the output')
    parser.add_argument('-s', '--strand', type=int, default=0,
                        help='Strand: this can be 0 (unstranded), 1 (secondstrand oriented) or 2 (firststrand oriented)')
    parser.add_argument('-a', '--append-file', action='store_true',
                        help='Appends results to file (and creates if not existing)')
    parser.add_argument('-r', '--reference', help='The reference FASTA file')
    parser.add_argument('-g', '--region', help='The region of the bam file to be analyzed')
    parser.add_argument('-m', '--omopolymeric-file', help='The file containing the omopolymeric positions')
    parser.add_argument('-c', '--create-omopolymeric-file', default=False,
                        help='Whether to create the omopolymeric span', action='store_true')
    parser.add_argument('-os', '--omopolymeric-span', type=int, default=5, help='The omopolymeric span')
    parser.add_argument('-sf', '--splicing-file', help='The file containing the splicing sites positions')
    parser.add_argument('-ss', '--splicing-span', type=int, default=4, help='The splicing span')
    parser.add_argument('-mrl', '--min-read-length', type=int, default=30,
                        help='The minimum read length. Reads whose length is below this value will be discarded.')
    parser.add_argument('-q', '--min-read-quality', type=int, default=20,
                        help='The minimum read quality. Reads whose mapping quality is below this value will be discarded.')
    parser.add_argument('-bq', '--min-base-quality', type=int, default=30,
                        help='The minimum base quality. Bases whose quality is below this value will not be included in the analysis.')
    parser.add_argument('-mbp', '--min-base-position', type=int, default=0,
                        help='The minimum base position. Bases which reside in a previous position (in the read) will not be included in the analysis.')
    parser.add_argument('-Mbp', '--max-base-position', type=int, default=0,
                        help='The maximum base position. Bases which reside in a further position (in the read) will not be included in the analysis.')
    parser.add_argument('-l', '--min-column-length', type=int, default=1,
                        help='The minimum length of editing column (per position). Positions whose columns have length below this value will not be included in the analysis.')
    parser.add_argument('-men', '--min-edits-per-nucleotide', type=int, default=1,
                        help='The minimum number of editing for events each nucleotide (per position). Positions whose columns have bases with less than min-edits-per-base edits will not be included in the analysis.')
    parser.add_argument('-me', '--min-edits', type=int, default=0,
                        help='The minimum number of editing events (per position). Positions whose columns have bases with less than \'min-edits-per-base edits\' will not be included in the analysis.')
    parser.add_argument('-Men', '--max-editing-nucleotides', type=int, default=100,
                        help='The maximum number of editing nucleotides, from 0 to 4 (per position). Positions whose columns have more than \'max-editing-nucleotides\' will not be included in the analysis.')
    parser.add_argument('-d', '--debug', default=False, help='REDItools is run in DEBUG mode.', action='store_true')
    parser.add_argument('-T', '--strand-confidence', default=1,
                        help='Strand inference type 1:maxValue 2:useConfidence [1]; maxValue: the most prominent strand count will be used; useConfidence: strand is assigned if over a prefixed frequency confidence (-TV option)')
    parser.add_argument('-C', '--strand-correction', default=False,
                        help='Strand correction. Once the strand has been inferred, only bases according to this strand will be selected.',
                        action='store_true')
    parser.add_argument('-Tv', '--strand-confidence-value', type=float, default=0.7, help='Strand confidence [0.70]')
    parser.add_argument('-V', '--verbose', default=False, help='Verbose information in stderr', action='store_true')
    parser.add_argument('-H', '--remove-header', default=False, help='Do not include header in output file',
                        action='store_true')
    parser.add_argument('-N', '--dna', default=False, help='Run REDItools 2.0 on DNA-Seq data', action='store_true')
    parser.add_argument('-B', '--bed_file', help='Path of BED file containing target regions')

    args = parser.parse_known_args()[0]

    global activate_debug
    activate_debug = args.debug

    logging.getLogger().setLevel(logging.DEBUG if args.verbose else logging.INFO)

    bamfile = args.file
    if bamfile is None:
        logging.error("An input baif VERBOSE: m file is mandatory. Please, provide one (-f|--file)")
        exit(1)

    omopolymeric_file = args.omopolymeric_file
    global OMOPOLYMERIC_SPAN
    OMOPOLYMERIC_SPAN = args.omopolymeric_span
    create_omopolymeric_file = args.create_omopolymeric_file

    reference_file = args.reference
    if reference_file is None:
        logging.error("An input reference file is mandatory. Please, provide one (-r|--reference)")
        exit(1)

    output = args.output_file
    append = args.append_file

    global strict_mode
    strict_mode = args.strict

    global strand
    strand = args.strand

    global strand_correction
    strand_correction = args.strand_correction

    global use_strand_confidence
    use_strand_confidence = bool(args.strand_confidence)

    global strand_confidence_value
    strand_confidence_value = float(args.strand_confidence_value)

    splicing_file = args.splicing_file
    global SPLICING_SPAN
    SPLICING_SPAN = args.splicing_span

    global MIN_READ_LENGTH
    MIN_READ_LENGTH = args.min_read_length

    global MIN_QUALITY
    MIN_QUALITY = args.min_read_quality

    global MIN_BASE_QUALITY
    MIN_BASE_QUALITY = args.min_base_quality

    global DEFAULT_BASE_QUALITY
    DEFAULT_BASE_QUALITY = 30

    global MIN_BASE_POSITION
    MIN_BASE_POSITION = args.min_base_position

    global MAX_BASE_POSITION
    MAX_BASE_POSITION = args.max_base_position

    global MIN_COLUMN_LENGTH
    MIN_COLUMN_LENGTH = args.min_column_length

    global MIN_EDITS_SINGLE
    MIN_EDITS_SINGLE = args.min_edits_per_nucleotide

    global MIN_EDITS_NO
    MIN_EDITS_NO = args.min_edits

    global MAX_CHANGES
    MAX_CHANGES = args.max_editing_nucleotides

    global IS_DNA
    IS_DNA = args.dna

    bed_file = args.bed_file

    if IS_DNA and bed_file is None:
        logging.error(
            "When analyzing DNA-Seq files it is mandatory to provide a BED file containing the positions of target regions (-B|--bed_file)")
        exit(1)

    region = None

    if args.region:
        region = re.split("[:-]", args.region)
        if not region or len(region) == 2 or (len(region) == 3 and region[1] == region[2]):
            logging.error(
                "Please provide a region of the form chrom:start-end (with end > start). Region provided: {}".format(
                    region))
            exit(1)
        if len(region) >= 2:
            region[1] = int(region[1])
            region[2] = int(region[2])

    options = {
        "bamfile": bamfile,
        "region": region,
        "reference": reference_file,
        "output": output,
        "append": append,
        "omopolymeric_file": omopolymeric_file,
        "create_omopolymeric_file": create_omopolymeric_file,
        "splicing_file": splicing_file,
        "remove_header": args.remove_header,
        "bed_file": bed_file
    }

    return options


if __name__ == '__main__':
    # Python 2.7.17 Done: 0:02:56.701005
    # Python 3.7.6 Done: 0:03:49.384387
    begin = datetime.datetime.now()
    options = parse_options()

    analyze(options)

    logging.info('Done: {}'.format(datetime.datetime.now() - begin))
