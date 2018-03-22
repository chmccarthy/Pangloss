# -*- coding: utf-8 -*-
"""
PanGuess: Gene prediction for PanGLOSS using homology, Hidden Markov Models
		  and PVM prediction.

PanGuess is a gene prediction pipeline used to generate protein and genomic
location data for pangenomic analysis of eukaryotes using PanGLOSS. PanGuess
is the sister program to PanGLOSS and can be used in conjuction
to determine pangenomic structure of species of interest based on microsynteny
using only genomic data and a reference protein set.

Requirements:
	- Python (written for 2.7.x)
		- BioPython (1.70)
	- Exonerate (>2.2)
	- GeneMark-ES (>4.30)
	- TransDecoder (>5.0.2)
	- MacOS (tested on MacOS 10.12) or Linux (tested on SLES 11)

Recent changes:

	v0.2.0 (March 2018)
	- Defined ExonerateGene as class, moved some functions to Tools module.
	- Removed other old functions.

	v0.1.0 (Winter 2017)
	- Initial version.

To-do:
	- Improve CLI interface: add options, arguments, &c.
	- Replace majority of path/string formatting with dedicated variable to make code easier to read.

Written by Charley McCarthy, Genome Evolution Lab, Department of Biology,
Maynooth University in 2017-2018 (Charley.McCarthy@nuim.ie).
"""

from __future__ import division

import cStringIO
import multiprocessing as mp
import os
import re
import shutil
import subprocess as sp
import sys
import time
from collections import OrderedDict as od
from csv import reader
from glob import glob

from Bio import SearchIO, SeqIO

from PanGLOSS.Tools import pairwise, get_gene_lengths, exoneratecmdline

logfile = open("PanGuess.log", "a", 0)


##### Functions for gene prediction via exonerate. #####
def buildrefset(fasta):
    """
	Given a reference protein FASTA file, build a reference protein bin.

	Splitting an reference set-vs.-genome exonerate search is (probably) much
	faster than searching every reference protein against a genome in one
	command, plus this lets us parallelize searching later on.
	"""
    if not os.path.isdir("{0}/reference_proteins".format(os.getcwd())):
        os.makedirs("{0}/reference_proteins".format(os.getcwd()))
    db = SeqIO.index(fasta, "fasta")
    for seq in db:
        SeqIO.write(db[seq], "reference_proteins/{0}.faa".format(db[seq].id),
                    "fasta")
    db.close()


def buildexontasks(genome, protein_dir):
    """
	Generate list of exonerate commands to run through multiprocessing.
	"""
    exon_cmds = []
    for prot in glob("{0}/*.faa".format(protein_dir)):
        exon_cmds.append(["exonerate", "--model", "protein2genome",
                          "-t", genome, "-q", prot, "--bestn", "1"])
    return exon_cmds


def check_overlap(gene, ref_lengths):
    if gene:
        longest = max(ref_lengths[gene.ref.split("=")[1]], len(gene.called))
        shortest = min(ref_lengths[gene.ref.split("=")[1]], len(gene.called))
        overlap = shortest / longest
        print overlap
        if overlap >= 0.5:
            return True
        else:
            return False
    else:
        return False


def run_exonerate(genome, protein_dir, len_dict=None, cores=None):
    """
	Farm list of exonerate commands to CPU threads using multiprocessing.

	Returns an unordered list of ExonerateGene instances. Default number of
	threads = (number of cores on computer - 1).
	"""
    exon_cmds = buildexontasks(genome, protein_dir)
    if not cores:
        cores = mp.cpu_count() - 1
    farm = mp.Pool(processes=cores)
    genes = farm.map(exoneratecmdline, exon_cmds)
    farm.close()
    farm.join()
    if len_dict:
        return [gene for gene in genes if check_overlap(gene, len_dict)]
    else:
        return [gene for gene in genes if gene]


def run_exonerate_for_transdecoder(genome, protein_dir, cores=None):
    """
	Farm list of exonerate commands to CPU threads using multiprocessing.

	Returns an unordered list of ExonerateGene instances. Default number of
	threads = (number of cores on computer - 1).
	"""
    exon_cmds = buildexontasks(genome, protein_dir)
    if not cores:
        cores = mp.cpu_count() - 1
    farm = mp.Pool(processes=cores)
    genes = farm.map(exoneratecmdline, exon_cmds)
    farm.close()
    farm.join()
    return [gene for gene in genes if gene]


##### Functions for gene prediction using GeneMark-ES. #####
def run_genemark(genome, cores=None):
    """
	Run GeneMark-ES gene prediction on a genome with multithreading.

	Prediction method uses self-training and a specific fungal model. Also runs
	(a customized version of) GeneMark's retrevial script to generate protein
	sequences from GeneMark's own GTF/GFF output. Returns a CSV reader object.
	Default number of threads = (number of cores on computer - 1).
	"""
    if not cores:
        cores = mp.cpu_count() - 1
    sp.call(["gmes_petap.pl", "--ES", "--fungus", "--cores", str(cores), "--sequence", genome])
    sp.call(["get_sequence_from_GTF.pl", "genemark.gtf", genome])
    return reader(open("genemark.gtf"), delimiter="\t")


def genemark_gtf_to_attributes(gtf, tag):
    """
	Convert a GeneMark-produced GTF/GFF file (which is NOT in a valid format)
	into an attributes file for easier merging with the predictions from
	exonerate and TransDecoder downstream.

	Uses a version of a FSM approach to get the correct locations of genes as
	predicted by GeneMark-ES as well as exon count. I'm not a fan of the way
	GeneMark produces GTF/GFF files or how its other Perl scripts parse them.
	"""
    attributes = []  # List for holding converted attribute info.
    locs = []
    exon_count = 0
    for row, next_row in pairwise(gtf):  # "gtf" will be a CSV.reader object.
        if next_row is not None:
            if row[8].split("\"") == next_row[8].split("\""):
                if row[2] == "exon":
                    exon_count = exon_count + 1
                locs = locs + [int(row[3]), int(row[4])]
            else:
                if row[2] == "exon":
                    exon_count = exon_count + 1
                locs = locs + [int(row[3]), int(row[4])]
                contig_id = row[0].split(" ")[0]
                gene_id = row[8].split("\"")[1]
                annotations = "GeneMark={0};IS=False;Introns={1}".format(gene_id, exon_count - 1)
                attributes.append([contig_id, gene_id, min(locs), max(locs),
                                   annotations, tag])
                locs = []
                exon_count = 0
        else:
            if row[2] == "exon":
                exon_count = exon_count + 1
            locs = locs + [int(row[3]), int(row[4])]
            contig_id = row[0].split(" ")[0]
            gene_id = row[8].split("\"")[1]
            annotations = "GeneMark={0};IS=False;Introns={1}".format(gene_id, exon_count - 1)
            attributes.append([contig_id, gene_id, min(locs), max(locs),
                               annotations, tag])
    return sorted(attributes, key=lambda x: (x[0], int(x[2])))


# For locs, we assume lowest value is start, highest is stop.


def genemark_folder_handler(genome, to_move):
    """
	Handles temporary folders/files created by GeneMark-ES.

	For first-time predictions, folders/files are moved to an new folder within
	a given genome's directory called "genemark_output" (which is created here
	if not extant beforehand). For subsequent predictions (i.e. after aborted
	runs), these temporary folders/files are deleted.
	"""
    try:
        os.makedirs("{0}/gene_calling/{1}/genemark_output".format(os.getcwd(), genome))
    except OSError as e:
        if e.errno != os.errno.EEXIST:
            raise
    for f in to_move:
        if os.path.isdir(f):
            if not os.path.isdir("{0}/gene_calling/{1}/genemark_output/{2}".format(os.getcwd(), genome, f)):
                shutil.move(f, "{0}/gene_calling/{1}/genemark_output".format(os.getcwd(), genome))
            else:
                shutil.rmtree(f)
        elif os.path.isfile(f):
            if not os.path.isfile("{0}/gene_calling/{1}/genemark_output/{2}".format(os.getcwd(), genome, f)):
                shutil.move(f, "{0}/gene_calling/{1}/genemark_output".format(os.getcwd(), genome))
            else:
                os.remove(f)


##### Functions for writing exonerate and GeneMark-ES output to files. ######
def write_gene_calls(ordered_exonerate_genes, genemark_genes, tag):
    """
	Write gene calls from exonerate and GeneMark-ES to files.
	"""
    with open("{0}/gene_calling/{1}/{1}_exonerate.txt".format(os.getcwd(), tag), "w") as outfile, open(
            "{0}/gene_calling/{1}/{1}_exonerate.faa".format(os.getcwd(), tag), "w") as outfaa:
        for gene in ordered_exonerate_genes:
            outfile.write("{0}\t{1}|{2}\t{3}\t{4}\t{5}\t{1}\n".format(gene.contig_id, tag,
                                                                      gene.id, gene.locs[0], gene.locs[1],
                                                                      ";".join([gene.ref, str(gene.internal_stop),
                                                                                str(gene.introns)])))
            outfaa.write(">{0}|{1}\n{2}\n".format(tag, gene.id, gene.called))
    with open("{0}/gene_calling/{1}/{1}_genemark.txt".format(os.getcwd(), tag), "w") as outfile:
        genemark_attributes = genemark_gtf_to_attributes(genemark_genes, tag)
        for line in genemark_attributes:
            outfile.write("\t".join([str(column) for column in line]) + "\n")


##### Functions for merging exonerate and GeneMark-ES gene calls. ######
def call_overlap(a, l):
    """
	Check overlapping co-ordinates for calls via exonerate vs. GeneMark-ES.
	"""
    for b in map(lambda x: (int(x[1]), int(x[2])), l):
        if int(a[0]) <= b[0] <= int(a[1]):
            return True
        elif int(a[0]) <= b[1] <= int(a[1]):
            return True
        else:
            for coord in a:
                if b[0] <= int(coord) <= b[1]:
                    return True
                elif b[0] <= int(coord) - 20 <= b[1]:
                    return True
                elif b[0] <= int(coord) + 20 <= b[1]:
                    return True


def get_unique_calls(exonerate_output, genemark_output):
    """
	Return genes called via GeneMark-ES that do not overlap with the
	co-ordinates of genes called via exonerate.
	"""
    unique_calls = []
    exonerate_csv = reader(open(exonerate_output), delimiter="\t")
    exonerate_dict = {}
    for row in exonerate_csv:
        if row[0] not in exonerate_dict.keys():
            exonerate_dict[row[0]] = [row[1:]]
        else:
            exonerate_dict[row[0]].append(row[1:])
    genemark_csv = reader(open(genemark_output), delimiter="\t")
    for row in genemark_csv:  # Safest to do this on a per-chromosome basis.
        if row[0] in exonerate_dict.keys():
            if call_overlap((row[2], row[3]), exonerate_dict[row[0]]):
                pass
            else:
                unique_calls.append(row)
    return unique_calls


def strip_duplicates(gtf):
    """
	Remove gene calls with duplicated locations.

	Tends to effect exonerate calls moreso than GeneMark-ES calls. Necessary
	because PanOCT can't handle duplicate locations (i.e. isoforms).
	"""
    stripped = []
    for line, next_line in pairwise(gtf):
        if next_line is not None:
            if int(line[2]) == int(next_line[2]):
                # If gene has identical start.
                pass
            elif int(line[3]) == int(next_line[3]):
                # If gene has identical end (does this happen?).
                pass
            elif int(next_line[2]) <= int(line[2]) <= int(next_line[3]):
                # If gene starts within next gene.
                pass
            elif int(next_line[2]) <= int(line[3]) <= int(next_line[3]):
                # If gene ends within next gene.
                pass
            elif int(next_line[3]) <= int(line[3]):
                # If gene overlaps with entirety of next gene.
                pass
            else:
                stripped.append(line)
        else:
            if line not in stripped:
                stripped.append(line)
    return stripped


##### Functions for predicting remaining genes using TransDecoder. ######
def get_noncoding_regions(seq_name, seq, list_of_coords):
    """
	Generate noncoding sequences from a genome by slicing around known coordinates.
	"""
    ncr = []
    for coord, next_coord in pairwise(list_of_coords):
        if list_of_coords.index(coord) == 0:  # first gene in a chromosome.
            if coord[0] != 0:  # Check that gene's co-ord isn't 0-to-n!
                ncr.append(">{0}\n{1}\n".format(seq_name + "_NCR_0_{0}".format(coord[0] - 1),
                                                seq.seq[0:coord[0] - 1]))
                if next_coord:  # For single-gene contigs/scaffolds (can happen!).
                    ncr.append(">{0}\n{1}\n".format(seq_name + "_NCR_{0}_{1}".format(coord[1] + 1,
                                                                                     next_coord[0] - 1),
                                                    seq.seq[coord[1] + 1:next_coord[0] - 1]))
        elif next_coord is None:  # (coord) is the last gene in a chromosome.
            ncr.append(">{0}\n{1}\n".format(seq_name + "_NCR_{0}_{1}".format(coord[1] + 1,
                                                                             len(seq)), seq.seq[coord[1] + 1:]))
        else:
            ncr.append(">{0}\n{1}\n".format(seq_name + "_NCR_{0}_{1}".format(coord[1] + 1,
                                                                             next_coord[0] - 1),
                                            seq.seq[coord[1] + 1:next_coord[0] - 1]))
    return ncr


def run_transdecoder(genome, combined_output, tag):
    """
	Predict potential protein-coding ORFs from non-coding regions (NCR) using TransDecoder.

	We generate regions by extracting (per chromosome/contig) subsequences
	not associated with any gene called either by exonerate or GeneMark. These
	regions are then run through TransDecoder, which predicts the longest "ORF"
	per region and then assesses whether it is coding or not.
	"""
    full_genome = SeqIO.index(genome, "fasta")
    combined_csv = reader(open(combined_output), delimiter="\t")
    combined_dict = {}
    noncoding = []
    for row in combined_csv:
        if row[0] not in combined_dict.keys():
            combined_dict[row[0]] = [row[1:]]
        else:
            combined_dict[row[0]].append(row[1:])
    for seq in full_genome:
        if seq in combined_dict.keys():
            known_coords = map(lambda x: (int(x[1]), int(x[2])), combined_dict[seq])
            noncoding = noncoding + get_noncoding_regions(seq, full_genome[seq], known_coords)
    with open("{0}/gene_calling/{1}/{1}_noncoding.fna".format(os.getcwd(), tag), "w") as outncr:
        for line in noncoding:
            outncr.write(line)
    sp.call(["TransDecoder.LongOrfs", "-t", "{0}/gene_calling/{1}/{1}_noncoding.fna".format(os.getcwd(), tag)])
    sp.call(["TransDecoder.Predict", "-t", "{0}/gene_calling/{1}/{1}_noncoding.fna".format(os.getcwd(), tag)])


def transdecoder_gtf_to_attributes(feature_file, tag):
    """
	Convert a TransDecoder-produced GTF/GFF file into an attributes list for
	merging with exonerate and GeneMark-ES attributes.
	"""
    attributes = []  # List for holding converted attribute info.
    exon_count = 0
    gtf = reader(open(feature_file), delimiter="\t")
    for row, next_row in pairwise(gtf):
        if len(row) == 9:
            if row[2] == "gene":
                gene_id = row[8].split(";")[0].split("~")[2]
                contig_id = row[0]
                locs = (row[3], row[4])
            if row[2] == "exon":
                exon_count = exon_count + 1
            if not next_row:
                annotations = "TransDecoder={0};IS=False;Introns={1}".format(gene_id, exon_count - 1)
                attributes.append([contig_id, gene_id, min(locs), max(locs),
                                   annotations, tag])
                exon_count = 0
        else:
            pass
    return sorted(attributes, key=lambda x: (x[0], int(x[2])))


# For locs, we assume lowest value is start, highest is stop.


def split_retained_orfs(tag):
    """
	Split retained ORFs file for realignment using exonerate.

	This way we can realign retained ORFs to a genome using the exact same
	methods we used to search reference homologs against the genome.
	"""
    try:
        os.makedirs("{0}/gene_calling/{1}/temp_retained_orfs".format(os.getcwd(), tag))
    except OSError as e:
        if e.errno != os.errno.EEXIST:
            raise
    for seq in SeqIO.parse("{0}/gene_calling/{1}/transdecoder_output/{1}_retained_orfs.faa".format(os.getcwd(), tag),
                           "fasta"):
        SeqIO.write(seq, open("{0}/gene_calling/{1}/temp_retained_orfs/{2}.faa".format(os.getcwd(), tag, seq.id), "w"),
                    "fasta")


def realign_orfs(genome, tag):
    """
	Re-align putative ORFs back to their genome and get their locations.

	Uses the same parallelized exonerate searches as our homology search does. If
	a putative ORF's top exonerate hit is not within the original ORF's co-ordinates,
	or on a different contig, the ORF is discarded (it's probably poor quality anyway).
	"""
    realigned_orfs = run_exonerate_for_transdecoder(genome, "{0}/gene_calling/{1}/temp_retained_orfs".format(os.getcwd(), tag))
    with open("{0}/gene_calling/{1}/{1}_transdecoder.txt".format(os.getcwd(), tag), "w") as outfile, open(
            "{0}/gene_calling/{1}/{1}_transdecoder.faa".format(os.getcwd(), tag), "w") as outfaa:
        for gene in realigned_orfs:
            original_contig = gene.ref.split("=")[1].split("_")[0]
            if gene.contig_id == original_contig:  # If realigned ORF is on same contig as TransDecoder's call.
                outfile.write("{0}\t{1}|{2}\t{3}\t{4}\t{5}\t{1}\n".format(gene.contig_id, tag,
                                                                          gene.id, gene.locs[0], gene.locs[1],
                                                                          ";".join(["TransDecoder={0}".format(
                                                                              gene.ref.split("=")[1]),
                                                                                    str(gene.internal_stop),
                                                                                    str(gene.introns)])))
                outfaa.write(">{0}|{1}\n{2}\n".format(tag, gene.id, gene.called))


def remove_dubious_orfs(predicted_orfs, dubious_orf_faa):
    """
	If dubious_orfs.faa is present, BLAST against TransDecoder ORFs.

	Such a file is available for Saccharomyces cerevisiae (via SGD), but
	IDK if that's the case for everything (e.g. it isn't for Aspergillus).
	"""
    remove = []
    keep = []
    sp.call(["makeblastdb", "-in", dubious_orf_faa, "-dbtype", "prot"])
    process = sp.check_output(["blastp", "-query", predicted_orfs, "-db",
                               dubious_orf_faa, "-evalue", "0.0001", "-num_alignments", "1"])
    for query in SearchIO.parse(cStringIO.StringIO(process), "blast-text"):
        for subject in query:
            lens = (subject.seq_len, query.seq_len)
            if (min(lens) / max(lens)) >= 0.7:
                remove.append(query.id.split("\n")[0])
    logfile.write("{0} dubious ORFS identified...".format(len(remove)))
    with open("{0}.not_dubious".format(predicted_orfs), "w") as outfile:
        for seq in SeqIO.parse(predicted_orfs, "fasta"):
            if seq.id in remove:
                pass
            else:
                keep.append(seq)
        SeqIO.write(keep, outfile, "fasta")
    # We need to have the seq.descriptions in this file for filter_transdecoder_calls to work!


def filter_transdecoder_calls(tag):
    """
	Align TransDecoder reads back to your genome and filter based on coding score/ORF length.
	"""
    attributes = transdecoder_gtf_to_attributes(
        "{0}/gene_calling/{1}/transdecoder_output/{1}_noncoding.fna.transdecoder.gff3".format(os.getcwd(), tag), tag)
    if os.path.isfile("{0}/gene_calling/{1}/transdecoder_output/{1}_noncoding.fna.transdecoder.pep.not_dubious".format(
            os.getcwd(), tag)):
        orf_index = SeqIO.index(
            "{0}/gene_calling/{1}/transdecoder_output/{1}_noncoding.fna.transdecoder.pep.not_dubious".format(
                os.getcwd(), tag), "fasta")
    else:
        orf_index = SeqIO.index(
            "{0}/gene_calling/{1}/transdecoder_output/{1}_noncoding.fna.transdecoder.pep".format(os.getcwd(), tag),
            "fasta")
    retained_attributes = []
    retained_seqs = []
    for attribute in attributes:
        if attribute[1] in orf_index.keys():
            score_match = re.search("(score=.*) ", orf_index[attribute[1]].description)
            seq_score = float(score_match.group().strip(" ").split("=")[1])
            if all([seq_score >= 100, len(orf_index[attribute[1]]) >= 200]):
                retained_seqs.append(orf_index[attribute[1]])
                retained_attributes.append(attribute)
    with open("{0}/gene_calling/{1}/transdecoder_output/{1}_retained_orfs.txt".format(os.getcwd(), tag), "w") as outatt:
        for line in retained_attributes:
            outatt.write("\t".join(col for col in line) + "\n")
    with open("{0}/gene_calling/{1}/transdecoder_output/{1}_retained_orfs.faa".format(os.getcwd(), tag), "w") as outfaa:
        for seq in retained_seqs:
            outfaa.write(">{0}\n{1}\n".format(seq.id, seq.seq))


def transdecoder_folder_handler(genome):
    """
	Handles temporary folders/files created by TransDecoder

	For first-time predictions, folders/files are moved to an new folder within
	a given genome's directory called "transdecoder_output" (which is created here
	if not extant beforehand). For subsequent predictions (i.e. after aborted
	runs), these temporary folders/files are deleted.
	"""
    try:
        os.makedirs("{0}/gene_calling/{1}/transdecoder_output".format(os.getcwd(), genome))
    except OSError as e:
        if e.errno != os.errno.EEXIST:
            raise
    to_move = glob("*transdecoder*") + glob("pipeliner*")
    for f in to_move:
        if os.path.isdir(f):
            if not os.path.isdir("{0}/gene_calling/{1}/transdecoder_output/{2}".format(os.getcwd(), genome, f)):
                shutil.move(f, "{0}/gene_calling/{1}/transdecoder_output".format(os.getcwd(), genome))
            else:
                shutil.rmtree(f)
        elif os.path.isfile(f):
            if not os.path.isfile("{0}/gene_calling/{1}/transdecoder_output/{2}".format(os.getcwd(), genome, f)):
                shutil.move(f, "{0}/gene_calling/{1}/transdecoder_output".format(os.getcwd(), genome))
            else:
                os.remove(f)


##### Final functions! #####
def remove_duplicates(calls, tag, out_suffix):
    """
	Remove duplicate proteins from exonerate and TransDecoder predictions.

	Basic reason for this is that SeqIO.index can't handle duplicated sequence
	IDs, and PanOCT won't like it either.
	"""
    ids = []
    unique = []
    for seq in SeqIO.parse(calls, "fasta"):
        if seq.id not in ids:
            ids.append(seq.id)
            unique.append(seq)
        else:
            pass
    with open("{0}/gene_calling/{1}/{1}_{2}".format(os.getcwd(), tag, out_suffix), "w") as outfaa:
        SeqIO.write(unique, outfaa, "fasta")


def merge_all_calls(tag):
    """
	Et voila! Kinda messy.

	Note: for TransDecoder calls, because their genomic locations get recalibrated
	in exonerating them back to the genome, the final gene IDs and locations may
	vary slightly from the original IDs and locations as assigned by TransDecoder.
	"""
    exonerate_index = SeqIO.index("{0}/gene_calling/{1}/{1}_exonerate_unique.faa".format(os.getcwd(), tag), "fasta")
    genemark_index = SeqIO.index("{0}/gene_calling/{1}/genemark_output/prot_seq.faa".format(os.getcwd(), tag), "fasta")
    transdecoder_index = SeqIO.index("{0}/gene_calling/{1}/{1}_transdecoder_unique.faa".format(os.getcwd(), tag),
                                     "fasta")
    final_faa = open("{0}/gene_calling/{1}/{1}.faa".format(os.getcwd(), tag), "w")
    final_attributes = open("{0}/gene_calling/{1}/{1}_attributes.txt".format(os.getcwd(), tag), "w")
    all_attributes = [line.strip("\n").split("\t") for line in
                      open("{0}/gene_calling/{1}/{1}_exon_gm.txt".format(os.getcwd(), tag))]
    for line in open("{0}/gene_calling/{1}/{1}_transdecoder.txt".format(os.getcwd(), tag)):
        all_attributes.append(line.strip("\n").split("\t"))
    sorted_attributes = sorted(all_attributes, key=lambda x: (x[0], int(x[2])))
    for line in sorted_attributes:
        if line[4].startswith("Exonerate"):
            seq = line[1]
            final_faa.write(">{0}\n{1}\n".format(exonerate_index[seq].id, exonerate_index[seq].seq))
            final_attributes.write("\t".join(row for row in line) + "\n")
        elif line[4].startswith("GeneMark"):
            seq = line[1]
            new_line = line
            new_line[1] = "{0}|{1}_{2}_{3}".format(tag, line[0], line[2], line[3])
            final_faa.write(">{0}\n{1}\n".format(new_line[1], genemark_index[seq].seq))
            final_attributes.write("\t".join(row for row in new_line) + "\n")
        elif line[4].startswith("TransDecoder"):
            seq = line[1]
            new_line = line
            new_line[1] = "{0}|{1}_{2}_{3}".format(tag, line[0], line[2], line[3])
            final_faa.write(">{0}\n{1}\n".format(new_line[1], transdecoder_index[seq].seq))
            final_attributes.write("\t".join(row for row in new_line) + "\n")


##### Main. #####
def main():
    """
	Main workflow of gene prediction per genome.

		1.  Create master gene_calling folder (outside of per genome loop).
		2.  Generate dictionary list of genes from genome/tag list file (ditto).
		3.  Call genes in genome based on homology to reference genes using exonerate.
		4.  Call genes using self-trained branching HMM analysis with GeneMark-ES.
		5.  Write output from exonerate and GeneMark-ES into PanOCT-compatible format.
		6.  Clean up output from GeneMark-ES (compress it in future maybe?).
		7.  Identify genes called by GeneMark that are non-overlapping relative to exonerate.
		8.  Remove gene calls from this combined prediction that have duplicated locations.
		9.  Call potential ORFs in (still) non-coding regions of genome using TransDecoder.
		10. Filter out dubious ORFs and filter remainder based on sequence length and coding potenital.
		10. Merge unique TransDecoder calls with exonerate/GeneMark-ES calls.
		11. Write unified protein set file and attributes file.
	"""
    try:
        os.makedirs("{0}/gene_calling".format(os.getcwd()))
    except OSError as e:
        if e.errno != os.errno.EEXIST:
            raise
    genomes = od()
    if not os.path.isdir("{0}/reference_proteins".format(os.getcwd())):
        os.makedirs("{0}/reference_proteins".format(os.getcwd()))
        buildrefset(proteins)
    else:
        pass
    ref_lengths = get_gene_lengths(proteins)
    for line in open(genomes_list):
        genomes[line.split("\t")[1].strip("\n")] = line.split("\t")[0]
    for genome in genomes.keys():
        genome_time = time.time()
        logfile.write("Predicting genes for {0}...\n".format(genome))
        genome_tag = genomes[genome]
        try:
            os.makedirs("{0}/gene_calling/{1}".format(os.getcwd(), genome_tag))
        except OSError as e:
            if e.errno != os.errno.EEXIST:
                raise
        logfile.write("Exonerating reference genes against {0}...\n".format(genome))
        exonerate_genes = run_exonerate(genome, "reference_proteins", ref_lengths)
        ordered_exonerate_genes = sorted(exonerate_genes, key=lambda x: (x.contig_id, x.locs[0]))
        logfile.write("Running GeneMark-ES for {0}...\n".format(genome))
        genemark_genes = run_genemark(genome)
        gm_temp_data = ["data", "info", "output", "run", "gmes.log", "run.cfg", "prot_seq.faa", "nuc_seq.fna"]
        genemark_folder_handler(genome_tag, gm_temp_data)
        write_gene_calls(ordered_exonerate_genes, genemark_genes, genome_tag)
        if not os.path.isfile(
                "{0}/gene_calling/{1}/genemark_output/{2}".format(os.getcwd(), genome_tag, "genemark.gtf")):
            shutil.move("genemark.gtf", "{0}/gene_calling/{1}/genemark_output".format(os.getcwd(), genome_tag))
        unique_genes = get_unique_calls("{0}/gene_calling/{1}/{1}_exonerate.txt".format(os.getcwd(), genome_tag),
                                        "{0}/gene_calling/{1}/{1}_genemark.txt".format(os.getcwd(), genome_tag))
        exon_coords = [line.strip("\n").split("\t") for line in
                       open("{0}/gene_calling/{1}/{1}_exonerate.txt".format(os.getcwd(), genome_tag)).readlines()]
        combined_coords = sorted(exon_coords + unique_genes, key=lambda x: (x[0], int(x[2])))
        corrected_coords = strip_duplicates(combined_coords)
        logfile.write("Combined and corrected Exonerate and GeneMark predictions...\n")
        with open("{0}/gene_calling/{1}/{1}_exon_gm.txt".format(os.getcwd(), genome_tag), "w") as outfile:
            for line in corrected_coords:
                outfile.write("\t".join(element for element in line) + "\n")
        run_transdecoder(genome, "{0}/gene_calling/{1}/{1}_exon_gm.txt".format(os.getcwd(), genome_tag), genome_tag)
        transdecoder_folder_handler(genome_tag)
        if os.path.isfile("{0}/dubious_orfs.faa".format(os.getcwd())):
            logfile.write("Checking for dubious ORFs...\n")
            remove_dubious_orfs(
                "{0}/gene_calling/{1}/transdecoder_output/{1}_noncoding.fna.transdecoder.pep".format(os.getcwd(),
                                                                                                     genome_tag),
                "./dubious_orfs.faa")
        filter_transdecoder_calls(genome_tag)
        split_retained_orfs(genome_tag)
        realign_orfs(genome, genome_tag)
        logfile.write("Unifying all calls...\n")
        remove_duplicates("{0}/gene_calling/{1}/{1}_exonerate.faa".format(os.getcwd(), genome_tag), genome_tag,
                          "exonerate_unique.faa")
        remove_duplicates("{0}/gene_calling/{1}/{1}_transdecoder.faa".format(os.getcwd(), genome_tag), genome_tag,
                          "transdecoder_unique.faa")
        merge_all_calls(genome_tag)
        logfile.write(
            "Gene prediction for {0} completed... ({1} seconds)\n".format(genome_tag, time.time() - genome_time))


##### Additional checks. ######
def check_dependencies():
    """
	Check to see that everything is installed or in PATH.
	"""
    deps = ["exonerate", "gmes_petap.pl", "get_sequence_from_GTF.pl",
            "TransDecoder.LongOrfs", "TransDecoder.Predict"]
    for prog in deps:
        try:
            sp.check_output(["which", prog])
        except sp.CalledProcessError as which_exec:
            if which_exec.returncode != 0:
                logfile.write("{0} either not installed or missing from your $PATH! Exiting!\n".format(prog))
                exit(0)


##### User input. #####
if __name__ == "__main__":
    start_time = time.time()
    proteins = sys.argv[1]
    genomes_list = sys.argv[2]
    check_dependencies()
    # Need absolute path for TransDecoder to run NCR realignment to genome!
    main()
    logfile.write(
        "Time taken: {0} seconds. Thank you for choosing PanGuess, the friendly gene prediction software.\n".format(
            time.time() - start_time))
