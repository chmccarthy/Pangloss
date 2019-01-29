# -*- coding: utf-8 -*-
"""
PanGuess: Gene prediction for PanGLOSS using homology, Hidden Markov Models
          and PVM prediction.

PanGuess is a gene prediction pipeline used to generate protein and genomic
location data for pangenomic analysis of eukaryotes using PanGLOSS. PanGuess
is a component of PanGLOSS and can be used in conjuction with PanGLOSS to determine
pangenomic structure of species of interest based on microsynteny using only genomic data
and a reference protein set.

Requirements:
    - Python (written for 2.7.x)
        - BioPython (1.70)
    - Exonerate (>2.2)
    - GeneMark-ES (>4.30)
    - TransDecoder (>5.0.2)
    - MacOS or Linux system.

Recent changes:
    v0.3.0 (January 2019)
    - Massive rewrite, improved code in a number of ways.
    - ExonerateGene now includes nucleotide sequence data.
    
    v0.2.0 (March 2018)
    - Defined ExonerateGene as class, moved some functions to Tools module.
    - Better integrated codebase with PanGLOSS.
    - Removed other old functions.

    v0.1.0 (Winter 2017)
    - Initial version.

Written by Charley McCarthy, Genome Evolution Lab, Department of Biology,
Maynooth University in 2017-2019 (Charley.McCarthy@nuim.ie).
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
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from Tools import ExonerateCmdLine, LocationOverlap, Pairwise, get_gene_lengths

def LengthOverlap(gene, ref_lengths):
    if gene:
        longest = max(ref_lengths[gene.ref.split("=")[1]], len(gene.called))
        shortest = min(ref_lengths[gene.ref.split("=")[1]], len(gene.called))
        overlap = shortest / longest
        if overlap >= 0.5:
            return True
        else:
            return False
    else:
        return False


def MakeWorkingDir(workdir):
    """
    Tries to make work directory if not already present.
    """
    
    # Don't rewrite work directory if already there.
    try:
        os.makedirs(workdir)
    except OSError as e:
        if e.errno != os.errno.EEXIST:
            raise


def BuildRefSet(workdir, ref):
    """
    Build temporary set of reference proteins. It's faster to run Exonerate by splitting
    up the dataset into individual files and running them as separate queries against
    the genome than as a full file.
    """
    
    # Make folder for reference proteins, if not already present.
    ref_folder = "{0}/ref".format(workdir)
    try:
       os.makedirs(ref_folder)
    except OSError as e:
        if e.errno != os.errno.EEXIST:
            raise
    
    # Split user-provided reference set into individual proteins (have to do this).
    ref_db = SeqIO.index(ref, "fasta")
    for seq in ref_db:
        SeqIO.write(ref_db[seq], "{0}/{1}.faa".format(ref_folder, ref_db[seq].id), "fasta")
    ref_db.close()


def BuildExonerateCmds(workdir, genome):
    """
    Generate list of exonerate commands to run through multiprocessing.
    """
    
    # List of commands.
    exon_cmds = []
    
    # Generate and return commands.
    for prot in glob("{0}/ref/*.faa".format(workdir)):
        exon_cmds.append(["exonerate", "--model", "protein2genome",
                          "-t", genome, "-q", prot, "--bestn", "1"])
    return exon_cmds


def RunExonerate(cmds, len_dict=None, cores=None):
    """
    Farm list of exonerate commands to CPU threads using multiprocessing.
    
    Returns an unordered list of ExonerateGene instances. Default number of
    threads = (number of cores on computer - 1).
    """
    
    # If user doesn't specify cores in command line, just leave them with one free.
    if not cores:
        cores = mp.cpu_count() - 1
    
    # Farm out Exonerate processes, wait for all to finish and merge together.
    farm = mp.Pool(processes=cores)
    genes = farm.map(ExonerateCmdLine, cmds)
    farm.close()
    farm.join()
    
    # Check overlap of predicted gene models against their query homolog.
    if len_dict:
        return [gene for gene in genes if check_overlap(gene, len_dict)]
    else:
        return [gene for gene in genes if gene]


def GetExonerateAttributes(exonerate_genes, tag):
    """
    Extract attributes from ExonerateGene data. Is somewhat redundant, but makes merging
    Exonerate and GeneMark-ES calls easier downstream.
    """
    # Master list of attributes.
    exonerate_attributes = []
    
    # Loop through called genes and extract info. Could be done in one line, obviously.
    for gene in exonerate_genes:
        gene_att = [gene.contig_id, gene.id, gene.locs[0], gene.locs[1]]
        gene_att.append("{0};{1};{2}".format(gene.ref, gene.internal_stop, gene.introns))
        gene_att.append(tag)
        
        # Add to master list.
        exonerate_attributes.append(gene_att)
    
    # Return separate Exonerate attributes object.
    return exonerate_attributes


def RunGeneMark(genome, gm_branch, cores=1):
    """
    Run GeneMark-ES on given genome, with optional arguments for fungal-specific
    prediction models and number of cores.
    """
    
    # If user doesn't specify cores in command line, just leave them with one free.
    if not cores:
        cores = mp.cpu_count() - 1
    
    # Run GeneMark-ES and extract data.
    if gm_branch:
        sp.call(["gmes_petap.pl", "--ES", "--fungus", "--cores", str(cores), "--sequence", genome])
    else:
        sp.call(["gmes_petap.pl", "--ES", "--cores", str(cores), "--sequence", genome])
    sp.call(["get_sequence_from_GTF.pl", "genemark.gtf", genome])
    
    # Return CSV object to convert into attribute data.
    return reader(open("genemark.gtf"), delimiter="\t")


def GeneMarkGTFConverter(gtf, tag):
    """
    Convert a GeneMark-produced GTF/GFF file (which is NOT in a valid GTF/GFF format)
    into an attributes object for easier merging with the predictions from
    exonerate and TransDecoder downstream.
    
    Uses a version of a FSM approach to get the correct locations of genes as
    predicted by GeneMark-ES as well as exon count.
    """
    attributes = []  # List for holding converted attribute info.
    locs = []
    exon_count = 0
    for row, next_row in Pairwise(gtf):  # "gtf" will be a CSV.reader object.
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


def MergeAttributes(tag, first_attributes, second_attributes):
    """
    Return called genes from two methods that do not overlap.
    """
    unique_calls = first_attributes + second_attributes
    unique_calls.sort(key=lambda x: (x[0], int(x[2])))
    to_remove = []
    for call, next_call in Pairwise(unique_calls):
        if next_call:
            overlap = LocationOverlap(call, next_call)
            if overlap:
                to_remove.append(overlap[1])
    return [call for call in unique_calls if call[1] not in to_remove]


def MoveGeneMarkFiles(workdir, genome):
    """
    Handles temporary folders/files created by GeneMark-ES.
    """
    to_move = ["data", "info", "output", "run", "gmes.log", "run.cfg",
               "prot_seq.faa", "nuc_seq.fna", "genemark.gtf"]
    
    gmes = "{0}/gmes/{1}/".format(workdir, genome)
    
    try:
        os.makedirs(gmes)
    except OSError as e:
        if e.errno != os.errno.EEXIST:
            raise
    
    for file in to_move:
        if os.path.isdir(file):
            if not os.path.isdir("{0}/{1}".format(gmes, file)):
                shutil.move(file, gmes)
            else:
                shutil.rmtree(file)
        elif os.path.isfile(file):
            if not os.path.isfile("{0}/{1}".format(gmes, file)):
                shutil.move(file, gmes)
            else:
                os.remove(file)


def ExtractNCR(attributes, genome):
    """
    Generate noncoding sequences from a genome by slicing around known coordinates.
    """
    # List of non-coding regions of genome.
    ncr = []
    
    # Strings for NCR IDs and sequence data.
    extract_id = ""
    extract = ""
    
    # Parse genome file.
    db = SeqIO.parse(open(genome), "fasta")
    
    # Loop over every contig/chromosome in the genome.
    for seq in db:
        coding = filter(lambda x: x[0] == seq.id, attributes)
        for gene, next_gene in Pairwise(coding):
            if coding.index(gene) == 0:
                if gene[2] != 0:
                    extract_id = seq.id + "_NCR_0_{0}".format(gene[2] - 1)
                    extract = seq.seq[0:gene[2] - 2]
                    ncr.append(">{0}\n{1}\n".format(extract_id, extract))
            elif next_gene is None:
                extract_id = seq.id + "_NCR_{0}_{1}".format(gene[3] + 1, len(seq) + 1)
                extract = seq.seq[gene[3]:]
                ncr.append(">{0}\n{1}\n".format(extract_id, extract))
            else:
                extract_id = seq.id + "_NCR_{0}_{1}".format(gene[3] + 1, next_gene[2] - 1)
                extract = seq.seq[gene[3]:next_gene[2] - 2]
                ncr.append(">{0}\n{1}\n".format(extract_id, extract))
    
    return ncr


def RunTransDecoder(ncr, workdir, genome):
    """
    Run the two TransDecoder commands via the command line.
    """
    
    # Path to TransDecoder directory.
    tdir = "{0}/td/{1}/".format(workdir, genome)
    
    # Try to make a directory for TransDecoder. Might as well do it now.
    try:
        os.makedirs(tdir)
    except OSError as e:
        if e.errno != os.errno.EEXIST:
            raise    
    
    # Write NCRs to FASTA file
    with open("{0}/NCR.fna".format(tdir), "w") as outfile:
        for line in ncr:
            outfile.write(line)

    # Run both TransDecoder processes sequentially
    sp.call(["TransDecoder.LongOrfs", "-t", "{0}/NCR.fna".format(tdir), "-m", "200"])
    sp.call(["TransDecoder.Predict", "-t", "{0}/NCR.fna".format(tdir), "--single_best_only"])

    # Return the TransDecoder directory for MoveTransDecoderFiles.
    return tdir


def MoveTransDecoderFiles(tdir):
    """
    Move all temporary TransDecoder files and folders to the TransDecoder directory.
    """
    
    to_move = glob("NCR*") + glob("pipeliner*")
    
    # Move all files named "NCR*".
    for file in to_move:
        if os.path.isdir(file):
            if not os.path.isdir("{0}/{1}".format(tdir, file)):
                shutil.move(file, tdir)
            else:
                shutil.rmtree(file)
        elif os.path.isfile(file):
            if not os.path.isfile("{0}/{1}".format(tdir, file)):
                shutil.move(file, tdir)
            else:
                os.remove(file)


def TransDecoderGTFToAttributes(tdir, tag):
    """
    Extract genomic attributes information from TransDecoder GTF file.
    """
    gtf = reader(open("{0}/NCR.fna.transdecoder.gff3".format(tdir)), delimiter="\t")
    attributes = []  # List for holding converted attribute info.
    locs = []
    exon_count = 0
    for row, next_row in Pairwise(gtf):
        if next_row is not None:
            if row:
                contig_id = row[0].split("_")[0]
                global_locs = map(int, row[0].split("_")[2:])
                if row[2] == "exon":
                      exon_count + 1
                if row[2] == "CDS":
                      relative_locs = map(int, row[3:5])
                      start = global_locs[0] + relative_locs[0] - 1
                      stop = global_locs[0] + relative_locs[1] - 1
                      locs = [start, stop]
                      gene_id = row[-1].split(";")[1].strip("Parent=")
                      annotations = "TransDecoder={0};IS=False;Introns={1}".format(
                                    gene_id, exon_count)
                      
            else:
                attributes.append([contig_id, gene_id, min(locs), max(locs), annotations, tag])
                locs = []
                exon_count = 0
    return sorted(attributes, key=lambda x: (x[0], int(x[2])))


def ConstructGeneModelSets(attributes, exonerate_genes, workdir, genome, tag):
    """
    """
    print workdir
    gm_prot_db = SeqIO.index("{0}/gmes/{1}/prot_seq.faa".format(workdir, genome), "fasta")
    gm_nucl_db = SeqIO.index("{0}/gmes/{1}/nuc_seq.fna".format(workdir, genome), "fasta")
    td_prot_db = SeqIO.index("{0}/td/{1}/NCR.fna.transdecoder.pep".format(workdir, genome), "fasta")
    td_nucl_db = SeqIO.index("{0}/td/{1}/NCR.fna.transdecoder.cds".format(workdir, genome), "fasta")
    
    prot_models = []
    nucl_models = []
    
    for gene in attributes:
        if gene[4].startswith("TransDecoder"):
            print gene[1]
            prot_seq = td_prot_db[gene[1]]
            nucl_seq = td_nucl_db[gene[1]]
            prot_seq.id = "{0}|{1}_{2}_{3}".format(tag, gene[0], gene[2], gene[3])
            nucl_seq.id = prot_seq.id
            gene[1] = prot_seq.id
            prot_models.append(prot_seq)
            nucl_models.append(nucl_seq)
        elif gene[4].startswith("GeneMark"):
            print gene[1]
            prot_seq = gm_prot_db[gene[1]]
            nucl_seq = gm_nucl_db[gene[1]]
            prot_seq.id = "{0}|{1}_{2}_{3}".format(tag, gene[0], gene[2], gene[3])
            nucl_seq.id = prot_seq.id
            gene[1] = prot_seq.id
            prot_models.append(prot_seq)
            nucl_models.append(nucl_seq)
        elif gene[4].startswith("Exonerate"):
            match = filter(lambda x: x.id == gene[1], exonerate_genes)
            prot_seq = SeqRecord(Seq(match[0].prot), id = match[0].id)
            nucl_seq = SeqRecord(Seq(match[0].nucl), id = match[0].id)
            prot_seq.id = "{0}|{1}".format(tag, prot_seq.id)
            nucl_seq.id = "{0}|{1}".format(tag, nucl_seq.id)
            gene[1] = prot_seq.id
            prot_models.append(prot_seq)
            nucl_models.append(nucl_seq)
    
    with open("{0}.faa".format(tag), "w") as outpro:
        SeqIO.write(prot_models, outpro, "fasta")
    
    with open("{0}.nucl".format(tag), "w") as outnuc:
        SeqIO.write(nucl_models, outnuc, "fasta")
        
    with open("{0}.attributes".format(tag), "w") as outatt:
         for line in attributes:
              outatt.write("\t".join(str(el) for el in line) + "\n")
