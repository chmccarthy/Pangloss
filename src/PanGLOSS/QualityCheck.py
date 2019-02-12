# -*- coding: utf-8 -*-
"""
Search a user-provided set of genes of dubious-quality (i.e. pseudogenes, transposable elements or
transposons &c.) against predicted gene model sets and filter out sufficiently similar genes in the latter.

Arguments:
    gene_sets   = List of strains in analysis (easy access to all files associated with a strain).
    queries     = Set of genes (protein sequences, in fact) to search against all gene model sets.
"""

from __future__ import division

import cStringIO
import logging
import multiprocessing as mp
import os
import shutil

from Bio import SeqIO, SearchIO
from csv import reader

from Tools import MakeBLASTDBCmdLine, QCBLASTCmdLine


def BuildMakeBLASTDBs(gene_sets, cores=None):
    """
    Builds BLAST binary database for each strain in a dataset.
    """
    # If user doesn't specify cores in command line, just leave them with one free.
    logging.info("QualityCheck: Constructing QCBLAST databases using makeblastdb.")
    if not cores:
        cores = mp.cpu_count() - 1

    # Holds makeblastdb commands.
    make_cmds = []

    # Generate commands for every strain.
    for strain in gene_sets:
        cmd = ["makeblastdb", "-in", "{0}.faa".format(strain), "-dbtype", "prot"]
        make_cmds.append(cmd)

    # Run simultaneous makeblastdb commands.
    logging.info("QualityCheck: Farming makeblastdb tasks to {0} threads.".format(cores))
    farm = mp.Pool(processes=int(cores))
    farm.map(MakeBLASTDBCmdLine, make_cmds)
    farm.close()
    logging.info("QualityCheck: QCBLAST databases constructed.")


def QCBLAST(queries, tags, cores=None):
    """
    BLASTs user-provided proteins against strains datasets and returns a list of parsed SearchIO objects. We parse the
    results of the BLAST(s) as SearchIO objects in the return statement rather than within QCBLASTCmdLine (see Tools.py)
    because mp.Pool can't handle lists of SearchIO objects (something to do with pickling non-standard data types?).
    """
    # If user doesn't specify cores in command line, just leave them with one free.
    logging.info("QualityCheck: BLASTing dubious genes against gene model sets.")
    if not cores:
        cores = mp.cpu_count() - 1

    # Holds BLASTp commands.
    blast_cmds = []

    # Generate commands for every strain.
    for strain in tags:
        cmd = ["blastp", "-query", queries, "-db", "{0}.faa".format(strain), "-outfmt", "5", "-evalue", "0.0001",
               "-num_alignments", "1"]
        blast_cmds.append(cmd)

    # Run simultaneous BLASTp commands.
    logging.info("QualityCheck: Farming BLASTp tasks using {0} threads.".format(cores))
    farm = mp.Pool(processes=int(cores))
    blasts = farm.map(QCBLASTCmdLine, blast_cmds)
    farm.close()
    farm.join()

    # Return list of parsed BLASTp results.
    logging.info("QualityCheck: Parsing QCBLAST results into SeqIO XML objects.")
    return [SearchIO.parse(cStringIO.StringIO(blast), "blast-xml") for blast in blasts]


def RemoveDubiousCalls(results, tags):
    """


    """
    logging.info("QualityCheck: Filtering gene model sets for dubious calls.")
    # Master list for calls to remove.
    to_remove = []

    # Loop through all QCBLAST results, flag top-hits that have >=70% sequence coverage with a dubious gene.
    for result in results:
        for query in result:
            if query.hits:
                query_len = query.seq_len
                subj_len = query.hits[0].seq_len
                ratio = min(query_len, subj_len) / max(query_len, subj_len)
                if ratio >= 0.7:
                    to_remove.append(query.hits[0].id)
                    logging.info("QualityCheck: {0} has >=70% length overlap with {1}, assigning {0} as a"
                                 " dubious call.".format(query.hits[0].id, query.id))

    # Remove flagged calls from nucleotide and protein sets, and genomic attributes file.
    for strain in tags:
        tr_strain = filter(lambda x: x.split("|")[0] == strain, to_remove)
        if tr_strain:
            current_prot = list(SeqIO.parse(open("{0}.faa".format(strain)), "fasta"))
            current_nucl = list(SeqIO.parse(open("{0}.nucl".format(strain)), "fasta"))
            current_att = list(reader(open("{0}.attributes".format(strain)), delimiter="\t"))

            new_prot = filter(lambda x: x.id not in tr_strain, current_prot)
            new_nucl = filter(lambda x: x.id not in tr_strain, current_nucl)
            new_att = filter(lambda x: x[1] not in tr_strain, current_att)

            logging.info("QualityCheck: Removed {0} dubious calls from {1},"
                         " writing remaining calls to new files.".format(len(tr_strain), strain))

            # Write protein sequences to file.
            with open("{0}.faa_new".format(strain), "w") as outpro:
                SeqIO.write(new_prot, outpro, "fasta")

            # Write nucleotide sequences to file.
            with open("{0}.nucl_new".format(strain), "w") as outnuc:
                SeqIO.write(new_nucl, outnuc, "fasta")

            # Write attributes to file.
            with open("{0}.attributes_new".format(strain), "w") as outatt:
                for line in new_att:
                    outatt.write("\t".join(str(el) for el in line) + "\n")

    logging.info("QualityCheck: Completed removal of dubious calls from all datasets.")
