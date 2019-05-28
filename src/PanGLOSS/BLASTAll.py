# -*- coding: utf-8 -*-
"""
BLASTAll: Module for handling parallelized all-vs.-all BLASTp searches, if enabled by user.
"""

import cStringIO
import logging
import multiprocessing as mp
import subprocess as sp

from Bio import SeqIO, SearchIO

from Tools import StringBLAST


def BLASTAll(cores=None):
    """
    Load all query files into memory as strings and BLAST them against all other gene models
    using mp.Pool and the StringBLAST function with n number of cores.
    """
    # If user doesn't specify cores in command line, just leave them with one free.
    if not cores:
        cores = mp.cpu_count() - 1

    logging.info("BLASTAll: Running makeblastdb command for protein sequence database.")
    sp.call(["makeblastdb", "-in", "./gm_pred/sets/allprot.db", "-dbtype", "prot"])

    # Extract FASTA header/sequence strings from database to memory.
    queries = []
    for seq in SeqIO.parse(open("./gm_pred/sets/allprot.db"), "fasta"):
        queries.append(">{0}\n{1}".format(seq.id, seq.seq))

    # Run individual StringBLAST tasks simultaneously.
    logging.info("BLASTAll: Running all-vs.-all BLASTp searches using {0} threads.".format(cores))
    farm = mp.Pool(processes=int(cores))
    results = farm.map(StringBLAST, queries)
    farm.close()
    farm.join()

    # Return the raw strings (could incorporate some filtering here, but I don't think it's necessary).
    logging.info("BLASTAll: All-vs.-all BLAST finished.")
    return results


def MergeBLASTsAndWrite(results):
    """
    Merge all individual BLASTp searches together and write to file in tabular format (without comments this time).
    One thing to note is we have to fool SearchIO into correctly parsing all the results as one big "file" by removing
    the last two lines ("# BLAST processed x queries" &c) from each result object while we're merging everything
    together (see join line).
    """
    # Filter last two lines of each BLASTp result and join remaining lines together, making one big SearchIO object.
    logging.info("BLASTAll: Merging all-vs.-all results together and parsing into tabular format.")
    merged = "\n".join((["\n".join(result.split("\n")[:-2]) for result in results if result]))
    parsed = SearchIO.parse(cStringIO.StringIO(merged), "blast-tab", comments=True)

    # Write merged BLASTp results to file for PanOCT.
    logging.info("BLASTAll: Writing BLASTp results to file panoct.blast.")
    SearchIO.write(parsed, "panoct.blast", "blast-tab")