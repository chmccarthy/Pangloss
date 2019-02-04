import cStringIO
import multiprocessing as mp
import os
import shutil
import subprocess as sp

from Bio import SearchIO

from Tools import MakeBLASTDBCmdLine, QCBLASTCmdLine

"""
Search a user-provided set of genes of dubious-quality (i.e. pseudogenes, transposable elements or
transposons &c.) against predicted gene model sets and filter out sufficiently similar genes in the latter.

Arguments:
    gene_sets   = List of strains in analysis (easy access to all files associated with a strain).
    queries     = Set of genes (protein sequences, in fact) to search against all gene model sets.
"""


def BuildMakeBLASTDBs(gene_sets, cores=None):
    """
    Builds BLAST binary database for each strain in a dataset.
    """
    # If user doesn't specify cores in command line, just leave them with one free.
    if not cores:
        cores = mp.cpu_count() - 1

    # Holds makeblastdb commands.
    make_cmds = []

    # Generate commands for every strain.
    for strain in gene_sets:
        cmd = ["makeblastdb", "-in", "{0}.faa".format(strain), "-dbtype", "prot"]
        make_cmds.append(cmd)

    # Run simultaneous makeblastdb commands.
    farm = mp.Pool(processes=cores)
    farm.map(MakeBLASTDBCmdLine, make_cmds)
    farm.close()


def QCBLAST(queries, gene_sets, cores=None):
    """
    BLASTs user-provided proteins against strains datasets and returns a list of parsed SearchIO objects. We parse the
    results of the BLAST(s) as SearchIO objects in the return statement rather than within QCBLASTCmdLine (see Tools.py)
    because mp.Pool can't handle lists of SearchIO objects (something to do with pickling non-standard data types?).
    """
    if not cores:
        cores = mp.cpu_count() - 1

    # Holds BLASTp commands.
    blast_cmds = []

    # Generate commands for every strain.
    for strain in gene_sets:
        cmd = ["blastp", "-query", queries, "-db", "{0}.faa".format(strain), "-outfmt", "5", "-evalue", "0.0001",
               "-num_alignments", "1"]
        blast_cmds.append(cmd)

    # Run simultaneous BLASTp commands.
    farm = mp.Pool(processes=cores)
    blasts = farm.map(QCBLASTCmdLine, blast_cmds)
    farm.close()
    farm.join()

    # Return list of parsed BLASTp results.
    return [SearchIO.parse(cStringIO.StringIO(blast), "blast-xml") for blast in blasts]


