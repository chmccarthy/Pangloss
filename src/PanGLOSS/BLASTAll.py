
import logging
import multiprocessing as mp
import os
import shutil
import subprocess as sp

from Bio import SeqIO

from Tools import StringBLAST


def ConcatenateDatasets(genomes):
    """
    Concatenate all datasets and construct BLASTp database for gene model set.
    """
    # Generate cat commands for the three full datasets we have.
    tags = [line.strip("\n").split(".")[0].split("/")[1] for line in open(genomes)]
    nucl_cmd = ["cat"] + [tag + ".nucl" for tag in tags]
    prot_cmd = ["cat"] + [tag + ".faa" for tag in tags]
    att_cmd = ["cat"] + [tag + ".attributes" for tag in tags]

    # Run commands.
    logging.info("BLASTAll: Concatenating sequence databases and attributes files.")
    with open("allnucl.db", "w") as f:
        sp.call(nucl_cmd, stdout=f)
    with open("allprot.db", "w") as f:
        sp.call(prot_cmd, stdout=f)
    with open("allatt.db", "w") as f:
        sp.call(att_cmd, stdout=f)

    # Run makeblastdb for gene model set.
    logging.info("BLASTAll: Running makeblastdb command for protein sequence database.")
    sp.call(["makeblastdb", "-in", "allprot.db", "-dbtype", "prot"])


def BLASTAll(evalue=0.0001, cores=None):
    """


    """
    # If user doesn't specify cores in command line, just leave them with one free.
    logging.info("QualityCheck: Constructing QCBLAST databases using makeblastdb.")
    if not cores:
        cores = mp.cpu_count() - 1

    queries = []
    for seq in SeqIO.parse(open("allprot.db"), "fasta"):
        queries.append(">{0}\n{1}".format(seq.id, seq.seq))

    farm = mp.Pool(processes=int(cores))
    results = farm.map(StringBLAST, queries)
    farm.close()
    farm.join()

    return results
