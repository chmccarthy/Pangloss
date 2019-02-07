

import multiprocessing as mp
import os
import shutil
import subprocess as sp


from Tools import MakeBLASTDBCmdLine


def ConcatenateDatasets(genomes):
    pass


def BLASTAll(genomes, evalue=0.0001, cores=None):
    # If user doesn't specify cores in command line, just leave them with one free.
    logging.info("QualityCheck: Constructing QCBLAST databases using makeblastdb.")
    if not cores:
        cores = mp.cpu_count() - 1