"""
Simple module for running UpSet analysis of accessory genomes using UpSetR.
"""
import os
import subprocess as sp
import sys


def UpSetR(tags, matchtable):
    """
    Parse matchtable and generate UpSetR plot of distribution of accessory orthologs within pangenome.
    """
    upsetpath = os.path.dirname(os.path.realpath(sys.argv[0])) + "/UpSet.R"
    sp.call(["Rscript", upsetpath, matchtable, tags])
