import os
import shutil
import subprocess as sp

from Tools import TryMkDirs


def RunBUSCO(buscopath, lineagepath, tags):
    """
    Runs BUSCO analysis on every protein set and writes output files to BUSCO folder.
    """
    bdir = "./busco"

    # Don't rewrite work directory if already there.
    TryMkDirs(bdir)

    genomes = [line.strip("\n") for line in open(tags)]

    for genome in genomes:
        pset = "./sets/" + genome + ".faa"
        cmd = ["python", buscopath, "-i", pset,  "-l", lineagepath, "-o", "{0}.busco".format(genome),"-m", "prot"]
        sp.call(cmd)

        shutil.move("run_{0}.busco".format(genome), bdir)