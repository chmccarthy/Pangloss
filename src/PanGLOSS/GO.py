# -*- coding: utf-8 -*-
"""

"""
import os
import subprocess as sp
from csv import reader
from Tools import Flatten, ParseMatchtable

def MakeWorkingDirs():
    """
    Tries to make work directory if not already present.
    """
    tdir = "go"
    try:
        os.makedirs(tdir)
    except OSError as e:
        if e.errno != os.errno.EEXIST:
            logging.info("GO: Program output directory already exists, using it instead.")
            raise


def GenerateAnnoDict(ips):
    """
    Load InterProScan output (must be tsv format), and generate dictionary of GO annotation data.
    """
    with open(ips) as infile:
        anno_dict = {}
        for line in reader(infile, delimiter="\t"):
            protein = line[0]
            if protein not in anno_dict:
                if len(line) == 14:
                    if line[13]:
                        anno_dict[protein] = [go for go in line[13].split("|") if go]
    return anno_dict


def GenerateAssociations(annos):
    """
    Given GO annotation dictionary, write associations file for use in GOATools.
    """
    with open("go/associations.txt", "w") as assocs:
        for gene in annos:
            if annos[gene]:
                assocs.write("{0}\t{1}\n".format(gene, ";".join(annos[gene])))


def GeneratePopulations(annos, matchtable):
    """
    Write out background (full) population and study (core, accessory) population files for use in GOATools.
    """
    core, acc = ParseMatchtable(matchtable)
    c_pop = [val for val in Flatten(core.values()) if val in annos]
    a_pop = [val for val in Flatten(acc.values()) if val in annos]
    full_pop = c_pop + a_pop
    with open("go/core_pop.txt", "w") as cp_file, open("go/acc_pop.txt", "w") as ap_file,\
         open("go/full_pop.txt", "w") as fp_file:
        cp_file.write("\n".join(c_pop))
        ap_file.write("\n".join(a_pop))
        fp_file.write("\n".join(full_pop))


def GenerateSlimData(assocs, go_obo, slim_obo):
    """

    """
    sp.call(["map_to_slim.py", "--association_file={0}".format(assocs), go_obo, slim_obo],
            stdout=open("go/pangenome_slim_temp.txt", "w"))
    with open("pangenome_slim.txt", "w") as slim:
        for line in open("pangenome_slim_temp.txt").readlines():
            if "|" in line:
                slim.write(line)

def CoreEnrichment():
    """
    """


def AccessoryEnrichment():
    """
    """