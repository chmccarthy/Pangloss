"""

"""


import os
import shutil
import subprocess as sp
import sys
from csv import reader
from glob import glob

from Bio import SeqIO

from Tools import Flatten, ParseMatchtable, TryMkDirs


def GenerateContigLengths(genomes):
    """
    Parse sequences in original genomes (contigs, chromosomes, &c.) and write their lengths to file.
    """
    lengths = []
    for genome in glob("{0}/*.fna".format(genomes)):
        tag = genome.split("/")[1].split(".")[0]
        gen = []
        for seq in SeqIO.parse(genome, "fasta"):
            gen.append([seq.id, "1", str(len(seq.seq)), tag])
        lengths.append(gen)

    with open("{0}/lengths.txt".format(genomes), "w") as out:
        for gen in lengths:
            for contig in gen:
                row = "\t".join(contig)
                row += "\n"
                out.write(row)

def GenerateKaryotypeFiles(attributes, matchtable):
    """
    Parse concatenated attributes file and PanOCT matchtable, and generate the input needed for Karyotype.R.
    """
    attread = reader(open(attributes), delimiter="\t")
    core, acc = ParseMatchtable(matchtable)
    karyotype = []

    for row in attread:
        karyo = [row[0], row[1], row[2], row[3]]
        if row[1] in Flatten(core.values()):
            karyo = karyo + ["core", row[5]]
        else:
            karyo = karyo + ["acc", row[5]]
            #for key in acc:
            #    if row[1] in acc[key]:
            #        clus_len = len([gene for gene in acc[key] if gene])
            #        karyo.append(str(clus_len))
        karyotype.append(karyo)

    with open("karyotypes.txt", "w") as out:
        for karyo in karyotype:
            line = "\t".join(karyo)
            line += "\n"
            out.write(line)


def KaryoPloteR(tags, karyotypes, lengths):
    """
    Run Karyoplot.R for all strains in a dataset and write the plots to the karyplots folder.
    """
    karyopath = os.path.dirname(os.path.realpath(sys.argv[0])) + "/Karyotype.R"
    sp.call(["Rscript", karyopath, tags, karyotypes, lengths])

    # Don't rewrite work directory if already there.
    kdir = "./karyoplots"
    TryMkDirs(kdir)

    genomes = [line.strip("\n") for line in open(genomelist)]

    for tag in genomes:
        shutil.copy("{0}.eps".format(tag), kdir)
        os.remove("{0}.eps".format(tag))


