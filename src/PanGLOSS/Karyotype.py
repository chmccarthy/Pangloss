"""

"""


import os
import shutil
import subprocess as sp
import sys
from Bio import SeqIO
from csv import reader
from glob import glob
from Tools import Flatten, ParseMatchtable

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
    Parse concatenated attributes file and PanOCT matchtable, and
    """
    attread = reader(open(attributes), delimiter="\t")
    core, acc = ParseMatchtable(matchtable)
    karyotype = []
    core_len = len(core[core.keys()[0]])

    for row in attread:
        karyo = [row[0], row[2], row[3]]
        if row[1] in Flatten(core.values()):
            karyo = karyo + ["core", row[5]]#, str(core_len)]
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
    :param tags:
    :param karyotypes:
    :param lengths:
    :return:
    """
    karyopath = os.path.dirname(os.path.realpath(sys.argv[0])) + "/Karyotype.R"
    sp.call(["Rscript", karyopath, tags, karyotypes, lengths])

    # Don't rewrite work directory if already there.
    try:
        os.makedirs("karyoplots")
    except OSError as e:
        if e.errno != os.errno.EEXIST:
            logging.info("Karyotype: Plot directory already exists, using it instead.")
            raise

    for f in glob("*.eps"):
        shutil.move(f, "karyoplots")


