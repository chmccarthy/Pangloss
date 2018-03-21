
from Bio import SeqIO
from csv import reader
import multiprocessing as mp
import subprocess as sp
import os


def subprocess_MUSCLE(cmd):
    sp.call(cmd)


def parallel_MUSCLE(muscle_cmds, cores=4):
    farm = mp.Pool(processes=cores)
    farm.map(subprocess_MUSCLE, muscle_cmds)
    farm.close()
    farm.join()


def generate_superalignment(panoct_clusters, fasta_handle):
    db = SeqIO.index(fasta_handle, "fasta")
    matchtable = reader(open(panoct_clusters), delimiter="\t")

    core = {}
    muscle_cmds = []

    for row in matchtable:
        if "----------" in row:
            pass
        else:
            core[row[0]] = row[1:]  # Populating our core dict.

    work_dir = os.getcwd()
    core_dir = "{0}/core_families".format(work_dir)
    if not os.path.isdir(core_dir):
        os.makedirs(core_dir)
    for cluster in core:
        with open("{0}/Cluster_{1}.faa".format(core_dir, cluster), "w") as seqfile:
            for protein in core[cluster]:
                SeqIO.write(db[protein], seqfile, "fasta")
        muscle_cmds.append(["muscle", "-in", "{0}/Cluster_{1}.faa".format(core_dir, cluster), "-out", "{0}/Cluster_{1}.mus".format(core_dir, cluster)])
    parallel_MUSCLE(muscle_cmds)


def main():
    generate_superalignment("matchtable.txt", "panoct_db.fasta")

if __name__ == "__main__":
    main()