# -*- coding: utf-8 -*-
import os
import shutil
import subprocess as sp

from Bio import SeqIO

from glob import glob

from Tools import ParseMatchtable


def RunPanOCT(fasta_db, attributes, blast, tags, **kwargs):
    """
    Run PanOCT analysis of gene model dataset. By default, PanGLOSS runs PanOCT with the default parameters
    without specifiying anything.
    """
    cmd = ["/Users/cmccarthy/Documents/GitHub/PanGLOSS/src/panoct.pl", "-t", blast, "-f",
           tags, "-g", attributes, "-P", fasta_db]
    if kwargs:
        pass
    else:
        pass
    sp.call(cmd)


def PanOCTOutputHandler():
    """
    Move expected PanOCT output (might differ from what user actually specifies) to dedicated
    PanOCT output directory.
    """
    to_move = glob("*pairwise*") + glob("*cluster*") + glob("*paralog*") \
              + glob("matchtable*") + ["centroids.fasta", "fragments_fusions.txt", "id.txt",
                                       "missing_blast_results.txt", "parameters.txt", "report.txt"]

    tdir = "panoct"
    try:
        os.makedirs(tdir)
    except OSError as e:
        if e.errno != os.errno.EEXIST:
            logging.info("PanOCT: Program output directory already exists, using it instead.")
            raise

    for f in to_move:
        if os.path.isdir(f):
            if not os.path.isdir("{0}/{1}".format(tdir, f)):
                shutil.move(f, tdir)
            else:
                shutil.rmtree(f)
        elif os.path.isfile(f):
            if not os.path.isfile("{0}/{1}".format(tdir, f)):
                shutil.move(f, tdir)
            else:
                os.remove(f)

def GenerateClusterFASTAs():
    """
    """
    nt_index = SeqIO.index("allnucl.db", "fasta")
    aa_index = SeqIO.index("allprot.db", "fasta")
    fdir = "panoct/clusters"
    matchtable = "panoct/matchtable.txt"
    try:
        os.makedirs(fdir)
        os.makedirs("{0}/core/faa".format(fdir))
        os.makedirs("{0}/core/fna".format(fdir))
        os.makedirs("{0}/acc/faa".format(fdir))
        os.makedirs("{0}/acc/fna".format(fdir))
    except OSError as e:
        if e.errno != os.errno.EEXIST:
            logging.info("PanOCT: Cluster directory already exists, using it instead.")
            raise

    core, acc = ParseMatchtable(matchtable)

    for cluster in core:
        nt_seqs = [nt_index[member] for member in core[cluster]]
        aa_seqs = [aa_index[member] for member in core[cluster]]
        with open("{0}/core/fna/Core_{1}.fna".format(fdir, cluster), "w") as aa_out:
            SeqIO.write(nt_seqs, aa_out, "fasta")

        with open("{0}/core/faa/Core_{1}.faa".format(fdir, cluster), "w") as aa_out:
            SeqIO.write(aa_seqs, aa_out, "fasta")

    for cluster in acc:
        nt_seqs = [nt_index[member] for member in acc[cluster] if member]
        aa_seqs = [aa_index[member] for member in acc[cluster] if member]
        with open("{0}/acc/fna/Acc_{1}.fna".format(fdir, cluster), "w") as aa_out:
            SeqIO.write(nt_seqs, aa_out, "fasta")

        with open("{0}/acc/faa/Acc_{1}.faa".format(fdir, cluster), "w") as aa_out:
            SeqIO.write(aa_seqs, aa_out, "fasta")


