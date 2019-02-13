# -*- coding: utf-8 -*-
import os
import shutil
import subprocess as sp

from glob import glob


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
            logging.info("PanGuess: Working directory already exists, using it instead.")
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
