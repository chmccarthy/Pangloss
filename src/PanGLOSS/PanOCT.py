# -*- coding: utf-8 -*-
import os
import subprocess as sp


def RunPanOCT(fasta_db, attributes, blast, tags, **kwargs):
    """
    Run PanOCT analysis of gene model dataset. By default, PanGLOSS runs PanOCT with the default parameters
    without specifiying anything, but these can be changed by
    """
    cmd = ["~/Documents/GitHub/PanGLOSS/src/panoct.pl", "-t", blast, "-f", tags, "-g", attributes, "-P", fasta_db]
    if kwargs:
        pass
    else:
        pass
    sp.call(cmd)


def PanOCTOutputHandler():
    pass