# -*- coding: utf-8 -*-
import logging
import os
import shutil
import sys
import subprocess as sp

from Bio import SeqIO, SearchIO
from glob import glob
from Tools import Flatten, ParseMatchtable, QueryClusterFirstHits, Reciprocal


def RunPanOCT(fasta_db, attributes, blast, tags, **kwargs):
    """
    Run PanOCT analysis of gene model dataset. By default, PanGLOSS runs PanOCT with the default parameters
    without specifiying anything.
    """
    panoct_path = os.path.dirname(os.path.realpath(sys.argv[0])) + "/panoct.pl"
    cmd = [panoct_path, "-t", blast, "-f", tags, "-g", attributes, "-P", fasta_db]
    if kwargs:
        pass
    else:
        pass
    logging.info("PanOCT: Running PanOCT on species dataset.")
    sp.call(cmd)


def FillGaps(blast, matchtable, seqs, tags):
    """
    Try to fill in gaps in syntenic clusters that might have arisen via genomic events and/or assembly artefacts.
    """
    # Load core and accessory cluster sets, BLAST+ data and sequence data.
    core, acc = ParseMatchtable(matchtable)
    searches = SearchIO.index(blast, "blast-tab")
    idx = SeqIO.index(seqs, "fasta")
    tags = [line.strip("\n") for line in open(tags)]

    # Get total number of strains in dataset by taking the length of any accessory cluster list (including Nones).
    total = len(acc[acc.keys()[0]])
    n = {}

    # Loop over every accessory cluster.
    og_acc = acc.keys()
    ignore = []
    for q_cluster_id in og_acc:
        current_acc = [key for key in acc.keys() if key not in ignore]
        if q_cluster_id in current_acc:
            q_cluster = acc[q_cluster_id]
            q_present = [gene.split("|")[0] for gene in q_cluster if gene]
            q_missing = filter(lambda tag: tag not in q_present, tags)
            q_blasts = QueryClusterFirstHits(q_cluster, searches, 30, q_missing)
            q_first_hits = set(Flatten(q_blasts.values()))
            merger = []
            if None not in q_first_hits:
                if len(q_first_hits) == len(q_missing):
                    for s_cluster_id in current_acc:
                        if bool(set(acc[s_cluster_id]).intersection(q_first_hits)):
                            s_cluster = acc[s_cluster_id]
                            s_present = [gene.split("|")[0] for gene in s_cluster if gene]
                            s_missing = filter(lambda tag: tag not in s_present, tags)
                            s_blasts = QueryClusterFirstHits(s_cluster, searches, 30, s_missing)
                            s_first_hits = set(Flatten(s_blasts.values()))
                            if Reciprocal(q_cluster, q_first_hits, s_cluster, s_first_hits):
                                if not merger:
                                    merger = [x or y for x, y in zip(q_cluster, s_cluster)]
                                elif len(merger) < total:
                                    intersect = set(s_cluster).intersection(q_first_hits)
                                if len(intersect) == len(q_first_hits):
                                    print q_cluster_id, s_cluster_id, intersect
                                    print str(q_cluster_id) + str(s_cluster_id) + " can be merged into a new core cluster."
                                    ignore = ignore + [q_cluster_id, s_cluster_id]
                                    n[q_cluster_id] = new
                                    del acc[q_cluster_id], acc[s_cluster_id]
                                    print len(n.keys()), len(acc.keys())
                                    break
                                else:
                                    print str(q_cluster_id) + str(s_cluster_id) + " can be merged into a new accessory cluster."
                                    print ""
                    if not intersect:
                        print str(q_cluster_id) + " doesn't have a potential subject cluster in the accessory genome."
                        ignore.append(q_cluster_id)
                        print ""
                else:
                    print str(q_cluster_id) + " doesn't have a possible subject cluster, inconsistent first hits."
                    ignore.append(q_cluster_id)
                    print ""
            else:
                print str(q_cluster_id) + " cannot be combined with another subject cluster."
                ignore.append(q_cluster_id)
                print ""


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
    Extract gene model clusters from full database and write out nucleotide and protein sequence families to file.
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