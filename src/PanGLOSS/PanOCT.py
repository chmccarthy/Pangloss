# -*- coding: utf-8 -*-
import logging
import os
import shutil
import subprocess as sp
import sys
from glob import glob

from Bio import SeqIO, SearchIO

from Tools import Flatten, ParseMatchtable, QueryClusterFirstHits, Reciprocal, TryMkDirs

def RunPanOCT(fasta_db, attributes, blast, genome_list, **kwargs):
    """
    Run PanOCT analysis of gene model dataset. By default, PanGLOSS runs PanOCT with the default parameters
    without specifiying anything.
    """
    panoct_path = os.path.dirname(os.path.realpath(sys.argv[0])) + "/panoct.pl"

    tag_list = [i.strip("\n").split(".")[0].split("/")[1] for i in open(genome_list).readlines()]
    print tag_list
    with open("./panoct_tags.txt", "w") as tag_file:
        tag_file.write("\n".join([str(tag) for tag in tag_list]))

    cmd = [panoct_path, "-t", blast, "-f", "./panoct_tags.txt", "-g", attributes, "-P", fasta_db]
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
    new_mergers = {}

    # Loop over every accessory cluster.
    og_acc = acc.keys()
    ignore = []
    for q_cluster_id in og_acc:
        ignore = []
        current_acc = [key for key in acc.keys() if key not in ignore]
        if q_cluster_id in current_acc:
            q_cluster = acc[q_cluster_id]
            q_pos = [pos for pos, gene in enumerate(q_cluster) if gene]
            q_present = set([tags[pos] for pos in q_pos])
            q_missing = set(filter(lambda tag: tag not in q_present, tags))
            q_blasts = QueryClusterFirstHits(q_cluster, searches, 30, q_missing)
            q_first_hits = set(filter(lambda x: x is not None, Flatten(q_blasts.values())))
            filtered_acc = [key for key in current_acc if not any(acc[key][pos] for pos in q_pos)]
            for s_cluster_id in filtered_acc:
                s_cluster = acc[s_cluster_id]
                s_members = [gene for gene in s_cluster if gene]
                s_present = set([gene.split("|")[0] for gene in s_members if gene])
                s_missing = set(filter(lambda tag: tag not in s_present, tags))
                s_blasts = QueryClusterFirstHits(s_cluster, searches, 30, s_missing)
                s_first_hits = set(filter(lambda x: x is not None, Flatten(s_blasts.values())))
                if len(s_first_hits) == len(q_present):
                    if sorted(q_cluster)
                    print "yes!"

    # Write new matchtable to file.
    with open("refined_matchtable.txt", "w") as out:
        for comp in [core, acc, new_mergers]:
            for cluster in comp:
                line = ["----------" if not a else str(a) for a in comp[cluster]]
                print line
                out.write("\t".join(line) + "\n")


def PanOCTOutputHandler():
    """
    Move expected PanOCT output (might differ from what user actually specifies) to dedicated
    PanOCT output directory.
    """
    to_move = glob("*pairwise*") + glob("*cluster*") + glob("*paralog*") \
              + glob("*matchtable*") + ["centroids.fasta", "fragments_fusions.txt", "id.txt",
                                       "missing_blast_results.txt", "parameters.txt", "report.txt"]

    tdir = "panoct"
    TryMkDirs(tdir)

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


def GenerateClusterFASTAs(matchtable):
    """
    Extract gene model clusters from full database and write out nucleotide and protein sequence families to file.
    """
    nt_index = SeqIO.index("./gm_pred/sets/allnucl.db", "fasta")
    aa_index = SeqIO.index("./gm_pred/sets/allprot.db", "fasta")
    fdir = "./panoct/clusters"
    matchtable = "./panoct/matchtable.txt"
    TryMkDirs(fdir)
    TryMkDirs("{0}/core/faa".format(fdir))
    TryMkDirs("{0}/core/fna".format(fdir))
    TryMkDirs("{0}/acc/faa".format(fdir))
    TryMkDirs("{0}/acc/fna".format(fdir))

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