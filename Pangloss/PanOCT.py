# -*- coding: utf-8 -*-
import logging
import os
import shutil
import subprocess as sp
import sys
from glob import glob

from Bio import SeqIO, SearchIO

from .Tools import ConcatenateDatasets, ClusterMerge, Flatten, MultipleInsert, ParseMatchtable, \
                  QueryClusterFirstHits, Reciprocal, TryMkDirs

def RunPanOCT(fasta_db, attributes, blast, genome_list, **kwargs):
    """
    Run PanOCT analysis of gene model dataset. By default, Pangloss runs PanOCT with the default parameters
    without specifiying anything.
    """
    panoct_path = os.path.dirname(os.path.realpath(sys.argv[0])) + "/panoct.pl"
    tag_list = []
    for genome in open(genome_list).readlines():
        if "/" in genome:
            tag_list.append(genome.split(".")[0].split("/")[1])
        else:
            tag_list.append(genome.split(".")[0])
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
    new_clusters = {}
    if acc:
        searches = SearchIO.index(blast, "blast-tab")
        tags = [line.strip("\n") for line in open(tags)]

        # Loop over every accessory cluster.
        og_acc = list(acc.keys())
        ignore = []
        for q_cluster_id in og_acc:
            print("{0} out of {1} clusters searched".format(og_acc.index(q_cluster_id), len(og_acc)))
            if q_cluster_id not in ignore:
                current_acc = [key for key in list(acc.keys()) if key not in ignore]
                if q_cluster_id in current_acc:
                    q_cluster = acc[q_cluster_id]
                    q_pos = [pos for pos, gene in enumerate(q_cluster) if gene]
                    q_present = set([tags[pos] for pos in q_pos])
                    q_members = set(sorted([x for x in q_cluster if x is not None]))
                    q_missing = set([tag for tag in tags if tag not in q_present])
                    q_blasts = QueryClusterFirstHits(q_cluster, searches, 30, q_missing)
                    q_first_hits = set([x for x in Flatten(list(q_blasts.values())) if x is not None])
                    q_query = MultipleInsert(list(q_first_hits), tags)
                    if q_query in list(acc.values()):
                        s_cluster_id = list(acc.keys())[list(acc.values()).index(q_query)]
                        if s_cluster_id not in ignore:
                            s_cluster = acc[s_cluster_id]
                            s_members = set(sorted([x for x in s_cluster if x is not None]))
                            if s_members == q_first_hits:
                                s_present = set([gene.split("|")[0] for gene in s_members])
                                s_missing = set([tag for tag in tags if tag not in s_present])
                                s_blasts = QueryClusterFirstHits(s_cluster, searches, 30, s_missing)
                                s_first_hits = set([x for x in Flatten(list(s_blasts.values())) if x is not None])
                                reciprocal = Reciprocal(q_members, q_first_hits, s_members, s_first_hits)
                                if reciprocal:
                                    new_cluster = ClusterMerge(q_cluster, s_cluster)
                                    new_clusters[q_cluster_id] = new_cluster
                                    acc.pop(q_cluster_id, "None")
                                    acc.pop(s_cluster_id, "None")
                                    print("clusters merged: {0} {1}\n".format(str(q_cluster_id), str(s_cluster_id)))
                                    print("size of clusters merged: {0} {1}\n".format(len(q_members), len(s_members)))
                                    ignore = ignore + [q_cluster_id, s_cluster_id]
            else:
                pass
    else:
        pass

    # Write new matchtable to file.
    with open("refined_matchtable.txt", "w") as out:
        if acc:
            dataset = [core, acc, new_clusters]
        else:
            dataset = [core]
        for comp in dataset:
            for cluster in comp:
                line = ["----------" if not a else str(a) for a in comp[cluster]]
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


def GenerateClusterFASTAs(genomes, refined=False):
    """
    Extract gene model clusters from full database and write out nucleotide and protein sequence families to file.
    """
    if not os.path.isfile("./gm_pred/sets/allnucl.db"):
        ConcatenateDatasets(genomes)
    elif not os.path.isfile("./gm_pred/sets/allprot.db"):
        ConcatenateDatasets(genomes)
    nt_index = SeqIO.index("./gm_pred/sets/allnucl.db", "fasta")
    aa_index = SeqIO.index("./gm_pred/sets/allprot.db", "fasta")
    fdir = "./panoct/clusters/"
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

    if refined:
        matchtable = "./panoct/refined_matchtable.txt"
        rdir = "./panoct/clusters/refined"
        TryMkDirs(rdir)
        TryMkDirs("{0}/core/faa".format(rdir))
        TryMkDirs("{0}/core/fna".format(rdir))
        TryMkDirs("{0}/acc/faa".format(rdir))
        TryMkDirs("{0}/acc/fna".format(rdir))

        core, acc = ParseMatchtable(matchtable)

        for cluster in core:
            nt_seqs = [nt_index[member] for member in core[cluster]]
            aa_seqs = [aa_index[member] for member in core[cluster]]
            with open("{0}/core/fna/Core_{1}.fna".format(rdir, cluster), "w") as aa_out:
                SeqIO.write(nt_seqs, aa_out, "fasta")

            with open("{0}/core/faa/Core_{1}.faa".format(rdir, cluster), "w") as aa_out:
                SeqIO.write(aa_seqs, aa_out, "fasta")

        for cluster in acc:
            nt_seqs = [nt_index[member] for member in acc[cluster] if member]
            aa_seqs = [aa_index[member] for member in acc[cluster] if member]
            with open("{0}/acc/fna/Acc_{1}.fna".format(rdir, cluster), "w") as aa_out:
                SeqIO.write(nt_seqs, aa_out, "fasta")

            with open("{0}/acc/faa/Acc_{1}.faa".format(rdir, cluster), "w") as aa_out:
                SeqIO.write(aa_seqs, aa_out, "fasta")