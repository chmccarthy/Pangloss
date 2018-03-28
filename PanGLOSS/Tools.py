"""
Short functions used throughout PanGLOSS and PanGuess.

Functions imported explictly via "from PanGLOSS.Tools import <name>".
"""

from __future__ import division

import cStringIO
import subprocess as sp

from difflib import SequenceMatcher
from itertools import chain, izip_longest, tee

from Bio import SeqIO
from ExonerateGene import ExonerateGene

def pairwise(iterable):
    """
    Enable pairwise iteration.

    Taken from the Python Standard Library.
    """
    a, b = tee(iterable)
    next(b, None)
    return izip_longest(a, b)  # Allows (line, None) for EOF.


def grouper(iterable, n):
    """
    Return a slice of size n from a iterable.

    Taken from the Python Standard Library.
    """
    args = [iter(iterable)] * n
    return izip_longest(*args, fillvalue=None)


def similar(a, b):
    """
    Compare similarity of two strings.
    """
    return SequenceMatcher(None, a, b).ratio()


def get_gene_lengths(fasta):
    """
    Generate dictionary of sequence length for a given SeqIO.index.
    """
    ref_lengths = {}
    db = SeqIO.index(fasta, "fasta")
    for seq in db:
        ref_lengths[db[seq].id] = len(db[seq].seq)
    return ref_lengths


def flatten(iterable):
    """
    Flatten a list of lists, essential for ClusterClean and GapFinder.

    Taken from the Python Standard Library.
    """
    return list(chain.from_iterable(iterable))


def merge_clusters(larger_cluster, smaller_cluster):
    for index, member in enumerate(larger_cluster):
        if member == "----------":
            larger_cluster[index] = smaller_cluster[index]
    return larger_cluster


def seq_ratio(seqindex, query, subject):
    """
    Return the ratio of the lengths of two sequences.

    Essential for GapFinder, and useful downstream too.
    """
    lengths = [len(seqindex[query].seq), len(seqindex[subject].seq)]
    longest = max(lengths)
    shortest = min(lengths)
    ratio = shortest / longest
    return ratio


def called_ratio(called_alignment, query_gene):
    """
    Return the ratio of lengths for a query sequence and a called gene.
    """
    longest = max(called_alignment, query_gene)
    shortest = min(called_alignment, query_gene)
    ratio = shortest / longest
    return ratio


def subject_top_hit(list_of_lists, gene_id, size, strain_cutoff):
    """
    Return boolean for whether a gene is the top BLASTp hit for its strain.

    This function loops through each set of BLASTp results for a given
    protein cluster as identified by PanOCT, and checks if the gene of
    interest (assuming it passes all prior criteria, see GapFinder for more)
    is the top BLASTp hit from that strain for each member. If this is the case
    for >cutoff of members, it returns the default value of True. If not, it
    returns False. Crucial for GapFinder!
    """
    count = 0
    for li in list_of_lists:
        if filter(lambda x: x.split("|")[0] == gene_id.split("|")[0], li):
            if not filter(lambda x: x.split("|")[0] == gene_id.split("|")[0], li)[0] == gene_id:
                pass  # Hit is not top hit for that strain.
            else:
                count = count + 1
        else:
            pass  # Hit's source strain is not represented in member protein's results.
    if count / size >= strain_cutoff:
        top = True
    else:
        top = False
    return top


def query_hit_dict(members, blast_results, min_id_cutoff):
    """
    Generate dictionary of all hits for all members of a query cluster >min_id_cutoff identity.
    """
    blast_hit_dict = {member: [hit.id for hit in blast_results[member].hits if hit.hsps[0].ident_pct
                               >= float(min_id_cutoff)] for member in members if member in blast_results}
    return blast_hit_dict

def subject_hit_dict(subject_cluster, blast_results, min_id_cutoff):
    """
    Generate dictionary of all hits for all members of a subject cluster >min_id_cutoff identity.
    """
    subjhits = {subj: [hit.id for hit in blast_results[subj].hits if
                hit.hsps[0].ident_pct >= float(min_id_cutoff)] for subj in
                subject_cluster if subj in blast_results}
    return subjhits


def query_top_hit(cluster_members, strain_list, blast_hits, size, strain_cutoff):
    """
    Return boolean for whether a set of genes are all top BLASTp strain hits.

    This function loops through the BLASTp results of a candidate homologous
    subject cluster that has passed all other critera (see GapFinder),
    and checks to see whether all members of the query cluster are the top
    BLASTp hits for their respective strains for every member of the subject
    cluster. If so, return the default True, if not (or if a strain is missing from the
    subject cluster's BLAST results) return False. In this way, we can determine
    reciprocality between query and subject clusters in terms of BLASTp hits. Crucial for GapFinder!
    """
    count = 0
    if any(isinstance(el, list) for el in blast_hits):
        for li in blast_hits:
            for strain in strain_list:
                if not filter(lambda x: x.split("|")[0] == strain, li):
                    pass
                elif filter(lambda x: x.split("|")[0] == strain, li)[0] not in cluster_members:
                    pass
                else:
                    count = count + 1
    else:
        for strain in strain_list:
            if not filter(lambda x: x.split("|")[0] == strain, blast_hits):
                pass
            elif filter(lambda x: x.split("|")[0] == strain, blast_hits)[0] not in cluster_members:
                pass
            else:
                count = count + 1
    if (count / size) >= strain_cutoff:
        top = True
    else:
        top = False
    return top


def exonerate_first_hits(instream):
    """
    Return first exonerate hit per protein-vs.-genome analysis.
    """
    first_align = False
    block = []
    for line in instream:
        if not first_align:
            block.append(line)
            if line.startswith("vulgar"):
                first_align = True
        if "completed exonerate analysis" in line:
            block.append(line)
            first_align = False
    return block


def gene_within(left_end, right_end, query_coords):
    """
    Check overlap of co-ordinates of exonerate gene within known gene.
    """
    overlap = False
    for coord_left, coord_right in pairwise(query_coords):
        if all([int(left_end), int(coord_left), int(right_end), int(left_end), int(coord_right), int(right_end)]):
            if all([int(left_end) <= int(coord_left) <= int(right_end), int(left_end) <= int(coord_right) <= int(right_end)]):
                overlap = True
    return overlap


def gene_overlap(left_gene, right_gene, query_coords, threshold=0):
    """
    Check overlap of co-ordinates of exonerate gene between known genes.
    """
    overlap = False
    if int(left_gene[1]) <= int(query_coords[0]) <= int(left_gene[2]):
        overlap = True
    elif int(right_gene[1]) <= int(query_coords[1]) <= int(right_gene[2]):
        overlap = True
    elif int(left_gene[1]) <= (query_coords[0] - threshold) <= int(left_gene[2]):
        overlap = True
    elif int(right_gene[1]) <= (query_coords[1] + threshold) <= int(right_gene[2]):
        overlap = True
    return overlap


def exoneratecmdline(cmd):
    """
    Carry out an exonerate command and return output as a ExonerateGene object.

    If an exonerate command does not find a suitable homolog to the query gene
    within the target genome (which is fine!), then the output will fail to be
    passed as a ExonerateGene object correctly (which makes sense, as there's
    no information to make an object from). As such, the contains check makes
    sure only full exonerate hits are returned.
    """
    print "Running {0}".format(" ".join(cmd))
    process = sp.check_output(cmd)
    if "C4 Alignment:" in process:  # Empty results don't contain this line!
        return ExonerateGene(cStringIO.StringIO(process))
    else:
        pass
