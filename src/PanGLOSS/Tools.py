# -*- coding: utf-8 -*-
"""
Short functions used throughout PanGLOSS and PanGuess.

Functions imported explictly via "from PanGLOSS.Tools import <name>".
"""

from __future__ import division

import cStringIO
import subprocess as sp
from Bio import SeqIO, SeqRecord
from csv import reader
from ExonerateGene import ExonerateGene
from difflib import SequenceMatcher
from itertools import chain, izip_longest, tee


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

def merge_clusters(larger_cluster, smaller_cluster):
    for index, member in enumerate(larger_cluster):
        if member == "----------":
            larger_cluster[index] = smaller_cluster[index]
    return larger_cluster

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

def subject_hit_dict(subject_cluster, blast_results, min_id_cutoff):
    """
    Generate dictionary of all hits for all members of a subject cluster >min_id_cutoff identity.

    FTR I think this is almost identical in function to query_hit_dict.
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

    Arguments:
        cluster_members = List of proteins in query cluster.
        strain_list = List of strains in query cluster.
        blast_hits = BLAST hit dictionary of subject cluster.
        size =

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


def Pairwise(iterable):
    """
    Enables pairwise iteration. Taken from the Python Standard Library.
    """
    a, b = tee(iterable)
    next(b, None)
    return izip_longest(a, b)  # Allows (line, None) for EOF.


def Flatten(iterable):
    """
    Flatten a list of lists, essential for ClusterClean and GapFinder.

    Taken from the Python Standard Library.
    """
    return list(chain.from_iterable(iterable))


def ExonerateCmdLine(cmd):
    """
    Carries out an exonerate command and return output as a ExonerateGene object.

    If an exonerate command does not find a suitable homolog to the query gene
    within the target genome (which is fine!), then the output will fail to be
    passed as a ExonerateGene object correctly (which makes sense, as there's
    no information to make an object from). As such, the contains check makes
    sure only full exonerate hits are returned.
    """
    process = sp.check_output(cmd)
    if "C4 Alignment:" in process:  # Empty results don't contain this line!
        return ExonerateGene(cStringIO.StringIO(process))
    else:
        pass


def LocationOverlap(call, next_call):
    """
    Check overlapping co-ordinates for calls via exonerate vs. GeneMark-ES.
    """
    overlap = False
    call_len = int(call[3]) - int(call[2])
    next_call_len = int(next_call[3]) - int(next_call[2])
    
    if (int(call[2]) <= int(next_call[2]) <= int(call[3])):
        overlap = True
    elif (int(call[2]) <= int(next_call[3]) <= int(call[3])):
        overlap = True
    elif (int(next_call[2]) - int(call[3]) <= 20):
        print call[3], next_call[2]
        overlap = True
    
    if overlap:
        if call_len > next_call_len:
            return next_call
        else:
            return call


def MakeBLASTDBCmdLine(cmd):
    """
    Generalized function for running makeBLASTDB. Nothing fancy.
    """
    sp.call(cmd)


def QCBLASTCmdLine(cmd):
    """
    BLAST pseudogenes/transposable elements/&c against a given gene model set.
    """
    process = sp.check_output(cmd)
    if "<BlastOutput>" in process:  # Empty results don't contain this line!
        return process
    else:
        pass


def StringBLAST(query):
    """
    Runs BLASTp against an intended all-vs.-all database given a valid FASTA gene model as
    a pipeable string. We run BLASTp with an output format set to tabular with comments to
    enable a check for empty results (see if line).
    """
    cmd = ['blastp', '-db', 'allprot.db', '-evalue', '0.0001', '-outfmt', '7', '-query', "-"]
    process = sp.Popen(cmd, stdin=sp.PIPE, stdout=sp.PIPE)
    output = process.communicate(query)
    if not "# 0 hits found" in output[0]:  # Empty results don't contain this line!
        return output[0]
    else:
        pass


def ParseMatchtable(matchtable):
    """
    """
    core = {}
    acc = {}
    count = 1
    clusters = reader(open(matchtable), delimiter="\t")

    for cluster in clusters:
        if "----------" in cluster:
            members = [gene if gene != "----------" else None for gene in cluster]
            acc[count] = members
        else:
            core[count] = cluster
        count = count + 1
    return core, acc


def StringMUSCLE(seqs):
    """
    Runs a MUSCLE alignment given a valid set of translated nucleotides as stdin, returns
    the alignment to stdout which is then processed within memory.
    """
    cmd = ["muscle", "-quiet"]
    process = sp.Popen(cmd, stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.PIPE)
    SeqIO.write(seqs, process.stdin, "fasta")
    process.stdin.close()
    return process.stdout.read()


def Untranslate(aseq, nseq):
    """
    Untranslates translated nucleotide alignments back to nucleotides but keeps codons together.
    """
    unseq = ""
    locs = [0, 3]
    stops = ["TAG", "TAA", "TGA"]
    for site in aseq:
        if site == "-":
            unseq += "---"
        else:
            if nseq[locs[0]:locs[1]] in stops:
                unseq += "AAA"
            else:
                unseq += nseq[locs[0]:locs[1]]
            locs[0] = locs[0] + 3
            locs[1] = locs[1] + 3
    return SeqRecord.SeqRecord(unseq)


def QueryClusterFirstHits(q_cluster, blast_idx, ident, tags):
    """
    Generate dictionary of all hits for all members of a query cluster >min_id_cutoff identity.
    """
    hit_dict = {member: [hit.id for hit in blast_idx[member].hits if hit.hsps[0].ident_pct
                               >= float(ident)] for member in q_cluster if member in blast_idx}
    for member in q_cluster:
        top_hits = []
        if member in hit_dict:
            for tag in tags:
                top_hits.append(next((hit for hit in hit_dict[member] if hit.startswith(tag)), None))
            hit_dict[member] = top_hits
    return hit_dict


def Reciprocal(q_cluster, q_first_hits, s_cluster, s_first_hits):
    """
    Return True if all first hits for the source strain of a cluster member are the member itself, otherwise False.
    """
    reciprocal = False
    q_members = set([member for member in q_cluster if member])
    q_value = bool(q_members <= s_first_hits)
    s_members = set([member for member in s_cluster if member])
    s_value = bool(s_members <= q_first_hits)
    if all([q_value, s_value]):
        reciprocal = True
    return reciprocal
