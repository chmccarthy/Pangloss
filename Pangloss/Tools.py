# -*- coding: utf-8 -*-
"""
Short functions used throughout Pangloss and PanGuess.

Functions imported explictly via "from Pangloss.Tools import <name>".
"""

from __future__ import division

import cStringIO
import os
import subprocess as sp
from collections import Counter, OrderedDict as od
from csv import reader
from itertools import chain, izip_longest, tee

from Bio import SeqIO, SeqRecord

from ExonerateGene import ExonerateGene


def TryMkDirs(path):
    """
    Tries to make directory at specified path.
    """
    try:
        os.makedirs(path)
    except OSError as e:
        if e.errno != os.errno.EEXIST:
            logging.info("TryMkDirs: {0} directory already exists, using it instead.".format(path))
            raise


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
        overlap = True
    
    if overlap:
        if call_len > next_call_len:
            return next_call
        else:
            return call


def ConcatenateDatasets(genomes):
    """
    Concatenate all datasets and construct BLASTp database for gene model set.
    """
    # Generate cat commands for the three full datasets we have.
    tags = [line.strip("\n").split(".")[0].split("/")[1] for line in open(genomes)]
    nucl_cmd = ["cat"] + ["./gm_pred/sets/" + tag + ".nucl" for tag in tags]
    prot_cmd = ["cat"] + ["./gm_pred/sets/" + tag + ".faa" for tag in tags]
    att_cmd = ["cat"] + ["./gm_pred/sets/" + tag + ".attributes" for tag in tags]

    # Run commands.
    with open("./gm_pred/sets/allnucl.db", "w") as f:
        sp.call(nucl_cmd, stdout=f)
    with open("./gm_pred/sets/allprot.db", "w") as f:
        sp.call(prot_cmd, stdout=f)
    with open("./gm_pred/sets/allatt.db", "w") as f:
        sp.call(att_cmd, stdout=f)


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
    cmd = ['blastp', '-db', './gm_pred/sets/allprot.db', '-evalue', '0.0001', '-outfmt', '7', '-query', "-"]
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


def UnparseMatchtable(components):
    """
    """
    with open("panoct/refined_matchtable.txt", "w") as outfile:
        for comp in components:
            for cluster in comp:
                members = [gene if gene != "None" else "----------" for gene in comp[cluster]]
                outfile.write("\t".join(members))
                outfile.write("\n")


def ParseKaryotypes(karyotypes):
    """
    """
    karyoreader = reader(open(karyotypes), delimiter="\t")
    karyodict = od()

    for row in karyoreader:
        karyodict[row[1]] = [row[0], row[2], row[3]]

    return karyodict


def StringMUSCLE(ml_path, seqs):
    """
    Runs a MUSCLE alignment given a valid set of translated nucleotides as stdin, returns
    the alignment to stdout which is then processed within memory.
    """
    cmd = [ml_path, "-quiet"]
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


def MultipleInsert(i_cluster, tags):
    """
    Add multiple nones into potential matching cluster for gap filling.
    """
    query_cluster = [None] * len(tags)
    for tag in tags:
        for gene in i_cluster:
            if gene.split("|")[0] == tag:
                query_cluster[tags.index(tag)] = gene
    return query_cluster


def Reciprocal(q_members, q_first_hits, s_members, s_first_hits):
    """
    Return True if all first hits for the source strain of a cluster member are the member itself, otherwise False.
    """
    reciprocal = False
    q_value = bool(q_members <= s_first_hits)
    s_value = bool(s_members <= q_first_hits)
    if all([q_value, s_value]):
        reciprocal = True
    return reciprocal


def ClusterMerge(q_cluster, s_cluster):
    for index, member in enumerate(q_cluster):
        if not member:
            q_cluster[index] = s_cluster[index]
    return q_cluster


def ClusterSizes(component):
    """
    Return counts of cluster sizes within a component.
    """
    clusters = component.values()
    counts = [len(filter(lambda x: x is not None, cluster)) for cluster in clusters]
    sizes = Counter(counts)
    return sizes
