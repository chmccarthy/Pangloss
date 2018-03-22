# -*- coding: utf-8 -*-
"""
PanGLOSS v0.1.2: PanGenomic Largescale Orthology & microSynteny Search.

Python software (with some future Cython and R components) for pangenomic
analysis of eukaryotes. In something of a provisional state at the moment.

PanGLOSS is the sister program to PanGuess, and can be used in conjuction
to determine pangenomic structure of species of interest using only
genomic data and a reference protein set through homology searches, HMM and
PVM gene prediction techniques.

Requirements:
    - Python (written for 2.7.x)
        - BioPython (1.70)
    - PanOCT (>3.23)
    - R (>3.0)
        - Cairo (>1.5)
        - UpSetR (>1.3)
    - MacOS (tested on MacOS 10.12) or Linux (tested on SLES 11).

Recent changes:

    v0.1.4 (March 2018)
    - Implementing full cluster merging (as opposed to n-1 + 1 merging).
    - Improved logging.

	v0.1.3 (March 2018)
	- Added R scripts.

    v0.1.2 (January 2018)
    - Added parallelized BLASTp searches via multiprocessing.Pool.
    - Made parallel_BLAST code easier to read.

    v0.1.1 (Winter 2017)
    - Completely rewritten to link in with PanGuess gene prediction pipeline.

    v0.1.0 (Autumn 2017)
    - Initial test version based on GFF files taken from NCBI, AUGUSTUS, &c.

To-do:

    - Integrate functions from PanGuess into PanGLOSS.
    - Improve CLI.
    - Write replacement orthology search.
    - Convert some of the more maths-y parts to Cython if feasible.

Written by Charley McCarthy, Genome Evolution Lab, Department of Biology,
Maynooth University in 2017-2018 (Charley.McCarthy@nuim.ie).
"""

from __future__ import division

import datetime
import multiprocessing as mp
import os
import subprocess as sp
import time
from csv import reader
from glob import glob

from Bio import SearchIO, SeqIO
from PanGLOSS.Tools import flatten, grouper, seq_ratio
from PanGLOSS.Tools import merge_clusters
from PanGLOSS.Tools import subject_top_hit, query_top_hit


##### Here is the function to handle individual BLASTp searches in parallel_BLAST. #####
def subprocess_BLAST(cmd):
    """
    Run an individual instance of a parallel_BLAST search.
    """
    subblastlog.write("Running {0}".format(" ".join(cmd)) + "\n")
    sp.call(cmd)


##### Here are the major functions used in cluster_clean. #####
def parallel_BLAST(list_of_genes, seqindex, split_by, out):
    """
    Run a BLASTp all-vs.-all search on a subset of genes from a database.

    Splits queries into separate files by a divisor split_by, BLASTs them
    simultaneously and then concatenates them using cat. Currently uses
    subprocess over BioPython's BLASTp wrapper for easy parallelization
    using mp.Pool and subprocess_BLAST. Returns a SearchIO BLAST tabular object.

    Requires cat, makeblastdb and blastp in your $PATH!

    Arguments:
        list_of_genes = List of (remaining) noncore protein IDs.
        seqindex      = SeqIO.index of all proteins.
        split_by      = Divisor to split files into.
        out           = File prefix for results files given current cluster size being investigated.
    """

    ##### Generate FASTA database of (remaining) noncore proteins. #####
    outfast = open(out, "w")
    for seq in list_of_genes:
        outfast.write(">{0}\n{1}\n".format(seqindex[seq].id, seqindex[seq].seq))
    outfast.close()  # Close file here, makeblastdb has problems otherwise.
    subblastlog.write("Wrote {0} sequences to {1}...\n".format(len(list_of_genes), out))

    ##### Run makeblastdb on FASTA database. #####
    sp.call(["makeblastdb", "-in", out, "-dbtype", "prot", "-out", "{0}.db".format(out)], stdout=subblastlog)

    ##### Split FASTA database by split_by and generate list of BLASTp commands. #####
    count = 0
    query_cmds = []  # Handle for simultaneous BLASTp commands.
    to_split = SeqIO.parse("{0}".format(out), "fasta")
    for part in grouper(to_split, int(round(len(list_of_genes) / split_by))):
        count = count + 1
        seqs = filter(lambda x: x is not None, part)  # Remove fill values.
        if glob("{0}.part{1}.faa".format(out, count)):  # Remove previous versions of file if present.
            os.remove("{0}.part{1}.faa".format(out, count))
        SeqIO.write(seqs, "{0}.part{1}.faa".format(out, count), "fasta")
        query_cmds.append(["blastp", "-query", "{0}.part{1}.faa".format(out, count),
                           "-db", "{0}.db".format(out), "-outfmt", "6 std qlen slen", "-evalue",
                           "0.0001", "-out",
                           "{0}.part{1}.subblast".format(out, count),
                           "-num_threads", "1"])
    subblastlog.write("Split original query file {0} ({1} sequences) into {2} files.\n".format(out, len(list_of_genes),
                                                                                               str(split_by)))

    ##### Run BLASTp processes simultaneously using mp.Pool. #####
    farm = mp.Pool(processes=split_by)
    farm.map(subprocess_BLAST, query_cmds)
    farm.close()
    subblastlog.write(
        "Finished BLAST+ searches for {0} ({1} sequences), split into {2} files.\n".format(out, len(list_of_genes),
                                                                                           str(split_by)))

    ##### Concatenate parallel_BLAST results together in the shell and remove other files. #####
    sp.call(["cat"] + glob("*subblast"), stdout=open("{0}.results".format(out), "wb"))
    for bin_file in glob("%s.db*" % out):
        os.remove(bin_file)
    for part_file in glob("%s.part*" % out):
        os.remove(part_file)
    os.remove("%s" % out)
    subblastlog.write(
        "Concatenated BLAST+ output for {0} ({1} sequences).\n".format(out, len(list_of_genes), str(split_by)))

    ##### Parse parallel_BLAST results and return them as a SearchIO.index instance to cluster_clean. #####
    blast = SearchIO.index("{0}.results".format(out), "blast-tab", fields=blast_fields)
    subblastlog.write("Loaded {0} as SearchIO.index for cluster_clean.\n".format(out))
    return blast


def gap_finder(blast_results, seqindex, noncore, total, current, min_id_cutoff, strain_cutoff):
    """
    Find potential "gaps" in noncore clusters arising from microsynteny loss.

    Workflow:
        1.  Loop through each noncore cluster.
        2.  Ignore clusters whose number of members is not the current size
            under investigation.
        3.  Generate dictionary of all BLAST hit IDs for every member of a cluster over the identity cutoff.
        4a. Loop through each "query" cluster in the BLAST hits dictionary.
        4b. If no homolog has already been found, loop through the hits of a given member of the cluster.
        4c. If one of the hits is from a strain missing from or not already added to the current cluster.
        4d. Calculate sequence length ratio between member protein and hit.
        4e. If ratio is greater than/equal to 60%, check hit in member protein is also the top hit for that strain for
            every other member protein (.
        4f. If it is, get that hit protein's associated cluster and check that its size is not greater than total - current.
        5a. If "hit" cluster is a singleton, check if every member of the "query" cluster is (first) a BLASTp hit and
            (second) the top BLASTp hit for their respective strains for the "hit" cluster protein.
        5b. If so, add that "hit" cluster to a list of suitable "fill" clusters for the "query" cluster.
        6a. If "hit" cluster is not a singleton, check that every member of the "hit" cluster are also subjects
            of every member of the "query" cluster (e.g. if "hit" cluster has two proteins, and "query" cluster
            has 4 proteins, each member of the "query" cluster must have the two "hit" cluster proteins as a BLASTp
            hit such that the number of times the "hit" cluster proteins are present in the BLAST hit dictionary
            for the "query" cluster is equal to 8).
        6b. Check that every strain represented in the "hit" cluster is missing from the "query" cluster.
        6c. Check that every member of the "query" cluster is (first) a BLASTp hit and
            (second) the top BLASTp hit for their respective strains for each "hit" cluster protein.
        6d. If so, add that "hit" cluster to a list of suitable "fill" clusters for the "query" cluster.


    Returns a dictionary of "query" clusters whose values are "subject" clusters which pass these criteria.

    Arguments:
        blast_results = All-vs.-all noncore proteins BLASTp results.
        seqindex      = SeqIO.index of all proteins.
        noncore       = Noncore cluster dictionary.
        total         = Total number of genomes.
        current       = Cluster size being queried.
        min_id_cutoff = Percentage identity of a BLASTp hit (default = 30).
        strain_cutoff = Cutoff fraction of reciprocal strain top hits between two clusters (default = 1).
    """
    homologs = {}
    for cluster in noncore.keys():  # Loop through non-core clusters.
        found = []  # Default for strains "added" to cluster.
        members = filter(lambda x: x != "----------", noncore[cluster])  # Get actual cluster.
        if len(members) != current:
            pass  # Ignore clusters outside of the current size.
        else:  # If cluster size = current.
            gap = False  # Default state of our cluster.
            blast_hit_dict = {
                member: [hit.id for hit in blast_results[member].hits if hit.hsps[0].ident_pct >= float(min_id_cutoff)]
            for member in members if member in blast_results}  # Generate dict of (cluster member: BLAST hit IDs with ID > 30%) for every cluster member.
            for key in blast_hit_dict:  # Loop through each protein in the "query cluster".
                if not gap:  # If no suitable homolog has been found.
                    for hit in blast_hit_dict[key]:  # Loop through every hit from a given query cluster protein.
                        if hit.split("|")[0] not in [i.split("|")[0] for i in
                                                     blast_hit_dict.keys()]:  # If the source strain of a given hit is missing from the query cluster.
                            if hit.split("|")[0] not in found:  # Ignore if missing strain has already been "added" to cluster.
                                if seq_ratio(seqindex, key,
                                             hit) >= 0.6:  # If the sequence length ratio (see seq_ratio function) between hit and query is greater than or = 60%.
                                    print key
                                    if len(filter(lambda x: x == hit, flatten(
                                            blast_hit_dict.values()))) / current == strain_cutoff:  # If that hit is found in the BLASTp results of every member of the query cluster.
                                        if subject_top_hit(blast_hit_dict.values(), hit, current,
                                                           strain_cutoff):  # If that hit is the top hit from the source strain of every member of the query cluster.
                                            for query_cluster in noncore.keys():  # Loop through noncore clusters AGAIN.
                                                if hit in noncore[query_cluster]:  # Locate hit's source cluster within noncore clusters, henceforth the "subject cluster".
                                                    print query_cluster
                                                    if (total - current) >= len(filter(lambda x: x != "----------",
                                                                                       noncore[
                                                                                           query_cluster])):  # If the size of the subject cluster is =< the number of missing strains from the query cluster.
                                                        subjhits = {subj: [hit.id for hit in blast_results[subj].hits if
                                                                           hit.hsps[0].ident_pct >= float(
                                                                               min_id_cutoff)] for subj in
                                                                    noncore[query_cluster] if subj in blast_results}
                                                        if len(filter(lambda x: x != "----------", noncore[query_cluster])) > 1:  # If the subject cluster is not a singleton cluster.
                                                            if len(filter(lambda x: x in noncore[query_cluster],
                                                                          flatten(
                                                                              blast_hit_dict.values()))) % current == 0:  # If all proteins in the subject cluster are hits of the query cluster (i.e. if all proteins from a 2-size subject cluster have reciprocality with a 10-size query cluster they should be present 20 times in the query clusters' BLASTp results.)
                                                                spec = [i.split("|")[0] for i in
                                                                        filter(lambda x: x != "----------", noncore[
                                                                            query_cluster])]  # Strains present in the subject cluster.
                                                                if not filter(lambda x: x in spec,
                                                                              [i.split("|")[0] for i in
                                                                               blast_hit_dict.keys()]):  # If all strains present in the subject cluster are missing from the query cluster.
                                                                    if len(filter(lambda x: x in blast_hit_dict.keys(),
                                                                                  set(flatten(
                                                                                      subjhits.values())))) == current:
                                                                        strains = [key.split("|")[0] for key in
                                                                                   blast_hit_dict.keys()]
                                                                        if query_top_hit(blast_hit_dict.keys(), strains,
                                                                                         subjhits.values(), current,
                                                                                         strain_cutoff):
                                                                            for s in spec:
                                                                                found.append(s)
                                                                            gap = True
                                                                            if cluster in homologs.keys():  # Allow more than one cluster to be associated.
                                                                                homologs[cluster].append(query_cluster)
                                                                            else:
                                                                                homologs[cluster] = [query_cluster]
                                                        else:
                                                            if len(filter(lambda x: x in blast_hit_dict.keys(), flatten(
                                                                    subjhits.values()))) / current == strain_cutoff:  # If every member of the query cluster is also a subject of the protein in the singleton subject cluster.
                                                                strains = [key.split("|")[0] for key in
                                                                           blast_hit_dict.keys()]
                                                                if query_top_hit(blast_hit_dict.keys(), strains,
                                                                                 subjhits.values(), current,
                                                                                 strain_cutoff):
                                                                    gap = True  # Suitable homolog has been found.
                                                                    found.append(hit.split("|")[
                                                                                     0])  # Shortcut: append hit ID substring because we're only looking at a singleton.
                                                                    if cluster in homologs.keys():  # Allow more than one subject cluster to be associated to a query cluster.
                                                                        homologs[cluster].append(query_cluster)
                                                                    else:
                                                                        homologs[cluster] = [query_cluster]
    return homologs


def run_PanOCT(blast_results, fasta_handle, attributes, tags):
    """
    Wrapper for running PanOCT within PanGLOSS.

    Expects panoct.pl in $PATH!
    """
    sp.call(["panotc.pl", "-t", blast_results, "-f", tags,
             "-g", attributes, "-P", fasta_handle])


def cluster_clean(panoct_clusters, fasta_handle, split_by=4, min_id_cutoff=30, strain_cutoff=1):
    """
    Tidy up non-core clusters found by PanOCT.

    Feeds into parallel_BLAST, which expects you to have BLAST+ installed.
    Also feeds into gap_finder, which doesn't require anything else.

    Workflow:
        1.  Load in FASTA database and clusters (matchtable.txt) from PanOCT.
        2a. Identify core clusters (i.e. size = total number of genomes).
        2b. Idenfify noncore clusters (i.e. size > total number of genomes).
        3a. Loop through noncore clusters of size = (total genomes - 1).
        3b. Get IDs of all proteins from noncore clusters.
        3c. Run all-vs.-all BLASTp search (parallel_BLAST) on all noncore proteins,
            parallelized by the split_by variable (default = 4).
        3d. Get concatenated parallel_BLAST results together.
        3e.
    """
    ##### Load in FASTA database and PanOCT results. #####
    db = SeqIO.index(fasta_handle, "fasta")
    matchtable = reader(open(panoct_clusters), delimiter="\t")

    ##### Initialize empty dictionaries for cluster types. #####
    core = {}
    noncore = {}
    softcore = {}

    ##### Initialize variables for total/starting number of genomes. #####
    total = 0
    start = 0

    ##### Populate core & noncore dictionaries by reading PanOCT results. #####
    for row in matchtable:
        if "----------" in row:
            noncore[row[0]] = row[1:]  # Populating our initial noncore dict.
        else:
            core[row[0]] = row[1:]  # Populating our core dict.
            if total == 0:
                total = len(row) - 1  # Total number of genomes.
                start = total - 2  # Max noncore cluster size.

    mainlogfile.write(
        "{0} core clusters and {1} noncore clusters identified...\n".format(len(core.keys()), len(noncore.keys())))
    mainlogfile.write("Creating ring chart in R...\n")
    sp.call(["Rscript", "{0}/PlotRingChart.R".format(dirname), str(len(core.keys())), str(len(noncore.keys()))])

    ##### Loop through noncore clusters from size (total -1) to 2. #####
    for size in range(start, 1, -1):
        filled_count = 0
        merged_count = 0
        ##### Get list of (remaining) noncore protein IDs. #####
        to_blast = filter(lambda x: x != "----------", flatten([noncore[key] for key in noncore]))
        mainlogfile.write("All-vs.-all BLAST of {0} proteins...\n".format(len(to_blast)))

        ##### Run parallel_BLAST. #####
        results = parallel_BLAST(to_blast, db, split_by, "ClusterBLAST_{0}.fasta".format(str(size)))
        mainlogfile.write("Finding potential homology gaps in clusters of size {0}...\n".format(str(size)))

        ##### Run gap_finder. #####
        gaps = gap_finder(results, db, noncore, total, size, min_id_cutoff, strain_cutoff)

        ##### Identify clusters that need to be merged and move merged clusters to appropriate dictionary. #####
        for cluster in gaps:
            if cluster in noncore.keys():
                for candidate in gaps[cluster]:
                    if candidate in noncore.keys():
                        merge_size = (len(filter(lambda x: x != "----------", noncore[cluster])) + len(
                            filter(lambda x: x != "----------", noncore[candidate])))
                        if merge_size == total:
                            mainlogfile.write("{0} (size: {1}) has a homologous cluster: {2} (size: {3})\n".format
                                              (cluster, len(filter(lambda x: x != "----------", noncore[cluster])),
                                               candidate, len(filter(lambda x: x != "----------", noncore[candidate]))))
                            mainlogfile.write(
                                "Merging smaller cluster {0} into larger cluster {1}...\n".format(candidate, cluster))
                            mainlogfile.write("Merged cluster {0} has size {1}.\n".format(cluster, merge_size))
                            filled = merge_clusters(noncore[cluster], noncore[candidate])
                            softcore[cluster] = filled
                            del noncore[cluster], noncore[candidate]
                            filled_count = filled_count + 2
                    else:
                        if merge_size > total:
                            mainlogfile.write("{0} (size: {1}) has a homologous cluster: {2} (size: {3})\n".format
                                              (cluster, len(filter(lambda x: x != "----------", noncore[cluster])),
                                               candidate, len(filter(lambda x: x != "----------", noncore[candidate]))))
                            mainlogfile.write(
                                "Merging smaller cluster {0} into larger cluster {1}...\n".format(candidate, cluster))
                            mainlogfile.write("Merged cluster {0} has size {1}.\n".format(cluster, merge_size))
                            merged = merge_clusters(noncore[cluster], noncore[candidate])
                            noncore[cluster] = merged
                            del noncore[candidate]
                            merged_count = merged_count + 2

        mainlogfile.write(
            "At cluster size (n = {0}): merged {1} homologous clusters into {2} softcore clusters.\n".format(size,
                                                                                                             filled_count,
                                                                                                             filled_count / 2))
        mainlogfile.write(
            "At cluster size (n = {0}): merged {1} homologous clusters into {2} noncore lusters.\n".format(size,
                                                                                                           filled_count,
                                                                                                           merged_count / 2))

    if not os.path.isdir("{0}/sub_BLASTs".format(os.getcwd())):
        os.makedirs("{0}/sub_BLASTs/faa".format(os.getcwd()))
        os.makedirs("{0}/sub_BLASTs/results".format(os.getcwd()))
    for sub_faa in glob("ClusterBLAST_*.fasta"):
        os.rename(sub_faa, "{0}/sub_BLASTs/faa/{1}".format(os.getcwd(), sub_faa))
    for sub_results in glob("*.results"):
        os.rename(sub_results, "{0}/sub_BLASTs/results/{1}".format(os.getcwd(), sub_results))
    mainlogfile.write("Remaining noncore clusters after sub_BLASTing: {0}\n".format(len(noncore.keys())))
    mainlogfile.write(
        "PanGLOSS analysis finished in {0} seconds. Thank you for choosing PanGLOSS, the friendly pangenome software.\n".format(
            time.time() - start_time))
    mainlogfile.write("=== Finished PanGLOSS job at {0}. ===\n".format(str(datetime.datetime.now())))


##### Here we define the workflow for PanGLOSS. #####


def main():
    """
    Main software workflow.
    """
    cluster_clean("matchtable.txt", "panoct_db.fasta", split_by=cores)


if __name__ == "__main__":
    ##### Open log files. #####
    mainlogfile = open("PanGLOSS.log", "a", 0)  # Flush to log immediately.
    mainlogfile.write("\n=== Started PanGLOSS job at {0}. ===\n".format(str(datetime.datetime.now())))
    start_time = time.time()
    subblastlog = open("parallel_BLAST.log", "a", 0)  # Log for parallel_BLAST.

    ##### Define custom columns for parallel_BLAST. #####
    blast_fields = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
                    'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore',
                    'qlen', 'slen']

    ##### Set default amount of cores used. #####
    cores = mp.cpu_count() - 1

    ##### Get absolute path of script, for running R commands. #####
    dirname = os.path.dirname(os.path.abspath(__file__))

    ##### Run PanGLOSS. #####
    main()
