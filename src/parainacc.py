from __future__ import division

import datetime
import multiprocessing as mp
import os
import subprocess as sp
import time
from collections import Counter
from csv import reader
from glob import glob

from Bio import SearchIO, SeqIO
from PanGLOSS.Tools import flatten, grouper, seq_ratio
from PanGLOSS.Tools import merge_clusters
from PanGLOSS.Tools import subject_top_hit, query_top_hit, query_hit_dict, subject_hit_dict

def paralog_finder(blast_results, seqindex, noncore, strain_cutoff):
    paralogs = {}
    count = 0
    for query_cluster in noncore:
        count = count + 1
        subblastlog.write("Query cluster {0}...".format(count))
        already_matched = []
        members = filter(lambda x: x != "----------", noncore[query_cluster])
        strains_in_query = [member.split("|")[0] for member in members]
        query_cluster_length = len(members)
        member_dict = query_hit_dict(members, blast_results, 30)
        for protein in member_dict:
            for hit in member_dict[protein]:
                if hit.split("|")[0] in strains_in_query:
                    if seq_ratio(seqindex, protein, hit) >= 0.6:
                        singlehitcount = len(filter(lambda x: x == hit, flatten(member_dict.values())))
                        if singlehitcount / query_cluster_length >= strain_cutoff:
                            if hit in flatten(noncore.values()):
                                for subject_cluster in noncore:
                                    if hit in noncore[subject_cluster]:
                                        if subject_cluster != query_cluster:
                                            if subject_cluster != paralogs:
                                                if subject_cluster != flatten(paralogs.values()):
													if subject_cluster not in already_matched:
														subject_members = filter(lambda x: x != "----------", noncore[subject_cluster])
														strains_in_subject = [subject.split("|")[0] for subject in subject_members]
														subject_cluster_length = len(subject_members)
														strain_match = set(set(strains_in_query) & set(strains_in_subject))
														print strains_in_query, strains_in_subject, strain_match
														if len(strain_match) / subject_cluster_length >= strain_cutoff:
															strain_dict = [member_dict[key] for key in member_dict if key.split("|")[0] in strain_match]
															subhitsinquery = len(filter(lambda x: x in subject_members, set(flatten(strain_dict))))
															if subhitsinquery / len(strain_match) >= strain_cutoff:
																subj_dict = subject_hit_dict(subject_members, blast_results, 30)
																strain_dict = [subj_dict[key] for key in subj_dict if key.split("|")[0] in strain_match]
																reciphitcount = len(filter(lambda x: x in members, set(flatten(strain_dict))))
																if reciphitcount / len(strain_match) >= strain_cutoff:
																	already_matched.append(subject_cluster)
																	print subject_cluster
																	if query_cluster in paralogs:
																		paralogs[query_cluster].append(subject_members)
																	else:
																		paralogs[query_cluster] = [subject_members]
    return paralogs
    
    
    
    
noncore = {}
matchtable = reader(open("../work/new_nontable.txt"), delimiter="\t")
for row in matchtable:
    noncore[row[0]] = row[1:]

seqindex = SeqIO.index("../work/panoct_db.fasta", "fasta")
blast_results = SearchIO.index("../blast_results.txt", "blast-tab")
strain_cutoff = 1

subblastlog = open("accparalog.txt", "w")

paralogs = paralog_finder(blast_results, seqindex, noncore, strain_cutoff)
with open("acc_para.txt", "w") as output:
    for res in paralogs:
        output.write("{0}\t{1}\n".format(res, ",".join(flatten(paralogs[res]))))
        print len(paralogs.values())