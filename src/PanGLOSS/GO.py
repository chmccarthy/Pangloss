# -*- coding: utf-8 -*-
import os
import subprocess as sp
from csv import reader


def MakeWorkingDirs():
    """"""
    tdir = "go"
    try:
        os.makedirs(tdir)
    except OSError as e:
        if e.errno != os.errno.EEXIST:
            logging.info("PanOCT: Program output directory already exists, using it instead.")
            raise


def GenerateAnnoDict(ips):
    """

    :param ips:
    :return:
    """
    with open(ips) as infile:
        anno_dict = {}
        for line in reader(infile, delimiter="\t"):
            if len(line) == 14:
                anno_dict[protein]["GO"] = anno[protein]["GO"] + [go for go in line[13].split("|") if go]
    return anno_dict


def get_complement_enrichments(core, noncore, annotations, go_terms, go_slim):
	with open("go/pangenome_associations.txt", "w") as passocs, open("go/pangenome_population.txt", "w") as panpopl:
		print "Getting association and background population for pangenome..."
		for gene in annotations:
			if annotations[gene]["GO"]:
				passocs.write("{0}\t{1}\n".format(gene, ";".join(annotations[gene]["GO"])))
				panpopl.write("{0}\n".format(gene))
	sp.call(["map_to_slim.py", "--association_file=go/pangenome_associations.txt", go_terms, go_slim], stdout=open("go/pangenome_slim_temp.txt", "w"))
	with open("pangenome_slim.txt", "w") as panslim:
		print "Tidying GO Slim file..."
		for line in open("pangenome_slim_temp.txt").readlines():
			if "|" in line:
				panslim.write(line)
	with open("core_population.txt", "w") as corepop:
		print "Getting core study population..."
		for cluster in core:
			for gene in core[cluster]:
				if gene in annotations:
					if annotations[gene]["GO"]:
						corepop.write("{0}\n".format(gene))
	with open("noncore_population.txt", "w") as noncorepop:
		print "Getting noncore study population..."
		for cluster in noncore:
			for gene in noncore[cluster]:
				if gene != "----------":
					if gene in annotations:
						if annotations[gene]["GO"]:
							noncorepop.write("{0}\n".format(gene))




def GenerateAssociations(annos):
    """
    """
    with open("go/associations.txt", "w") as assocs:
        for gene in annos:
            if annos[gene]["GO"]:
                assocs.write("{0}\t{1}\n".format(gene, ";".join(annotations[gene]["GO"])))


def GeneratePopulations(annos, core, acc):
    """
    """
    c_pop = [key for key in core if key in annos]
    a_pop = [key for key in acc if key in annos]
    full_pop = c_pop + a_pop
    with open("go/core_pop.txt", "w") as cp_file, open("go/acc_pop.txt", "w") as ap_file,\
         open("go/full_pop.txt", "w") as fp_file:
        cp_file.write("\n".join(c_pop))
        ap_file.write("\n".join(a_pop))
        fp_file.write("\n".join(full_pop))


def GenerateSlimData(assocs, go_obo, slim_obo):
    """
    """
    sp.call(["map_to_slim.py", "--association_file={0}".format(assocs), go_obo, slim_obo],
            stdout=open("go/pangenome_slim_temp.txt", "w"))
    with open("pangenome_slim.txt", "w") as panslim:
        print "Tidying GO Slim file..."
        for line in open("pangenome_slim_temp.txt").readlines():
            if "|" in line:
                panslim.write(line)

def CoreEnrichment():
    """
    """


def AccessoryEnrichment():
    """
    """