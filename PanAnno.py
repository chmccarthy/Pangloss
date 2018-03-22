import re
import subprocess as sp
from csv import reader
core = {}
noncore = {}

with open("../matchtable.txt") as infile:
	matchtable = reader(infile, delimiter="\t")
	for row in matchtable:
		if "----------" in row:
			noncore[row[0]] = row[1:]  # Populating our initial noncore dict.
		else:
			core[row[0]] = row[1:]  # Populating our core dict

def generate_annotation_dict(ips_file):
	with open(ips_file) as infile:
		annotations_dict = {}
		for line in reader(infile, delimiter="\t"):
			protein = line[0]
			if protein not in annotations_dict.keys():
				annotations_dict[protein] = {}
				annotations_dict[protein]["PFAM"] = {}
				annotations_dict[protein]["IPR"] = {}
				annotations_dict[protein]["GO"] = []
			if line[4]:
				annotations_dict[protein]["PFAM"][line[4]] = line[5]
			if len(line) > 11:
				annotations_dict[protein]["IPR"][line[11]] = line[12]
			if len(line) == 14:
				annotations_dict[protein]["GO"] = annotations_dict[protein]["GO"] + [go for go in line[13].split("|") if go]
	return annotations_dict


# def generate_mapping_dict(tnt_output):
# 	node_match = re.compile(r"([ ]{3,4})[\S]")
# 	char_match = re.compile(r"([ ]{6})[\S]")
# 	states = {}
# 	apo = False
# 	current_node = ""
# 	with open(tnt_output) as infile:
# 		for line in infile.readlines():
# 			if line.startswith("Tree 0 :"):
# 				apo = True
# 			if node_match.match(line) is not None:
# 				current_node = line.lstrip().split(":")[0].strip()
# 				states[current_node] = []
# 			if apo:
# 				if char_match.match(line) is not None:
# 					if not line.lstrip().startswith("No "):
# 						data = line.lstrip().split(":")
# 						cluster = int(data[0].split(" ")[1]) + 1
# 						if data[1].lstrip().strip() == "1 --> 0":
# 							gain = False
# 						elif data[1].lstrip().strip() == "0 --> 1":
# 							gain = True
# 						states[current_node].append((str(cluster), gain))
# 	return states


def get_strains(core):
	key = core.keys()[0]
	strains = [gene.split("|")[0] for gene in core[key]]
	return strains


def get_genome_enrichments(strain, states, annotations, go_terms, go_slim):
	with open("{0}_associations.txt".format(strain), "w") as assocs, open("{0}_population.txt".format(strain), "w") as popl:
		for gene in annotations.keys():
			if gene.startswith(strain):
				if annotations[gene]["GO"]:
					assocs.write("{0}\t{1}\n".format(gene, ";".join(annotations[gene]["GO"])))
					popl.write("{0}\n".format(gene))
	sp.call(["map_to_slim.py", "--association_file={0}_associations.txt".format(strain), go_terms, go_slim], stdout=open("{0}_slim_temp.txt".format(strain), "w"))
	with open("{0}_slim.txt".format(strain), "w") as slim:
		for line in open("{0}_slim_temp.txt".format(strain)).readlines():
			if line.startswith(strain):
				slim.write(line)
	with open("{0}_gained_popl.txt".format(strain), "w") as gained:
		for cluster in states[strain]:
			if cluster[1] == True:
					member = filter(lambda x: x.startswith(strain), noncore[cluster[0]])[0]
					if member in annotations.keys():
						if annotations[member]["GO"]:
							gained.write("{0}\n".format(member))
	sp.call(["find_enrichment.py", "--pval=0.05", "--method=fdr", "--obo", go_terms, "{0}_gained_popl.txt".format(strain), "{0}_population.txt".format(strain), "{0}_slim.txt".format(strain), "--outfile={0}_enrichment.tsv".format(strain)])


def get_complement_enrichments(core, noncore, annotations, go_terms, go_slim):
	with open("pangenome_associations.txt", "w") as passocs, open("pangenome_population.txt", "w") as panpopl:
		print "Getting association and background population for pangenome..."
		for gene in annotations.keys():
			if annotations[gene]["GO"]:
				passocs.write("{0}\t{1}\n".format(gene, ";".join(annotations[gene]["GO"])))
				panpopl.write("{0}\n".format(gene))
	sp.call(["map_to_slim.py", "--association_file=pangenome_associations.txt", go_terms, go_slim], stdout=open("pangenome_slim_temp.txt", "w"))
	with open("pangenome_slim.txt", "w") as panslim:
		print "Tidying GO Slim file..."
		for line in open("pangenome_slim_temp.txt").readlines():
			if "|" in line:
				panslim.write(line)
	with open("core_population.txt", "w") as corepop:
		print "Getting core study population..."
		for cluster in core:
			for gene in core[cluster]:
				if gene in annotations.keys():
					if annotations[gene]["GO"]:
						corepop.write("{0}\n".format(gene))
	with open("noncore_population.txt", "w") as noncorepop:
		print "Getting noncore study population..."
		for cluster in noncore:
			for gene in noncore[cluster]:
				if gene != "----------":
					if gene in annotations.keys():
						if annotations[gene]["GO"]:
							noncorepop.write("{0}\n".format(gene))
	sp.call(["find_enrichment.py", "--pval=0.05", "--method=fdr", "--obo", go_terms, "core_population.txt", "pangenome_population.txt", "pangenome_slim.txt", "--outfile=core_enrichment.tsv"])
	sp.call(["find_enrichment.py", "--pval=0.05", "--method=fdr", "--obo", go_terms, "noncore_population.txt", "pangenome_population.txt", "pangenome_slim.txt", "--outfile=noncore_enrichment.tsv"])


def main():
	print "Loading IPS annotations file..."
	anno_dict = generate_annotation_dict("interpro_results.txt")
	print "Loaded annotations.\nLoading TNT output..."
	#cluster_states = generate_mapping_dict("tnt.output")
	print "Loaded TNT output."
	get_complement_enrichments(core, noncore, anno_dict, "/Users/cmccarthy/Desktop/go.obo", "/Users/cmccarthy/Desktop/goslim_generic.obo")


if __name__ == "__main__":
	main()