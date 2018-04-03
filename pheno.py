from csv import reader
from Bio import SeqIO
from itertools import chain, izip_longest, tee

def flatten(iterable):
	"""
	Flatten a list of lists, essential for ClusterClean and GapFinder.
	
	Taken from the Python Standard Library.
	"""
	return list(chain.from_iterable(iterable))

def matchdict(matchtable):
	core = {}
	noncore = {}
	for row in matchtable:
		if "----------" in row:
			noncore[row[0]] = row[1:]  # Populating our initial noncore dict.
		else:
			core[row[0]] = row[1:]  # Populating our core dict.
	return core, noncore

phen = reader(open("AF293_phenogram.txt"), delimiter="\t")
core, noncore = matchdict(reader(open("matchtable.txt"), delimiter="\t"))

with open("phen_input.txt", "w") as outfile:
	for row in phen:
		if row[1] in flatten(core.values()):
			col = "3"
		elif row[1] in flatten(noncore.values()):
			col = "2"
		outfile.write("\t".join([row[0], row[2], row[3], col]))
		outfile.write("\n")
