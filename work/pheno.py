from csv import reader
from Bio import SeqIO
from itertools import chain, izip_longest, tee

def flatten(iterable):
	"""
	Flatten a list of lists, essential for ClusterClean and GapFinder.
	
	Taken from the Python Standard Library.
	"""
	return list(chain.from_iterable(iterable))

def matchdict(matchtable, core=True):
	dict = {}
	d2 = {}
	for row in matchtable:
		if core:
			if "----------" in row:
				d2[row[0]] = row[1:]
			else:
				dict[row[0]] = row[1:]  # Populating our core dict.
	return dict, d2


phen = reader(open("S288C_attributes.txt"), delimiter="\t")
core, non = matchdict(reader(open("new_matchtable.txt"), delimiter="\t"))

with open("phen_input.txt", "w") as outfile:
	chrn = 0
	chr = ""
	for row in phen:
		if row[0] != chr:
			chr = row[0]
			chrn = chrn + 1
		if row[1] in flatten(core.values()):
			col = "3"
			outfile.write("\t".join([row[0], str(row[2]), str(row[3]), col]))
			outfile.write("\n")
		elif row[1] in flatten(non.values()):
			col = "2"
			outfile.write("\t".join([row[0], str(row[2]), str(row[3]), col]))
			outfile.write("\n")
		else:
			print row[1], "not in pg"
