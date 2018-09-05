import sys
from Bio import SeqIO
from csv import reader
from random import choice

core = {}
acc = {}

idx = SeqIO.index("all.fasta", "fasta")
m = reader(open("new_matchtable.txt"), delimiter="\t")
for row in m:
	if len(filter(lambda x: x.startswith("-"), row[1:])) > 0:
		acc[row[0]] = row[1:]
	else:
		core[row[0]] = row[1:]

for key in acc:
	genes = filter(lambda x: "|" in x, acc[key])
	get = choice(genes)
	#print ">{0}\n{1}\n".format(idx[get].id, idx[get].seq)

for key in core:
	for gene in core[key]:
		print ">{0}\n{1}\n".format(idx[gene].id, idx[gene].seq)
