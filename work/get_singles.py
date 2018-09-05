import sys
from Bio import SeqIO
from csv import reader

m = reader(open("new_matchtable.txt"), delimiter="\t")
idx = SeqIO.index(sys.argv[1], "fasta")
l = []


for row in m:
	if len(filter(lambda x: x.startswith("-"), row[1:])) == 11:
		if "|" in row[1]:
 			l.append(row[1])

for i in l:
	print ">{0}\n{1}\n".format(idx[i].id, idx[i].seq)
