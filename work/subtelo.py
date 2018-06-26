from csv import reader


sub_regions = {}
core = [0, 0]
acc = [0, 0]


with open("chr_sizes.txt") as sizes:
	for line in sizes.readlines():
		if line.startswith("ID"):
			pass
		else:
			fields = line.strip("\n").split("\t")
			first = round(int(fields[1]) * 0.1)
			second = round(int(fields[1]) * 0.9)
			sub_regions[fields[0]] = (first, second)

with open("phen_input.txt") as phen:
	chrs = reader(phen, delimiter="\t")
	for row in chrs:
		chr = row[0]
		start = int(row[1])
		stop = int(row[2])
		type = row[3]
		if sub_regions[chr][0] <= start <= sub_regions[chr][1]:
			if type == "3":
				core[0] = core[0] + 1
			else:
				acc[0] = acc[0] + 1
		else:
			if type == "3":
				core[1] = core[1] + 1
			else:
				acc[1] = acc[1] + 1
print core, acc
print sub_regions