import re
import numpy as np
import subprocess as sp


def generate_mapping_dict(tnt_output):
    node_match = re.compile(r"([ ]{3,4})[\S]")
    char_match = re.compile(r"([ ]{6})[\S]")
    states = {}
    apo = False
    current_node = ""
    with open(tnt_output) as infile:
        for line in infile.readlines():
            if line.startswith("Tree 0 :"):
                apo = True
            if node_match.match(line) is not None:
                current_node = line.lstrip().split(":")[0].strip()
                states[current_node] = []
            if apo:
                if char_match.match(line) is not None:
                    if not line.lstrip().startswith("No "):
                        print line
                        data = line.lstrip().split(":")
                        cluster = int(data[0].split(" ")[1]) + 1
                        if data[1].lstrip().strip() == "1 --> 0":
                            gain = False
                        elif data[1].lstrip().strip() == "0 --> 1":
                            gain = True
                        states[current_node].append((str(cluster), gain))
    return states


def main():
    match_matrix = np.loadtxt(open("noncore_matchtable.txt"), dtype=int, delimiter="\t")
    outgroup_row = [0 for i in range(0, len(match_matrix))]
    outgroup_array = np.array([outgroup_row, outgroup_row, outgroup_row])

    tnt_matrix = np.concatenate((match_matrix.T, outgroup_array), axis = 0)
    n_tax = tnt_matrix.shape[0]
    n_char = tnt_matrix.shape[1]

    with open("tnt.input", "w") as outfile:
        outfile.write("#NEXUS \n\nBegin data;\n\tDimensions ntax={0} nchar={1};\n\tFormat datatype=standard gap=- missing=? matchchar=.;\n\tMatrix\n".format(n_tax, n_char))
        for row in tnt_matrix:
            outfile.write("\t\tSpecies name\n\t\t{0}\n".format("".join([str(i) for i in row.tolist()])))
        outfile.write("\t;\nEnd;\n")

    with open("tnt.run", "w") as outrun:
        outrun.write("mxram 50000;\nproc tnt.input;\nproc tnt.tree;\nlog tnt.output;\nnaked -;\ntplot;\napo -;\nquit")

    sp.call(["tnt.command", "run", "tnt.run", ","])
    tnt_dict = generate_mapping_dict("tnt.output")
    for key in tnt_dict:
        print key, len(filter(lambda x: x[1] == True, tnt_dict[key])), len(filter(lambda x: x[1] == False, tnt_dict[key]))


