import os
import sys
import subprocess as sp
from Tools import ParseMatchtable, ClusterSizes


def GenerateSizeNumbers(matchtable):
    """

    :param matchtable:
    :return:
    """
    core, acc = ParseMatchtable(matchtable)

    with open("cluster_sizes.txt", "w") as outfile:
        outfile.write("Size\tCount\n")
        for component in [acc, core]:
            sizes = ClusterSizes(component)
            for size in sizes:
                outfile.write(str(size) + "\t" + str(sizes[size]) + "\n")

def GenerateBarChart(sizes):
    """

    :return:
    """
    barpath = os.path.dirname(os.path.realpath(sys.argv[0])) + "/BarChart.R"
    sp.call(["Rscript", barpath, sizes])
