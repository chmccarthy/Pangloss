from __future__ import division

"""
Module for generating bar and ring charts.
"""
import os
import sys
import subprocess as sp
from Tools import ParseMatchtable, ClusterSizes

def GenerateRingChart(matchtable):
    """
    Generate ring chart using RingChart.R.
    """
    ringchart = os.path.dirname(os.path.realpath(sys.argv[0])) + "/RingChart.R"
    core, acc = ParseMatchtable(matchtable)
    core_sizes = ClusterSizes(core)
    acc_sizes = ClusterSizes(acc)
    core_total = sum(core_sizes.values())
    acc_total = sum(acc_sizes.values())
    sp.call(["Rscript", ringchart, str(core_total), str(acc_total)])


def GenerateSizeNumbers(matchtable):
    """
    Generate cluster sizes file for bar and ring charts.
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
    Generate bar chart using BarChart.R.
    """
    barpath = os.path.dirname(os.path.realpath(sys.argv[0])) + "/BarChart.R"
    sp.call(["Rscript", barpath, sizes])