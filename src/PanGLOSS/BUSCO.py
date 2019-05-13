import subprocess as sp


def RunBUSCO(tags):
    """

    """
    buscopath = "/bin/busco-master/scripts/run_BUSCO.py"

    for tag in tags:
        pset = "./sets" + tag + ".faa"
        output = pset + ".busco"
        cmd = ["python", buscopath, "-i", pset,  "-l", "saccharomycetales_odb9", "-o", output ,"-m", "prot"]
        sp.call(cmd)