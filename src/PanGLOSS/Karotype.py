from Bio import SeqIO

from glob import glob

def GenerateContigLengths(genomes):
    """
    """
    lengths = []
    for genome in glob("{0}/*.fna".format(genomes)):
        tag = genome.split("/")[1].split(".")[0]
        gen = []
        for seq in SeqIO.parse(genome, "fasta"):
            gen.append([seq.id, "1", str(len(seq.seq)), tag])
        lengths.append(gen)

    with open("{0}/lengths.txt".format(genomes), "w") as out:
        for gen in lengths:
            for contig in gen:
                row = "\t".join(contig)
                row += "\n"
                out.write(row)

