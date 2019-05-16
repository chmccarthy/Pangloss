# -*- coding: utf-8 -*-
"""
PAML: Module for handling yn00 selection analysis (and maybe CodeML in the future), if enabled by user.
"""

import cStringIO

from Bio import AlignIO, SeqIO
from Bio.Phylo.PAML import yn00
from glob import glob

from Tools import StringMUSCLE, Untranslate


def TranslateCDS(seqs):
    """
    Given a nucleotide FASTA file, translate each sequence.
    """
    tn = []
    for seq in SeqIO.parse(seqs, "fasta"):
        t = seq.translate()
        t.id = seq.id
        tn.append(t)

    return tn


def MUSCLEAlign(seqs):
    """
    Align translated nucleotides in StringMUSCLE, return parsed alignment.
    """
    output = StringMUSCLE(seqs)
    return AlignIO.parse(cStringIO.StringIO(output), "fasta")


def PutGaps(alignment, cluster):
    """
    Given amino acid alignment, generate a matching nucleotide alignment which respects codons. Note: the alignment
    file that this function outputs contains sequence IDs without location data, as yn00 limits sequence IDs to 30
    characters in length in its output files and more importantly Biopython's yn00 output parser has problems
    parsing sequence IDs that contains dots (which will be fixed in 1.40) and/or sequence IDs that contains dots
    and are >30 characters in length.
    """
    nucl = SeqIO.index(cluster, "fasta")
    nucl_aln = ""
    for aln in alignment:
        for seq in aln._records:
            nseq = nucl[seq.id].seq
            aseq = seq.seq
            unseq = Untranslate(aseq, nseq)
            unseq.id = seq.id.split("|")[0]
            nucl_aln += (">{0}\n{1}\n".format(unseq.id, unseq.seq))

    fas_aln = AlignIO.read(cStringIO.StringIO(nucl_aln), "fasta")
    AlignIO.write(fas_aln, "{0}.aln".format(cluster), "phylip-sequential")

    return "{0}.aln".format(cluster)


def RunYn00(alignment):
    """
    Run yn00 on untranslated alignment with default parameters, output to file and then parse file for
    the results of Yang and Nielsen (2000) analysis.
    """
    yn = yn00.Yn00(alignment=alignment, out_file="{0}.yn00".format(alignment))
    yn.set_options(verbose=0, icode=0, weighting=0, commonf3x4=0)
    yn.run(ctl_file=None, command="yn00", parse=True)


def SummarizeYn00():
    """
    Summarize yn00 results for core and accessory genomes.
    """
    results = {}
    clusters = glob("./panoct/clusters/core/fna/Core*.fna.aln.yn00") + glob("./panoct/clusters/acc/fna/Acc*.fna.aln.yn00")
    for cluster in clusters:
        cl_number = cluster.split("_")[1].split(".")[0]
        results[cl_number] = {"Component": "", "Kappa": 0, "Alignments with omega > 1":}
        yn = yn00.read(cluster)
        kappa = 0
        for gene in yn:
            for alignment in yn[gene]:
                pairwise_result = yn[gene][alignment][yn00]

        print yn.items()



