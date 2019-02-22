import cStringIO

from Bio import AlignIO, SeqIO
from Bio.Phylo.PAML import yn00
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


def PutGaps(seqs, alignment):
    """
    Given amino acid alignment,
    """
    nucl = SeqIO.index(seqs, "fasta")
    nucl_aln = ""
    for aln in alignment:
        for seq in aln._records:
            nseq = nucl[seq.id].seq
            aseq = seq.seq
            unseq = Untranslate(aseq, nseq)
            unseq.id = seq.id
            nucl_aln += (">{0}\n{1}\n".format(unseq.id, unseq.seq))

    with open("{0}.aln".format(seqs), "w") as out:
        for seq in nucl_aln:
            out.write(seq)

    return seqs


def RunYn00(alignment):
    yn = yn00.Yn00(alignment=alignment, out_file="{0}.yn00".format(alignment))
    yn.set_options(verbose=0, icode=0, weighting=0, commonf3x4=0)
    yn.run(ctl_file=None, command="yn00", parse=True)



