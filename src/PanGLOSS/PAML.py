import cStringIO

from Bio import AlignIO, SeqIO
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





