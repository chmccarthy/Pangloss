import cStringIO

from Bio import AlignIO, SeqIO, SeqRecord
from Bio.Phylo.PAML import codeml

from Tools import StringMUSCLE, Untranslate


def TranslateCDS():
    nucl = SeqIO.parse("test.faa", "fasta")
    tn = []
    for seq in nucl:
        t = seq.translate()
        t.id = seq.id
        tn.append(t)

    return nucl, tn


def MUSCLEAlign(seqs):
    output = StringMUSCLE(seqs)
    return AlignIO.parse(cStringIO.StringIO(output), "fasta")


def PutGaps(alignment, nucl):
    nucl = SeqIO.index("test.faa", "fasta")
    nucl_aln = []
    for aln in alignment:
        for seq in aln._records:
            nseq = nucl[seq.id].seq
            aseq = seq.seq
            unseq = Untranslate(aseq, nseq)
            nucl_aln.append(SeqRecord.SeqRecord(unseq, id=seq.id))
    return nucl_aln
