import cStringIO

from Bio import AlignIO, SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.PAML import yn00

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
    nucl_aln = ""
    for aln in alignment:
        for seq in aln._records:
            nseq = nucl[seq.id].seq
            aseq = seq.seq
            unseq = Untranslate(aseq, nseq)
            unseq.id = seq.id
            nucl_aln += (">{0}\n{1}\n".format(unseq.id, unseq.seq))
    with open("test.phylip", "w") as h:
        for seq in nucl_aln:
            h.write(seq)


def RunYn00(alignment):
    yn = yn00.Yn00(alignment=alignment, out_file="test.yn00")
    yn.set_options(verbose=0, icode=0, weighting=0, commonf3x4=0)
    results = yn.run(ctl_file=None, command="yn00", parse=True)
    print results

