from Bio.Align.Applications import ClustalwCommandline
cline = ClustalwCommandline("clustalw2", infile="test.fasta")
from Bio import AlignIO
align = AlignIO.read("data/opuntia.aln", "clustal")
print(align)