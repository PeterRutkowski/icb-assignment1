from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio.Phylo import draw_ascii, draw
import numpy as np

labels = ['SARS-CoV-2','SARS', 'bat coronavirus', 'MERS', 'influenza A', 'hepatitis A']

spike_proteins = AlignIO.read("s_out.fasta", "fasta")

ids = []
for i in range(len(spike_proteins)):
    spike_proteins[i].id = labels[i]



calculator = DistanceCalculator('blosum62')
dm = calculator.get_distance(spike_proteins)
constructor = DistanceTreeConstructor()
njtree = constructor.nj(dm)
print('\n\n')
draw_ascii(njtree)