from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio.Phylo import draw_ascii, draw
import numpy as np

labels = ['SARS-CoV-2','SARS', 'bat coronavirus', 'MERS', 'influenza A', 'hepatitis A']

spike_protein_dna = AlignIO.read("s_out.fasta", "fasta")
genome_dna = AlignIO.read("g_out.fasta", "fasta")

for i in range(len(spike_protein_dna)):
    spike_protein_dna[i].id = labels[i]
    genome_dna[i].id = labels[i]

calculator = DistanceCalculator('blosum62')

spike_protein_dm = calculator.get_distance(spike_protein_dna)
genome_dm = calculator.get_distance(genome_dna)

constructor = DistanceTreeConstructor()

spike_protein_njtree = constructor.nj(spike_protein_dm)
genome_njtree = constructor.nj(genome_dm)

print('\n\n')
print('Spike protein NJ tree')
draw_ascii(spike_protein_njtree)
print('\n\n')
print('Genome NJ tree')
draw_ascii(genome_njtree)
print('\n\n')