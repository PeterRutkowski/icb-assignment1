import os
from Bio import Entrez, SeqIO, AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Phylo import draw
import matplotlib.pyplot as plt

Entrez.email = "pr386097@students.mimuw.edu.pl"

def import_dna_from_db(genome_db_codes, spike_protein_db_codes):
    # write given database entries' dna codes to fasta files
    print('\nImporting following genomes:\n')

    genome_dna = []
    for i in range(len(genome_db_codes)):
        genome_handle = Entrez.efetch(db="nucleotide",id=genome_db_codes[i],rettype="fasta")
        genome = SeqIO.read(genome_handle,"fasta")
        genome_dna.append(genome)
        print(genome.description)

    print('\nImporting following spike proteins:\n')

    spike_protein_dna = []
    for i in range(len(spike_protein_db_codes)):
        spike_protein_handle = Entrez.efetch(db="protein",id=spike_protein_db_codes[i],rettype="fasta")
        spike_protein = SeqIO.read(spike_protein_handle,"fasta")
        spike_protein_dna.append(spike_protein)
        print(spike_protein.description)

    # write imported dna codes to fasta files
    SeqIO.write(genome_dna, 'data/genome_dna.fasta', 'fasta')
    SeqIO.write(spike_protein_dna, 'data/spike_protein_dna.fasta', 'fasta')

# respective database codes for: SARS-CoV-2, Bat SARS CoV HKU3-4,
# Bat CoV, MERS-CoV, Influenza A, Duck hepatitis A DHAV-3
genome_db_codes = [1821109035, 292660135, 1180422623, 1386872249, 1820140354, 1812620187]
spike_protein_db_codes = [1812779093, 292660137, 1179780473, 1386872252, 1820506005, 1812620188]

import_dna_from_db(genome_db_codes, spike_protein_db_codes)

# multiple sequence alignment using clustalo
print('\nExecuting multiple sequence alignment...\n')
input_files = ['data/genome_dna.fasta', 'data/spike_protein_dna.fasta']
output_files = ['data/genome_dna_aligned.fasta', 'data/spike_protein_dna_aligned.fasta']

for i in range(len(input_files)):
    os.system('clustalo -i ' + input_files[i] + ' -o ' + output_files[i] + ' --force')

# import aligned dna sequences
spike_protein_dna = AlignIO.read("data/spike_protein_dna_aligned.fasta", "fasta")
genome_dna = AlignIO.read("data/genome_dna_aligned.fasta", "fasta")

# virus names for tree labels
viruses = ['SARS-CoV-2','Bat SARS CoV\n      HKU3-4', 'Bat CoV', 'MERS-CoV',
          'Influenza A', 'Duck hepatitis A\n        DHAV-3']

# change virus codes to virus labels in the data
for i in range(len(spike_protein_dna)):
    spike_protein_dna[i].id = viruses[i]
    genome_dna[i].id = viruses[i]

# calculate distance matrices
print('Calculating distance matrices...\n')
calculator = DistanceCalculator('blosum62')
spike_protein_dm = calculator.get_distance(spike_protein_dna)
genome_dm = calculator.get_distance(genome_dna)

# build trees
print('Building trees...\n')
constructor = DistanceTreeConstructor()
spike_protein_nj_tree = constructor.nj(spike_protein_dm)
spike_protein_upgma_tree = constructor.upgma(spike_protein_dm)
genome_nj_tree = constructor.nj(genome_dm)
genome_upgma_tree = constructor.upgma(genome_dm)

def plot_tree(tree, title, output_file):
    # build the plot

    # clear names of nodes that are not leaves
    for clade in tree.find_clades():
        if clade.name[0:5] == 'Inner':
            clade.name = ''

    fig, ax = plt.subplots()
    axes = fig.add_subplot(1, 1, 1)
    fig = draw(tree, axes = axes, do_show=0)

    ax.get_xaxis().set_ticks([])
    ax.get_yaxis().set_ticks([])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)

    plt.gcf().subplots_adjust(right=0.99, left=0.01, top=0.9, bottom=0.1)
    plt.suptitle(title)
    plt.yticks([])
    plt.ylabel('')

    plt.savefig(output_file, dpi=450)
    plt.close(fig)

# plot trees
print('Saving trees...\n')
trees = [genome_nj_tree, genome_upgma_tree, spike_protein_nj_tree, spike_protein_upgma_tree]
plot_titles = ['genome phylogenetic tree (NJ algorithm)', 'genome phylogenetic tree (UPGMA algorithm)',
               'spike protein phylogenetic tree (NJ algorithm)', 'spike protein phylogenetic tree (UPGMA algorithm)']
plot_output_files = ['genome_tree_nj.png', 'genome_tree_upgma.png', 'spike_protein_tree_nj.png',
                     'spike_protein_tree_upgma.png']

for i in range(len(trees)):
    plot_tree(trees[i], 'Virus ' + plot_titles[i], 'plots/' + plot_output_files[i])
