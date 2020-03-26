import os
from Bio import Entrez, SeqIO, AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Phylo import draw
import matplotlib.pyplot as plt

Entrez.email = "pr386097@students.mimuw.edu.pl"

def entrez_import():
    # write given database entries' dna codes to fasta files

    # respective database codes for: SARS-CoV-2, Bat SARS CoV HKU3-4,
    # Bat CoV, MERS-CoV, Murine Hepatitis Virus
    genome_db_codes = [1821109035, 292660135, 1189488873, 1386872249, 1268318180]
    spike_protein_db_codes = [1812779093, 292660137, 1023043585, 510937295, 1825483]

    print('\nImporting following DNA codes...\n')

    genome_dna = []
    for i in range(len(genome_db_codes)):
        genome_handle = Entrez.efetch(db="nucleotide",id=genome_db_codes[i],rettype="fasta")
        genome = SeqIO.read(genome_handle,"fasta")
        genome_dna.append(genome)

    # import additional Influenza A genome which is stored in 8 segments
    influenza_segments_codes = [310699718, 310699715, 310699713, 310699701, 310699708, 310699706, 310699703, 310699710]

    for i in range(len(influenza_segments_codes)):
        influenza_handle = Entrez.efetch(db="nucleotide", id=influenza_segments_codes[i], rettype="fasta")
        influenza_segment = SeqIO.read(influenza_handle, "fasta")
        if i == 0:
            genome = influenza_segment
            genome.id = 'Influenza A Virus'
        else:
            genome.seq = genome.seq + influenza_segment.seq

    genome_dna.append(genome)

    spike_protein_dna = []
    for i in range(len(spike_protein_db_codes)):
        spike_protein_handle = Entrez.efetch(db="protein",id=spike_protein_db_codes[i],rettype="fasta")
        spike_protein = SeqIO.read(spike_protein_handle,"fasta")
        spike_protein_dna.append(spike_protein)

    # write imported dna codes to fasta files
    SeqIO.write(genome_dna, 'data/genome_dna.fasta', 'fasta')
    SeqIO.write(spike_protein_dna, 'data/spike_protein_dna.fasta', 'fasta')

entrez_import()

# multiple sequence alignment using clustalo
print('Executing multiple sequence alignment...\n')
input_files = ['data/genome_dna.fasta', 'data/spike_protein_dna.fasta']
output_files = ['data/genome_dna_aligned.fasta', 'data/spike_protein_dna_aligned.fasta']

for i in range(len(input_files)):
    os.system('clustalo -i ' + input_files[i] + ' -o ' + output_files[i] + ' --force')

# import aligned dna sequences
spike_protein_dna = AlignIO.read("data/spike_protein_dna_aligned.fasta", "fasta")
genome_dna = AlignIO.read("data/genome_dna_aligned.fasta", "fasta")

# virus names for tree labels
viruses = ['SARS-CoV-2','Bat SARS CoV\n      HKU3-4', 'Bat CoV', 'MERS-CoV',
          'Murine Hepatitis\n          Virus', 'Influenza A']

# change virus codes to virus labels in the data
for i in range(0,5):
    genome_dna[i].id = viruses[i]
    spike_protein_dna[i].id = viruses[i]

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

    plt.savefig(output_file, dpi=200)
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
