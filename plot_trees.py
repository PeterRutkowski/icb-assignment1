from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio.Phylo import draw
import matplotlib.pyplot as plt

# import aligned dna sequences
spike_protein_dna = AlignIO.read("spike_protein_dna_aligned.fasta", "fasta")
genome_dna = AlignIO.read("genome_dna_aligned.fasta", "fasta")

# virus names for tree labels
viruses = ['SARS-CoV-2','Bat SARS CoV\n      HKU3-4', 'Bat CoV', 'MERS-CoV',
          'Influenza A', 'Duck hepatitis A\n        DHAV-3']

# change virus codes to virus names in the data
for i in range(len(spike_protein_dna)):
    spike_protein_dna[i].id = viruses[i]
    genome_dna[i].id = viruses[i]

# calculate distance matrices
calculator = DistanceCalculator('blosum62')
spike_protein_dm = calculator.get_distance(spike_protein_dna)
genome_dm = calculator.get_distance(genome_dna)

# build trees
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
trees = [genome_nj_tree, genome_upgma_tree, spike_protein_nj_tree, spike_protein_upgma_tree]
plot_titles = ['genome phylogenetic tree (NJ algorithm)', 'genome phylogenetic tree (UPGMA algorithm)',
               'spike protein phylogenetic tree (NJ algorithm)', 'spike protein phylogenetic tree (UPGMA algorithm)']
plot_output_files = ['genome_tree_nj.png', 'genome_tree_upgma.png', 'spike_protein_tree_nj.png',
                     'spike_protein_tree_upgma.png']

for i in range(len(trees)):
    plot_tree(trees[i], 'Virus ' + plot_titles[i], plot_output_files[i])
