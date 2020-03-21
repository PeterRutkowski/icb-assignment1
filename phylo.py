from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio.Phylo import draw_ascii, draw
import matplotlib.pyplot as plt

labels = ['SARS-CoV-2','Bat SARS CoV HKU3-4', 'Bat CoV', 'MERS-CoV',
          'Influenza A', 'Duck hepatitis A\n        DHAV-3']

spike_protein_dna = AlignIO.read("s_out.fasta", "fasta")
genome_dna = AlignIO.read("g_out.fasta", "fasta")

for i in range(len(spike_protein_dna)):
    spike_protein_dna[i].id = labels[i]
    genome_dna[i].id = labels[i]

calculator = DistanceCalculator('blosum62')

spike_protein_dm = calculator.get_distance(spike_protein_dna)
genome_dm = calculator.get_distance(genome_dna)

constructor = DistanceTreeConstructor()

spike_protein_nj_tree = constructor.nj(spike_protein_dm)
spike_protein_upgma_tree = constructor.upgma(spike_protein_dm)
genome_nj_tree = constructor.nj(genome_dm)
genome_upgma_tree = constructor.upgma(genome_dm)

def plot_tree(tree, title, output_file):
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

plot_tree(spike_protein_nj_tree, 'Spike protein NJ tree','spike_protein_nj_tree.png')
