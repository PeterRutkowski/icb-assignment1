from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio.Phylo import draw_ascii, draw
import matplotlib.pyplot as plt

labels = ['SARS-CoV-2','Bat SARS CoV HKU3-4', 'Bat CoV', 'MERS-CoV', 'Influenza A', 'Duck hepatitis A DHAV-3']

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

#print('\n\n')
#print('Spike protein NJ tree')
#draw_ascii(spike_protein_nj_tree)
#print('\n\n')
#print('Genome NJ tree')
#draw_ascii(genome_nj_tree)
#print('\n\n')
#print('Spike protein UPGMA tree')
#draw_ascii(spike_protein_upgma_tree)
#print('\n\n')
#print('Genome UPGMA tree')
#draw_ascii(genome_upgma_tree)
#print('\n\n')

def plot_tree(tree, title, output_file):
    fig, ax = plt.subplots()

    axes = fig.add_subplot(1, 1, 1)
    fig = draw(tree, axes = axes, do_show=0)
    plt.tight_layout()

    ax.get_xaxis().set_ticks([])
    ax.get_yaxis().set_ticks([])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)

    plt.suptitle(title)
    plt.axis('off')
    plt.savefig(output_file, dpi=450)

    plt.close(fig)

plot_tree(spike_protein_nj_tree, 'Spike protein NJ tree','spike_protein_nj_tree.png')
