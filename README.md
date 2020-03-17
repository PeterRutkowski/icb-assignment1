# ICB1920AS1
MIMUW 19/20 Introduction to Computational Biology Assignment 1

As we are spending our time in quarantine, I thought we might try and consider reconstructing some phylogenetic trees that will tell us something about the evolutionary history of the coronavirus 2019-nCov. In particular, we will try and replicate (not completely, but the general idea) the panels b and c from Figure 2 in this paper:Identification_of_a_novel_coronavirus . The idea is to create and compare two phylogenetic trees created from the whole genome sequence of the coronavirus and from its spike protein (you might read more about spike proteins in this article, back from the previous SARS epidemic).

Our goal will be to write a script that creates these two plots using the following steps:

Downloading the sequences from Entrez using Bio.Entrez (2 pts)
Creating multiple alignments of both DNA (whole genome) and protein (spike protein) sequences (3 pts)
Creating phylogenetic trees for both protein and genome alignments using the UPGMA and neighbor joining algorithms (3pts)
Visualizing the trees (1 pt)
Now let us describe all these steps in more detail

1. We need to download whole genome sequences and spike protein sequences from the Genbank Entrez databases (“nucleotide” and “protein” databases, respectively), together with some other viruses. We will only look at 6 species:

2019-nCov virus (current coronavirus)
SARS – another coronavirus back from 2002
bat coronavirus
MERS – middle eastern respiratory syndrome coronavirus
influenza A -  flu virus, to see how far are coronaviruses from the flu
hepatitis A virus – less related virus
You can use the Bio.Entrez interface. An example of how this works can be seen in a jupyter notebook here (please remember to use your e-mail and not mine in your code…). We will need to do this for all 6 species. You should provide a single script that does that for all species, using different query terms for each genome and protein. You will probably need some trial and error work here, but I want to see just the end result – a script that works.

2. Here we will want to run the ClustalW algorithm on the downloaded sequences. You will need to have the clustalw executable installed (it is the package clustalw in ubuntu). Your script should use the Bio.Align.Applications.Clustalw interface. in the result you should get two multiple sequence alignments (for genomes and for spike proteins).

3. Here you should use the two multiple alignments to create distance matrices and from these distance matrices to create phylogenetic trees (both with upgma and nj methods).

4. You should create images of the 4 trees (upgma or nj x protein or genome) using the draw function.

A full solution is a script (in a .py file sent as a real attachment -not some usos attachment, not an ipynb, not a link) together with the images for part 4 (in png format). All of this should be sent to me (bartek[at]mimuw.edu.pl) with the tag [WBO] in the e-mail subject by April 6th.
