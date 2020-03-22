import os
from Bio import Entrez
from Bio import SeqIO

Entrez.email = "piotrrutkowski97@gmail.com"

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
    SeqIO.write(genome_dna, 'genome_dna.fasta', 'fasta')
    SeqIO.write(spike_protein_dna, 'spike_protein_dna.fasta', 'fasta')

# respective database codes for: SARS-CoV-2, Bat SARS CoV HKU3-4,
# Bat CoV, MERS-CoV, Influenza A, Duck hepatitis A DHAV-3
genome_db_codes = [1821109035, 292660135, 1180422623, 1386872249, 1820140354, 1812620187]
spike_protein_db_codes = [1812779093, 292660137, 1179780473, 1386872252, 1820506005, 1812620188]

import_dna_from_db(genome_db_codes, spike_protein_db_codes)

# multiple sequence alignment using clustalo
input_files = ['genome_dna.fasta', 'spike_protein_dna.fasta']
output_files = ['genome_multiple_seq_alignment.fasta', 'spike_protein_multiple_seq_alignment.fasta']

for i in range(len(input_files)):
    os.system('clustalo -i ' + input_files[i] + ' -o ' + output_files[i] + ' --force')