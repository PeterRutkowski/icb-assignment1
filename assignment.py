# identifying ourselves to the Entrez database
from Bio import Entrez
from Bio import SeqIO
from Bio.Align.Applications import ClustalwCommandline
import numpy as np
Entrez.email = "piotrrutkowski97@gmail.com"

def genome_sequence_import(query):
    #searching for coronavirus sequence
    handle=Entrez.esearch(db="nuccore",term=query)
    rec=Entrez.read(handle)
    # now we have the search results in a dictionary, we can take the list of sequence IDs
    rec["IdList"]
    #print(rec["IdList"])
    print()
    for i in range(len(rec["IdList"])):
        #we can fetch the first one - the whole genome sequence
        rec_handle=Entrez.efetch(db="nucleotide",id=rec["IdList"][i],rettype="fasta")

        # now we have the handle, we need to read the fasta file from it

        ncov19=SeqIO.read(rec_handle,"fasta")

        print(ncov19.description, len(ncov19.seq), rec["IdList"][i])
    #print(ncov19.seq)

    # now we will do the same for the spike protein
    handle=Entrez.esearch(db="protein",term=query)
    spike_rec=Entrez.read(handle)
    print()
    for i in range(len(spike_rec["IdList"])):
        spike_handle=Entrez.efetch(db="protein",id=spike_rec["IdList"][i],rettype="fasta")
        spike_ncov19=SeqIO.read(spike_handle,"fasta")

        print(spike_ncov19.description, len(spike_ncov19.seq), spike_rec["IdList"][i])
    #print(spike_ncov19.seq)
    print()
    return [ncov19.seq, spike_ncov19.seq]

def genome_sequence_import2(genomes, spikes):
    #searching for coronavirus sequence
    print()
    for i in range(len(genomes)):
        handle=Entrez.efetch(db="nucleotide",id=genomes[i],rettype="fasta")
        genome=SeqIO.read(handle,"fasta")
        print(genome.description, len(genome.seq))
        #print(genome.seq)

    print()
    for i in range(len(spikes)):
        spike_handle=Entrez.efetch(db="protein",id=spikes[i],rettype="fasta")
        spike=SeqIO.read(spike_handle,"fasta")

        print(spike.description, len(spike.seq))
    #print(spike_ncov19.seq)
    print()
    #return [ncov19.seq, spike_ncov19.seq]

query_list = ['2019-nCov','SARS 2002', 'bat coronavirus', 'MERS', 'influenza A', 'hepatitis A']
genomes = [1821109035, 292660135, 1180422623, 1386872249, 1820140354, 1812620187]
spikes = [1812779093, 292660137, 1179780473, 1386872252, 1820506005, 1812620188] # last two may be wrong

genetic_codes = []

#genetic_codes.append(genome_sequence_import('hepatitis a'))

genome_sequence_import2(genomes, spikes)

#for query in query_list:
#    genetic_codes.append(genome_sequence_import(query))


