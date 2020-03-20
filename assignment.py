# identifying ourselves to the Entrez database
from Bio import Entrez
from Bio import SeqIO
import numpy as np
Entrez.email = "piotrrutkowski97@gmail.com"

def genome_sequence_import(query):
    #searching for coronavirus sequence
    handle=Entrez.esearch(db="nuccore",term=query)
    rec=Entrez.read(handle)
    # now we have the search results in a dictionary, we can take the list of sequence IDs
    rec["IdList"]
    print(rec["IdList"])

    #we can fetch the first one - the whole genome sequence
    rec_handle=Entrez.efetch(db="nucleotide",id=rec["IdList"][2],rettype="fasta")

    # now we have the handle, we need to read the fasta file from it

    ncov19=SeqIO.read(rec_handle,"fasta")

    print(ncov19.description)
    #print(ncov19.seq)

    # now we will do the same for the spike protein
    handle=Entrez.esearch(db="protein",term=query)
    spike_rec=Entrez.read(handle)
    spike_handle=Entrez.efetch(db="protein",id=spike_rec["IdList"][0],rettype="fasta")
    spike_ncov19=SeqIO.read(spike_handle,"fasta")

    print(spike_ncov19.description)
    #print(spike_ncov19.seq)

    return ncov19.seq, spike_ncov19.seq

a, b = genome_sequence_import('bat coronavirus')
#print(a)
#print()
#print(b)

