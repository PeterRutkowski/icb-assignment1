from Bio import Entrez
import numpy as np
Entrez.email = "piotrrutkowski97@gmail.com"

def genome_search(query, database, query_results=10):
    handle = Entrez.esearch(db=database, retmax=query_results, term=query)
    records = Entrez.read(handle)
    handle.close()
    #print(records['IdList'])
    
    id_records_text = ''
    for i in range(len(records['IdList'])):
        if len(id_records_text) > 0:
            id_records_text = id_records_text + ',' + records['IdList'][i]
        else:
            id_records_text = records['IdList'][i]
    #print(id_records_text)

    handle = Entrez.esummary(db=database, id=id_records_text, retmode="xml")
    records = Entrez.parse(handle)
    for record in records:
        print(record['Title'])
        print()
    handle.close()

genome_search(query='influenza A', database='nucleotide', query_results=15)