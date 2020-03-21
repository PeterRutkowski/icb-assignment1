from Bio.Align.Applications import ClustalwCommandline
in_file = "test.fasta"
clustalw_cline = ClustalwCommandline("clustalw2", infile=in_file, align=1)
print(clustalw_cline)