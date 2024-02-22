#Import data and select the basics from the blast hit tables.
import pandas as pd
from Bio import Entrez
seqs1 = pd.read_csv('Sequences1_Hittable.csv',header=None)
seqs2 = pd.read_csv('Sequences2_Hittable.csv',header=None)
seqs2 = seqs2[[0,1,2,10]]
seqs1 = seqs1[[0,1,2,10]]
seqs1 = seqs1.groupby(0).head(5)
seqs2 = seqs2.groupby(0).head(5)
Entrez.email = "myemail@domain.com"
myl_seqs1,myl_seqs2 = [],[]

#Use Entrez.efetch to query the NT database for our accessions from the hittable.
#Parse the return from efetch, return scientific name as out.
for acc in seqs1[1]: 
    ret = Entrez.efetch(db = "nucleotide", id = acc, retmode="text")
    for ele in ret:
        if "taxname" in ele:
            start,end = ele.index("\""),ele.rindex(",")
            out = (ele[start+1:end-1])
            myl_seqs1.append(out)
            break
            
for acc in seqs2[1]: 
    ret = Entrez.efetch(db = "nucleotide", id = acc, retmode="text")
    for ele in ret:
        if "taxname" in ele:
            start,end = ele.index("\""),ele.rindex(",")
            out = (ele[start+1:end-1])
            myl_seqs2.append(out)
            break

#Rewrite NCBI accesions to Human readable name, columns names from numbers
seqs1[1]=myl_seqs1
seqs2[1]=myl_seqs2
seqs1=seqs1.set_axis(["Seq","Blast_Name","Pct.Identity","E.val"],axis=1)
seqs2=seqs2.set_axis(["Seq","Blast_Name","Pct.Identity","E.val"],axis=1)
seqs1.to_csv("seqs1_blastreformat.csv")
seqs2.to_csv("seqs2_blastreformat.csv")