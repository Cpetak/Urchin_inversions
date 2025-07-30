import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm

"""**Question**: do uniprotIDs mapped to the same LOC share similar GO terms?"""

locs=pd.read_csv("supp_files/idmapping_2024_08_21.tsv", sep="\t")
#locs=pd.read_csv("uniprot-pcangsd_all.tab", sep="\t", names=("LOC","uni","1","2","3","4","5","6"), header=None)
#locs = locs.iloc[1: , :]

locs = locs.rename(columns={'From': 'LOC'})
locs = locs.rename(columns={'Entry': 'uni'})

mapping=pd.read_csv("supp_files/GO_mapping_topGO", sep="\t", names=("uni","go"))

mydic={}
for index, row in locs.iterrows():
  loc = row[0]
  same_locs_list=locs[locs["LOC"]==loc]
  if len(same_locs_list) > 1:
    nlist=[]
    for i,l in same_locs_list.iterrows():
      goterms=list(mapping[mapping["uni"]==l[1]].go)
      if len(goterms) != 0:
        g=goterms[0].split(",")
        nlist.append(g)
    nlist=[item for sublist in nlist for item in sublist]
    nnlist=[x for x in nlist if nlist.count(x)==1]
    if len(nlist) != 0:
      if len(nnlist)/len(nlist) > 0:
        mydic[loc]=(len(nnlist)/len(nlist)) # Number of unique GO terms over all GO terms, over all different matches of LOC to UniprotID

plt.hist(mydic.values())

"""
I am just gonna choose the uniprotID with the most GO terms
"""

newdic={}
for loc in tqdm(locs.LOC.unique()):
  same_locs_list=locs[locs["LOC"]==loc]
  tempdic={}
  for i,l in same_locs_list.iterrows():
    goterms=list(mapping[mapping["uni"]==l[1]].go) # for each of the rows with the same loc, get GO term corresponding to uniprotID
    if len(goterms) != 0:
      g=goterms[0].split(",")
      tempdic[l[1]] = g

  if len(tempdic) == 0: # no uniprotIDs mapped to any goterms, select first uniprotID
    newdic[loc]=list(same_locs_list.uni)[0]

  else:
    m=0
    k=list(tempdic.keys())[0]
    for key in tempdic:
      num_go=len(tempdic[key])
      if num_go > m:
        m=num_go
        k=key
    newdic[loc]=k

textfile = open("all_locs_to_uniprotIDs.txt", "w")

for key, val in newdic.items():
    textfile.write(key + "," + val + "\n")
textfile.close()
