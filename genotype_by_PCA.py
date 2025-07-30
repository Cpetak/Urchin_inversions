#Import packages

from collections import Counter
import matplotlib.pyplot as plt
from collections import defaultdict
import numpy as np
import pandas as pd
import argparse
from scipy import stats
from scipy.stats import chisquare
import csv

from sklearn.cluster import KMeans
from sklearn.cluster import AgglomerativeClustering
from mpl_toolkits.basemap import Basemap as Basemap


#K-means clustering
def get_KMeans(values, k):
    X = np.array(values).reshape(-1, 1)

    # Perform K-means clustering
    kmeans = KMeans(n_clusters=k,n_init="auto",random_state=42)
    kmeans.fit(X)

    cluster_labels = kmeans.labels_
    clusters = [[] for _ in range(k)]
    for value, label in zip(values, cluster_labels):
        clusters[label].append(value)

    return sorted(clusters, key=lambda x: x[0])

#Agglomerative clustering
def get_agglo(values, k):
  X = np.array(values).reshape(-1, 1)

  agg_clustering = AgglomerativeClustering(n_clusters=k, linkage='ward')
  agg_clustering.fit(X)

  labels = agg_clustering.labels_

  return labels

#Agglomerative clustering, 2D
def get_agglo_2d(x,y,k):
  X = np.array(x).reshape(-1, 1)
  Y = np.array(y).reshape(-1, 1)

  agg_clustering = AgglomerativeClustering(n_clusters=k, linkage='ward')
  agg_clustering.fit(np.concatenate((X, Y), axis=1))

  labels = agg_clustering.labels_

  return labels

#Elbow plot
def find_elbow(x,y,dt):
  X = np.array(x).reshape(-1, 1)
  Y = np.array(y).reshape(-1, 1)

  agg_clustering = AgglomerativeClustering(n_clusters=None, linkage='ward',distance_threshold=dt)
  agg_clustering.fit(np.concatenate((X, Y), axis=1))

  return agg_clustering.n_clusters_


#Initialize variables
parser = argparse.ArgumentParser(description="Description of your script")
parser.add_argument("df1", type=str, help="PC1 csv")
parser.add_argument("df2", type=str, help="PC2 csv")
parser.add_argument("k", type=int, help="Number of clusters")
parser.add_argument("perc_explained", type=str, help="File for percent explained for PC 1 and 2")

args = parser.parse_args()

df1=pd.read_csv(args.df1)
df2=pd.read_csv(args.df2)
dim1=df1["dim1"].tolist()
dim2=df2["dim2"].tolist()

pe=pd.read_csv(args.perc_explained,names=["pes"])
pes_list=pe.pes.to_list()

#Population labels
dim=np.array(list(zip(dim1,dim2)))
pops=["BOD", "CAP", "FOG", "KIB", "LOM", "SAN", "TER"]
tpop=[]
for p in pops:
  if p != "CAP" and p!= "FOG":
    tpop.append([p] * 20)
  else:
    if p == "CAP":
      tpop.append([p] * 19)
    if p == "FOG":
      tpop.append([p] * 18)

pops=np.array([x for xs in tpop for x in xs])
ind_id=np.array(list(range(1,len(pops)+1)))

NUM_CLUST=args.k

#Make elbow plot, a.k.a. how many groups are there?
f_clusts=[]
for i in np.arange(0.05,1,0.01):
  f_clusts.append(find_elbow(dim1,dim2,i))

plt.plot(np.arange(0.05,1,0.01), f_clusts)
plt.axhline(3,color="red",linestyle=":")
plt.axhline(5,color="red",linestyle=":")
plt.ylim(0,max(f_clusts)+1)
plt.yticks(np.arange(0, max(f_clusts), step=1))
plt.ylabel("Number of clusters")
plt.xlabel("The linkage distance threshold (above which clusters will not be merged)")
plt.savefig(args.df1[:-9]+"_elbow.pdf")


#Cluster individuals based on how many groups there are
if NUM_CLUST == 3:
  clusts=get_agglo(dim1, NUM_CLUST)
else:
  clusts=get_agglo_2d(dim1,dim2,NUM_CLUST)

#Assign colors based on clustering
#order of colors: homoq, homop, hetero
all_cols=["C0","C1","C2",'red', 'blue', 'orange', 'cyan', 'magenta', 'brown', 'lime']
p_colors= all_cols[:NUM_CLUST]

f, axs = plt.subplots(figsize=(10, 10))

coloring_dic={}
for c in range(NUM_CLUST):
  t=np.where(np.array(clusts) == c)[0]
  x=[i[0] for i in dim[t]]
  if min(x) == min(dim1):
    coloring_dic["left"] = c
  elif max(x) == max(dim1):
    coloring_dic["right"] = c
  else:
    coloring_dic["middle"] = c

coloring_dic2={}
coloring_dic2[coloring_dic["middle"]] = all_cols[2]
if Counter(clusts)[coloring_dic["left"]] > Counter(clusts)[coloring_dic["right"]]:
  coloring_dic2[coloring_dic["left"]] = all_cols[1]
  coloring_dic2[coloring_dic["right"]] = all_cols[0]
else:
  coloring_dic2[coloring_dic["left"]] = all_cols[0]
  coloring_dic2[coloring_dic["right"]] = all_cols[1]

pop_dic={}
id_dic={}
lens_for_props = []
xs=[]
ys=[]
cs=[]
for c in range(NUM_CLUST):
  t=np.where(np.array(clusts) == c)[0]
  print(len(t))
  lens_for_props.append(len(t))
  x=[i[0] for i in dim[t]]
  y=[i[1] for i in dim[t]]
  plt.scatter(x,y,c=coloring_dic2[c],s=100)

  #saving to file for separate recreation of the figure:
  for data_idx in range(len(x)):
    xs.append(x[data_idx])
    ys.append(y[data_idx])
    cs.append(coloring_dic2[c])

  avx=sum(x)/len(x)
  ids=pops[t]
  pop_dic[c]=[avx,ids]
  id_dic[c]=[avx,ind_id[t]]
plt.ylabel("PC 2, "+str(round(pes_list[1], 2))+" %",fontsize=64)
plt.xlabel("PC 1, "+str(round(pes_list[0], 2))+" %",fontsize=64)
plt.xticks([])
plt.yticks([])
print("made first figure")
plt.savefig(args.df1[:-9]+"_PCA.pdf")

with open(args.df1[:-9]+'_pca_data.csv', 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(xs)
    writer.writerow(ys)
    writer.writerow(cs)

print("Hetero freq: ", lens_for_props[0]/(sum(lens_for_props)))
print("Homop freq: ", lens_for_props[1]/(sum(lens_for_props)))
print("Homoq freq: ", lens_for_props[2]/(sum(lens_for_props)))

if NUM_CLUST == 3: #code to find "homozygote minor allele, heterozygote, homozygote major allele"
  sorted_values = sorted(pop_dic.values(), key=lambda x: x[0]) #sort based on average value along PC1
  sorted_values2 = sorted(id_dic.values(), key=lambda x: x[0])

  thomoq=sorted_values[0][1]
  hetero=sorted_values[1][1]
  thomop=sorted_values[2][1]

  thomoq2=sorted_values2[0][1]
  hetero2=sorted_values2[1][1]
  thomop2=sorted_values2[2][1]

  if len(thomop) > len(thomoq): #making sure q is assigned to minor allele
    homop=thomop
    homoq=thomoq
    homop2=thomop2
    homoq2=thomoq2
  else:
    homoq=thomop
    homop=thomoq
    homoq2=thomop2
    homop2=thomoq2

  print("saving genotype ids")
  np.savetxt(args.df1[:-9]+"_homoq_nums.txt",homoq2,fmt="%d")
  np.savetxt(args.df1[:-9]+"_homop_nums.txt",homop2,fmt="%d")

  q_dic=Counter(homoq)
  p_dic=Counter(homop)
  pq_dic={}
  pops=["FOG","CAP","KIB","BOD","TER","LOM","SAN"]
  per_pops=defaultdict(dict)
  for pop in pops:
    one=q_dic[pop]
    two=p_dic[pop]
    if pop == "FOG":
      three=18-(one+two)
    elif pop == "CAP":
      three=19-(one+two)
    else:
      three=20-(one+two)
    pq_dic[pop]=three
    per_pops[pop]["one"]= one
    per_pops[pop]["two"]= two
    per_pops[pop]["three"]= three

else:
  for idx,i in enumerate(id_dic.items()):
    g=i[1][1]
    print("saving genotype ids")
    np.savetxt(args.df1[:-9]+"_clust"+str(idx)+".txt",g,fmt="%d")

  per_pops=defaultdict(dict)
  for i in pop_dic.items():
    temp_dic=Counter(i[1][1])
    for p in ["FOG","CAP","KIB","BOD","TER","LOM","SAN"]:
      per_pops[p][i[0]]=temp_dic[p]

f, axs = plt.subplots(1, 7, sharey=True, figsize=(12, 3), subplot_kw=dict(aspect="equal"))


#Correlation with latitude
pops=["FOG","CAP","KIB","BOD","TER","LOM","SAN"]

temp_nums=[]
for i,k in enumerate(per_pops.keys()):
  labels = []
  sizes = []
  for x, y in per_pops[k].items():
      labels.append(x)
      sizes.append(y)
  temp_nums.append(sizes)
  # Plot
  axs[i].pie(sizes, autopct='%1.0f%%', colors=p_colors)
  axs[i].title.set_text(pops[i])
plt.savefig(args.df1[:-9]+"_pies.pdf")

NS_coordinates = [44.840000, 42.840000, 39.60412, 38.3182, 36.94841, 34.718893, 32.651389]

f, axs = plt.subplots(figsize=(10, 4))
pop_lens = np.array(temp_nums).sum(axis = 1)
x = np.arange(7)
for i,l in enumerate(np.array(temp_nums).T):
  print(l/pop_lens)
  if i == 2:
    het_freqs = l/pop_lens
  if i == 1:
    homop_freqs = l/pop_lens
  if i == 0:
    homoq_freqs = l/pop_lens
  plt.plot(NS_coordinates, l/pop_lens, ".",linestyle="-",markersize=10,c=p_colors[i])
  plt.ylabel("Genotype frequency")
  plt.xticks(NS_coordinates, pops)
  plt.xlabel("Populations")
plt.savefig(args.df1[:-9]+"_lines.pdf")


from os import P_PGID
llcrnrlon=-130 #lower left corner westness
llcrnrlat=30 #lower left corner northernness
urcrnrlon=-110 #upper right corner westness
urcrnrlat=50 #upper right corner northernness

plt.figure(figsize=(8,8))

m = Basemap(
        projection='cyl',
        llcrnrlon=llcrnrlon,
        llcrnrlat=llcrnrlat,
        urcrnrlon=urcrnrlon,
        urcrnrlat=urcrnrlat,
        resolution='l')

m.drawcountries()
m.drawcoastlines()

pop_coords=[(-124,44.8),(-124.5,42.8),(-123.8,39.6),(-123,38.3),(-122,36.9),(-120.6,34.7),(-117.25,32.6)]

for i,k in enumerate(per_pops.keys()):
  labels = []
  sizes = []
  for x, y in per_pops[k].items():
      labels.append(x)
      sizes.append(y)
  # Plot
  plt.pie(sizes, center=(m(pop_coords[i][0],pop_coords[i][1])), colors=p_colors,radius=0.8,textprops={'fontsize': 8},startangle = -90,wedgeprops={"edgecolor":"k"})#,autopct='%1.0f%%')


from matplotlib.patches import Patch

axis = plt.gca()
axis.set_xlim([llcrnrlon, urcrnrlon])
axis.set_ylim([llcrnrlat, urcrnrlat])
plt.savefig(args.df1[:-9]+"_map.pdf")

if NUM_CLUST == 3:

  observed_counts = [len(homoq), len(hetero), len(homop)]
  tot=sum(observed_counts)

  p=(len(homop)*2 + len(hetero))/(tot*2)
  q=(len(homoq)*2 + len(hetero))/(tot*2)
  ehomoq = q*q*tot #24, actual 16
  ehomop = p*p*tot #49, actual 41
  ehete = 2*q*p*tot #68, actual 83

  print(ehomoq/len(homoq)) #expected minor is 1.5 more
  print(ehomop/len(homop)) #expected major is 1.2 more
  print(ehete/len(hetero)) #actual hete is 1.2 more

  result = chisquare(f_obs=observed_counts, f_exp=[ehomoq,ehete,ehomop])
  print(result)
  print(observed_counts)
  print([ehomoq,ehete,ehomop])

  #per pop
  print(pop_lens)
  hets=het_freqs*pop_lens
  homps=homop_freqs*pop_lens
  homqs=homoq_freqs*pop_lens
  for i in range(7):
    pop_p=(homps[i]*2 + hets[i])/(pop_lens[i]*2)
    pop_q=(homqs[i]*2 + hets[i])/(pop_lens[i]*2)
    ehomq = pop_q*pop_q*pop_lens[i]
    ehomp = pop_p*pop_p*pop_lens[i]
    ehet = 2*pop_p*pop_q*pop_lens[i]

    result = chisquare(f_obs=[homqs[i], hets[i],homps[i]], f_exp=[ehomq,ehet,ehomp])
    print(result)


  with open(args.df1[:-9]+"_chi2.txt", "w") as file:
      # Write the variable to the file
      file.write(str(result.pvalue) + "\n")
      file.write(str(ehomoq/len(homoq)) + "\n")
      file.write(str(ehomop/len(homop)) + "\n")
      file.write(str(ehete/len(hetero)) + "\n")


slope, intercept, r_value, p_value, std_err = stats.linregress(NS_coordinates, het_freqs)
print("hetero correlation")
print(slope, r_value, p_value)

slope, intercept, r_value, p_value, std_err = stats.linregress(NS_coordinates, homop_freqs)
print("homop correlation")
print(slope, r_value, p_value)

slope, intercept, r_value, p_value, std_err = stats.linregress(NS_coordinates, homoq_freqs)
print("homoq correlation")
print(slope, r_value, p_value)