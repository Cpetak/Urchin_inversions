import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
import numpy as np
import argparse

parser = argparse.ArgumentParser(description="Description of your script")
parser.add_argument("chrom", type=str, help="Chromosome name")
parser.add_argument("mtype", type=str, help="snp or bp")
parser.add_argument("--thr", type=float, default=0.1, help="threshold for MDS outliers, default=0.1")

args = parser.parse_args()
chrom=args.chrom
mtype=args.mtype
thr=args.thr

def get_groups(values, k, df):
    # Convert the list to a numpy array
    X = np.array(values).reshape(-1, 1)

    # Define the number of clusters (2 in this case, for the two groups)
    n_clusters = k

    # Perform K-means clustering
    kmeans = KMeans(n_clusters=n_clusters,n_init="auto",random_state=42)
    kmeans.fit(X)

    # Get the cluster labels
    cluster_labels = kmeans.labels_

    # Initialize lists for each cluster
    clusters = [[] for _ in range(n_clusters)]

    # Assign each value to its corresponding cluster
    for value, label in zip(values, cluster_labels):
        clusters[label].append(value)

    # Print the clusters
    outs=[]
    for i, cluster in enumerate(clusters):
        #print(f"Cluster {i + 1}: {cluster}")
        out_start=df[df["pos"]==min(cluster)]["start"].tolist()[0]
        out_end=df[df["pos"]==max(cluster)]["end"].tolist()[0]
        outs.append([out_start, out_end])
        
    return sorted(outs, key=lambda x: x[0])

def get_figs(size, m, ax, mythr):
    
    route="~/WGS/inversion_results/lostruct_results/type_" + mtype + "_size_" + str(size) + "_chromosome_" + chrom
    coords = pd.read_csv(route + "/" + chrom + ".regions.csv")
    print(coords.head())
    mds = pd.read_csv(route + "/mds_coords.csv")
    lpca = coords.join(mds["MDS1"])
    lpca = lpca.join(mds["MDS2"])
    lpca["pos"] = (lpca["end"]+lpca["start"])/2
    lpca["pos"] = lpca["pos"].astype(int)
    
    if abs(lpca[m].max()) > abs(lpca[m].min()):
        outs = lpca[lpca[m]>=mythr]
    else:
        mythr = -1*mythr
        outs = lpca[lpca[m]<=mythr]

    ax.plot(lpca["pos"], lpca[m], ".") #plotting window middle points
    ax.plot(outs["pos"], outs[m], ".")
    ax.axhline(mythr, color="red")
    
    all_regions=[]
    for k in range(2):
        #int_regions is actual genomic coordinates, not middle points!
        int_regions=get_groups(outs["pos"].tolist(), k+1, lpca)
        
        if k == 0:
            linestyle="--"
            ax.axvline(int_regions[0][0], linestyle=linestyle, color="black")
            ax.axvline(int_regions[0][1], linestyle=linestyle, color="black")
        else:
            linestyle=":"
            ax.axvline(int_regions[0][0], linestyle=linestyle, color="black")
            ax.axvline(int_regions[0][1], linestyle=linestyle, color="black")
            ax.axvline(int_regions[1][0], linestyle=linestyle, color="black")
            ax.axvline(int_regions[1][1], linestyle=linestyle, color="black")
               
        all_regions.append(int_regions)
            
    ax.text(1, 1, all_regions,
        verticalalignment='bottom', horizontalalignment='right',
        transform=ax.transAxes,
        color='black', fontsize=8) 
    
    ax.set_title(m+" Window size: "+str(size),fontsize=8)
    
            
        
f, axs = plt.subplots(8,1,figsize=(14, 24),sharex=True)

lens=[500,1000,5000,10000]
c=0
for m in ["MDS2","MDS1"]:
    for l in lens:
        get_figs(l,m,axs[c],thr)
        c+=1
plt.subplots_adjust(hspace=0.3)
plt.xlabel("Genomic position")
fname="outlier_regions_"
plt.savefig(fname+ mtype + "_chromosome_" + chrom+".pdf")
