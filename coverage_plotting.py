import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

#COORDINATES

#linkage_outlier = 12785121 #largest distance, r2>0.8
linkage_outlier = 16266898

#a = 12780817
#b = 12783216

a=16274350
b=16277490

#PLOTTING

path_list=["homoq_list_cov","hetero_list_cov","homop_list_cov"]
side="right_bp"

all_covs=[]
for p in path_list:
    plt.figure(figsize=(26, 5))
    covs=[]
    file_path='/users/c/p/cpetak/WGS/local_pca_pipe/'+p
    with open(file_path, 'r') as file:
        # Loop through each line in the file
        for line in file:
            #print(line)
            df=pd.read_csv("/users/c/p/cpetak/WGS/inversion_results/"+side+"/"+line.rstrip(),sep='\t',names=["chr",'pos','cov'])
            plt.plot(df["pos"], df['cov'],color="b",alpha=0.1)
            plt.axvline(linkage_outlier, c="red")
            plt.axvline(a, c="blue")
            plt.axvline(b, c="blue")
            plt.xlim(min(df["pos"]), max(df["pos"]))
            plt.xlim(linkage_outlier-100_000, linkage_outlier+100_000)
            plt.ylim(0,70)

            ndf=df[(df["pos"]<=b)&(df["pos"]>=a)]
            covs.append(np.sum(np.array(ndf['cov'].tolist())))
    
    plot_name="/users/c/p/cpetak/WGS/inversion_results/"+p+"_"+side+".png"
    print(plot_name)
    plt.savefig(plot_name)

    all_covs.append(covs)

fig, ax = plt.subplots(figsize=(12, 10))
ax.violinplot(all_covs, showmedians=True)
ax.set_xticks([1,2,3], labels=['Homoq', 'Heterozygote','Homop'])
mylabel="Sum coverage "+str(a)+" - "+str(b)
plt.ylabel(mylabel)
plot_name="/users/c/p/cpetak/WGS/inversion_results/"+"cov_boxplot_"+side+".png"
plt.savefig(plot_name)

h_statistic, p_value = stats.kruskal(all_covs[0], all_covs[1], all_covs[2])

# Print the results
print("homop vs hetero")
print("t-statistic:", h_statistic)
print("p-value:", p_value)

# Save data
#np.save("homop_right2", np.array(covs))
