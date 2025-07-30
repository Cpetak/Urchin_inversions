import pandas as pd
import matplotlib.pyplot as plt

#Part 1: make PCA file ---------------------------------------

# Read csv file

df = pd.read_csv('intermediary_files/inv9/NW_022145610.1_30779143_31460853_pca_data.csv', header=None)
df = df.T #transpose
df.columns = ['x', 'y', 'color']
df["x"] = df["x"].astype(float)
df["y"] = df["y"].astype(float)

# Scatter plot
plt.scatter(df['x'], df['y'], c=df['color'])
plt.xlabel('X')
plt.ylabel('Y')
plt.savefig("name_your_file")

#Part 2: make MDS file ------------------------------------------

all_chrs=["NW_022145594.1" ,"NW_022145594.1" ,"NW_022145597.1" ,"NW_022145600.1" ,"NW_022145601.1" ,"NW_022145603.1" ,"NW_022145606.1" ,"NW_022145609.1" ,"NW_022145610.1"]
mtype="snp"

f, axs = plt.subplots(3,3,figsize=(24, 14),sharex=False, sharey=False)
axs=axs.flatten()

counter=0
for c in all_chrs:
    if counter == 8:
        route="~/WGS/Urchin_inversions/lostruct_results/type_" + mtype + "_size_" + str(1000) + "_chromosome_" + c
    else:
        route="~/WGS/Urchin_inversions/lostruct_results/type_" + mtype + "_size_" + str(500) + "_chromosome_" + c
    coords = pd.read_csv(route + "/" + c + ".regions.csv")
    mds = pd.read_csv(route + "/mds_coords.csv")
    lpca = coords.join(mds["MDS1"])
    lpca = lpca.join(mds["MDS2"])
    lpca["pos"] = (lpca["end"]+lpca["start"])/2
    lpca["pos"] = lpca["pos"].astype(int)

    if counter in [0,2,8]:
        axs[counter].plot(lpca["pos"], lpca["MDS2"], ".", color="black")
    else:
        axs[counter].plot(lpca["pos"], lpca["MDS1"], ".", color="black")

    axs[counter].set_title(f'Chromosome: {c}', fontsize=18)
    axs[counter].tick_params(axis='both', which='major', labelsize=16)

    counter+=1

plt.text(-0.2, 0.5, 'Local PCA MDS values', va='center', rotation='vertical', fontsize=22, transform=axs[3].transAxes)
plt.text(0.2, -0.3, 'Genomic position', va='center', rotation='horizontal', fontsize=22, transform=axs[7].transAxes)

plt.subplots_adjust(wspace=0.2)
plt.subplots_adjust(hspace=0.3)

plt.savefig(f'temp_all_mdss.pdf', format='pdf')
