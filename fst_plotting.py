import pandas as pd
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser(description="Description of your script")
parser.add_argument("chr", type=str, help="Chromosome name")
parser.add_argument("start", type=int, help="start region")
parser.add_argument("stop", type=int, default=0.1, help="stop region")
args = parser.parse_args()

filename="~/WGS/inversion_results/" + args.chr + "_" + str(args.start) + "_" + str(args.stop) + "_fst.weir.fst"
df=pd.read_csv(filename, sep="\t")

fig, ax1 = plt.subplots(figsize=(16, 8))
plt.plot(df["POS"], df["WEIR_AND_COCKERHAM_FST"], ".")
plt.ylabel("Fst, homo1 vs homo2")
plt.xlabel("Genomic position")

plt.savefig(args.chr + "_" + str(args.start) + "_" + str(args.stop) + "fst_plot.jpg")
