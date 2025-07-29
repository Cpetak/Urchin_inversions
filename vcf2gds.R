library(SeqArray)
library(SNPRelate)

seqParallelSetup(cluster=10, verbose=TRUE)

args <- commandArgs(trailingOnly = TRUE)
vcf_file <- args[1]
out_file <- args[2]

print(args)
#print(paste(args[1],"sometime"))

snpgdsVCF2GDS(vcf_file,paste(out_file,".gds",sep=""),verbose=T)
