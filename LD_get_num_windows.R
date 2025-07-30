args <- commandArgs(trailingOnly=TRUE)
myFileName <- args[1]
myWindowSize <- as.integer(args[2])

pattern <- "NW_(.*?)\\.1"
myChrom <- regmatches(myFileName, regexpr(pattern, myFileName, perl=TRUE))

### libraries
library(SNPRelate)
library(SeqArray)
library(data.table)
library(foreach)
library(doMC)
registerDoMC(20)

### convert to snprelate GDS
#vcf.fn <- "/project/berglandlab/Dmel_Single_Individuals/Phased_Whatshap_Shapeit4_final/CM_pops.AllChrs.Whatshap.shapeit.annot.wSNPids.vcf.gz"
gds.fn <- myFileName

### open
genofile <- snpgdsOpen(gds.fn, allow.fork=TRUE)

### make snp table
snp.dt <- as.data.table(snpgdsSNPList(genofile, sample.id=NULL))
snp.dt <- snp.dt[chromosome==myChrom][afreq>.05 & afreq<.95]
setnames(snp.dt, c("chromosome", "position"), c("chr", "pos"))
setkey(snp.dt, chr, pos)

### bring in environmental model
#load("/project/berglandlab/alan/environmental_ombibus_global/temp.max;2;5.Cville/temp.max;2;5.Cville.glmRNP.Rdata")
#glm.out <- glm.out[perm==0]
#setkey(glm.out, chr, pos)

### merge
#snp.dt <- merge(snp.dt, glm.out, all.x=T)
#table(is.na(snp.dt$perm))

### make windows
win.size <- myWindowSize
step.size <- myWindowSize

wins <- data.table(chr=myChrom,
                   start=seq(from=min(snp.dt$pos), to=max(snp.dt$pos)-win.size, by=step.size),
                   stop=seq(from=min(snp.dt$pos), to=max(snp.dt$pos)-win.size, by=step.size) + win.size)
print(dim(wins))
