library("SNPRelate")

args <- commandArgs(trailingOnly = TRUE)
filename <- args[1]
genofile <- snpgdsOpen(filename)

ccm_pca<-snpgdsPCA(genofile, autosome.only=FALSE)
dim1<-ccm_pca$eigenvect[,1]
dim2<-ccm_pca$eigenvect[,2]

eigenval_list <- ccm_pca$eigenval
eigenval_list[is.nan(eigenval_list)] <- 0
esum <- sum(eigenval_list)
e1 <- eigenval_list[1] / esum * 100
e2 <- eigenval_list[2] / esum * 100

my_df1 <- as.data.frame(dim1)
my_df2 <- as.data.frame(dim2)

special_character <- "\\."
split_string <- strsplit(filename, special_character)[[1]]
result1 <- paste(paste(head(split_string, -1), collapse = "."),"_dim1.csv",sep="")
result2 <- paste(paste(head(split_string, -1), collapse= "."),"_dim2.csv",sep="")

eresult <- paste(paste(head(split_string, -1), collapse= "."),"_perc_explained.csv",sep="")

writeLines(c(as.character(e1), as.character(e2)), eresult)

write.csv(my_df1, result1, row.names=FALSE)
write.csv(my_df2, result2, row.names=FALSE)

