library(data.table)
library(foreach)
library(stringr)

#### Obtain results

args <- commandArgs(trailingOnly=TRUE)

root_folder <- paste0("~/WGS/local_pca_pipe/makegrid_",args[1],"_",args[2])

files_r <- system(paste("ls -v",root_folder, sep = " "), inter = T)

#Checking missing files
c<-1
for (f in files_r){
  first<-strsplit(f, split = "job")[[1]][2]
  second<-strsplit(first, split = "_")[[1]][1]
  if (strtoi(second) != c){
    print(second)
    c<-c+1
  }
  c<-c+1
}

winsize <- strsplit(files_r[1], split = "_")[[1]][2]
stepsize <- strsplit(files_r[1], split = "_")[[1]][3]
stepsize <- strsplit(stepsize, split = "[.]")[[1]][1]

# Now load all results using a foreach loop with rbind as .combine
o <- foreach(i=files_r, .combine = "rbind")%do%{
  
  message(i)
  tmp <- get(load(paste(root_folder,i, sep = "/")))
  return(tmp)
  
}
setDT(o)
save(o, file=paste("~/WGS/local_pca_pipe/makegrid_",args[1],"_",args[2], "/combined_", args[1], "_", args[2], ".Rdata", sep = ""))
