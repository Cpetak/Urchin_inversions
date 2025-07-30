library(topGO)

args <- commandArgs(trailingOnly=TRUE)
chr <- args[1]
mystart <- args[2]
mystop <- args[3]

geneID2GO <- readMappings("~/WGS/Urchin_inversions/supp_files/GO_mapping_topGO") # uniprot to GO mapping
geneNames <- names(geneID2GO)

get_go_res <- function(ont_type,geneList,fname_out){
  GOdata <- new("topGOdata", 
                ontology = ont_type, # ontology of interest (BP, MF or CC)
                allGenes = geneList,
                annot = annFUN.gene2GO, 
                gene2GO = geneID2GO)
  resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
  print(resultFisher)
  num_sig<-length(resultFisher@score[resultFisher@score < 0.01])
  print("num significant")
  print(num_sig)
  allRes <- GenTable(GOdata, classicFisher = resultFisher, topNodes = num_sig)
  print("printing allRes")
  print(allRes)
  if (num_sig > 0) {
  write.table(allRes, file = fname_out, sep = ",", 
              append = TRUE, quote = FALSE, 
              col.names = FALSE, row.names = FALSE) 
  }
}

for (i in 1:3){
  fname_in <- paste0(chr,"_",mystart,"_",mystop,"_","uniprot",i,".txt")
  print(fname_in)
  fname_out <- paste0(chr,"_",mystart,"_",mystop,"_","GoEn",i,".csv")
  
  myInterestingGenes <- read.csv(fname_in, header = FALSE) # list of interesting genes, output of LOC to uniprot mapping
  intgenes <- myInterestingGenes[, "V1"]
  geneList <- factor(as.integer(geneNames %in% intgenes)) # mask of 0 and 1 if geneName is interesting
  names(geneList) <- geneNames # geneList but annotated with the gene names
  
  for (t in list("BP","MF","CC")){
    get_go_res(ont_type=t,geneList=geneList,fname_out=fname_out)
  }
  
}




