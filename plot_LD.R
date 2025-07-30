library(ggplot2)
library(data.table)
library(viridis)

### load data
#load("~/WGS/LD_chr/makegrid_outs/job1_50000_10000.Rdata")

args <- commandArgs(trailingOnly=TRUE)

if (interactive()) {
  # Define defaults for interactive mode
  args <- c("NW_022145597.1", "500","all")
}

#file <- paste("~/WGS/Urchin_inversions/ld_data", "/combined_", args[1], "_", args[2],args[3],".Rdata", sep = "")
file <- paste("~/WGS/Urchin_inversions/ld_data", "/combined_", args[1], "_", args[2],".Rdata", sep = "")
print(file)
#outfile <- paste(args[1], "_", args[2], "_LD_", args[3], ".pdf", sep="")
outfile <- paste(args[1], "_", args[2], "_LD.pdf", sep="")

load(file)

### pad empty spaces
grid <- data.table(expand.grid(1:max(c(o$win1, o$win2)), 1:max(c(o$win1, o$win2))))
setnames(grid, names(grid), c("win1", "win2"))

grid <- grid[win1<win2]
setkey(o, win1, win2)
setkey(grid, win1, win2)

#o2 <- merge(o[poolOnly==T][rnp.thr==1], grid, all=T)
#o2[is.na(meanR2), meanR2:=-.01]

table(o$win1>o$win2)
#table(o2$win1>o2$win2)

#mstart<-14219351
#mstop<-14298524
#line1 <- o[o$start1 <= mstart & o$stop1 >= mstart, ]$win1[1]
#line2 <- o[o$start1 <= mstop & o$stop1 >= mstop, ]$win1[1]
#line3 <- o[o$start2 <= mstop & o$stop2 >= mstop, ]$win2[1]

new_o <- o#[o$stop1 < upper & o$stop2 < upper, ]
#new_o <- new_o[new_o$stop1 > lower & new_o$stop2 > lower]

#testing
#abs(win1-win2)>0 makes it so that the diagonal is not shown. 
#this changes the r2 ranges shown because the diagonal will have high values

p1 <- ggplot(data=new_o[abs(win1-win2)>0][meanR2>0], aes(x=win1, y=win2, fill=meanR2)) +
  geom_tile() + 
  theme_minimal() + coord_fixed(ratio = 1)+
  #geom_tile(data=new_o[!complete.cases(new_o), ], aes(x=win1, y=win2), fill="black", alpha=0.4) +
  #scale_fill_gradient2(low = "blue", mid = "white", high = "red")+
  scale_fill_viridis(option="H") #limits = c(0, 0.2), trans = "sqrt")
  #theme(
    #panel.background = element_rect(fill = "lightblue")  # change panel background
  #)
  #geom_hline(yintercept = line3, color = "red", linetype = "dashed", size = 1) +
  #geom_vline(xintercept = line2, color = "red", linetype = "dashed", size = 1)

p1

ggsave(outfile, plot = p1, width = 7, height = 5)

