# R script to produce Figure 1. 

#Read packages
library(readr)
library(ggplot2)

#Read csv

df <- read_csv(
  "intermediary_files/inv9/NW_022145610.1_30779143_31460853_pca_data.csv",
  col_names = FALSE
)

df <- as.data.frame(t(df))
colnames(df) <- c("x", "y", "color")

df$x <- as.numeric(df$x)
df$y <- as.numeric(df$y)

#Make PCA

p <- ggplot(df, aes(x = x, y = y, color = color)) +
  geom_point() +
  xlab("X") +
  ylab("Y")

ggsave("test_R_PCA.png", plot = p)