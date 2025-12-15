# R script to produce Figure 1. 

#Read packages
library(readr)
library(ggplot2)

#Subfigure B

#Find all files
dir_path <- "intermediary_files/inv9"

files <- list.files(
  path = dir_path,
  pattern = "_pca_data\\.csv$",
  full.names = TRUE
)

#Read csvs
# read all files into one data frame
all_df <- lapply(files, function(f) {
  df <- read_csv(f, col_names = FALSE)
  df <- as.data.frame(t(df))
  colnames(df) <- c("x", "y", "color")
  df$x <- as.numeric(df$x)
  df$y <- as.numeric(df$y)
  df$file <- basename(f)
  df
})

all_df <- do.call(rbind, all_df)

#Make PCAs

p <- ggplot(all_df, aes(x = x, y = y, color = color)) +
  geom_point() +
  facet_wrap(~ file, nrow = 1) +
  xlab("X") +
  ylab("Y")

ggsave("all_pca_plots.png", p, width = 4 * length(files), height = 4)