# R script to produce Figure 1. 
# module load R/4.5.1

#Read packages
library(readr)
library(ggplot2)
library(patchwork)
library(cowplot)

#Subfigure A

#Read data
all_chrs <- c(
  "NW_022145594.1","NW_022145594.1","NW_022145597.1",
  "NW_022145600.1","NW_022145601.1","NW_022145603.1",
  "NW_022145606.1","NW_022145609.1","NW_022145610.1"
)
mtype <- "snp"
all_data <- list()

for (i in seq_along(all_chrs)) {
  c <- all_chrs[i]
  
  if (i == 9) { #For the last chromosome, 1000 snp window size was better
    route <- paste0(
      "~/WGS/Urchin_inversions/lostruct_results/type_",
      mtype, "_size_1000_chromosome_", c
    )
  } else {
    route <- paste0(
      "~/WGS/Urchin_inversions/lostruct_results/type_",
      mtype, "_size_500_chromosome_", c
    )
  }
  
  coords <- read.csv(file.path(route, paste0(c, ".regions.csv")))
  mds    <- read.csv(file.path(route, "mds_coords.csv"))
  
  df <- cbind(coords, MDS1 = mds$MDS1, MDS2 = mds$MDS2)
  df$pos <- as.integer((df$start + df$end) / 2)
  
  if ((i - 1) %in% c(0, 2, 8)) {
    print(i)
    df$value <- df$MDS2
  } else {
    df$value <- df$MDS1
  }
  
  #to avoid collapsing data due to inversions 1 and 2 being on the same chromosome
  if (i == 2){
    df$chromosome <- paste0(c,"2") 
  } else {
    df$chromosome <- c
  }
  df$panel <- i
  all_data[[i]] <- df
}

plot_df <- do.call(rbind, all_data)
plot_df <- plot_df[!is.na(plot_df$value), ] #drop NA
plot_df$sign <- ifelse(plot_df$value < 0, "Negative", "Positive") #remember which is negative

#Plot
p <- ggplot(plot_df, aes(x = pos, y = abs(value), color = sign)) +
  geom_point(size = 2, alpha = 0.8) + #HERE you can change the size of the dots
  facet_wrap(~ panel, nrow = 1, scales = "free") +
  scale_y_continuous(limits = c(0, 0.8)) + #HERE you can change the y axis limits
  labs(
    x = "Genomic position",
    y = "Local PCA |MDS values|"
  ) +
  scale_color_manual(values = c(Negative = "red", Positive = "black"), name = NULL) +
  theme_bw() +
  theme(
    strip.text = element_blank(),
    axis.text  = element_text(size = 12),
    axis.title = element_text(size = 16)
  ) +
  geom_text( #HERE is where I label them 1-9
    data = data.frame(panel = unique(plot_df$panel), idx = unique(plot_df$panel)),
    aes(x = Inf, y = Inf, label = idx),
    inherit.aes = FALSE,
    hjust = 2.5,   # change hjust and vjust to adjust location of labels
    vjust = 1.2,
    size = 6
  )

#Moving around legend so that it's not outside the whole plot
# 1.Remove legend from main plot
p_no_legend <- p + theme(legend.position = "none")

# 2.Extract legend from a plot
legend <- get_legend(p + theme(legend.position = "right", legend.text = element_text(size = 14))) #HERE you can change size of legend labels

# 3.Combine plot and legend
p <- ggdraw(p_no_legend) +
  draw_plot(legend, x = 0.89, y = 0.6, width = 0.15, height = 0.3)  #HERE adjust position of legend

ggsave("Fig1_A.pdf", p, width = 24, height =4)

#-----------------------------------------------------------------------------------------------

#Subfigure B

#Find all files
dir_path <- "intermediary_files"

files <- list.files(
  path = dir_path,
  pattern = "_pca_data\\.csv$",
  full.names = TRUE,
  recursive = TRUE
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

#Get PCA percent explained labels
get_axis_labels <- function(file_name) {
  label_name <- sub("_pca_data\\.csv$", "_perc_explained.csv", file_name)

  label_path <- list.files(
    path = dir_path,
    pattern = paste0("^", label_name, "$"),
    full.names = TRUE,
    recursive = TRUE
  )

  vals <- read_csv(label_path[1],
                   col_names = FALSE,
                   show_col_types = FALSE)[[1]]

  vals <- round(as.numeric(vals), 2) #HERE you can adjust level of rounding

  list(
    x = as.character(vals[1]),
    y = as.character(vals[2])
  )
}

#Two rows, where last "subplot" is legend

#Function to create plots
make_plot <- function(file_name) {
  df <- subset(all_df, file == file_name)
  labs_xy <- get_axis_labels(file_name)
  idx <- match(file_name, c(first_row_files, second_row_files)) #label 1-9
  
  ggplot(df, aes(x, y, color = color)) +
    geom_point(size = 3, alpha = 0.8) + #HERE to adjust size of dots! also adjust in legend! also the transparency
    scale_color_manual(values = c(C0 = "#1f77b4", #COLORs where specifically chosen to match the rest of the plots! (if we are still using them)
                                C1 = "#ff7f0e",
                                C2 = "#2ca02c")) +
    labs(x = paste0("PC 1, ", labs_xy$x, "%"), y = paste0("PC 2, ", labs_xy$y, "%")) +
    annotate( #labels 1-9
      "text",
      x = max(df$x, na.rm = TRUE), y = max(df$y, na.rm = TRUE),
      label = idx,
      hjust = 1.1, vjust = 2, #HERE you can adjust location of labels
      size = 6
    ) +
    theme_minimal() +
    theme(legend.position = "none",
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"))
}

#BIT OF MAGIC to make the custom layout with the last "plot" being the legend
# Make plots
n_files <- length(files)
first_row_files <- basename(files)[1:5]
second_row_files <- basename(files)[6:9]

plots_row1 <- lapply(first_row_files, make_plot)
plots_row2 <- lapply(second_row_files, make_plot)

# Add empty placeholder for legend in last cell of second row
plots_row2[[5]] <- ggplot() + theme_void()

# Extract legend from one plot
legend_plot <- ggplot(subset(all_df, file == first_row_files[1]),
                      aes(x, y, color = color)) +
  geom_point(size = 3) + #HERE to adjust size of the dots!
  scale_color_manual(
    values = c(C0 = "#1f77b4", C1 = "#ff7f0e", C2 = "#2ca02c"),
    labels = c(C0 = "Less frequent homokaryotypes", C1 = "More frequent homokaryotypes", C2 = "Heterokaryotypes"),
    name = NULL
  ) +
  theme_minimal() +
  theme(legend.position = "right",
  legend.text = element_text(size = 14))
legend <- get_legend(legend_plot)


# Add legend to second row
plots_row2[[5]] <- legend

# Combine everything
final_plot <- plot_grid(
  plotlist = c(plots_row1, plots_row2),
  ncol = 5,
  align = 'hv'
)

ggsave("Fig1_B.png", final_plot,
       width = 18, height = 8, dpi = 300, bg = "white")

