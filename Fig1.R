# R script to produce Figure 1. 
# module load R/4.5.1

#Read packages
library(readr)
library(ggplot2)
library(patchwork)
library(cowplot)

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

  vals <- round(as.numeric(vals), 2)

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
  idx <- match(file_name, c(first_row_files, second_row_files))
  ggplot(df, aes(x, y, color = color)) +
    geom_point(size = 3, alpha = 0.8) + #HERE to adjust size of dots! also adjust in legend! also the transparency
    scale_color_manual(values = c(C0 = "#1f77b4",
                                C1 = "#ff7f0e",
                                C2 = "#2ca02c")) +
    labs(x = paste0("PC 1, ", labs_xy$x, "%"), y = paste0("PC 2, ", labs_xy$y, "%")) +
    annotate(
      "text",
      x = max(df$x, na.rm = TRUE), y = max(df$y, na.rm = TRUE),
      label = idx,
      hjust = 1.1, vjust = 2,
      size = 6
    ) +
    theme_minimal() +
    theme(legend.position = "none",
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"))
}

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

# save
ggsave("all_pca_plots_custom_layout.png", final_plot,
       width = 18, height = 8, dpi = 300, bg = "white")


#Subfigure A

