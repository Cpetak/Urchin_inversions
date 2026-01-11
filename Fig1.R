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

#Plot (individual panels so legend can occupy last grid cell)
make_panel_plot <- function(panel_num) {
  dfp <- subset(plot_df, panel == panel_num)
  ggplot(dfp, aes(x = pos, y = abs(value), color = sign)) +
    geom_point(size = 2, alpha = 0.8) +
    scale_x_continuous(labels = function(x) paste0(x / 1e6, "M")) +
    scale_y_continuous(limits = c(0, 0.8)) +
    scale_color_manual(values = c(Negative = "red", Positive = "black"), name = NULL) +
    theme_minimal() +
    theme(
      strip.text = element_blank(),
      axis.text  = element_text(size = 12),
      axis.title = element_blank(),
      legend.position = "none",
      panel.grid = element_blank(),
      axis.line = element_line(color = "black")
    ) +
    annotate( #labels 1-9
      "text",
      x = Inf, y = Inf, label = panel_num,
      hjust = 2.5, vjust = 1.2,
      size = 6
    )
}

# Create panels 1-9
panels <- lapply(1:9, make_panel_plot)

# Arrange into two rows: first row 1-5, second row 6-9 + placeholder for legend
panels_row1 <- panels[1:5]
panels_row2 <- panels[6:9]
panels_row2[[5]] <- ggplot() + theme_void()

# Extract legend from a plot (showing legend)
legend_plot <- ggplot(plot_df, aes(x = pos, y = abs(value), color = sign)) +
  geom_point(size =4) +
  scale_color_manual(values = c(Negative = "red", Positive = "black"), name = NULL) +
  theme_bw() +
  theme(legend.position = "right", legend.text = element_text(size = 18))
legend <- get_legend(legend_plot)

# Place legend in the empty slot
panels_row2[[5]] <- legend

final_plot_A <- plot_grid(
  plotlist = c(panels_row1, panels_row2),
  ncol = 5,
  align = 'hv'
)

# Draw the grid inset to create larger left/bottom margins, then add shared labels in those margins
canvas <- ggdraw() +
  draw_plot(final_plot_A, x = 0.06, y = 0.06, width = 0.88, height = 0.88)

final_plot_A_labeled <- canvas +
  draw_label("Genomic position", x = 0.5, y = 0.02, vjust = 0, hjust = 0.5, size = 18) +
  draw_label("Local PCA |MDS values|", x = 0.05, y = 0.5, angle = 90, vjust = 1, hjust = 0.5, size = 18)

ggsave("Fig1_A.pdf", final_plot_A_labeled, width = 34, height = 12)

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

