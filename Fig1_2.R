library(ggplot2)

all_chrs <- c(
  "NW_022145594.1","NW_022145594.1","NW_022145597.1",
  "NW_022145600.1","NW_022145601.1","NW_022145603.1",
  "NW_022145606.1","NW_022145609.1","NW_022145610.1"
)

mtype <- "snp"

all_data <- list()

for (i in seq_along(all_chrs)) {
  c <- all_chrs[i]

  if (i == 9) {
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
  print(route)

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

  if (i == 2){
    df$chromosome <- paste0(c,"2")
  } else {
    df$chromosome <- c
  }
  df$panel <- i

  all_data[[i]] <- df
}

plot_df <- do.call(rbind, all_data)

plot_df$sign <- ifelse(plot_df$value < 0, "negative", "positive")

p <- ggplot(plot_df, aes(x = pos, y = abs(value), color = sign)) +
  geom_point(size = 0.6) +
  facet_wrap(~ chromosome, nrow = 1, scales = "free") +
  labs(
    x = "Genomic position",
    y = "Local PCA MDS values"
  ) +
  scale_color_manual(values = c(negative = "red", positive = "black")) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 14),
    axis.text  = element_text(size = 12),
    axis.title = element_text(size = 16)
  )

ggsave("temp_all_mdss_fromR.pdf", p, width = 24, height =7)

