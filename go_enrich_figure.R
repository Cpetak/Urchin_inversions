# Load required libraries
library(ggplot2)
library(dplyr)
library(readr)
library(stringr)

# Read topGO output CSV without column names and assign names explicitly
topgo_data <- read_csv("go_enrich_results/NW_022145594.1_39429440_42445994_GoEn1.csv", col_names = FALSE) # Replace with your file path

#Specifying GO terms to keep due to p-value adjustments and geneontology.com hierarchical structure
go_to_keep <- c("GO:1900037", "GO:0032205", "GO:2001243", "GO:0008544", "GO:0009913", "GO:0030855", "GO:0043124", "GO:2001242", "GO:1902532","GO:0010521","GO:0140313","GO:0140311","GO:0008395", "GO:1990879", "GO:0034663", "GO:0000782", "GO:0140445", "GO:0030126", "GO:0030663", "GO:0032993")

# Assign column names manually
colnames(topgo_data) <- c("GO", "GO_term", "Annotated", "Significant", "Expected", "classicFisher")

print(topgo_data)

# Preprocess data: calculate GeneRatio, assign sizes, and filter significant terms
topgo_filtered <- topgo_data %>%
  mutate(GeneRatio = as.numeric(Significant) / as.numeric(Annotated),       # Calculate GeneRatio
         Count = as.numeric(Annotated),                       # Use Significant as Count
         p.adjust = as.numeric(classicFisher)) %>%  # Ensure classicFisher is numeric
  #filter(p.adjust < 0.01) %>%                       # Filter significant GO terms
  #filter(Count < 500 & Count > 10) %>%
  arrange(p.adjust)                                 # Sort by p-value

topgo_filtered <- topgo_filtered[topgo_filtered$GO %in% go_to_keep, ]

print(nrow(topgo_filtered))

#order by GO in geneontology.com, BF, MF CC.
gos_kept <- rev(go_to_keep[go_to_keep %in% topgo_filtered$GO])

topgo_filtered <- topgo_filtered %>%
  mutate(GO = factor(GO, levels = topgo_filtered$GO[match(gos_kept, topgo_filtered$GO)])) %>%
  arrange(factor(GO, levels = gos_kept), desc(GeneRatio)) # Arrange by go_to_keep order, then by GeneRatio

#topgo_filtered$GO_term <- str_wrap(topgo_filtered$GO_term, width = 20)

# Make the dot plot
p <- ggplot(topgo_filtered, aes(x = GeneRatio,
                           y = GO, #reorder(GO_term, GeneRatio), # Reorder terms by GeneRatio
                           size = Count,
                           color = p.adjust)) +
  geom_point() +
  scale_y_discrete(labels = topgo_filtered$GO_term) +
  scale_color_gradient(low = "red", high = "blue",
                       name = "p (adj.)") +
  scale_size(range = c(1, 8), name = "Count") +
  labs(x = "Gene Ratio",
       y = NULL) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20),   # Title font size
    axis.title.x = element_text(size = 16), # X-axis label font size
    axis.title.y = element_text(size = 16), # Y-axis label font size
    axis.text = element_text(size = 14),     # Axis tick label font size
    legend.title = element_text(size = 16),     # Legend title font size
    legend.text = element_text(size = 14),
    panel.grid.major = element_line(color = "gray70"), # Darker major grid lines
    panel.grid.minor = element_line(color = "gray70")
  )

#NOTE: some labels are cut off on the left side of the figure. This is because the csv files already have them cut off!
#Meaning that in order to make them look nice, you make to manually replace the entries in the csv file with the longer version of the label.
#E.g. hypo... to hypoxia!
#I've done this for inversion 2.

ggsave(
  filename = "inv2_GO_plot.pdf", # File name
  plot = p)

