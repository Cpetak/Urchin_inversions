library(ggplot2)
library(data.table)
library(RColorBrewer)
library(dplyr)
library(ggplot2)
library(fuzzyjoin)
library(broom)
library(readr)
library(rstatix)
##localscore code https://forge-dga.jouy.inra.fr/projects/local-score/
source("/projappl/project_2003522/localscore/scorelocalfunctions.R")
##snp list with positions and chromosome for whole genome
snps<-read.table("snp_pos_chr.txt", header=T)
##read baypass results
xtx<-read.table("baypass_summary_pi_xtx.out", header=T)
##combine files
xtx_pos_dt<-cbind(snps,xtx)
xtx_pos_dt <- as.data.table(xtx_pos_dt)


# calculate lindley score peaks -------------------------------------------

setkey(xtx_pos_dt, chr)
xtx_pos_dt$pval[xtx_pos_dt$pval==0]=1e-16#min(mydata$pval[-which(mydata$pval ==0)])
Nchr=length(xtx_pos_dt[,unique(chr)])

### Computation of absolute position in the genome. 
#This is useful for doing genomewide plots. 
chrInfo=xtx_pos_dt[,.(L=.N,cor=autocor(log10.1.pval.)),chr]
setkey(chrInfo,chr)
data.table(chr=xtx_pos_dt[,unique(chr),], S=cumsum(c(0,chrInfo$L[-Nchr])))

### Choice of $\xi$ (1,2,3 or 4)

xi=2
xtx_pos_dt[,score:= log10.1.pval.-xi]
xtx_pos_dt[,lindley:=lindley(score),chr]

mean(xtx_pos_dt$score)

# verify mean is negative and max is positive
mean(xtx_pos_dt$log10.1.pval.);max(xtx_pos_dt$log10.1.pval.)

#hist(mydata$score)
max(xtx_pos_dt$lindley)

# Compute significance threshold for each chromosome
## Uniform distribution of p-values
chrInfo[,th:=thresUnif(L, cor, xi),chr]
(xtx_pos_dt=xtx_pos_dt[chrInfo])

(sigZones=xtx_pos_dt[chrInfo,sig_sl(lindley, pos, unique(th)),chr])



# plot results ------------------------------------------------------------


invplot <- xtx_pos_dt %>%
  mutate(
    Inversion_flag = case_when(
      chr == "NW_022145594.1" & pos >= 12702886 & pos <= 16793794 ~ "1",
      chr == "NW_022145594.1" & pos >= 39429440 & pos <= 42445994 ~ "2",
      chr == "NW_022145597.1" & pos >= 14219351 & pos <= 14298524 ~ "3",
      chr == "NW_022145600.1" & pos >= 33842803 & pos <= 34582517 ~ "4",
      chr == "NW_022145601.1" & pos >= 28950395 & pos <= 29566247 ~ "5",
      chr == "NW_022145603.1" & pos >= 8291260 & pos <= 9094120     ~ "6",
      chr == "NW_022145606.1" & pos >= 16610900 & pos <= 16737625   ~ "7",
      chr == "NW_022145609.1" & pos >= 29427741 & pos <= 29776196   ~ "8",
      chr == "NW_022145610.1" & pos >= 30779143 & pos <= 31460853   ~ "9"
    ),
    Sig_Inversion_flag = case_when(
      log10.1.pval. >= 2 & !is.na(Inversion_flag) ~ "sig",
      TRUE ~ "non-sig"
    )
  )

###using sig snps that are above 0.01 (-log10 2)
invplot_sig<-invplot %>% filter(log10.1.pval. >= 2)


##Get chromosomes with at least one inversion
chr_with_inversions <- test_sig %>%
  filter(!is.na(Inversion_flag)) %>%
  distinct(chr) %>%
  pull(chr)

##Filter to those chromosomes, and create a grouping variable
plot_data <- test_sig %>%
  filter(chr %in% chr_with_inversions) %>%
  mutate(
    Inversion_group = ifelse(is.na(Inversion_flag), "Outside inversion", paste0("Inv", Inversion_flag))
  )


invplot <- invplot %>%
  mutate(label_x = (xmin + xmax) / 2, label_y = max(xtx_pos_dt$lindley, na.rm = TRUE) * 1.05)

p1 <- ggplot(xtx_pos_dt, aes(x = pos_cum, y = lindley, colour = as.factor(chr))) +
  # Highlight inversion regions
  geom_rect(data = invplot, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE, fill = "red4", alpha = 0.8) +
    geom_point(size = 1, alpha = 0.7) +
  geom_hline(yintercept = min(xtx_pos_dt$th), linetype = "dashed", color = "red", linewidth = 0.8) +
  geom_hline(yintercept = max(xtx_pos_dt$th), linetype = "solid", color = "red", linewidth = 0.8) +
  scale_x_continuous(label = chr_labels$chr, breaks = chr_labels$center) +
  scale_color_manual(values = rep(c("green4", "purple3"), length(unique(xtx_pos_dt$chr)))) +
  theme_classic(base_size = 14) +
  theme(
    panel.grid.major.y = element_line(color = "grey90", linetype = "dotted"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 12),
    legend.position = "none"
  ) +
  labs(x = "Chromosome", y = "Lindley Score")+ 
  geom_text(data = invplot, aes(x = label_x, y = label_y, label = Inversion),
            inherit.aes = FALSE, size = 4, vjust = 0, fontface = "bold",
            hjust= -0.5)
output_dir <- "/gpfs2/scratch/dsadler1/SeaUrchinEnvPop/GDS/BayPass/Plots"  # Change to your desired path

ggsave(
    filename = file.path(output_dir, "pvalue_xtx_null_invhighlight_v2.png"),
    plot = p1,
    dpi = 300,
    width = 20,
    height = 8,
    units = "in"
  )



# Fishers Exact test ------------------------------------------------------



setwd("/gpfs2/scratch/dsadler1/SeaUrchinEnvPop/GDS/sigoutliersstrict/FEanalyses")


load("annotationfileallsnps_combined.RData")

## Clean inversion peak table
invpeaks <- read.csv("inversionID.csv") %>%
  mutate(across(c(Chr), ~ gsub("[^0-9A-Za-z_.]", "", as.character(.))),
         startpos = as.numeric(gsub("[^0-9]", "", as.character(startpos))),
         stopos  = as.numeric(gsub("[^0-9]", "", as.character(stopos))))

## Prepare SNP annotations (get highest-impact annotation per SNP)
annot_highest_impact <- combined_annots %>%
  group_by(SNP_id) %>%
  arrange(Annotation_Impact) %>%
  slice(1) %>%
  ungroup()

allsnps <- annot_highest_impact[, "SNP_id"]

##Extract chr and pos from SNP_id
snp_positions <- annot_highest_impact %>%
  mutate(
    chr = sub("^(.*)_(\\d+)$", "\\1", SNP_id),
    pos = as.numeric(sub("^.*_(\\d+)$", "\\1", SNP_id))
  )

##Find SNPs that fall within any inversion region using fuzzy join
snp_in_inversion <- fuzzy_inner_join(
  snp_positions,
  invpeaks,
  by = c("chr" = "Chr", "pos" = "startpos", "pos" = "stopos"),
  match_fun = list(`==`, `>=`, `<=`)
)
save(snp_in_inversion, file="snp_in_inversion.RData")

inv_snps <- unique(snp_in_inversion$SNP_id)

## Load XTX-significant SNPs
xtx <- read.csv("/gpfs2/scratch/dsadler1/SeaUrchinEnvPop/GDS/BayPass/random_chunks/nullmodel/sigsnpspeakhighestimpact.csv")

xtx_prepped <- xtx %>%
  mutate(SNP_id = paste(chr, pos, sep = "_"))

sig_xtx <- xtx_prepped$SNP_id

## Fisher test comparing XTX and inversion overlaps
fisher_test_for_inversion <- function(var_name = "Inversion") {
  all_snp_ids <- allsnps$SNP_id
  
  A <- sum(all_snp_ids %in% intersect(inv_snps, sig_xtx))     # In both inversion and XTX
  B <- sum(all_snp_ids %in% setdiff(inv_snps, sig_xtx))       # In inversion only
  C <- sum(all_snp_ids %in% setdiff(sig_xtx, inv_snps))       # In XTX only
  D <- sum(!(all_snp_ids %in% union(inv_snps, sig_xtx)))      # In neither
  
  mat <- matrix(c(A, B, C, D), nrow = 2, byrow = TRUE,
                dimnames = list("Inversion" = c("Sig", "NotSig"),
                                "XTX" = c("Sig", "NotSig")))
  
  test_result <- fisher.test(mat)
  tidy_result <- broom::tidy(test_result)
  
  tidy_result$variable <- var_name
  tidy_result$A <- A
  tidy_result$B <- B
  tidy_result$C <- C
  tidy_result$D <- D
  tidy_result
}

fisher_test_for_inversion()


# XtX comparison ----------------------------------------------------------

# Get chromosomes with at least one inversion
chr_with_inversions <- test_sig %>%
  filter(!is.na(Inversion_flag)) %>%
  distinct(chr) %>%
  pull(chr)

# Filter to those chromosomes, and create a grouping variable
plot_data <- test_sig %>%
  filter(chr %in% chr_with_inversions) %>%
  mutate(
    Inversion_group = ifelse(is.na(Inversion_flag), "Outside inversion", paste0("Inv", Inversion_flag))
  )


ggplot(plot_data, aes(x = chr, y = log10(XtXst), fill = Inversion_group)) +
  geom_boxplot(outliers= F, alpha = 0.8, position = position_dodge(width = 0.75)) +
  scale_fill_brewer(palette = "Set1") +
  scale_colour_brewer(palette = "Set1") +
  theme_minimal(base_size = 13) +
  labs(
    title = "XtXst Distribution by Chromosome and Inversion Region",
    x = "Chromosome",
    y = "log10(XtXst)",
    fill = "Region"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


results <- plot_data %>%
  filter(chr %in% inversion_chrs) %>%
  mutate(inv_status = ifelse(Inversion_group == "Outside inversion", "outside", "inside")) %>%
  group_by(chr) %>%
  summarise(
    # Perform the Wilcoxon test manually so we can extract more details
    test = list(wilcox.test(XtXst ~ inv_status)),
    p_value = test[[1]]$p.value,
    w_stat = test[[1]]$statistic,
    
    # Sample sizes
    n_inside = sum(inv_status == "inside"),
    n_outside = sum(inv_status == "outside"),
    
    # Medians
    median_inside = median(XtXst[inv_status == "inside"]),
    median_outside = median(XtXst[inv_status == "outside"]),
    
    # Rank-biserial correlation (effect size)
    effect_size_r = (2 * test[[1]]$statistic) / (n_inside * n_outside) - 1
  ) %>%
  arrange(p_value)

