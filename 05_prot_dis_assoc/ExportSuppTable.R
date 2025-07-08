#!/usr/bin/env Rscript

# Load necessary libraries
library(data.table)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggsci)
library(magrittr)
library(forcats)
library(patchwork)
library(ggrepel)
library(writexl)

rm(list=ls())

## set wd
path <- "/sc-projects/sc-proj-computational-medicine/people/Kamil/projects/03_UKB_Protein_Variance/01_Cox_prediction/output"

##########################################
#### Supplementary Table Prot-Disease ####
##########################################
## load non-adj
result_files_non_adj <- list.files(
  path,
  pattern = "^res\\.surv\\.[0-9].*\\.txt$",
  full.names = TRUE
)

## rbind
combined_data_non_adj <- rbindlist(
  lapply(result_files_non_adj, fread),
  fill = TRUE
)

## excl 5 phecodes with cases < 200
## 438 phecodes
combined_data_non_adj %<>% filter(!is.na(beta))
## excl non-significant
## rename combined non-adj data
final_data <- combined_data_non_adj %>% filter(pval<0.05/(2919*438))

## define the patterns for adjusted data
prefixes <- c(
  "res.surv.adj.",
  "res.surv.adj.top3.",
  "res.surv.adj.top5.",
  "res.surv.adj.top10.",
  "res.surv.adj.top15.",
  "res.surv.cis.",
  "res.surv.trans."
)

## suffixes to label columns accordingly
suffixes <- c(
  "adj", "top3", "top5", "top10", "top15", "cis", "trans"
)

## loop
for (i in seq_along(prefixes)) {
  prefix <- prefixes[i]
  suffix <- suffixes[i]
  
  ## find matching
  files <- list.files(
    path,
    pattern = paste0("^", prefix, "[0-9].*\\.txt$"),
    full.names = TRUE
  )
  
  if (length(files) == 0) next
  
  ## load and rbind for each adj set
  tmp_data <- rbindlist(
    lapply(files, fread),
    fill = TRUE
  )
  
  ## select columns to keep, for cis & trans also keep n all & event bc these differ
  cols_to_keep <- c("id", "phecode", "beta", "se", "pval")
  if (suffix %in% c("cis", "trans")) {
    cols_to_keep <- c(cols_to_keep, "nall", "nevent")
  }
  
  tmp_data <- tmp_data[, ..cols_to_keep]
  
  ## rename columns accordingly
  new_names <- setNames(
    paste0(names(tmp_data)[-(1:2)], "_", suffix),
    names(tmp_data)[-(1:2)])
  
  setnames(tmp_data, old=names(new_names), new=new_names)
  
  ## left join with significant non-adj
  final_data <- left_join(
    final_data,
    tmp_data,
    by = c("id", "phecode")
  )
}

rm(tmp_data, combined_data_non_adj)

## add the phecode labels
lab.phe <- fread("/sc-projects/sc-proj-computational-medicine/data/UK_biobank/phecodes/data/Comparison.phecodes.cases.UKBB.UCL.Charite.20210811.txt", sep="\t", header=T)
lab.phe <- lab.phe[,c(1:3)]
lab.phe$phecode <- as.factor(lab.phe$phecode)
final_data$phecode <- as.factor(final_data$phecode)
final_data %<>% left_join(., lab.phe, by="phecode")

## reorder columns
final_data <- final_data %>%
  select(phenotype, category, phecode, everything())

## write to csv
fwrite(final_data, "output/SupplementaryTableProtDisAssoc.csv", row.names = F, na = NA)

##########################################
#### AAA MMP12 non-adj vs full adjust ####
##########################################

## load non-adjusted
file_non_adj <- file.path(path, "res.surv.442.11.txt")
dat_non_adj <- fread(file_non_adj) %>%
  select(id, phecode, beta_nonadj = beta, pval_nonadj = pval)

## load adjusted
file_adj <- file.path(path, "res.surv.adj.442.11.txt")
dat_adj <- fread(file_adj) %>%
  select(id, phecode, beta_adj = beta, pval_adj = pval)

## merge
volcano_AAA_data <- left_join(dat_non_adj, dat_adj, by = c("id", "phecode"))

volcano_AAA_minimal_top10 <- volcano_AAA_data %>% 
  filter(pval_nonadj < 0.05/(2919*438)) %>% 
  arrange(desc(beta_nonadj)) %>% slice(1:10)

a <- volcano_AAA_data %>%
  ggplot(aes(x=beta_nonadj, y=-log10(pval_nonadj)))+
  geom_hline(yintercept = -log10(0.05/(2919*438)), linetype = "dashed", color = "gray30")+
  geom_point(fill="#4DBBD5FF",
             color = "black", shape = 21,alpha=.3) +
  geom_label_repel(
    data = volcano_AAA_minimal_top10, aes(label = id),
    size = 2, color = "grey20",
    segment.color = "grey60",
    max.overlaps = Inf,
    nudge_y = .5,
    fill = alpha(c("white"),0.5)
  ) +
  ylim(0,100) + xlim(-0.5, 1.5)+
  labs(x = "Beta coefficient for abdominal aortic aneurysm",
       y = "-log10(P value)",
       title = "Age & sex adjustment") +
  theme_light() +
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        panel.grid = element_blank(),
        strip.text = element_text(color="white"),
        strip.background = element_rect(
          color="black", fill="black", linetype="solid"))

volcano_AAA_extended_top10 <- volcano_AAA_data %>% 
  filter(pval_adj < 0.05/(2919*438)) %>% 
  arrange(desc(beta_adj)) %>% slice(1:10)

b <- volcano_AAA_data %>% 
  ggplot(aes(x=beta_adj, y=-log10(pval_adj)))+
  geom_hline(yintercept = -log10(0.05/(2919*438)), linetype = "dashed", color = "gray30")+
  geom_point(fill="#4DBBD5FF",
             color = "black", shape = 21,alpha=.3) +
  geom_label_repel(
    data = volcano_AAA_extended_top10, aes(label = id),
    size = 2, color = "grey20",
    segment.color = "grey60",
    max.overlaps = Inf,
    nudge_y = .5,
    fill = alpha(c("white"),0.5)
  ) +
  ylim(0,100) + xlim(-0.5, 1.5)+
  labs(x = "Beta coefficient for abdominal aortic aneurysm",
       y = "-log10(P value)",
       title = "Comprehensive adjustment") +
  theme_light() +
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        panel.grid = element_blank(),
        strip.text = element_text(color="white"),
        strip.background = element_rect(
          color="black", fill="black", linetype="solid"))

jpeg("graphics/Supplementary_Figure_AAA.jpg",units = "in",width=7, height=3.5, res=1500)
a+b + plot_layout(guides="collect") & theme(legend.position="bottom", text=element_text(size=8))
dev.off()

pdf("graphics/Supplementary_Figure_AAA.pdf", width=7, height=3.5)
a+b + plot_layout(guides="collect") & theme(legend.position="bottom", text=element_text(size=8))
dev.off()
