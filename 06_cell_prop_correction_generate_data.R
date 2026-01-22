# Author: Alec Eames
# Script: generate_data.R
# Description: Load and process data, prepare cell proportions and CpG regression data for plotting

### Load libraries
library(data.table)
library(dplyr)
library(tidyverse)
library(assertthat)

### Set working directory
setwd(file.path("/Users", "aleceames", "Desktop_no_icloud_sync", 
                "gladyshev_lab", "mahdi_multi_omics_deconv_CONFIDENTIAL", 
                "cell_composition_code_v2_revision"))

### Load methylation data and metadata
mgb_meth <- setDF(fread(file.path("mgb_data", "dnam.csv")))
rownames(mgb_meth) <- mgb_meth$cg
mgb_meth <- mgb_meth[, -1]

mgb_meta <- read.csv(file.path("mgb_data", "sex_age.csv"))
top_cpgs <- read.csv(file.path("mgb_data", "top_cpgs_alec.csv"))
top_cpgs <- distinct(top_cpgs, Name, .keep_all = TRUE)

### Load Salas et al. 2022 deconvolution results
mgb_prop_salas_22 <- read.csv(file.path("mgb_data", "Aging_samples_Salas_12model.csv"), check.names = FALSE)
colnames(mgb_prop_salas_22)[1] <- "new_id"

# Map to old IDs
salas_22_mapping <- read.csv(file.path("mgb_data", "GSE246337_series_matrix_mapping_metadata.csv"))
salas_22_mapping <- salas_22_mapping[, 2:ncol(salas_22_mapping)] %>% t() %>% as.data.frame()
salas_22_mapping <- salas_22_mapping[, c("V1", "V23")]
colnames(salas_22_mapping) <- c("new_id", "old_id")

mgb_prop_salas_22 <- inner_join(mgb_prop_salas_22, salas_22_mapping, by = "new_id")
colnames(mgb_meta)[1] <- "old_id"
mgb_prop_salas_22 <- inner_join(mgb_prop_salas_22, mgb_meta, by = "old_id")
colnames(mgb_prop_salas_22)[16] <- "age"

# Match with methylation data
mgb_prop_salas_22 <- mgb_prop_salas_22[mgb_prop_salas_22$old_id %in% colnames(mgb_meth), ]
idx_match <- match(mgb_prop_salas_22$old_id, colnames(mgb_meth))
assert_that(all(!is.na(idx_match)))
mgb_meth <- mgb_meth[, idx_match]

### Prepare data for Figure 5a
cell_prop_df <- mgb_prop_salas_22[, c(2:13, 16)]
colnames(cell_prop_df) <- c("Basophil", "Memory B cell", "Naive B cell", "Memory CD4+ T cell", 
                            "Naive CD4+ T cell", "Memory CD8+ T cell", "Naive CD8+ T cell", 
                            "Eosinophil", "Monocyte", "Neutrophil", "Natural killer cell", 
                            "Regulatory T cell", "age")
cell_prop_df <- reshape2::melt(cell_prop_df, id = "age")
write.csv(cell_prop_df, file = file.path("plots_data", "fig_5a_data_cell_proportions.csv"), row.names = FALSE)

### Cell composition correction and multivariate regression for each CpG
mgb_prop_salas_22_sub <- mgb_prop_salas_22[, c(2:13, 16)]

# Initial model to check multicollinearity
model_collinear_check <- lm(age ~ Bas + Bmem + Bnv + CD4mem + CD4nv + CD8mem + CD8nv + Eos + Mono + Neu + NK + Treg,
                            data = mgb_prop_salas_22_sub)
ols_vif_tol(model_collinear_check)

alias(model_collinear_check)
# Remove alias (Treg)
model_collinear_check <- lm(age ~ Bas + Bmem + Bnv + CD4mem + CD4nv + CD8mem + CD8nv + Eos + Mono + NK,
                            data = mgb_prop_salas_22_sub)

ols_vif_tol(model_collinear_check)

# Remove cell type with largest VIF until all VIFs under 5

model_collinear_check <- lm(age ~ Bas + Bmem + Bnv + CD4mem + CD4nv + CD8mem + CD8nv + Eos + Mono + NK,
                            data = mgb_prop_salas_22_sub)

ols_vif_tol(model_collinear_check)

# Final model used: remove Treg + Neutrophils to address collinearity

### Run regression for each CpG
age_pval <- numeric()
age_coeff <- numeric()
cpgs <- character()

for (i in 1:nrow(top_cpgs)) {
  print(i)
  cpg_i <- top_cpgs$Name[i]
  cpg_prop_df <- data.frame()
  cpg_prop_df <- cbind(mgb_prop_salas_22_sub, meth = as.numeric(mgb_meth[cpg_i, ]))
  cpg_prop_df <- as.data.frame(scale(cpg_prop_df))
  model <- lm(meth ~ age + Bas + Bmem + Bnv + CD4mem + CD4nv + CD8mem + CD8nv + Eos + Mono + NK, data = cpg_prop_df)
  age_pval[i] <- summary(model)$coefficients["age", 4]
  age_coeff[i] <- summary(model)$coefficients["age", 1]
  cpgs[i] <- cpg_i
}

age_qval <- p.adjust(age_pval, method = "BH")

age_pval_df <- data.frame(
  cpg = cpgs,
  gene = top_cpgs$gene,
  chr = top_cpgs$ch,
  age_coeff = age_coeff,
  age_pval = age_pval,
  age_qval = age_qval
)

age_pval_df$is_sig <- ifelse(age_pval_df$age_qval < 0.05, "Cell-composition-independent", "Cell-composition-dependent")
age_pval_df$age_qval_log10 <- -log10(age_pval_df$age_qval)

highlight_genes <- c("SATB1", "CD248", "CD27", "CD28")
age_pval_df$idx_highlight <- (age_pval_df$age_qval_log10 > 14) |
  (age_pval_df$age_qval < 0.05 & age_pval_df$gene %in% highlight_genes)

write.csv(age_pval_df, file = file.path("plots_data", "fig_5b_cell_composition_correction.csv"), row.names = FALSE)
