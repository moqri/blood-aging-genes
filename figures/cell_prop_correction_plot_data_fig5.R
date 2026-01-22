### Author: Alec Eames
### Script: plot_data.R
### Description: Generate plots for cell proportions vs age and cell-composition-corrected CpG effects

# Load libraries
library(tidyverse)
library(ggrepel)
library(ggbreak)
library(scales)

# Set working directory
setwd(file.path("/Users", "aleceames", "Desktop_no_icloud_sync", 
                "gladyshev_lab", "mahdi_multi_omics_deconv_CONFIDENTIAL", 
                "cell_composition_code_v2_revision"))

# Load preprocessed data
cell_prop_df <- read.csv(file.path("plots_data", "fig_5a_data_cell_proportions.csv"))
age_pval_df <- read.csv(file.path("plots_data", "fig_5b_cell_composition_correction.csv"))

# Plot: Cell proportion vs Age

palette <- c('#8c564b', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#1f77b4', 
             '#e377c2', '#7f7f7f', '#bcbd22', '#17becf', 'grey30', 'orange3')

cell_prop_age <- ggplot(cell_prop_df, aes(x = age, y = value, fill = variable, color = variable)) +
  geom_point(alpha = 0.03, stroke = NA) +
  geom_smooth(method = "lm") +
  theme_classic() +
  scale_color_manual(values = palette, name = "Cell type") +
  scale_fill_manual(values = palette, name = "Cell type") +
  labs(x = "Age", y = "Cell proportion", title = "Cell proportion vs age") +
  scale_x_continuous(minor_breaks = scales::breaks_width(4)) +
  scale_y_continuous(minor_breaks = scales::breaks_width(0.01)) +
  guides(
    x = guide_axis(minor.ticks = TRUE),
    y = guide_axis(minor.ticks = TRUE)
  ) +
  theme(
    axis.ticks.length = unit(5, "pt"),
    axis.minor.ticks.length = rel(0.6),
    axis.line = element_line(linewidth = 0.3, color = "grey60"),
    axis.minor.ticks.x.bottom = element_line(linewidth = 0.3, color = "grey60"),
    axis.minor.ticks.y.left = element_line(linewidth = 0.3, color = "grey60"),
    axis.ticks = element_line(linewidth = 0.45, color = "grey20"),
    axis.line.y.right = element_blank(),
    axis.ticks.y.right = element_blank(),
    axis.text.y.right = element_blank()) +
  scale_y_break(c(0.16, 0.5), scales = 0.8, space = 0.55, ticklabels = NULL,
                expand = c(0, 0))
cell_prop_age


ggsave(file.path("plots", "fig_5a_salas_12_mgb_cell_prop_vs_age.png"), plot = cell_prop_age, 
       dpi = 400, width = 9, height = 7, scale = 0.6)

# Plot: Age effect on methylation corrected for cell composition
age_pval_df$is_sig <- factor(age_pval_df$is_sig, 
                             levels = c("Cell-composition-independent", "Cell-composition-dependent"))

cell_prop_corrected_point_plot <- ggplot(age_pval_df, aes(x = age_coeff, y = age_qval_log10, color = is_sig, label = gene)) +
  geom_point(stroke = NA, size = 2.2, alpha = 0.8) +
  theme_classic() +
  scale_color_manual(name = "", values = c('#1f77b4', '#ff7f0e')) +
  labs(x = "Age-CpG coefficient", y = expression(-log[10]~("adjusted p-value")), title = "Cell-composition-corrected aging CpGs") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey60", linewidth = 0.3) +
  geom_text_repel(aes(label = ifelse(idx_highlight, as.character(gene), '')), 
                  size = 2,
                  show.legend = F,
                  max.overlaps = 30,
                  color = "black") +
  scale_x_continuous(minor_breaks = scales::breaks_width(0.01)) +
  scale_y_continuous(minor_breaks = scales::breaks_width(1)) +
  guides(
    x = guide_axis(minor.ticks = TRUE),
    y = guide_axis(minor.ticks = TRUE)
  ) +
  theme(
    legend.spacing.x = unit(0, 'cm'),
    legend.position = "top",
    legend.justification = "left",
    legend.text = element_text(size = 7, margin = ggplot2::margin(l = -2)),
    legend.margin = ggplot2::margin(t = 0, l = 0, unit='cm'),
    axis.ticks.length = unit(5, "pt"),
    axis.minor.ticks.length = rel(0.6),
    axis.ticks = element_line(linewidth = 0.45, color = "grey20"),
    axis.line = element_line(linewidth = 0.3, color = "grey60"),
    axis.minor.ticks.x.bottom = element_line(linewidth = 0.3, color = "grey60"),
    axis.minor.ticks.y.left = element_line(linewidth = 0.3, color = "grey60"))
cell_prop_corrected_point_plot


ggsave(file.path("plots", "fig_5b_cell_prop_corrected_independent_point_plot_lrm_clock_foundation_additional_cell_types.png"),
       plot = cell_prop_corrected_point_plot, dpi = 400, width = 3.7, height = 3.5, scale = 1.2)
