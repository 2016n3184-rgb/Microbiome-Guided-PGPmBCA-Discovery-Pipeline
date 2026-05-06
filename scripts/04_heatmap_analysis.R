# 04_heatmap_analysis.R
# Genus-level abundance heatmap

library(tidyverse)
library(phyloseq)
library(pheatmap)

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

ps <- readRDS("data/processed/tomato_microbiome_phyloseq.rds")

ps_rel <- transform_sample_counts(
  ps,
  function(x) x / sum(x)
)

ps_genus <- tax_glom(
  ps_rel,
  taxrank = "Genus",
  NArm = TRUE
)

genus_df <- as.data.frame(psmelt(ps_genus))
genus_df$Abundance <- as.numeric(genus_df$Abundance)

top_genera <- genus_df %>%
  group_by(Genus) %>%
  summarise(
    TotalAbundance = sum(Abundance),
    .groups = "drop"
  ) %>%
  arrange(desc(TotalAbundance)) %>%
  head(20)

write.csv(
  top_genera,
  "results/top20_abundant_genera.csv",
  row.names = FALSE
)

heatmap_df <- genus_df %>%
  filter(Genus %in% top_genera$Genus)

heatmap_summary <- heatmap_df %>%
  group_by(Genus, Zone) %>%
  summarise(
    MeanAbundance = mean(Abundance),
    .groups = "drop"
  )

heatmap_matrix <- heatmap_summary %>%
  pivot_wider(
    names_from = Zone,
    values_from = MeanAbundance
  ) %>%
  column_to_rownames("Genus") %>%
  as.matrix()

heatmap_matrix[is.na(heatmap_matrix)] <- 0

heatmap_matrix_log <- log10(
  heatmap_matrix + 1e-6
)

saveRDS(
  heatmap_matrix_log,
  "data/processed/top_genera_heatmap_matrix.rds"
)

png(
  "results/figures/top_genera_heatmap.png",
  width = 1800,
  height = 1400,
  res = 220
)

pheatmap(
  heatmap_matrix_log,
  scale = "row",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  fontsize_row = 10,
  fontsize_col = 12,
  border_color = NA,
  main = "Top Genera Across Tomato Root Compartments"
)

dev.off()

cat("Heatmap analysis complete\n")