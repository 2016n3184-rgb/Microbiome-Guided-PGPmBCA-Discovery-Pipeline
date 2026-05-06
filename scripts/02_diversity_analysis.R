# 02_diversity_analysis.R
# Alpha diversity and PCA analysis

library(tidyverse)
library(phyloseq)
library(vegan)
library(ggplot2)

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

# Load processed phyloseq object
ps <- readRDS(
  "data/processed/tomato_microbiome_phyloseq.rds"
)

# -----------------------------
# Alpha diversity
# -----------------------------

alpha_div <- estimate_richness(
  ps,
  measures = c("Shannon")
)

alpha_div$sample_id <- rownames(alpha_div)

metadata <- as.data.frame(sample_data(ps))

metadata$sample_id <- rownames(metadata)

alpha_df <- left_join(
  alpha_div,
  metadata,
  by = "sample_id"
)

write.csv(
  alpha_df,
  "results/alpha_diversity_results.csv",
  row.names = FALSE
)

shannon_plot <- ggplot(
  alpha_df,
  aes(
    x = Zone,
    y = Shannon,
    fill = Zone
  )
) +
  
  geom_boxplot(
    alpha = 0.8,
    outlier.shape = NA
  ) +
  
  geom_jitter(
    width = 0.15,
    size = 1.8,
    alpha = 0.7
  ) +
  
  theme_classic(base_size = 14) +
  
  labs(
    title = "Shannon Diversity Across Tomato Compartments",
    x = "Compartment",
    y = "Shannon Diversity"
  ) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      size = 18
    ),
    
    legend.position = "none"
  )

ggsave(
  "results/figures/shannon_diversity.png",
  shannon_plot,
  width = 8,
  height = 6,
  dpi = 600
)

# -----------------------------
# PCA analysis
# -----------------------------

otu <- as(otu_table(ps), "matrix")

if(!taxa_are_rows(ps)) {
  otu <- t(otu)
}

otu_hellinger <- decostand(
  t(otu),
  method = "hellinger"
)

pca <- prcomp(
  otu_hellinger,
  scale. = FALSE
)

pca_df <- as.data.frame(pca$x)

pca_df$sample_id <- rownames(pca_df)

pca_df <- left_join(
  pca_df,
  metadata,
  by = "sample_id"
)

write.csv(
  pca_df,
  "results/pca_coordinates.csv",
  row.names = FALSE
)

pca_plot <- ggplot(
  pca_df,
  aes(
    x = PC1,
    y = PC2,
    color = Zone
  )
) +
  
  geom_point(
    size = 3,
    alpha = 0.85
  ) +
  
  theme_classic(base_size = 14) +
  
  labs(
    title = "PCA of Tomato Root-Associated Microbiomes",
    subtitle = "Hellinger-transformed ASV abundances",
    x = "PC1",
    y = "PC2"
  ) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      size = 18
    )
  )

ggsave(
  "results/figures/pca_hellinger_final.png",
  pca_plot,
  width = 10,
  height = 7,
  dpi = 600
)

cat("Diversity analysis complete\n")