# 03_differential_abundance.R
# Differential abundance analysis using DESeq2

library(tidyverse)
library(phyloseq)
library(DESeq2)
library(ggplot2)

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

# Load processed phyloseq object
ps <- readRDS("data/processed/tomato_microbiome_phyloseq.rds")

# Make Zone a clean factor
sample_data(ps)$Zone <- as.factor(sample_data(ps)$Zone)

# Check available Zone levels
print(levels(sample_data(ps)$Zone))

# Build DESeq2 object
dds <- phyloseq_to_deseq2(
  ps,
  ~ Zone
)

# Run DESeq2 with poscounts for sparse microbiome data
dds <- DESeq(
  dds,
  fitType = "parametric",
  sfType = "poscounts"
)

# Print available result names
print(resultsNames(dds))

# Use exact coefficient name from DESeq2
res_rhizo <- results(
  dds,
  name = "Zone_Rhizosphere_vs_Bulk.soil"
)

# Convert to normal dataframe
res_rhizo_df <- as.data.frame(res_rhizo)
res_rhizo_df$ASV <- rownames(res_rhizo_df)

res_rhizo_df$padj <- as.numeric(res_rhizo_df$padj)
res_rhizo_df$log2FoldChange <- as.numeric(res_rhizo_df$log2FoldChange)

# Add taxonomy
tax_table_df <- as.data.frame(tax_table(ps))
tax_table_df$ASV <- rownames(tax_table_df)

res_rhizo_df <- left_join(
  res_rhizo_df,
  tax_table_df,
  by = "ASV"
)

# Save full DESeq2 results
write.csv(
  res_rhizo_df,
  "results/rhizosphere_vs_bulksoil_deseq2.csv",
  row.names = FALSE
)

# Top 20 taxa by adjusted p-value
top_taxa <- res_rhizo_df[
  !is.na(res_rhizo_df$padj),
]

top_taxa <- top_taxa[
  order(top_taxa$padj),
]

top_taxa <- head(top_taxa, 20)

write.csv(
  top_taxa,
  "results/top20_rhizosphere_vs_bulksoil_taxa.csv",
  row.names = FALSE
)

# Remove missing/unclassified genera for cleaner plot
top_taxa_clean <- top_taxa %>%
  filter(
    !is.na(Genus),
    Genus != "",
    Genus != "NA"
  )

# Differential abundance plot
diff_plot <- ggplot(
  top_taxa_clean,
  aes(
    x = reorder(Genus, log2FoldChange),
    y = log2FoldChange,
    fill = log2FoldChange
  )
) +
  geom_col(width = 0.72) +
  scale_fill_gradient2(
    low = "#2563EB",
    mid = "#F8FAFC",
    high = "#DC2626",
    midpoint = 0
  ) +
  coord_flip() +
  theme_classic(base_size = 14) +
  labs(
    title = "Differentially Abundant Genera",
    subtitle = "Rhizosphere vs Bulk Soil",
    x = NULL,
    y = "Log2 Fold Change"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 18),
    plot.subtitle = element_text(size = 11),
    axis.text.y = element_text(face = "italic", size = 11),
    axis.text.x = element_text(color = "black"),
    legend.position = "none"
  )

ggsave(
  "results/figures/differential_genera_bright.png",
  diff_plot,
  width = 10,
  height = 6,
  dpi = 600
)

cat("Differential abundance analysis complete\n")