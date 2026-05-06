# 01_import_and_qc.R
# Import, clean, filter and prepare tomato 16S microbiome data

library(tidyverse)
library(phyloseq)

metadata <- read_tsv("data/raw/16S_rar_mctools_meta.txt")

asv_table <- read.delim(
  "data/raw/16S_rar_mctools_data.txt",
  sep = "\t",
  header = TRUE,
  skip = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

taxonomy <- asv_table %>%
  select(`#OTU ID`, taxonomy)

counts <- asv_table %>%
  select(-taxonomy)

colnames(counts)[1] <- "ASV"
colnames(taxonomy)[1] <- "ASV"

counts_matrix <- counts %>%
  column_to_rownames("ASV")

asv_sums <- rowSums(counts_matrix)

filtered_counts <- counts_matrix[
  asv_sums >= 10,
]

colnames(filtered_counts) <- paste0("TR_", colnames(filtered_counts))

metadata$sample_id <- as.character(metadata$sample_id)

metadata <- metadata %>%
  filter(sample_id %in% colnames(filtered_counts))

metadata <- metadata[
  match(colnames(filtered_counts), metadata$sample_id),
]

metadata <- as.data.frame(metadata)
rownames(metadata) <- metadata$sample_id

taxonomy_clean <- taxonomy %>%
  filter(ASV %in% rownames(filtered_counts)) %>%
  separate(
    taxonomy,
    into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
    sep = "; ",
    fill = "right"
  )

tax_matrix <- taxonomy_clean %>%
  column_to_rownames("ASV") %>%
  as.matrix()

OTU <- otu_table(
  as.matrix(filtered_counts),
  taxa_are_rows = TRUE
)

META <- sample_data(metadata)
TAX <- tax_table(tax_matrix)

ps <- phyloseq(OTU, META, TAX)

saveRDS(ps, "data/processed/tomato_microbiome_phyloseq.rds")

write.csv(filtered_counts, "data/processed/filtered_asv_table.csv")
write.csv(metadata, "data/processed/metadata_clean.csv", row.names = FALSE)
write.csv(taxonomy_clean, "data/processed/taxonomy_table.csv", row.names = FALSE)

cat("Import and QC complete\n")
cat("Filtered taxa:", ntaxa(ps), "\n")
cat("Samples:", nsamples(ps), "\n")