# 05_candidate_scoring.R
# Trait-informed candidate prioritization framework

library(tidyverse)
library(phyloseq)
library(ggplot2)

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

# Load phyloseq object
ps <- readRDS(
  "data/processed/tomato_microbiome_phyloseq.rds"
)

# Relative abundance transformation
ps_rel <- transform_sample_counts(
  ps,
  function(x) x / sum(x)
)

# Genus-level aggregation
ps_genus <- tax_glom(
  ps_rel,
  taxrank = "Genus",
  NArm = TRUE
)

genus_df <- as.data.frame(
  psmelt(ps_genus)
)

genus_df$Abundance <- as.numeric(
  genus_df$Abundance
)

# --------------------------------
# Literature-informed trait table
# --------------------------------

pgp_traits <- data.frame(
  
  Genus = c(
    "Pseudomonas",
    "Devosia",
    "Massilia",
    "Sphingomonas",
    "Sphingobium",
    "Flavobacterium",
    "Cellvibrio",
    "Pedobacter",
    "Brevundimonas",
    "Rheinheimera"
  ),
  
  Nitrogen_Fixation = c(
    1,1,0,0,0,0,0,0,0,0
  ),
  
  Phosphate_Solubilization = c(
    1,1,1,0,0,1,1,1,0,0
  ),
  
  Siderophore_Production = c(
    1,0,0,1,1,0,0,0,0,0
  ),
  
  Biocontrol = c(
    1,0,0,0,0,0,0,0,0,0
  ),
  
  Stress_Tolerance = c(
    1,0,1,1,1,0,0,0,0,1
  ),
  
  Biofilm_Formation = c(
    1,0,1,1,1,0,0,0,0,0
  ),
  
  Root_Colonization = c(
    1,1,1,1,1,1,1,0,0,0
  )
)

# --------------------------------
# Prevalence + abundance
# --------------------------------

genus_prevalence <- genus_df %>%
  
  group_by(Genus) %>%
  
  summarise(
    
    Prevalence = mean(Abundance > 0),
    
    MeanAbundance = mean(Abundance),
    
    .groups = "drop"
  )

# --------------------------------
# Merge
# --------------------------------

pgp_candidates <- pgp_traits %>%
  
  left_join(
    genus_prevalence,
    by = "Genus"
  )

pgp_candidates[is.na(pgp_candidates)] <- 0

# --------------------------------
# Weighted scoring system
# --------------------------------

pgp_candidates <- pgp_candidates %>%
  
  mutate(
    
    Trait_Score =
      
      Nitrogen_Fixation * 3 +
      
      Phosphate_Solubilization * 2 +
      
      Siderophore_Production * 2 +
      
      Biocontrol * 3 +
      
      Stress_Tolerance * 1 +
      
      Biofilm_Formation * 2 +
      
      Root_Colonization * 3,
    
    
    Abundance_Score = as.numeric(
      scale(log10(MeanAbundance + 1e-6))
    ),
    
    Prevalence_Score = as.numeric(
      scale(Prevalence)
    ),
    
    
    Final_PGP_Score =
      
      Trait_Score +
      
      Abundance_Score +
      
      Prevalence_Score
  )

# --------------------------------
# Validation tiers
# --------------------------------

pgp_candidates <- pgp_candidates %>%
  
  mutate(
    
    Validation_Level =
      case_when(
        
        Final_PGP_Score >= 12 ~
          "Tier 1: High-priority validation",
        
        Final_PGP_Score >= 7 ~
          "Tier 2: Secondary validation",
        
        Final_PGP_Score >= 3 ~
          "Tier 3: Exploratory screening",
        
        TRUE ~
          "Low-priority candidate"
      )
  )

# --------------------------------
# Sort candidates
# --------------------------------

pgp_candidates <- pgp_candidates %>%
  
  arrange(
    desc(Final_PGP_Score)
  )

# Save results
write.csv(
  pgp_candidates,
  "results/pgp_candidate_ranking_enhanced.csv",
  row.names = FALSE
)

# --------------------------------
# Bubble plot
# --------------------------------

candidate_bubble <- ggplot(
  pgp_candidates,
  aes(
    x = Trait_Score,
    y = Abundance_Score,
    size = Final_PGP_Score,
    color = Validation_Level
  )
) +
  
  geom_point(
    alpha = 0.85
  ) +
  
  geom_text(
    aes(label = Genus),
    vjust = -1,
    size = 4
  ) +
  
  theme_classic(base_size = 14) +
  
  labs(
    title = "Trait-Informed Prioritization of PGP Candidates",
    subtitle = "Integrated ecological and functional scoring framework",
    x = "Trait-Based Functional Score",
    y = "Normalized Abundance Score"
  ) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      size = 18
    )
  )

ggsave(
  "results/figures/pgp_candidate_bubble_plot.png",
  candidate_bubble,
  width = 11,
  height = 8,
  dpi = 600
)

cat("Candidate scoring analysis complete\n")