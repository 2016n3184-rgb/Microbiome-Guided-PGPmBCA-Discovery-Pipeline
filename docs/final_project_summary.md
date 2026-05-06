
# Final Project Summary

## Project Title
Microbiome-Guided PGP-mBCA Discovery Pipeline

## Aim
This project develops a reproducible R-based workflow for analyzing public tomato root-associated microbiome data and prioritizing microbial genera with potential relevance as plant growth-promoting microbial biocontrol agent candidates.

## Dataset
The project uses real public 16S rRNA amplicon data from a tomato belowground microbiome study. The analysis focuses on bacterial community patterns across bulk soil, rhizosphere, rhizoplane and endosphere compartments.

## Workflow
1. Data import and QC
2. Phyloseq object construction
3. Alpha diversity analysis
4. Hellinger-transformed PCA
5. DESeq2 differential abundance
6. Genus-level heatmap
7. Trait-informed candidate scoring
8. Validation and limitations framework

## Main Finding
Tomato-associated belowground microbiomes show compartment-associated ecological structuring. Candidate genera including Pseudomonas, Devosia, Massilia, Sphingomonas and Sphingobium were prioritized using trait evidence, abundance and prevalence.

## Candidate Prioritization Logic
The final score integrates:
- weighted PGP trait score
- mean genus-level abundance
- prevalence across samples
- validation confidence tier

## Interpretation
This project does not claim confirmed PGP-mBCA activity. It provides a reproducible prioritization framework for selecting microbial candidates for downstream culturing, WGS, antiSMASH analysis and greenhouse validation.

## Key Limitations
- 16S data cannot resolve strain-level function.
- Literature-informed traits require experimental validation.
- Relative abundance data are compositional.
- Candidate taxa require culturing and greenhouse testing.

## Future Work
Future work should include isolate recovery, functional assays, WGS, AMR/virulence screening, antiSMASH biosynthetic gene cluster mining, greenhouse validation and synthetic community testing.

