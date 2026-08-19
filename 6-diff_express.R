
# Script function: Differential Gene Expression analysis using DESeq2
# Version: 1.0
# Date created: 18-AUG-2026
# 
# Dependencies:
#  - DESeq2 (Differential Gene Expression analysis)
#    
#  - ggplot2 (Plotting figures)
#    
#  = pheatmap (Plotting heatmaps)
#    
#  - RColorBrewer (Colour palettes)
#    
#  - ggrepel (Non-overlapping plot labels)
#    
#  - PoiClaClu (Poisson distance calculation)
#    
#  - ashr (LFC shrinkage for composite contrasts)
#    

#===============================================================================
# SETUP
#===============================================================================

# Load packages
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(ggrepel)
library(PoiClaClu)
library(ashr)

#===============================================================================
# DESEQ2
#===============================================================================

# Wald Test --------------------------------------------------------------------

# The Wald Test finds the significance of specific contrasts

# Create DESeq object with interaction design: ~ sex + timepoint + sex:timepoint
# This models the main effect of sex, the main effect of timepoint, and whether
# the time course differs between sexes (the interaction)
dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = sample_info,
    design = ~ sex + timepoint + sex:timepoint
)

# Remove genes with low counts and update DESeq object
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

# Set 2h post-emergence as the reference level
# All timepoint contrasts will be relative to this baseline
dds$timepoint <- relevel(dds$timepoint, ref = "02")

# Run Wald test (default DESeq2 statistical test)
dds <- DESeq(dds)

# Set 2h post-emergence as the reference level
# All timepoint contrasts will be relative to this baseline
dds$timepoint <- relevel(dds$timepoint, ref = "02")

# Run Wald test (default DESeq2 statistical test)
dds <- DESeq(dds)

# Print available named coefficients for reference
# TODO: Remove if unecessary
cat("\nAvailable contrasts:\n")
print(resultsNames(dds))

# Likelihood Ratio Test (LRT) --------------------------------------------------

# The LRT tests whether the interaction term as a whole (across all timepoints)
# is significant, providing more statistical power than testing the interaction
# at each timepoint individually.
#
# This type of statistical test is particularly useful for timecourse
# experimental designs.

# Create new DESeq object with interaction design as before
ddsTC <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = sample_info,
    design = ~ sex + timepoint + sex:timepoint
)

# Filter out low reads and set baseline
keep_tc <- rowSums(counts(ddsTC)) >= 10
ddsTC <- ddsTC[keep_tc, ]
ddsTC$timepoint <- relevel(ddsTC$timepoint, ref = "02")

# Compare full model vs reduced model without interaction term
# This determines the significance of the interaction
ddsTC <- DESeq(ddsTC, test = "LRT", reduced = ~ sex + timepoint)

#===============================================================================
# RESULTS
#===============================================================================

# Wald Test --------------------------------------------------------------------
