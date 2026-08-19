
# Script function: Performing initial DESeq analysis, then removing potentially
#   potentially mislabelled samples after QC
# Version: 1.0
# Date created: 18-AUG-2026
# 
# Dependencies:
#  - DESeq2 (Differential Gene Expression analysis)
#    
#  - ggplot2 (Plotting figures)
#    
#  - RColorBrewer (Colour palettes)
#    
#  - ggrepel (Non-overlapping plot labels)
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

#===============================================================================
# DESEQ2
#===============================================================================

# Create DESeq object with interaction design: ~ sex + timepoint + sex:timepoint
# This models the main effect of sex, the main effect of timepoint, and whether
# the time course differs between sexes (the interaction)

dds_og <- DESeqDataSetFromMatrix(
    countData = counts_og,
    colData = sample_info_og,
    design = ~ sex + timepoint +sex:timepoint
)

# Remove genes with low counts and update DESeq object
keep <- rowSums(counts(dds_og)) >= 10
dds_og <- dds_og[keep, ]

# Set 2h post-emergence as the reference level
# All timepoint contrasts will be relative to this baseline
dds_og$timepoint <- relevel(dds_og$timepoint, ref = "02")

# Run Wald test (default DESeq2 statistical test)
dds_og <- DESeq(dds_og)

#===============================================================================
# VARIANCE STABILISING TRANSFORMATION (VST)
#===============================================================================

# VST normalises count data for visualisation (PCA)
# blind = FALSE uses the design formula to inform the transformation
vsd_og <- vst(dds_og, blind = FALSE)

#===============================================================================
# PRINCIPAL COMPONENT ANALYSIS (PCA)
#===============================================================================

# Extract PCA coordinates from the VST data
pca_data_og <- plotPCA(
  vsd_og,                                         # Input dataset
  intgroup = c("sex", "timepoint", "replicate"),  # Grouping variables
  returnData = TRUE                               # PCA data.frame for plotting
)

# Calculate percentage variance
percent_var_og <- round(100 * attr(pca_data_og, "percentVar"))

# Define colourblind-friendly palettes
colour_blind_timepoint <- c(
    "02" = "#D55E00",
    "12" = "#CC79A7",
    "24" = "#0072B2",
    "48" = "#F0E442",
    "96" = "#009E73")

# PCA plot
ggplot(pca_data_og, aes(x = PC1, y = PC2)) +

    # Females (filled)
    geom_point(
        data = subset(pca_data_og, sex == "female"),
        aes(colour = timepoint, fill = timepoint, shape = replicate),
        size = 2.5,
        stroke = 1
    ) +
    
    # Males (open)
    geom_point(
        data = subset(pca_data_og, sex == "male"),
        aes(colour = timepoint, shape = replicate),
        fill = "white",
        size = 2.5,
        stroke = 1
    ) +
    
    # Add labels to each sample point
    geom_text_repel(aes(label = name, colour = timepoint), size = 2.5) +
    
    scale_colour_manual(
        # Add previously defined colourblind palette
        values = colour_blind_timepoint,
        # Change timepoint labels in figure legend
        labels = c(
            "02" = "2h",
            "12" = "12h",
            "24" = "24h",
            "48" = "48h",
            "96" = "96h"
        )
    ) +
    
    scale_fill_manual(
        # Add previously defined colourblind palette
        values = colour_blind_timepoint,
        # Change timepoint labels in figure legend
        labels = c(
            "02" = "2h",
            "12" = "12h",
            "24" = "24h",
            "48" = "48h",
            "96" = "96h"
        )
    ) +
    
    # Sample shape is defined by replicate
    scale_shape_manual(values = c(
        "1" = 21,  # Circle
        "2" = 22,  # Square
        "3" = 24   # Triangle
    )) +
    
    # Add percentage vairance to axes labels
    xlab(paste0("PC1: ", percent_var_og[1], "% variance")) +
    ylab(paste0("PC2: ", percent_var_og[2], "% variance")) +

    # Figure legend customisation
    labs(
        colour = "Timepoint",
        fill = "Timepoint",
        shape = "Replicate"
    ) +
    
    theme_classic()

# Save plot
ggsave(
    filename = "PCA_mislabelled_samples.png",
    path = file.path(working_dir, "7-output", "figures"),
    width = 8, height = 5, dpi = 300
)
