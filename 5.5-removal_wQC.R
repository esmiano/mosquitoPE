#===============================================================================
# SAMPLE REMOVAL
#===============================================================================
# Author: esmiano
# Version: 1.0
# Date created: 18-AUG-2026
# 
# Description:
# This script performs initial DESeq analysis and QC, as well as removal of 
# potentially mislabelled samples.

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

dds_og <- DESeqDataSetFromMatrix(
    countData = counts_og,
    colData = sample_info_og,
    design = ~ 1  # Placeholder intercept design (overwritten below)
)

# Remove genes with low counts and update DESeq object
keep <- rowSums(counts(dds_og)) >= 10
dds_og <- dds_og[keep, ]

# Define sex x timepoint condition coefficients
dds_og$condition <- factor(paste0(dds_og$sex, "_", dds_og$timepoint))

# Set level order explicitly
dds_og$condition <- factor(dds_og$condition, levels = c(
  "Female_02", "Female_12", "Female_24", "Female_48", "Female_96",
  "Male_02",   "Male_12",   "Male_24",   "Male_48",   "Male_96"
))

# Update model design
# No intercept design means each coefficient is a mean, rather than comparing
# with a reference level
design(dds_og) <- ~ 0 + condition

# Run standard DESeq2 for Wald Test
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

# PCA plot
ggplot(pca_data_og, aes(x = PC1, y = PC2)) +
  
  # Females (filled)
  geom_point(
    data = subset(pca_data_og, sex == "Female"),
    aes(colour = timepoint, fill = timepoint, shape = replicate),
    size = 2.5,
    stroke = 1
  ) +
  
  # Males (empty)
  geom_point(
    data = subset(pca_data_og, sex == "Male"),
    aes(colour = timepoint, shape = replicate),
    fill = "white",
    size = 2.5,
    stroke = 1
  ) +
  
  # Add labels to each sample point
  geom_text_repel(aes(label = name, colour = timepoint), size = 2.5) +
  
  # Stroke colour
  scale_colour_manual(
    # Add previously defined colourblind palette
    values = colourblind_timepoint,
    # Change timepoint labels in figure legend
    labels = c(
        "02" = "2h",
        "12" = "12h",
        "24" = "24h",
        "48" = "48h",
        "96" = "96h"
    )
  ) +
  
  # Fill colour
  scale_fill_manual(
    # Add previously defined colourblind palette
    values = colourblind_timepoint,
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
  
  # Add percentage variance to axes labels
  xlab(paste0("PC1: ", percent_var_og[1], "% variance")) +
  ylab(paste0("PC2: ", percent_var_og[2], "% variance")) +

  # Figure legen customisation
  labs(
    colour = "Timepoint",
    fill = "Timepoint",
    shape = "Replicate"
  ) +
  
  theme_classic() +
  
  # Remove figure legend
  theme(legend.position = "none")

# Save plot
ggsave(
    filename = "pca_before.png",
    path = file.path(working_dir, "7-output", "figures"),
    width = 8, height = 5, dpi = 300
)
