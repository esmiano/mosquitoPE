
# Script function: Differential Gene Expression analysis using DESeq2
# Version: 1.0
# Date created: 18-AUG-2026
# 
# Dependencies:
#  - DESeq2 (Differential Gene Expression analysis)
#    
#  - ggplot2 (Plotting figures)
#    
#  - pheatmap (Plotting heatmaps)
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

# DESeq object setup -----------------------------------------------------------

# Initialise DESeq objecT
dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = sample_info,
    design = ~ 1  # Placeholder intercept design (overwritten below)
)

# Remove genes with low counts and update DESeq object
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

# Define sex x timepoint condition coefficients
dds$condition <- factor(paste0(dds$sex, "_", dds$timepoint))

# Set level order explicitly
dds$condition <- factor(dds$condition, levels = c(
  "female_02", "female_12", "female_24", "female_48", "female_96",
  "male_02",   "male_12",   "male_24",   "male_48",   "male_96"
))

# Update model design
# No intercept design means each coefficient is a mean, rather than comparing
# with a reference level
design(dds) <- ~ 0 + condition

# Wald Test --------------------------------------------------------------------

# The Wald Test estimates and tests the significance of specific contrasts 
# between explicitly named coefficients.

# Run standard DESeq2 for Wald Test
dds <- DESeq(dds)

# Define contrasts
all_contrasts <- list(
    # Sex effect at each timepoint (female used as reference)
    list(label = "sex_at_02h", contrast = c("condition", "male_02", "female_02")),
    list(label = "sex_at_12h", contrast = c("condition", "male_12", "female_12")),
    list(label = "sex_at_24h", contrast = c("condition", "male_24", "female_24")),
    list(label = "sex_at_48h", contrast = c("condition", "male_48", "female_48")),
    list(label = "sex_at_96h", contrast = c("condition", "male_96", "female_96")),

    # Timepoint effect in females (female at 2hPE used as reference)
    list(label = "female_12v02", contrast = c("condition", "female_12", "female_02")),
    list(label = "female_24v02", contrast = c("condition", "female_24", "female_02")),
    list(label = "female_48v02", contrast = c("condition", "female_48", "female_02")),
    list(label = "female_96v02", contrast = c("condition", "female_96", "female_02")),

    # Timepoint effect in males (male 2hPE used as reference)
    list(label = "male_12v02", contrast = c("condition", "male_12", "male_02")),
    list(label = "male_24v02", contrast = c("condition", "male_24", "male_02")),
    list(label = "male_48v02", contrast = c("condition", "male_48", "male_02")),
    list(label = "male_96v02", contrast = c("condition", "male_96", "male_02"))
)

# Create empty results data.frames 
res_wald_df <- data.frame()
res_wald_shrunken_df <- data.frame()

# Extract results and apply Log Fold Change Shrinkage
for (i in 1:length(all_contrasts)) {
  # Extract i-th contrast
  contrast <- all_contrasts[[i]]
  # Extract results for contrast
  res <- results(dds, contrast = contrast$contrast)
  # Apply ashr shrinkage
  res_shrunken <- lfcShrink(dds, contrast = contrast$contrast, type = "ashr")
  
  # Build temporary results data.frame (overwritten by each i-th iteration)
  tmp_res_df <- as.data.frame(res) %>%
    # Add gene to gene ID column
    rownames_to_column("gene_id") %>%
    # Add contrast label
    mutate(contrast = contrast$label)

  # Build temporary shrunken results data.frame (as above)
  tmp_shrunken_df <- as.data.frame(res_shrunken) %>%
    # Add gene to gene ID column
    rownames_to_column("gene_id") %>%
    # Add contrast label
    mutate(contrast = contrast$label)

  # Append i-th contrast results to results data.frames
  res_wald_df <- rbind(res_wald_df, tmp_res_df)
  res_wald_shrunken_df <- rbind(res_wald_shrunken_df, tmp_shrunken_df)
}

res_wald_df_ord <- res_wald_df[order(res_wald_df$padj), ]
res_wald_shrunk_df_ord <- res_wald_shrunken_df[order(res_wald_shrunken_df), ]


# Likelihood Ratio Test (LRT) -------------------------------------------------

# The LRT compares the full model (which can represent different time courses
# per sex) against a reduced model (which forces the same time course for both
# sexes). All timepoints are effectively pooled into a single test per gene, 
# making it possible to test single variables (such as the interaction between
# sex and timepoint) across the experiment timecourse.

# Run DESeq2 comparing full vs reduced model for LRT
# Log2FC here defaults to arbitrary contrast and can't be interpreted as 
# interaction effect size
dds_lrt <- dds
design(dds_lrt) <- ~ sex + timepoint + sex:timepoint        # Full model design
dds_lrt <- nbinomLRT(dds_lrt, reduced = ~ sex + timepoint)  # Refit with LRT

# Extract and order LRT results
res_lrt <- results(dds_lrt)
res_lrt_ord <- res_lrt[order(res_lrt$padj), ]


# Interaction effect sizes -----------------------------------------------------



#===============================================================================
# VARIANCE STABILISING TRANSFORMATION (VST)
#===============================================================================

# VST normalises count data for visualising (PCA, heatmaps, distances)
# blind = FALSE uses the design formula to inform the transformation
vsd <- vst(dds, blind = FALSE)

#===============================================================================
# PRINCIPAL COMPONENT ANALYSIS (PCA)
#===============================================================================

# Extract PCA coordinated from the VST data
pca_data <- plotPCA(
    vsd,
    intgroup = c("sex", "timepoint", "replicate"),  # Grouping variables
    returnData = TRUE                               # Return PCA data.frame
)

# Calculate percentage variance
percent_var <- round(100 * attr(pca_data, "percentVar"))

# Define colourblind-friendly palette
colourblind_timepoint <- c(
    "02" = "#D55E00",
    "12" = "#CC79A7",
    "24" = "#0072B2",
    "48" = "#F0E442",
    "96" = "#009E73"
)

# PCA plot
ggplot(pca_data, aes(x = PC1, y = PC2)) +
  
  # Females (filled)
  geom_point(
    data = subset(pca_data, sex == "female"),
    aes(colour = timepoint, fill = timepoint, shape = replicate),
    size = 2.5,
    stroke = 1
  ) +
  
  # Males (empty)
  geom_point(
    data = subset(pca_data, sex == "male"),
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
  xlab(paste0("PC1: ", percent_var[1], "% variance")) +
  ylab(paste0("PC2: ", percent_var[2], "% variance")) +

  # Figure legen customisation
  labs(
    colour = "Timepoint",
    fill = "Timepoint",
    shape = "Replicate"
  ) +
  
  theme_classic()

    
#===============================================================================
# SAMPLE DISTANCES
#===============================================================================

# TODO: Decide whether to keep

# Euclidean --------------------------------------------------------------------

# Compute pairwise Euclidean distances between samples (VST data)
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

pheatmap(
  sampleDistMatrix,
  clustering_distance_rows = sampleDists,
  clustering_distance_cols = sampleDists,
  col = colors,
  main = "Sample distances (Euclidean)"
)

# Poisson ----------------------------------------------------------------------

# Compute pairwise Poisson distances between samples (raw counts)
poisd <- PoissonDistance(t(counts(dds)))
samplePoisDistMatrix <- as.matrix(poisd$dd)
rownames(samplePoisDistMatrix) <- paste(vsd$condition)
colnames(samplePoisDistMatrix) <- NULL

pheatmap(
  samplePoisDistMatrix,
  clustering_distance_rows = poisd$dd,
  clustering_distance_cols = poisd$dd,
  col = colors,
  main = "Sample distances (Poisson)"
)

#===============================================================================
# MINUS VS. AVERAGE (MA)
#===============================================================================

# MA plots show the relationship between mean expression (x-axis) and log2 fold
# change (y-axis).

# Iterate over each contrast
for (i in 1:length(all_contrasts)) {
  # Extract i-th contrast
  contrast <- all_contrasts[[i]]
  # Open PNG device
  png(file.path(figures_dir, "ma_plots", paste0("ma_", contrast$label, ".png")),
      width = 800, height = 600)
  # Plot MA for each contrast
  plotMA(results(dds, contrast = contrast$contrast),
         main = paste("MA plot:", contrast$label))
  # Close PNG device to write file
  dev.off()
}

#===============================================================================
# EXPRESSION PLOTS
#===============================================================================

# Define and create output directory for expression plots
expression_dir <- file.path(figures_dir, "expression_plots")
dir.create(expression_dir, showWarnings = FALSE, recursive = TRUE)



#===============================================================================
# NORMALISED COUNTS TABLES
#===============================================================================

# Full table -------------------------------------------------------------------

# Retrieve normalised counts
# normalized = TRUE adjusts for differences in sequencing depth between samples
norm_counts <- counts(dds, normalized = TRUE)

# Convert to dataframe and add Gene ID column
norm_counts_df <- as.data.frame(norm_counts) %>% rownames_to_column("gene_id")

# Join results and normalised counts data.frames according to "gene_id"
master_final <- left_join(res_wald_df_ord, norm_counts_df, by = "gene_id")

# Export CSV
# row.names = FALSE revents duplication of the gene ID column
write.csv(master_final, "DEG_results.csv", row.names = FALSE)

# Genes of interest only -------------------------------------------------------

goi_final <- filter(master_final, gene_id %in% genes_of_interest)
write.csv(goi_final, file.path(working_dir, "GOI_results.csv"), row.names = FALSE)