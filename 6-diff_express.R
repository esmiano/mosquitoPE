#===============================================================================
# DIFFERENTIAL GENE EXPRESSION ANALYSIS
#===============================================================================
# Author: esmiano
# Version: 1.0
# Date created: 18-AUG-2026
# 
# Description:
# This script uses DESeq to identify Differentially Expressed Genes (DEGs) 
# within the previously defined genes of interest.

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
  "Female_02", "Female_12", "Female_24", "Female_48", "Female_96",
  "Male_02",   "Male_12",   "Male_24",   "Male_48",   "Male_96"
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
    list(label = "sex_at_02h", contrast = c("condition", "Male_02", "Female_02")),
    list(label = "sex_at_12h", contrast = c("condition", "Male_12", "Female_12")),
    list(label = "sex_at_24h", contrast = c("condition", "Male_24", "Female_24")),
    list(label = "sex_at_48h", contrast = c("condition", "Male_48", "Female_48")),
    list(label = "sex_at_96h", contrast = c("condition", "Male_96", "Female_96")),

    # Timepoint effect in females (female at 2hPE used as reference)
    list(label = "female_12v02", contrast = c("condition", "Female_12", "Female_02")),
    list(label = "female_24v02", contrast = c("condition", "Female_24", "Female_02")),
    list(label = "female_48v02", contrast = c("condition", "Female_48", "Female_02")),
    list(label = "female_96v02", contrast = c("condition", "Female_96", "Female_02")),

    # Timepoint effect in males (male 2hPE used as reference)
    list(label = "male_12v02", contrast = c("condition", "Male_12", "Male_02")),
    list(label = "male_24v02", contrast = c("condition", "Male_24", "Male_02")),
    list(label = "male_48v02", contrast = c("condition", "Male_48", "Male_02")),
    list(label = "male_96v02", contrast = c("condition", "Male_96", "Male_02"))
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

# Order results in orrder of padj value
res_wald_df_ord <- arrange(res_wald_df, padj)
res_wald_shrunk_df_ord <- arrange(res_wald_shrunken_df, padj)


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

# Extract LRT results
res_lrt <- results(dds_lrt)
# Convert to data.frame
res_lrt_ord_df <- as.data.frame(res_lrt) %>%
  # Create Gene ID column
  rownames_to_column("gene_id") %>%
  # Order according to padj value
  arrange(padj)

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
p_pca <- ggplot(pca_data, aes(x = PC1, y = PC2)) +
  
  # Females (filled)
  geom_point(
    data = subset(pca_data, sex == "Female"),
    aes(colour = timepoint, fill = timepoint, shape = replicate),
    size = 2.5,
    stroke = 1
  ) +
  
  # Males (empty)
  geom_point(
    data = subset(pca_data, sex == "Male"),
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
    fill   = "Timepoint",
    shape  = "Replicate"
  ) +
  
  theme_classic() +
  
  # Remove figure legend
  theme(legend.position = "none")

# Save pca plot without a legend
ggsave(
    filename = "pca_after.png", plot = p_pca,
    path = file.path(working_dir, "7-output", "figures"),
    width = 8, height = 5, dpi = 300
)

# PCA plot with legend
pca_with_legend <- p_pca + theme(legend.position = "bottom")

# Save PCA legend separately
ggsave("pca_legend.png", plot = ggdraw(get_legend(pca_with_legend)),
       path = figures_dir, width = 7, height = 0.5, dpi = 300)

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

# Define position dodge
dodge <- position_dodge(width=0.2)

# Colourblind-friendly palette
colourblind_sex <- c("Female" = "#E69F00", "Male" = "#56B4E9")

# Gene expression  -------------------------------------------------------------

# Define normalised gene expression plotting function
plot_expression <- function(gene_name) {

  # Extract normalised counts for gene from DESeq object
  gene_counts <- plotCounts(dds, gene = gene_name,
                            intgroup = c("sex", "timepoint"),
                            returnData = TRUE)
  
  # Initialise ggplot object
  p <- ggplot(gene_counts, aes(x = timepoint, y = count,
                                   colour = sex, group = sex)) +
    
    # Add line passing through mean values
    stat_summary(
      geom      = "line",
      fun       = mean,
      aes(group = sex),
      linewidth = 0.6,
      position  = dodge
    ) +
    
    # Add standard error bar
    stat_summary(
      geom      = "errorbar",
      fun.data  = mean_se,
      linewidth = 0.5,
      width     = 0.2,
      position  = dodge
    ) +
    
    # Add point for mean values
    stat_summary(
      geom      = "point",
      fun       = mean,
      aes(shape = sex),
      size      = 3,
      position  = dodge
    ) +

    # Set colours and shapes + capitalise labels
    scale_colour_manual(values = colourblind_sex) +
    scale_shape_manual(values = c(19, 17)) +

    # Add labels
    labs(
      # Axes labels
      x = "Hours post-emergence",
      y = "Normalised counts",
      # Legend labels
      colour = "Sex",
      shape = "Sex"
    ) +
    
    # Remove through line in legend shapes
    guides(color = guide_legend(override.aes = list(linetype = 0))) +
    
    # Set theme
    theme_classic() +

    # Remove figure legend
    theme(legend.position = "none")
}

# Plot and save figures with no legend
for (gene in no_legend) {
  ggsave(paste0(gene, "_expression.png"),
  # Plot to save and path to output directory
  plot = plot_expression(gene), path = figures_dir,
  # Plot dimensions
  width = 6, height = 4, dpi =300)
}

# Generate plots without legends
for (gene in c(hygroreceptors, co2_receptors)) {
  p <- plot_expression(gene)
  ggsave(paste0(gene, "_expression.png"), plot = plot_expression(gene),
         path = figures_dir, width = 5, height = 4, dpi = 300)
}

# Plot with legend
p_with_legend <- plot_expression("Ir93a") + theme(legend.position = "bottom")

# Save legend separately
ggsave("expression_legend.png", plot = ggdraw(get_legend(p_with_legend)),
       path = figures_dir, width = 3, height = 0.5, dpi = 300)


#===============================================================================
# NORMALISED COUNTS TABLES
#===============================================================================

# Full results -----------------------------------------------------------------

# Retrieve normalised counts
# normalized = TRUE adjusts for differences in sequencing depth between samples
norm_counts <- counts(dds, normalized = TRUE)

# Convert to dataframe and add Gene ID column
norm_counts_df <- as.data.frame(norm_counts) %>% rownames_to_column("gene_id")

# Extract LRT results columns
lrt_cols <- res_lrt_ord_df %>%
  dplyr::select(gene_id,
                lrt_stat   = stat,
                lrt_pvalue = pvalue,
                lrt_padj   = padj)


# Combine Wald results, LRT columns & normalised counts
master_final <- res_wald_df_ord %>%
  left_join(lrt_cols, by = "gene_id") %>%
  left_join(norm_counts_df, by = "gene_id") %>%
  # Rename Wald Test columns to be explicit
  dplyr::rename(wald_stat = "stat", wald_pvalue = "pvalue", wald_padj = "padj") %>%
  # Move contrast column to position after Gene ID
  relocate(contrast, .after = gene_id)

# Export CSV
# row.names = FALSE revents duplication of the gene ID column
write.csv(master_final, file.path(results_dir, "DEG_results.csv"), row.names = FALSE)

# Genes of interest only -------------------------------------------------------

# Filter gor genes of interest
goi_final <- filter(master_final, gene_id %in% genes_of_interest) %>%
  # Arrange in order of Gene ID first, then contrast
  arrange(gene_id, contrast)
write.csv(goi_final, file.path(results_dir, "GOI_results.csv"),
          row.names = FALSE)

# Filter for significant (padj<=0.05) Wald Test contrasts
sig_wald_goi_final <- filter(goi_final, wald_padj<=0.05) %>%
  dplyr::select(-lrt_stat, -lrt_pvalue, -lrt_padj)
write.csv(sig_wald_goi_final, file.path(results_dir, "Sig_Wald_GOI_results.csv"),
          row.names = FALSE)

# Filter for significant (padj<=0.05) LRT genes
sig_lrt_goi_final <- filter(goi_final, lrt_padj<=0.05) %>%
  dplyr::select(-contrast, -log2FoldChange, -lfcSE, -wald_stat, -wald_pvalue, -wald_padj) %>%
  distinct()
write.csv(sig_lrt_goi_final, file.path(results_dir, "Sig_LRT_GOI_results.csv"),
          row.names = FALSE)