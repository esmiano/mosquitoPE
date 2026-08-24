
# Script function: Pairwise correlation and co-expression network analysis 
#   using WGCNA
# Version: 1.0
# Date created: 23-AUG-2026
# 
# Dependencies:
#  - tidyverse ()
#    
#  - ggplot2 (Plotting figures)
#    
#  - pheatmap (Plotting heatmaps)
#    
#  - RColorBrewer (Colour palettes)
#    
#  - WGCNA ()
#


#===============================================================================
# SETUP
#===============================================================================

# Load packages
library(tidyverse)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(WGCNA)
library(flashClust)

#===============================================================================
# PAIRWISE CORRELATION ANALYSIS
#===============================================================================

# Spearman correlation ---------------------------------------------------------

# Filter normalised counts
matrix_counts <- norm_counts[rownames(norm_counts) %in% genes_of_interest, ]

# Perform Spearman correlation
# t() transposes counts so cor() correlates genes instead of samples
cor_matrix <- cor(t(matrix_counts), method = "spearman")

# Correlation heatmap 
pheatmap(
    cor_matrix,
    color = colorRampPalette(c("blue", "white", "red"))(100)
)

# Genome-wide co-expression with Ir68a ----------------------------------------

# TODO: Potentially remove this section
# Correlate every gene with Ir68a
ir68a_counts <- norm_counts["Ir68a", ]
ir68a_cors <- apply(norm_counts, 1, function(x) cor(x, ir68a_counts, method = "spearman"))

# Sort by absolute correlation
ir68a_cors_sorted <- sort(abs(ir68a_cors), decreasing = TRUE)

# Top 200 co-expressed genes (excluding Ir68a itself)
top_coexpressed <- names(ir68a_cors_sorted)
top_coexpressed <- top_coexpressed[top_coexpressed != "Ir68a"]
top_coexpressed <- head(top_coexpressed, 200)

# Heatmap of top 50 co-expressed genes
top50 <- c("Ir68a", head(top_coexpressed, 50))
pheatmap(norm_counts[top50, ],
         scale = "row",
         show_colnames = FALSE,
         annotation_col = sample_info[, c("sex", "timepoint")],
         main = "Top 50 genes co-expressed with Ir68a",
         color = colorRampPalette(c("blue", "white", "red"))(50))

# Save for GO enrichment
saveRDS(top_coexpressed, file.path(results_dir, "ir68a_coexpressed_genes.rds"))

#===============================================================================
# WGCNA PREPARATION
#===============================================================================

# Retrieve VST data
vsd_matrix <- assay(vsd)

# Check for NA values if any are present, keep only complete cases
if (sum(is.na(vsd_matrix)) > 0) {
    vsd_matrix <- vsd_matrix[complete.cases(vsd_matrix), ]
}

# Removing low-variance genes --------------------------------------------------

# Low-variance genes provide little co-expression information while adding noise

# Calculate variance of each gene across samples
gene_var <- apply(vsd_matrix, 1, var)

# Inspect distribution before choosing a cutoff
hist(log10(gene_var), breaks = 50,
     main = "Distribution of gene variance",
     xlab = "log10(variance)")

# Remove lowest-variance genes
var_threshold <- quantile(gene_var, 0.3)
vsd_matrix_filt <- vsd_matrix[gene_var > var_threshold, ]
cat("genes retained:", nrow(vsd_matrix_filt), "\n")

# Transpose so rows = samples and columns = genes
vsd_matrix_tp <- t(vsd_matrix_filt)

# Outlier detection ------------------------------------------------------------

# Cluster samples by overall expression distance
sample_tree <- hclust(dist(vsd_matrix_tp), method = "average")

# Plot clustered samples to detect outliers
# Height represents dissimilarity between samples
plot(sample_tree,
     main = "Sample clustering to detect outliers",
     xlab = "", cex = 0.7)

# Plot data to see distribution
# Should be approximately normally distributed
hist(as.matrix(vsd_matrix_tp), breaks = 100,
     main = "Distribution of expression values",
     xlab = "Expression value")

#===============================================================================
# BUILDING THE WGCNA NETWORK
#===============================================================================

# Choosing the soft-thresholding power -----------------------------------------

# Test a range of candidate powers
powers <- c(seq(1, 10, by = 1), seq(12, 20, by =2))

sft <- pickSoftThreshold(
    vsd_matrix_tp,
    powerVector = powers,
    networkType = "signed",
    verbose = 3
)

par(mfrow = c(1, 2))

# Plot scale-free topology fit against power
plot(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft threshold (power)",
     ylab = "Scale-free topology model fit, signed R^2",
     main = "Scale independence", type = "n")
text(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = 0.9, col = "red")
abline(h = 0.85, col = "red", lty = 2)

# Plot mean connectivity against power
plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft threshold (power)",
     ylab = "Mean connectivity",
     main = "Mean connectivity", type = "n")
text(sft$fitIndices[, 1], sft$fitIndices[, 5],
     labels = powers, cex = 0.9, col = "red")
abline(h = 10, col = "red", lty = 2)
par(mfrow = c(1, 1))

# View as table
sft$fitIndices

# Assign soft-thresholding power
# 12 achieves scale independence ~0.83 with mean connectivity of 140
soft_power <- 12

# Constructing the network and detecting modules -------------------------------

# Use WGCNA cor() function for duration of construction
cor <- WGCNA::cor

# Construct network and detect modules
network <- blockwiseModules(
    vsd_matrix_tp,
    power = soft_power,
    networkType = "signed",
    TOMType = "signed",
    minModuleSize = 30,
    reassignThreshold = 0,
    mergeCutHeight = 0.25,
    numericLabels = TRUE,       # Returns numeric labels, converted to colours below
    pamRespectsDendro = FALSE,
    saveTOMs = TRUE,            # Writes TOM to disk
    saveTOMFileBase = "TOM",
    verbose = 3
)

# Restore standard cor() function
cor <- stats::cor

# Convert numeric module labels to colour names for visualisation
module_colr <- labels2colors(network$colors)

# Count genes assigned to each module
table(module_colr)
cat("Modules identified:", length(unique(module_colr)) - 1, "\n")

# Plot dendrogram
plotDendroAndColors(
    network$dendrograms[[1]],
    module_colr[network$blockGenes[[1]]],
    "Module colours",
    dendroLabels = FALSE,
    hang = 0.03,
    addGuide = TRUE,
    guideHang = 0.05,
    main = "Gene dendrogram and module colours"
)

# Module eigengenes ----

# Extract module eigengenes
module_eigen <- network$MEs

# Reorder modules so that similar ones are adjacent
module_eigen <- orderMEs(module_eigen)

# Plot heatmap of module behaviour across samples, annotated by condition
# TODO: Fix annotation column
pheatmap(
    t(module_eigen),
    scale = "row",
    clustering_distance_cols = "euclidean",
    clustering_method = "average",
    annotation_col = "????",
    show_colnames = FALSE,
    main = "Module eigengenes across samples",
    color = colorRampPalette(c("blue", "white", "red"))(50)
)

#===============================================================================
# CONNECTING WGCNA MODULES TO BIOLOGY
#===============================================================================

