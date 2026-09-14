#===============================================================================
# WGCNA: GENE CO-EXPRESSION ANALYSIS
#===============================================================================
# Author: esmiano
# Version: 1.0
# Date created: 23-AUG-2026
# 
# Description:
# This script performs a pairwise correlation analysis of the genes of interest
# as well as a more in-depth co-expression network analysis of the full gene set
# using WGCNA.

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
library(viridis)

#===============================================================================
# INPUT PREPARATION
#===============================================================================

# Retrieve VST expression data
expr_data <- assay(vsd)

# Low-variance gene removal ----------------------------------------------------

# Low-variance genes provide little co-expression information while adding noise

# Calculate variance of each gene across samples
gene_var <- apply(expr_data, 1, var)

# Inspect distribution before choosing a cutoff
hist(log10(gene_var), breaks = 50,
     main = "Distribution of gene variance",
     xlab = "log10(variance)")

# Remove lowest-variance genes
var_threshold <- quantile(gene_var, 0.3)
expr_data_filt <- expr_data[gene_var > var_threshold, ]

# Check genes before and after filtering
message("\nInitial genes: ", nrow(expr_data), "\n")
message("\nGenes after filtering: ", nrow(expr_data_filt), "\n")

# Transpose so rows = samples and columns = genes (required by WGCNA)
datExpr <- t(expr_data_filt)

# Quality control --------------------------------------------------------------

# Perform QC
# Flags samples/genes with too many missing values or zero variance
gsg <- goodSamplesGenes(datExpr)

# If any samples or genes are bad or missing, then add okay ones to datExpr
if (!gsg$allOK) {datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]}

# Check genes remaining after QC
message("\nGenes after QC: ", ncol(datExpr), "\n")

# Align sample metadata to dataset
sample_info_wgcna <- sample_info[rownames(datExpr), ]

# Outlier detection ------------------------------------------------------------

# Cluster samples by overall expression distance
sample_tree <- hclust(dist(datExpr), method = "average")

# Add sex/timepoint information to the sample tree
tree_labels <- paste0(sample_info_wgcna$sex, "_", sample_info_wgcna$timepoint, "_", rownames(sample_info_wgcna))

# Assign labels to the dendrogram
sample_tree$labels <- tree_labels

# Plot clustered samples to detect outliers
# Height represents dissimilarity between samples
plot(
    sample_tree,
    main = "Sample clustering to detect outliers",
    xlab = "",
    sub  = "",
    cex  = 0.7
)

# Plot data to see distribution
hist(as.matrix(datExpr), breaks = 100,
     main = "Distribution of expression values",
     xlab = "Expression value")

# Trait data matrix construction -----------------------------------------------

# WGCNA correlates module eigengenes against numeric traits, so it requires 
# sample information to be numerically encoded in datTrait matrix

# Create trait data matrix
datTrait <- data.frame(
    # as.numeric() labels a cell with "1" if the condition is true for that 
    # sample, otherwise cell is labelled as "0"
    "Female" = as.numeric(sample_info_wgcna$sex == "Female"),
    "Male"   = as.numeric(sample_info_wgcna$sex == "Male"),
    "2h"     = as.numeric(sample_info_wgcna$timepoint == "02"),
    "12h"    = as.numeric(sample_info_wgcna$timepoint == "12"),
    "24h"    = as.numeric(sample_info_wgcna$timepoint == "24"),
    "48h"    = as.numeric(sample_info_wgcna$timepoint == "48"),
    "96h"    = as.numeric(sample_info_wgcna$timepoint == "96"),
    row.names = rownames(sample_info_wgcna),
    check.names = FALSE   # Keeps column names as inputted
)

# Soft-thresholding power choice -----------------------------------------------

# Test a range of candidate powers
powers <- c(seq(1, 10, by = 1), seq(12, 20, by =2))

# Initialise soft threshold exploration object
sft <- pickSoftThreshold(
    datExpr,
    powerVector = powers,
    networkType = "signed",
    verbose     = 3
)

# Set up multi-plot (2 plots side by side)
par(mfrow = c(1, 2))

# Plot scale-free topology fit against power
plot(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft threshold (power)",
     ylab = "Scale-free topology model fit, signed R^2",
     main = "Scale independence",
     type = "n")
text(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers,
     cex    = 0.9,
     col    = "red")
abline(h = 0.85, col = "red", lty = 2)

# Plot mean connectivity against power
plot(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     xlab = "Soft threshold (power)",
     ylab = "Mean connectivity",
     main = "Mean connectivity",
     type = "n")
text(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     labels = powers,
     cex    = 0.9,
     col    = "red")
abline(h = 10, col = "red", lty = 2)
par(mfrow = c(1, 1))

# View results as table
sft$fitIndices

# Assign soft-thresholding power
# 14 achieves scale independence ~0.86 with mean connectivity of 100
soft_power <- 14

#===============================================================================
# NETWORK CONSTRUCTION AND MODULE DETECTION
#===============================================================================

# Network construction and module detection ------------------------------------

# Use WGCNA cor() function for duration of construction
cor <- WGCNA::cor

# Construct network and detect modules
#   blockwiseModules builds the correlation network, computes the topological
#   overlap matrix (TOM), clusters genes, and cuts the tree into modules
network <- blockwiseModules(
    datExpr,
    maxBlockSize      = 46000,         # Force gene block size (max WGCNA can handle)
    power             = soft_power,    # Previously determined soft power
    networkType       = "signed",      # +ve and -ve corr genes not grouped together
    TOMType           = "signed",      # TOM follows same logic as above
    minModuleSize     = 30,            # Minimum genes per module
    reassignThreshold = 0,             # Reassignment of genes disabled
    mergeCutHeight    = 0.25,          # Merge modules with eigengene corr > 0.75
    numericLabels     = TRUE,          # Retrieve numeric module labels
    pamRespectsDendro = FALSE,         # Dendrogram ignored during PAM stage
    saveTOMs          = TRUE,
    saveTOMFileBase   = "TOM",
    verbose           = 3
)


# Restore standard cor() function
cor <- stats::cor

# Convert numeric module labels to colour names for visualisation
module_colours <- labels2colors(network$colors)

# Count genes assigned to each module
module_sizes <- table(module_colours)

# Open png device
png(file.path(figures_dir, "WGCNA_dendrogram.png"),
    width = 10, height = 4, units = "in", res = 300)
# Set figure margins
par(mar = c(8, 10, 3, 3))
# Plot gene dendrogram with module colours
plotDendroAndColors(
    dendro = network$dendrograms[[1]],                # View first (only) gene block
    colors = module_colours[network$blockGenes[[1]]], # Module colour per gene
    groupLabels  = "Module colours",                  # Label for the colour bar row
    dendroLabels = FALSE,    # Hide individual gene labels
    hang         = 0.03,     # How far leaves hang below branches
    addGuide     = TRUE,     # Guide lines from dendrogram to colour bar
    guideHang    = 0.05,     # Gap between guide lines and dendrogram
    main         = ""        # Plot title (blank for export)
)
# Turn off png device
dev.off()

# Module eigengenes ------------------------------------------------------------

# Extract module eigengenes
module_eigengenes <- network$MEs

# Reorder modules so that similar ones are adjacent
module_eigengenes <- orderMEs(module_eigengenes)

# Map eigengene number column names to colour names
eigengene_numbers <- as.numeric(gsub("ME", "", colnames(module_eigengenes)))
colnames(module_eigengenes) <- paste0("ME", labels2colors(eigengene_numbers))

# Build figure annotation dataframe
# Ensures numeric traits grouped together into two annotation legends 
# (sex/timepoint) rather than separate legends for each trait
fig_annotation <- sample_info_wgcna %>%
    # Keep only sex and timepoint columns from sample info dataframe
    select(sex, timepoint) %>%
    # Convert to factors
    mutate(across(c(sex, timepoint), factor)) %>%
    # Rename columns for cleaner legends
    dplyr::rename(Sex = sex, Timepoint = timepoint)

# Link colour palettes to biological traits
annotation_colours <- list(Sex = sex_palette, Timepoint = timepoint_palette)

# Heatmap of module eigengenes across all samples, annotated by trait
pheatmap(
    t(module_eigengenes), # Transpose so modules are rows and samples are columns
    scale                    = "row",
    clustering_distance_cols = "euclidean",
    clustering_method        = "average",
    annotation_col           = fig_annotation,
    annotation_colors        = annotation_colours,
    cluster_cols             = TRUE,
    show_colnames            = FALSE,
    main   = "Module eigengenes across samples",
    color  = colorRampPalette(c("blue", "white", "red"))(50)
)

#===============================================================================
# CONNECTING MODULES TO BIOLOGY
#===============================================================================

# Module-trait correlation -----------------------------------------------------

# Correlate module eigengenes with datTrait
module_trait_corr <- cor(module_eigengenes, datTrait, use = "pairwise.complete.obs")

# Derive p-values from correlations and sample count
nSamples <- nrow(datExpr)
module_trait_pvalue <- corPvalueStudent(module_trait_corr, nSamples)

# Build cell labels including correlation and p-value
textMatrix <- paste(signif(module_trait_corr, 2), "\n(",
                    signif(module_trait_pvalue, 2), ")", sep = "")
dim(textMatrix) <- dim(module_trait_corr)

# Open png device
png(file.path(figures_dir, "WGCNA_module_trait_heatmap.png"),
    width = 10, height = 11, units = "in", res = 300)
# Set figure margins
par(mar = c(5, 10, 2, 3))
# Plot module-trait relationships
labeledHeatmap(
    Matrix        = module_trait_corr,
    xLabels       = names(datTrait),
    yLabels       = names(module_eigengenes),
    ySymbols      = names(module_eigengenes),
    colorLabels   = FALSE,
    colors        = blueWhiteRed(50),
    textMatrix    = textMatrix,
    setStdMargins = FALSE,
    zlim          = c(-1, 1),
    cex.text      = 0.8,
)
# Force png device to close
# While current device is greater than 1 switch off
while (dev.cur() > 1) dev.off()

# Ir68a module identification --------------------------------------------------

# Identify Ir68a module and convert to colour
ir68a_module <- labels2colors(network$colors[match("Ir68a", colnames(datExpr))])

# Define Ir68a module eigengene
ir68a_ME <- paste0("ME", ir68a_module)

# Determine module membership (genes correlated with Ir68a module eigengene)
module_membership <- cor(
    datExpr,
    module_eigengenes[, ir68a_ME],
    use = "pairwise.complete.obs"
)

# Convert to numeric vector and add corresponding Gene IDs to each element
module_membership <- as.numeric(module_membership)
names(module_membership) <- colnames(datExpr)

# Identification of genes co-expressed with Ir68a ------------------------------

# Correlate Ir68a with every other gene in the WGCNA dataset
ir68a_wgcna_corr <- cor(datExpr, datExpr[, "Ir68a"], use = "pairwise.complete.obs")

# Convert to numeric vector and add corresponding Gene IDs to each element
ir68a_wgcna_corr <- as.numeric(ir68a_wgcna_corr)
names(ir68a_wgcna_corr) <- colnames(datExpr)

# Create results table
ir68a_wgcna_res <- data.frame(
    gene                  = colnames(datExpr),
    module                = module_colours,
    corr_with_Ir68a       = ir68a_wgcna_corr,
    abs_corr_with_Ir68a   = abs(ir68a_wgcna_corr),
    module_membership     = module_membership,
    abs_module_membership = abs(module_membership)
)

# Filter and order results
ir68a_wgcna_res <- ir68a_wgcna_res %>%
  # Remove Ir68a itself and restrict remaining remaining genes to those sharing 
  # the same module
  filter(gene != "Ir68a", module == ir68a_module) %>%
  # Rank genes by direct correlation with Ir68a (descending)
  arrange(desc(abs_corr_with_Ir68a))

# Save results
write.csv(
    ir68a_wgcna_res,
    file.path(results_dir, "Ir68a_WGCNA_module_coexpression.csv"),
    row.names = FALSE
)

# Retrieve Ir68a module correlation with traits
ir68a_module_corr <- module_trait_corr[ir68a_ME, ]
ir68a_module_p    <- module_trait_pvalue[ir68a_ME, ]
data.frame(
    trait = names(ir68a_module_corr),
    r     = round(ir68a_module_corr, 3),
    p     = signif(ir68a_module_p, 3)
)

# Ir68a module gene expression -------------------------------------------------

# Include Ir68a itself (previously filtered out of ir68a_wgcna_res)
module_genes_w_ir68a <- c("Ir68a", ir68a_wgcna_res$gene)

# Gene expression heatmap
pheatmap(
    t(scale(datExpr[, module_genes_w_ir68a])),
    annotation_col = fig_annotation,
    annotation_colors = annotation_colours,
    cluster_cols   = TRUE,
    show_rownames  = FALSE,
    show_colnames  = FALSE,
    treeheight_row = 0,
    filename = file.path(figures_dir, paste0("WGCNA_module_expr_heatmap.png")),
    width = 10,
    height = 12,
    color  = colorRampPalette(c("black", rev(brewer.pal(11, "RdYlBu")), "black"))(100)
)
