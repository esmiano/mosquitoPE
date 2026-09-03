#===============================================================================
# GENE-ONTOLOGY (GO) ENRICHMENT ANALYSIS
#===============================================================================
# Author: esmiano
# Version: 1.0
# Date created: 24-AUG-2026
# 
# Description:
# This script retrieves VectorBase IDs for each gene, where available, and 
# performs GO enrichment analysis on the genes in Ir68a's WGCNA module. The GO enrichment 
# tools used in the following analysis require VectorBase IDs for each gene as input.
# 
# The Adavi annotation file uses mixed Gene IDs:
#  1. Chemoreceptors are indicated by Ae. aegypti symbol
#  2. Genes not identified here are indicated by VectorBase (AAEL) IDs only
#  3. Genes with D. melanogaster homologues identified here are indicated by 
#     AAEL_Drosophilasymbol
# 
# No other previously identified genes of interest fall in the same WGCNA module 
# as Ir68a, so their lack of VectorBase IDs isn't an issue.

#===============================================================================
# SETUP AND PREPARATION
#===============================================================================

# Packages ---------------------------------------------------------------------

# Load packages
library(tidyverse)
library(gprofiler2)
library(clusterProfiler)
library(enrichplot)
library(AnnotationHub)
library(AnnotationDbi)
library(GO.db)

# Gene list extraction ---------------------------------------------------------

# Strip Drosophila name suffixes from VectorBase IDs where present
background <- gsub("_.*", "", colnames(datExpr))          # All genes tested
query      <- gsub("_.*", "", rownames(ir68a_wgcna_res))  # Ir68a module genes

#===============================================================================
# GPROFILER2
#===============================================================================

# g:Profiler takes mixed gene ID types as input so can accept both VectorBase 
# IDs and gene names. It retrieves GO annotations from VectorBase.

# g:GOSt functional profiling ----------------------------------------------------

# Run g:GOSt
gostRes <- gost(
  query             = query,          # Module gene list
  organism          = "aalvpagwg",    # Ae. aegypti strain LVP_AGWG
  custom_bg         = background,     # Background (all tested genes) 
  domain_scope      = "custom",       # Enforces custom background
  correction_method = "g_SCS",        # g:Profiler's own calibrated correction
  sources           = c("GO:BP", "GO:MF", "GO:CC"),
  significant       = TRUE,
  highlight         = TRUE
)

# Visualisation ----------------------------------------------------------------

# Plot "Manhatten-like" gost plot
publish_gostplot(
  gostplot(gostRes, capped = FALSE, interactive = FALSE),
  filename = file.path(figures_dir, "GO_gost_plot.png"),
  width = 10, height = 6
)

# Generate results table
publish_gosttable(
  gostRes,
  use_colors = TRUE,
  filename = file.path(figures_dir, "GO_gost_table.png")
)

# Export results
gost_export <- gostRes$result
gost_export$parents <- sapply(gost_export$parents, paste, collapse = ",")
write.csv(gost_export, file.path(results_dir, "GO_gost_res.csv"), 
    row.names = FALSE)


#===============================================================================
# CLUSTERPROFILER
#===============================================================================

# clusterProfiler has a wider variety of plotting options for more effective 
# visualisation of the results. However, while it does have options for 
# inputting a custom set GO annotations (e.g. retrieved from VectorBase using 
# g:GOSt as above), it works best using an OrgDb built off NCBI annotations.

# Retrieve OrgDb from AnnotationHub --------------------------------------------

# There isn't an OrgDb package on Bioconductor but AnnotationHub hosts a number 
# of such databases for non-typical organisms.

ah <- AnnotationHub()
orgdb <- ah[["AH121343"]]   # Ae. aegypti OrgDb record

# Convert VectorBase IDs to Entrez ---------------------------------------------------

# The Ae. aegypti OrgDb is indexed by Entrez ID so need to convert gene lists so 
# they are compatible

# Convert query gene list
query_conv <- gconvert(
  query = query,
  organism = "aalvpagwg",
  target = "ENTREZGENE_ACC",
  filter_na = TRUE
)

# Convert background gene list
bg_conv <- gconvert(
  query = background,
  organism = "aalvpagwg",
  target = "ENTREZGENE_ACC",
  filter_na = TRUE
)

# Assign converted IDs to new gene lists, only keeping unique IDs (no duplicates)
query_entrez <- unique(query_conv$target)
bg_entrez    <- unique(bg_conv$target)

# Over-Representation Analysis (ORA) -------------------------------------------

# ORA tests whether genes in the Ir68a WGCNA module are over-represented in any
# given GO terms compared to the background of all tested genes

# Run enrichGO on Ir68a module genes
eGO <- enrichGO(
  gene          = query_entrez,
  universe      = bg_entrez,
  OrgDb         = orgdb,
  keytype       = "ENTREZID",
  ont           = "ALL",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  minGSSize     = 10,
  maxGSSize     = 500,
  readable      = TRUE
)

# Reduce redundnacy by computing semantic similarity and collapsing redundant
# terms (keeps most significant representative of each cluster)
simGO <- simplify(eGO, cutoff = 0.7, by = "p.adjust", select_fun = min)

# Plot ORA dot plot
p_dot <- dotplot(simGO, showCategory = 10, split = "ONTOLOGY", label_format = 60) +
  facet_grid(ONTOLOGY ~ ., scales = "free")
p_dot
ggsave(file.path(figures_dir, ""))