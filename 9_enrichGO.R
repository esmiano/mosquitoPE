#===============================================================================
# GENE-ONTOLOGY (GO) ENRICHMENT ANALYSIS
#===============================================================================
# Author: esmiano
# Version: 1.0
# Date created: 24-AUG-2026
# 
# Description:
# This script performs GO enrichment of the genes in Ir68a's WGCNA module using 
# two separate annotation sources, these being VectorBase via g:Profiler (GOSt) 
# and an NCBI-derived OrgDb via clusterProfiler (enrichGO / gseGO).
#
# Types of GO enrichment analyses:
#  1. Over-Representation Analysis (ORA)
#  2. Gene Set Enrichment Analysis (GSEA)
# 
# The Adavi annotation file uses mixed Gene IDs:
#  1. Ae. aegypti symbol alone
#  2. VectorBase (AAEL) IDs alone
#  3. VectorBase ID with D. melanogaster suffix (AAEL_Dmelname format)

#===============================================================================
# SETUP AND PREPARATION
#===============================================================================

# Packages ---------------------------------------------------------------------

# Load packages
library(tidyverse)
library(gprofiler2)
library(ltc)
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
  filename = file.path(figures_dir, "GOSt_plot.png"),
  width = 10, height = 6
)

# Generate results table
publish_gosttable(
  gostRes,
  use_colors = TRUE,
  filename = file.path(figures_dir, "GOSt_table.png")
)

# Export results
gost_export <- gostRes$result
gost_export$parents <- sapply(gost_export$parents, paste, collapse = ",")
write.csv(gost_export, file.path(results_dir, "GOSt_results.csv"), 
    row.names = FALSE)


#===============================================================================
# CLUSTERPROFILER
#===============================================================================

# clusterProfiler has a wider variety of plotting options for more effective 
# visualisation of the results. However, while it does have options for 
# inputting a custom set GO annotations (e.g. retrieved from VectorBase using 
# g:GOSt as above), it works best using an OrgDb built off NCBI annotations.

# Retrieve OrgDb from AnnotationHub --------------------------------------------

# There isn't an OrgDb package on Bioconductor but AnnotationHub hosts a more
# extensive list of databases that include non-typical organisms.

ah <- AnnotationHub()
orgdb <- ah[["AH121343"]]   # Ae. aegypti OrgDb record

# Convert VectorBase IDs to Entrez ---------------------------------------------------

# The Ae. aegypti OrgDb is indexed by Entrez ID so need to convert gene lists so 
# they are compatible

# Convert query gene list
query_conv <- gconvert(
  query     = query,
  organism  = "aalvpagwg",
  target    = "ENTREZGENE_ACC",
  filter_na = TRUE
)

# Convert background gene list
bg_conv <- gconvert(
  query     = background,
  organism  = "aalvpagwg",
  target    = "ENTREZGENE_ACC",
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
  gene          = query_entrez, # Ir68a module genes
  universe      = bg_entrez,    # All tested genes
  OrgDb         = orgdb,        # Aedes aegypti GO database
  keyType       = "ENTREZID",   # Gene identification key
  ont           = "ALL",        # BP, CC, and MF
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  pAdjustMethod = "BH",         # Benjamini-Hochberg p value adjustment method
  minGSSize     = 10,
  maxGSSize     = 500,
  readable      = TRUE          # Retrieve readable gene names
)

# Reduce redundnacy by computing semantic similarity and collapsing redundant
# terms (keeps most significant representative of each cluster)
simGO <- simplify(eGO, cutoff = 0.7, by = "p.adjust", select_fun = min)

# ORA dot plot (include up to 20 GO terms)
p_dot <- dotplot(simGO, showCategory = 20, split = "ONTOLOGY", label_format = 60) +
  # Facet by ontology
  facet_grid(ONTOLOGY ~ ., scales = "free") +
  
  # Customise appearance
  theme(
    strip.text       = element_blank(),
    # Force blank background
    strip.background = element_blank(),
    panel.background = element_blank(),
    panel.grid       = element_blank(),
    # Create more space on left margins of plot for GO term labels
    plot.margin      = margin(5, 5, 5, 25),
    # Increase GO term label text size
    axis.text        = element_text(size=13),
    # Create thick border around each facetted plot
    panel.border     = element_rect(colour = "black", fill = NA, linewidth = 1),
  )

# Save dot plot
ggsave(file.path(figures_dir, "ORA_dot.png"), plot = p_dot,
       width = 12, height = 10, dpi = 300)


# Manhattan plot (GO enrichment landscape)
p_man <- manhattanplot(eGO,
  color        = "p.adjust",  # Colour according to adjusted p value
  showCategory = 12,          # Show 12 GO terms
  size         = "Count",     # Size indicates number of associated genes
  split        = "ONTOLOGY") +
  theme_classic()

# Save manhattan plot
ggsave(file.path(figures_dir, "ORA_manhattan.png"), plot = p_man,
       width = 9, height = 10, dpi = 300)


# Perform pairwise similarity on simplified terms (required for ernichment map)
pwsimGO <- pairwise_termsim(simGO)

# ORA enrichment map
p_emap <- emapplot(pwsimGO, node_label_size = 3, size_edge = 0.5, showCategory = 25) +
  # Scale fill map node colour according to adjusted p value
  scale_fill_continuous(low = "#e06663", high = "#327eba", name = "p.adjust",
                # Add legend
                guide = guide_colorbar(reverse = TRUE, order = 1), trans='log10')

# Save enrichment map
ggsave(file.path(figures_dir, "ORA_emap.png"), plot = p_emap,
       width = 9, height = 10, dpi = 300)


# Plot concept-gene network (connects which genes drive which terms)
p_cnet <- cnetplot(
  simGO,
  node_label     = "all",
  color_category = "firebrick",
  color_item     = "steelblue"
)

# Save concept-gene network
ggsave(file.path(figures_dir, "ORA_cnet.png"), plot = p_cnet,
       width = 15, height = 8, dpi = 300)


# Gene Set Enrichment Analysis (GSEA) ------------------------------------------

# GSEA reveals which GO terms are enriched towards the top or bottom of a 
# pre-ranked list of all tested genes. Genes are ranked according to shrunken
# Log2FoldChange and the 96h vs 2h contrast is used here to examine 
# post-emergence maturation processes occuring as Ir68a is downregulated.

# Contrasts
# Male   = "male_96v02"
# Female = "female_96v02"

# Define function for performing GSEA
performGSEA <- function(GSEA_contrast) {

  # Retrieve genes to be tested in GSEA
  ranks_df <- shrunk_final %>%
    # Filter for specific contrasst
    filter(contrast == GSEA_contrast) %>%
    # Rename gene ID column for downstream table joining
    dplyr::rename(input = gene_id)

  # Convert filtered gene IDs to Entrez format
  ranks_conv <- gconvert(
    query     = ranks_df$input,
    organism  = "aalvpagwg",
    target    = "ENTREZGENE_ACC",
    filter_na = TRUE
  )

  # Filter and rank genes
  ranks_df <- ranks_df %>%
    # Merge Entrez ID conversion and ranking dataframe, joining by Gene ID
    dplyr::inner_join(ranks_conv, by = "input") %>%
    # Remove columns not required for ranking
    dplyr::select(target, log2FoldChange) %>%
    # Keep only the max log2fc value for each Entrez ID
    slice_max(abs(log2FoldChange), by = target, n = 1, na_rm = TRUE, with_ties = FALSE) %>%
    # Explicitly ungroup
    ungroup() %>%
    # Arrange in order of descending Wald Stat
    arrange(desc(log2FoldChange))
    
  # Convert to named vector, ready for gseGO
  ranks <- setNames(ranks_df$log2FoldChange, ranks_df$target)

  # Run gseGO
  gseaGO <- gseGO(
    geneList     = ranks,        # Ranked list of genes
    OrgDb        = orgdb,        # Aedes aegypti GO database
    ont          = "ALL",        # BP, CC, and MF
    keyType      = "ENTREZID",   # Gene identification key
    minGSSize    = 10,
    maxGSSize    = 500,
    pvalueCutoff = 0.05,
    verbose      = FALSE
  )

  # Simplify GSEA GO terms
  simgGO <- simplify(gseaGO, cutoff = 0.7, by = "p.adjust", select_fun = min)
  
  # Plot ridge plot
  p_ridge <- ridgeplot(simgGO, showCategory = 15, label_format = 60) +
    # Widen the left margin to give terms more space
    theme(plot.margin = margin(5, 5, 5, 25)) +
    theme_classic()
  
  # Save figure
  ggsave(file.path(figures_dir, paste0("GSEA_", GSEA_contrast, ".png")),
    plot = p_ridge, width = 10, height = 7, dpi = 300)
}

# Perform GSEA on 96h vs 2h samples
gsea_male   <- performGSEA("male_96v02")
gsea_female <- performGSEA("female_96v02")

#===============================================================================
# RESULTS EXPORT
#===============================================================================

# Save ORA results
write.csv(as.data.frame(eGO),
          file.path(results_dir, "ORA_results.csv"), row.names = FALSE)

# Save male GSEA results
write.csv(as.data.frame(gsea_male),
          file.path(results_dir, "GSEA_male_results.csv"), row.names = FALSE)

# Save female GSEA results
write.csv(as.data.frame(gsea_female),
          file.path(results_dir, "GSEA_female_results.csv"), row.names = FALSE)
