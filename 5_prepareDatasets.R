#===============================================================================
# DATASET PREPARATION
#===============================================================================
# Author: esmiano
# Version: 1.1
# Date created: 09-AUG-2026
# 
# Description:
# This script imports gene expression counts, metadata, and gene annotations
# into R and prepares them for downstream analysis.

#===============================================================================
# CONFIGURATION & SETUP
#===============================================================================

# Clear environment to remove any artefacts
rm(list = ls())

# Configure path to working directory (requires user input)
working_dir <- "/Users/mspensley/Documents/university/research_project/analysis"

# Set working directory
setwd(file.path(working_dir))

# Configure annotations file name
annotations <- "Adavi2024_annotations.gtf"

# Configure paths to output directories
results_dir <- file.path(working_dir, "7-output", "results")
figures_dir <- file.path(working_dir, "7-output", "figures")

# Load packages
library(vroom)     # Fast data imports
library(tidyverse) # Data handling
library(ltc)       # Colour palettes

# Create separate output directories for results and figures
#   showWarnings set to FALSE to prevent warnings in case of pre-existing directory
#   recursive set to TRUE to create parent directories, if necessary
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

# Define timepoint palette
timepoint_palette <- ltc(kiss)[1:5]
names(timepoint_palette) <- c("02", "12", "24", "48", "96")

# Define sex colour palette
sex_palette <- c(
    "Female" = "#EE6C4D",
    "Male"   = "#3D5A80"
)

#===============================================================================
# COUNT DATA PREPARATION
#===============================================================================

# DESeq2 takes only numerical data as input so the output from featureCounts
# needs to be prepared accordingly.

# Import featureCounts output
count_data <- read.table(file.path(working_dir, "6-quant", "counts.txt"),
    # header = TRUE to indicate first line is a header
    header = TRUE,
    # Columns tab separated
    sep = "\t",
    # Replace quotes with the escaped equivalent
    quote = ""
)

# Set gene IDs as row names, ready for DESeq2
rownames(count_data) <- count_data[, 1]

# Remove Geneid column so only numeric counts remain
count_data <- dplyr::select(count_data, -Geneid)

# Extract just the SRR accession IDs from column names
# Replace column names with accessions for readability
# Regex pattern matches "SRR" followed by any digit then 1 or more of the preceding
# characters
colnames(count_data) <- str_extract(colnames(count_data), "SRR\\d+")

#===============================================================================
# METADATA PREPARATION
#===============================================================================

# Extracting condition info from the sample metadata is necessary to draw any 
# useful insights from the count data.

# Import the ENA metadata file mapping SRR accession IDs to library names
metadata <- vroom(file.path(working_dir, "reference", "PRJNA659517.tsv"))

# Parse sample information from the library name encoding and define explicitly:
#   Position 1     = sex (F/M)
#   Positions 2-3  = hours post-emergence (02, 12, 24, 48, 96)
#   For 12h only   = D (day) or N (night)
#   Next letter    = tissue (H = head, B = body)
#   Final position = biological replicate (1, 2, 3)
meta_parsed <- metadata %>%
  # Keep only head samples (ends in H followed by replicate number)
  filter(grepl("H[0-9]$", library_name)) %>%
  mutate(
    # Position 1: If "F", then label as "female", else "male"
    sex = ifelse(str_sub(library_name, 1, 1) == "F", "Female", "Male"),
    # Positions 2-3: Label timepoint
    timepoint = str_sub(library_name, 2, 3),
    # Final position: Label replicate
    replicate = str_sub(library_name, -1, -1)
  )

# Build sample info dataframe with accession IDs as row names
# DESeq2 uses this to reference which experimental group each sample belongs to
sample_info <- data.frame(
    sex = factor(meta_parsed$sex),
    timepoint = factor(meta_parsed$timepoint),
    replicate = factor(meta_parsed$replicate),
    row.names = meta_parsed$run_accession
)

# Reorder count matrix columns to match sample_info row order
counts <- count_data[, rownames(sample_info)]

# Verify alignment (stops script if counts and metadata don't match)
stopifnot(all(colnames(counts) == rownames(sample_info)))

#===============================================================================
# SAMPLE EXCLUSION
#===============================================================================

# Downstream PCA analysis (see 6_IDmislabelledSamples script) indicated that 2 
# pairs of samples may have been mislablled and effectively swapped. The 
# following were excluded from the final analysis:
#   SRR12522935 (Female, 2h PE, rep 2) and SRR12522925 (Female 24h PE, rep 2)
#   SRR12522892 (Male 12h (day), rep 2) and SRR12522931 (Female 12h (day), rep 2)

# Save original info and counts separately
sample_info_og <- sample_info
counts_og <- counts

# Remove potentially mislabelled samples
samples_to_remove <- c("SRR12522935", "SRR12522925", "SRR12522892", "SRR12522931")
sample_info <- sample_info[!rownames(sample_info) %in% samples_to_remove, ]
counts <- count_data[, rownames(sample_info)]

#===============================================================================
# GENES OF INTEREST
#===============================================================================

# Genes involved in the detection of humidity, CO2, and various Volatile Organic
# Compounds (VOCs) were identified from the literature. These are defined as
# objects so that plots can be easily made.

# Define genes of interest
genes_of_interest <- c(
    # Hygroreceptors (humidity sensing)
    "Ir40a",  # Dry tuning receptor
    "Ir68a",  # Moist tuning receptor
    "Ir93a",  # Hygrosensation co-receptor
    
    # CO2 receptors
    "Gr1",    # CO2 receptor subunit
    "Gr2",    # CO2 receptor subunit
    "Gr3",    # CO2 receptor subunit
        
    # VOC receptors
    "Ir8a",   # Acidic volatile detection gating co-receptor
    "Ir25a",  # Amine detection gating co-receptor
    "Ir41a",  # Amine detecting receptor
    "Ir41c",  # Amine detecting receptor
    "Ir75k",  # Carboxylic acid detecting receptor
    "Ir76b",  # Amine detection gating co-receptor
    "Or4",    # Sulcatone detecting receptor
    "Or8",    # (R)-1-octen-3-ol detecting receptor
    "Or11",   # Cis-decahydro-1-naphthol detecting receptor
    "Or49",   # Multiple mammalian/plant volatiles
    "Or71",   # Indole detecting receptor
    "Orco"    # Obligate odorant co-receptor
)

# Define gene subsets
hygroreceptors <- c("Ir93a", "Ir40a", "Ir68a")
co2_receptors <- c("Gr1", "Gr2", "Gr3")
voc_receptors <- c("Orco", "Ir8a", "Ir25a", "Ir76b", "Ir41a", "Ir41c", "Ir75k", 
                   "Or2", "Or4", "Or7", "Or8", "Or11", "Or49", "Or71")

# Save counts of genes of interest to separate object
goi_data <- genes_of_interest[genes_of_interest %in% rownames(count_data)]

# Identify any missing genes of interest 
goi_missing <- genes_of_interest[!genes_of_interest %in% goi_data]