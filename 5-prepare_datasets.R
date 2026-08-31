#===============================================================================
# DATASET PREPARATION
#===============================================================================
# Author: esmiano
# Version: 1.1
# Date created: 09-AUG-2026
# Date modified: 
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
library(vroom)
library(tidyverse)

# Create separate output directories for results and figures
#   showWarnings set to FALSE to prevent warnings in case of pre-existing directory
#   recursive set to TRUE to create parent directories, if necessary
dir.create(file.path(results_dir), showWarnings = FALSE, recursive = TRUE)

dir.create(file.path(figures_dir, "ma_plots"),
                     showWarnings = FALSE, recursive = TRUE)

dir.create(file.path(figures_dir, "expression_plots"),
                     showWarnings = FALSE, recursive = TRUE)

#===============================================================================
# COUNT DATA PREPARATION
#===============================================================================

# DESeq2 takes only numerical data as input so the output from featureCounts
# needs to be prepared accordingly.

# Import featureCounts output
count_data <- read.table(
    file.path(working_dir, "6-quant", "counts.txt"),
        # header set to TRUE to indicate first line is a header
        header = TRUE,
        # Columns tab separated
        sep = "\t",
        # Replace quotes with the escaped equivalent
        quote = ""
    )

# Remove AAEL accession IDs when gene name present in Geneid column
# Iterate over each row in Geneid column
# TODO: Finish this
#for(gene in 1:nrow(count_data[, 1])) {
  #if (str_extract(count_data[, 1], "AAEL\\d+") == F) {
    #test <- data.frame(
        #geneid = c(if)
    #)
  #}, 
#}

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
metadata <- vroom(
    file.path(working_dir, "reference", "PRJNA659517.tsv"))

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

# PCA analysis indicated that 2 pairs of samples may have been mislablled and
# effectively swapped. The following were excluded from the final analysis:
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
# GENE ANNOTATION
#===============================================================================



#===============================================================================
# GENES OF INTEREST
#===============================================================================

# Genes involved in the detection of humidity, CO2, and various Volatile Organic
# Compounds (VOCs) were identified from the literature. These are defined as
# objects so that plots can be easily made.

# Define genes of interest
genes_of_interest <- c(
    # Hygroreceptors (humidity sensing)
    "Ir40a",  # Expressed in dry cells
    "Ir68a",  # Expressed in moist cells
    "Ir93a",  # Co-receptor, expressed in both moist and dry cells
    
    # CO2 receptors
    "Gr1",    # CO2 receptor subunit
    "Gr2",    # CO2 receptor subunit
    "Gr3",    # CO2 receptor subunit
    
    # Volatile/olfactory co-receptors
    "Orco",   # Obligate odorant co-receptor (DeGennaro et al. 2013; Larsson 2004)
    "Ir8a",   # Detects acidic volatiles incl. lactic acid (Raji et al. 2019)
    "Ir25a",  # Highly conserved co-receptor (Abuin et al. 2011)
    "Ir76b",  # Co-receptor (Ye et al. 2022; Goldman et al. 2025)
    
    # VOC receptors
    "Ir21a",  # Hill et al. 2021
    "Ir41a",  # Herre et al. 2022
    "Ir41c",  # Raji et al. 2023
    "Ir41j",  # Hill et al. 2021
    "Ir75d",  # Hill et al. 2021
    "Ir75g",  # Herre et al. 2022
    "Ir75l",  # Raji et al. 2023; Hill et al. 2021
    "Ir100a", # Herre et al. 2022
    "Ir101",  # Hill et al. 2021
    "Ir161",  # Herre et al. 2022
    "Or2",    # Goldman et al. 2025; Xiong et al. 2025
    "Or4",    # Tallon et al. 2019
    "Or8",    # Recognises (R)-1-octen-3-ol (Bohbot & Dickens 2009)
    "Or10",   # Xiong et al. 2025
    "Or11",   # Hill et al. 2021; Tallon et al. 2019; Xiong et al. 2025
    "Or23",   # Tallon et al. 2019
    "Or28",   # Tallon et al. 2019
    "Or41",   # Tallon et al. 2019
    "Or47",   # Tallon et al. 2019
    "Or49",   # Herre et al. 2022
    "Or52",   # Hill et al. 2021; Tallon et al. 2019
    "Or57",   # Hill et al. 2021
    "Or59",   # Tallon et al. 2019
    "Or66",   # Tallon et al. 2019
    "Or69",   # Tallon et al. 2019
    "Or70",   # Tallon et al. 2019
    "Or71",   # Herre et al. 2022; Tallon et al. 2019
    "Or72",   # Hill et al. 2021
    "Or81",   # Hill et al. 2021
    "Or82",   # Goldman et al. 2025
    "Or84",   # Tallon et al. 2019
    "Or87",   # Tallon et al. 2019
    "Or88",   # Tallon et al. 2019
    "Or91",   # Tallon et al. 2019
    "Or94",   # Tallon et al. 2019
    "Or100",  # Tallon et al. 2019
    "Or103",  # Tallon et al. 2019
    "Or104",  # Tallon et al. 2019
    "Or105",  # Hill et al. 2021; Tallon et al. 2019
    "Or112",  # Tallon et al. 2019
    "Or113"   # Tallon et al. 2019
)

# Define gene subsets
hygroreceptors <- c("Ir93a", "Ir40a", "Ir68a")
co2_receptors <- c("Gr1", "Gr2", "Gr3")
coreceptors <- c("Orco", "Ir8a", "Ir25a", "Ir76b")
voc_receptors <- c("Ir21a", "Ir41a", "Ir41c", "Ir41j", "Ir75d", "Ir75g", "Ir75l", 
    "Ir100a", "Ir101", "Ir161", "Or2", "Or4", "Or8", "Or10", "Or11", "Or23", 
    "Or28", "Or41", "Or47", "Or49", "Or52", "Or57", "Or59", "Or66", "Or69", 
    "Or70", "Or71", "Or72", "Or81", "Or82", "Or84", "Or87", "Or88", "Or91", 
    "Or94", "Or100", "Or103", "Or104", "Or105", "Or112", "Or113"
)

# Save counts of genes of interest to separate object
goi_data <- genes_of_interest[genes_of_interest %in% rownames(count_data)]

# Identify any missing genes of interest 
goi_missing <- genes_of_interest[!genes_of_interest %in% goi_data]
