#!/bin/zsh

# Script function: Configure environmental variables for directories and
#   reference genome files, ready for downstream analysis
# Version: 1.0
# Date created: 09-AUG-2026
#
# Dependencies:
#  - sra-tools (downloading RNA-Seq data from the Sequence Read Archive (SRA))
#    https://github.com/ncbi/sra-tools
#  - wget (downloading reference genome, annotations, and sample metadata)
#    https://github.com/rockdaboot/wget2
#  - HISAT2 (building HISAT2 index for alignment to reference genome)
#    https://github.com/DaehwanKimLab/hisat2

#===============================================================================
# CONFIGURATION
#===============================================================================

# User must configure the path to the working directory, all other variables
#   are configured automatically.

# Configure working directory (requires user input)
export WORK=/path/to/working/directory

# Configure sub-directories
export SRA=${WORK}/1-sra
export FASTQ=${WORK}/2-fastq
export FASTQC_PRE=${WORK}/3-fastqc/3.1-pre_trim
export FASTQC_POST=${WORK}/3-fastqc/3.2-post_trim
export TRIM=${WORK}/4-trim
export ALIGN=${WORK}/5-align
export QUANT=${WORK}/6-quant
export OUTPUT=${WORK}/7-output
export REFERENCE=${WORK}/reference

# Configure reference genome files
export GENOME=${REFERENCE}/AaegL5_reference_genome.fa
export ANNOTATIONS=${REFERENCE}/Adavi2024_annotations.gtf
export SPLICE_SITES=${REFERENCE}/splice_sites.txt
export HISAT2_INDEX=${REFERENCE}/index_hisat2/Aedes_aegypti_lvpagwg.AaegL5.dna.toplevel.fa

# Configure threads 
export THREADS=8

#===============================================================================
# SETUP
#===============================================================================

# Create sub-directories
mkdir -p ${SRA} ${FASTQ} ${FASTQC_PRE} ${FASTQC_POST} ${TRIM} \
    ${ALIGN} ${QUANT} ${OUTPUT} ${REFERENCE}/index_hisat2

# Transcriptomic dataset -------------------------------------------------------

# Pre-fetch (download) RNA-Seq data from SRA
if [ ! -d "${SRA}/PRJNA659517" ]; then
    echo "Downloading transcriptomic dataset..."
    prefetch PRJNA659517 --output-directory ${SRA}
fi

# Reference genome and annotations ---------------------------------------------

# Download reference genome if not present
if [ ! -f "${GENOME}" ]; then
    echo "Downloading reference genome..."
    wget -q -O ${GENOME}.gz \
        https://ftp.ensemblgenomes.ebi.ac.uk/pub/metazoa/release-63/fasta/aedes_aegypti_lvpagwg/dna/Aedes_aegypti_lvpagwg.AaegL5.dna.toplevel.fa.gz
    gunzip ${GENOME}.gz
fi
 
# Download updated annotation file if not present
if [ ! -f "${ANNOTATIONS}" ]; then
    echo "Downloading updated annotation file..."
    wget -q -O ${ANNOTATIONS} \
        "https://zenodo.org/records/12801833/files/AaegyptiLVP_AGWG_ThreePrimeUTRextended_Adavi2024.gtf?download=1"
fi

# Strip AaegL5_ prefix from chromosome names in annotations if still present
if grep -q '^AaegL5_' ${ANNOTATIONS}; then
    echo "Stripping AaegL5_ prefix from GTF chromosome names..."
    # Zsh version of sed command is used below
    sed -i '' 's/^AaegL5_//' ${ANNOTATIONS}
    # For bash, please use the following
    # sed -i 's/^AaegL5_//' ${ANNOTATIONS}
fi
 
# Extract splice sites for splice-aware alignment
if [ ! -f "${SPLICE_SITES}" ]; then
    echo "Extracting splice sites..."
    hisat2_extract_splice_sites.py ${ANNOTATIONS} > ${SPLICE_SITES}
fi
 
# Build HISAT2 index if not present
if [ ! -f "${HISAT2_INDEX}.1.ht2" ]; then
    echo "Building HISAT2 index..."
    hisat2-build -p ${THREADS} ${GENOME} ${HISAT2_INDEX}
fi

# Sample metadata --------------------------------------------------------------

# Download dataset report containing sample metadata if not present
if [ ! -f "${REFERENCE}/PRJNA659517.tsv" ]; then
    echo "Downloading sample metadata..."
    wget -q -O "${REFERENCE}/PRJNA659517.tsv" \
        "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA659517&result=read_run&fields=run_accession,library_name&format=tsv&download=true&limit=0"
fi

# Process finish message
echo "Ended on: $(date)"