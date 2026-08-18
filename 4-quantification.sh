#!/bin/zsh

# Script function: Quantify gene expression from aligned reads
# Version: 1.0
# Date created: 09-AUG-2026
#
# Dependencies:
#  - featureCounts (quantifying gene expression)
#    https://subread.sourceforge.net/

#===============================================================================
# CONFIGURATION
#===============================================================================

# Retrieve configured environmental variables
source 1-config_setup.sh

#===============================================================================
# QUANTIFICATION
#===============================================================================

# Quantify all samples
if [ -f "${QUANT}/gene_counts.txt" ]; then

	# Skip file
    echo "Already quantified - skipping"

else
	# Run feature counts
    echo "Running featureCounts..."
    featureCounts -T ${THREADS} -t exon -g gene_id -a ${ANNOTATIONS} \
        -o ${QUANT}/gene_counts.txt ${ALIGN}/SRR*_sorted.bam
        
fi

#===============================================================================
# CLEANING UP OF COUNTS FILE
#===============================================================================

# Remove irrelevant columns from counts file not needed for downstream analysis
if [ ! -f "${QUANT}/counts.txt" ]; then
    grep -v "^#" ${QUANT}/gene_counts.txt | cut -f1,7- > ${QUANT}/counts.txt
fi

# Process finish message
echo "Ended on: $(date)"