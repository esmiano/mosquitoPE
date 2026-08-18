#!/bin/zsh

# Script function: Align reads by mapping to reference genome
# Version: 1.0
# Date created: 09-AUG-2026
#
# Dependencies:
#  - HISAT2 (building HISAT2 index for alignment to reference genome)
#    https://github.com/DaehwanKimLab/hisat2
#  - samtools (converting SAM files to sorted BAM files)
#    https://github.com/samtools/samtools

#===============================================================================
# CONFIGURATION
#===============================================================================

# Retrieve configured environmental variables
source 1-config_setup.sh

#===============================================================================
# ALIGNMENT
#===============================================================================

# Loop over trimmed fastq files
for trm in ${TRIM}/SRR*.fastq.gz; do

    # Extract SRR accession
    srr=$(basename "$trm" _trimmed.fastq.gz)
    echo "$(date) - Starting ${srr}..."

    # Check if read already mapped
    if [ -f "${ALIGN}/${srr}.sam" ] \
        || [ -f "${ALIGN}/${srr}_sorted.bam" ]; then

        # Skip file
		echo "${srr} already mapped - Skipping..."

    else
		# Map reads
        hisat2 -x ${HISAT2_INDEX} -U "$trm" \
            --known-splicesite-infile ${SPLICE_SITES} \
            -p ${THREADS} -S ${ALIGN}/${srr}.sam \
            --summary-file ${ALIGN}/${srr}_summary.txt
        
        # Convert to sorted BAM file
        samtools sort ${ALIGN}/${srr}.sam \
         -o ${ALIGN}/${srr}_sorted.bam

        # Index used for some downstream tools
        # TODO: Uncomment if needed
		# samtools index ${ALIGN}/${srr}_sorted.bam

        # Remove original SAM file to save space
        rm ${ALIGN}/${srr}.sam

    fi

    echo "$(date) - Finished mapping ${srr}."

done

# Process finish message
echo "Ended on: $(date)"