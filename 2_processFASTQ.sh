#!/bin/zsh
#===============================================================================
# RAW DATA PRE-PROCESSING AND QUALITY CONTROL
#===============================================================================
# Author: esmiano
# Version: 1.0
# Date created: 09-AUG-2026
# 
# Description:
# This script pre-processes and QCs raw RNA-Seq data, ready for alignment and
# quantification.

#===============================================================================
# CONFIGURATION
#===============================================================================

# Retrieve configured environmental variables
source 1-config_setup.sh

#===============================================================================
# DATA PRE-PROCESSING AND QC
#===============================================================================

# Loop over each downloaded SRA file
for sra in ${SRA}/PRJNA659517/SRR*.sra; do

    # Extract SRR accession from the filename
    srr=$(basename "$sra" .sra)	
	
	#===========================================================================
	# SRA TO COMPRESSED FASTQ CONVERSION
	#===========================================================================
	
	# Check if SRA file already converted
	if [ -f "${FASTQ}/${srr}.fastq.gz" ]; then
		
		# Skip file
		echo "${srr} already converted - Skipping..."
		
	else
		# Convert to FASTQ and output to new directory
		echo "Converting to FASTQ..."
		fasterq-dump "$sra" --outdir ${FASTQ}
		
		# Compress output FASTQ files
		echo "Compressing..."
		gzip ${FASTQ}/${srr}.fastq
		
	fi
	
	#===========================================================================
	# PRE-TRIM FASTQC
	#===========================================================================
	
	# Check if compressed FASTQ file already checked
	if [ -f "${FASTQC_PRE}/${srr}_fastqc.html" ]; then
		
		# Skip file
		echo "${srr} already checked - Skipping..."
		
	else
		# Run FastQC on compressed file and output to new directory
		echo "Running FastQC..."
		fastqc -o ${FASTQC_PRE} ${FASTQ}/${srr}.fastq.gz
		
	fi
	
	#===========================================================================
	# TRIMMING OF COMPRESSED FASTQ
	#===========================================================================
	
	# Check if compressed FASTQ file already trimmed
	if [ -f "${TRIM}/${srr}_trimmed.fastq.gz" ]; then
		
		# Skip file
		echo "${srr} already trimmed - Skipping..."
		
	else
		# Trim compressed files with Fastp and output to new directory
		echo "Running Fastp..."
		fastp -i ${FASTQ}/${srr}.fastq.gz -o ${TRIM}/${srr}_trimmed.fastq.gz \
    	-q 30 -l 40 -j ${TRIM}/${srr}_fastp.json -h ${TRIM}/${srr}_fastp.html
		
	fi
	
	#===========================================================================
	# POST-TRIM FASTQC
	#===========================================================================
	
	# Check if trimmed FASTQ file already checked
	if [ -f "${FASTQC_POST}/${srr}_trimmed_fastqc.html" ]; then
		
		# Skip file
		echo "${srr} already checked - Skipping..."
		
	else
		# Run FastQC on trimmed file and output to new directory
		echo "Running FastQC..."
		fastqc -o ${FASTQC_POST} ${TRIM}/${srr}_trimmed.fastq.gz
		
	fi
	
	echo "$(date) - Finished pre-processing ${srr}."
	
done

# Process finish message
echo "Ended on: $(date)"