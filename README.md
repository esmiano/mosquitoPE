# mosquitoPE

A RNA-Seq differential expression and functional enrichment pipeline using the _Aedes aegypti_ mosquito transcriptome post-emergence.

## Description

This pipeline downloads raw RNA-Seq reads from the NCBI Sequence Read Archive (SRA), performs QC and trimming on these reads before aligning them to a reference genome and obtaining gene abundances in each sample. Gene counts are then used to identify Differentially Expressed Genes (DEGs), with WGCNA being used to identify genes co-expressed with genes of interest and functional those that share a module.

The project, as it is currently, has been built to identify DEGs within a set of genes of interest (for more information see background folder) in the PRJNA659517 dataset, and then identifies genes sharing a co-expression module with Ir68a, the only gene of interest to have both significant contrasts with the Wald Test and a significant sex x timepoint interaction with the Likelihood-Ratio Test.

Future version of the pipeline will be more configurable/adaptable to different datasets and usage cases.


## Dependencies

More details to come here.

### Command-Line tools

SRA tools

FastQC

fastp

hisat2

samtools

featureCounts


### R packages

tidyverse

DESeq2

WGCNA

gprofiler2

AnnotationHub

clusterProfiler

ashr

ggplot2

pheatmap

enrichplot

## Installation

For the SRA toolkit, follow the instructions on the project [GitHub](https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit) to install.

FastQC can be downloaded from the Babraham Bioinformatics Projects [webpage](https://www.bioinformatics.babraham.ac.uk/projects/download.html#fastqc).

fastp can be installed with bioconda as below or compiled from source, for which the instructions can be found on the project [GitHub](https://github.com/OpenGene/fastp#install-with-bioconda).

```bash
conda install -c bioconda fastp
```

For HISAT2, download the binary from the project [website](https://daehwankimlab.github.io/hisat2/download/). samtools can also be downloaded directly from the project [website](https://www.htslib.org/download/).

featureCounts is part of the subread package, which can be installed from SourceForge following the instructions [here](https://subread.sourceforge.net/subread-package.html).

For R packages available through CRAN, please install them using the following command:

```R
install.packages("tidyverse", "gprofiler2", "ashr", "ggplot2", "pheatmap")
```

For packages from the Bioconductor ecosystem, first install Bioconductor:

```R
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
```

Once Bioconductor is installed, run the following command:

```R
BiocManager::install("DESeq2", "WGCNA", "AnnotationHub", "clusterProfiler", "enrichplot")
```

### Pipeline

Clone the repository in the location where you wish to perform the analysis. This folder will be the working directory. 

## Usage

To run the pipeline, run each script in ascending order from 1 to 9. Ensure that all dependencies are installed before running.

Running the pipeline will download data and reference files from the NCBI. This takes a long time and requires a significant amount of storage so it is recommended that the initial 4 scripts are run on a High-Performance Computing (HPC) system.

All datasets, referece files, directories, and output files will be created within the working directory as each file is run.

## Contact

If you’d like to get in contact with any comments or issues when running the pipeline, please email me at spensleymj@gmail.com.