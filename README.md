# mosquitoPE

A RNA-Seq differential expression and functional enrichment pipeline using the _Aedes aegypti_ mosquito transcriptome post-emergence.

## Description

This pipeline downloads raw RNA-Seq reads from the NCBI Sequence Read Archive (SRA), performs QC and trimming on these reads before aligning them to a reference genome and obtaining gene abundances in each sample. Gene counts are then used to identify Differentially Expressed Genes (DEGs), with WGCNA being used to identify genes co-expressed with genes of interest and functional those that share a module.

The project, as it is currently, has been built to identify DEGs within a set of genes of interest (for more information see background folder) in the PRJNA659517 dataset, and then identifies genes sharing a co-expression module with Ir68a, the only gene of interest to have both significant contrasts with the Wald Test and a significant sex x timepoint interaction with the Likelihood-Ratio Test.

Future version of the pipeline will be more configurable/adaptable to different datasets and usage cases.

## Dependencies

### Command-Line tools

[SRA toolkit](https://github.com/esmiano/mosquitoPE/edit/main/README.md): Retrieval of RNA-Seq reads from the NCBI Sequence Read Archive (SRA). SRA files are prefetched and then converted to FASTQ.

[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/): Quality Control, assessing quality of sequencing reads.

[fastp](https://github.com/OpenGene/fastp): Trimming of low quality reads.

[hisat2](https://github.com/esmiano/mosquitoPE/edit/main/README.md): Aligning sequencing reads to reference genome.

[samtools](https://www.htslib.org/): Converting FASTQ files to sorted BAM files

[featureCounts](https://subread.sourceforge.net/featureCounts.html): Quantifying how many reads map to each gene according to GTF annotation file.

### R packages

[vroom](https://vroom.tidyverse.org/) and [tidyverse]() are required for fast importing of data and data handling, respectively. The tidyverse is a suite of packages that provide a number of quality of life improvements primarily utilising pipes and dataframes to make code easier to generate and comprehend.

[ltc](), [RColorBrewer](), and [viridis]() are provide colour palettes, although the exact palettes used is a purely aesthetic choice and shouldn't affect the overall function of the pipeline if edited. [ggrepel]() makes labels easier to read on figures, of particular relevance to the Principal Component Analysis (PCA) plots in this analysis. [cowplot]() is required for plotting figures separately, so isn't integral to pipeline function but is required to produce ready to publish figures.

[ggplot2](), [pheatmap](), and [enrichplot]() are required for plotting.

[DESeq2]() performs the differential gene expression analysis, [WGCNA]() identifies co-expression modules along with [flashClust](https://cran.r-project.org/web/packages/flashClust/index.html), and [gprofiler2]() and [clusterProfiler]() functionally annotate co-expressed genes with GO enrichment.

[PoiClaClu]() is required for Poisson similarity matrices (DESeq2 results QC) and [ashr]() was the algorithm used for shrinkage of those results (used when ranking genes as part of the Gene Set Enrichment Analysis).

[AnnotationHub](), [AnnotationDbi](), and [GO.db]() are required for GO term database handling and standardisation.

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
install.packages("tidyverse", "ggplot2", "pheatmap", "ashr", "flashClust", "gprofiler2")
```

For R packages from the Bioconductor ecosystem, first install Bioconductor:

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
