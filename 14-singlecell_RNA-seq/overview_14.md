# Single cell RNA-seq -- principles and processing

## Overview

Now that you're comfortable with bulk RNA-seq data analysis, we'll shift our focus to the rapidly developing landscape of single cell RNA-seq (scRNA-seq). In this lecture, you'll learn about the underlying technology and we'll demonstrate how to process raw single cell data directly on your laptop (!) for importing into R/bioconductor.

With your data already preprocessed with Kallisto-Bustools, you're now ready to import into R and use a variety of packages to filter, plot and analyze your data.

## Learning objectives

* Understand droplet-based scRNA-seq technology
* Be able to compare and contrast single cell and bulk RNA-seq methods
* Understand cost and experimental design considerations for scRNA-seq experiments.
* Familiarity with multiplexed single cell assays (CITE-seq, 'multiome', TEA-seq)
* Be able to define common terms and concepts in single cell genomics
* Use Kallisto-BUStools to preprocess raw scRNA-seq data (via [kb-python](https://www.kallistobus.tools/))

In the second part of the course:

* Be able to import preprocessed data into R and create a Seurat object
* Carry out filtering using [DropUtils](https://bioconductor.org/packages/release/bioc/html/DropletUtils.html)
* Generate a Quality Control report of your scRNA-seq data directly within R
* Use standard QC metrics and plots to filter your data
* Generate clusters and visualize via UMAP dimensional reduction
* Find cluster-specific marker genes with Seurat
* Annotate unknown clusters using public databases and [CellAssign](https://www.rdocumentation.org/packages/cellassign/) and [SingleR](https://bioconductor.org/packages/release/bioc/html/SingleR.html)
* Integrate multiple samples and use sample details to analyze integrated data, e.g., compared different conditions in different samples

## What you need to do

[Download raw files](https://drive.google.com/drive/folders/1DbLRO4kv-y3W06adFR26RdSaDPmfB4UA?usp=sharing). You will need about 5Gb of storage space on your harddrive to accomodate this download. Please do not uncompress these files (leave them as .gz files). This is data from 1000 peripheral blood mononuclear cells (PBMCs) and is one of the sample datasets provided by 10X Genomics [here](https://bit.ly/10xPBMC_small). I merged the separate lane files to make this simpler to work with for the course.
Of note, the _sample index read_ is not used here as there is a single type of sample, see [more information](https://kb.10xgenomics.com/hc/en-us/articles/115004553003-How-do-I-demultiplex-a-run-with-no-index-sequences-).

Download the [human genome data from Ensembl](http://ftp.ensembl.org/pub/release-108/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz), and the [corresponding annotations](http://ftp.ensembl.org/pub/release-108/gtf/homo_sapiens/Homo_sapiens.GRCh38.108.gtf.gz). The latest genome and annotation files for every species on ensembl can be found with the gget command-line tool `gget ref -w dna,gtf homo_sapiens`.
Another option is to use `kb ref -d human` to download a prebuild index file, as described later.

[t2g.txt](http://diytranscriptomics.github.io/Code/files/t2g.txt) - this is a human transcript-to-gene mapping file that we will use with Kallisto-Bustools to preprocess our data. This file is easy to generate with kb ref (see later), but downloading it now will save you some time.

## Code and files

* [pre-processed data for 1000 PBMCs](https://drive.google.com/drive/folders/1RO45z5DEVpuaq5qwlF5QNdhc0tbGVK7l?usp=sharing) - You only need to download this if you were unable to use kb-python in the last lecture to process raw scRNA-seq data. This ensures that everyone can follow along with this lecture, regardless of whether you were able to install or use Kb-python.
* [DIY_scRNAseq.R](http://diytranscriptomics.github.io/Code/files/DIY_scRNAseq.R) - this is the R script that we'll use for this lecture.
* [functions.R](http://diytranscriptomics.github.io/Code/files/functions.R) - this is the custom R function we'll use for generating a QC report with our scRNA-seq data (see Reading material below for source).
* [Seurat objects](https://drive.google.com/drive/folders/1SEEr70W6D9itvVaLfEXuRVWqyj_Gribc?usp=sharing) - this folder contains two Seurat objects from an unpublished mouse experiment (courtesy of Chris Hunter's lab). One sample is from a naive control mouse, while the second is from a mouse infected with the protozoan parasite, Toxoplasma gondii (14 days post-infection). We'll use these data in the second 1/2 of the lecture to practice integration and differential gene testing between conditions.

[kb-python](https://github.com/pachterlab/kb_python) - You will need to have this software installed in a Conda environment on your laptop. We did this way back in [lecture 1](https://diytranscriptomics.com/project/lecture-01). If you are unable to install or use kb-python, just follow along with the lecture so you understand the concepts.

### Note for installation

kb-python will only work on WSL, see the [explanations on stackoverflow](https://stackoverflow.com/questions/60197890/why-does-installing-pysam-python-package-fail). Install the package using `sudo pip3 install kb_python`.

```bash
# prior to the installation of kb_python, install PySAM
sudo pip3 install pysam
# there may be some error messages pointing to missing packages, e.g.,
sudo apt install liblzma-dev
sudo apt install libbz2-dev

# finally check the list of supported technologies
kb --list
```

### How to prepare the reference files scRNA-seq

#### Without genome and annotations files

```bash
kb ref \
-d human \ #download a prebuilt index
-i Homo_sapiens.GRCh38.cdna.all.index \
-g t2g.txt
```

#### Build locally from genome and annotation files [as described elsewhere](https://github.com/pachterlab/kb_python)

```bash
kb ref \
-i Homo_sapiens.GRCh38.cdna.all.index \
-g t2g.txt \
-f1 transcriptome.fa \
Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz \
Homo_sapiens.GRCh38.108.gtf.gz
```

The _de novo_ building of the index can take a couple of minutes.
Note -- See the [issue on generating t2g.txt and index](https://github.com/pachterlab/kb_python/issues/180)

### Preprocess our raw fastq files

```bash
kb count \
pbmc_1k_v3_S1_mergedLanes_R1.fastq.gz pbmc_1k_v3_S1_mergedLanes_R2.fastq.gz \
-i Homo_sapiens.GRCh38.cdna.all.index \
-x 10XV3 \ # technology used
-g t2g.txt \
--overwrite \ #overwrite existing output.bus file
#-t 8 \ # default is 8 threads
--cellranger # convert counts to cellranger-compatible format
```

In addition to unfiltered counts, different diagnostic files are bieng generated and will be used later for the QC.

### Note for analyzing data in VSCode

The graphs and figures are more complex, so it's better to use an optimized backend/grpahics server such as [`httpgd`](https://cran.rstudio.com/web/packages/httpgd/index.html). Once installed, select the following options in the VSCode preferences:

- R > Plot > Defaults: Full Window Mode: set on
- R > Plot > Use Httpgd: set on

## Reading

* [Modular and efficient pre-processing of single-cell RNA-seq](https://doi.org/10.1101/673285) - describes the full Kallisto-Bustools workflow for memory efficient processing of scRNA-seq data.
* [The barcode, UMI, set format and BUStools, Bioinformatics](https://doi.org/10.1093/bioinformatics/btz279) - Describes the BUS format as an efficient and platform-independent way to store information from scRNA-seq data.
 *[A curated database reveals trends in single-cell transcriptomics](https://doi.org/10.1093/database/baaa073) - describes the growing collection of scRNA-seq experiments found [here](https://diytranscriptomics.com/project/www.nxn.se/single-cell-studies/gui), which I used to produce two of the plots in the slides for this lecture.
* [EmptyDrops: distinguishing cells from empty droplets in droplet-based single-cell RNA sequencing data..](https://doi.org/10.1186/s13059-019-1662-y) - This is the paper describing the DropletUtils package that we use in this lecture to identify empty drops.
* [Sarah Ennis' Github repo for preprocessing scRNA-seq data](https://github.com/Sarah145/scRNA_pre_process) - This is the source of the custom script we use to generate the CellRanger-esque html QC report.
* [Probabilistic cell-type assignment of single-cell RNA-seq for tumor microenvironment profiling](https://pubmed.ncbi.nlm.nih.gov/31501550/) - Describes the CellAssign algorithm and R package that we use to identify clusters.
* [Reference-based analysis of lung single-cell sequencing reveals a transitional profibrotic macrophage](https://doi.org/10.1038/s41590-018-0276-y) - describes the SingleR and celldex packages that allow us to leverage bulk RNA-seq data in public repositories to curate clusters in our scRNA-seq.
* [Comprehensive Integration of Single-Cell Data](https://doi.org/10.1016/j.cell.2019.05.031) - This 2019 paper describes the underlying statistical approach for data integration in Seurat.
