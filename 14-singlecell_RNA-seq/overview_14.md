# Single cell RNA-seq -- principles and processing

## Overview

Now that you're comfortable with bulk RNA-seq data analysis, we'll shift our focus to the rapidly developing landscape of single cell RNA-seq (scRNA-seq). In this lecture, you'll learn about the underlying technology and we'll demonstrate how to process raw single cell data directly on your laptop (!) for importing into R/bioconductor.

## Learning objectives

* Understand droplet-based scRNA-seq technology
* Be able to compare and contrast single cell and bulk RNA-seq methods
* Understand cost and experimental design considerations for scRNA-seq experiments.
* Familiarity with multiplexed single cell assays (CITE-seq, 'multiome', TEA-seq)
* Be able to define common terms and concepts in single cell genomics
* Use Kallisto-BUStools to preprocess raw scRNA-seq data (via [kb-python](https://www.kallistobus.tools/))

## What you need to do

[Download raw files](https://drive.google.com/drive/folders/1DbLRO4kv-y3W06adFR26RdSaDPmfB4UA?usp=sharing). You will need about 5Gb of storage space on your harddrive to accomodate this download. Please do not uncompress these files (leave them as .gz files). This is data from 1000 peripheral blood mononuclear cells (PBMCs) and is one of the sample datasets provided by 10X Genomics [here](https://bit.ly/10xPBMC_small). I merged the separate lane files to make this simpler to work with for the course.

[Human transcriptome reference index file](https://diytranscriptomics.com/project/lecture-02) - this is the index you created using Kallisto way back in lecture 2. If you don't have this, remember it's easy to create using `kallisto` index.

[t2g.txt](http://diytranscriptomics.github.io/Code/files/t2g.txt) - this is a human transcript-to-gene mapping file that we will use with Kallisto-Bustools to preprocess our data. This file is easy to generate with kb ref, but downloading it now will save you some time.

### How to prepare the reference files scRNA-seq

```bash
kb ref \
-d human \ #download a prebuilt index
-i Homo_sapiens.GRCh38.cdna.all.index \
-g t2g.txt
```

[kb-python](https://github.com/pachterlab/kb_python) - You will need to have this software installed in a Conda environment on your laptop. We did this way back in [lecture 1](https://diytranscriptomics.com/project/lecture-01). If you are unable to install or use kb-python, just follow along with the lecture so you understand the concepts.

### Note for installation

kb-python will only work on WSL, see the [explanations on stackoverflow](https://stackoverflow.com/questions/60197890/why-does-installing-pysam-python-package-fail). Install the package using `pip3 install kb_python`.

### Preprocess our raw fastq files

```bash
kb count \
pbmc_1k_v3_S1_mergedLanes_R1.fastq.gz pbmc_1k_v3_S1_mergedLanes_R2.fastq.gz \
-i Homo_sapiens.GRCh38.cdna.all.index \
-x 10XV3 \ # technology used
-g t2g.txt \
-t 8 \
--cellranger # convert counts to cellranger-compatible format
```

## Reading

[Modular and efficient pre-processing of single-cell RNA-seq](https://doi.org/10.1101/673285) - describes the full Kallisto-Bustools workflow for memory efficient processing of scRNA-seq data.

[The barcode, UMI, set format and BUStools, Bioinformatics](https://doi.org/10.1093/bioinformatics/btz279) - Describes the BUS format as an efficient and platform-independent way to store information from scRNA-seq data.

[A curated database reveals trends in single-cell transcriptomics](https://doi.org/10.1093/database/baaa073) - describes the growing collection of scRNA-seq experiments found [here](https://diytranscriptomics.com/project/www.nxn.se/single-cell-studies/gui), which I used to produce two of the plots in the slides for this lecture.
