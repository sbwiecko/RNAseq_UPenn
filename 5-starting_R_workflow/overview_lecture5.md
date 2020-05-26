# Overview

We’ll begin this class by reviewing how to access R packages and help documentation, as well as understanding the basic structure of an R script and RStudio project. We’ll then access annotation data before reading our Kallisto results into R.

## Learning objectives

- Review basic elements of an R script
- Learn how to access R packages, and their help documentation
- Understand the importance of a study design file
- Access annotation information for transcripts using Bioconductor
- Read Kallisto transcript abundance measurements into the R environment using TxImport

## Reading

- [Differential analysis of RNA-seq incorporating quantification uncertainty](http://diytranscriptomics.github.io/Reading/files/sleuth.pdf). Nature Methods, June, 2017 - Original paper describing Sleuth
- [Lior Pachter’s blog post on Sleuth](https://liorpachter.wordpress.com/2015/08/17/a-sleuth-for-rna-seq/)
- [vignette for the Tximport package](https://bioconductor.org/packages/devel/bioc/vignettes/tximport/inst/doc/tximport.html) - the R package we’ll use to read the Kallisto mapping results into R.
- [Differential analyses for RNA-seq: transcript-level estimates improve gene-level inferences*](http://f1000research.com/articles/4-1521/v2) F1000Research, Dec 2015. This paper describes the Tximport package and its application for handling transcript-level expression measurments from lightweight aligners (Salmon, Sailfish, Kallisto)
