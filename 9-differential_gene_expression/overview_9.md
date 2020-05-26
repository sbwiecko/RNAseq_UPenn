# Differential gene expression

## Overview

The ultimate goal of most transcriptional profiling experiments is to identify differentially expressed genes or transcripts. In this class, we’ll dig into differential expression using the popular and venerable Limma package in R, while continuing to explore options for producing compelling plots from your differential expression results. Finally, we’ll discuss a workflow for going beyond DGE analysis to look at differentail transcript (isoform) usage (DTU).

## Learning objectives

* Talk about the _model.matrix()_ function
* Talk about how to set your pairwise comparisons using the _contrast.matrix()_ function
* Use the limma package to identify differentially expressed genes
* Produce static and interactive volcano plots
* Use the [isoformSwitchAnalyzeR](https://bioconductor.org/packages/release/bioc/vignettes/IsoformSwitchAnalyzeR/inst/doc/IsoformSwitchAnalyzeR.html) package to carry out an analysis of differential transcript usage (DTU).

## Reading

* [VOOM: precision weights unlock linear model analysis tools for RNA-seq read counts. Genome Biology, Feb, 2014](http://diytranscriptomics.github.io/Reading/files/voom.pdf) - Describes one of the approaches for adjusting RNAseq count data based on the mean-variance relationship
* [Differential analyses for RNA-seq: transcript-level estimates improve gene-level inferences](https://doi.org/10.12688/f1000research.7563.2) – describes the TxImport package with extensive discussion/consideration for DGE vs DTE vs DTU analysis. A must read for this lecture
* [Harold Pimentel’s talk on differential expression with RNAseq](https://www.youtube.com/watch?v=BRWj6re9iGc)
* [how to set-up a design matrix](http://genomicsclass.github.io/book/pages/expressing_design_formula.html)
* [Count-based differential expression analysis of RNA sequencing data using R and Bioconductor. Nature Protocols, Aug 22, 2013](http://diytranscriptomics.github.io/Reading/files/nprot.2013.099.pdf) - This is a great overview of the edgeR and DESeq packages, their use, and explains how each one approaches differential gene expression.
* [Limma user’s guide](http://www.bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf)
* [EdgeR user’s guide](https://www.bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf). See section 3.4 and 3.5 for details about how to modify your model.matrix function for a ‘blocking’ design

## Other videos

* [Josh Starmer from StatQuest describing False Discovery Rates for RNAseq](https://youtu.be/K8LQSvtjcEo)
* [StatQuest for linear regression and least squares](https://youtu.be/PaFPbb66DxQ)
* [More StatQuest for Linear models] (https://youtu.be/nk2CQITm_eo)
