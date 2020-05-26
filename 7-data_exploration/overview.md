# Overview

In this class you’ll learn about a variety of approaches exploring your data. You’ll use multivariate statistical approaches such as Principal Component Analysis (PCA) to understand sources of variance in our data, while continuing to build your plotting skills by using ggplot2 to graph the results of PCA analyses. You’ll also learn how to use the dplyr package to take control over our gene expression dataframes, allowing us to change, sort, filter, arrange and summarize large data sets quickly and easily using simple commands in R. We’ll discuss common missteps and how to identify sources of bias in transcriptional data sets.

## Learning objectives

* Start and finish Step 3 script
* Discuss basics of multivariate statistical analysis
* Carry out hierarchical clustering of samples
* Discuss and perform principal component analyses (PCA)
* Produce ‘small multiples’ plot
* Use standard dplyr ‘verbs’ to quickly query our data
* Produce interactive graphics using the plotly package
* Produce interactive tables with DT package

## Reading

* [Ten quick tips for effective dimensionality reduction](https://doi.org/10.1371/journal.pcbi.1006907) - a absolute must-read for understanding data exploration methods.
* [Blog post describing T-SNE](http://distill.pub/2016/misread-tsne/) - I mentioned various unsupervised linear methods for dimensional reduction of your data (PCA, MDS). T-SNE and UMAP are non-linear unsupervised methods that have become popular for representing single-cell RNAseq data and flow cytometry data.
* [Original T-SNE paper](http://diytranscriptomics.github.io/Reading/files/TSNE.pdf).
* [UMAP paper](https://www.nature.com/articles/nbt.4314) - A new algorithm, called uniform manifold approximation and projection (UMAP) has been recently published and is gaining popularity in single cell RNAseq and flow cytometry analysis. UMAP is proposed to preserve as much of the local and more of the global data structure than t-SNE, with a shorter run time.