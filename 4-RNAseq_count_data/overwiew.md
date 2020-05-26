# Overview

Now that we've aligned our reads, it time to discuss units for measuring gene expression. We'll discuss differences between RPKM and TPM, and how these units relate to basic properties of your reference file and data. We'll also discuss normalization within and between samples. To conclude this class, we'll fire up RStudio and take a look at our first script.

## Learning objectives

* Review steps from last class (using Kallisto).
* Discuss output from Kallisto and units of measurement for RNAseq and 'normalization'
* Start an RStudio Project directory that we'll use for the rest of the course.
* Open and discuss our first script, including installation of packages

## Reading

* [The RNA-seq abundance zoo](http://robpatro.com/blog/?p=235) - Blog post by Rob Patro (developer of Salfish and Salmon software) that describes units for RNAseq, and has a nice description of 'effective length' for transcripts.

* [What the FPKM?](https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/) - Blog post by Harold Pimentel discussing within sample normalization and the meaning of RNAseq expression units.

* [Measurement of mRNA abundance using RNA-seq data: RPKM measure is inconsistent among samples. Theory in Biosciences, Dec 2012](http://diytranscriptomics.github.io/Reading/files/wagnerTPM.pdf)

* [Between sample normalization in RNAseq](https://haroldpimentel.wordpress.com/2014/12/08/in-rna-seq-2-2-between-sample-normalization/) - another great blog post from Harold Pimentel on between-sample normalization.
