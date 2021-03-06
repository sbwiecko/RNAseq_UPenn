# Functional Enrichment Analysis

Now that you’ve identified differentially expressed genes, what do they mean and how do you begin to elucidate the biological pathways governed by these genes? Toward this end, you will learn how to carry out functional enichment analyses using [Gene Ontology](http://geneontology.org/) and [Gene Set Enrichment Analysis (GSEA)](http://software.broadinstitute.org/gsea/index.jsp). We’ll also explore different options for how to present your functional enrichment results graphically.

## Learning objectives

* Carry out Gene Ontology (GO) enrichment analysis using modules identified in the previous script
* Carry out a GSEA enrichment analysis using our full dataset
* Understand the differences between GO and GSEA
* Understand the MSigDB resource and how to access signature collections
* Discuss the basic building blocks for assembling a figure

## Reading

* [The what, where, how and why of gene ontology – a primer for bioinformaticians](http://diytranscriptomics.github.io/Reading/files/GO.pdf). Briefings in Bioinformatics, Feb 2011
* [A nice blog post on the hypergeometric test and Fisher’s exact test](http://diytranscriptomics.com/project/httP;//mengnote.blogspot.com/2012/12/calculate-correct-hypergeometric-p.html) - these statistical tests are at the core of many functional enrichment approaches
* [Analyzing gene expression data in terms of gene sets: methodological issues](https://doi.org/10.1093/bioinformatics/btm051) - A seminal paper on the statistics of enrichment analysis in gene expression studies
* [Toward a gold standard for benchmarking gene set enrichment analysis](https://doi.org/10.1093/bib/bbz158) - A excellent and recent benchmarking study for enrichment tools
* [original 2003 Nat. Methods paper describing Gene Set Enrichment Analysis (GSEA)](http://diytranscriptomics.github.io/Reading/files/Mootha2003_GSEA.pdf), and the [2005 PNAS paper](http://mootha.med.harvard.edu/PubPDFs/Subramanian2005.pdf) that formally detailed its usage.
* You can carry out self-contained and competitive GSEA in R using [ROAST](http://diytranscriptomics.github.io/Reading/files/ROAST.pdf) and [CAMERA](http://diytranscriptomics.github.io/Reading/files/CAMERA.pdf), respectively
* [Gene Set VARIATION Analysis (GSVA)](http://diytranscriptomics.github.io/Reading/files/GSVA.pdf) - I find GSVA useful for producing GSEA-type results across a heterogeneous dataset (like a large cohort of patients)
* [The Molecular Signatures Database (MSigDB)](http://software.broadinstitute.org/gsea/msigdb)
* [2016 Immunity Paper describing the creation of Immunological Signatures’ collection (C7)](http://diytranscriptomics.github.io/Reading/files/ImmuneSigDB.pdf)
* You know how I feel about Venn diagrams, so if you’re interested in exploring interactions between many groups of genes, have a look at this [Nature Methods paper](http://diytranscriptomics.github.io/Reading/files/upSet_plot.pdf), the accompanying R package, [UpSetR](https://cran.r-project.org/web/packages/UpSetR/README.html), as well as the [UpSet website](http://caleydo.org/tools/upset/). Note, there’s a shiny app for this as well!

