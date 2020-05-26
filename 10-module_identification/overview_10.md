# 10 - Module Identification

Lists of differentially expressed transcripts often include different patterns or modules of genes that are coordinately regulated across treatments or conditions, and these patterns can provide powerful insight into biology. In this class you’ll use correlation-based clustering methods and heatmap visualization to interrogate DEGs to reveal modules of co-regulated genes.

## Learning objectives

* Be able to interpret and construct heatmaps
* Understand color choice in R
* Understand how clustering methods are used to identify coordinately expressed genes (a.k.a modules)
* Learn to use the command-line [clust](https://github.com/BaselAbujamous/clust) program

### Downloads

If you want to try _Clust_ on your own, you’ll need to install the program first (see [github page](https://github.com/BaselAbujamous/clust) for instructions). You can test it out using the Schisto hackdash dataset to reproduce what I showed you in class. There are two files you’ll need for this:

1. [this text file of DEGs](https://drive.google.com/file/d/1OgrR7YbSuhbFxvwvGdHLRViJixWr2bGZ/view) in female LE-strain worms; and
2. [a reps file](https://drive.google.com/file/d/1qv5x-MHqg-bh9OllAf-JsE1eB-AkEgIC/view) that maps the columns (samples) to groups (conditions)

You’ll need to read these two files into your R environment before running _Clust_. Run _clust_ in terminal:

```bash
clust DiffGenes.txt -n 101 4 -r reps.txt
```

You can view an example of the output [here](https://drive.google.com/drive/folders/1BWVl42rhzC1Kd7GA5JF0OwX73bAUd_mm).

### Reading

* [Clust: automatic extraction of optimal co-expressed gene clusters from gene expression data. Genome Biology, 2018.](https://doi.org/10.1186/s13059-018-1536-8)
* Colors in R - Colors palettes are an often underappreciated aspect of making beautiful and informative plots in R. You can access a suite of color palettes using the [RColorBrewer package](http://colorbrewer2.org/). These palettes can be viewed in this [cheatsheet](http://diytranscriptomics.github.io/Reading/files/colorbrewerPalettes.pdf). Unfortunately, these standard palettes often don’t cut it, and you’ll need custom palettes. For this, I love using [Sip](https://sipapp.io/) to pick, organize and access color palettes
