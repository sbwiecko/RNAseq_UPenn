# Ultra-fast read mapping with Kallisto

## Overview

In this class, we’ll finally get down to the business of using [Kallisto](https://pachterlab.github.io/kallisto/about) for memory-efficient mapping of raw reads to a reference transcriptome. You’ll carry out this mapping in class, right on your laptop, while we discuss what’s happening ‘under the hood’ with Kallisto and how this compares to more traditional alignment methods. You’ll be introduced to using command line software and will learn about automation and reproducibility through shell scripts.

## Learning objectives

* Discuss the course dataset.
* Download and examine a reference transcriptome from Ensembl.
* Use Kallisto to construct an index from this reference file.
* Use Kallisto to map our raw reads to this index
* Talk a bit about how an index is built and facilitates read alignment

## Reading

### Papers and blogs posts on Kallisto

* [2016 Nature Biotech paper](http://diytranscriptomics.github.io/Reading/files/Kallisto.pdf) from Lior Pachter’s lab describing Kallisto

* [2017 Nature Methods paper](http://diytranscriptomics.github.io/Reading/files/sleuth.pdf) from Lior Pachter’s lab describing Sleuth

* [Lior Pachter’s blog post on Kallisto](https://liorpachter.wordpress.com/2015/05/10/near-optimal-rna-seq-quantification-with-kallisto/)


* [blog post on pseudoalignments](http://tinyheero.github.io/2015/09/02/pseudoalignments-kallisto.html) - helps understand how Kallisto maps reads to transcripts. Did you notice that Kallisto is using ‘Expectation Maximization (EM)’ during the alignment? You can read more about what this is [here](http://diytranscriptomics.github.io/Reading/files/EM.pdf)

* [Kallisto discussions/questions](https://groups.google.com/forum/#!forum/kallisto-sleuth-users) and [Kallisto announcements](https://groups.google.com/forum/#!forum/kallisto-sleuth-announcements) are available on Google groups

### General info about ultra lightweight methods for transcript quantification

* [2014 Nature Biotech paper](http://diytranscriptomics.github.io/Reading/files/Sailfish.pdf) - describes Sailfish, which implimented the first lightweight method for quantifying transcript expression.

* [Not quite alignments](http://robpatro.com/blog/?p=248) - Rob Patro, the first author of the Sailfish paper, wrote a nice blog post comparing and contrasting alignment-free methods used by Sailfish, Salmon and Kallisto.

* [2018 Nature Methods paper describing Salmon](https://www.nature.com/articles/nmeth.4197) - A lightweight aligment tool from Rob Patro and Carl Kinsford. [Check out the website too](https://combine-lab.github.io/salmon/).

* [2011 Nature Biotechnology](http://diytranscriptomics.github.io/Reading/files/deBruijn.pdf) - Great primer to better understand what de Bruijn graph is.

Greg Grant’s recent paper comparing different aligners. This should be a helpful guide in choosing alignment software outside of what we used in class.

## Mapping single-end data

```powershell
kallisto quant `
-i Homo_sapiens.GRCh38.cdna.all.index `
-o test `
-t 4 `
--single -l 250 -s 30 `
SRR8668755_1M_subsample.fastq `
*> test.log
```