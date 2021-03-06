# The COVID collaborative challenge

## Background on the dataset

You will be given raw counts and a study design file for an extensive RNA-seq exploration of the host response to SARS-CoV-2 and related viruses in different cell lines, tissues and time points, and which include human primary samples as well as _in vivo_ studies in ferrets - an incredible dataset from Benjamin TenOever’s lab at Mt Sinai! You should check out their recent BioRxiv preprint: [SARS-CoV-2 launches a unique transcriptional signature from in vitro, ex vivo, and in vivo systems.](https://doi.org/10.1101/2020.03.24.004655)

## The challenge

This challenge certainly falls into the category of 'easier said than done'.  
Your challenge is to look at the study design file and decide on a question YOU would like to ask using this highly multivariate dataset, and then carry out your analysis using code from the course.

## Logistics

You have 1 week until we meet for the scheduled hackdash on the last day of the course. Use that time to work through the dataset either alone, or together with your classmates (don’t forget that you can use the call/videoconference function directly within the course Slack page).  
Our hackdash will follow the same format as usual, with random assignment to breakout rooms of 3-4 people. Once in your room, you will need to discuss the various approaches/questions that each person pursued, then bring your different perspectives together to produce a final figure or set of figures that illustrates what you found.  
Your group’s final figure(s) must be turned in at the conclusion of the hackdash, along with a description of your key findings.

## Key concepts

* Thinking about what you want to ask with a large dataset, priortizing questions when there are many possible things you _could ask_, and then constructing your analysis around one or a few key questions
* Flying solo on a data analysis project, then coming together with collaborators at different points in the project to see how colleagues differ in their approach and perspective, then incorporating these different perspectives in a final product is a key part of the research process when large datasets are involved

## What you’ll need to get started

* [Human expression data](http://diytranscriptomics.github.io/Data/files/GSE147507_RawReadCounts_Human.tsv) - raw counts obtained as a table, obtained directly from the public GEO record
* [Ferret expression data](http://diytranscriptomics.github.io/Data/files/GSE147507_RawReadCounts_Ferret.tsv) - raw counts obtained as a table, obtained directly from the public GEO record
* [Study design file](http://diytranscriptomics.github.io/Data/files/covid_metadata.txt) - Manually assembled and curated from sample descriptions in GEO
* [Rscript](http://diytranscriptomics.github.io/Data/files/loadData.R) - to get you over the hump of importing the three files above into R

## Tips

* You are not required to produce an Rmarkdown document for this challenge, but you can certainly do so if you prefer this over making a figure as a way to show how you worked through the data
* You can ignore the step 1 script for this challenge, since we didn’t align the data, but rather got raw counts directly from the GEO entry
* Since the data is already in the form of a count table (genes as rows and samples as columns) you don’t need to worry about annotations either. Go straight to creating a DGEList object
* dplyr is going to be critical in this challenge, as you will need to wrangle the study design file and the raw count tables to get what you need to address your question(s)
