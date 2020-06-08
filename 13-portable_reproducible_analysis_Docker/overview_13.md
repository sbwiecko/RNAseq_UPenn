# 13 - Making your analysis portable and reproducible with Docker

## Overview

Reproducing an analysis requires more than just code. You need the original raw data, access to the appropriate programming languages, and application specific packages (and often specific versions of these packages). This poses a major impediment to reproducibility, even for researchers with a background in bioinformatics. To address this challenge, you’ll learn how to ‘containerize’ your data, scripts and software into a code timecapsule, making it easy to share and rerun an entire analysis with the push of a button.

## Learning objectives

* Learn how to make your research analyses reproducible
* Create a reproducible package environment with [renv](https://rstudio.github.io/renv/articles/renv.html)
* Share your project via [GitHub](https://github.com/) and git
* Understand how to streamline your code by making your own R functions
* Share your work as an R package
* Discuss the basics of Docker and containerized software
* Use CodeOcean to access the entire DIY course as a reproducible computing environment

### What you need to do

* Sign-up for a an account on [Code Ocean](https://codeocean.com/) (be sure to use the same email address from your DataCamp login). You’ll get 15 compute hrs/month for free
* Sign-up for a free [GitHub](https://github.com/) account (doesn’t matter which email you use)
* Download [this gitignore file](https://github.com/DIYtranscriptomics/DIYtranscriptomics.github.io/tree/master/Code/files/gitignore.txt) - useful for updating your own .gitignore file in a project repo

## Code

We’ll use Code Ocean to interact with a dockerized container that packs all the code, data and software from the course into one reproducible and web-accessible environment. Simply login (or set-up a free account if you don’t already have one) and you’ll be able to re-run the entire course in a matter of minutes, without any software installation or data download. Your first run may take ~15min, since the full computing environment must be being built, but subsequent runs will be much faster. Note that this capsule includes raw fastq files, kallisto outputs, all of MSigDB for running GSEA, and the entire ARCHS4 database for interrogating ~700,000 publically available mouse and human RNAseq datasets. Have fun adapting this capsule for your own analyses!

## Reading

* [Happy Git and GitHub with RStudio](https://happygitwithr.com/) - Jenny Bryan and team walk through every step of how to install git, connect to GitHub and access version control from within RStudio
* [Code Ocean whitepaper](http://diytranscriptomics.github.io/Reading/files/codeOcean_whitepaper.pdf) describing the need for better tools for reproducible research and introducting their cloud-based computational platform for addressing this need
* [Intro to Docker](https://docker-curriculum.com/) - CodeOcean is based on Docker, a free and open-source tool that allows you to ‘build, share and run applications anywhere’
* [On-boarding document](https://docker-curriculum.com/) for gettign started with Code Ocean
* [Our recent paper](https://stm.sciencemag.org/content/11/519/eaax4204), showing a code capsule embedded directly in the joural webpage (a first for any AAAS journal)
* [Google Collaboratory](https://colab.research.google.com/) - Write, edit and share Python code directly in your browser
