# 12 - Documenting your analyses with Rmarkdown

## Overview

In order to document the code and results from a bioinformatics workflow you’ll be introduced to [R Markdown](http://rmarkdown.rstudio.com/) and [Knitr](http://yihui.name/knitr/), and will use these tools to wrap up all your code and outputs (figures, tables, etc) together in a dyanmic document that can be placed in your lab notebook or published as a supplementary file in your manuscript.

## Learning objectives

* Discuss the basic building blocks for assembling a figure
* Understand the basics of Rmarkdown and dynamic documents
* Learn how to construct an Rmarkdown doc using ‘essential’ code chunks from the course
* Make a flexdashboard the incorporates graphics and code

## Code

* [Rmarkdown_template.Rmd](https://www.dropbox.com/s/f1c610tx0360cbd/Rmarkdown_template.Rmd?dl=0) - a skeletal template that you can fill out with the code ‘essentials’ that we covered in class
* [Rmarkdown_essentials.Rmd](https://www.dropbox.com/s/elui0r4x3uoxvzw/Rmarkdown_essentials.Rmd?dl=0) - same template as above, but already filled in with essentials
* [flexdashboard_essentials.Rmd](https://www.dropbox.com/s/pp1oow3u8wnfp2i/flexdashboard_essentials.Rmd?dl=0) - a different style of Rmarkdown that puts the focus on creating a dashboard of graphics. Particularly powerful when combined with the interactive graphs we produced during the course

## Reading

[Blog post on Knitr by Karl Broman](http://kbroman.org/knitr_knutshell/)

### Other examples of markdown in action

Supplementary code files from some of our recent papers:

* [Eosinophils and IL-4 Support Nematode Growth Coincident with an Innate Response to Tissue Injury](https://journals.plos.org/plospathogens/article?id=10.1371/journal.ppat.1005347)
* [Life cycle progression and sexual development of the apicomplexan parasite Cryptosporidium parvum](https://www.nature.com/articles/s41564-019-0539-x)
* [Variable gene expression and parasite load predict treatment outcome in cutaneous leishmaniasis](https://stm.sciencemag.org/content/11/519/eaax4204)
* [Make your own CV in markdown](https://stm.sciencemag.org/content/11/519/eaax4204) - uses the Vitae and Scholar packages in R. Here’s mine on github
* [This course website!](https://github.com/DIYtranscriptomics/DIYtranscriptomics.github.io) - Each page on the course website is just a simple markdown document.
* [My lab website](http://hostmicrobe.org/) is just a bunch of markdown files served from github (repo [here](https://github.com/hostmicrobe/hostmicrobe.github.io))