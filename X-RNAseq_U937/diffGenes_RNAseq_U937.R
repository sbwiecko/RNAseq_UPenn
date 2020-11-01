# INTRODUCTION ----
# The goal of this script is to identify differentially expressed genes (DEGs)
# in Undifferentiated vs. Differentiated U937 using public RNAseq dataset
# from Haney et al.,"Identification of Phagocytosis Regulators Using Magnetic
# Genome-Wide CRISPR Screens" with accession id GSE107566.


# LOAD PACKAGES ----
library(tidyverse) # provides access to R packages for data science
library(rhdf5) # to access data from ARCHS4 (BiocManager::install(rhdf5))
library(edgeR) # well known package for differential expression analysis, but
        # we only use for the DGEList object and for normalization methods
library(matrixStats) # calculate stats on rows or columns of a data matrix


# LOAD ARCHS4 DATABASE ----
# we should have already downloaded the most recent versions of human
# RNAseq data from ARCHS4 database in hdf5 format
archs4_human <- "./human_matrix.h5"

# use the h5 list (h5ls) function from the rhdf5 package to look at the
# contents of these databases
h5ls(archs4_human)
# 238,522 samples from human
all_samples_human <- h5read(archs4_human, name = "meta/Sample_geo_accession")


# QUERY ARCHS4 DATABASE ----
# choose your samples based on GEO or SRA ID
# see https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE107566
samples <- c("GSM2310941", # Undifferentiated Ctrl sgRNA 1
             "GSM2871539", # Undifferentiated Ctrl sgRNA 2
             "GSM2871540", # Undifferentiated Ctrl sgRNA 3
             "GSM2871541", # Undifferentiated Ctrl sgRNA 4
             "GSM2871546", # Differentiated Ctrl sgRNA 1
             "GSM2871547", # Differentiated Ctrl sgRNA 2
             "GSM2871548", # Differentiated Ctrl sgRNA 3
             "GSM2871549"  # Differentiated Ctrl sgRNA 4
            )

# identify columns to be extracted from ARCHS4 database
sample_locations <- which(all_samples_human %in% samples)

#####
# !!! apparently our samples are not in the database :-(
# so the next steps are not usefull in our current project
#####

# extract gene symbols from the metadata
genes <- h5read(archs4_human, "meta/genes")
# extract expression data from ARCHS4
expression <- h5read(archs4_human, "data/expression",
                     index = list(1 : length(genes),
                     sample_locations)
                    )
rownames(expression) <- genes
colnames(expression) <- all_samples_human[sample_locations]
colSums(expression) # this shows the sequencing depth for each of the samples
archs4_dgelist <- DGEList(expression)
archs4_cpm <- cpm(archs4_dgelist)
colSums(archs4_cpm)


# SET UP OF THE DESIGN MATRIX ----
targets <- read_tsv("studydesign.txt")
sample_labels <- targets$sample
# gene annotation is already provided in the data file as follow:
# "RNA-seq data was mapped to the human genome annotation hg38
# resulting in 8 to 12 million mapped reads per sample.
# CBMH2_Counts_GEO.xlsx contains counts from **DESeq2**"
# Therefore, no need to annotate the read counts (EnsemblDB and TxImport)
# and no need to use edgeR or Limma
# note -the data are raw counts with 4 different sample types

data <- read_excel("GSE107566_CBMH2_Counts_GEO.xlsx")
data <- data[sample_labels]

# Make a DGElist from your counts, and plot ####
data_DGEList <- DGEList(data)
