library(tidyverse)
library(edgeR)
library(tximport)
library(ensembldb)
library(EnsDb.Hsapiens.v86)

# read in your study design
targets <- read_tsv("../3-read_mapping_Kallisto/test/studydesign.txt")
sampleLabels <- targets$sample

# set file paths to your mapped data
path <- file.path("../3-read_mapping_Kallisto/test", targets$sample, "abundance.tsv")
Tx <- transcripts(EnsDb.Hsapiens.v86,
                  columns = c("tx_id", "gene_name"))
Tx <- as_tibble(Tx)
Tx <- dplyr::rename(Tx, target_id = tx_id)
Tx <- dplyr::select(Tx, "target_id", "gene_name")
Txi_gene <- tximport(path,
                     type = "kallisto",
                     tx2gene = Tx,
                     txOut = FALSE,
                     countsFromAbundance = "lengthScaledTPM",
                     ignoreTxVersion = TRUE)


myDGEList <- DGEList(Txi_gene$counts)
# first run the function in profile_function.R
profile_function(data = myDGEList,
                 samples = sampleLabels)

cpm <- cpm(myDGEList)
keepers <- rowSums(cpm > 1) >= 5 #user defined
myDGEList_filtered <- myDGEList[keepers, ]
profile_function(data = myDGEList_filtered,
                 samples = sampleLabels)

myDGEList_filtered_norm <- calcNormFactors(myDGEList_filtered,
                                           method = "TMM")
profile_function(data = myDGEList_filtered_norm,
                 samples = sampleLabels)