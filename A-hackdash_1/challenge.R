# This is my script for the worm challenge
# Sébastien Wieckowski, 2020-05-14

# A colleague has asked for your help to mine data produced from a very
# large RNAseq study of the parasitc worm, _Schistosoma mansoni_. In
# this experiment, male (M), female (F), juvenile (J) and mixed sex (X)
# worms were recovered from infected mice at various timepoints (control,
# 3hr, 12hr, and 24hr) following _in vivo_ treatment with a low dose of the
# frontline anti-parasitic drug, praziquantel. Experiments were carried out
# with three different strains of worms: NMRI, LE, and LEPZQ. In total, 144
# samples were sequenced! To begin this challenge, you’ll need to download
# the Kallisto read mapping results and study design file [here]
# (https://www.dropbox.com/s/k0959o6luvda69p/schistosoma.zip?dl=0) and read
# these data into a new R project.

# * don’t try to annotate the data, just import at the transcript level

# very first step, import of all libraries
library(tidyverse)
library(tximport)
library(edgeR)
#library(matrixStats)
library(gt)

# get in the data
# read in the study design
targets <- read_tsv("./schistosoma/studydesign.txt")

# we only keep the columns necessary for the study as rep, treatment,
# dataType, sampleAccession, runAccession, internalRunID are not necessary
feats <- c("sample", "sex", "strain", "drugSensitivity", "timpoint")
targets <- targets %>%
    dplyr::select(feats)

# now we build the sample set based on targets$sample
path <- file.path("./schistosoma", targets$sample, "abundance.tsv")
all(file.exists(path)) # worked great

# final step for import
txi_transcripts <- tximport(path,
                            type = "kallisto",
                            txOut = TRUE, # we keep transcripts
                            countsFromAbundance = "lengthScaledTPM")
                            #ignoreTxVersion = TRUE

### Question 1
# Use plots to show the impact of filtering and normalization on the data
sample_labels <- targets$sample

# we 1st prepare a df of unfiltered, non-normalized data
# it will always be log2.cpm.df.pivot from now
dgelist <- DGEList(txi_transcripts$counts) # appropriate format
data_unfilt_non_norm <- cpm(dgelist, log = TRUE) %>%
    as_tibble(, rownames = "transcriptID") # EdgeR utilizes counts for stats
# reorder colnames
colnames(data_unfilt_non_norm) <- c("transcriptID", sample_labels)

# create the final pivoted dataframe for unfiltered/non-normalized data
data_unfilt_non_norm <- pivot_longer(data_unfilt_non_norm,
                                     cols = M1h_LEPZQ_1:X24h_LEPZQ_3, # or cols = 2:145
                                     # looked at names(data_unfilt_non_norm)
                                     names_to = "sample",
                                     values_to = "expression")

# we now prepare a df of filtered, non-normalized data
# a quick look at the distribution of prev data for chosing the threshold value
ggplot(data_unfilt_non_norm) +
    aes(x = sample, y = expression, fill = sample) +
    geom_violin(trim = FALSE, show.legend = FALSE)

# many data points below log2=0, therefore we choose to filter > 2^0 i.e. 1
# and because we deal with triplicate we keep data from 3 or more samples
cpm <- cpm(dgelist, log = TRUE)
keepers <- rowSums(cpm > 1) >= 3
dgelist_filt <- dgelist[keepers, ]
# check dim(dgelist_filt) and dim(dgelist) to get number of genes filtered
data_filt_non_norm <- cpm(dgelist_filt, log = TRUE) %>%
    as_tibble(, rownames = "transcriptID") # EdgeR utilizes counts for stats
colnames(data_filt_non_norm) <- c("transcriptID", sample_labels)

# create the final pivoted dataframe for filtered/non-normalized data
data_filt_non_norm <- pivot_longer(data_filt_non_norm,
                                   cols = M1h_LEPZQ_1:X24h_LEPZQ_3,
                                   names_to = "sample",
                                   values_to = "expression")
# quick check on the new data
ggplot(data_filt_non_norm) +
    aes(x = sample, y = expression) +
    geom_violin(trim = FALSE, show.legend = FALSE)

# finally we prepare the df for filterd and normalized data
dgelist_filt_norm <- calcNormFactors(dgelist_filt, method = "TMM")

data_filt_norm <- cpm(dgelist_filt_norm, log = TRUE) %>%
    as_tibble(, rownames = "transcriptID") # EdgeR utilizes counts for stats
colnames(data_filt_norm) <- c("transcriptID", sample_labels)
data_filt_norm <- pivot_longer(data_filt_norm,
                               cols = M1h_LEPZQ_1:X24h_LEPZQ_3,
                               names_to = "sample",
                               values_to = "expression")

# quick check on the new data
ggplot(data_filt_norm) +
    aes(x = sample, y = expression) +
    geom_violin(trim = FALSE, show.legend = FALSE)

# we create one PDF for each plot
p1 <- ggplot(data_unfilt_non_norm) +
    aes(x = sample, y = expression, fill = sample) +
    geom_violin(trim = FALSE, show.legend = FALSE) +
    stat_summary(fun = "median",
                 geom = "point",
                 shape = 95,
                 size = 5,
                 color = "black",
                 show.legend = FALSE) +
    labs(y = "log2 expression", x = "sample",
         title = "Log2 Counts per Million (CPM)",
         subtitle = "unfiltered, non-normalized",
         caption = paste0("produced on ", Sys.time())) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, size = 4))

p2 <- ggplot(data_filt_non_norm) +
    aes(x = sample, y = expression, fill = sample) +
    geom_violin(trim = FALSE, show.legend = FALSE) +
    stat_summary(fun = "median",
                 geom = "point",
                 shape = 95,
                 size = 5,
                 color = "black",
                 show.legend = FALSE) +
    labs(y = "log2 expression", x = "sample",
         title = "Log2 Counts per Million (CPM)",
         subtitle = "filtered, non-normalized",
         caption = paste0("produced on ", Sys.time())) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, size = 4))


p3 <- ggplot(data_filt_norm) +
    aes(x = sample, y = expression, fill = sample) +
    geom_violin(trim = FALSE, show.legend = FALSE) +
    stat_summary(fun = "median",
                 geom = "point",
                 shape = 95,
                 size = 5,
                 color = "black",
                 show.legend = FALSE) +
    labs(y = "log2 expression", x = "sample",
         title = "Log2 Counts per Million (CPM)",
         subtitle = "filtered, normalized",
         caption = paste0("produced on ", Sys.time())) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, size = 4))


ggsave(filename = "cpm_unfilt_non_norm.pdf", plot = p1,
       width = 25, height = 15, units = "cm", dpi = 150)
ggsave(filename = "cpm_filt_non_norm.pdf", plot = p2,
       width = 25, height = 15, units = "cm", dpi = 150)
ggsave(filename = "cpm_filt_norm.pdf", plot = p3,
       width = 25, height = 15, units = "cm", dpi = 150)
# the median homogenized in all samples after normalization to ca. log2=5

### Question 2
# We look at the variables sex, strain and timpoint,
# not sure what/how to do with variable drugSensitivity

# capture experimental variables as factors from this study design
# sex abreviated
sex <- factor(targets$sex,
              levels = c("juvenile", "female", "male", "mixed"),
              labels = c("J", "F", "M", "X"))
strain <- factor(targets$strain)
# for time we keep levels in order and give a float label
# be carefull there are more unique values than in the instructions
# >unique(targets$timpoint)
time <- factor(x = targets$timpoint,
               levels = c("control", "30min", "1hr", "3hr", "12hr", "24hr"),
               labels = c(0, 0.5, 1, 3, 12, 24))
sample_names <- targets$sample

# Principal component analysis (PCA) on filtered and normalized cpm DGElist data
pca_res <- prcomp(t(cpm(dgelist_filt_norm, log = TRUE)))

# look at pca_res in environment
summary(pca_res) # Prints variance summary for all principal components
# we observe cumul prop var 72.85% at PC4
screeplot(pca_res) # elbow observed around PC2

# variance values needed for the plots
pc_var <- pca_res$sdev^2 # sdev^2 captures these eigenvalues from the PCA result
pc_per <- round(pc_var / sum(pc_var) * 100, 1) # %variance explained by each PC

# lets then plot PC2 vs PC1
pca_res_df <- as_tibble(pca_res$x)
ggplot(pca_res_df) +
    aes(x = PC1, y = PC2, label = sample_names,
        size = time,
        color = sex,
        shape = strain) +
    geom_point(alpha = .8) +
    xlab(paste0("PC1 (", pc_per[1], "%", ")")) +
    ylab(paste0("PC2 (", pc_per[2], "%", ")")) +
    labs(title = "PCA plot",
         caption = paste0("produced on ", Sys.time())) +
    coord_fixed() +
    theme_light()
# also ploting PC4 vs PC2 gives a nice separation based on strain

ggsave(filename = "PCA.pdf", dpi = 150)
# clearly sex explains the majority of the variance in this dataset
# shorter timepoints have mainly higher PC2
# X are F+M mixed together; J ~ M

# PCA 'small multiples' chart to understand impact of each sample on each PC
# we first prepare a tibble with data and variables
pca_res_df <- pca_res$x[, 1:4] %>%
    as_tibble() %>%
    add_column(sample = sample_names,
               time = time,
               sex = sex,
               strain = strain)
  
pca_pivot <- pivot_longer(pca_res_df,
                          cols = PC1:PC4,
                          names_to = "PC",
                          values_to = "loadings")

ggplot(pca_pivot) +
    aes(x = sample, y = loadings, fill = sex) +
    # we can iteratively 'paint' different covariates using the 'fill' aes
    geom_bar(stat = "identity") +
    facet_wrap(~PC) +
    labs(title = "PCA 'small multiples' plot",
         caption = paste0("produced on ", Sys.time())) +
    theme_bw() +
    theme(axis.text.y = element_blank()) +
    coord_flip()

ggsave(filename = "PCA_multi_sex.pdf", dpi = 150)
# this plot confirms the implication of sex on PC1 and PC2

ggplot(pca_pivot) +
    aes(x = sample, y = loadings, fill = strain) +
    # we can iteratively 'paint' different covariates using the 'fill' aes
    geom_bar(stat = "identity") +
    facet_wrap(~PC) +
    labs(title = "PCA 'small multiples' plot",
         caption = paste0("produced on ", Sys.time())) +
    theme_bw() +
    theme(axis.text.y = element_blank()) +
    coord_flip()

ggsave(filename = "PCA_multi_strain.pdf", dpi = 150)
# here we detect implication of strain on PC4

ggplot(pca_pivot) +
    aes(x = sample, y = loadings, fill = time) +
    # we can iteratively 'paint' different covariates using the 'fill' aes
    geom_bar(stat = "identity") +
    facet_wrap(~PC) +
    labs(title = "PCA 'small multiples' plot",
         caption = paste0("produced on ", Sys.time())) +
    theme_bw() +
    theme(axis.text.y = element_blank()) +
    coord_flip()

ggsave(filename = "PCA_multi_time.pdf", dpi = 150)
# some trend of longer timepoints on PC3


### Question 3
# identify the top 10 parasite genes induced by praziquantel treatment
# in female LE strain worms at the 24hr timepoint compared to control
# we first select data for F/LE at 0 and 24hr treatment from CPM dataset
data_filt_norm <- cpm(dgelist_filt_norm, log = TRUE) %>%
    as_tibble(, rownames = "transcriptID") # EdgeR utilizes counts for stats
colnames(data_filt_norm) <- c("transcriptID", sample_labels)
# data are in triplicate subset=
F_LE_samples = c("transcriptID", "FCtl_LE_1", "FCtl_LE_2", "FCtl_LE_3", "F24h_LE_1", "F24h_LE_2", "F24h_LE_3")
# better to do %>% dplyr::select(geneID, starts_with("F24h_LE_"), starts_with("FCtl_LE_"))
# we then compute the average expression and finally fold change for each gene
data_fc <- data_filt_norm[, F_LE_samples] %>%
    mutate(ctr_avg = (FCtl_LE_1 + FCtl_LE_2 + FCtl_LE_3) / 3,
           trt_avg = (F24h_LE_1 + F24h_LE_2 + F24h_LE_3) / 3,
           LogFC = (trt_avg - ctr_avg)) %>%
    mutate_if(is.numeric, round, 2) %>% # all data rounded
# we then arrange by ascending absolute logFC
    dplyr::arrange(desc(abs(LogFC)))
# and we select the top 10 gens
top10 <- data_fc %>%
    dplyr::top_n(10) # 12 values because if ties

print(top10[, c("transcriptID", "LogFC")])
# Produce publication-quality tables using the gt package
gt(top10[, c("transcriptID", "LogFC")]) %>% gtsave(filename = "top10_genes.png")

### Question bonus
# Smp_138080 : G4LZT9; MEG-3 (Grail) family
# Smp_138070 : G4LZU0; MEG-3 (Grail) family
# Smp_049300.3:G4V8Y7; Putative major egg antigen
# Smp_049300.1:G4V8Y2; Putative heat shock protein
# Smp_186020.1:G4V8Y2; Putative heat shock protein
# as described in UniProtKB