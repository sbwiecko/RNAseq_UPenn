# The Malaria challenge
# script written by Sébastien Wieckowski, 2020

# Let's first get the data into R

library(tidyverse)
library(tximport)
library(edgeR)
library(matrixStats)
library(cowplot)
library(gplots)

targets <- read_tsv("./malaria/studyDesign.txt")
path <- file.path("./malaria", targets$sample, "abundance.tsv")
all(file.exists(path)) # worked great

txi <- tximport(path,
                type = "kallisto",
                txOut = TRUE,      # we import data at the transcript level
                countsFromAbundance = "lengthScaledTPM",
                ignoreTxVersion = TRUE)

sample_labels <- targets$sample
# just need to delete the redundant prefix on each label
sample_labels <- substr(sample_labels, 10, nchar(sample_labels))

# we 1st prepare a df of unfiltered, non-normalized data
# it will always be log2.cpm.df.pivot from now
dgelist <- DGEList(txi$counts) # appropriate format
data_unfilt_non_norm <- cpm(dgelist, log = TRUE) %>%
    as_tibble(, rownames = "transcriptID") # EdgeR utilizes counts for stats
# reorder colnames
colnames(data_unfilt_non_norm) <- c("transcriptID", sample_labels)

# create the final pivoted dataframe for unfiltered/non-normalized data
data_unfilt_non_norm <- pivot_longer(data_unfilt_non_norm,
                                     cols = 2:15,
                                     # looked at names(data_unfilt_non_norm)
                                     names_to = "sample",
                                     values_to = "expression")

# we now prepare a df of filtered, non-normalized data
# a quick look at the distribution of prev data for chosing the threshold value
ggplot(data_unfilt_non_norm) +
    aes(x = sample, y = expression) +
    geom_violin(trim = FALSE, show.legend = FALSE)
# it look very strange, the second replicate (rep2) for each sample has
# lower values systematically
# nevertheless we select threshold at 10^-1 based on that graph

# and because we deal with duplicate we keep data from 2 or more samples
cpm <- cpm(dgelist, log = TRUE)
keepers <- rowSums(cpm > -1) >= 2
dgelist_filt <- dgelist[keepers, ]
# check dim(dgelist_filt) and dim(dgelist) to get number of genes filtered
data_filt_non_norm <- cpm(dgelist_filt, log = TRUE) %>%
    as_tibble(, rownames = "transcriptID") # EdgeR utilizes counts for stats
colnames(data_filt_non_norm) <- c("transcriptID", sample_labels)

# create the final pivoted dataframe for filtered/non-normalized data
data_filt_non_norm <- pivot_longer(data_filt_non_norm,
                                   cols = 2:15,
                                   names_to = "sample",
                                   values_to = "expression")
# quick check on the new data
ggplot(data_filt_non_norm) +
    aes(x = sample, y = expression) +
    geom_violin(trim = FALSE, show.legend = FALSE)
# a bit better but not stellar...

# finally we prepare the df for filterd and normalized data
dgelist_filt_norm <- calcNormFactors(dgelist_filt, method = "TMM")

data_filt_norm <- cpm(dgelist_filt_norm, log = TRUE) %>%
    as_tibble(, rownames = "transcriptID") # EdgeR utilizes counts for stats
colnames(data_filt_norm) <- c("transcriptID", sample_labels)
data_filt_norm <- pivot_longer(data_filt_norm,
                               cols = 2:15,
                               names_to = "sample",
                               values_to = "expression")

# quick check on the new data
ggplot(data_filt_norm) +
    aes(x = sample, y = expression) +
    geom_violin(trim = FALSE, show.legend = FALSE)

# we create one PDF for the grid plot
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

plot <- plot_grid(p1, p2, p3,
                  labels = c("A", "B", "C"),
                  label_size = 12)
ggsave(filename = "cpm_all_grid.pdf", plot = plot)

# time for data analysis
# PCA
pca_res <- prcomp(t(cpm(dgelist_filt_norm, log = TRUE)))
# look at pca_res in environment
summary(pca_res) # Prints variance summary for all principal components
# we observe cumul prop var 72.85% at PC4
screeplot(pca_res) # elbow clearly observed between PC2 and PC3
# variance values needed for the plots
pc_var <- pca_res$sdev^2 # sdev^2 captures these eigenvalues from the PCA result
pc_per <- round(pc_var / sum(pc_var) * 100, 1) # %variance explained by each PC

# lets then plot PC2 vs PC1
pca_res_df <- as_tibble(pca_res$x)
ggplot(pca_res_df) +
    aes(x = PC1, y = PC2, label = sample_labels,
        color = group,
        size = 10) +
    geom_point(alpha = .8) +
    geom_label() +
    xlab(paste0("PC1 (", pc_per[1], "%", ")")) +
    ylab(paste0("PC2 (", pc_per[2], "%", ")")) +
    labs(title = "PCA plot",
         caption = paste0("produced on ", Sys.time())) +
    coord_fixed() +
    theme_light()

# looks great and makes a kind of cycle from 0 to 48hrs with 48hrs similar to 0
ggsave(filename = "PCA.pdf", dpi = 150)

# There is only one variable in this study
group <- targets$group
group <- factor(group,
                levels = c("0hrs", "8hrs", "16hrs",
                           "24hrs", "32hrs", "40hrs", "48hrs"),
                labels = c('A', 'B', 'C', 'D', 'E', 'F', 'G'))
# let's do some stat using time as a variable
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

v.DEGList.filtered.norm <- voom(dgelist_filt_norm, design, plot = TRUE)
# the graph show that we could take a threshold a little bit higher
fit <- lmFit(v.DEGList.filtered.norm, design)
# I excluded the 48hrs - 0hrs test because these are similar
contrast.matrix <- makeContrasts(B - A, C - A, D - A, E - A, F - A,
                                 levels = design)

fits <- contrasts.fit(fit, contrast.matrix)
ebFit <- eBayes(fits)

# Volcano Plots ----
# in topTable function above, set 'number=40000' to capture all genes
Hits <- topTable(ebFit, adjust = "BH", coef = 1, number = 40000, sort.by = "logFC") %>%
  as_tibble(rownames = "geneID")

# we choose less restrictive threshold because of stat power (duplicate)
ggplot(Hits) +
  aes(y = -log10(adj.P.Val), x = logFC, text = paste("Symbol:", geneID)) +
  geom_point(size = 2) +
  annotate("rect", xmin = 2.5, xmax = 7, ymin = -log10(0.05), ymax = 2.5, alpha=.2, fill="#BE684D") +
  annotate("rect", xmin = -2.5, xmax = -7, ymin = -log10(0.05), ymax = 2.5, alpha=.2, fill="#2C467A") +
  labs(title = "Volcano plot",
       subtitle = "Malaria challenge",
       caption = paste0("produced on ", Sys.time())) +
  theme_bw()

# Venn diagram with those parameters
# NOTE - not more than 5 sets usable with this tool...
results <- decideTests(ebFit, method="global", adjust.method="BH", p.value=0.05, lfc=2.5)
vennDiagram(results,
            include = "both",
            lwd = 2,
            circle.col = c("blue","green","red","yellow","black"))

# Let's build a heatmap
colnames(v.DEGList.filtered.norm$E) <- sample_labels
diffGenes <- v.DEGList.filtered.norm$E[results[, 1] != 0, ]
clustRows <- hclust(as.dist(1 - cor(t(diffGenes),
                            method = "pearson")
                   ), method = "complete") 
clustColumns <- hclust(as.dist(1 - cor(diffGenes,
                                       method="spearman")
                              ), method = "complete") #cluster columns by spearman correlation
heatmap.2(diffGenes, 
          Rowv = as.dendrogram(clustRows),
          Colv = as.dendrogram(clustColumns),
          scale = 'row', labRow = NA,
          density.info = "none", trace = "none",
          cexRow = 1, cexCol = 1, margins = c(8,20))

