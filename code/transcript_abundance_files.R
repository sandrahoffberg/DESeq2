library(DESeq2)
library(ggplot2)
library(tximport)

source("config.R", local=TRUE)

# Read in the samples and the experimental condition
samples <- read.csv(metadata)

# Read in the table that says which transcript belong to which gene
tx2gene <- read.csv(transcripts_genes)

txi <- tximport(transcripts, type=analysis, tx2gene=tx2gene, ignoreAfterBar = TRUE)

# Create the DEseq2DataSet object
deseq2Data <- DESeqDataSetFromTximport(txi,
                                   colData = samples,
                                   design = ~ condition)


condition <- condition_name
#condition

source("filter_plot.R", local=TRUE)
