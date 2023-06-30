library(DESeq2)
library(ggplot2)
library(tximport)

source("config.R", local=TRUE)

# Read in the samples and the experimental condition
samples <- read.csv(metadata)
samples

# Read in the table that says which transcript belong to which gene
tx2gene <- read.csv("../data/transcript_abundance/tx2gencode.v27.csv")

txi <- tximport(transcripts, type="kallisto", tx2gene=tx2gene, ignoreAfterBar = TRUE)

# Create the DEseq2DataSet object
deseq2Data <- DESeqDataSetFromTximport(txi,
                                   colData = samples,
                                   design = ~ condition)

deseq2Data

condition <- "condition"

source("filter_plot.R", local=TRUE)
