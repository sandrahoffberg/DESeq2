library(DESeq2)
library(ggplot2)

source("config.R", local=TRUE)

# Read in the samples and the experimental condition
samples <- read.csv(metadata)

sampleTable <- data.frame(sampleName = samples$file,
                          fileName = transcripts,
                          condition = samples$condition)

deseq2Data <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = "../data",
                                       design= ~ condition)

condition <- "condition"

source("filter_plot.R", local=TRUE)
