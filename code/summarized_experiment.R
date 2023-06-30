library(DESeq2)
library(ggplot2)

source("config.R", local=TRUE)

deseq_object <- readRDS(counts)

deseq2Data <- DESeqDataSet(deseq_object, design= ~ cell + dex)


condition <- condition_name
#dex

source("filter_plot.R", local=TRUE)
