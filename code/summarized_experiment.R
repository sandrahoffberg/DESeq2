library(DESeq2)
library(ggplot2)

source("config.R", local=TRUE)

deseq_object <- readRDS(counts)

deseq2Data <- DESeqDataSet(deseq_object, design= formula(design))
# for sample data, the design should be "~ cell + dex"

condition <- colData(deseq_object)[[condition_name]]
# for sample data, the condition should be "dex"


source("filter_plot.R", local=TRUE)
