library(DESeq2)
library(ggplot2)

source("config.R", local=TRUE)

# Read in the raw read counts
rawCounts <- read.csv(counts)

# Read in the sample mappings
sampleData <- read.csv(metadata)

### Construct DESEQDataSet Object

# Convert count data to a matrix of appropriate form that DEseq2 can read
geneID <- rawCounts$Gene.ID
sampleIndex <- setdiff(colnames(rawCounts), c("Gene.ID", "Gene.Name"))
rawCounts <- as.matrix(rawCounts[,sampleIndex])
rownames(rawCounts) <- geneID
head(rawCounts)

# Convert sample variable mappings to an appropriate form that DESeq2 can read
head(sampleData)
rownames(sampleData) <- sampleData$Run
keep <- setdiff(strsplit(design, " ")[[1]], c('~', '+'))
sampleData <- sampleData[,keep]
sampleData <- apply(sampleData, 2, factor)
head(sampleData)

# Put the columns of the count data in the same order as rows names of the sample mapping, then make sure it worked
rawCounts <- rawCounts[,unique(rownames(sampleData))]
all(colnames(rawCounts) == rownames(sampleData))

# Create the DEseq2DataSet object
deseq2Data <- DESeqDataSetFromMatrix(countData=rawCounts, colData=sampleData, design= formula(design))
# for sample data, the design should be "~ Sample.Characteristic.individual. + Sample.Characteristic.biopsy.site."

condition = condition_name
# for sample data, the condition should be "Sample.Characteristic.biopsy.site."

source("filter_plot.R", local=TRUE)

