library(DESeq2)
library(ggplot2)

source("config.R", local=TRUE)

# Read in the raw read counts
rawCounts <- read.csv(raw_exp)
head(rawCounts)

# Read in the sample mappings
sampleData <- read.csv(metadata)
head(sampleData)


### Construct DESEQDataSet Object

# Convert count data to a matrix of appropriate form that DEseq2 can read
geneID <- rawCounts$Gene.ID
sampleIndex <- grepl("SRR\\d+", colnames(rawCounts))
rawCounts <- as.matrix(rawCounts[,sampleIndex])
rownames(rawCounts) <- geneID
head(rawCounts)

# Convert sample variable mappings to an appropriate form that DESeq2 can read
head(sampleData)
rownames(sampleData) <- sampleData$Run
keep <- c("Sample.Characteristic.biopsy.site.", "Sample.Characteristic.individual.")
sampleData <- sampleData[,keep]
colnames(sampleData) <- c("tissueType", "individualID")
sampleData$individualID <- factor(sampleData$individualID)
head(sampleData)

# Put the columns of the count data in the same order as rows names of the sample mapping, then make sure it worked
rawCounts <- rawCounts[,unique(rownames(sampleData))]
all(colnames(rawCounts) == rownames(sampleData))

# Create the DEseq2DataSet object
deseq2Data <- DESeqDataSetFromMatrix(countData=rawCounts, colData=sampleData, design= ~ individualID + tissueType)
deseq2Data

### Pre-filtering of data to reduce dataset size

# See what affect this filter will have.
dim(deseq2Data)
dim(deseq2Data[rowSums(counts(deseq2Data)) > filter, ])

# Perform pre-filtering of the genes that have less than 5 reads in all samples
deseq2Data <- deseq2Data[rowSums(counts(deseq2Data)) > filter, ]


### Run the DESEQ function and look at results table

deseq2Data <- DESeq(deseq2Data)

# Extract differential expression results
# For "tissueType" perform primary vs normal comparison
deseq2Results <- results(deseq2Data)#, contrast=c("tissueType", "primary colorectal cancer", "normal colonic epithelium"))

### Summary of differential gene expression

# View summary of results
summary(deseq2Results)

### Sort summary list by p-value

res <- deseq2Results[order(deseq2Results$padj),]
head(res)


### Plot the MA plot: log ratio (M) vs an average (A)

output_ma_file <- '../results/MA_plot.png'

png(output_ma_file)
plotMA(deseq2Results, alpha=alpha, main=paste("MA plot with alpha of", alpha, sep=" "))
dev.off()
include_graphics(output_ma_file)

### Plot the counts

a <- plotCounts(deseq2Data, gene=res@rownames[1], intgroup="tissueType", returnData=TRUE)
b <- plotCounts(deseq2Data, gene=res@rownames[2], intgroup="tissueType", returnData=TRUE)
c <- plotCounts(deseq2Data, gene=res@rownames[3], intgroup="tissueType", returnData=TRUE)
d <- plotCounts(deseq2Data, gene=res@rownames[4], intgroup="tissueType", returnData=TRUE)
e <- plotCounts(deseq2Data, gene=res@rownames[5], intgroup="tissueType", returnData=TRUE)
f <- plotCounts(deseq2Data, gene=res@rownames[6], intgroup="tissueType", returnData=TRUE)


g <- ggplot(a, aes(x=tissueType, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

h <- ggplot(b, aes(x=tissueType, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

i <- ggplot(c, aes(x=tissueType, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

j <- ggplot(d, aes(x=tissueType, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

k <- ggplot(e, aes(x=tissueType, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

l <- ggplot(f, aes(x=tissueType, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

par(mfrow=c(2,3))
g
h
i
j
k
l

dir.create("../results/plots_by_gene/", recursive=TRUE)
ggsave(file="../results/plots_by_gene/ENSG00000080910.svg", plot=g, width=10, height=8)
ggsave(file="../results/plots_by_gene/ENSG00000142959.svg", plot=h, width=10, height=8)
ggsave(file="../results/plots_by_gene/ENSG00000196611.svg", plot=i, width=10, height=8)
ggsave(file="../results/plots_by_gene/ENSG00000122641.svg", plot=j, width=10, height=8)
ggsave(file="../results/plots_by_gene/ENSG00000168748.svg", plot=k, width=10, height=8)
ggsave(file="../results/plots_by_gene/ENSG00000167767.svg", plot=l, width=10, height=8)


### Volcano Plot

#reset par
par(mfrow=c(1,1))

output_volcano_file <- '../results/volcano_plot.png'

png(output_volcano_file)
# Make a basic volcano plot
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))

# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
dev.off()
include_graphics(output_volcano_file)
ggsave(file="../results/PCA.svg", plot=PCA, width=10, height=10)


### PCA

#First we need to transform the raw count data
#vst function will perform variance stabilizing transformation

vsdata <- vst(deseq2Data, blind=FALSE)
PCA <- plotPCA(vsdata, intgroup="tissueType") + scale_color_manual(values=c("red","forestgreen","blue"))
PCA
ggsave(file="../results/PCA.svg", plot=PCA, width=10, height=10)
