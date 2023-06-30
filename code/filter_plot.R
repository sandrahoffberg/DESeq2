### Pre-filtering of data to reduce dataset size

# See what affect this filter will have.
dim(deseq2Data)
dim(deseq2Data[rowSums(counts(deseq2Data)) > filter, ])

# Perform pre-filtering of the genes that have less than 5 reads in all samples
deseq2Data <- deseq2Data[rowSums(counts(deseq2Data)) > filter, ]


# Run the DESEQ function
deseq2Data <- DESeq(deseq2Data)

# Extract differential expression results
deseq2Results <- results(deseq2Data)#, contrast=c("condition", "A", "B"))

# View summary of results
summary(deseq2Results)

# Sort summary list by p-value
res <- deseq2Results[order(deseq2Results$padj),]


### Plot the MA plot: log ratio (M) vs an average (A)
output_ma_file <- '../results/MA_plot.png'
png(output_ma_file)

MAplot <- plotMA(deseq2Results, alpha=alpha, main=paste("MA plot with alpha of", alpha, sep=" "))
MAplot <- plotMA(deseq2Results)
dev.off()

### Plot the counts

dir.create("../results/plots_by_gene/", recursive=TRUE)

for (i in 1:plots) {
  a <- plotCounts(deseq2Data, gene=res@rownames[i], intgroup=condition, returnData=TRUE)

  b <- ggplot(a, aes(x=condition, y=count)) + 
      geom_point(position=position_jitter(w=0.1,h=0)) + 
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

  ggsave(file=paste("../results/plots_by_gene/", res@rownames[i], ".png", sep=""), plot=b, width=10, height=8)
}


### Volcano Plot

output_volcano_file <- '../results/volcano_plot.png'
png(output_volcano_file)

# Make a basic volcano plot
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))

# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
dev.off()


### PCA

#First we need to transform the raw count data
#vst function will perform variance stabilizing transformation

vsdata <- vst(deseq2Data, blind=FALSE)
PCA <- plotPCA(vsdata, intgroup=condition) + scale_color_manual(values=c("red","forestgreen","blue"))
ggsave(file="../results/PCA.png", plot=PCA, width=10, height=10)

