#!/usr/bin/env Rscript
system("set -ex")

source("config.R", local=TRUE)

if (run_capsule == "yes") {

library(DESeq2)
    
#Load expression data. Set sample ID as row names.
raw.exp <- read.csv(raw_exp)
raw.exp <- raw.exp[!(raw.exp$Treatment == 'none'),]

#Load targets and controls.
targets <- read.csv(target_file)
controls <- read.csv(control_genes)
control.condition <- control_condition

#Grab metadata.
metadata <- raw.exp[,c("SampleID","TreatmentDosageTreatTime")]
rownames(metadata) <- metadata$SampleID
metadata$SampleID <- NULL

#Format expression data.
rownames(raw.exp) <- raw.exp$SampleID
raw.exp <- t(raw.exp[ , gsub('-','.',targets[,1])])

#Load DESeq2 dataset. Specify user-provided controls.
dds <- DESeqDataSetFromMatrix(countData = raw.exp, colData = metadata,
                              design = ~ TreatmentDosageTreatTime)
dds <- estimateSizeFactors(dds, type = "poscount",
                          controlGenes=rownames(dds) %in% controls[,1])
dds <- DESeq(dds, fitType = 'mean')

#Output results sorted by adjusted p-value.
for (condition in sort(unique(metadata$TreatmentDosageTreatTime))) {
    if (condition != control.condition) {
        res <- results(dds, contrast = c("TreatmentDosageTreatTime",
                                     condition,control.condition))
        write.csv(as.data.frame(res[order(res$padj),]),
                  file=paste("../results/",base_name,"_",condition,'.csv',sep=""))
    }
}

} else {
    dummy_file<-file(paste("../results/",base_name,"_deseq_not_run.txt", sep=""))
    writeLines("DESeq2 capsule was not run", dummy_file)
    close(dummy_file)
}