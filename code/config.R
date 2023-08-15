#!/usr/bin/env Rscript
system("set -ex")

## Passing Args
args <- commandArgs(trailingOnly=TRUE)


if (length(args) == 0 | nchar(args[2])==0) {
    counts <- list.files(path = "../data", pattern = ".*raw.*", full.names = TRUE, recursive=TRUE)
} else {
    counts <- args[2]
}
print(paste("Single data file:", counts))


if (length(args) == 0 | nchar(args[3])==0) {
    transcripts <- list.files(path = "../data", pattern = "abundance", full.names = TRUE, recursive=TRUE)

} else {
    file_name <- paste(".*", args[3], ".*", sep ="")
    transcripts <- list.files(path = "../data", pattern = file_name, full.names = TRUE, recursive=TRUE)

}
print(paste("Multiple data files:", transcripts))


if (length(args) == 0 | nchar(args[4])==0) {
    transcripts_genes <- list.files(path = "../data", pattern = ".*gene.*", full.names = TRUE, recursive=TRUE)
} else {
    transcripts_genes <- args[4]
}
print(paste("File that associates genes with transcripts:", transcripts_genes))


if (length(args) == 0 | nchar(args[5])==0) {
    analysis <- "kallisto"
} else {
    analysis <- args[5]
}
print(paste("Type of transcript abundance analysis:", analysis))


if (length(args) == 0 | nchar(args[6])==0) {
    metadata <- list.files(path = "../data", pattern = ".*etadata.*", full.names = TRUE, recursive=TRUE)
} else {
    metadata <- args[6]
}
print(paste("Sample info file:", metadata))


if (length(args) == 0 | nchar(args[7])==0) {
    condition_name <- "condition"
} else {
    condition_name <- args[7]
}
print(paste("Column Name of condition in metadata file:", condition_name))


if (length(args) == 0 | nchar(args[8])==0) {
    filter <- "5"
} else {
    filter <- args[8]
}
print(paste("Minimum number of reads:", filter))



if (length(args) == 0 | nchar(args[9])==0) {
    alpha <- "0.1"
} else {
    alpha <- args[9]
}
print(paste("Alpha signifiicance level:", alpha))


if (length(args) == 0 | nchar(args[10])==0) {
    plots <- 6
} else {
    plots <- args[10]
}
print(paste("Number of genes to plot:", plots))


if (length(args) == 0 | nchar(args[11])==0) {
    outfile <- "DESeq2_results"
} else {
    outfile <- args[11]
}
print(paste("Name of output file:", outfile))