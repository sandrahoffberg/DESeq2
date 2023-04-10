#!/usr/bin/env Rscript

## Passing Args
args <- commandArgs(trailingOnly=TRUE)



if (length(str(args[1]))==0) {
    raw_exp <- list.files(path = "../", pattern = "raw_exp", recursive=TRUE)
} else {
    raw_exp <- args[1]
}
gene_dir <- dirname(gene_file)
print(gene_dir)


if (length(str(args[2]))==0) {
    target_list <- list.files(path = "../", pattern = "target", recursive=TRUE)
} else {
    target_list <- args[2]
}
print(target_list)


if (length(str(args[3]))==0) {
    control_genes <- list.files(path = "../", pattern = "control_genes", recursive=TRUE)
} else {
    control_genes <- args[3]
}
print(target_list)


if (length(str(args[4]))==0) {
    control_condition <- list.files(path = "../", pattern = "control_condition", recursive=TRUE)
} else {
    control_condition <- args[4]
}

base_name <- sample