#!/usr/bin/env Rscript
system("set -ex")

## Passing Args
args <- commandArgs(trailingOnly=TRUE)


if (args[1]=="yes") {
    run_capsule <- "yes"
} else {
    run_capsule <- "no"
}
print(run_capsule)


if (length(args) == 0 | nchar(args[2])==0) {
#if (length(str(args[1]))==0) {
    raw_exp <- list.files(path = "../data", pattern = ".*raw.exp.csv", full.names = TRUE, recursive=TRUE)
} else {
    raw_exp <- args[2]
}
print(raw_exp)


if (length(args) == 0 | nchar(args[3])==0) {
#if (length(str(args[2]))==0) {
    target_list <- list.files(path = "../data", pattern = ".*arget.*", full.names = TRUE, recursive=TRUE)
} else {
    target_list <- args[3]
}
print(target_list)


if (length(args) == 0 | nchar(args[4])==0) {
#if (length(str(args[3]))==0) {
    control_genes <- list.files(path = "../data", pattern = ".*Control.*", full.names = TRUE, recursive=TRUE)
} else {
    control_genes <- args[4]
}
print(control_genes)


if (length(args) == 0 | nchar(args[5])==0) {
#if (length(str(args[4]))==0) {
    print("Please add a Control Condition in the App Panel.")
} else {
    control_condition <- args[5]
}

# When running in CW
#control_condition <- 'DMSO_0_144'
#control_condition <- 'DMSO_0.00_264'


sample_name <- basename(raw_exp)
base_name <- substr(sample_name, 1, nchar(sample_name) - 12)
