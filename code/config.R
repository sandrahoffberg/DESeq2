#!/usr/bin/env Rscript
system("set -ex")

## Passing Args
args <- commandArgs(trailingOnly=TRUE)


if (length(args) == 0 | nchar(args[2])==0) {
#if (length(str(args[1]))==0) {
    raw_exp <- list.files(path = "../data", pattern = ".*raw.*", full.names = TRUE, recursive=TRUE)
} else {
    raw_exp <- args[2]
}
print(raw_exp)


if (length(args) == 0 | nchar(args[3])==0) {
#if (length(str(args[2]))==0) {
    metadata <- list.files(path = "../data", pattern = ".*etadata.*", full.names = TRUE, recursive=TRUE)
} else {
    metadata <- args[3]
}
print(metadata)


if (length(args) == 0 | nchar(args[4])==0) {
#if (length(str(args[3]))==0) {
    filter <- "5"
} else {
    filter <- args[4]
}
print(filter)


if (length(args) == 0 | nchar(args[5])==0) {
#if (length(str(args[3]))==0) {
    alpha <- "0.05"
} else {
    alpha <- args[5]
}
print(alpha)
