#!/usr/bin/env Rscript
system("set -ex")

## Passing Args
args <- commandArgs(trailingOnly=TRUE)

if (length(args) == 0 | nchar(args[1])==0) {
    alpha = 0.05
} else {
    alpha <- args[1]
}

if (length(args) == 0 | nchar(args[2])==0) {
    colNonSig <- "darkgray"
} else {
    colNonSig <- args[2]
}

if (length(args) == 0 | nchar(args[3])==0) {
    colSig <- "blue"
} else {
    colSig <- args[3]
}
