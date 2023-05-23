#!/usr/bin/env Rscript
system("set -ex")

## Passing Args
args <- commandArgs(trailingOnly=TRUE)


raw_exp <- list.files(path = "../data", pattern = ".*raw.exp.csv", full.names = TRUE, recursive=TRUE)

print(raw_exp)


target_list <- list.files(path = "../data", pattern = ".*arget.*", full.names = TRUE, recursive=TRUE)

print(target_list)


control_genes <- list.files(path = "../data", pattern = ".*Control.*", full.names = TRUE, recursive=TRUE)

print(control_genes)


if (length(args) == 0 | nchar(args[1])==0) {
    print("Please add a Control Condition in the App Panel.")
} else {
    control_condition <- args[4]
}

# When running in CW
#control_condition <- 'DMSO_0_144'
#control_condition <- 'DMSO_0.00_264'


sample_name <- basename(raw_exp)
base_name <- substr(sample_name, 1, nchar(sample_name) - 12)
