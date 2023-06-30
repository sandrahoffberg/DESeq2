#!/usr/bin/env bash

echo "${1}"

if [ "${1}" == "Transcript_abundance" ]; then
    Rscript transcript_abundance_files.R "$@"

elif [ "${1}" == "Count_matrix" ]; then 
    Rscript count_matrix.R "$@"

elif [ "${1}" == "htseq-count" ]; then
    Rscript htseq-count_files.R "$@"

elif [ "${1}" == "SummarizedExperiment" ]; then
    Rscript summarized_experiment.R "$@"

else 
    echo "The input data format must be specified."
    exit 1
    
fi