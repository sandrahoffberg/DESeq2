#!/usr/bin/env bash

echo "${1}"

if [ "${1}" == "Transcript_abundance" ]; then
    echo "Running /data/transcript_abundance.R"
    Rscript transcript_abundance.R "$@"

elif [ "${1}" == "Counts_data" ]; then
    echo "Running /data/counts_data.R"
    Rscript counts_data.R "$@"

elif [ "${1}" == "HTseq_data" ]; then
    echo "Running /data/htseq_data.R"
    Rscript htseq_data.R "$@"

elif [ "${1}" == "Summarized_experiment" ]; then
    echo "Running /data/summarized_experiment.R"
    Rscript summarized_experiment.R "$@"

else 
    echo "The input data format must be specified."
    exit 1
    
fi