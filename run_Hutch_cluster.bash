#!/bin/bash
#SBATCH -c 8 -t 1-0

# https://vaneyckt.io/posts/safer_bash_scripts_with_set_euxo_pipefail/
set -euo pipefail

printf "Running snakemake...\n"
snakemake \
    -j 8 \
    --software-deployment-method conda \
    --keep-going \
    --rerun-incomplete
printf "Run of snakemake complete.\n"
