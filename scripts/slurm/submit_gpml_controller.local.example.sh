#!/bin/bash
# Copy to submit_gpml_controller.local.sh and edit it for your cluster.
# The local copy is git-ignored so no site-private settings are committed.
set -euo pipefail

export MINIFORGE_PATH="/path/to/miniforge"
export SNAKEMAKE_ENV="snakemake"
export SNAKEMAKE_CONDA_PREFIX="/path/to/shared-snakemake-conda"
export GPML_SHORT_PARTITION="short"
export GPML_LONG_PARTITION="long"
export GPML_SHORT_MAX_MINUTES="720"
export GPML_SLURM_ACCOUNT=""
export GPML_SLURM_NO_ACCOUNT="0"
export GPML_EXCLUDE_NODES=""
export GPML_MAX_JOBS="10"
export GPML_MAX_THREADS_PER_JOB="4"
export GPML_UNLOCK="0"

mkdir -p logs
sbatch --partition="$GPML_SHORT_PARTITION" --time=1-00:00:00 --nodes=1 \
    --cpus-per-task=1 --mem=4GB --job-name=gpml-controller \
    --output=logs/%j_%u_%N_gpml_controller.out \
    --error=logs/%j_%u_%N_gpml_controller.err \
    scripts/slurm/run_gpml_controller.sbatch "$@"
