#!/bin/bash
# Validate GPML's SLURM executor setup and perform a no-submit DAG dry run.
set -euo pipefail

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
    echo "Usage: bash scripts/slurm/check_slurm_profile.sh [targets or Snakemake args]"
    exit 0
fi

PROFILE="${GPML_WORKFLOW_PROFILE:-profiles/slurm}"
[[ -f "Snakefile" ]] || { echo "ERROR: run this command from the GPML repository root." >&2; exit 1; }
[[ -f "$PROFILE/config.yaml" ]] || { echo "ERROR: profile not found: $PROFILE/config.yaml" >&2; exit 1; }

for command in snakemake python sbatch sacct; do
    command -v "$command" >/dev/null 2>&1 || { echo "ERROR: required command not found: $command" >&2; exit 1; }
done
python -c 'import snakemake_executor_plugin_slurm' >/dev/null 2>&1 || {
    echo "ERROR: snakemake-executor-plugin-slurm is not installed in this environment." >&2
    exit 1
}

CONDA_PREFIX_ARGS=()
[[ -z "${SNAKEMAKE_CONDA_PREFIX:-}" ]] || CONDA_PREFIX_ARGS=(--conda-prefix "$SNAKEMAKE_CONDA_PREFIX")
CONFIG_ARGS=()
if [[ -n "${GPML_SHORT_PARTITION:-}" ]]; then
    [[ -n "${GPML_LONG_PARTITION:-}" ]] || {
        echo "ERROR: GPML_LONG_PARTITION is required with GPML_SHORT_PARTITION." >&2
        exit 2
    }
    CONFIG_ARGS=(--config "slurm_short_partition=${GPML_SHORT_PARTITION}" "slurm_long_partition=${GPML_LONG_PARTITION}" "slurm_short_max_minutes=${GPML_SHORT_MAX_MINUTES:-720}")
fi

snakemake --snakefile Snakefile --configfile config/config.yaml "${CONFIG_ARGS[@]}" \
    --workflow-profile "$PROFILE" --jobs 1 --cores 16 --dry-run \
    "${CONDA_PREFIX_ARGS[@]}" "$@"
echo "SLURM profile preflight passed; no jobs were submitted."
