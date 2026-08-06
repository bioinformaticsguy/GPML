#!/bin/bash
# GPML SLURM controller.  This lightweight job runs Snakemake while the SLURM
# executor plugin submits each workflow rule as a separate child job.
#SBATCH --job-name=gpml-controller
#SBATCH --output=logs/%j_%u_%N_gpml_controller.out
#SBATCH --error=logs/%j_%u_%N_gpml_controller.err
set -euo pipefail

usage() {
    cat <<'EOF'
Usage: sbatch scripts/slurm/run_gpml_controller.sbatch [options] [targets]
Options: -c/--configfile FILE, -j/--jobs N, --workflow-profile DIR, -h/--help.
All other arguments are passed to Snakemake.
Required environment: MINIFORGE_PATH
Optional environment: SNAKEMAKE_ENV, SNAKEMAKE_CONDA_PREFIX,
GPML_SHORT_PARTITION, GPML_LONG_PARTITION, GPML_SHORT_MAX_MINUTES,
GPML_SLURM_ACCOUNT, GPML_SLURM_NO_ACCOUNT, GPML_EXCLUDE_NODES,
GPML_MAX_JOBS, GPML_MAX_THREADS_PER_JOB, GPML_UNLOCK.
EOF
}

CONFIG_FILE="config/config.yaml"
WORKFLOW_PROFILE="profiles/slurm"
MAX_JOBS="${GPML_MAX_JOBS:-10}"
MAX_THREADS_PER_JOB="${GPML_MAX_THREADS_PER_JOB:-4}"
SHORT_PARTITION="${GPML_SHORT_PARTITION:-}"
LONG_PARTITION="${GPML_LONG_PARTITION:-}"
SHORT_MAX_MINUTES="${GPML_SHORT_MAX_MINUTES:-720}"
EXCLUDE_NODES="${GPML_EXCLUDE_NODES:-}"
EXTRA_SNAKEMAKE_ARGS=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        -c|--configfile) CONFIG_FILE="$2"; shift 2 ;;
        -j|--jobs) MAX_JOBS="$2"; shift 2 ;;
        --workflow-profile) WORKFLOW_PROFILE="$2"; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) EXTRA_SNAKEMAKE_ARGS+=("$1"); shift ;;
    esac
done

for value in "$MAX_JOBS" "$MAX_THREADS_PER_JOB" "$SHORT_MAX_MINUTES"; do
    [[ "$value" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: expected a positive integer, got: $value" >&2; exit 2; }
done
if [[ -n "$SHORT_PARTITION" && -z "$LONG_PARTITION" ]]; then
    echo "ERROR: GPML_LONG_PARTITION is required when GPML_SHORT_PARTITION is set." >&2
    exit 2
fi
if [[ -n "$EXCLUDE_NODES" && ! "$EXCLUDE_NODES" =~ ^[][A-Za-z0-9._,-]+$ ]]; then
    echo "ERROR: GPML_EXCLUDE_NODES is not a valid SLURM hostlist: $EXCLUDE_NODES" >&2
    exit 2
fi
if [[ -z "${MINIFORGE_PATH:-}" || ! -f "${MINIFORGE_PATH:-}/etc/profile.d/conda.sh" ]]; then
    echo "ERROR: set MINIFORGE_PATH to a Miniforge/Conda installation." >&2
    exit 1
fi
[[ -f "$CONFIG_FILE" ]] || { echo "ERROR: config file not found: $CONFIG_FILE" >&2; exit 1; }
[[ -f "$WORKFLOW_PROFILE/config.yaml" ]] || { echo "ERROR: profile not found: $WORKFLOW_PROFILE/config.yaml" >&2; exit 1; }

mkdir -p logs logs/slurm .snakemake/locks
TIMESTAMP=$(date +%Y%m%dT%H%M%S)
RUN_LOG="logs/controller_${TIMESTAMP}_${SLURM_JOB_ID:-interactive}.log"
exec > >(tee -a "$RUN_LOG") 2>&1

# shellcheck disable=SC1091
source "${MINIFORGE_PATH}/etc/profile.d/conda.sh"
conda activate "${SNAKEMAKE_ENV:-snakemake}"
python -c 'import snakemake_executor_plugin_slurm' >/dev/null 2>&1 || {
    echo "ERROR: snakemake-executor-plugin-slurm is missing." >&2
    exit 1
}

CONFIG_ARGS=()
if [[ -n "$SHORT_PARTITION" ]]; then
    CONFIG_ARGS=(--config "slurm_short_partition=${SHORT_PARTITION}" "slurm_long_partition=${LONG_PARTITION}" "slurm_short_max_minutes=${SHORT_MAX_MINUTES}")
fi
CONDA_PREFIX_ARGS=()
[[ -z "${SNAKEMAKE_CONDA_PREFIX:-}" ]] || CONDA_PREFIX_ARGS=(--conda-prefix "$SNAKEMAKE_CONDA_PREFIX")
DEFAULT_RESOURCE_ARGS=()
if [[ -n "$SHORT_PARTITION" || -n "${GPML_SLURM_ACCOUNT:-}" || -n "$EXCLUDE_NODES" ]]; then
    DEFAULT_RESOURCE_ARGS=(--default-resources mem_mb=8000 runtime=60)
    [[ -z "$SHORT_PARTITION" ]] || DEFAULT_RESOURCE_ARGS+=("slurm_partition=${SHORT_PARTITION}")
    [[ -z "${GPML_SLURM_ACCOUNT:-}" ]] || DEFAULT_RESOURCE_ARGS+=("slurm_account=${GPML_SLURM_ACCOUNT}")
    [[ -z "$EXCLUDE_NODES" ]] || DEFAULT_RESOURCE_ARGS+=("slurm_extra=--exclude=${EXCLUDE_NODES}")
fi
SLURM_PLUGIN_ARGS=()
if [[ "${GPML_SLURM_NO_ACCOUNT:-0}" == "1" ]]; then
    [[ -z "${GPML_SLURM_ACCOUNT:-}" ]] || { echo "ERROR: configure an account or GPML_SLURM_NO_ACCOUNT=1, not both." >&2; exit 2; }
    SLURM_PLUGIN_ARGS+=(--slurm-no-account)
fi

MAX_CORES=$((MAX_JOBS * MAX_THREADS_PER_JOB))
echo "GPML SLURM controller: ${SLURM_JOB_ID:-interactive}"
echo "Config: $CONFIG_FILE | profile: $WORKFLOW_PROFILE | child jobs: $MAX_JOBS | CPU cap: $MAX_CORES"
echo "Controller log: $RUN_LOG | child logs: logs/slurm"
if [[ "${GPML_UNLOCK:-0}" == "1" ]]; then
    echo "WARNING: removing Snakemake lock because GPML_UNLOCK=1."
    snakemake --snakefile Snakefile --configfile "$CONFIG_FILE" "${CONFIG_ARGS[@]}" --unlock || true
fi

snakemake --snakefile Snakefile --configfile "$CONFIG_FILE" \
    "${CONFIG_ARGS[@]}" --workflow-profile "$WORKFLOW_PROFILE" \
    --jobs "$MAX_JOBS" --cores "$MAX_CORES" --use-conda --keep-going \
    "${CONDA_PREFIX_ARGS[@]}" "${SLURM_PLUGIN_ARGS[@]}" \
    "${EXTRA_SNAKEMAKE_ARGS[@]}" "${DEFAULT_RESOURCE_ARGS[@]}"
