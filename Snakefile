# Snakefile — GPML pipeline entry point
#
# Usage:
#   snakemake --use-conda --cores <N>
#
# For HPC (SLURM), add a cluster profile or use:
#   snakemake --use-conda --cores <N> --executor slurm

configfile: "config/config.yaml"

# ---------------------------------------------------------------------------
# Resource policy for optional distributed SLURM execution
# ---------------------------------------------------------------------------
# These values affect scheduler requests only.  They do not change the
# workflow DAG, analysis parameters, commands, or output paths.  Cluster
# partition names are deliberately supplied at launch time by the ignored
# local controller launcher, so normal local Snakemake runs remain portable.
RESOURCE_POLICY = config.get("resource_policy") or {}
RESOURCE_RULE_OVERRIDES = RESOURCE_POLICY.get("rules") or {}
RESOURCE_MEMORY_MULTIPLIER = float(RESOURCE_POLICY.get("memory_retry_multiplier", 1.5))
RESOURCE_RUNTIME_MULTIPLIER = float(RESOURCE_POLICY.get("runtime_retry_multiplier", 1.5))
RESOURCE_MAX_MEMORY_MB = int(RESOURCE_POLICY.get("max_memory_mb", 128000))
RESOURCE_MAX_RUNTIME_MINUTES = int(RESOURCE_POLICY.get("max_runtime_minutes", 4320))
SLURM_SHORT_PARTITION = config.get("slurm_short_partition")
SLURM_LONG_PARTITION = config.get("slurm_long_partition")
SLURM_SHORT_MAX_MINUTES = int(config.get("slurm_short_max_minutes", 720))


def _rule_resource(rule_name, key, default):
    return RESOURCE_RULE_OVERRIDES.get(rule_name, {}).get(key, default)


def _rule_threads(rule_name, default):
    return int(_rule_resource(rule_name, "threads", default))


def _scaled_resource(base, attempt, multiplier, maximum):
    value = int(round(int(base) * (float(multiplier) ** (int(attempt) - 1))))
    return min(value, int(maximum))


def _rule_mem_mb(rule_name, default, attempt):
    base = _rule_resource(rule_name, "mem_mb", default)
    return _scaled_resource(base, attempt, RESOURCE_MEMORY_MULTIPLIER, RESOURCE_MAX_MEMORY_MB)


def _rule_runtime(rule_name, default, attempt):
    base = _rule_resource(rule_name, "runtime", default)
    return _scaled_resource(base, attempt, RESOURCE_RUNTIME_MULTIPLIER, RESOURCE_MAX_RUNTIME_MINUTES)


def _rule_slurm_partition(rule_name, default_runtime, attempt):
    runtime = _rule_runtime(rule_name, default_runtime, attempt)
    if runtime <= SLURM_SHORT_MAX_MINUTES:
        return SLURM_SHORT_PARTITION
    return SLURM_LONG_PARTITION or SLURM_SHORT_PARTITION

# ---------------------------------------------------------------------------
# Include rule modules
# ---------------------------------------------------------------------------

include: "workflow/rules/preprocess.smk"
include: "workflow/rules/baseline.smk"
include: "workflow/rules/correlation.smk"
include: "workflow/rules/plots.smk"
include: "workflow/rules/tables.smk"

# ---------------------------------------------------------------------------
# Convenience variables
# ---------------------------------------------------------------------------

PICKLED_DIR = config["pickled_dir"]
PLOTS_DIR   = config["plots_dir"]
TABLES_DIR  = config["tables_dir"]
FMT         = config["plot_format"]

# Map config plot names to their output file paths
PLOT_OUTPUT_MAP = {
    "pie_chart":  f"{PLOTS_DIR}/pie_plot.{FMT}",
    "bar_strict": f"{PLOTS_DIR}/bar_strict.{FMT}",
    "bar_all":    f"{PLOTS_DIR}/bar_all.{FMT}",
    "mutpred_comparison": f"{PLOTS_DIR}/MutPred.{FMT}",
    "deogen2_comparison": f"{PLOTS_DIR}/DEOGEN2.{FMT}",
    "clinpred_comparison": f"{PLOTS_DIR}/ClinPred.{FMT}",
    "primateai_comparison": f"{PLOTS_DIR}/PrimateAI.{FMT}",
    "fathmm_comparison": f"{PLOTS_DIR}/FATHMM.{FMT}",
    "mutationtaster_comparison": f"{PLOTS_DIR}/MutationTaster.{FMT}",
    "mutpred_strict": f"{PLOTS_DIR}/MutPred_strict_cor.{FMT}",
    "deogen2_strict": f"{PLOTS_DIR}/DEOGEN2_strict_cor.{FMT}",
    "clinpred_strict": f"{PLOTS_DIR}/ClinPred_strict_cor.{FMT}",
    "primateai_strict": f"{PLOTS_DIR}/PrimateAI_strict_cor.{FMT}",
    "fathmm_strict": f"{PLOTS_DIR}/FATHMM_strict_cor.{FMT}",
    "mutationtaster_strict": f"{PLOTS_DIR}/MutationTaster_strict_cor.{FMT}",
    "all_tools_all_exclude": f"{PLOTS_DIR}/all_tools_all_exclude.{FMT}",
    "all_tools_non_exclude": f"{PLOTS_DIR}/all_tools_non_exclude.{FMT}",
    "normal_corr": f"{PLOTS_DIR}/normal_corr.{FMT}",
    "mean_bar": f"{PLOTS_DIR}/mean_bar.{FMT}",
}

# ---------------------------------------------------------------------------
# Rule all — defines the final outputs that drive the whole pipeline
# ---------------------------------------------------------------------------

rule all:
    input:
        f"{TABLES_DIR}/{config['human_protein_table']}.csv",
        [PLOT_OUTPUT_MAP[p] for p in config["enabled_plots"]],


# ---------------------------------------------------------------------------
# Create output directories before any rules run
# ---------------------------------------------------------------------------

onstart:
    shell(
        "mkdir -p logs "
        "{config[pickled_dir]} "
        "{config[plots_dir]} "
        "{config[tables_dir]}"
    )
