# Snakefile — GPML pipeline entry point
#
# Usage:
#   snakemake --use-conda --cores <N>
#
# For HPC (SLURM), add a cluster profile or use:
#   snakemake --use-conda --cores <N> --executor slurm

configfile: "config/config.yaml"

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

# ---------------------------------------------------------------------------
# Rule all — defines the final outputs that drive the whole pipeline
# ---------------------------------------------------------------------------

rule all:
    input:
        f"{TABLES_DIR}/{config['human_protein_table']}.csv",
        f"{PLOTS_DIR}/pie_plot.{FMT}",
        f"{PLOTS_DIR}/MutPredDEOGEN2ClinPredPrimateAIall_exclude.{FMT}",


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
