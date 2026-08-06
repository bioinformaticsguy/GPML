# GPML

A benchmarking pipeline for evaluating protein variant effect prediction tools against MAVE (Multiplex Assay of Variant Effect) gold standard data.

## Overview

GPML computes Spearman rank correlations between experimental MAVE measurements and scores from 9 computational prediction tools. While the MAVE gold standard includes 36 proteins spanning multiple organisms (human, bacteria, virus, yeast, E. coli), the study is limited to human proteins — dbNSFP, the source of tool scores, only contains annotations for human variants. It accounts for training set overlap to enable unbiased performance comparisons and generates publication-ready plots and summary tables.

**Prediction tools evaluated:** MutPred, REVEL, EVE, AlphaMissense, DEOGEN2, ClinPred, PrimateAI, FATHMM, MutationTaster

## Project Structure

```
GPML/
├── Snakefile                               # Pipeline entry point
├── install.sh                              # Database download and setup script
├── gpml_env.yml                            # Full monolithic conda environment
│
├── config/
│   └── config.yaml                         # Pipeline paths and parameters
│
├── src/                                    # Core library modules
│   ├── constants.py                        # Global config, column names, protein mappings
│   ├── utils.py                            # Shared utility functions
│   ├── dataframe_preprocessor.py           # MAVE data processing and tool score integration
│   ├── baseline_calculation.py             # PSSM baseline computation (LOPO cross-validation)
│   ├── pssm_baseline.py                    # PSSM matrix operations
│   ├── corelation_calculator.py            # Spearman correlation computations
│   ├── dbNSFP_preprocessor.py             # dbNSFP database query handling
│   ├── plot_graphs.py                      # Visualization generation
│   ├── tables_generator.py                 # Summary table generation
│   └── parallel_script.py                  # Parallel processing utilities
│
├── workflow/
│   ├── rules/                              # Snakemake rule modules
│   │   ├── preprocess.smk                  # Step 1: Data integration
│   │   ├── baseline.smk                    # Step 2: PSSM baseline setup + calculation
│   │   ├── correlation.smk                 # Step 3: Spearman correlations
│   │   ├── plots.smk                       # Step 4: Visualization
│   │   └── tables.smk                      # Step 5: Summary tables
│   ├── envs/                               # Per-step conda environments
│   │   ├── preprocess.yaml
│   │   ├── baseline.yaml
│   │   ├── correlation.yaml
│   │   ├── plots.yaml
│   │   └── tables.yaml
│   └── scripts/                            # Python scripts called by Snakemake rules
│       ├── main.py                         # Path definitions and data loading
│       ├── main_dataframe_preprocessor.py  # Step 1: Build master DataFrame
│       ├── main_pssm_baseline.py           # Step 2a: Add PSSM baseline columns
│       ├── main_baseline_calculation.py    # Step 2b: Compute PSSM baseline scores
│       ├── main_corelation_calculator.py   # Step 3: Calculate Spearman correlations
│       ├── main_plot_graphs.py             # Step 4: Generate comparison plots
│       ├── main_tables.py                  # Step 5: Generate summary tables
│       ├── main_dbNSFP.py                  # dbNSFP score integration and analysis
│       ├── main_mutepred_overlap.py        # MutPred training overlap analysis
│       ├── parallel_main.py                # Parallel dbNSFP processing
│       ├── all_runner.py                   # Legacy sequential runner (superseded by Snakefile)
│       └── testing_code.py                 # Development and testing scripts
│
├── Data/                                   # Input and cached data (git-ignored)
│   ├── mave_gs_data/                       # MAVE gold standard FASTA files
│   ├── dbNSFP/                             # dbNSFP database files (from install.sh)
│   ├── dbNSFP_input_dir/                   # dbNSFP query inputs
│   ├── dbNSFP_output_dir/                  # dbNSFP tool score outputs
│   ├── mutepred_scores/                    # MutPred prediction outputs
│   ├── pickled_dataframes/                 # Cached intermediate DataFrames
│   └── HandPickedMSAs/                     # Multiple sequence alignments
│
├── Plots/                                  # Generated figures (git-ignored)
└── Tables/                                 # Generated summary tables (git-ignored)
```

## Pipeline

The pipeline runs in sequential steps, each producing a cached DataFrame used by the next step:

1. **Data integration** (`main_dataframe_preprocessor.py`) — loads MAVE gold standard data, integrates dbNSFP scores for all 9 tools, and flags mutations that appear in tool training sets (DEOGEN2, ClinVar, MutPred). Outputs `gold_std_df.pkl`.

2. **PSSM baseline setup** (`main_pssm_baseline.py`) — filters down to human proteins only (non-human proteins are excluded as dbNSFP has no annotations for them) and adds PSSM columns. Outputs `pssm_base.pkl`.

3. **Baseline calculation** (`main_baseline_calculation.py`) — computes PSSM baseline scores using Leave-One-Protein-Out (LOPO) cross-validation.

4. **Correlation analysis** (`main_corelation_calculator.py`) — calculates Spearman correlations for each tool, both including and excluding training set mutations (`correlation_all` vs `correlation_strict`). Outputs `gold_std_df_only_human_with_baseline_corelation.pkl`.

5. **Visualization** (`main_plot_graphs.py`) — generates bar plots comparing tool performance, sorted by PSSM baseline correlation. Outputs to `Plots/`.

6. **Tables** (`main_tables.py`) — generates summary tables of human proteins with SNP counts and UniProt IDs. Outputs to `Tables/`.

To run the full pipeline:

```bash
python all_runner.py
```

To run dbNSFP queries in parallel:

```bash
python parallel_main.py
```

## Setup

```bash
conda env create -f gpml_env.yml
conda activate gpml_env
```

**Key dependencies:** Python 3.9, pandas 1.4, numpy 1.23, scipy 1.9, scikit-learn 1.0, matplotlib 3.5, tensorflow 2.12

## SLURM execution

GPML can submit each Snakemake rule as a separate SLURM job. The controller is
lightweight: it builds the DAG and submits child jobs; it does not process the
large data itself. This keeps the existing workflow, commands, parameters, and
outputs unchanged.

Install the SLURM executor plugin in the environment that runs Snakemake:

```bash
conda install -n snakemake -c conda-forge -c bioconda snakemake-executor-plugin-slurm
```

On the cluster, make a site-local launcher and set the placeholders for your
Conda installation and scheduler allocation. The local launcher is ignored by
Git, so account, partition, and node settings remain private.

```bash
cp scripts/slurm/submit_gpml_controller.local.example.sh \
   scripts/slurm/submit_gpml_controller.local.sh
chmod +x scripts/slurm/submit_gpml_controller.local.sh
# Edit scripts/slurm/submit_gpml_controller.local.sh for your cluster.
```

Validate the installation without submitting jobs:

```bash
conda activate snakemake
bash scripts/slurm/check_slurm_profile.sh
```

Submit the controller after the preflight succeeds:

```bash
bash scripts/slurm/submit_gpml_controller.local.sh
```

Controller logs are written to `logs/controller/`, stable per-rule logs to
`logs/rules/`, and SLURM child-job logs to `logs/slurm/`. Active jobs can be
viewed with:

```bash
bash scripts/slurm/show_gpml_jobs.sh
```
