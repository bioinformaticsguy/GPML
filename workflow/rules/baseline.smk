# baseline.smk
#
# Step 2a — PSSM setup: filter to human proteins and initialise PSSM columns.
# Step 2b — PSSM calculation: compute LOPO (Leave-One-Protein-Out) PSSM
#            baseline scores.
#
# Both rules depend on the preprocessed DataFrame and can run in parallel
# once preprocess completes.
#
# Outputs:
#   Data/pickled_dataframes/gold_std_df_only_human_with_baseline.pkl  (step 2a)
#   Data/pickled_dataframes/pssm_base.pkl                             (step 2b)


rule pssm_setup:
    input:
        pkl = f"{config['pickled_dir']}/{config['gold_std_pkl']}",
    output:
        pkl = f"{config['pickled_dir']}/{config['human_baseline_pkl']}",
    log:
        "logs/pssm_setup.log",
    threads: _rule_threads("pssm_setup", 1)
    resources:
        mem_mb=lambda wc, attempt: _rule_mem_mb("pssm_setup", 8000, attempt),
        runtime=lambda wc, attempt: _rule_runtime("pssm_setup", 60, attempt),
        slurm_partition=lambda wc, attempt: _rule_slurm_partition("pssm_setup", 60, attempt),
    conda:
        "../envs/baseline.yaml"
    shell:
        "PYTHONPATH=. python workflow/scripts/pssm_setup.py > {log} 2>&1"


rule pssm_baseline_calc:
    input:
        pkl = f"{config['pickled_dir']}/{config['gold_std_pkl']}",
    output:
        pkl = f"{config['pickled_dir']}/{config['pssm_base_pkl']}",
    log:
        "logs/pssm_baseline_calc.log",
    threads: _rule_threads("pssm_baseline_calc", 4)
    resources:
        mem_mb=lambda wc, attempt: _rule_mem_mb("pssm_baseline_calc", 16000, attempt),
        runtime=lambda wc, attempt: _rule_runtime("pssm_baseline_calc", 480, attempt),
        slurm_partition=lambda wc, attempt: _rule_slurm_partition("pssm_baseline_calc", 480, attempt),
    conda:
        "../envs/baseline.yaml"
    shell:
        "PYTHONPATH=. python workflow/scripts/pssm_baseline_calc.py > {log} 2>&1"
