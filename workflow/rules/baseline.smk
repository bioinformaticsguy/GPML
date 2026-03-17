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
    conda:
        "../envs/baseline.yaml"
    shell:
        "PYTHONPATH=. python workflow/scripts/pssm_baseline_calc.py > {log} 2>&1"
