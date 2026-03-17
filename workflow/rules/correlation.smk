# correlation.smk
#
# Step 3: Calculate Spearman correlations between MAVE measurements and each
# prediction tool, both including and excluding training set mutations
# (correlation_all vs correlation_strict).
#
# Input  : gold_std_df_only_human_with_baseline.pkl, pssm_base.pkl
# Output : Data/pickled_dataframes/gold_std_df_only_human_with_baseline_corelation.pkl


rule correlation:
    input:
        human_baseline = f"{config['pickled_dir']}/{config['human_baseline_pkl']}",
        pssm_base      = f"{config['pickled_dir']}/{config['pssm_base_pkl']}",
    output:
        pkl = f"{config['pickled_dir']}/{config['correlation_pkl']}",
    log:
        "logs/correlation.log",
    conda:
        "../envs/correlation.yaml"
    shell:
        "PYTHONPATH=. python workflow/scripts/correlation.py > {log} 2>&1"
