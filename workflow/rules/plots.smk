# plots.smk
#
# Step 4: Generate publication-ready plots comparing prediction tool
# performance, sorted by PSSM baseline correlation.
#
# Input  : gold_std_df_only_human_with_baseline_corelation.pkl, gold_std_df.pkl
# Output : Plots/pie_plot.png
#          Plots/MutPredDEOGEN2ClinPredPrimateAIall_exclude.png


rule plots:
    input:
        correlation_pkl = f"{config['pickled_dir']}/{config['correlation_pkl']}",
        gold_std_pkl    = f"{config['pickled_dir']}/{config['gold_std_pkl']}",
    output:
        pie  = f"{config['plots_dir']}/pie_plot.{config['plot_format']}",
        bars = f"{config['plots_dir']}/MutPredDEOGEN2ClinPredPrimateAIall_exclude.{config['plot_format']}",
    log:
        "logs/plots.log",
    conda:
        "../envs/plots.yaml"
    shell:
        "PYTHONPATH=. python workflow/scripts/plots.py > {log} 2>&1"
