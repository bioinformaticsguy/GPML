# tables.smk
#
# Step 5: Generate summary tables of human proteins including protein name,
# UniProt ID, and SNP counts.
#
# Input  : gold_std_df_only_human_with_baseline_corelation.pkl
# Output : Tables/human_protein_table.csv


rule tables:
    input:
        pkl = f"{config['pickled_dir']}/{config['correlation_pkl']}",
    output:
        csv = f"{config['tables_dir']}/{config['human_protein_table']}.csv",
    log:
        "logs/tables.log",
    conda:
        "../envs/tables.yaml"
    shell:
        "PYTHONPATH=. python workflow/scripts/tables.py > {log} 2>&1"
