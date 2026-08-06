# preprocess.smk
#
# Step 1: Build the master DataFrame from MAVE gold standard data, integrate
# dbNSFP scores for all 9 prediction tools, and flag mutations that appear in
# tool training sets (DEOGEN2, ClinVar, MutPred).
#
# Input  : MAVE FASTA file, dbNSFP output CSVs, training data files
# Output : Data/pickled_dataframes/gold_std_df.pkl


rule preprocess:
    input:
        mave_gs   = config["mave_gs_file"],
        dbnsfp    = config["dbnsfp_output_dir"],
        mutpred   = config["mutpred_training"],
        deogen2   = config["deogen2_training"],
        clinvar   = config["clinvar_data"],
    output:
        pkl = f"{config['pickled_dir']}/{config['gold_std_pkl']}",
    log:
        "logs/preprocess.log",
    threads: _rule_threads("preprocess", 2)
    resources:
        mem_mb=lambda wc, attempt: _rule_mem_mb("preprocess", 32000, attempt),
        runtime=lambda wc, attempt: _rule_runtime("preprocess", 360, attempt),
        slurm_partition=lambda wc, attempt: _rule_slurm_partition("preprocess", 360, attempt),
    conda:
        "../envs/preprocess.yaml"
    shell:
        "PYTHONPATH=. python workflow/scripts/preprocess.py > {log} 2>&1"
