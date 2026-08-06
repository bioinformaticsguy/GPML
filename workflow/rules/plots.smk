# plots.smk
#
# One rule per plot type. The Snakefile's `rule all` uses config["enabled_plots"]
# to decide which of these are required, so commenting a plot out of
# config/config.yaml is all that's needed to skip it.


rule plot_pie:
    input:
        pkl = f"{config['pickled_dir']}/{config['gold_std_pkl']}",
    output:
        f"{config['plots_dir']}/pie_plot.{config['plot_format']}",
    log:
        "logs/plot_pie.log",
    threads: _rule_threads("plot_pie", 1)
    resources:
        mem_mb=lambda wc, attempt: _rule_mem_mb("plot_pie", 4000, attempt),
        runtime=lambda wc, attempt: _rule_runtime("plot_pie", 60, attempt),
        slurm_partition=lambda wc, attempt: _rule_slurm_partition("plot_pie", 60, attempt),
    conda:
        "../envs/plots.yaml"
    shell:
        "PYTHONPATH=. python workflow/scripts/plots.py --type pie > {log} 2>&1"


rule plot_bar_strict:
    input:
        correlation_pkl = f"{config['pickled_dir']}/{config['correlation_pkl']}",
        gold_std_pkl    = f"{config['pickled_dir']}/{config['gold_std_pkl']}",
    output:
        f"{config['plots_dir']}/bar_strict.{config['plot_format']}",
    log:
        "logs/plot_bar_strict.log",
    threads: _rule_threads("plot_bar_strict", 1)
    resources:
        mem_mb=lambda wc, attempt: _rule_mem_mb("plot_bar_strict", 4000, attempt),
        runtime=lambda wc, attempt: _rule_runtime("plot_bar_strict", 60, attempt),
        slurm_partition=lambda wc, attempt: _rule_slurm_partition("plot_bar_strict", 60, attempt),
    conda:
        "../envs/plots.yaml"
    shell:
        "PYTHONPATH=. python workflow/scripts/plots.py --type bar-strict > {log} 2>&1"


rule plot_bar_all:
    input:
        correlation_pkl = f"{config['pickled_dir']}/{config['correlation_pkl']}",
        gold_std_pkl    = f"{config['pickled_dir']}/{config['gold_std_pkl']}",
    output:
        f"{config['plots_dir']}/bar_all.{config['plot_format']}",
    log:
        "logs/plot_bar_all.log",
    threads: _rule_threads("plot_bar_all", 1)
    resources:
        mem_mb=lambda wc, attempt: _rule_mem_mb("plot_bar_all", 4000, attempt),
        runtime=lambda wc, attempt: _rule_runtime("plot_bar_all", 60, attempt),
        slurm_partition=lambda wc, attempt: _rule_slurm_partition("plot_bar_all", 60, attempt),
    conda:
        "../envs/plots.yaml"
    shell:
        "PYTHONPATH=. python workflow/scripts/plots.py --type bar-all > {log} 2>&1"
