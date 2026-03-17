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
    conda:
        "../envs/plots.yaml"
    shell:
        "PYTHONPATH=. python workflow/scripts/plots.py --type bar-all > {log} 2>&1"
