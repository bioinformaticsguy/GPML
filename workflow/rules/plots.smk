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
        "logs/rules/plot_pie.log",
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
        "logs/rules/plot_bar_strict.log",
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
        "logs/rules/plot_bar_all.log",
    threads: _rule_threads("plot_bar_all", 1)
    resources:
        mem_mb=lambda wc, attempt: _rule_mem_mb("plot_bar_all", 4000, attempt),
        runtime=lambda wc, attempt: _rule_runtime("plot_bar_all", 60, attempt),
        slurm_partition=lambda wc, attempt: _rule_slurm_partition("plot_bar_all", 60, attempt),
    conda:
        "../envs/plots.yaml"
    shell:
        "PYTHONPATH=. python workflow/scripts/plots.py --type bar-all > {log} 2>&1"


# Legacy per-tool comparisons.  These retain the historical output filenames
# while making all six thesis-era tool plots reproducible workflow targets.
rule plot_tool_comparison:
    input:
        correlation_pkl = f"{config['pickled_dir']}/{config['correlation_pkl']}",
        gold_std_pkl    = f"{config['pickled_dir']}/{config['gold_std_pkl']}",
    output:
        f"{config['plots_dir']}/{{tool}}.{config['plot_format']}",
    log:
        "logs/rules/plot_tool_comparison_{tool}.log",
    wildcard_constraints:
        tool="MutPred|DEOGEN2|ClinPred|PrimateAI|FATHMM|MutationTaster",
    threads: _rule_threads("plot_tool_comparison", 1)
    resources:
        mem_mb=lambda wc, attempt: _rule_mem_mb("plot_tool_comparison", 4000, attempt),
        runtime=lambda wc, attempt: _rule_runtime("plot_tool_comparison", 60, attempt),
        slurm_partition=lambda wc, attempt: _rule_slurm_partition("plot_tool_comparison", 60, attempt),
    conda:
        "../envs/plots.yaml"
    shell:
        "PYTHONPATH=. python workflow/scripts/plots.py --type tool-comparison --tool {wildcards.tool} > {log} 2>&1"


rule plot_tool_strict:
    input:
        correlation_pkl = f"{config['pickled_dir']}/{config['correlation_pkl']}",
        gold_std_pkl    = f"{config['pickled_dir']}/{config['gold_std_pkl']}",
    output:
        f"{config['plots_dir']}/{{tool}}_strict_cor.{config['plot_format']}",
    log:
        "logs/rules/plot_tool_strict_{tool}.log",
    wildcard_constraints:
        tool="MutPred|DEOGEN2|ClinPred|PrimateAI|FATHMM|MutationTaster",
    threads: _rule_threads("plot_tool_strict", 1)
    resources:
        mem_mb=lambda wc, attempt: _rule_mem_mb("plot_tool_strict", 4000, attempt),
        runtime=lambda wc, attempt: _rule_runtime("plot_tool_strict", 60, attempt),
        slurm_partition=lambda wc, attempt: _rule_slurm_partition("plot_tool_strict", 60, attempt),
    conda:
        "../envs/plots.yaml"
    shell:
        "PYTHONPATH=. python workflow/scripts/plots.py --type tool-strict --tool {wildcards.tool} > {log} 2>&1"


rule plot_all_tools_all_exclude:
    input:
        correlation_pkl = f"{config['pickled_dir']}/{config['correlation_pkl']}",
        gold_std_pkl    = f"{config['pickled_dir']}/{config['gold_std_pkl']}",
    output:
        f"{config['plots_dir']}/all_tools_all_exclude.{config['plot_format']}",
    log:
        "logs/rules/plot_all_tools_all_exclude.log",
    threads: _rule_threads("plot_all_tools_all_exclude", 1)
    resources:
        mem_mb=lambda wc, attempt: _rule_mem_mb("plot_all_tools_all_exclude", 8000, attempt),
        runtime=lambda wc, attempt: _rule_runtime("plot_all_tools_all_exclude", 120, attempt),
        slurm_partition=lambda wc, attempt: _rule_slurm_partition("plot_all_tools_all_exclude", 120, attempt),
    conda:
        "../envs/plots.yaml"
    shell:
        "PYTHONPATH=. python workflow/scripts/plots.py --type all-tools-all-exclude > {log} 2>&1"


rule plot_all_tools_non_exclude:
    input:
        correlation_pkl = f"{config['pickled_dir']}/{config['correlation_pkl']}",
        gold_std_pkl    = f"{config['pickled_dir']}/{config['gold_std_pkl']}",
    output:
        f"{config['plots_dir']}/all_tools_non_exclude.{config['plot_format']}",
    log:
        "logs/rules/plot_all_tools_non_exclude.log",
    threads: _rule_threads("plot_all_tools_non_exclude", 1)
    resources:
        mem_mb=lambda wc, attempt: _rule_mem_mb("plot_all_tools_non_exclude", 4000, attempt),
        runtime=lambda wc, attempt: _rule_runtime("plot_all_tools_non_exclude", 60, attempt),
        slurm_partition=lambda wc, attempt: _rule_slurm_partition("plot_all_tools_non_exclude", 60, attempt),
    conda:
        "../envs/plots.yaml"
    shell:
        "PYTHONPATH=. python workflow/scripts/plots.py --type all-tools-non-exclude > {log} 2>&1"


rule plot_normal_corr:
    input:
        correlation_pkl = f"{config['pickled_dir']}/{config['correlation_pkl']}",
        gold_std_pkl    = f"{config['pickled_dir']}/{config['gold_std_pkl']}",
    output:
        f"{config['plots_dir']}/normal_corr.{config['plot_format']}",
    log:
        "logs/rules/plot_normal_corr.log",
    threads: _rule_threads("plot_normal_corr", 1)
    resources:
        mem_mb=lambda wc, attempt: _rule_mem_mb("plot_normal_corr", 4000, attempt),
        runtime=lambda wc, attempt: _rule_runtime("plot_normal_corr", 60, attempt),
        slurm_partition=lambda wc, attempt: _rule_slurm_partition("plot_normal_corr", 60, attempt),
    conda:
        "../envs/plots.yaml"
    shell:
        "PYTHONPATH=. python workflow/scripts/plots.py --type normal-corr > {log} 2>&1"


rule plot_mean_bar:
    input:
        correlation_pkl = f"{config['pickled_dir']}/{config['correlation_pkl']}",
    output:
        f"{config['plots_dir']}/mean_bar.{config['plot_format']}",
    log:
        "logs/rules/plot_mean_bar.log",
    threads: _rule_threads("plot_mean_bar", 1)
    resources:
        mem_mb=lambda wc, attempt: _rule_mem_mb("plot_mean_bar", 4000, attempt),
        runtime=lambda wc, attempt: _rule_runtime("plot_mean_bar", 60, attempt),
        slurm_partition=lambda wc, attempt: _rule_slurm_partition("plot_mean_bar", 60, attempt),
    conda:
        "../envs/plots.yaml"
    shell:
        "PYTHONPATH=. python workflow/scripts/plots.py --type mean-bar > {log} 2>&1"


rule plot_averaged_over_savs:
    input:
        correlation_pkl = f"{config['pickled_dir']}/{config['correlation_pkl']}",
    output:
        f"{config['plots_dir']}/averaged_over-savs.{config['plot_format']}",
    log:
        "logs/rules/plot_averaged_over_savs.log",
    threads: _rule_threads("plot_averaged_over_savs", 1)
    resources:
        mem_mb=lambda wc, attempt: _rule_mem_mb("plot_averaged_over_savs", 4000, attempt),
        runtime=lambda wc, attempt: _rule_runtime("plot_averaged_over_savs", 60, attempt),
        slurm_partition=lambda wc, attempt: _rule_slurm_partition("plot_averaged_over_savs", 60, attempt),
    conda:
        "../envs/plots.yaml"
    shell:
        "PYTHONPATH=. python workflow/scripts/plots.py --type averaged-over-savs > {log} 2>&1"
