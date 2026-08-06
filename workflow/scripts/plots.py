"""
plots.py — Step 4

Generates one plot per invocation, selected via --type:
  pie        : species distribution pie chart
  bar-strict : bar plot for 6 tools with training SNPs excluded (strict correlations)
  bar-all    : bar plot across all 9 tools (Spearman correlations, no exclusion)
  tool-comparison : a legacy all-SNP versus training-SNP-excluded tool plot
  tool-strict : a legacy strict-correlation plot for one tool
  all-tools-all-exclude : legacy aggregate of all/excluded correlations
  all-tools-non-exclude : legacy aggregate of strict correlations
  normal-corr : legacy aggregate of unfiltered correlations
  mean-bar : legacy mean absolute correlation summary

Usage:
  python workflow/scripts/plots.py --type pie
  python workflow/scripts/plots.py --type bar-strict
  python workflow/scripts/plots.py --type bar-all
  python workflow/scripts/plots.py --type tool-comparison --tool MutPred
  python workflow/scripts/plots.py --type tool-strict --tool MutPred
"""

import argparse

import matplotlib.pyplot as plt

from src.constants import (
    MAVE_DATAFRAME_PICKLE_FILE_NAME,
    MAVE_DATAFRAME_ONLY_HUMAN_WITH_BASELINE_CORELATION_PICKLE_FILE_NAME,
    PLOTS_DIRECTORY_PATH,
    PLOT_FORMAT,
    COLUMN_NAME_OF_MAVE_GOLD_STANDARD_ID,
    COLUMN_NAME_OF_BASELINE_SCORES_DICTIONARY,
    PROTEIN_SHORT_DICTMAP,
    TOOLS_LIST,
    MUTEPRED_TOOL_NAME,
    DEOGEN_TOOL_NAME,
    CLINPRED_TOOL_NAME,
    PRIMATEAI_TOOL_NAME,
    FATHMM_TOOL_NAME,
    MUTATION_TASTER,
    SPEAR_COR_SUFFIX,
    STRICT_COR_SUFFIX,
    EXCLUDE_TRAINING_SAV_SUFFIX,
    PICKLED_DATAFRAMES_DIRECTORY_PATH,
)
from src.utils import load_dataframe
from src.plot_graphs import PlotGeneroator

TOOLS_WITH_TRAINING_DATA = [
    MUTEPRED_TOOL_NAME,
    DEOGEN_TOOL_NAME,
    CLINPRED_TOOL_NAME,
    PRIMATEAI_TOOL_NAME,
    FATHMM_TOOL_NAME,
    MUTATION_TASTER,
]

ID_COL       = COLUMN_NAME_OF_MAVE_GOLD_STANDARD_ID
BASELINE_COL = COLUMN_NAME_OF_BASELINE_SCORES_DICTIONARY + SPEAR_COR_SUFFIX


def plot_pie():
    full_df = load_dataframe(
        file_path=PICKLED_DATAFRAMES_DIRECTORY_PATH,
        file_name=MAVE_DATAFRAME_PICKLE_FILE_NAME,
    )
    PlotGeneroator.plot_pie_with_counts(full_df, column_name="species")


def _load_correlation_dataframe():
    return load_dataframe(
        file_path=PICKLED_DATAFRAMES_DIRECTORY_PATH,
        file_name=MAVE_DATAFRAME_ONLY_HUMAN_WITH_BASELINE_CORELATION_PICKLE_FILE_NAME,
    )


def _generate_bar_plot(corr_df, column_names, output_name, **plot_options):
    sorted_df = corr_df.loc[:, column_names].sort_values(by=BASELINE_COL, ascending=True)
    protein_names = sorted_df.pop(ID_COL).tolist()
    data_dict = {
        key.replace(SPEAR_COR_SUFFIX, ""): value
        for key, value in sorted_df.abs().to_dict("list").items()
    }
    PlotGeneroator.generate_bar_plot(
        protein_names,
        data_dict,
        PLOTS_DIRECTORY_PATH / f"{output_name}.{PLOT_FORMAT}",
        corr_df,
        **plot_options,
    )


def plot_bar_strict(output_name="bar_strict"):
    """
    Bar plot for the 6 tools that have training data, using strict correlations
    (training SNPs excluded). Matches legacy active plot logic.
    """
    corr_df = _load_correlation_dataframe()
    column_names = [ID_COL, BASELINE_COL] + [tool + STRICT_COR_SUFFIX for tool in TOOLS_WITH_TRAINING_DATA]
    _generate_bar_plot(
        corr_df,
        column_names,
        output_name,
        height=0.123,
        fig_height=12,
        legend_font_size="small",
        barlabel_font_size=4,
        padding=0.035,
        barlabel_flag=True,
        removed_snp_flag_value=False,
    )


def plot_bar_all():
    """
    Bar plot for all 9 tools using Spearman correlations (no training SNP exclusion).
    Matches legacy commented-out 'All Tools Plot' logic.
    """
    corr_df = _load_correlation_dataframe()
    column_names = [ID_COL, BASELINE_COL] + [tool + SPEAR_COR_SUFFIX for tool in TOOLS_LIST]
    _generate_bar_plot(
        corr_df,
        column_names,
        "bar_all",
        height=0.13,
        fig_height=8,
        barlabel_font_size=4,
        legend_font_size="large",
    )


def plot_legacy_tool_comparison(tool_name):
    """Reproduce the legacy per-tool all-SNP versus excluded-training-SNP plot."""
    corr_df = _load_correlation_dataframe()
    all_col = tool_name + SPEAR_COR_SUFFIX
    excluded_col = all_col + EXCLUDE_TRAINING_SAV_SUFFIX
    column_names = [ID_COL, BASELINE_COL, all_col, excluded_col]
    _generate_bar_plot(corr_df, column_names, tool_name)


def plot_legacy_tool_strict(tool_name):
    """Reproduce the legacy PSSM-baseline versus strict-correlation plot."""
    corr_df = _load_correlation_dataframe()
    column_names = [ID_COL, BASELINE_COL, tool_name + STRICT_COR_SUFFIX]
    _generate_bar_plot(
        corr_df,
        column_names,
        f"{tool_name}_strict_cor",
        removed_snp_flag_value=False,
    )


def plot_all_tools_all_exclude():
    """Reproduce the legacy six-tool aggregate, including excluded SNP scores."""
    corr_df = _load_correlation_dataframe()
    tool_columns = [
        column
        for tool in TOOLS_WITH_TRAINING_DATA
        for column in (tool + SPEAR_COR_SUFFIX + EXCLUDE_TRAINING_SAV_SUFFIX, tool + SPEAR_COR_SUFFIX)
    ]
    _generate_bar_plot(
        corr_df,
        [ID_COL, BASELINE_COL] + tool_columns,
        "all_tools_all_exclude",
        height=0.073,
        fig_height=12,
        legend_font_size="small",
        barlabel_font_size=4,
        padding=0.035,
        barlabel_flag=True,
        removed_snp_flag_value=False,
    )


def plot_all_tools_non_exclude():
    """Reproduce the legacy aggregate strict-correlation figure."""
    plot_bar_strict(output_name="all_tools_non_exclude")


def plot_normal_corr():
    """Reproduce the legacy aggregate of unfiltered correlations for six tools."""
    corr_df = _load_correlation_dataframe()
    column_names = [ID_COL, BASELINE_COL] + [tool + SPEAR_COR_SUFFIX for tool in TOOLS_WITH_TRAINING_DATA]
    _generate_bar_plot(
        corr_df,
        column_names,
        "normal_corr",
        height=0.13,
        fig_height=12,
        legend_font_size="large",
        barlabel_font_size=8,
        barlabel_flag=True,
        removed_snp_flag_value=False,
    )


def plot_mean_bar():
    """Reproduce the legacy mean absolute Spearman-correlation summary."""
    corr_df = _load_correlation_dataframe()
    columns = [BASELINE_COL] + [tool + SPEAR_COR_SUFFIX for tool in TOOLS_WITH_TRAINING_DATA]
    means = {
        column.replace(SPEAR_COR_SUFFIX, ""): corr_df[column].abs().mean()
        for column in columns
    }
    labels, values = zip(*sorted(means.items(), key=lambda item: item[1]))

    fig, ax = plt.subplots()
    bars = ax.bar(labels, values, color="skyblue")
    ax.set_title("Mean values of Spearman's correlation")
    ax.set_xlabel("Tool Names")
    ax.set_ylabel("Mean Value")
    ax.tick_params(axis="x", rotation=15)
    ax.bar_label(bars, fmt="%.3f", padding=3, fontsize=8)
    fig.tight_layout()
    fig.savefig(PLOTS_DIRECTORY_PATH / f"mean_bar.{PLOT_FORMAT}")
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--type",
        required=True,
        choices=[
            "pie", "bar-strict", "bar-all", "tool-comparison", "tool-strict",
            "all-tools-all-exclude", "all-tools-non-exclude", "normal-corr", "mean-bar",
        ],
        help="Which plot to generate",
    )
    parser.add_argument("--tool", choices=TOOLS_WITH_TRAINING_DATA, help="Tool for per-tool legacy plots")
    args = parser.parse_args()

    if args.type == "pie":
        plot_pie()
    elif args.type == "bar-strict":
        plot_bar_strict()
    elif args.type == "bar-all":
        plot_bar_all()
    elif args.type == "tool-comparison":
        if args.tool is None:
            parser.error("--tool is required with --type tool-comparison")
        plot_legacy_tool_comparison(args.tool)
    elif args.type == "tool-strict":
        if args.tool is None:
            parser.error("--tool is required with --type tool-strict")
        plot_legacy_tool_strict(args.tool)
    elif args.type == "all-tools-all-exclude":
        plot_all_tools_all_exclude()
    elif args.type == "all-tools-non-exclude":
        plot_all_tools_non_exclude()
    elif args.type == "normal-corr":
        plot_normal_corr()
    elif args.type == "mean-bar":
        plot_mean_bar()
