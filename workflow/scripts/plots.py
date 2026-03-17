"""
plots.py — Step 4

Generates one plot per invocation, selected via --type:
  pie        : species distribution pie chart
  bar-strict : bar plot for 6 tools with training SNPs excluded (strict correlations)
  bar-all    : bar plot across all 9 tools (Spearman correlations, no exclusion)

Usage:
  python workflow/scripts/plots.py --type pie
  python workflow/scripts/plots.py --type bar-strict
  python workflow/scripts/plots.py --type bar-all
"""

import argparse

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


def plot_bar_strict():
    """
    Bar plot for the 6 tools that have training data, using strict correlations
    (training SNPs excluded). Matches legacy active plot logic.
    """
    corr_df = load_dataframe(
        file_path=PICKLED_DATAFRAMES_DIRECTORY_PATH,
        file_name=MAVE_DATAFRAME_ONLY_HUMAN_WITH_BASELINE_CORELATION_PICKLE_FILE_NAME,
    )

    column_names = [ID_COL, BASELINE_COL] + [tool + STRICT_COR_SUFFIX for tool in TOOLS_WITH_TRAINING_DATA]

    sorted_df     = corr_df.loc[:, column_names].sort_values(by=BASELINE_COL, ascending=True)
    protein_names = sorted_df.pop(ID_COL).tolist()
    dict_df       = {k.replace(SPEAR_COR_SUFFIX, ""): v for k, v in sorted_df.abs().to_dict("list").items()}

    PlotGeneroator.generate_bar_plot(
        protein_names,
        dict_df,
        PLOTS_DIRECTORY_PATH / f"bar_strict.{PLOT_FORMAT}",
        corr_df,
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
    corr_df = load_dataframe(
        file_path=PICKLED_DATAFRAMES_DIRECTORY_PATH,
        file_name=MAVE_DATAFRAME_ONLY_HUMAN_WITH_BASELINE_CORELATION_PICKLE_FILE_NAME,
    )

    column_names = [ID_COL, BASELINE_COL] + [tool + SPEAR_COR_SUFFIX for tool in TOOLS_LIST]

    sorted_df     = corr_df.loc[:, column_names].sort_values(by=BASELINE_COL, ascending=True)
    protein_names = [PROTEIN_SHORT_DICTMAP[name] for name in sorted_df.pop(ID_COL).tolist()]
    dict_df       = {k.replace(SPEAR_COR_SUFFIX, ""): v for k, v in sorted_df.abs().to_dict("list").items()}

    PlotGeneroator.generate_bar_plot(
        protein_names,
        dict_df,
        PLOTS_DIRECTORY_PATH / f"bar_all.{PLOT_FORMAT}",
        corr_df,
        height=0.13,
        fig_height=8,
        barlabel_font_size=4,
        legend_font_size="large",
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--type",
        required=True,
        choices=["pie", "bar-strict", "bar-all"],
        help="Which plot to generate",
    )
    args = parser.parse_args()

    if args.type == "pie":
        plot_pie()
    elif args.type == "bar-strict":
        plot_bar_strict()
    elif args.type == "bar-all":
        plot_bar_all()
