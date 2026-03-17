"""
plots.py — Step 4

Generates one plot per invocation, selected via --type:
  pie        : species distribution pie chart
  bar-strict : bar plot with training SNPs excluded (strict correlations)
  bar-all    : bar plot across all SNPs (no exclusion)

Usage:
  python workflow/scripts/plots.py --type pie
  python workflow/scripts/plots.py --type bar-strict
  python workflow/scripts/plots.py --type bar-all
"""

import argparse
from pathlib import Path

from src.constants import (
    PICKLED_DIR,
    MAVE_DATAFRAME_PICKLE_FILE_NAME,
    MAVE_DATAFRAME_ONLY_HUMAN_WITH_BASELINE_CORELATION_PICKLE_FILE_NAME,
    PLOTS_DIRECTORY_PATH,
    PLOT_FORMAT,
    COLUMN_NAME_OF_MAVE_GOLD_STANDARD_ID,
    COLUMN_NAME_OF_BASELINE_SCORES_DICTIONARY,
    TOOLS_LIST,
    MUTEPRED_TOOL_NAME,
    DEOGEN_TOOL_NAME,
    CLINPRED_TOOL_NAME,
    PRIMATEAI_TOOL_NAME,
    FATHMM_TOOL_NAME,
    MUTATION_TASTER,
    SPEAR_COR_SUFFIX,
    STRICT_COR_SUFFIX,
)
from src.utils import load_dataframe
from src.plot_graphs import PlotGeneroator

# Shorthand
PICKLED_DATAFRAMES_DIRECTORY_PATH = Path("Data/pickled_dataframes")

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


def plot_bar(use_strict: bool):
    corr_df = load_dataframe(
        file_path=PICKLED_DATAFRAMES_DIRECTORY_PATH,
        file_name=MAVE_DATAFRAME_ONLY_HUMAN_WITH_BASELINE_CORELATION_PICKLE_FILE_NAME,
    )

    if use_strict:
        tools      = TOOLS_WITH_TRAINING_DATA
        cor_suffix = STRICT_COR_SUFFIX
        file_name  = PLOTS_DIRECTORY_PATH / f"bar_strict.{PLOT_FORMAT}"
    else:
        tools      = TOOLS_LIST
        cor_suffix = SPEAR_COR_SUFFIX
        file_name  = PLOTS_DIRECTORY_PATH / f"bar_all.{PLOT_FORMAT}"

    cor_cols  = [tool + cor_suffix for tool in tools]
    columns   = [ID_COL, BASELINE_COL] + cor_cols

    sorted_df     = corr_df.loc[:, columns].sort_values(by=BASELINE_COL, ascending=True)
    protein_names = sorted_df.pop(ID_COL).tolist()
    dict_df       = {k.replace(cor_suffix, ""): v for k, v in sorted_df.abs().to_dict("list").items()}

    PlotGeneroator.generate_bar_plot(
        protein_names,
        dict_df,
        file_name,
        corr_df,
        height=0.123,
        fig_height=12,
        legend_font_size="small",
        barlabel_font_size=4,
        padding=0.035,
        barlabel_flag=True,
        removed_snp_flag_value=False,
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
        plot_bar(use_strict=True)
    elif args.type == "bar-all":
        plot_bar(use_strict=False)
