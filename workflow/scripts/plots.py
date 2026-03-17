"""
plots.py — Step 4

Generates a species pie chart and a bar plot comparing strict Spearman
correlations for all tools with known training data, sorted by PSSM baseline.

Input  : Data/pickled_dataframes/gold_std_df_only_human_with_baseline_corelation.pkl
         Data/pickled_dataframes/gold_std_df.pkl
Output : Plots/pie_plot.png
         Plots/MutPredDEOGEN2ClinPredPrimateAIall_exclude.png
"""

from pathlib import Path

from src.constants import (
    PICKLED_DATAFRAMES_DIRECTORY_PATH,
    MAVE_DATAFRAME_PICKLE_FILE_NAME,
    MAVE_DATAFRAME_ONLY_HUMAN_WITH_BASELINE_CORELATION_PICKLE_FILE_NAME,
    PLOTS_DIRECTORY_PATH,
    PLOT_FORMAT,
    COLUMN_NAME_OF_MAVE_GOLD_STANDARD_ID,
    COLUMN_NAME_OF_BASELINE_SCORES_DICTIONARY,
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


if __name__ == "__main__":

    full_df = load_dataframe(
        file_path=PICKLED_DATAFRAMES_DIRECTORY_PATH,
        file_name=MAVE_DATAFRAME_PICKLE_FILE_NAME,
    )
    corr_df = load_dataframe(
        file_path=PICKLED_DATAFRAMES_DIRECTORY_PATH,
        file_name=MAVE_DATAFRAME_ONLY_HUMAN_WITH_BASELINE_CORELATION_PICKLE_FILE_NAME,
    )

    # --- Pie chart: species breakdown across the full MAVE database ---
    PlotGeneroator.plot_pie_with_counts(full_df, column_name="species")

    # --- Bar plot: strict correlations for all tools with training data ---
    strict_cols   = [tool + STRICT_COR_SUFFIX for tool in TOOLS_WITH_TRAINING_DATA]
    column_names  = [ID_COL, BASELINE_COL] + strict_cols

    sorted_df     = corr_df.loc[:, column_names].sort_values(by=BASELINE_COL, ascending=True)
    protein_names = sorted_df.pop(ID_COL).tolist()
    dict_df       = {k.replace(SPEAR_COR_SUFFIX, ""): v for k, v in sorted_df.abs().to_dict("list").items()}

    bar_plot_name = "".join(TOOLS_WITH_TRAINING_DATA[:4]) + "all_exclude." + PLOT_FORMAT
    bar_plot_path = PLOTS_DIRECTORY_PATH / Path(bar_plot_name)

    PlotGeneroator.generate_bar_plot(
        protein_names,
        dict_df,
        bar_plot_path,
        corr_df,
        height=0.123,
        fig_height=12,
        legend_font_size="small",
        barlabel_font_size=4,
        padding=0.035,
        barlabel_flag=True,
        removed_snp_flag_value=False,
    )
