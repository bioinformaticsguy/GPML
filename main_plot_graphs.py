from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

from src.constants import PICKLED_DATAFRAMES_DIRECTORY_PATH, \
    MAVE_DATAFRAME_ONLY_HUMAN_WITH_BASELINE_CORELATION_PICKLE_FILE_NAME, PROTEIN_SHORT_DICTMAP, MUTEPRED_TOOL_NAME, \
    COLUMN_NAME_OF_MAVE_GOLD_STANDARD_ID, COLUMN_NAME_OF_BASELINE_SCORES_DICTIONARY, SPEAR_COR_SUFFIX, \
    DEOGEN_TOOL_NAME, EXCLUDE_TRAINING_SAV_SUFFIX, TOOLS_LIST, PLOTS_DIRECTORY_PATH
from src.utils import load_dataframe
from src.plot_graphs import PlotGeneroator


def generate_plot(tool_name, file_name, exclude_training_snps=False):
    column1 = tool_name + SPEAR_COR_SUFFIX
    if exclude_training_snps:
        column1 += EXCLUDE_TRAINING_SAV_SUFFIX
    column2 = column1 + EXCLUDE_TRAINING_SAV_SUFFIX

    column_names = [id_column, baseline_column, column1, column2]
    sorted_df = LOADED_MAVE_DF.loc[:, column_names].sort_values(by=baseline_column, ascending=True)

    protein_names = [PROTEIN_SHORT_DICTMAP[name] for name in sorted_df.pop(id_column).tolist()]
    dict_df = {k.replace(SPEAR_COR_SUFFIX, ""): v for k, v in sorted_df.abs().to_dict('list').items()}
    PlotGeneroator.generate_bar_plot(protein_names, dict_df, file_name,)

if __name__ == '__main__':
    LOADED_MAVE_DF = load_dataframe(file_path=PICKLED_DATAFRAMES_DIRECTORY_PATH,
                                    file_name=MAVE_DATAFRAME_ONLY_HUMAN_WITH_BASELINE_CORELATION_PICKLE_FILE_NAME)


    # PlotGeneroator.plot_correlations(LOADED_MAVE_DF)
    # PlotGeneroator.plot_pie_with_counts(LOADED_MAVE_DF, column_name="species")


    id_column = "protein_name"
    baseline_column  = COLUMN_NAME_OF_BASELINE_SCORES_DICTIONARY + SPEAR_COR_SUFFIX

    # # Deogen Plot
    file_name = PLOTS_DIRECTORY_PATH / Path(DEOGEN_TOOL_NAME + ".svg")
    generate_plot(DEOGEN_TOOL_NAME, file_name)

    # # # MutPred Plot
    file_name = PLOTS_DIRECTORY_PATH / Path(MUTEPRED_TOOL_NAME + ".svg")
    generate_plot(MUTEPRED_TOOL_NAME, file_name)

    # # # Both Tools Plot Exclude Training SNPs
    file_name = PLOTS_DIRECTORY_PATH / Path(MUTEPRED_TOOL_NAME + DEOGEN_TOOL_NAME + ".svg")
    column1 = MUTEPRED_TOOL_NAME + SPEAR_COR_SUFFIX + EXCLUDE_TRAINING_SAV_SUFFIX
    column2 = DEOGEN_TOOL_NAME + SPEAR_COR_SUFFIX + EXCLUDE_TRAINING_SAV_SUFFIX

    column_names = [id_column, baseline_column, column1, column2]
    sorted_df = LOADED_MAVE_DF.loc[:, column_names].sort_values(by=baseline_column, ascending=True)

    protein_names = [PROTEIN_SHORT_DICTMAP[name] for name in sorted_df.pop(id_column).tolist()]
    dict_df = {k.replace(SPEAR_COR_SUFFIX, ""): v for k, v in sorted_df.abs().to_dict('list').items()}
    PlotGeneroator.generate_bar_plot(protein_names, dict_df, file_name)
    #
    # # # All Tools Plot Exclude Training SNPs
    file_name = PLOTS_DIRECTORY_PATH / Path("all_tools" + ".svg")
    column_names = [id_column, baseline_column] + [tool + SPEAR_COR_SUFFIX for tool in TOOLS_LIST]
    sorted_df = LOADED_MAVE_DF.loc[:, column_names].sort_values(by=baseline_column, ascending=True)

    protein_names = [PROTEIN_SHORT_DICTMAP[name] for name in sorted_df.pop(id_column).tolist()]
    dict_df = {k.replace(SPEAR_COR_SUFFIX, ""): v for k, v in sorted_df.abs().to_dict('list').items()}
    PlotGeneroator.generate_bar_plot(protein_names, dict_df, file_name, height = 0.14)

    print("Debug Pause")