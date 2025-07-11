from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

from src.constants import PICKLED_DATAFRAMES_DIRECTORY_PATH, \
    MAVE_DATAFRAME_ONLY_HUMAN_WITH_BASELINE_CORELATION_PICKLE_FILE_NAME, PROTEIN_SHORT_DICTMAP, MUTEPRED_TOOL_NAME, \
    COLUMN_NAME_OF_MAVE_GOLD_STANDARD_ID, COLUMN_NAME_OF_BASELINE_SCORES_DICTIONARY, SPEAR_COR_SUFFIX, \
    DEOGEN_TOOL_NAME, EXCLUDE_TRAINING_SAV_SUFFIX, TOOLS_LIST, PLOTS_DIRECTORY_PATH, PLOT_FORMAT, CLINPRED_TOOL_NAME, \
    PRIMATEAI_TOOL_NAME, FATHMM_TOOL_NAME, MUTATION_TASTER, STRICT_COR_SUFFIX, MAVE_DATAFRAME_PICKLE_FILE_NAME
from src.utils import load_dataframe, extract_value
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
    FULL_MAVE_DB = load_dataframe(file_path=PICKLED_DATAFRAMES_DIRECTORY_PATH,
                                    file_name=MAVE_DATAFRAME_PICKLE_FILE_NAME)

    LOADED_MAVE_DF = load_dataframe(file_path=PICKLED_DATAFRAMES_DIRECTORY_PATH,
                                    file_name=MAVE_DATAFRAME_ONLY_HUMAN_WITH_BASELINE_CORELATION_PICKLE_FILE_NAME)

    AAAAA = extract_value(df=LOADED_MAVE_DF,
                          id_column=COLUMN_NAME_OF_MAVE_GOLD_STANDARD_ID,
                          row_value="CBS_urn:mavedb:00000005-a",
                          target_col_name="DEOGEN2_training_savs")

    # PlotGeneroator.plot_correlations(LOADED_MAVE_DF)
    PlotGeneroator.plot_pie_with_counts(FULL_MAVE_DB, column_name="species")


    id_column = "protein_name"
    baseline_column  = COLUMN_NAME_OF_BASELINE_SCORES_DICTIONARY + SPEAR_COR_SUFFIX

    # # # Deogen Plot
    # file_name = PLOTS_DIRECTORY_PATH / Path(DEOGEN_TOOL_NAME + "." + PLOT_FORMAT)
    # generate_plot(DEOGEN_TOOL_NAME, file_name)
    #
    # # # # MutPred Plot
    # file_name = PLOTS_DIRECTORY_PATH / Path(MUTEPRED_TOOL_NAME + "." + PLOT_FORMAT)
    # generate_plot(MUTEPRED_TOOL_NAME, file_name)
    #
    # # # # Both Tools Plot Exclude Training SNPs
    # file_name = PLOTS_DIRECTORY_PATH / Path(MUTEPRED_TOOL_NAME + DEOGEN_TOOL_NAME + "." + PLOT_FORMAT)
    # column1 = MUTEPRED_TOOL_NAME + SPEAR_COR_SUFFIX + EXCLUDE_TRAINING_SAV_SUFFIX
    # column2 = DEOGEN_TOOL_NAME + SPEAR_COR_SUFFIX + EXCLUDE_TRAINING_SAV_SUFFIX
    #
    # column_names = [id_column, baseline_column, column1, column2]
    # sorted_df = LOADED_MAVE_DF.loc[:, column_names].sort_values(by=baseline_column, ascending=True)
    #
    # protein_names = [PROTEIN_SHORT_DICTMAP[name] for name in sorted_df.pop(id_column).tolist()]
    # dict_df = {k.replace(SPEAR_COR_SUFFIX, ""): v for k, v in sorted_df.abs().to_dict('list').items()}
    # PlotGeneroator.generate_bar_plot(protein_names, dict_df, file_name)
    #
    # # # All Tools Plot
    # file_name = PLOTS_DIRECTORY_PATH / Path("all_tools" + "." + PLOT_FORMAT)
    # column_names = [id_column, baseline_column] + [tool + SPEAR_COR_SUFFIX for tool in TOOLS_LIST]
    # sorted_df = LOADED_MAVE_DF.loc[:, column_names].sort_values(by=baseline_column, ascending=True)
    #
    # protein_names = [PROTEIN_SHORT_DICTMAP[name] for name in sorted_df.pop(id_column).tolist()]
    # dict_df = {k.replace(SPEAR_COR_SUFFIX, ""): v for k, v in sorted_df.abs().to_dict('list').items()}
    # PlotGeneroator.generate_bar_plot(protein_names, dict_df, file_name,
    #                                  height = 0.13 ,
    #                                  fig_height=8,
    #                                  barlabel_font_size=4,
    #                                  legend_font_size="large")


    # # Both Tools Plot Exclude Training SAVs and Normal SAVs
    tool_names_list = [MUTEPRED_TOOL_NAME,
                       DEOGEN_TOOL_NAME,
                       CLINPRED_TOOL_NAME,
                       PRIMATEAI_TOOL_NAME,
                       FATHMM_TOOL_NAME,
                       MUTATION_TASTER]

    file_name = PLOTS_DIRECTORY_PATH / Path(tool_names_list[0] +
                                            tool_names_list[1] +
                                            tool_names_list[2] +
                                            tool_names_list[3] +
                                            "all_exclude" + "." + PLOT_FORMAT)

    # column_names = [id_column, baseline_column,
    #                 tool_names_list[0] + SPEAR_COR_SUFFIX + EXCLUDE_TRAINING_SAV_SUFFIX,
    #                 tool_names_list[0] + SPEAR_COR_SUFFIX,
    #                 tool_names_list[1] + SPEAR_COR_SUFFIX + EXCLUDE_TRAINING_SAV_SUFFIX,
    #                 tool_names_list[1] + SPEAR_COR_SUFFIX,
    #                 tool_names_list[2] + SPEAR_COR_SUFFIX + EXCLUDE_TRAINING_SAV_SUFFIX,
    #                 tool_names_list[2] + SPEAR_COR_SUFFIX,
    #                 tool_names_list[3] + SPEAR_COR_SUFFIX + EXCLUDE_TRAINING_SAV_SUFFIX,
    #                 tool_names_list[3] + SPEAR_COR_SUFFIX,
    #                 tool_names_list[4] + SPEAR_COR_SUFFIX + EXCLUDE_TRAINING_SAV_SUFFIX,
    #                 tool_names_list[4] + SPEAR_COR_SUFFIX,
    #                 tool_names_list[5] + SPEAR_COR_SUFFIX + EXCLUDE_TRAINING_SAV_SUFFIX,
    #                 tool_names_list[5] + SPEAR_COR_SUFFIX,]

    column_names = [id_column, baseline_column,
                    MUTEPRED_TOOL_NAME + STRICT_COR_SUFFIX,
                    DEOGEN_TOOL_NAME + STRICT_COR_SUFFIX,
                    CLINPRED_TOOL_NAME + STRICT_COR_SUFFIX,
                    PRIMATEAI_TOOL_NAME + STRICT_COR_SUFFIX,
                    FATHMM_TOOL_NAME + STRICT_COR_SUFFIX,
                    MUTATION_TASTER + STRICT_COR_SUFFIX,]



    sorted_df = LOADED_MAVE_DF.loc[:, column_names].sort_values(by=baseline_column, ascending=True)

    # protein_names = [PROTEIN_SHORT_DICTMAP[name] for name in sorted_df.pop(id_column).tolist()]
    protein_names = sorted_df.pop(id_column).tolist()
    dict_df = {k.replace(SPEAR_COR_SUFFIX, ""): v for k, v in sorted_df.abs().to_dict('list').items()}
    # PlotGeneroator.generate_bar_plot(protein_names, dict_df, file_name, LOADED_MAVE_DF,
    #                                  height = 0.073,
    #                                  fig_height=12,
    #                                  legend_font_size="small",
    #                                  barlabel_font_size=4,
    #                                  padding=0.035,
    #                                  barlabel_flag=True,
    #                                  removed_snp_flag_value=False)

    PlotGeneroator.generate_bar_plot(protein_names, dict_df, file_name, LOADED_MAVE_DF,
                                     height = 0.123,
                                     fig_height=12,
                                     legend_font_size="small",
                                     barlabel_font_size=4,
                                     padding=0.035,
                                     barlabel_flag=True,
                                     removed_snp_flag_value=False)



    print("Debug Pause")