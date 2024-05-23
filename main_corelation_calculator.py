from src.constants import MAVE_DATAFRAME_PICKLE_FILE_NAME, TOOLS_LIST, PICKLED_DATAFRAMES_DIRECTORY_PATH, \
    MUTEPRED_TOOL_NAME, MAVE_DATAFRAME_ONLY_HUMAN_WITH_BASELINE_PICKLE_FILE_NAME, \
    COLUMN_NAME_OF_BASELINE_SCORES_DICTIONARY, MAVE_DATAFRAME_ONLY_HUMAN_WITH_BASELINE_CORELATION_PICKLE_FILE_NAME
from src.utils import load_dataframe, pickle_dataframe
from src.corelation_calculator import CorelationUpdator


if __name__ == '__main__':
    MAVE_GS_DATAFRAME_HUMAN_WITH_BASELINE = load_dataframe(file_path=PICKLED_DATAFRAMES_DIRECTORY_PATH,
                                                           file_name=MAVE_DATAFRAME_ONLY_HUMAN_WITH_BASELINE_PICKLE_FILE_NAME)


    MAVE_GS_DATAFRAME_HUMAN_WITH_BASELINE = CorelationUpdator. \
        add_tool_data_for_multiple_tools(MAVE_GS_DATAFRAME_HUMAN_WITH_BASELINE,
                                         tools_names_list=TOOLS_LIST)


    LOADED_MAVE_DF = CorelationUpdator. \
        add_tool_correlation_and_snp_percentage_column(mave_goldstandard_df=MAVE_GS_DATAFRAME_HUMAN_WITH_BASELINE,
                                                       tool_name=MUTEPRED_TOOL_NAME,
                                                       exclude_tool_training_snps_flag=True)

    LOADED_MAVE_DF = CorelationUpdator. \
        add_tool_correlation_and_snp_percentage_column(mave_goldstandard_df=MAVE_GS_DATAFRAME_HUMAN_WITH_BASELINE,
                                                       tool_name=COLUMN_NAME_OF_BASELINE_SCORES_DICTIONARY,)




    mean_normal = LOADED_MAVE_DF["MutPred_pearson_correlation"].mean()

    mean_no_bias = LOADED_MAVE_DF["MutPred_pearson_correlation_excluded_training_snps"].mean()


    strict = CorelationUpdator.calculate_tool_bias(df_with_spearman_scores=LOADED_MAVE_DF,
                        tool_name=MUTEPRED_TOOL_NAME)

    relax = abs(mean_no_bias - mean_normal)

    #
    # mutepred_tool_bias = CorelationUpdator. \
    #     calculate_tool_bias(df_with_spearman_scores=LOADED_MAVE_DF,
    #                         tool_name=MUTEPRED_TOOL_NAME)
    #
    pickle_dataframe(dataframe=LOADED_MAVE_DF,
                     file_path=PICKLED_DATAFRAMES_DIRECTORY_PATH,
                     file_name=MAVE_DATAFRAME_ONLY_HUMAN_WITH_BASELINE_CORELATION_PICKLE_FILE_NAME)


    print("Debug Pause")
