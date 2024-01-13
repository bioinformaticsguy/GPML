
from src.constants import MAVE_DATAFRAME_PICKLE_FILE_NAME, TOOLS_LIST, PICKLED_DATAFRAMES_DIRECTORY_PATH, \
    MUTEPRED_TOOL_NAME
from src.utils import load_dataframe, pickle_dataframe, add_missing_columns
from src.corelation_calculator import CorelationUpdator



if __name__ == '__main__':
    LOADED_MAVE_DF = load_dataframe(file_path=PICKLED_DATAFRAMES_DIRECTORY_PATH,
                                        file_name=MAVE_DATAFRAME_PICKLE_FILE_NAME)

    LOADED_MAVE_DF = CorelationUpdator. \
        add_tool_data_for_multiple_tools(mave_goldstandard_df=LOADED_MAVE_DF,
                                         tools_names_list=TOOLS_LIST)

    mutepred_tool_bias = CorelationUpdator. \
        calculate_tool_bias(df_with_spearman_scores=LOADED_MAVE_DF,
                            tool_name=MUTEPRED_TOOL_NAME)

    pickle_dataframe(dataframe=LOADED_MAVE_DF,
                     file_path=PICKLED_DATAFRAMES_DIRECTORY_PATH,
                     file_name=MAVE_DATAFRAME_PICKLE_FILE_NAME)


    print("Debug Pause")
