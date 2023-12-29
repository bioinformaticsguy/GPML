from main_dataframe_preprocessor import PICKLED_DATAFRAMES_DIRECTORY_PATH, MUTEPRED_TOOL_NAME
from src.constants import MAVE_DATAFRAME_PICKLE_FILE_NAME
from src.utils import load_dataframe, filter_dataframe_by_species, pickle_dataframe
from src.corelation_calculator import CorelationUpdator

def add_missing_columns(dataframe1, dataframe2):
    """
    Dynamically adds missing columns from dataframe2 to dataframe1 and copies values.

    Parameters:
    - dataframe1 (pd.DataFrame): Target DataFrame.
    - dataframe2 (pd.DataFrame): DataFrame with additional columns.

    Returns:
    pd.DataFrame: DataFrame1 with missing columns added and values copied.
    """
    missing_columns = set(dataframe2.columns) - set(dataframe1.columns)

    for column in missing_columns:
        dataframe1[column] = dataframe2[column]

    return dataframe1


if __name__ == '__main__':
    LOADED_MAVE_DF = load_dataframe(file_path=PICKLED_DATAFRAMES_DIRECTORY_PATH,
                                        file_name=MAVE_DATAFRAME_PICKLE_FILE_NAME)

    ONLY_HUMAN_MAVE_DF = filter_dataframe_by_species(LOADED_MAVE_DF)

    ONLY_HUMAN_MAVE_DF_WITH_CORELATIONS = CorelationUpdator.add_tool_correlation_and_snp_percentage_column(mave_goldstandard_df=ONLY_HUMAN_MAVE_DF,
                                                                                         tool_name=MUTEPRED_TOOL_NAME)

    mutepred_tool_bias = CorelationUpdator.calculate_tool_bias(df_with_spearman_scores=ONLY_HUMAN_MAVE_DF_WITH_CORELATIONS,
                                             tool_name=MUTEPRED_TOOL_NAME)

    appended_df = add_missing_columns(dataframe1=LOADED_MAVE_DF,
                                      dataframe2=ONLY_HUMAN_MAVE_DF_WITH_CORELATIONS)

    pickle_dataframe(dataframe=appended_df,
                     file_path=PICKLED_DATAFRAMES_DIRECTORY_PATH,
                     file_name=MAVE_DATAFRAME_PICKLE_FILE_NAME)


    print("Debug Pause")
