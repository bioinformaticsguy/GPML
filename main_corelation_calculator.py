from main_dataframe_preprocessor import PICKLED_DATAFRAMES_DIRECTORY_PATH, MAVE_DATAFRAME_PICKLE_FILE_NAME, \
    MUTEPRED_TOOL_NAME
from src.utils import load_dataframe, filter_dataframe_by_species
from src.corelation_calculator import CorelationUpdator


if __name__ == '__main__':
    LOADED_MAVE_DF = load_dataframe(file_path=PICKLED_DATAFRAMES_DIRECTORY_PATH,
                                        file_name=MAVE_DATAFRAME_PICKLE_FILE_NAME)

    ONLY_HUMAN_MAVE_DF = filter_dataframe_by_species(LOADED_MAVE_DF)

    ONLY_HUMAN_MAVE_DF = CorelationUpdator.add_tool_correlation_and_snp_percentage_column(mave_goldstandard_df=ONLY_HUMAN_MAVE_DF,
                                                                                         tool_name=MUTEPRED_TOOL_NAME)

    overall_score = ONLY_HUMAN_MAVE_DF["MutPred_pearson_correlation"].mean()
    un_bias_score = ONLY_HUMAN_MAVE_DF[ONLY_HUMAN_MAVE_DF['in_mutepred_training'] == 0][
        'MutPred_pearson_correlation'].mean()
    bias = abs(un_bias_score - overall_score)

    print("Debug Pause")