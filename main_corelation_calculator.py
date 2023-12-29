from main_dataframe_preprocessor import PICKLED_DATAFRAMES_DIRECTORY_PATH, MAVE_DATAFRAME_PICKLE_FILE_NAME, \
    MUTEPRED_TOOL_NAME
from src.dataframe_preprocessor import COLUMN_NAME_OF_MAVE_GOLD_STANDARD_SPECIES
from src.utils import load_dataframe
from src.corelation_calculator import CorelationUpdator

NAME_OF_SPECIES_TO_FILTER = "Human"



if __name__ == '__main__':
    LOADED_MAVE_DF = load_dataframe(file_path=PICKLED_DATAFRAMES_DIRECTORY_PATH,
                                        file_name=MAVE_DATAFRAME_PICKLE_FILE_NAME)

    ONLY_HUMAN_MAVE_DF = LOADED_MAVE_DF[LOADED_MAVE_DF[COLUMN_NAME_OF_MAVE_GOLD_STANDARD_SPECIES] == NAME_OF_SPECIES_TO_FILTER]

    ONLY_HUMAN_MAVE_DF = ONLY_HUMAN_MAVE_DF.dropna()

    ONLY_HUMAN_MAVE_DF = CorelationUpdator.add_tool_corelation_and_snp_percentage_column(mave_goldstandard_df=ONLY_HUMAN_MAVE_DF,
                                                                                         tool_name=MUTEPRED_TOOL_NAME)

    overall_score = ONLY_HUMAN_MAVE_DF["MutPred_pearson_correlation"].mean()
    un_bias_score = ONLY_HUMAN_MAVE_DF[ONLY_HUMAN_MAVE_DF['in_mutepred_training'] == 0][
        'MutPred_pearson_correlation'].mean()
    bias = abs(un_bias_score - overall_score)


    print("Debug Pause")
