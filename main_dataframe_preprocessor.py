from pathlib import Path
from src.constants import COLUMN_NAMES_OF_MAVE_GS_DATAFRAME_LIST, TRAINING_FLAG_SUFFIX, MUTEPRED_TOOL_NAME, \
    TOOL_SCORE_COLUMN_SUFFIX, AMINO_ACID_SEQUENCE_COLUMN_NAME, SNP_COLUMN_NAME, MAVE_DATAFRAME_PICKLE_FILE_NAME, \
    TOOLS_LIST
from src.dataframe_preprocessor import MaveGoldStandard, MutepredTrainingProcessor, dbNSFPProcessor
from main import MAVE_GS_FILE_PATH
from src.utils import pickle_dataframe

## Path and Strings
MUTEPRED_SCORE_COLUMN_NAME = MUTEPRED_TOOL_NAME + TOOL_SCORE_COLUMN_SUFFIX
TRAINING_DATA_FILE_PATH = Path("Data/mutepred_training_data/wo_exclusive_hgmd_mp2_training_data_MavedbData.csv")
OUTPUT_DIR_DB_NSFP = Path("Data/dbNSFP_output_dir")
PICKLED_DATAFRAMES_DIRECTORY_PATH = Path("Data/pickled_dataframes")

MUTEPRED_TRAINING_FLAG_COLUMN_NAME = MUTEPRED_TOOL_NAME + TRAINING_FLAG_SUFFIX

def add_data_from_list_of_tools(mave_gs_dataframe,
                                db_nsfp_output_dir_path=OUTPUT_DIR_DB_NSFP,
                                tool_list=TOOLS_LIST,
                                snp_column_name=SNP_COLUMN_NAME):
    for tool in tool_list:
        dbNSFPProcessor.add_tool_score_column(mave_gs_dataframe,
                                              db_nsfp_output_dir_path,
                                              tool,
                                              snp_column_name)

    return mave_gs_dataframe

if __name__ == '__main__':
    ## Dataframes
    MAVE_GS_DATAFRAME = MaveGoldStandard. \
                            get_dataframe_for_mave_gs_data(MAVE_GS_FILE_PATH,
                                                           column_names=COLUMN_NAMES_OF_MAVE_GS_DATAFRAME_LIST)

    MUTEPRED_DATAFRAME = MutepredTrainingProcessor.get_mutepred_df(TRAINING_DATA_FILE_PATH)

    FILTERED_MAVE_GOLDSTANDARD_MUTEPRED_TRAINING = MutepredTrainingProcessor. \
        filter_data_with_common_sequences(dataframe_to_filter=MAVE_GS_DATAFRAME,
                                          dataframe_with_common_values=MUTEPRED_DATAFRAME,
                                          df_sequence_column_name=AMINO_ACID_SEQUENCE_COLUMN_NAME)

    MAVE_GS_DATAFRAME = MaveGoldStandard. \
                        mark_rows_present_in_subset(superset_df=MAVE_GS_DATAFRAME,
                                                    subset_df=FILTERED_MAVE_GOLDSTANDARD_MUTEPRED_TRAINING,
                                                    new_column_name=MUTEPRED_TRAINING_FLAG_COLUMN_NAME)


    MAVE_GS_DATAFRAME = add_data_from_list_of_tools(MAVE_GS_DATAFRAME,
                                                    db_nsfp_output_dir_path=OUTPUT_DIR_DB_NSFP,
                                                    tool_list=TOOLS_LIST,
                                                    snp_column_name=SNP_COLUMN_NAME)

    pickle_dataframe(dataframe=MAVE_GS_DATAFRAME,
                     file_path=PICKLED_DATAFRAMES_DIRECTORY_PATH,
                     file_name=MAVE_DATAFRAME_PICKLE_FILE_NAME)

