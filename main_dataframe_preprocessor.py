from pathlib import Path
from src.dataframe_preprocessor import MaveGoldStandard, MutepredTrainingProcessor, dbNSFPProcessor, \
    COLUMN_NAMES_OF_MAVE_GS_DATAFRAME_LIST
from main import MAVE_GS_FILE_PATH
from src.utils import pickle_dataframe

## Path and Strings
AMINO_ACID_SEQUENCE_COLUMN_NAME = 'Prot_sequence'
SNP_COLUMN_NAME = "HGVSp_ANNOVAR"
MAVE_SCORE_DICTIONARY_COLUMN_NAME = "SNP_dict"
# MAVE_SCORE_COLUMN_NAME = "Mave_score"
MUTEPRED_TOOL_NAME = "MutPred"
MUTEPRED_SCORE_COLUMN_NAME = MUTEPRED_TOOL_NAME + "_score"
TRAINING_DATA_FILE_PATH = Path("Data/mutepred_training_data/wo_exclusive_hgmd_mp2_training_data_MavedbData.csv")
OUTPUT_DIR_DB_NSFP = Path("Data/dbNSFP_output_dir")
PICKLED_DATAFRAMES_DIRECTORY_PATH = Path("Data/pickled_dataframes")
MAVE_DATAFRAME_PICKLE_FILE_NAME = "MAVE_DATAFRAME.pkl"


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
                                                    new_column_name="in_mutepred_training")


    MAVE_GS_DATAFRAME = dbNSFPProcessor.add_tool_score_column(mave_gs_dataframe=MAVE_GS_DATAFRAME,
                                                              db_nsfp_output_dir_path=OUTPUT_DIR_DB_NSFP,
                                                              tool_name=MUTEPRED_TOOL_NAME,
                                                              snp_column_name=SNP_COLUMN_NAME)

    pickle_dataframe(dataframe=MAVE_GS_DATAFRAME,
                     file_path=PICKLED_DATAFRAMES_DIRECTORY_PATH,
                     file_name=MAVE_DATAFRAME_PICKLE_FILE_NAME)

