from pathlib import Path
from src.dataframe_preprocessor import MaveGoldStandard, MutepredTrainingProcessor, dbNSFPProcessor
from main import MAVE_GS_FILE_PATH
from src.utils import COLUMN_NAMES_OF_MAVE_GS_DATAFRAME_LIST

## Path and Strings
AMINO_ACID_SEQUENCE_COLUMN_NAME = 'Prot_sequence'
TRAINING_DATA_FILE_PATH = Path("Data/mutepred_training_data/wo_exclusive_hgmd_mp2_training_data_MavedbData.csv")
TEMP_dbNSFP_PROTEIN_PATH = Path("Data/dbNSFP_output_dir/CBS_urn_mavedb_00000005-a_output.csv")
OUTPUT_DIR_MUTEPRED = Path("Data/dbNSFP_output_dir")

## Dataframes
MAVE_GS_DATAFRAME = MaveGoldStandard.get_dataframe_for_mave_gs_data(MAVE_GS_FILE_PATH,
                                                                    column_names=COLUMN_NAMES_OF_MAVE_GS_DATAFRAME_LIST)

MUTEPRED_DATAFRAME = MutepredTrainingProcessor.get_mutepred_df(TRAINING_DATA_FILE_PATH)


FILTERED_MAVE_GOLDSTANDARD_MUTEPRED_TRAINING = MutepredTrainingProcessor. \
    filter_data_with_common_sequences(dataframe_to_filter=MAVE_GS_DATAFRAME,
                                      dataframe_with_common_values=MUTEPRED_DATAFRAME,
                                      df_sequence_column_name=AMINO_ACID_SEQUENCE_COLUMN_NAME)



MAVE_GS_DATAFRAME = MaveGoldStandard.mark_rows_present_in_subset(superset_df=MAVE_GS_DATAFRAME,
                                                                 subset_df=FILTERED_MAVE_GOLDSTANDARD_MUTEPRED_TRAINING,
                                                                 new_column_name="in_mutepred")

dataframe = dbNSFPProcessor.get_dbNSFP_df(TEMP_dbNSFP_PROTEIN_PATH)

check_list = dbNSFPProcessor.check_if_all_dbNSFP_file_have_same_names_as_mave_goldstandard(MAVE_GS_DATAFRAME, OUTPUT_DIR_MUTEPRED)

print("Break Point Stoper")