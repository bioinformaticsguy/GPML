from pathlib import Path

import pandas
import pandas as pd

from src.dataframe_preprocessor import MaveGoldStandard, MutepredTrainingProcessor, dbNSFPProcessor
from main import MAVE_GS_FILE_PATH
from src.utils import COLUMN_NAMES_OF_MAVE_GS_DATAFRAME_LIST, get_mave_tool_scores_dataframe
import copy

## Path and Strings
AMINO_ACID_SEQUENCE_COLUMN_NAME = 'Prot_sequence'
SNP_COLUMN_NAME = "HGVSp_ANNOVAR"
MUTEPRED_SCORE_COLUMN_NAME = "MutPred_score"
MAVE_SCORE_DICTIONARY_COLUMN_NAME = "SNP_dict"
MAVE_SCORE_COLUMN_NAME = "Mave_score"
MUTEPRED_TOOL_NAME = "MutPred"
TRAINING_DATA_FILE_PATH = Path("Data/mutepred_training_data/wo_exclusive_hgmd_mp2_training_data_MavedbData.csv")
OUTPUT_DIR_DB_NSFP = Path("Data/dbNSFP_output_dir")

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





MAVE_GS_DATAFRAME = dbNSFPProcessor.add_tool_score_column(mave_gs_dataframe=MAVE_GS_DATAFRAME,
                                                    db_nsfp_output_dir_path=OUTPUT_DIR_DB_NSFP,
                                                    tool_name=MUTEPRED_TOOL_NAME,
                                                    snp_column_name=SNP_COLUMN_NAME)




your_column_values = MAVE_GS_DATAFRAME['protein_name'].tolist()

protein_name = your_column_values[1]

mave_score_dict = MAVE_GS_DATAFRAME.loc[MAVE_GS_DATAFRAME['protein_name'] == protein_name, MAVE_SCORE_DICTIONARY_COLUMN_NAME].values[0]
tool_score_dict = MAVE_GS_DATAFRAME.loc[MAVE_GS_DATAFRAME['protein_name'] == protein_name, MUTEPRED_SCORE_COLUMN_NAME].values[0]


mave_mutepred_scores_df = get_mave_tool_scores_dataframe(mave_score_dict,
                                                         tool_score_dict,
                                                         mave_score_column_name=MAVE_SCORE_COLUMN_NAME,
                                                         tool_score_column_name=MUTEPRED_SCORE_COLUMN_NAME)





result = get_correlation_and_percentage_used(mave_mutepred_scores_df, MAVE_SCORE_COLUMN_NAME, MUTEPRED_SCORE_COLUMN_NAME)


print(f"Percentage of Rows Used: {result[0]:.2f}%")
print(f"Pearson Correlation: {result[1]:.2f}")


print("pause")


