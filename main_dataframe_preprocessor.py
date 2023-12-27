from pathlib import Path

import pandas
import pandas as pd

from src.dataframe_preprocessor import MaveGoldStandard, MutepredTrainingProcessor, dbNSFPProcessor
from main import MAVE_GS_FILE_PATH
from src.utils import COLUMN_NAMES_OF_MAVE_GS_DATAFRAME_LIST
import copy

## Path and Strings
AMINO_ACID_SEQUENCE_COLUMN_NAME = 'Prot_sequence'
SNP_COLUMN_NAME = "HGVSp_ANNOVAR"
MUTEPRED_SCORE_COLUMN_NAME = "MutPred_score"
MAVE_SCORE_COLUMN_NAME = "SNP_dict"
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



def create_mave_tool_scores_dataframe(mave_scores_dict: dict, tool_scores_dict: dict, mave_score_column_name: str,
                                      tool_score_column_name: str) -> pandas.DataFrame:
    """
    Create a DataFrame with MAVE and tool scores.

    Parameters:
    - mave_scores_dict (dict): Dictionary containing MAVE scores.
    - tool_scores_dict (dict): Dictionary containing tool scores.
    - mave_score_column_name (str): Name of the column for MAVE scores.
    - tool_score_column_name (str): Name of the column for tool scores.

    Returns:
    - pd.DataFrame: DataFrame with 'SNPs', 'MAVE_Scores', and 'Tool_Scores' columns.
    """

    # Create a set of all keys from both dictionaries
    all_keys = set(mave_scores_dict.keys()).union(tool_scores_dict.keys())

    # Create a DataFrame with None as the default value
    df = pd.DataFrame({key: [mave_scores_dict.get(key), tool_scores_dict.get(key)] for key in all_keys},
                      index=[mave_score_column_name, tool_score_column_name]).T.reset_index()

    # Rename the columns
    df.columns = ['SNPs', mave_score_column_name, tool_score_column_name]

    return df

your_column_values = MAVE_GS_DATAFRAME['protein_name'].tolist()

protein_name = your_column_values[1]

mave_score_dict = MAVE_GS_DATAFRAME.loc[MAVE_GS_DATAFRAME['protein_name'] == protein_name, MAVE_SCORE_COLUMN_NAME].values[0]
tool_score_dict = MAVE_GS_DATAFRAME.loc[MAVE_GS_DATAFRAME['protein_name'] == protein_name, MUTEPRED_SCORE_COLUMN_NAME].values[0]


df = create_mave_tool_scores_dataframe(mave_score_dict,
                                       tool_score_dict,
                                       mave_score_column_name=MAVE_SCORE_COLUMN_NAME,
                                       tool_score_column_name=MUTEPRED_SCORE_COLUMN_NAME)

print("pause")
