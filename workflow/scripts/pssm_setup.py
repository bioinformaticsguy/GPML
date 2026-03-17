"""
pssm_setup.py — Step 2a

Filters the master DataFrame to human proteins only (dbNSFP is human-only)
and initialises the PSSM baseline column structure.

Input  : Data/pickled_dataframes/gold_std_df.pkl
Output : Data/pickled_dataframes/gold_std_df_only_human_with_baseline.pkl
"""

import copy

from src.constants import (
    PICKLED_DATAFRAMES_DIRECTORY_PATH,
    MAVE_DATAFRAME_PICKLE_FILE_NAME,
    MAVE_DATAFRAME_ONLY_HUMAN_WITH_BASELINE_PICKLE_FILE_NAME,
    SPECIE_NAME_HUMAN,
    COLUMN_NAME_OF_MAVE_GOLD_STANDARD_SPECIES,
    COLUMN_NAME_OF_MAVE_GOLD_STANDARD_SAV_DICTIONARY,
    COLUMN_NAME_OF_MAVE_GOLD_STANDARD_ID,
    COLUMN_NAME_OF_BASELINE_SCORES_DICTIONARY,
)
from src.utils import load_dataframe, filter_dataframe_by_species, get_protein_list, pickle_dataframe
from src.pssm_baseline import pssmBaseline


if __name__ == "__main__":

    df = load_dataframe(
        file_path=PICKLED_DATAFRAMES_DIRECTORY_PATH,
        file_name=MAVE_DATAFRAME_PICKLE_FILE_NAME,
    )

    df = filter_dataframe_by_species(
        df=df,
        target_species=SPECIE_NAME_HUMAN,
        species_column=COLUMN_NAME_OF_MAVE_GOLD_STANDARD_SPECIES,
    )

    # Initialise PSSM column as a deep copy of the MAVE SAV dictionary
    df[COLUMN_NAME_OF_BASELINE_SCORES_DICTIONARY] = df[
        COLUMN_NAME_OF_MAVE_GOLD_STANDARD_SAV_DICTIONARY
    ].apply(lambda x: copy.deepcopy(x))

    for protein_name in get_protein_list(dataframe=df):
        pssmBaseline.update_all_values_in_snp_scores_dict(
            dataframe=df,
            id_column=COLUMN_NAME_OF_MAVE_GOLD_STANDARD_ID,
            protein_name=protein_name,
            dict_col_name=COLUMN_NAME_OF_BASELINE_SCORES_DICTIONARY,
            default_value=False,
        )

    pickle_dataframe(
        dataframe=df,
        file_path=PICKLED_DATAFRAMES_DIRECTORY_PATH,
        file_name=MAVE_DATAFRAME_ONLY_HUMAN_WITH_BASELINE_PICKLE_FILE_NAME,
    )
