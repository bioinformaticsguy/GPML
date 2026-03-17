"""
pssm_baseline_calc.py — Step 2b

Computes PSSM baseline scores using Leave-One-Protein-Out (LOPO)
cross-validation on human proteins.

Input  : Data/pickled_dataframes/gold_std_df.pkl
Output : Data/pickled_dataframes/pssm_base.pkl
"""

from src.constants import (
    PICKLED_DATAFRAMES_DIRECTORY_PATH,
    MAVE_DATAFRAME_PICKLE_FILE_NAME,
    SPECIE_NAME_HUMAN,
    COLUMN_NAME_OF_MAVE_GOLD_STANDARD_SPECIES,
    COLUMN_NAME_OF_MAVE_GOLD_STANDARD_ID,
)
from src.baseline_calculation import LopoBaseline
from src.utils import load_dataframe, filter_dataframe_by_species, pickle_dataframe

PSSM_BASE_PICKLE_FILE_NAME = "pssm_base.pkl"


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

    df = LopoBaseline.add_pssm_column_to_df(
        df=df,
        id_column_name=COLUMN_NAME_OF_MAVE_GOLD_STANDARD_ID,
    )

    df = LopoBaseline.add_pssm_predictions_to_df(
        df, COLUMN_NAME_OF_MAVE_GOLD_STANDARD_ID
    )

    pickle_dataframe(
        dataframe=df,
        file_path=PICKLED_DATAFRAMES_DIRECTORY_PATH,
        file_name=PSSM_BASE_PICKLE_FILE_NAME,
    )
