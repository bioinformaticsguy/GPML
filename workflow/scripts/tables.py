"""
tables.py — Step 5

Generates a summary CSV table of human proteins with their UniProt IDs
and SNP counts.

Input  : Data/pickled_dataframes/gold_std_df_only_human_with_baseline_corelation.pkl
Output : Tables/human_protein_table.csv
"""

from src.constants import (
    PICKLED_DATAFRAMES_DIRECTORY_PATH,
    MAVE_DATAFRAME_ONLY_HUMAN_WITH_BASELINE_CORELATION_PICKLE_FILE_NAME,
    COLUMN_NAME_OF_MAVE_GOLD_STANDARD_ID,
    COLUMN_NAME_OF_MAVE_GOLD_STANDARD_UNIPROT_ID,
    COLUMN_NAME_OF_MAVE_SNPS,
    HUMAN_PROTEIN_TABLE,
)
from src.tables_generator import TableGenerator
from src.utils import load_dataframe


if __name__ == "__main__":

    df = load_dataframe(
        file_path=PICKLED_DATAFRAMES_DIRECTORY_PATH,
        file_name=MAVE_DATAFRAME_ONLY_HUMAN_WITH_BASELINE_CORELATION_PICKLE_FILE_NAME,
    )

    TableGenerator.generate_table(
        df,
        HUMAN_PROTEIN_TABLE,
        columns=[
            COLUMN_NAME_OF_MAVE_GOLD_STANDARD_ID,
            COLUMN_NAME_OF_MAVE_GOLD_STANDARD_UNIPROT_ID,
            COLUMN_NAME_OF_MAVE_SNPS,
        ],
    )
