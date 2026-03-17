from src.constants import PICKLED_DATAFRAMES_DIRECTORY_PATH, \
    COLUMN_NAME_OF_MAVE_GOLD_STANDARD_ID, COLUMN_NAME_OF_MAVE_GOLD_STANDARD_UNIPROT_ID, COLUMN_NAME_OF_MAVE_SNPS, \
    HUMAN_PROTEIN_TABLE, MAVE_DATAFRAME_ONLY_HUMAN_WITH_BASELINE_CORELATION_PICKLE_FILE_NAME
from src.tables_generator import TableGenerator
from src.utils import load_dataframe

if __name__ == '__main__':
    FULL_MAVE_DB = load_dataframe(file_path=PICKLED_DATAFRAMES_DIRECTORY_PATH,
                                  file_name=MAVE_DATAFRAME_ONLY_HUMAN_WITH_BASELINE_CORELATION_PICKLE_FILE_NAME)

    cols = [COLUMN_NAME_OF_MAVE_GOLD_STANDARD_ID,
               COLUMN_NAME_OF_MAVE_GOLD_STANDARD_UNIPROT_ID,
               COLUMN_NAME_OF_MAVE_SNPS]

    TableGenerator.generate_table(FULL_MAVE_DB, HUMAN_PROTEIN_TABLE, columns=cols)

    print('debug pausse')