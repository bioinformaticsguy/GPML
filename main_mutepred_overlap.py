from pathlib import Path
from main import MAVE_GS_FILE_PATH
from src.mutepred_overlaps import MutepredDataPreprocessor
from src.utils import COLUMN_NAMES_OF_MAVE_GS_DATAFRAME_LIST
from src.dataframe_preprocessor import MAVE_GOLD_STANDARD

TRAINING_DATA_FILE_PATH = Path("Data/mutepred_training_data/wo_exclusive_hgmd_mp2_training_data_MavedbData.csv")
AMINO_ACID_SEQUENCE_SUBSTITUTIONS_COLUMN_NAME = 'Amino_acid_substitutions'
AMINO_ACID_SEQUENCE_COLUMN_NAME = 'Prot_sequence'

MUTEPRED_DATAFRAME = MutepredDataPreprocessor.get_mutepred_df(TRAINING_DATA_FILE_PATH)
MAVE_GS_DATAFRAME = MAVE_GOLD_STANDARD.get_dataframe_for_mave_gs_data(MAVE_GS_FILE_PATH,
                                                   column_names=COLUMN_NAMES_OF_MAVE_GS_DATAFRAME_LIST)



filtered_data = MutepredDataPreprocessor. \
    filter_data_with_common_sequences(dataframe_to_filter=MAVE_GS_DATAFRAME,
                                      dataframe_with_common_values=MUTEPRED_DATAFRAME,
                                      df_sequence_column_name=AMINO_ACID_SEQUENCE_COLUMN_NAME)


disjunction_df = MutepredDataPreprocessor. \
                get_disjunction(superset_df=MAVE_GS_DATAFRAME,
                                subset_df=filtered_data,
                                on_column="protein_name")




print(1)
print(1)
