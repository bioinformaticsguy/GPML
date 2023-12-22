from src.dataframe_preprocessor import MAVE_GOLD_STANDARD
from main import MAVE_GS_FILE_PATH
from src.utils import COLUMN_NAMES_OF_MAVE_GS_DATAFRAME_LIST


MAVE_GS_DATAFRAME = MAVE_GOLD_STANDARD.get_dataframe_for_mave_gs_data(MAVE_GS_FILE_PATH,
                                                   column_names=COLUMN_NAMES_OF_MAVE_GS_DATAFRAME_LIST)

ONLY_HUMAN_MAVE_GS_DATAFRAME = MAVE_GS_DATAFRAME[MAVE_GS_DATAFRAME['species'] == "Human"]