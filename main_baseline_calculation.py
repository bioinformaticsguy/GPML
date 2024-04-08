from src.constants import PICKLED_DATAFRAMES_DIRECTORY_PATH, MAVE_DATAFRAME_PICKLE_FILE_NAME, SPECIE_NAME_HUMAN, \
    MAVE_DATAFRAME_ONLY_HUMAN_PICKLE_FILE_NAME
from src.utils import load_dataframe, filter_dataframe_by_species, pickle_dataframe

if __name__ == '__main__':
    GOLD_STD_DF_ONLY_HUMAN = load_dataframe(file_path=PICKLED_DATAFRAMES_DIRECTORY_PATH,
                                 file_name=MAVE_DATAFRAME_ONLY_HUMAN_PICKLE_FILE_NAME)


    print("Debug Pause")