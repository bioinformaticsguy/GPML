from main_dataframe_preprocessor import PICKLED_DATAFRAMES_DIRECTORY_PATH, MAVE_DATAFRAME_PICKLE_FILE_NAME
from src.utils import load_dataframe

if __name__ == '__main__':
    LOADED_MAVE_DF = load_dataframe(file_path=PICKLED_DATAFRAMES_DIRECTORY_PATH, file_name=MAVE_DATAFRAME_PICKLE_FILE_NAME)

    print("Debug Pause")