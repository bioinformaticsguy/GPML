from src.constants import PICKLED_DATAFRAMES_DIRECTORY_PATH, MAVE_DATAFRAME_PICKLE_FILE_NAME, SPECIE_NAME_HUMAN
from src.utils import load_dataframe, filter_dataframe_by_species

if __name__ == '__main__':
    LOADED_MAVE_DF = load_dataframe(file_path=PICKLED_DATAFRAMES_DIRECTORY_PATH,
                                        file_name=MAVE_DATAFRAME_PICKLE_FILE_NAME)

    ONLY_HUMAN_DF = filter_dataframe_by_species(df=LOADED_MAVE_DF,
                                                target_species=SPECIE_NAME_HUMAN)

    print("Debug Pause")