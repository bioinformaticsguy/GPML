from src.baseline_calculation import LopoBaseline
from src.constants import PICKLED_DATAFRAMES_DIRECTORY_PATH, MAVE_DATAFRAME_PICKLE_FILE_NAME, SPECIE_NAME_HUMAN, \
    COLUMN_NAME_OF_MAVE_GOLD_STANDARD_SPECIES, COLUMN_NAME_OF_MAVE_GOLD_STANDARD_ID
from src.utils import load_dataframe, filter_dataframe_by_species

if __name__ == '__main__':
    GOLD_STD_DF = load_dataframe(file_path=PICKLED_DATAFRAMES_DIRECTORY_PATH,
                                 file_name=MAVE_DATAFRAME_PICKLE_FILE_NAME)

    GOLD_STD_DF_HUMAN = filter_dataframe_by_species(df=GOLD_STD_DF,
                                                   target_species=SPECIE_NAME_HUMAN,
                                                   species_column=COLUMN_NAME_OF_MAVE_GOLD_STANDARD_SPECIES)

    GOLD_STD_DF_HUMAN = LopoBaseline.add_pssm_column_to_df(df=GOLD_STD_DF_HUMAN,
                                                           id_column_name=COLUMN_NAME_OF_MAVE_GOLD_STANDARD_ID)


    GOLD_STD_DF_HUMAN = LopoBaseline.add_pssm_predictions_to_df(GOLD_STD_DF_HUMAN, COLUMN_NAME_OF_MAVE_GOLD_STANDARD_ID)

    print("Debug Pause")