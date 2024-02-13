from src.constants import MAVE_DATAFRAME_PICKLE_FILE_NAME, PICKLED_DATAFRAMES_DIRECTORY_PATH, PEARSON_CORELATION_SUFFIX, \
    PROTEIN_SHORT_MAPPING
from src.utils import load_dataframe, filter_dataframe_by_species
from src.plot_graphs import PlotGeneroator

if __name__ == '__main__':
    LOADED_MAVE_DF = load_dataframe(file_path=PICKLED_DATAFRAMES_DIRECTORY_PATH,
                                        file_name=MAVE_DATAFRAME_PICKLE_FILE_NAME)

    ONLY_HUMAN_DATABASE = filter_dataframe_by_species(LOADED_MAVE_DF)

    PlotGeneroator.plot_correlations(ONLY_HUMAN_DATABASE, PROTEIN_SHORT_MAPPING, PEARSON_CORELATION_SUFFIX)


    print("Debug Pause")