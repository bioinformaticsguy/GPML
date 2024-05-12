from src.constants import MAVE_DATAFRAME_PICKLE_FILE_NAME, PICKLED_DATAFRAMES_DIRECTORY_PATH, PEARSON_CORELATION_SUFFIX, \
    PROTEIN_SHORT_MAPPING
from src.utils import load_dataframe, filter_dataframe_by_species
from src.plot_graphs import PlotGeneroator

if __name__ == '__main__':
    LOADED_MAVE_DF = load_dataframe(file_path=PICKLED_DATAFRAMES_DIRECTORY_PATH,
                                        file_name="pssm_base.pkl")

    # ONLY_HUMAN_DATABASE = filter_dataframe_by_species(LOADED_MAVE_DF)

    # ONLY_HUMAN_DATABASE = ONLY_HUMAN_DATABASE.loc[ONLY_HUMAN_DATABASE['protein_name'] != 'CYP2C9_urn:mavedb:00000095-a']

    LOADED_MAVE_DF.fillna(0, inplace=True)

    PlotGeneroator.plot_correlations(LOADED_MAVE_DF, PROTEIN_SHORT_MAPPING, PEARSON_CORELATION_SUFFIX)


    # PlotGeneroator.plot_pie_with_counts(LOADED_MAVE_DF, "species")

    print("Debug Pause")