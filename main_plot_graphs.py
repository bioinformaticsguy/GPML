from src.constants import PICKLED_DATAFRAMES_DIRECTORY_PATH, \
    MAVE_DATAFRAME_ONLY_HUMAN_WITH_BASELINE_CORELATION_PICKLE_FILE_NAME
from src.utils import load_dataframe
from src.plot_graphs import PlotGeneroator

if __name__ == '__main__':
    LOADED_MAVE_DF = load_dataframe(file_path=PICKLED_DATAFRAMES_DIRECTORY_PATH,
                                    file_name=MAVE_DATAFRAME_ONLY_HUMAN_WITH_BASELINE_CORELATION_PICKLE_FILE_NAME)


    LOADED_MAVE_DF.drop("pssmBaseline_score_pearson_correlation", axis=1, inplace=True)


    PlotGeneroator.plot_correlations(LOADED_MAVE_DF)
    # PlotGeneroator.plot_pie_with_counts(LOADED_MAVE_DF, column_name="species")


    print("Debug Pause")