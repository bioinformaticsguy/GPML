import numpy as np
import matplotlib.pyplot as plt

from src.constants import PICKLED_DATAFRAMES_DIRECTORY_PATH, \
    MAVE_DATAFRAME_ONLY_HUMAN_WITH_BASELINE_CORELATION_PICKLE_FILE_NAME, PROTEIN_SHORT_DICTMAP, MUTEPRED_TOOL_NAME, \
    COLUMN_NAME_OF_MAVE_GOLD_STANDARD_ID, COLUMN_NAME_OF_BASELINE_SCORES_DICTIONARY, PEARSON_CORELATION_SUFFIX
from src.utils import load_dataframe
from src.plot_graphs import PlotGeneroator

if __name__ == '__main__':
    LOADED_MAVE_DF = load_dataframe(file_path=PICKLED_DATAFRAMES_DIRECTORY_PATH,
                                    file_name=MAVE_DATAFRAME_ONLY_HUMAN_WITH_BASELINE_CORELATION_PICKLE_FILE_NAME)


    # PlotGeneroator.plot_correlations(LOADED_MAVE_DF)
    # PlotGeneroator.plot_pie_with_counts(LOADED_MAVE_DF, column_name="species")


    id_column = "protein_name"
    baseline_column  = "pssmBaseline_pearson_correlation"
    column1 = "MutPred_pearson_correlation"
    column2 = "MutPred_pearson_correlation_excluded_training_snps"

    column_names = [id_column, baseline_column, column1, column2]
    sorted_df = LOADED_MAVE_DF.loc[:, column_names].sort_values(by=baseline_column, ascending=True)

    protein_names = [PROTEIN_SHORT_DICTMAP[name] for name in sorted_df.pop(id_column).tolist()]
    dict_df = {k.replace("_pearson_correlation", ""): v for k, v in sorted_df.abs().to_dict('list').items()}
    PlotGeneroator.generate_bar_plot(protein_names, dict_df)

    print("Debug Pause")