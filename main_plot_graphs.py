from src.constants import MAVE_DATAFRAME_PICKLE_FILE_NAME, PICKLED_DATAFRAMES_DIRECTORY_PATH, PEARSON_CORELATION_SUFFIX, \
    PROTEIN_SHORT_MAPPING
from src.utils import load_dataframe, filter_dataframe_by_species
from itertools import cycle
import matplotlib.pyplot as plt

if __name__ == '__main__':
    LOADED_MAVE_DF = load_dataframe(file_path=PICKLED_DATAFRAMES_DIRECTORY_PATH,
                                        file_name=MAVE_DATAFRAME_PICKLE_FILE_NAME)

    LOADED_MAVE_DF = filter_dataframe_by_species(LOADED_MAVE_DF)

    # Sample data
    df = LOADED_MAVE_DF

    columns_list = df.columns.tolist()
    filtered_columns_list = [column for column in columns_list if column.endswith(PEARSON_CORELATION_SUFFIX)]

    legend = [value.rstrip(PEARSON_CORELATION_SUFFIX) for value in filtered_columns_list]

    x = [PROTEIN_SHORT_MAPPING[protein_name] for protein_name in df['protein_name'].tolist()]

    # Generating dynamic colors based on the number of legend entries
    num_colors = len(legend)
    color_cycle = cycle(plt.cm.get_cmap('tab10').colors)
    colors = [next(color_cycle) for _ in range(num_colors)]

    fig, ax = plt.subplots()

    for column, label, color in zip(filtered_columns_list, legend, colors):
        y = df[column].abs().tolist()
        ax.scatter(x, y, label=label, color=color)

    ax.set_xlabel('Protein Names', fontsize=14)
    ax.set_ylabel('Correlations', fontsize=14)
    ax.legend(title='Correlation Types', fontsize=12, labels=legend)
    ax.set_xticks(x)
    ax.set_xticklabels(x, rotation=45, fontsize=14)
    fig.set_size_inches(20, 20)
    fig.subplots_adjust(bottom=0.1)

    plt.show()

    print("Debug Pause")