from src.constants import MAVE_DATAFRAME_PICKLE_FILE_NAME, PICKLED_DATAFRAMES_DIRECTORY_PATH, PEARSON_CORELATION_SUFFIX, \
    PROTEIN_SHORT_MAPPING
from src.utils import load_dataframe, filter_dataframe_by_species
from itertools import cycle

import matplotlib.pyplot as plt

if __name__ == '__main__':
    LOADED_MAVE_DF = load_dataframe(file_path=PICKLED_DATAFRAMES_DIRECTORY_PATH,
                                        file_name=MAVE_DATAFRAME_PICKLE_FILE_NAME)

    LOADED_MAVE_DF = filter_dataframe_by_species(LOADED_MAVE_DF)

    columns_list = LOADED_MAVE_DF.columns.tolist()
    filtered_columns_list = [column
                             for column in columns_list
                             if column.endswith(PEARSON_CORELATION_SUFFIX)]

    legend = [value.rstrip(PEARSON_CORELATION_SUFFIX) for value in filtered_columns_list]

    df = LOADED_MAVE_DF

    columns_list = filtered_columns_list
    x = [PROTEIN_SHORT_MAPPING[key] for key in df['protein_name'].tolist()]  # Assuming the first column is on the x-axis

    # Generating dynamic colors based on the number of legend entries
    num_colors = len(legend)
    color_cycle = cycle(plt.cm.get_cmap('tab10').colors)  # You can replace 'tab10' with other colormaps
    colors = [next(color_cycle) for _ in range(num_colors)]

    for column, label, color in zip(columns_list, legend, colors):
        y = df[column].abs()  # Calculate absolute values
        plt.scatter(x, y, label=label, color=color)

    fig_width = max(6, len(x) * 0.5)  # Minimum width of 6 inches, adjust as needed
    fig_height = max(8, len(columns_list) * 0.5)  # Minimum height of 8 inches, adjust as needed




    plt.xlabel('Protein Names', fontsize=14)  # Replace with your x-axis label
    plt.ylabel('Corelations', fontsize=14)  # Replace with your y-axis label
    plt.legend(title='Correlation Types', bbox_to_anchor=(1, 1), fontsize=12,  labels=legend)  # Adjust the legend title and position
    plt.xticks(rotation=45, fontsize=14)
    plt.gcf().set_size_inches(12, 10)  # Adjust the size as needed

    # plt.gcf().set_size_inches(fig_width, fig_height)
    plt.subplots_adjust(bottom=0.1)  # Adjust the bottom margin as needed
    plt.tight_layout()

    plt.show()

    print("Debug Pause")