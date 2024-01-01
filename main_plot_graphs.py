from src.constants import MAVE_DATAFRAME_PICKLE_FILE_NAME, PICKLED_DATAFRAMES_DIRECTORY_PATH, PEARSON_CORELATION_SUFFIX, \
    PROTEIN_SHORT_MAPPING
from src.utils import load_dataframe, filter_dataframe_by_species

import matplotlib.pyplot as plt

if __name__ == '__main__':
    LOADED_MAVE_DF = load_dataframe(file_path=PICKLED_DATAFRAMES_DIRECTORY_PATH,
                                        file_name=MAVE_DATAFRAME_PICKLE_FILE_NAME)

    LOADED_MAVE_DF = filter_dataframe_by_species(LOADED_MAVE_DF)

    import pandas as pd
    from pandas.plotting import scatter_matrix

    columns_list = LOADED_MAVE_DF.columns.tolist()
    filtered_columns_list = [column for column in columns_list if column.endswith(PEARSON_CORELATION_SUFFIX)]

    # Create a sample DataFrame
    data = {'ID_Colu': ["A", "B", "C", "D"],
            'Column2': [5, 6, 7, 8],
            'Column3': [9, 10, 11, 12]}
    df = LOADED_MAVE_DF


    columns_list = filtered_columns_list
    x = df.iloc[:, 0]  # Assuming the first column is on the x-axis


    for column in columns_list:
        y = df[column]
        plt.scatter(x, y, label=column)

    plt.xlabel('Protein Names')  # Replace with your x-axis label
    plt.ylabel('Corelations')  # Replace with your y-axis label
    plt.legend()
    plt.xticks(rotation=45)
    plt.show()

    print("Debug Pause")