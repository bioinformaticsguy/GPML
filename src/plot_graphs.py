from itertools import cycle
from operator import itemgetter

import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import UnivariateSpline

from src.constants import PROTEIN_SHORT_DICTMAP, SPEAR_COR_SUFFIX


class PlotGeneroator:
    @staticmethod
    def plot_correlations(dataframe,
                          protein_short_mapping=PROTEIN_SHORT_DICTMAP,
                          pearson_corelation_suffix=SPEAR_COR_SUFFIX):
        dataframe.fillna(0, inplace=True)
        columns_list = dataframe.columns.tolist()
        filtered_columns_list = [column for column in columns_list if column.endswith(pearson_corelation_suffix)]

        legend = [value.rstrip(pearson_corelation_suffix) for value in filtered_columns_list]

        x = [protein_short_mapping[protein_name] for protein_name in dataframe['protein_name'].tolist()]

        # Generating dynamic colors based on the number of legend entries
        num_colors = len(legend)
        color_cycle = cycle(plt.cm.get_cmap('tab10').colors)
        colors = [next(color_cycle) for _ in range(num_colors)]

        fig, ax = plt.subplots()

        for column, label, color in zip(filtered_columns_list, legend, colors):
            y = dataframe[column].abs().tolist()

            # Create a UnivariateSpline object with your data
            x_values = np.arange(len(y))
            spline = UnivariateSpline(x_values, y, s=0.00)

            # Generate new, smoother y values
            ynew = spline(x_values)

            if column == "pssmBaseline_pearson_correlation" and label == "pssmB":
                ax.plot(x, ynew, label=label, color='black', linewidth=15, linestyle='--')
            else:
                ax.plot(x, ynew, label=label, color=color, linewidth=6)

        ax.set_xlabel('Protein Names', fontsize=30)
        ax.set_ylabel('Correlations', fontsize=30)
        ax.legend(title='Tool Name', title_fontsize=30, fontsize=25, labels=legend)

        # ax.set_xticks(x)
        fig.subplots_adjust(top=0.9)  # Adjust the top border
        fig.subplots_adjust(bottom=0.2)  # Adjust the bottom border
        plt.tick_params(axis='y', labelsize=30)
        ax.set_xticklabels(x, rotation=45, fontsize=20, weight='bold')
        fig.set_size_inches(22, 22)


        plt.gca().invert_xaxis()

        plt.show()

    @staticmethod
    def plot_pie_with_counts(df, column_name):
        """
        Generate a pie chart with percentage and count labels for unique entries in the specified column of a DataFrame.

        Parameters:
            df (DataFrame): The DataFrame containing the data.
            column_name (str): The name of the column to create the pie chart from.

        Returns:
            None
        """

        # Count the frequency of each unique entry in the specified column
        value_counts = df[column_name].value_counts()

        # Define a custom function to display both percentage and count

        # Plot the pie chart
        labels = [f"{label} ({count})" for label, count in zip(value_counts.index, value_counts.values)]

        # plt.pie(value_counts.values, labels=labels, startangle=140)
        plt.pie(value_counts.values, labels=labels, startangle=140, autopct='%1.1f%%')

        # value_counts.plot.pie(autopct=lambda pct: func(pct, value_counts), startangle=100)
        # value_counts.plot.pie()

        # Add a title
        plt.title(f'Distribution of Values: Counts and Percentages (Species)')

        # Show the plot
        plt.show()

    @staticmethod
    def generate_bar_plot(species_tuple, data_dict, file_name, height = 0.25):
        y = np.arange(len(species_tuple))  # the label locations
        multiplier = 0

        # fig, ax = plt.subplots(layout='constrained')
        fig, ax = plt.subplots(figsize=(10, 8), dpi=500, layout='constrained')

        for attribute, measurement in data_dict.items():
            offset = height * multiplier
            rects = ax.barh(y + offset, measurement, height, label=attribute)
            ax.bar_label(rects, padding=3, fontsize=4)
            multiplier += 1

        # Add some text for labels, title and custom y-axis tick labels, etc.
        ax.set_xlabel("Absolute Spearman's correlation")
        ax.set_ylabel('Protein Names')
        # ax.set_title('Attributes by species')
        ax.set_yticks(y + height)
        ax.set_yticklabels(species_tuple, rotation=0)
        ax.legend(loc='lower right', bbox_to_anchor=(1, 1), ncol=3, fontsize='x-small')
        # ax.legend(loc='upper right', ncols=1)
        # all_values = [value for values in data_dict.values() for value in values]
        # x_min = min(all_values)
        # x_max = max(all_values)
        #
        # ax.set_xlim(x_min, x_max+0.05*x_max)
        ax.set_xlim(0, 0.8)
        fig.savefig(file_name, format='svg')
        # plt.show()
