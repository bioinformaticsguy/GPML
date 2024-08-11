from itertools import cycle
from operator import itemgetter
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import UnivariateSpline

from src.constants import PROTEIN_SHORT_DICTMAP, SPEAR_COR_SUFFIX, PLOT_FORMAT, EXCLUDE_TRAINING_SAV_SUFFIX, \
    TOOL_SCORE_COLUMN_SUFFIX, COLUMN_NAME_OF_MAVE_GOLD_STANDARD_ID, TRAINING_SAVS_COLUMN_SIFFIX, PIE_PLOT_FILE_NAME, \
    PLOTS_DIRECTORY_PATH
from src.utils import extract_value


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
    def plot_pie_with_counts(df, column_name,
                             plot_path=PLOTS_DIRECTORY_PATH,
                             file_name=PIE_PLOT_FILE_NAME,
                             plot_format=PLOT_FORMAT,
                             title=False):
        """
        Generate a pie chart with percentage and count labels for unique entries in the specified column of a DataFrame.

        Parameters:
            df (DataFrame): The DataFrame containing the data.
            column_name (str): The name of the column to create the pie chart from.

        Returns:
            None
        """

        plot_path = Path(plot_path) / f"{file_name}.{plot_format}"

        # Count the frequency of each unique entry in the specified column
        value_counts = df[column_name].value_counts()

        # Define a custom function to display both percentage and count

        # Plot the pie chart
        labels = [f"{label} ({count})" for label, count in zip(value_counts.index, value_counts.values)]

        # Create a figure object
        fig = plt.figure()

        # plt.pie(value_counts.values, labels=labels, startangle=140)
        plt.pie(value_counts.values, labels=labels, startangle=140, autopct='%1.1f%%')

        # value_counts.plot.pie(autopct=lambda pct: func(pct, value_counts), startangle=100)
        # value_counts.plot.pie()

        # Add a title
        if title:
            plt.title(f'Distribution of Values: Counts and Percentages (Species)')

        fig.savefig(plot_path, format=PLOT_FORMAT)


        # Show the plot
        plt.show()

    @staticmethod
    def generate_bar_plot(protein_names_list, data_dict, file_name, main_dataframe,
                          padding=0.01,  # adjust this value as needed
                          height = 0.25,
                          fig_height=12,
                          legend_font_size="xx-large",
                          barlabel_font_size=8,
                          barlabel_flag=False,
                          removed_snp_flag_value=True,):

        short_protein_names = [PROTEIN_SHORT_DICTMAP[name] for name in protein_names_list]

        y = np.arange(len(short_protein_names))  # the label locations
        multiplier = 0

        # fig, ax = plt.subplots(layout='constrained')
        fig, ax = plt.subplots(figsize=(10, fig_height), dpi=500, layout='constrained')

        for attribute, measurement in data_dict.items():
            offset = height * multiplier
            rects = ax.barh(y + offset, measurement, height, label=attribute)

            for i, rect in enumerate(rects):
                if attribute.replace("_excluded_training_savs", "") != "pssmBaseline" and "_excluded_training_savs" in attribute:
                    exc_savs = extract_value(df=main_dataframe,
                                             id_column=COLUMN_NAME_OF_MAVE_GOLD_STANDARD_ID,
                                             row_value=protein_names_list[i],
                                             target_col_name=attribute.replace("_excluded_training_savs", "")+TRAINING_SAVS_COLUMN_SIFFIX)

                    all_savs = extract_value(df=main_dataframe,
                                             id_column=COLUMN_NAME_OF_MAVE_GOLD_STANDARD_ID,
                                             row_value=protein_names_list[i],
                                             target_col_name=attribute.replace("_excluded_training_savs", "")+TOOL_SCORE_COLUMN_SUFFIX)



                    # mave_sav_num = "mave_savs_scores_dictinary"

                    devided = round(exc_savs / all_savs, 3)

                    perc = "{:.2f}%".format(devided * 100)
                    #
                    # print(protein_names_list[i],
                    #       attribute.replace("_excluded_training_savs", "")+EXCLUDE_TRAINING_SAV_SUFFIX,
                    #       attribute.replace("_excluded_training_savs", "")+TOOL_SCORE_COLUMN_SUFFIX,
                    #       exc_savs, all_savs, devided, perc, mave_sav_num)
                    # print(f"Width: {rect.get_width()}, Y: {rect.get_y()}, Height: {rect.get_height()}")

                    # To put in center: rect.get_width() / 2 or -0.00001
                    ax.text(rect.get_width() + padding, rect.get_y() + rect.get_height() / 2,
                            '(' + str(exc_savs) + ' / ' + str(all_savs) + " = " + str(perc) + ')',
                            ha='left', va='center', fontsize=barlabel_font_size, color='black', weight='normal')

                if attribute.replace("_excluded_training_savs", "") != "pssmBaseline" and "_excluded_training_savs" not in attribute and removed_snp_flag_value:
                    mave_sav_num = extract_value(df=main_dataframe,
                                             id_column=COLUMN_NAME_OF_MAVE_GOLD_STANDARD_ID,
                                             row_value=protein_names_list[i],
                                             target_col_name="mave_savs_scores_dictinary")

                    all_savs = extract_value(df=main_dataframe,
                                             id_column=COLUMN_NAME_OF_MAVE_GOLD_STANDARD_ID,
                                             row_value=protein_names_list[i],
                                             target_col_name=attribute.replace("_excluded_training_savs", "") + TOOL_SCORE_COLUMN_SUFFIX )

                    devided = round(all_savs / mave_sav_num  , 3)

                    perc = "{:.2f}%".format(devided * 100)


                    ax.text(rect.get_width() + padding, rect.get_y() + rect.get_height() / 2,
                             '(' + str(all_savs) + ' / ' + str(mave_sav_num) + " = " + str(perc) + ')' ,
                            ha='left', va='center', fontsize=barlabel_font_size, color='black', weight='normal')

                    print(protein_names_list[i], all_savs, mave_sav_num, attribute)

            if barlabel_flag:
                ax.bar_label(rects, padding=3, fontsize=barlabel_font_size, weight='normal')

            multiplier += 1

        # Add some text for labels, title and custom y-axis tick labels, etc.
        ax.set_xlabel("Absolute Spearman's correlation")
        ax.set_ylabel('Protein Names')
        # ax.set_title('Attributes by species')
        ax.set_yticks(y + height)
        ax.set_yticklabels(short_protein_names, rotation=0)
        ax.legend(loc='lower right', bbox_to_anchor=(1, 1), ncol=3, fontsize=legend_font_size)
        # ax.legend(loc='upper right', ncols=1)
        # all_values = [value for values in data_dict.values() for value in values]
        # x_min = min(all_values)
        # x_max = max(all_values)

        # ax.set_xlim(x_min, x_max+0.05*x_max)
        ax.set_xlim(0, 0.75)
        fig.savefig(file_name, format=PLOT_FORMAT)
        plt.show()
