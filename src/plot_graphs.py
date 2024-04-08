from itertools import cycle
import matplotlib.pyplot as plt

class PlotGeneroator:
    @staticmethod
    def plot_correlations(ONLY_HUMAN_DATABASE, PROTEIN_SHORT_MAPPING, PEARSON_CORELATION_SUFFIX):
        columns_list = ONLY_HUMAN_DATABASE.columns.tolist()
        filtered_columns_list = [column for column in columns_list if column.endswith(PEARSON_CORELATION_SUFFIX)]

        legend = [value.rstrip(PEARSON_CORELATION_SUFFIX) for value in filtered_columns_list]

        x = [PROTEIN_SHORT_MAPPING[protein_name] for protein_name in ONLY_HUMAN_DATABASE['protein_name'].tolist()]

        # Generating dynamic colors based on the number of legend entries
        num_colors = len(legend)
        color_cycle = cycle(plt.cm.get_cmap('tab10').colors)
        colors = [next(color_cycle) for _ in range(num_colors)]

        fig, ax = plt.subplots()

        for column, label, color in zip(filtered_columns_list, legend, colors):
            y = ONLY_HUMAN_DATABASE[column].abs().tolist()

            # # Add lines to join the points
            ax.plot(x, y, label=label, color=color, linewidth=4)
            # ax.scatter(x, y, label=label, color=color, linewidth=4, s=100)


        ax.set_xlabel('Protein Names', fontsize=20)
        ax.set_ylabel('Correlations', fontsize=30)
        ax.legend(title='Tool Name', fontsize=25, labels=legend)

        ax.set_xticks(x)
        plt.tick_params(axis='y', labelsize=30)
        ax.set_xticklabels(x, rotation=45, fontsize=20, weight='bold')
        fig.set_size_inches(21, 21)
        fig.subplots_adjust(bottom=0.1)

        plt.gca().invert_xaxis()
        # plt.gca().invert_yaxis()

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