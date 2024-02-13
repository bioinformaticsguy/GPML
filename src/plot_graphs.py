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
            ax.scatter(x, y, label=label, color=color, linewidth=4)

        ax.set_xlabel('Protein Names', fontsize=14)
        ax.set_ylabel('Correlations', fontsize=14)
        ax.legend(title='Tool Name', fontsize=12, labels=legend)
        ax.set_xticks(x)
        ax.set_xticklabels(x, rotation=45, fontsize=14)
        fig.set_size_inches(20, 20)
        fig.subplots_adjust(bottom=0.1)

        plt.show()

