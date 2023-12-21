import pandas as pd

from src.utils import MAVE_DB_GOLD_STANDARD_SEQUENCE_ONLY_FILE_PATH, read_fasta_iteration


class MutepredDataPreprocessor:
    @staticmethod
    def get_mutepred_df(file_path):
        """
        Input: str file path
        Output: returns a pandas dataframe of whole data
        """
        df = pd.read_csv(file_path, skiprows=8, sep='\t')

        return df

    @staticmethod
    def filter_data_with_common_sequences(dataframe_to_filter,
                                          dataframe_with_common_values,
                                          df_sequence_column_name):
        """
        This function takes tow dataframes and a string which is the name
        of the columns, then it filters the dataset on the basis of the common
         column in both dataframes.
        Input:
            dataframe_to_filter -> pandas.dataframe
            dataframe_with_common_values -> pandas.dataframe
            sequence_list -> python list with strings of sequences
            df_sequence_column_name -> name of column containing sequences

        Output:
            filtered_df ->  pandas dataframe (filtered)
        """
        return dataframe_to_filter[dataframe_to_filter[df_sequence_column_name]. \
                isin(dataframe_with_common_values[df_sequence_column_name])]

    @staticmethod
    def get_disjunction(superset_df, subset_df, on_column='ID'):
        """
            Get the rows from the superset DataFrame that are not present in the subset DataFrame.

            Parameters:
            - superset_df (pd.DataFrame): The superset DataFrame.
            - subset_df (pd.DataFrame): The subset DataFrame.
            - on_column (str, optional): The column used for merging the DataFrames. Defaults to 'ID'.

            Returns:
            - pd.DataFrame: A DataFrame containing rows from the superset that are not present in the subset.
            """

        result_df = superset_df.merge(subset_df, on=on_column, how='left', indicator=True)
        result_df = result_df[result_df['_merge'] == 'left_only'].drop(columns=['_merge']).rename(
            columns={f"{on_column}_x": on_column})

        result_df = result_df[[col for col in result_df.columns if not col.endswith('_y')]]
        result_df.columns = [col.replace('_x', '') for col in result_df.columns]

        return result_df



class MutepredDataAnalyzer:
    @staticmethod
    def count_overlaping_seqences(gold_standard_sequences, mutepred_sequences):
        count = 0
        for gold_standard_sequence in gold_standard_sequences:
            if gold_standard_sequence in mutepred_sequences:
                count += 1
        return count
