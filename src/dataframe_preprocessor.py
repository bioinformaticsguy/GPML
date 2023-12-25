import pathlib

import pandas
import pandas as pd
from src.utils import get_dictonary_of_scores_maveDB, get_list_to_add_in_dataframe, get_single_letter_point_mutation

DICTIONARY_PROTEINS_SPECIES = {
    'VIM-2_with_p.Met1_Phe2insGly_urn:mavedb:00000073-c': 'Human',
    'CBS_urn:mavedb:00000005-a': 'Bacteria',
    'NUDT15_urn:mavedb:00000055-0': 'E-Coli',
    'VIM-2_with_p.Met1_Phe2insGly_urn:mavedb:00000073-f': 'Virus',
    'CcdB_urn:mavedb:00000084-a': 'Bacteria',
    'VKOR_urn:mavedb:00000078-a': 'Yeast',
    'VIM-2_with_p.Met1_Phe2insGly_urn:mavedb:00000073-g': 'Virus',
    'TEM-1_beta-lactamase_urn:mavedb:00000086-d': 'Human',
    'VIM-2_with_p.Met1_Phe2insGly_urn:mavedb:00000073-d': 'Virus',
    'TEM-1_beta-lactamase_urn:mavedb:00000070-a': 'Human',
    'CCR5_urn:mavedb:00000047-c': 'Bacteria',
    'TEM-1_beta-lactamase_urn:mavedb:00000086-c': 'Human',
    'PTEN_urn:mavedb:00000013-a': 'Escherichia coli',
    'PTEN_urn:mavedb:00000054-a': 'Escherichia coli',
    'p53_urn:mavedb:00000059-a': 'Escherichia coli',
    'VKOR_urn:mavedb:00000078-b': 'Yeast',
    'TEM-1_beta-lactamase_urn:mavedb:00000086-e': 'Human',
    'TP53_(P72R)_urn:mavedb:00000068-c': 'Human',
    'SUMO1_urn:mavedb:00000001-b': 'Human',
    'NUDT15_urn:mavedb:00000055-a': 'E-Coli',
    'SARS-CoV-2_receptor_binding_domain_urn:mavedb:00000044-b': 'Human',
    'VIM-2_with_p.Met1_Phe2insGly_urn:mavedb:00000073-a': 'Human',
    'TEM-17_beta-lactamase_urn:mavedb:00000085-b': 'Human',
    'CYP2C9_urn:mavedb:00000095-a': 'Bacteria',
    'SARS-CoV-2_receptor_binding_domain_urn:mavedb:00000044-a': 'Human',
    'VIM-2_with_p.Met1_Phe2insGly_urn:mavedb:00000073-h': 'Virus',
    'CYP2C9_urn:mavedb:00000095-b': 'E-Coli',
    'A0A2Z5U3Z0_9INFA_A0A2Z5U3Z0_9INFA_Doud_2016': 'Bacteria',
    'BLAT_ECOLX_BLAT_ECOLX_Deng_2012': 'Bacteria',
    'BLAT_ECOLX_BLAT_ECOLX_Jacquier_2013': 'Bacteria',
    'CCDB_ECOLI_CCDB_ECOLI_Tripathi_2016': 'Bacteria',
    'IF1_ECOLI_IF1_ECOLI_Kelsic_2016': 'E-Coli',
    'KKA2_KLEPN_KKA2_KLEPN_Melnikov_2014': 'E-Coli',
    'Q2N0S5_9HIV1_Q2N0S5_9HIV1_Haddox_2018': 'Escherichia coli',
    'R1AB_SARS2_R1AB_SARS2_Flynn_growth_2022': 'Human',
    'RL401_YEAST_RL401_YEAST_Mavor_2016': 'Human',
}


class MaveGoldStandard:

    @staticmethod
    def get_dataframe_for_mave_gs_data(mave_gs_file_path, column_names):
        """
        Returns a pandas DataFrame from the MAVE DB gold standard data.

        Parameters:
        - mave_gs_file_path (str): The file path to the MAVE DB gold standard data.
        - column_names (list): List of column names for the resulting DataFrame.

        Returns:
        - pd.DataFrame: A DataFrame constructed from the MAVE DB gold standard data.
        """

        def _add_species_column(mave_gold_standard_df, dictionary_species, column_name="species"):
            """
                Adds a new 'species' column to the provided DataFrame based on a mapping from protein names to species.

                Parameters:
                - mave_gold_standard_df (pd.DataFrame): The DataFrame to which the 'species' column will be added.
                - dictionary_species (dict): A dictionary mapping protein names to their corresponding species.
                - column_name (str, optional): The name of the new column. Defaults to "species".

                Returns:
                - pd.DataFrame: The DataFrame with the added 'species' column.
                """

            mave_gold_standard_df[column_name] = mave_gold_standard_df['protein_name'].map(dictionary_species)
            return mave_gold_standard_df

        def _add_SNP_dict_column(mave_gold_standard_df: pd.DataFrame) -> pd.DataFrame:
            """
            Adds a new 'SNP_dict' column to the provided DataFrame based on 'SNPs' and 'scaled_effect' columns.

            Parameters:
            - df (pd.DataFrame): The DataFrame to which the 'SNP_dict' column will be added.

            Returns:
            - pd.DataFrame: The DataFrame with the added 'SNP_dict' column.
            """
            # Split 'SNPs' and 'scaled_effect' columns, then create a dictionary for each row
            mave_gold_standard_df['SNP_dict'] = mave_gold_standard_df.apply(
                lambda row: dict(zip(map(get_single_letter_point_mutation, row['SNPs'].split(';')),
                                     map(float, row['scaled_effect'].split(';')))),
                axis=1
            )

            return mave_gold_standard_df


        dictionary_of_data = get_dictonary_of_scores_maveDB(mave_gs_file_path)

        rows = []
        for protein_name, value in dictionary_of_data.items():
            row = [protein_name] + get_list_to_add_in_dataframe(value)
            rows.append(row)
        mave_goldstandard_dataframe = pd.DataFrame(rows, columns=column_names)

        mave_goldstandard_dataframe = _add_species_column(mave_goldstandard_dataframe,
                            dictionary_species=DICTIONARY_PROTEINS_SPECIES,
                            column_name="species")

        mave_goldstandard_dataframe = _add_SNP_dict_column(mave_goldstandard_dataframe)

        return mave_goldstandard_dataframe


    @staticmethod
    def mark_rows_present_in_subset(
            superset_df, subset_df, new_column_name='in_subset_dataframe',
            id_columns=('protein_name', 'protein_name')
    ):
        """
        Add a binary column to the superset dataframe indicating whether each row is present in the subset dataframe.

        Parameters:
        - superset_df (pd.DataFrame): The dataframe that serves as the superset.
        - subset_df (pd.DataFrame): The dataframe that serves as the subset.
        - new_column_name (str): The name of the new binary column to be added. Default is 'in_subset_dataframe'.
        - id_columns (tuple): Tuple containing the column names used as identifiers in superset and subset dataframes,
                            respectively. Default is ('id', 'id').

        Returns:
        pd.DataFrame: The superset dataframe with the new binary column added.

        Example:
        superset_df = pd.DataFrame({'id': [1, 2, 3, 4], 'other_column': ['A', 'B', 'C', 'D']})
        subset_df = pd.DataFrame({'id': [2, 4], 'other_column': ['B', 'D']})
        result_df = DataFrameUtility.mark_rows_present_in_subset(superset_df, subset_df)

        Output:
        >> result_df
           id other_column  in_subset_dataframe
        0   1            A                   0
        1   2            B                   1
        2   3            C                   0
        3   4            D                   1
        """
        superset_id_col, subset_id_col = id_columns

        if (
                superset_id_col not in superset_df.columns or
                subset_id_col not in subset_df.columns
        ):
            raise ValueError(
                f"Column '{superset_id_col}' or '{subset_id_col}' not found in respective dataframes."
            )

        superset_df[new_column_name] = superset_df[superset_id_col].isin(subset_df[subset_id_col]).astype(int)
        return superset_df


class MutepredTrainingProcessor:
    @staticmethod
    def get_mutepred_df(file_path: pathlib.Path) -> pandas.DataFrame:
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

class dbNSFPProcessor:
    @staticmethod
    def get_dbNSFP_df(file_path: pathlib.Path) -> pandas.DataFrame:
        """
        Input: str file path
        Output: returns a pandas dataframe of whole data for a single protein.
        """
        df = pd.read_csv(file_path, sep='\t')

        return df

    @staticmethod
    def check_if_all_dbNSFP_file_have_same_names_as_mave_goldstandard(mave_goldstandard_df, outputdirectory):
        # List all file names in the directory
        csv_file_names = [file.name for file in outputdirectory.iterdir() if
                          file.is_file() and file.name.endswith('.csv')]
        # file_names = [file.name for file in outputdirectory.iterdir() if file.is_file()]
        modified_list = [filename.replace('_output.csv', '') for filename in csv_file_names]

        first_column_values = mave_goldstandard_df.iloc[:, 0].tolist()

        unique_to_list1 = list(set(first_column_values) - set(modified_list))

        # Values unique to list2
        unique_to_list2 = list(set(modified_list) - set(first_column_values))

        return unique_to_list1, unique_to_list2, modified_list, first_column_values


