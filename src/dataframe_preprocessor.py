import copy
import pathlib
import re
import pandas
import pandas as pd

from src.constants import COLUMN_NAME_OF_MAVE_GOLD_STANDARD_ID, COLUMN_NAME_OF_MAVE_GOLD_STANDARD_SNP, \
    COLUMN_NAME_OF_MAVE_GOLD_STANDARD_SCALED_EFFECT, COLUMN_NAME_OF_MAVE_GOLD_STANDARD_SPECIES, \
    COLUMN_NAME_OF_MAVE_GOLD_STANDARD_SNP_DICTIONARY, PROTEIN_SPECIES_MAPPING
from src.utils import get_dictonary_of_scores_maveDB, get_list_to_add_in_dataframe, get_single_letter_point_mutation, \
    get_protein_names_from_db_nsfp_output_directory

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

        def _add_species_column(mave_gold_standard_df,
                                dictionary_species,
                                column_name=COLUMN_NAME_OF_MAVE_GOLD_STANDARD_SPECIES):
            """
                Adds a new 'species' column to the provided DataFrame based on a mapping from protein names to species.

                Parameters:
                - mave_gold_standard_df (pd.DataFrame): The DataFrame to which the 'species' column will be added.
                - dictionary_species (dict): A dictionary mapping protein names to their corresponding species.
                - column_name (str, optional): The name of the new column. Defaults to "species".

                Returns:
                - pd.DataFrame: The DataFrame with the added 'species' column.
                """

            mave_gold_standard_df[column_name] = mave_gold_standard_df[COLUMN_NAME_OF_MAVE_GOLD_STANDARD_ID]. \
                                                    map(dictionary_species)

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
            mave_gold_standard_df[COLUMN_NAME_OF_MAVE_GOLD_STANDARD_SNP_DICTIONARY] = mave_gold_standard_df.apply(
                lambda row: dict(zip(map(get_single_letter_point_mutation, row[COLUMN_NAME_OF_MAVE_GOLD_STANDARD_SNP].split(';')),
                                     map(float, row[COLUMN_NAME_OF_MAVE_GOLD_STANDARD_SCALED_EFFECT].split(';')))),
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
                                                          dictionary_species=PROTEIN_SPECIES_MAPPING,
                                                          column_name=COLUMN_NAME_OF_MAVE_GOLD_STANDARD_SPECIES)

        mave_goldstandard_dataframe = _add_SNP_dict_column(mave_goldstandard_dataframe)

        return mave_goldstandard_dataframe

    @staticmethod
    def mark_rows_present_in_subset(
            superset_df, subset_df, new_column_name,
            id_columns=(COLUMN_NAME_OF_MAVE_GOLD_STANDARD_ID, COLUMN_NAME_OF_MAVE_GOLD_STANDARD_ID)
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
    def get_disjunction(superset_df, subset_df, on_column=COLUMN_NAME_OF_MAVE_GOLD_STANDARD_ID):
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
    def get_dbNSFP_df(protein_csv_file_path: pathlib.Path) -> pandas.DataFrame:
        """
        Input: str file path
        Output: returns a pandas dataframe of whole data for a single protein.
        """
        df = pd.read_csv(protein_csv_file_path, sep='\t')

        return df

    @staticmethod
    def file_names_are_similar_to_mave_goldstandard_protein_names(mave_goldstandard_df, db_nsfp_output_path):
        """
        Checks if file names in the specified directory are similar to protein names in the provided DataFrame.

        Parameters:
        - mave_goldstandard_df (pd.DataFrame): DataFrame containing protein names.
        - dbNSFP_output_path (Path): Path object representing the directory with CSV files to be checked.

        Returns:
        - bool: True if all file names are similar to protein names, False otherwise.
        """

        protein_names_from_csv_files = get_protein_names_from_db_nsfp_output_directory(directory_path=db_nsfp_output_path)

        protein_names_in_mave_gold_standard_df = mave_goldstandard_df.iloc[:, 0].tolist()

        # Values unique to list2
        unique_names_in_output_directory = list(
            set(protein_names_from_csv_files) - set(protein_names_in_mave_gold_standard_df)
        )

        return True if len(unique_names_in_output_directory) == 0 else False

    @staticmethod
    def get_protein_dictionary_snps_scores(dbNSFP_protein_dataframe, snp_column_name, tool_score_column_name):

        def _extract_snp_and_score(snp_str, score_str):
            """
            Extract SNPs and their corresponding scores from input strings.

            Parameters:
            - snp_str (str): String containing SNP information.
            - score_str (str): String containing score information.

            Returns:
            - list: A list of tuples where each tuple contains an extracted SNP and its corresponding score.
              Returns an empty list if no valid SNPs are found.
            """
            # Extracting SNPs and their score
            snps = re.findall(r'p\.([A-Z0-9]+)', snp_str)
            scores = str(score_str).split(';')[:len(snps)]

            # Return the extracted SNP and its score or None if not valid
            return list(set(list(zip(snps, scores))))

        values_list = []

        for _, row in dbNSFP_protein_dataframe.iterrows():
            for tuple_snp_score in _extract_snp_and_score(snp_str=row[snp_column_name],
                                                          score_str=row[tool_score_column_name]):
                values_list.append(tuple_snp_score)

        filtered_tuples = [(snp, float(score)) for snp, score in values_list if score != '.']

        return dict(filtered_tuples)

    @staticmethod
    def add_tool_score_column(mave_gs_dataframe, db_nsfp_output_dir_path, tool_name, snp_column_name):
        """
        Add a tool score column to a copy of the MAVE GoldStandard DataFrame.

        Parameters:
        - mave_gs_dataframe (pd.DataFrame): The MAVE GoldStandard DataFrame.
        - db_nsfp_output_dir_path (pathlib.Path): The path to the directory containing dbNSFP output files.
        - tool_name (str): The name of the tool for which the score is being added.
        - snp_column_name (str): The name of the column containing SNPs in dbNSFP output files.

        Returns:
        - pd.DataFrame: A deep copy of the MAVE GoldStandard DataFrame with the added tool score column.
        """

        # Create a deep copy of the MAVE GoldStandard DataFrame
        mave_gs_df_deep_copy = copy.deepcopy(mave_gs_dataframe)

        # Create the column name for the tool score
        column_name = tool_name + TOOL_SCORE_COLUMN_SUFFIX

        # Add a new column to the deep copy DataFrame
        mave_gs_df_deep_copy[column_name] = None

        # Get a list of file paths for dbNSFP output files in the specified directory
        csv_file_names = [file for file in db_nsfp_output_dir_path.iterdir() if
                          file.is_file() and file.name.endswith('.csv')]

        # Iterate over each dbNSFP output file
        for csv_file_path in csv_file_names:
            # Extract protein name from the file path
            protein_name = csv_file_path.stem.replace("_output", "").replace("_urn_", "_urn:")

            # Get the DataFrame for the current protein from the dbNSFP output file
            protein_dataframe = dbNSFPProcessor.get_dbNSFP_df(protein_csv_file_path=csv_file_path)

            # Get the dictionary of SNPs and scores for the current tool
            protein_snp_score_dictionary = dbNSFPProcessor.get_protein_dictionary_snps_scores(
                dbNSFP_protein_dataframe=protein_dataframe,
                snp_column_name=snp_column_name,
                tool_score_column_name=column_name
            )

            # Update the tool score column for the corresponding protein in the deep copy DataFrame
            mave_gs_df_deep_copy.loc[
                mave_gs_df_deep_copy[COLUMN_NAME_OF_MAVE_GOLD_STANDARD_ID] == protein_name, column_name] = [
                protein_snp_score_dictionary]

        return mave_gs_df_deep_copy