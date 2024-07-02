import copy
import pathlib
import re
import pandas
import pandas as pd
import numpy as np
import ast


from src.constants import COLUMN_NAME_OF_MAVE_GOLD_STANDARD_ID, COLUMN_NAME_OF_MAVE_GOLD_STANDARD_SAVS, \
    COLUMN_NAME_OF_MAVE_GOLD_STANDARD_SCALED_EFFECT, COLUMN_NAME_OF_MAVE_GOLD_STANDARD_SPECIES, \
    COLUMN_NAME_OF_MAVE_GOLD_STANDARD_SAV_DICTIONARY, PROTEIN_SPECIES_DICTMAP, TOOL_SCORE_COLUMN_SUFFIX, \
    UNIPROT_ID_DICTMAP, COLUMN_NAME_OF_MAVE_GOLD_STANDARD_UNIPROT_ID, OUTPUT_DIR_DB_NSFP, TOOLS_LIST, \
    DBNSFP_SAV_COLUMN_NAME, DEOGEN_COLUMN_NAME_TO_FILTER, DEOGEN_VALUES_TO_FILTER, DEOGEN_UNIPROT_ID_COLUMN, \
    DEOGEN_AMINO_ACID_CHANGE_COLUMN, CLINVAR_DATA_FILE_PATH

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

        def _add_meta_data_column(mave_gold_standard_df,
                                  dictionary_meta,
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
                                                    map(dictionary_meta)

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
            mave_gold_standard_df[COLUMN_NAME_OF_MAVE_GOLD_STANDARD_SAV_DICTIONARY] = mave_gold_standard_df.apply(
                lambda row: dict(zip(map(get_single_letter_point_mutation, row[COLUMN_NAME_OF_MAVE_GOLD_STANDARD_SAVS].split(';')),
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
        mave_goldstandard_dataframe = _add_meta_data_column(mave_goldstandard_dataframe,
                                                            dictionary_meta=PROTEIN_SPECIES_DICTMAP,
                                                            column_name=COLUMN_NAME_OF_MAVE_GOLD_STANDARD_SPECIES)

        mave_goldstandard_dataframe = _add_meta_data_column(mave_goldstandard_dataframe,
                                                            dictionary_meta=UNIPROT_ID_DICTMAP,
                                                            column_name=COLUMN_NAME_OF_MAVE_GOLD_STANDARD_UNIPROT_ID)




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


class ClinPredTrainingProcessor:
    @staticmethod
    def get_clinvar_df(file_path=CLINVAR_DATA_FILE_PATH):
        """
        this function reads the clinvar data file and returns a pandas dataframe
        Input: str file path
        Output: returns a pandas dataframe of whole data
        """
        header_list = [COLUMN_NAME_OF_MAVE_GOLD_STANDARD_UNIPROT_ID, "snp", "val"]
        df = pd.read_csv(file_path, sep='\t')
        df.columns = header_list

        return df

    @staticmethod
    def add_savs_of_one_protein(mave_df,
                                clinphred_df,
                                protein_name,
                                clinphred_sav_column_name,
                                protein_name_id_col_name=COLUMN_NAME_OF_MAVE_GOLD_STANDARD_ID,
                                uniprot_id_col_name=COLUMN_NAME_OF_MAVE_GOLD_STANDARD_UNIPROT_ID):

        uniprot_id = mave_df.loc[mave_df[protein_name_id_col_name] == protein_name, uniprot_id_col_name].values[0]
        filtered_df = clinphred_df[clinphred_df[COLUMN_NAME_OF_MAVE_GOLD_STANDARD_UNIPROT_ID] == uniprot_id]
        sav_list = filtered_df["snp"].tolist()

        if clinphred_sav_column_name not in mave_df.columns:
            mave_df[clinphred_sav_column_name] = np.nan

        if len(sav_list) > 1:
            mave_df.loc[mave_df[protein_name_id_col_name] == protein_name, clinphred_sav_column_name] = str(sav_list)

        return mave_df

    @staticmethod
    def add_training_col_for_all_proteins(mave_gs_df,
                                          clinphred_df,
                                          clinphred_sav_column_name,
                                          protein_name_id_col_name=COLUMN_NAME_OF_MAVE_GOLD_STANDARD_ID,):

        protein_names = mave_gs_df[protein_name_id_col_name].tolist()
        for protein_name in protein_names:
            mave_gs_df = ClinPredTrainingProcessor.add_savs_of_one_protein(mave_df=mave_gs_df,
                                                                           clinphred_df=clinphred_df,
                                                                           protein_name=protein_name,
                                                                           clinphred_sav_column_name=clinphred_sav_column_name)
        for i in range(len(mave_gs_df)):
            training_snps = mave_gs_df.at[i, clinphred_sav_column_name]
            if isinstance(training_snps, float):
                mave_gs_df.at[i, clinphred_sav_column_name] = []
            if isinstance(training_snps, str):
                mave_gs_df.at[i, clinphred_sav_column_name] = ast.literal_eval(mave_gs_df.at[i, clinphred_sav_column_name])

        return  mave_gs_df







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


class Deogen2TrainingProcessor:

    @staticmethod
    def get_deogen2_training_df(file_path, column_names):
        """
        Input: str file path
        Output: returns a pandas dataframe of whole data
        """
        list_of_lines = []
        with open(file_path) as file:
            lines = file.readlines()
            for line in lines[30:-6]:
                # list_of_lines.append(line.split())
                curr_list = line[:74].split() + [line[74:].strip()]
                list_of_lines.append(curr_list)

        df = pd.DataFrame(list_of_lines, columns=column_names)
        return df

    @staticmethod
    def filter_unwanted_rows(deogen_training_df,
                             column_name=DEOGEN_COLUMN_NAME_TO_FILTER,
                             row_values_to_filter=DEOGEN_VALUES_TO_FILTER):
        """
         Filter unwanted rows from a dataframe.
         Input: df: pandas dataframe.
         Output: pandas dataframe.
        """

        for value in row_values_to_filter:
            deogen_training_df = deogen_training_df.loc[deogen_training_df[column_name] != value]
        return deogen_training_df

    @staticmethod
    def get_overlaping_snp_list(deogen_training_df,
                                uniprot_id,
                                uniPprot_id_column_name=DEOGEN_UNIPROT_ID_COLUMN,
                                amino_acid_change_column_name=DEOGEN_AMINO_ACID_CHANGE_COLUMN, ):

        """
        This function returns a list of SNPs for a given uniprot_id from the deogen2 training data.
        Input:
            deogen_training_df: pd.DataFrame
            uniprot_id: str
            uniPprot_id_column_name: str
            amino_acid_change_column_name: str
        Output:
            snp_list: list
        """

        uniprot_id_df = deogen_training_df.loc[deogen_training_df[uniPprot_id_column_name] == uniprot_id]
        snp_list = [get_single_letter_point_mutation(snp[2:]) for snp in
                    uniprot_id_df[amino_acid_change_column_name].tolist()]

        return snp_list

    @staticmethod
    def add_training_col_for_one_protein(mave_gs_df, deogen_training_df, protein_name, deogen_traininig_snp_col,
                                         id_col=COLUMN_NAME_OF_MAVE_GOLD_STANDARD_ID,
                                         uniprot_id_col=COLUMN_NAME_OF_MAVE_GOLD_STANDARD_UNIPROT_ID,):

        """
        This function adds the training column for one protein to the mave_gs_df
        Input: mave_gs_df: DataFrame, deogen_training_df: DataFrame, protein_name: str
        Output: mave_gs_df: DataFrame
        """

        if deogen_traininig_snp_col not in mave_gs_df.columns:
            mave_gs_df[deogen_traininig_snp_col] = np.nan


        curr_uni_prot_id = str(mave_gs_df.loc[mave_gs_df[id_col] == protein_name, uniprot_id_col].values[0])
        list_of_snps = Deogen2TrainingProcessor.get_overlaping_snp_list(deogen_training_df, curr_uni_prot_id, )

        if len(list_of_snps) > 1:
            mave_gs_df.loc[mave_gs_df[id_col] == protein_name, deogen_traininig_snp_col] = str(list_of_snps)


        return mave_gs_df

    @staticmethod
    def add_training_col_for_all_proteins(mave_gs_df, deogen_training_df, deogen_traininig_snp_col,
                                          id_col=COLUMN_NAME_OF_MAVE_GOLD_STANDARD_ID,):

        protein_names = mave_gs_df[id_col].tolist()
        for protein_name in protein_names:
            mave_gs_df = Deogen2TrainingProcessor.add_training_col_for_one_protein(mave_gs_df,
                                                                                   deogen_training_df,
                                                                                   protein_name,
                                                                                   deogen_traininig_snp_col, )
        for i in range(len(mave_gs_df)):
            training_snps = mave_gs_df.at[i, deogen_traininig_snp_col]
            if isinstance(training_snps, float):
                mave_gs_df.at[i, deogen_traininig_snp_col] = []
            if isinstance(training_snps, str):
                mave_gs_df.at[i, deogen_traininig_snp_col] = ast.literal_eval(mave_gs_df.at[i, deogen_traininig_snp_col])

        return  mave_gs_df


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
    def add_tool_score_column(mave_gs_df, db_nsfp_output_dir_path, tool_name, snp_column_name):
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
        # mave_gs_df_deep_copy = copy.deepcopy(mave_gs_dataframe)
        # # mave_gs_df_deep_copy = mave_gs_dataframe

        # Create the column name for the tool score
        column_name = tool_name + TOOL_SCORE_COLUMN_SUFFIX

        # Add a new column to the deep copy DataFrame
        mave_gs_df[column_name] = None

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
            mave_gs_df.loc[
                mave_gs_df[COLUMN_NAME_OF_MAVE_GOLD_STANDARD_ID] == protein_name, column_name] = [
                protein_snp_score_dictionary]

        return mave_gs_df

    @staticmethod
    def add_data_from_list_of_tools(mave_gs_dataframe,
                                    db_nsfp_output_dir_path=OUTPUT_DIR_DB_NSFP,
                                    tool_list=TOOLS_LIST,
                                    snp_column_name=DBNSFP_SAV_COLUMN_NAME):
        """
        Add tool scores to the MAVE Gold Standard DataFrame
        Input: MAVE Gold Standard DataFrame, dbNSFP output directory path, list of tools, SNP column name
        Output: MAVE Gold Standard DataFrame with all tool scores added.
        """

        for tool in tool_list:
            dbNSFPProcessor.add_tool_score_column(mave_gs_dataframe,
                                                  db_nsfp_output_dir_path,
                                                  tool,
                                                  snp_column_name)

        return mave_gs_dataframe


if __name__ == '__main__':
    pass