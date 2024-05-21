import pandas as pd
from src.constants import AMINO_ACIDS_SINGLE_LETTER, COLUMN_NAME_OF_MAVE_GOLD_STANDARD_ID, \
    COLUMN_NAME_OF_MAVE_GOLD_STANDARD_SNP_DICTIONARY
from src.utils import get_value_from_dataframe, get_protein_list, remove_digits_from_key


class pssmBaseline:
    @staticmethod
    def get_empty_pssm_dataframe(amino_acid_list=AMINO_ACIDS_SINGLE_LETTER):
        """
        This function creates an empty DataFrame with amino_acids as both row and column names.
        Each cell is filled with the tuple (0, 0).
        Input: amino_acid_list
        Output: df with amino_acids as both row and column names with (0, 0) in each cell.
        """
        return pd.DataFrame(data=[[(0, 0) for _ in amino_acid_list] for _ in amino_acid_list],
                            index=amino_acid_list, columns=amino_acid_list)

    @staticmethod
    def get_lopo_dict(protein_name_list):
        """
        Generate a dictionary where each key is a protein name and the value is a list of proteins with that protein left out.

        Parameters:
        - protein_name_list (list): The list of protein names.

        Returns:
        - dict: A dictionary where each key is a protein name and the value is a list of proteins with that protein left out.
        """
        lopo_dict = {}
        for i in range(len(protein_name_list)):
            lopo_list = protein_name_list[:i] + protein_name_list[i + 1:]
            lopo_dict[protein_name_list[i]] = lopo_list
        return lopo_dict

    @staticmethod
    def get_protein_pssm(dataframe, lopo_list):
        """
        This function calculates the mean of the PSSM values for each amino acid for a protein.
        Input: protein_name, lopo_list, df
        Output: amino_acid_pssm_dict.
        """

        pssm_df = pssmBaseline.get_empty_pssm_dataframe()
        for protein in lopo_list:
            protein_dict = get_value_from_dataframe(df=dataframe,
                                                    id_column=COLUMN_NAME_OF_MAVE_GOLD_STANDARD_ID,
                                                    id_value=protein,
                                                    value_column=COLUMN_NAME_OF_MAVE_GOLD_STANDARD_SNP_DICTIONARY)
            for key, value in protein_dict.items():
                if key[0] != key[-1]:
                    cur_tuple = pssm_df.at[key[0], key[-1]]
                    pssm_df.at[key[0], key[-1]] = (cur_tuple[0] + value, cur_tuple[1] + 1)
        return pssm_df


    @staticmethod
    def get_mean_from_pssm(key, pssm_df):
        """
        This function returns the mean value from the PSSM DataFrame for a given key.
        input: key, pssm_df
        output: mean_value
        """
        if pssm_df.at[key[0], key[-1]][1] == 0:
            return 0
        return pssm_df.at[key[0], key[-1]][0] / pssm_df.at[key[0], key[-1]][1]

    @staticmethod
    def update_value_in_snp_scores_dict(dataframe,
                                        id_column,
                                        protein_name,
                                        dict_col_name,
                                        key,
                                        value):

        dataframe.loc[dataframe[id_column] == protein_name, dict_col_name].values[0][key] = value
        return dataframe

    @staticmethod
    def update_all_values_in_snp_scores_dict(dataframe,
                                             id_column,
                                             protein_name,
                                             dict_col_name,
                                             default_value):
        temp_keys = dataframe.loc[dataframe[id_column] == protein_name, dict_col_name].values[0].keys()
        if default_value == None:
            for key in temp_keys:
                pssmBaseline.update_value_in_snp_scores_dict(dataframe,
                                                id_column,
                                                protein_name,
                                                dict_col_name,
                                                key,
                                                default_value)
        else:
            lopo_dict = pssmBaseline.get_lopo_dict(get_protein_list(dataframe=dataframe))
            protein_pssm = pssmBaseline.get_protein_pssm(dataframe=dataframe, lopo_list=lopo_dict[protein_name])
            for key in temp_keys:
                mean_value = pssmBaseline.get_mean_from_pssm(remove_digits_from_key(key), protein_pssm)
                pssmBaseline.update_value_in_snp_scores_dict(dataframe,
                                                id_column,
                                                protein_name,
                                                dict_col_name,
                                                key,
                                                mean_value)
        return dataframe




