import json

from src.constants import AMINO_ACIDS_SINGLE_LETTER, COLUMN_NAME_OF_MAVE_GOLD_STANDARD_ID, \
    COLUMN_NAME_OF_MAVE_GOLD_STANDARD_SNP_DICTIONARY
from src.utils import update_value_based_on_protein_name, get_value_from_dataframe, generate_amino_pssm_dict, \
    remove_digits_from_key, calculate_mean, get_protein_name_list


class LopoBaseline:

    @staticmethod
    def calculate_protein_mean(protein_snp_dict):
        """
        Calculates the mean of SNP values for a protein.

        Parameters:
        - protein_snp_dict (dict): Dictionary containing SNP values for a protein.

        Returns:
        float: Mean SNP value for the protein.
        """
        return sum(protein_snp_dict.values()) / len(protein_snp_dict)

    @staticmethod
    def calculate_and_update_mean(df, protein_name_column, dictionary_column, mean_column):
        """
        Calculate the mean of SNP values for each protein and update the mean column.

        Parameters:
        - df (pd.DataFrame): The DataFrame to update.
        - protein_name_column (str): The name of the column containing protein names.
        - dictionary_column (str): The name of the column containing SNP dictionaries.
        - mean_column (str): The name of the column to update with the calculated mean.

        Returns:
        - pd.DataFrame: The updated DataFrame.
        """
        protein_name_list = df[protein_name_column].tolist()

        for protein_name in protein_name_list:
            snp_scores_dictionary = df.loc[df[protein_name_column] == protein_name, dictionary_column].values[0]
            mean = LopoBaseline.calculate_protein_mean(protein_snp_dict=snp_scores_dictionary)
            df = update_value_based_on_protein_name(df=df, protein_name=protein_name, column_name=mean_column,
                                                    value=mean)

        return df

    @staticmethod
    def leave_one_protein_out_lists(protein_name_list):
        """
        Generate lists of proteins with one protein left out for each iteration.

        Parameters:
        - protein_name_list (list): The list of protein names.

        Returns:
        - list of lists: A list containing lists of proteins with one protein left out.
        """
        lopo_lists = []
        for i in range(len(protein_name_list)):
            lopo_list = protein_name_list[:i] + protein_name_list[i + 1:]
            lopo_lists.append(lopo_list)
        return lopo_lists

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
    def get_lopo_mean(df, id_column, value_column, lopo_list):
        """
        Calculate the mean of SNP values for each protein and update the mean column.
        input: df (pd.DataFrame): The DataFrame to update.
        output: pd.DataFrame: The updated DataFrame.
        """
        mean_values = []
        for protein in lopo_list:
            protein_mean = get_value_from_dataframe(df, id_column, protein, value_column)
            mean_values.append(protein_mean)

        return sum(mean_values) / len(mean_values)
    @staticmethod
    def get_PSSM_dict_for_protein(lopo_list, df,
                                  single_letter_amino_acids=AMINO_ACIDS_SINGLE_LETTER,
                                  mave_gs_id_column=COLUMN_NAME_OF_MAVE_GOLD_STANDARD_ID,
                                  mave_gs_snp_dict_column_name=COLUMN_NAME_OF_MAVE_GOLD_STANDARD_SNP_DICTIONARY):
        """
        This function calculates the mean of the PSSM values for each amino acid for a protein.
        Input: protein_name, lopo_list, df
        Output: amino_acid_pssm_dict
        """
        amino_acid_pssm_dict = generate_amino_pssm_dict(single_letter_amino_acids)
        for lopo_protein in lopo_list:
            snp_dict = get_value_from_dataframe(df=df,
                                     id_column=mave_gs_id_column,
                                     id_value=lopo_protein,
                                     value_column=mave_gs_snp_dict_column_name)
            for key, value in snp_dict.items():
                pssm_key = remove_digits_from_key(key)
                amino_acid_pssm_dict[pssm_key].append(value)

        for key in amino_acid_pssm_dict:
            mean_value = calculate_mean(amino_acid_pssm_dict[key])
            amino_acid_pssm_dict[key] = mean_value

        return amino_acid_pssm_dict

    @staticmethod
    def add_pssm_column_to_df(df, id_column_name):
        """
        This function adds the PSSM column to the DataFrame.
        input: df, protein_name, amino_acid_pssm_dict, id_column_name
        output: df
        """

        protein_list = get_protein_name_list(df, id_column_name)
        lopo_dict = LopoBaseline.get_lopo_dict(protein_list)

        for protein_name, lopo_list in lopo_dict.items():
            amino_acid_pssm_dict = LopoBaseline.get_PSSM_dict_for_protein(lopo_list=lopo_list, df=df)
            df = update_value_based_on_protein_name(df=df,
                                                    column_name="PSSM",
                                                    value=amino_acid_pssm_dict,
                                                    protein_name=protein_name,
                                                    id_column_name=id_column_name)

        return df


    @staticmethod
    def get_pssm_predictions(snp, pssm_dict):
        """
        This function returns the PSSM value for a given SNP.
        input: snp, pssm_dict
        output: pssm_value
        """
        # if isinstance(pssm_dict, str):
        #     pssm_dict = json.loads(pssm_dict)
        pssm_key = remove_digits_from_key(snp)
        pssm_value = pssm_dict[remove_digits_from_key(pssm_key)]
        return pssm_value


    def get_predicted_pssm_snp_scores(df, id_column_name, protein_name):
        """
        This function adds the PSSM predictions to the DataFrame.

        """
        pssm_dict = get_value_from_dataframe(df,
                                                 id_column=id_column_name,
                                                 id_value=protein_name,
                                                 value_column="PSSM")



        snp_dict = get_value_from_dataframe(df,
                                             id_column=COLUMN_NAME_OF_MAVE_GOLD_STANDARD_ID,
                                             id_value=protein_name,
                                             value_column=COLUMN_NAME_OF_MAVE_GOLD_STANDARD_SNP_DICTIONARY)




        predicted_snp_dict = {key: LopoBaseline.get_pssm_predictions(snp=key, pssm_dict=pssm_dict)
                       for key in snp_dict.keys()}

        return predicted_snp_dict


    @staticmethod
    def add_pssm_predictions_to_df(df, id_column_name):
        """
        This function adds the PSSM predictions to the DataFrame.
        """
        protein_list = get_protein_name_list(df, id_column_name)

        for protein_name in protein_list:
            predic_dict = LopoBaseline.get_predicted_pssm_snp_scores(df,
                                                                 COLUMN_NAME_OF_MAVE_GOLD_STANDARD_ID, protein_name)


            df = update_value_based_on_protein_name(df=df,
                                                    column_name="baseline_score",
                                                    value=predic_dict,
                                                    protein_name=protein_name,
                                                    id_column_name=id_column_name)

        return df


