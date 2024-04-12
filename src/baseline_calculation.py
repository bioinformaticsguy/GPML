from src.utils import update_value_based_on_protein_name, get_value_from_dataframe


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