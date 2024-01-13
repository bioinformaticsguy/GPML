from main_dataframe_preprocessor import MUTEPRED_SCORE_COLUMN_NAME
from src.constants import PEARSON_CORELATION_SUFFIX, USED_SNP_PERCENTAGE_SUFFIX, TRAINING_FLAG_SUFFIX, TOOLS_LIST
from src.dataframe_preprocessor import COLUMN_NAME_OF_MAVE_GOLD_STANDARD_ID, \
    COLUMN_NAME_OF_MAVE_GOLD_STANDARD_SNP_DICTIONARY
from src.utils import get_mave_tool_scores_dataframe, get_correlation_and_percentage_used


class CorelationUpdator:
    @staticmethod
    def add_tool_correlation_and_snp_percentage_column(mave_goldstandard_df,
                                                       tool_name,
                                                       mave_df_id_column_name=COLUMN_NAME_OF_MAVE_GOLD_STANDARD_ID,
                                                       mave_score_dict_column_name=COLUMN_NAME_OF_MAVE_GOLD_STANDARD_SNP_DICTIONARY):
        """
        Adds columns for tool-specific Pearson correlation and SNP usage percentage to the MAVE gold standard DataFrame.

        Parameters:
        - mave_goldstandard_df (pd.DataFrame): MAVE gold standard DataFrame.
        - tool_name (str): Name of the tool for which scores are being calculated.
        - mave_df_id_column_name (str, optional): Column name containing protein identifiers in the MAVE DataFrame.
          Defaults to COLUMN_NAME_OF_MAVE_GOLD_STANDARD_ID.
        - mave_score_dict_column_name (str, optional): Column name containing MAVE SNP dictionaries in the MAVE DataFrame.
          Defaults to COLUMN_NAME_OF_MAVE_GOLD_STANDARD_SNP_DICTIONARY.

        Returns:
        - pd.DataFrame: MAVE gold standard DataFrame with added tool-specific columns.
        """
        # Define column names for tool-specific scores
        tool_pearson_score_column = tool_name + PEARSON_CORELATION_SUFFIX
        tool_snps_percentage_column = tool_name + USED_SNP_PERCENTAGE_SUFFIX

        # Initialize new columns with None
        mave_goldstandard_df[tool_pearson_score_column] = None
        mave_goldstandard_df[tool_snps_percentage_column] = None

        # Extract protein names from the DataFrame
        protein_names = mave_goldstandard_df[mave_df_id_column_name].tolist()

        # Iterate through each protein to calculate tool-specific scores
        for protein_name in protein_names:
            # Extract MAVE and tool-specific scores dictionaries for the current protein
            mave_score_dict = mave_goldstandard_df.loc[
                mave_goldstandard_df[mave_df_id_column_name] == protein_name, mave_score_dict_column_name].values[0]
            tool_score_dict = mave_goldstandard_df.loc[
                mave_goldstandard_df[mave_df_id_column_name] == protein_name, MUTEPRED_SCORE_COLUMN_NAME].values[0]

            if tool_score_dict == None:
                continue

            # Create a DataFrame for MAVE and tool scores
            mave_tool_scores_df = get_mave_tool_scores_dataframe(mave_score_dict,
                                                                 tool_score_dict,
                                                                 mave_score_dictionary_column_name=mave_score_dict_column_name,
                                                                 tool_score_dictionary_column_name=MUTEPRED_SCORE_COLUMN_NAME)

            # Calculate Pearson correlation and percentage of used rows
            used_rows_percentage, pearson_correlation = get_correlation_and_percentage_used(
                mave_tool_scores_df,
                COLUMN_NAME_OF_MAVE_GOLD_STANDARD_SNP_DICTIONARY,
                MUTEPRED_SCORE_COLUMN_NAME)

            # Assign calculated values to the tool-specific columns based on the protein name
            mave_goldstandard_df.loc[mave_goldstandard_df[mave_df_id_column_name] == \
                                     protein_name, tool_pearson_score_column] = pearson_correlation
            mave_goldstandard_df.loc[mave_goldstandard_df[mave_df_id_column_name] == \
                                     protein_name, tool_snps_percentage_column] = used_rows_percentage

        return mave_goldstandard_df

    @staticmethod
    def calculate_tool_bias(df_with_spearman_scores,
                            tool_name,
                            pearson_corelation_suffix=PEARSON_CORELATION_SUFFIX,
                            training_flag_suffix=TRAINING_FLAG_SUFFIX):
        """
        Calculates the bias of a tool based on Pearson correlation scores.

        Parameters:
        - df_with_spearman_scores (pd.DataFrame): DataFrame containing Spearman correlation scores.
        - tool_name (str): Name of the tool for which bias is calculated.
        - pearson_corelation_suffix (str): Suffix for the Pearson correlation column names.
        - training_flag_suffix (str): Suffix for the training flag column names.

        Returns:
        float: Bias score representing the absolute difference between overall and unbiased tool scores.
        """
        overall_score = df_with_spearman_scores[tool_name + pearson_corelation_suffix].mean()
        unbiased_score = df_with_spearman_scores[df_with_spearman_scores[tool_name + training_flag_suffix] == 0][
            tool_name + pearson_corelation_suffix].mean()
        bias = abs(unbiased_score - overall_score)
        return bias

    @staticmethod
    def add_tool_data_for_multiple_tools(mave_goldstandard_df,
                                         tools_names_list=TOOLS_LIST,
                                         mave_df_id_column_name=COLUMN_NAME_OF_MAVE_GOLD_STANDARD_ID,
                                         mave_score_dict_column_name=COLUMN_NAME_OF_MAVE_GOLD_STANDARD_SNP_DICTIONARY):
        """
        Adds columns for tool-specific Pearson correlation and SNP usage percentage to the MAVE gold standard DataFrame
        for multiple tools.

        Parameters:
        - mave_goldstandard_df (pd.DataFrame): MAVE gold standard DataFrame.
        - tool_names (list): List of tool names for which scores are being calculated.
        - mave_df_id_column_name (str, optional): Column name containing protein identifiers in the MAVE DataFrame.
          Defaults to COLUMN_NAME_OF_MAVE_GOLD_STANDARD_ID.
        - mave_score_dict_column_name (str, optional): Column name containing MAVE SNP dictionaries in the MAVE DataFrame.
          Defaults to COLUMN_NAME_OF_MAVE_GOLD_STANDARD_SNP_DICTIONARY.

        Returns:
        - pd.DataFrame: MAVE gold standard DataFrame with added tool-specific columns for all tools.
        """
        for tool_name in tools_names_list:
            mave_goldstandard_df = CorelationUpdator.add_tool_correlation_and_snp_percentage_column(
                mave_goldstandard_df=mave_goldstandard_df,
                tool_name=tool_name,
                mave_df_id_column_name=mave_df_id_column_name,
                mave_score_dict_column_name=mave_score_dict_column_name)

        return mave_goldstandard_df
