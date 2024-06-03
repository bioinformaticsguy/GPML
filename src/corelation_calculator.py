import numpy as np
from src.constants import PEARSON_CORELATION_SUFFIX, TRAINING_FLAG_SUFFIX, TOOLS_LIST, \
    COLUMN_NAME_OF_BASELINE_SCORES_DICTIONARY
from src.dataframe_preprocessor import COLUMN_NAME_OF_MAVE_GOLD_STANDARD_ID, \
    COLUMN_NAME_OF_MAVE_GOLD_STANDARD_SNP_DICTIONARY
from src.utils import get_mave_tool_scores_dataframe, get_correlation_and_percentage_used, exclude_snps, \
    generate_tool_columns, get_training_snps_column_name


class CorelationUpdator:
    @staticmethod
    def add_tool_correlation_and_snp_percentage_column(mave_goldstandard_df,
                                                       tool_name,
                                                       exclude_tool_training_snps_flag=False,
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
        tool_pearson_score_column, tool_snps_percentage_column, tool_score_column_name = \
            generate_tool_columns(tool_name=tool_name, exclude_tool_training_snps_flag=exclude_tool_training_snps_flag)

        # Initialize new columns with None
        mave_goldstandard_df[tool_pearson_score_column] = np.NaN
        mave_goldstandard_df[tool_snps_percentage_column] = np.NaN

        # Extract protein names from the DataFrame
        protein_names = mave_goldstandard_df[mave_df_id_column_name].tolist()

        # Iterate through each protein to calculate tool-specific scores
        for protein_name in protein_names:
            # Extract MAVE and tool-specific scores dictionaries for the current protein
            mave_score_dict = mave_goldstandard_df.loc[
                mave_goldstandard_df[mave_df_id_column_name] == protein_name, mave_score_dict_column_name].values[0]
            tool_score_dict = mave_goldstandard_df.loc[
                mave_goldstandard_df[mave_df_id_column_name] == protein_name, tool_score_column_name].values[0]

            if tool_score_dict == None:
                continue

            # Create a DataFrame for MAVE and tool scores
            mave_tool_scores_df = get_mave_tool_scores_dataframe(mave_score_dict,
                                                                 tool_score_dict,
                                                                 mave_score_dictionary_column_name=mave_score_dict_column_name,
                                                                 tool_score_dictionary_column_name=tool_score_column_name)

            if exclude_tool_training_snps_flag:
                training_snps_column_name = get_training_snps_column_name(tool_name)
                training_snps_list = mave_goldstandard_df.loc[mave_goldstandard_df[mave_df_id_column_name] == \
                                                              protein_name, training_snps_column_name].iloc[0]

                mave_tool_scores_df = exclude_snps(df=mave_tool_scores_df, exclude_snps_list=training_snps_list)

            # Calculate Pearson correlation and percentage of used rows
            used_rows_percentage, pearson_correlation = get_correlation_and_percentage_used(
                df=mave_tool_scores_df,
                column1_name=mave_score_dict_column_name,
                column2_name=tool_score_column_name)

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
                                         exclude_tool_training_snps_flag=False,
                                         tools_names_list=TOOLS_LIST):
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
        if exclude_tool_training_snps_flag:
            for tool_name in tools_names_list:
                mave_goldstandard_df = CorelationUpdator. \
                    add_tool_correlation_and_snp_percentage_column(mave_goldstandard_df=mave_goldstandard_df,
                                                                   tool_name=tool_name,
                                                                   exclude_tool_training_snps_flag=True)
        else:
            for tool_name in tools_names_list:
                mave_goldstandard_df = CorelationUpdator. \
                    add_tool_correlation_and_snp_percentage_column(mave_goldstandard_df=mave_goldstandard_df,
                                                                   tool_name=tool_name)

        return mave_goldstandard_df


class DeogenCorelation:
    @staticmethod
    def add_deogen_baseline_corelation(mave_gs_df,
                                       baseline_score_col_name=COLUMN_NAME_OF_BASELINE_SCORES_DICTIONARY,
                                       mave_df_id_column_name=COLUMN_NAME_OF_MAVE_GOLD_STANDARD_ID,
                                       mave_score_dict_column_name=COLUMN_NAME_OF_MAVE_GOLD_STANDARD_SNP_DICTIONARY,
                                       pearson_corelation_suffix=PEARSON_CORELATION_SUFFIX, ):


        baseline_pearson_score_col = baseline_score_col_name + pearson_corelation_suffix
        mave_gs_df[baseline_pearson_score_col] = np.NaN

        protein_names = mave_gs_df[mave_df_id_column_name].tolist()

        for protein_name in protein_names:
            # Extract MAVE and tool-specific scores dictionaries for the current protein
            mave_score_dict = mave_gs_df.loc[
                mave_gs_df[mave_df_id_column_name] == protein_name, mave_score_dict_column_name].values[0]
            baeline_score_dict = mave_gs_df.loc[
                mave_gs_df[mave_df_id_column_name] == protein_name, baseline_score_col_name].values[0]

            mave_tool_scores_df = get_mave_tool_scores_dataframe(mave_score_dict,
                                                                 baeline_score_dict,
                                                                 mave_score_dictionary_column_name=mave_score_dict_column_name,
                                                                 tool_score_dictionary_column_name=baseline_score_col_name)

            used_rows_percentage, pearson_correlation = get_correlation_and_percentage_used(
                df=mave_tool_scores_df,
                column1_name=mave_score_dict_column_name,
                column2_name=baseline_score_col_name)

            mave_gs_df.loc[mave_gs_df[mave_df_id_column_name] == \
                                    protein_name, baseline_pearson_score_col] = pearson_correlation

        return mave_gs_df
