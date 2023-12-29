from main_dataframe_preprocessor import MUTEPRED_SCORE_COLUMN_NAME
from src.dataframe_preprocessor import COLUMN_NAME_OF_MAVE_GOLD_STANDARD_ID, \
    COLUMN_NAME_OF_MAVE_GOLD_STANDARD_SNP_DICTIONARY
from src.utils import get_mave_tool_scores_dataframe, get_correlation_and_percentage_used

PEARSON_CORELATION_SUFFIX = "_pearson_correlation"
USED_SNP_PERCENTAGE_SUFFIX = "_used_snp_percentage"

class CorelationUpdator:
    @staticmethod
    def add_tool_corelation_and_snp_percentage_column(mave_goldstandard_df,
                                                      tool_name,
                                                      mave_df_id_column_name=COLUMN_NAME_OF_MAVE_GOLD_STANDARD_ID,
                                                      mave_score_dictionary_column_name=COLUMN_NAME_OF_MAVE_GOLD_STANDARD_SNP_DICTIONARY):

        tool_pearson_score_column_name = tool_name + PEARSON_CORELATION_SUFFIX
        tool_snps_percentage_column_name = tool_name + USED_SNP_PERCENTAGE_SUFFIX

        mave_goldstandard_df[tool_pearson_score_column_name] = None
        mave_goldstandard_df[tool_snps_percentage_column_name] = None

        protein_names = mave_goldstandard_df[mave_df_id_column_name].tolist()

        for protein_name in protein_names:
            mave_score_dict = \
                mave_goldstandard_df.loc[
                    mave_goldstandard_df[mave_df_id_column_name] == protein_name, COLUMN_NAME_OF_MAVE_GOLD_STANDARD_SNP_DICTIONARY].values[0]

            tool_score_dict = \
                mave_goldstandard_df.loc[
                    mave_goldstandard_df[mave_df_id_column_name] == protein_name, MUTEPRED_SCORE_COLUMN_NAME].values[0]

            mave_tool_scores_df = get_mave_tool_scores_dataframe(mave_score_dict,
                                                                     tool_score_dict,
                                                                     mave_score_dictionary_column_name=mave_score_dictionary_column_name,
                                                                     tool_score_dictionary_column_name=MUTEPRED_SCORE_COLUMN_NAME)

            used_rows_percentage, pearson_correlation = get_correlation_and_percentage_used(mave_tool_scores_df,
                                                                                           COLUMN_NAME_OF_MAVE_GOLD_STANDARD_SNP_DICTIONARY,
                                                                                           MUTEPRED_SCORE_COLUMN_NAME)

            # Conditional indexing to assign values based on the specific name
            mave_goldstandard_df.loc[mave_goldstandard_df[mave_df_id_column_name] == protein_name, tool_pearson_score_column_name] = pearson_correlation
            mave_goldstandard_df.loc[mave_goldstandard_df[mave_df_id_column_name] == protein_name, tool_snps_percentage_column_name] = used_rows_percentage

        return mave_goldstandard_df