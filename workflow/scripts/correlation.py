"""
correlation.py — Step 3

Calculates Spearman correlations between MAVE measurements and each
prediction tool. Produces both correlation_all (all SNPs) and
correlation_strict (training SNPs excluded) for the 6 tools that have
known training data.

Input  : Data/pickled_dataframes/gold_std_df_only_human_with_baseline.pkl
         Data/pickled_dataframes/pssm_base.pkl  (consumed implicitly via baseline column)
Output : Data/pickled_dataframes/gold_std_df_only_human_with_baseline_corelation.pkl
"""

import numpy as np

from src.constants import (
    PICKLED_DATAFRAMES_DIRECTORY_PATH,
    MAVE_DATAFRAME_ONLY_HUMAN_WITH_BASELINE_PICKLE_FILE_NAME,
    MAVE_DATAFRAME_ONLY_HUMAN_WITH_BASELINE_CORELATION_PICKLE_FILE_NAME,
    TOOLS_LIST,
    MUTEPRED_TOOL_NAME,
    DEOGEN_TOOL_NAME,
    CLINPRED_TOOL_NAME,
    PRIMATEAI_TOOL_NAME,
    FATHMM_TOOL_NAME,
    MUTATION_TASTER,
    SPEAR_COR_SUFFIX,
    STRICT_COR_SUFFIX,
    TRAINING_FLAG_SUFFIX,
)
from src.corelation_calculator import CorelationUpdator, DeogenCorelation
from src.utils import load_dataframe, pickle_dataframe

# Tools that have known training data and therefore need strict correlation
TOOLS_WITH_TRAINING_DATA = [
    MUTEPRED_TOOL_NAME,
    DEOGEN_TOOL_NAME,
    CLINPRED_TOOL_NAME,
    PRIMATEAI_TOOL_NAME,
    FATHMM_TOOL_NAME,
    MUTATION_TASTER,
]


def _add_strict_correlation(df, tool_name):
    """Set strict correlation to NaN for rows where the tool's training flag is set."""
    flag_col = tool_name + TRAINING_FLAG_SUFFIX
    cor_col  = tool_name + SPEAR_COR_SUFFIX
    strict_col = tool_name + STRICT_COR_SUFFIX
    df[strict_col] = df.apply(
        lambda row: row[cor_col] if row[flag_col] == 0 else np.nan, axis=1
    )
    return df


if __name__ == "__main__":

    df = load_dataframe(
        file_path=PICKLED_DATAFRAMES_DIRECTORY_PATH,
        file_name=MAVE_DATAFRAME_ONLY_HUMAN_WITH_BASELINE_PICKLE_FILE_NAME,
    )

    # --- Correlations for all 9 tools (all SNPs) ---
    df = CorelationUpdator.add_tool_data_for_multiple_tools(
        df, tools_names_list=TOOLS_LIST
    )

    # --- Correlations excluding training SNPs for tools with known training data ---
    for tool in TOOLS_WITH_TRAINING_DATA:
        df = CorelationUpdator.add_tool_correlation_and_snp_percentage_column(
            mave_goldstandard_df=df,
            tool_name=tool,
            exclude_tool_training_snps_flag=True,
        )

    # --- PSSM baseline correlation ---
    df = DeogenCorelation.add_deogen_baseline_corelation(df)

    # --- Strict correlation columns (NaN where training flag is set) ---
    for tool in TOOLS_WITH_TRAINING_DATA:
        df = _add_strict_correlation(df, tool)

    pickle_dataframe(
        dataframe=df,
        file_path=PICKLED_DATAFRAMES_DIRECTORY_PATH,
        file_name=MAVE_DATAFRAME_ONLY_HUMAN_WITH_BASELINE_CORELATION_PICKLE_FILE_NAME,
    )
