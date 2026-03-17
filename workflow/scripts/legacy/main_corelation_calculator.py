import numpy as np

from main_dataframe_preprocessor import DEOGEN_TRAINING_FLAG_COLUMN_NAME, CLINPRED_TRAINING_FLAG_COLUMN_NAME, \
    PRIMATEAI_TRAINING_FLAG_COLUMN_NAME, FATHMM_TRAINING_FLAG_COLUMN_NAME, MUTATION_TASTER_TRAINING_FLAG_COLUMN_NAME
from src.constants import MAVE_DATAFRAME_PICKLE_FILE_NAME, TOOLS_LIST, PICKLED_DATAFRAMES_DIRECTORY_PATH, \
    MUTEPRED_TOOL_NAME, MAVE_DATAFRAME_ONLY_HUMAN_WITH_BASELINE_PICKLE_FILE_NAME, \
    COLUMN_NAME_OF_BASELINE_SCORES_DICTIONARY, MAVE_DATAFRAME_ONLY_HUMAN_WITH_BASELINE_CORELATION_PICKLE_FILE_NAME, \
    DEOGEN_TOOL_NAME, SPEAR_COR_SUFFIX, EXCLUDE_TRAINING_SAV_SUFFIX, CLINPRED_TOOL_NAME, PRIMATEAI_TOOL_NAME, \
    FATHMM_TOOL_NAME, MUTATION_TASTER, TRAINING_FLAG_SUFFIX, STRICT_COR_SUFFIX
from src.utils import load_dataframe, pickle_dataframe, filter_dataframe_by_species
from src.corelation_calculator import CorelationUpdator, DeogenCorelation

if __name__ == '__main__':
    MAVE_GS_DATAFRAME_HUMAN_WITH_BASELINE = load_dataframe(file_path=PICKLED_DATAFRAMES_DIRECTORY_PATH,
                                                           file_name=MAVE_DATAFRAME_ONLY_HUMAN_WITH_BASELINE_PICKLE_FILE_NAME)

    MAVE_GS_DATAFRAME_HUMAN_WITH_BASELINE = CorelationUpdator. \
        add_tool_data_for_multiple_tools(MAVE_GS_DATAFRAME_HUMAN_WITH_BASELINE,
                                         tools_names_list=TOOLS_LIST)

    MAVE_GS_DATAFRAME_HUMAN_WITH_BASELINE = CorelationUpdator. \
        add_tool_correlation_and_snp_percentage_column(mave_goldstandard_df=MAVE_GS_DATAFRAME_HUMAN_WITH_BASELINE,
                                                       tool_name=DEOGEN_TOOL_NAME,
                                                       exclude_tool_training_snps_flag=True)

    MAVE_GS_DATAFRAME_HUMAN_WITH_BASELINE = CorelationUpdator. \
        add_tool_correlation_and_snp_percentage_column(mave_goldstandard_df=MAVE_GS_DATAFRAME_HUMAN_WITH_BASELINE,
                                                       tool_name=MUTEPRED_TOOL_NAME,
                                                       exclude_tool_training_snps_flag=True)

    MAVE_GS_DATAFRAME_HUMAN_WITH_BASELINE = CorelationUpdator. \
        add_tool_correlation_and_snp_percentage_column(mave_goldstandard_df=MAVE_GS_DATAFRAME_HUMAN_WITH_BASELINE,
                                                       tool_name=CLINPRED_TOOL_NAME,
                                                       exclude_tool_training_snps_flag=True)

    MAVE_GS_DATAFRAME_HUMAN_WITH_BASELINE = CorelationUpdator. \
        add_tool_correlation_and_snp_percentage_column(mave_goldstandard_df=MAVE_GS_DATAFRAME_HUMAN_WITH_BASELINE,
                                                       tool_name=PRIMATEAI_TOOL_NAME,
                                                       exclude_tool_training_snps_flag=True)

    MAVE_GS_DATAFRAME_HUMAN_WITH_BASELINE = CorelationUpdator. \
        add_tool_correlation_and_snp_percentage_column(mave_goldstandard_df=MAVE_GS_DATAFRAME_HUMAN_WITH_BASELINE,
                                                       tool_name=FATHMM_TOOL_NAME,
                                                       exclude_tool_training_snps_flag=True)

    MAVE_GS_DATAFRAME_HUMAN_WITH_BASELINE = CorelationUpdator. \
        add_tool_correlation_and_snp_percentage_column(mave_goldstandard_df=MAVE_GS_DATAFRAME_HUMAN_WITH_BASELINE,
                                                       tool_name=MUTATION_TASTER,
                                                       exclude_tool_training_snps_flag=True)

    MAVE_GS_DATAFRAME_HUMAN_WITH_BASELINE = DeogenCorelation.add_deogen_baseline_corelation(MAVE_GS_DATAFRAME_HUMAN_WITH_BASELINE,)

    correlation_column = MUTEPRED_TOOL_NAME + SPEAR_COR_SUFFIX
    MUTEPRED_TRAINING_FLAG_COLUMN_NAME = MUTEPRED_TOOL_NAME + TRAINING_FLAG_SUFFIX


    def calculate_correlation(row, correlation_column, mask_column):
        if row[mask_column] == 0:
            return row[correlation_column]
        else:
            return np.nan


    # MAVE_GS_DATAFRAME_HUMAN_WITH_BASELINE['new_column_rankscore'] = MAVE_GS_DATAFRAME_HUMAN_WITH_BASELINE.apply(calculate_correlation, args=(correlation_column, MUTEPRED_TRAINING_FLAG_COLUMN_NAME), axis=1)
    MAVE_GS_DATAFRAME_HUMAN_WITH_BASELINE[MUTEPRED_TOOL_NAME + STRICT_COR_SUFFIX] = MAVE_GS_DATAFRAME_HUMAN_WITH_BASELINE.apply(calculate_correlation, args=(MUTEPRED_TOOL_NAME + SPEAR_COR_SUFFIX, MUTEPRED_TRAINING_FLAG_COLUMN_NAME), axis=1)
    MAVE_GS_DATAFRAME_HUMAN_WITH_BASELINE[DEOGEN_TOOL_NAME + STRICT_COR_SUFFIX] = MAVE_GS_DATAFRAME_HUMAN_WITH_BASELINE.apply(calculate_correlation, args=(DEOGEN_TOOL_NAME + SPEAR_COR_SUFFIX, DEOGEN_TRAINING_FLAG_COLUMN_NAME), axis=1)
    MAVE_GS_DATAFRAME_HUMAN_WITH_BASELINE[CLINPRED_TOOL_NAME + STRICT_COR_SUFFIX] = MAVE_GS_DATAFRAME_HUMAN_WITH_BASELINE.apply(calculate_correlation, args=(CLINPRED_TOOL_NAME + SPEAR_COR_SUFFIX, CLINPRED_TRAINING_FLAG_COLUMN_NAME), axis=1)
    MAVE_GS_DATAFRAME_HUMAN_WITH_BASELINE[PRIMATEAI_TOOL_NAME + STRICT_COR_SUFFIX] = MAVE_GS_DATAFRAME_HUMAN_WITH_BASELINE.apply(calculate_correlation, args=(PRIMATEAI_TOOL_NAME + SPEAR_COR_SUFFIX, PRIMATEAI_TRAINING_FLAG_COLUMN_NAME), axis=1)
    MAVE_GS_DATAFRAME_HUMAN_WITH_BASELINE[FATHMM_TOOL_NAME + STRICT_COR_SUFFIX] = MAVE_GS_DATAFRAME_HUMAN_WITH_BASELINE.apply(calculate_correlation, args=(FATHMM_TOOL_NAME + SPEAR_COR_SUFFIX, FATHMM_TRAINING_FLAG_COLUMN_NAME), axis=1)
    MAVE_GS_DATAFRAME_HUMAN_WITH_BASELINE[MUTATION_TASTER + STRICT_COR_SUFFIX] = MAVE_GS_DATAFRAME_HUMAN_WITH_BASELINE.apply(calculate_correlation, args=(MUTATION_TASTER + SPEAR_COR_SUFFIX, MUTATION_TASTER_TRAINING_FLAG_COLUMN_NAME), axis=1)



    # mean_normal = MAVE_GS_DATAFRAME_HUMAN_WITH_BASELINE[MUTEPRED_TOOL_NAME + SPEAR_COR_SUFFIX].mean()
    #
    # mean_no_bias = MAVE_GS_DATAFRAME_HUMAN_WITH_BASELINE[MUTEPRED_TOOL_NAME + SPEAR_COR_SUFFIX + EXCLUDE_TRAINING_SAV_SUFFIX].mean()
    #
    #
    # strict = CorelationUpdator.calculate_tool_bias(df_with_spearman_scores=MAVE_GS_DATAFRAME_HUMAN_WITH_BASELINE,
    #                     tool_name=MUTEPRED_TOOL_NAME)
    #
    # relax = abs(mean_no_bias - mean_normal)

    #
    # mutepred_tool_bias = CorelationUpdator. \
    #     calculate_tool_bias(df_with_spearman_scores=LOADED_MAVE_DF,
    #                         tool_name=MUTEPRED_TOOL_NAME)
    #
    pickle_dataframe(dataframe=MAVE_GS_DATAFRAME_HUMAN_WITH_BASELINE,
                     file_path=PICKLED_DATAFRAMES_DIRECTORY_PATH,
                     file_name=MAVE_DATAFRAME_ONLY_HUMAN_WITH_BASELINE_CORELATION_PICKLE_FILE_NAME)



    print("Debug Pause")
