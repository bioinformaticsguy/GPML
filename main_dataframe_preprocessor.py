from pathlib import Path
from src.constants import COLUMN_NAMES_OF_MAVE_GS_DATAFRAME_LIST, TRAINING_FLAG_SUFFIX, MUTEPRED_TOOL_NAME, \
    TOOL_SCORE_COLUMN_SUFFIX, AMINO_ACID_SEQUENCE_COLUMN_NAME, SNP_COLUMN_NAME, MAVE_DATAFRAME_PICKLE_FILE_NAME, \
    TOOLS_LIST, OUTPUT_DIR_DB_NSFP, TRAINING_DATA_FILE_PATH, PICKLED_DATAFRAMES_DIRECTORY_PATH, \
    COLUMN_NAME_OF_MAVE_GOLD_STANDARD_ID

from src.dataframe_preprocessor import MaveGoldStandard, MutepredTrainingProcessor, dbNSFPProcessor
from main import MAVE_GS_FILE_PATH
from src.utils import pickle_dataframe, merge_and_add_column

## Path and Strings
MUTEPRED_SCORE_COLUMN_NAME = MUTEPRED_TOOL_NAME + TOOL_SCORE_COLUMN_SUFFIX
MUTEPRED_TRAINING_FLAG_COLUMN_NAME = MUTEPRED_TOOL_NAME + TRAINING_FLAG_SUFFIX

def add_data_from_list_of_tools(mave_gs_dataframe,
                                db_nsfp_output_dir_path=OUTPUT_DIR_DB_NSFP,
                                tool_list=TOOLS_LIST,
                                snp_column_name=SNP_COLUMN_NAME):
    for tool in tool_list:
        dbNSFPProcessor.add_tool_score_column(mave_gs_dataframe,
                                              db_nsfp_output_dir_path,
                                              tool,
                                              snp_column_name)

    return mave_gs_dataframe

if __name__ == '__main__':
    MAVE_GS_DATAFRAME = MaveGoldStandard. \
                            get_dataframe_for_mave_gs_data(MAVE_GS_FILE_PATH,
                                                           column_names=COLUMN_NAMES_OF_MAVE_GS_DATAFRAME_LIST)

    MUTEPRED_DATAFRAME = MutepredTrainingProcessor.get_mutepred_df(TRAINING_DATA_FILE_PATH)

    FILTERED_MAVE_GOLDSTANDARD_MUTEPRED_TRAINING = MutepredTrainingProcessor. \
        filter_data_with_common_sequences(dataframe_to_filter=MAVE_GS_DATAFRAME,
                                          dataframe_with_common_values=MUTEPRED_DATAFRAME,
                                          df_sequence_column_name=AMINO_ACID_SEQUENCE_COLUMN_NAME)

    MAVE_GS_DATAFRAME = MaveGoldStandard. \
                        mark_rows_present_in_subset(superset_df=MAVE_GS_DATAFRAME,
                                                    subset_df=FILTERED_MAVE_GOLDSTANDARD_MUTEPRED_TRAINING,
                                                    new_column_name=MUTEPRED_TRAINING_FLAG_COLUMN_NAME)


    MAVE_GS_DATAFRAME = add_data_from_list_of_tools(MAVE_GS_DATAFRAME,
                                                    db_nsfp_output_dir_path=OUTPUT_DIR_DB_NSFP,
                                                    tool_list=TOOLS_LIST,
                                                    snp_column_name=SNP_COLUMN_NAME)

    import pandas as pd



    def merge_and_add_column(main_df, other_df, main_id_col, other_id_col, new_col_name):
        merged_df = pd.merge(main_df, other_df, left_on=main_id_col, right_on=other_id_col, how='left')
        main_df[new_col_name] = merged_df['SNPS']
        return main_df


    # Example usage:
    main_df = pd.DataFrame({'ID_main': [1, 2, 3, 4],
                            'Other_Column': ['A', 'B', 'C', 'D']})

    other_df = pd.DataFrame({'ID_other': [2, 4],
                             'SNPS': ['X', 'Y']})

    main_df = merge_and_add_column(main_df, other_df, 'ID', 'ID_in_other_df', 'Merged_SNPS')

    print(main_df)


    MAVE_GS_DATAFRAME = merge_and_add_column(main_df=MAVE_GS_DATAFRAME,
                                             other_df=MUTEPRED_DATAFRAME,
                                             main_id_column=COLUMN_NAME_OF_MAVE_GOLD_STANDARD_ID,
                                             other_id_column=AMINO_ACID_SEQUENCE_COLUMN_NAME,
                                             new_column="testing_SNPs")

    pickle_dataframe(dataframe=MAVE_GS_DATAFRAME,
                     file_path=PICKLED_DATAFRAMES_DIRECTORY_PATH,
                     file_name=MAVE_DATAFRAME_PICKLE_FILE_NAME)

