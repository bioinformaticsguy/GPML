from src.constants import COLUMN_NAMES_OF_MAVE_GS_DATAFRAME_LIST, TRAINING_FLAG_SUFFIX, MUTEPRED_TOOL_NAME, \
    TOOL_SCORE_COLUMN_SUFFIX, AMINO_ACID_SEQUENCE_COLUMN_NAME, SNP_COLUMN_NAME, MAVE_DATAFRAME_PICKLE_FILE_NAME, \
    TOOLS_LIST, OUTPUT_DIR_DB_NSFP, TRAINING_DATA_FILE_PATH, PICKLED_DATAFRAMES_DIRECTORY_PATH, \
    COLUMN_NAME_OF_MAVE_GOLD_STANDARD_PROTEIN_SEQUENCE, TRAINING_SNPS_COLUMN_SIFFIX, \
    MUTEPRED_AMINO_ACID_SUBSTITUTIONS_COLUMN_NAME, MAVE_DATAFRAME_HUMAN_LOPO_MEAN_PICKLE_FILE_NAME, \
    MAVE_DATAFRAME_HUMAN_LOPO_MEAN_PICKLE_FILE_NAME_WITH_TOOL_SCORES

from src.dataframe_preprocessor import MaveGoldStandard, MutepredTrainingProcessor, dbNSFPProcessor
from main import MAVE_GS_FILE_PATH
from src.utils import pickle_dataframe, add_column_from_tool_df_to_mave_df, add_flag_column, convert_column_to_list, \
    load_dataframe

## Path and Strings
MUTEPRED_SCORE_COLUMN_NAME = MUTEPRED_TOOL_NAME + TOOL_SCORE_COLUMN_SUFFIX
MUTEPRED_TRAINING_FLAG_COLUMN_NAME = MUTEPRED_TOOL_NAME + TRAINING_FLAG_SUFFIX
MUTEPRED_TRAINING_SNPS_COLUMN_NAME = MUTEPRED_TOOL_NAME + TRAINING_SNPS_COLUMN_SIFFIX
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
    MAVE_GS_DATAFRAME_HUMAN_LOPO_MEAN = load_dataframe(file_path=PICKLED_DATAFRAMES_DIRECTORY_PATH,
                                        file_name=MAVE_DATAFRAME_HUMAN_LOPO_MEAN_PICKLE_FILE_NAME)

    MUTEPRED_DATAFRAME = MutepredTrainingProcessor.get_mutepred_df(TRAINING_DATA_FILE_PATH)

    MUTEPRED_DATAFRAME = convert_column_to_list(MUTEPRED_DATAFRAME, MUTEPRED_AMINO_ACID_SUBSTITUTIONS_COLUMN_NAME)

    MAVE_GS_DATAFRAME_HUMAN_LOPO_MEAN = add_column_from_tool_df_to_mave_df(mave_df=MAVE_GS_DATAFRAME_HUMAN_LOPO_MEAN,
                                                           tool_df=MUTEPRED_DATAFRAME,
                                                           mave_df_prot_seq_col_name=COLUMN_NAME_OF_MAVE_GOLD_STANDARD_PROTEIN_SEQUENCE,
                                                           tool_df_prot_seq_col_name=AMINO_ACID_SEQUENCE_COLUMN_NAME,
                                                           tool_col_to_add=MUTEPRED_AMINO_ACID_SUBSTITUTIONS_COLUMN_NAME,
                                                           name_of_new_col=MUTEPRED_TRAINING_SNPS_COLUMN_NAME)

    MAVE_GS_DATAFRAME = add_flag_column(df=MAVE_GS_DATAFRAME_HUMAN_LOPO_MEAN,
                                        target_column=MUTEPRED_TRAINING_SNPS_COLUMN_NAME,
                                        flag_column_name=MUTEPRED_TRAINING_FLAG_COLUMN_NAME)


    MAVE_GS_DATAFRAME = add_data_from_list_of_tools(MAVE_GS_DATAFRAME,
                                                    db_nsfp_output_dir_path=OUTPUT_DIR_DB_NSFP,
                                                    tool_list=TOOLS_LIST,
                                                    snp_column_name=SNP_COLUMN_NAME)

    pickle_dataframe(dataframe=MAVE_GS_DATAFRAME,
                     file_path=PICKLED_DATAFRAMES_DIRECTORY_PATH,
                     file_name=MAVE_DATAFRAME_HUMAN_LOPO_MEAN_PICKLE_FILE_NAME_WITH_TOOL_SCORES)




    print("Debug Pause")


    # MAVE_GS_DATAFRAME = MaveGoldStandard. \
    #                         get_dataframe_for_mave_gs_data(MAVE_GS_FILE_PATH,
    #                                                        column_names=COLUMN_NAMES_OF_MAVE_GS_DATAFRAME_LIST)
    #
    # pickle_dataframe(dataframe=MAVE_GS_DATAFRAME,
    #                  file_path=PICKLED_DATAFRAMES_DIRECTORY_PATH,
    #                  file_name=MAVE_DATAFRAME_PICKLE_FILE_NAME)
    #
    # MUTEPRED_DATAFRAME = MutepredTrainingProcessor.get_mutepred_df(TRAINING_DATA_FILE_PATH)
    #
    # MUTEPRED_DATAFRAME = convert_column_to_list(MUTEPRED_DATAFRAME, MUTEPRED_AMINO_ACID_SUBSTITUTIONS_COLUMN_NAME)
    #
    #
    # MAVE_GS_DATAFRAME = add_column_from_tool_df_to_mave_df(mave_df=MAVE_GS_DATAFRAME,
    #                                                        tool_df=MUTEPRED_DATAFRAME,
    #                                                        mave_df_prot_seq_col_name=COLUMN_NAME_OF_MAVE_GOLD_STANDARD_PROTEIN_SEQUENCE,
    #                                                        tool_df_prot_seq_col_name=AMINO_ACID_SEQUENCE_COLUMN_NAME,
    #                                                        tool_col_to_add=MUTEPRED_AMINO_ACID_SUBSTITUTIONS_COLUMN_NAME,
    #                                                        name_of_new_col=MUTEPRED_TRAINING_SNPS_COLUMN_NAME)
    #
    # MAVE_GS_DATAFRAME = add_flag_column(df=MAVE_GS_DATAFRAME,
    #                                     target_column=MUTEPRED_TRAINING_SNPS_COLUMN_NAME,
    #                                     flag_column_name=MUTEPRED_TRAINING_FLAG_COLUMN_NAME)
    #
    # MAVE_GS_DATAFRAME = add_data_from_list_of_tools(MAVE_GS_DATAFRAME,
    #                                                 db_nsfp_output_dir_path=OUTPUT_DIR_DB_NSFP,
    #                                                 tool_list=TOOLS_LIST,
    #                                                 snp_column_name=SNP_COLUMN_NAME)
    #



