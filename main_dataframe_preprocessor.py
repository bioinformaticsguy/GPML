import pandas as pd

from src.constants import COLUMN_NAMES_OF_MAVE_GS_DATAFRAME_LIST, TRAINING_FLAG_SUFFIX, MUTEPRED_TOOL_NAME, \
    TOOL_SCORE_COLUMN_SUFFIX, AMINO_ACID_SEQUENCE_COLUMN_NAME, DBNSFP_SNP_COLUMN_NAME, MAVE_DATAFRAME_PICKLE_FILE_NAME, \
    TOOLS_LIST, OUTPUT_DIR_DB_NSFP, MUTPRED_TRAINING_DATA_FILE_PATH, PICKLED_DATAFRAMES_DIRECTORY_PATH, \
    COLUMN_NAME_OF_MAVE_GOLD_STANDARD_PROTEIN_SEQUENCE, TRAINING_SNPS_COLUMN_SIFFIX, \
    MUTEPRED_AMINO_ACID_SUBSTITUTIONS_COLUMN_NAME, MAVE_DATAFRAME_HUMAN_LOPO_MEAN_PICKLE_FILE_NAME, \
    MAVE_DATAFRAME_HUMAN_LOPO_MEAN_PICKLE_FILE_NAME_WITH_TOOL_SCORES, \
    MAVE_DATAFRAME_ONLY_HUMAN_WITH_BASELINE_PICKLE_FILE_NAME, DEOGEN2_TRAINING_DATA_FILE_PATH, \
    DEOGEN_TRAINING_DF_COLUMNS, DEOGEN_COLUMN_NAME_TO_FILTER, DEOGEN_VALUES_TO_FILTER, DEOGEN_UNIPROT_ID_COLUMN, \
    DEOGEN_AMINO_ACID_CHANGE_COLUMN

from src.dataframe_preprocessor import MaveGoldStandard, MutepredTrainingProcessor, dbNSFPProcessor, deogen2TrainingProcessor
from main import MAVE_GS_FILE_PATH
from src.utils import pickle_dataframe, add_column_from_tool_df_to_mave_df, add_flag_column, convert_column_to_list, \
    load_dataframe, get_single_letter_point_mutation

## Path and Strings
MUTEPRED_SCORE_COLUMN_NAME = MUTEPRED_TOOL_NAME + TOOL_SCORE_COLUMN_SUFFIX
MUTEPRED_TRAINING_FLAG_COLUMN_NAME = MUTEPRED_TOOL_NAME + TRAINING_FLAG_SUFFIX
MUTEPRED_TRAINING_SNPS_COLUMN_NAME = MUTEPRED_TOOL_NAME + TRAINING_SNPS_COLUMN_SIFFIX


## Temp variables
uniprot_id = "P35520"

def get_overlaping_snp_list(deogen_training_df,
                            uniprot_id,
                            uniPprot_id_column_name=DEOGEN_UNIPROT_ID_COLUMN,
                            amino_acid_change_column_name=DEOGEN_AMINO_ACID_CHANGE_COLUMN,):

    uniprot_id_df = deogen_training_df.loc[deogen_training_df[uniPprot_id_column_name] == uniprot_id]
    snp_list = [get_single_letter_point_mutation(snp[2:]) for snp in uniprot_id_df[amino_acid_change_column_name].tolist()]

    return snp_list



if __name__ == '__main__':
    AAA = deogen2TrainingProcessor.get_deogen2_training__df(DEOGEN2_TRAINING_DATA_FILE_PATH, DEOGEN_TRAINING_DF_COLUMNS)
    AAA = deogen2TrainingProcessor.filter_unwanted_rows(AAA)

    AAAAAAspecific_uniprot_df = get_overlaping_snp_list(AAA, uniprot_id,)




    MAVE_GS_DATAFRAME = MaveGoldStandard.get_dataframe_for_mave_gs_data(mave_gs_file_path=MAVE_GS_FILE_PATH,
                                                                        column_names=COLUMN_NAMES_OF_MAVE_GS_DATAFRAME_LIST)

    MUTEPRED_DATAFRAME = MutepredTrainingProcessor.get_mutepred_df(MUTPRED_TRAINING_DATA_FILE_PATH)
    MUTEPRED_DATAFRAME = convert_column_to_list(MUTEPRED_DATAFRAME, MUTEPRED_AMINO_ACID_SUBSTITUTIONS_COLUMN_NAME)

    MAVE_GS_DATAFRAME = add_column_from_tool_df_to_mave_df(mave_df=MAVE_GS_DATAFRAME,
                                                           tool_df=MUTEPRED_DATAFRAME,
                                                           mave_df_prot_seq_col_name=COLUMN_NAME_OF_MAVE_GOLD_STANDARD_PROTEIN_SEQUENCE,
                                                           tool_df_prot_seq_col_name=AMINO_ACID_SEQUENCE_COLUMN_NAME,
                                                           tool_col_to_add=MUTEPRED_AMINO_ACID_SUBSTITUTIONS_COLUMN_NAME,
                                                           name_of_new_col=MUTEPRED_TRAINING_SNPS_COLUMN_NAME)

    # MAVE_GS_DATAFRAME = add_flag_column(df=MAVE_GS_DATAFRAME,
    #                                     target_column=MUTEPRED_TRAINING_SNPS_COLUMN_NAME,
    #                                     flag_column_name=MUTEPRED_TRAINING_FLAG_COLUMN_NAME)
    #

    # MAVE_GS_DATAFRAME = dbNSFPProcessor.add_data_from_list_of_tools(MAVE_GS_DATAFRAME,
    #                                                                 db_nsfp_output_dir_path=OUTPUT_DIR_DB_NSFP,
    #                                                                 tool_list=TOOLS_LIST,
    #                                                                 snp_column_name=DBNSFP_SNP_COLUMN_NAME)

    # pickle_dataframe(dataframe=MAVE_GS_DATAFRAME_HUMAN_WITH_BASELINE,
    #                  file_path=PICKLED_DATAFRAMES_DIRECTORY_PATH,
    #                  file_name=MAVE_DATAFRAME_ONLY_HUMAN_WITH_BASELINE_PICKLE_FILE_NAME)

    print("Debug Pause")