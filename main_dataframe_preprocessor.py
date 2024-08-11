from src.constants import LIST_OF_COL_NAMES_OF_MAVE_GS_DF, TRAINING_FLAG_SUFFIX, MUTEPRED_TOOL_NAME, \
    AMINO_ACID_SEQUENCE_COLUMN_NAME, DBNSFP_SAV_COLUMN_NAME, \
    TOOLS_LIST, OUTPUT_DIR_DB_NSFP, MUTPRED_TRAINING_DATA_FILE_PATH, \
    COL_NAME_OF_MAVE_GS_PROTEIN_SEQ, TRAINING_SAVS_COLUMN_SIFFIX, \
    MUTEPRED_AMINO_ACID_SUBSTITUTIONS_COLUMN_NAME, DEOGEN2_TRAINING_DATA_FILE_PATH, \
    DEOGEN_TRAINING_DF_COLUMNS, DEOGEN_TOOL_NAME, PICKLED_DATAFRAMES_DIRECTORY_PATH, MAVE_DATAFRAME_PICKLE_FILE_NAME, \
    CLINPRED_TOOL_NAME, PRIMATEAI_TOOL_NAME, FATHMM_TOOL_NAME, MUTATION_TASTER, COLUMN_NAME_OF_MAVE_SNPS, \
    COLUMN_NAME_OF_MAVE_GOLD_STANDARD_SAV_DICTIONARY

from src.dataframe_preprocessor import MaveGoldStandard, MutepredTrainingProcessor, dbNSFPProcessor, \
    Deogen2TrainingProcessor, ClinPredTrainingProcessor

from main import MAVE_GS_FILE_PATH
from src.utils import add_column_from_tool_df_to_mave_df, add_flag_column, convert_column_to_list, pickle_dataframe

## Path and Strings
MUTEPRED_TRAINING_FLAG_COLUMN_NAME = MUTEPRED_TOOL_NAME + TRAINING_FLAG_SUFFIX
MUTEPRED_TRAINING_SNPS_COLUMN_NAME = MUTEPRED_TOOL_NAME + TRAINING_SAVS_COLUMN_SIFFIX

DEOGEN_TRAINING_SNPS_COLUMN_NAME = DEOGEN_TOOL_NAME + TRAINING_SAVS_COLUMN_SIFFIX
DEOGEN_TRAINING_FLAG_COLUMN_NAME = DEOGEN_TOOL_NAME + TRAINING_FLAG_SUFFIX

CLINPRED_TRAINING_SNPS_COLUMN_NAME = CLINPRED_TOOL_NAME + TRAINING_SAVS_COLUMN_SIFFIX
CLINPRED_TRAINING_FLAG_COLUMN_NAME = CLINPRED_TOOL_NAME + TRAINING_FLAG_SUFFIX

PRIMATEAI_TRAINING_SNPS_COLUMN_NAME = PRIMATEAI_TOOL_NAME + TRAINING_SAVS_COLUMN_SIFFIX
PRIMATEAI_TRAINING_FLAG_COLUMN_NAME = PRIMATEAI_TOOL_NAME + TRAINING_FLAG_SUFFIX

FATHMM_TARINING_SNPS_COLUMN_NAME = FATHMM_TOOL_NAME + TRAINING_SAVS_COLUMN_SIFFIX
FATHMM_TRAINING_FLAG_COLUMN_NAME = FATHMM_TOOL_NAME + TRAINING_FLAG_SUFFIX

MUTATION_TASTER_SNPS_COLUMN_NAME = MUTATION_TASTER + TRAINING_SAVS_COLUMN_SIFFIX
MUTATION_TASTER_TRAINING_FLAG_COLUMN_NAME = MUTATION_TASTER + TRAINING_FLAG_SUFFIX

if __name__ == '__main__':


    MAVE_GS_DATAFRAME = MaveGoldStandard.get_dataframe_for_mave_gs_data(mave_gs_file_path=MAVE_GS_FILE_PATH,
                                                                        column_names=LIST_OF_COL_NAMES_OF_MAVE_GS_DF)

    MAVE_GS_DATAFRAME[COLUMN_NAME_OF_MAVE_SNPS] = MAVE_GS_DATAFRAME[COLUMN_NAME_OF_MAVE_GOLD_STANDARD_SAV_DICTIONARY].apply(lambda x: len(x) if isinstance(x, dict) else 0)


    MAVE_GS_DATAFRAME = dbNSFPProcessor.add_data_from_list_of_tools(MAVE_GS_DATAFRAME,
                                                                    db_nsfp_output_dir_path=OUTPUT_DIR_DB_NSFP,
                                                                    tool_list=TOOLS_LIST,
                                                                    snp_column_name=DBNSFP_SAV_COLUMN_NAME)


    DEOGEN2_TRAINING_DF = Deogen2TrainingProcessor.get_deogen2_training_df(DEOGEN2_TRAINING_DATA_FILE_PATH, DEOGEN_TRAINING_DF_COLUMNS)
    DEOGEN2_TRAINING_DF = Deogen2TrainingProcessor.filter_unwanted_rows(DEOGEN2_TRAINING_DF)

    MAVE_GS_DATAFRAME = Deogen2TrainingProcessor.add_training_col_for_all_proteins(MAVE_GS_DATAFRAME,
                                                                                   DEOGEN2_TRAINING_DF,
                                                                                   DEOGEN_TRAINING_SNPS_COLUMN_NAME,)

    MAVE_GS_DATAFRAME = Deogen2TrainingProcessor.add_training_col_for_all_proteins(MAVE_GS_DATAFRAME,
                                                                                   DEOGEN2_TRAINING_DF,
                                                                                   FATHMM_TARINING_SNPS_COLUMN_NAME,)


    CLINVAR_TRAINING_DF = ClinPredTrainingProcessor.get_clinvar_df()

    MAVE_GS_DATAFRAME = ClinPredTrainingProcessor.add_training_col_for_all_proteins(mave_gs_df=MAVE_GS_DATAFRAME,
                                          clinphred_df=CLINVAR_TRAINING_DF,
                                          clinphred_sav_column_name=CLINPRED_TRAINING_SNPS_COLUMN_NAME,)

    MAVE_GS_DATAFRAME = ClinPredTrainingProcessor.add_training_col_for_all_proteins(mave_gs_df=MAVE_GS_DATAFRAME,
                                          clinphred_df=CLINVAR_TRAINING_DF,
                                          clinphred_sav_column_name=PRIMATEAI_TRAINING_SNPS_COLUMN_NAME,)


    MAVE_GS_DATAFRAME = ClinPredTrainingProcessor.add_training_col_for_all_proteins(mave_gs_df=MAVE_GS_DATAFRAME,
                                          clinphred_df=CLINVAR_TRAINING_DF,
                                          clinphred_sav_column_name=MUTATION_TASTER_SNPS_COLUMN_NAME,)



    MUTPRED_TRAINING_DF = MutepredTrainingProcessor.get_mutepred_df(MUTPRED_TRAINING_DATA_FILE_PATH)
    MUTPRED_TRAINING_DF = convert_column_to_list(MUTPRED_TRAINING_DF, MUTEPRED_AMINO_ACID_SUBSTITUTIONS_COLUMN_NAME)

    MAVE_GS_DATAFRAME = add_column_from_tool_df_to_mave_df(mave_df=MAVE_GS_DATAFRAME,
                                                           tool_df=MUTPRED_TRAINING_DF,
                                                           mave_df_prot_seq_col_name=COL_NAME_OF_MAVE_GS_PROTEIN_SEQ,
                                                           tool_df_prot_seq_col_name=AMINO_ACID_SEQUENCE_COLUMN_NAME,
                                                           tool_col_to_add=MUTEPRED_AMINO_ACID_SUBSTITUTIONS_COLUMN_NAME,
                                                           name_of_new_col=MUTEPRED_TRAINING_SNPS_COLUMN_NAME)


    MAVE_GS_DATAFRAME = add_flag_column(df=MAVE_GS_DATAFRAME,
                                 target_column=MUTEPRED_TRAINING_SNPS_COLUMN_NAME,
                                 flag_column_name=MUTEPRED_TRAINING_FLAG_COLUMN_NAME)

    MAVE_GS_DATAFRAME = add_flag_column(df=MAVE_GS_DATAFRAME,
                                 target_column=CLINPRED_TRAINING_SNPS_COLUMN_NAME,
                                 flag_column_name=CLINPRED_TRAINING_FLAG_COLUMN_NAME)

    MAVE_GS_DATAFRAME = add_flag_column(df=MAVE_GS_DATAFRAME,
                                 target_column=DEOGEN_TRAINING_SNPS_COLUMN_NAME,
                                 flag_column_name=DEOGEN_TRAINING_FLAG_COLUMN_NAME)

    MAVE_GS_DATAFRAME = add_flag_column(df=MAVE_GS_DATAFRAME,
                                 target_column=PRIMATEAI_TRAINING_SNPS_COLUMN_NAME,
                                 flag_column_name=PRIMATEAI_TRAINING_FLAG_COLUMN_NAME)

    MAVE_GS_DATAFRAME = add_flag_column(df=MAVE_GS_DATAFRAME,
                                 target_column=FATHMM_TARINING_SNPS_COLUMN_NAME,
                                 flag_column_name=FATHMM_TRAINING_FLAG_COLUMN_NAME)

    MAVE_GS_DATAFRAME = add_flag_column(df=MAVE_GS_DATAFRAME,
                                 target_column=MUTATION_TASTER_SNPS_COLUMN_NAME,
                                 flag_column_name=MUTATION_TASTER_TRAINING_FLAG_COLUMN_NAME)



    pickle_dataframe(dataframe=MAVE_GS_DATAFRAME,
                     file_path=PICKLED_DATAFRAMES_DIRECTORY_PATH,
                     file_name=MAVE_DATAFRAME_PICKLE_FILE_NAME)

    print("Debug Pause")