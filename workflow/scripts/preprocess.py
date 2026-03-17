"""
preprocess.py — Step 1

Builds the master DataFrame from the MAVE gold standard FASTA, integrates
dbNSFP scores for all 9 prediction tools, and flags mutations that appear in
the training sets of DEOGEN2, ClinVar, and MutPred.

Output: Data/pickled_dataframes/gold_std_df.pkl
"""

from pathlib import Path

from src.constants import (
    LIST_OF_COL_NAMES_OF_MAVE_GS_DF,
    TOOLS_LIST,
    OUTPUT_DIR_DB_NSFP,
    DBNSFP_SAV_COLUMN_NAME,
    MUTPRED_TRAINING_DATA_FILE_PATH,
    DEOGEN2_TRAINING_DATA_FILE_PATH,
    DEOGEN_TRAINING_DF_COLUMNS,
    COL_NAME_OF_MAVE_GS_PROTEIN_SEQ,
    AMINO_ACID_SEQUENCE_COLUMN_NAME,
    MUTEPRED_TOOL_NAME,
    DEOGEN_TOOL_NAME,
    CLINPRED_TOOL_NAME,
    PRIMATEAI_TOOL_NAME,
    FATHMM_TOOL_NAME,
    MUTATION_TASTER,
    TRAINING_FLAG_SUFFIX,
    TRAINING_SAVS_COLUMN_SIFFIX,
    MUTEPRED_AMINO_ACID_SUBSTITUTIONS_COLUMN_NAME,
    COLUMN_NAME_OF_MAVE_SNPS,
    COLUMN_NAME_OF_MAVE_GOLD_STANDARD_SAV_DICTIONARY,
    PICKLED_DATAFRAMES_DIRECTORY_PATH,
    MAVE_DATAFRAME_PICKLE_FILE_NAME,
)
from src.dataframe_preprocessor import (
    MaveGoldStandard,
    MutepredTrainingProcessor,
    dbNSFPProcessor,
    Deogen2TrainingProcessor,
    ClinPredTrainingProcessor,
)
from src.utils import (
    add_column_from_tool_df_to_mave_df,
    add_flag_column,
    convert_column_to_list,
    pickle_dataframe,
)

MAVE_GS_FILE_PATH = Path("Data/mave_gs_data/mave_db_gold_standard.fasta")

# Training column names derived from constants — no cross-script imports needed
MUTEPRED_TRAINING_SAVS_COL  = MUTEPRED_TOOL_NAME  + TRAINING_SAVS_COLUMN_SIFFIX
MUTEPRED_TRAINING_FLAG_COL  = MUTEPRED_TOOL_NAME  + TRAINING_FLAG_SUFFIX
DEOGEN_TRAINING_SAVS_COL    = DEOGEN_TOOL_NAME    + TRAINING_SAVS_COLUMN_SIFFIX
DEOGEN_TRAINING_FLAG_COL    = DEOGEN_TOOL_NAME    + TRAINING_FLAG_SUFFIX
CLINPRED_TRAINING_SAVS_COL  = CLINPRED_TOOL_NAME  + TRAINING_SAVS_COLUMN_SIFFIX
CLINPRED_TRAINING_FLAG_COL  = CLINPRED_TOOL_NAME  + TRAINING_FLAG_SUFFIX
PRIMATEAI_TRAINING_SAVS_COL = PRIMATEAI_TOOL_NAME + TRAINING_SAVS_COLUMN_SIFFIX
PRIMATEAI_TRAINING_FLAG_COL = PRIMATEAI_TOOL_NAME + TRAINING_FLAG_SUFFIX
FATHMM_TRAINING_SAVS_COL    = FATHMM_TOOL_NAME    + TRAINING_SAVS_COLUMN_SIFFIX
FATHMM_TRAINING_FLAG_COL    = FATHMM_TOOL_NAME    + TRAINING_FLAG_SUFFIX
MUTATION_TASTER_TRAINING_SAVS_COL = MUTATION_TASTER + TRAINING_SAVS_COLUMN_SIFFIX
MUTATION_TASTER_TRAINING_FLAG_COL = MUTATION_TASTER + TRAINING_FLAG_SUFFIX


if __name__ == "__main__":

    # --- Build base DataFrame from MAVE gold standard ---
    df = MaveGoldStandard.get_dataframe_for_mave_gs_data(
        mave_gs_file_path=MAVE_GS_FILE_PATH,
        column_names=LIST_OF_COL_NAMES_OF_MAVE_GS_DF,
    )
    df[COLUMN_NAME_OF_MAVE_SNPS] = df[COLUMN_NAME_OF_MAVE_GOLD_STANDARD_SAV_DICTIONARY].apply(
        lambda x: len(x) if isinstance(x, dict) else 0
    )

    # --- Integrate dbNSFP scores for all 9 tools ---
    df = dbNSFPProcessor.add_data_from_list_of_tools(
        df,
        db_nsfp_output_dir_path=OUTPUT_DIR_DB_NSFP,
        tool_list=TOOLS_LIST,
        snp_column_name=DBNSFP_SAV_COLUMN_NAME,
    )

    # --- DEOGEN2 training flags (also used for FATHMM — shared training source) ---
    deogen2_training_df = Deogen2TrainingProcessor.get_deogen2_training_df(
        DEOGEN2_TRAINING_DATA_FILE_PATH, DEOGEN_TRAINING_DF_COLUMNS
    )
    deogen2_training_df = Deogen2TrainingProcessor.filter_unwanted_rows(deogen2_training_df)

    df = Deogen2TrainingProcessor.add_training_col_for_all_proteins(
        df, deogen2_training_df, DEOGEN_TRAINING_SAVS_COL
    )
    df = Deogen2TrainingProcessor.add_training_col_for_all_proteins(
        df, deogen2_training_df, FATHMM_TRAINING_SAVS_COL
    )

    # --- ClinVar training flags (shared by ClinPred, PrimateAI, MutationTaster) ---
    clinvar_training_df = ClinPredTrainingProcessor.get_clinvar_df()

    df = ClinPredTrainingProcessor.add_training_col_for_all_proteins(
        mave_gs_df=df,
        clinphred_df=clinvar_training_df,
        clinphred_sav_column_name=CLINPRED_TRAINING_SAVS_COL,
    )
    df = ClinPredTrainingProcessor.add_training_col_for_all_proteins(
        mave_gs_df=df,
        clinphred_df=clinvar_training_df,
        clinphred_sav_column_name=PRIMATEAI_TRAINING_SAVS_COL,
    )
    df = ClinPredTrainingProcessor.add_training_col_for_all_proteins(
        mave_gs_df=df,
        clinphred_df=clinvar_training_df,
        clinphred_sav_column_name=MUTATION_TASTER_TRAINING_SAVS_COL,
    )

    # --- MutPred training flags ---
    mutpred_training_df = MutepredTrainingProcessor.get_mutepred_df(MUTPRED_TRAINING_DATA_FILE_PATH)
    mutpred_training_df = convert_column_to_list(
        mutpred_training_df, MUTEPRED_AMINO_ACID_SUBSTITUTIONS_COLUMN_NAME
    )
    df = add_column_from_tool_df_to_mave_df(
        mave_df=df,
        tool_df=mutpred_training_df,
        mave_df_prot_seq_col_name=COL_NAME_OF_MAVE_GS_PROTEIN_SEQ,
        tool_df_prot_seq_col_name=AMINO_ACID_SEQUENCE_COLUMN_NAME,
        tool_col_to_add=MUTEPRED_AMINO_ACID_SUBSTITUTIONS_COLUMN_NAME,
        name_of_new_col=MUTEPRED_TRAINING_SAVS_COL,
    )

    # --- Binary training flags for all tools ---
    for savs_col, flag_col in [
        (MUTEPRED_TRAINING_SAVS_COL,        MUTEPRED_TRAINING_FLAG_COL),
        (DEOGEN_TRAINING_SAVS_COL,          DEOGEN_TRAINING_FLAG_COL),
        (CLINPRED_TRAINING_SAVS_COL,        CLINPRED_TRAINING_FLAG_COL),
        (PRIMATEAI_TRAINING_SAVS_COL,       PRIMATEAI_TRAINING_FLAG_COL),
        (FATHMM_TRAINING_SAVS_COL,          FATHMM_TRAINING_FLAG_COL),
        (MUTATION_TASTER_TRAINING_SAVS_COL, MUTATION_TASTER_TRAINING_FLAG_COL),
    ]:
        df = add_flag_column(df=df, target_column=savs_col, flag_column_name=flag_col)

    # --- Save ---
    pickle_dataframe(
        dataframe=df,
        file_path=PICKLED_DATAFRAMES_DIRECTORY_PATH,
        file_name=MAVE_DATAFRAME_PICKLE_FILE_NAME,
    )
