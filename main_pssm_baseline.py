import copy

from src.constants import PICKLED_DATAFRAMES_DIRECTORY_PATH, MAVE_DATAFRAME_PICKLE_FILE_NAME, SPECIE_NAME_HUMAN, \
    COLUMN_NAME_OF_MAVE_GOLD_STANDARD_SPECIES, COLUMN_NAME_OF_MAVE_GOLD_STANDARD_SNP_DICTIONARY, \
    COLUMN_NAME_OF_MAVE_GOLD_STANDARD_ID, COLUMN_NAME_OF_BASELINE_SCORES_DICTIONARY, \
    MAVE_DATAFRAME_ONLY_HUMAN_WITH_BASELINE_PICKLE_FILE_NAME
from src.utils import load_dataframe, filter_dataframe_by_species, get_protein_list, pickle_dataframe
from src.pssm_baseline import pssmBaseline

if __name__ == '__main__':
    GOLD_STD_DF = load_dataframe(file_path=PICKLED_DATAFRAMES_DIRECTORY_PATH,
                                 file_name=MAVE_DATAFRAME_PICKLE_FILE_NAME)

    GOLD_STD_DF = filter_dataframe_by_species(df=GOLD_STD_DF,
                                              target_species=SPECIE_NAME_HUMAN,
                                              species_column=COLUMN_NAME_OF_MAVE_GOLD_STANDARD_SPECIES)

    GOLD_STD_DF[COLUMN_NAME_OF_BASELINE_SCORES_DICTIONARY] = GOLD_STD_DF[COLUMN_NAME_OF_MAVE_GOLD_STANDARD_SNP_DICTIONARY].apply(lambda x: copy.deepcopy(x))
    protein_names = get_protein_list(dataframe=GOLD_STD_DF)
    for protein_name in protein_names:
        pssmBaseline.update_all_values_in_snp_scores_dict(dataframe=GOLD_STD_DF,
                                                          id_column=COLUMN_NAME_OF_MAVE_GOLD_STANDARD_ID,
                                                          protein_name=protein_name,
                                                          dict_col_name=COLUMN_NAME_OF_BASELINE_SCORES_DICTIONARY,
                                                          default_value=False)

    pickle_dataframe(dataframe=GOLD_STD_DF,
                     file_path=PICKLED_DATAFRAMES_DIRECTORY_PATH,
                     file_name=MAVE_DATAFRAME_ONLY_HUMAN_WITH_BASELINE_PICKLE_FILE_NAME)

    print('Debug Finish Poiting Here!')
