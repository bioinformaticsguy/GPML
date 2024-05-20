import copy

from src.constants import PICKLED_DATAFRAMES_DIRECTORY_PATH, MAVE_DATAFRAME_PICKLE_FILE_NAME, SPECIE_NAME_HUMAN, \
    COLUMN_NAME_OF_MAVE_GOLD_STANDARD_SPECIES, COLUMN_NAME_OF_MAVE_GOLD_STANDARD_SNP_DICTIONARY, \
    COLUMN_NAME_OF_MAVE_GOLD_STANDARD_ID
from src.utils import load_dataframe, filter_dataframe_by_species, get_value_from_dataframe, get_protein_list, \
    update_value_based_on_protein_name, remove_digits_from_key
from src.pssm_baseline import pssmBaseline

if __name__ == '__main__':
    GOLD_STD_DF = load_dataframe(file_path=PICKLED_DATAFRAMES_DIRECTORY_PATH,
                                 file_name=MAVE_DATAFRAME_PICKLE_FILE_NAME)

    GOLD_STD_DF = filter_dataframe_by_species(df=GOLD_STD_DF,
                                              target_species=SPECIE_NAME_HUMAN,
                                              species_column=COLUMN_NAME_OF_MAVE_GOLD_STANDARD_SPECIES)

    column_data = get_protein_list(dataframe=GOLD_STD_DF)
    lopo_dict = pssmBaseline.get_lopo_dict(column_data)

    protein_name = 'CBS_urn:mavedb:00000005-a'
    lopo_list = lopo_dict[protein_name]

    protein_pssm = pssmBaseline.get_protein_pssm(dataframe=GOLD_STD_DF, lopo_list=lopo_list)

    GOLD_STD_DF['baseline'] = GOLD_STD_DF[COLUMN_NAME_OF_MAVE_GOLD_STANDARD_SNP_DICTIONARY].apply(lambda x: copy.deepcopy(x))

    baseline_dict = dict.fromkeys(get_value_from_dataframe(df=GOLD_STD_DF,
                                             id_column=COLUMN_NAME_OF_MAVE_GOLD_STANDARD_ID,
                                             id_value=protein_name,
                                             value_column='baseline'), 0)




    for key in baseline_dict.keys():
        # print(key[0], key[-1], get_mean_from_pssm(key, protein_pssm))
        # print(GOLD_STD_DF.loc[GOLD_STD_DF[COLUMN_NAME_OF_MAVE_GOLD_STANDARD_ID] == protein_name, 'baseline'].values[0][key])
        print(GOLD_STD_DF.loc[GOLD_STD_DF[COLUMN_NAME_OF_MAVE_GOLD_STANDARD_ID] == protein_name, 'baseline'].values[0][
                  key])
        print(key, get_mean_from_pssm(remove_digits_from_key(key), protein_pssm))
        GOLD_STD_DF.loc[GOLD_STD_DF[COLUMN_NAME_OF_MAVE_GOLD_STANDARD_ID] == protein_name, 'baseline'].values[0][
            key] =  get_mean_from_pssm(remove_digits_from_key(key), protein_pssm)

        # print(GOLD_STD_DF.loc[GOLD_STD_DF[COLUMN_NAME_OF_MAVE_GOLD_STANDARD_ID] == protein_name, 'baseline'].values[0][key])





    # accessed_value = \
    # GOLD_STD_DF.loc[GOLD_STD_DF[COLUMN_NAME_OF_MAVE_GOLD_STANDARD_ID] == protein_name, 'baseline'].values[0]

    # GOLD_STD_DF.loc[GOLD_STD_DF[COLUMN_NAME_OF_MAVE_GOLD_STANDARD_ID] == protein_name, 'baseline'].values[0]["V90Y"] = 10000





    # aa = get_mean_from_pssm("AA", protein_pssm)
    # av = get_mean_from_pssm("AV", protein_pssm)



    # GOLD_STD_DF = update_value_based_on_protein_name(df=GOLD_STD_DF,
    #                                    column_name="baseline",
    #                                    value=baseline_dict,
    #                                    protein_name=protein_name,
    #                                    id_column_name=COLUMN_NAME_OF_MAVE_GOLD_STANDARD_ID)


    print('Debug Finish Poiting Here!')
