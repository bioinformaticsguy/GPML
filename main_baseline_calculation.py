from src.baseline_calculation import LopoBaseline
from src.constants import PICKLED_DATAFRAMES_DIRECTORY_PATH, MAVE_DATAFRAME_PICKLE_FILE_NAME, SPECIE_NAME_HUMAN, \
    COLUMN_NAME_OF_MAVE_GOLD_STANDARD_SNP_DICTIONARY, COLUMN_NAME_OF_MAVE_GOLD_STANDARD_SPECIES, \
    COLUMN_NAME_OF_MAVE_GOLD_STANDARD_ID, MAVE_DATAFRAME_HUMAN_LOPO_MEAN_PICKLE_FILE_NAME, AMINO_ACIDS_SINGLE_LETTER
from src.utils import load_dataframe, filter_dataframe_by_species, get_value_from_dataframe, \
    update_value_based_on_protein_name, pickle_dataframe, generate_amino_pssm_dict, remove_digits_from_key, \
    calculate_mean

if __name__ == '__main__':
    GOLD_STD_DF = load_dataframe(file_path=PICKLED_DATAFRAMES_DIRECTORY_PATH,
                                 file_name=MAVE_DATAFRAME_PICKLE_FILE_NAME)

    GOLD_STD_DF_HUMAN= filter_dataframe_by_species(df=GOLD_STD_DF,
                                target_species=SPECIE_NAME_HUMAN,
                                species_column=COLUMN_NAME_OF_MAVE_GOLD_STANDARD_SPECIES)

    protein_name_list = GOLD_STD_DF_HUMAN[COLUMN_NAME_OF_MAVE_GOLD_STANDARD_ID].tolist()
    lopo_dict = LopoBaseline.get_lopo_dict(protein_name_list=protein_name_list)

    def get_PSSM_dict_for_protein(protein_name, lopo_list, df):
        """
        This function calculates the mean of the PSSM values for each amino acid for a protein.
        Input: protein_name, lopo_list, df
        Output: amino_acid_pssm_dict
        """

        amino_acid_pssm_dict = generate_amino_pssm_dict(AMINO_ACIDS_SINGLE_LETTER)
        for lopo_protein in lopo_list:
            snp_dict = get_value_from_dataframe(df=df,
                                     id_column=COLUMN_NAME_OF_MAVE_GOLD_STANDARD_ID,
                                     id_value=lopo_protein,
                                     value_column=COLUMN_NAME_OF_MAVE_GOLD_STANDARD_SNP_DICTIONARY)
            for key, value in snp_dict.items():
                pssm_key = remove_digits_from_key(key)
                amino_acid_pssm_dict[pssm_key].append(value)

        for key in amino_acid_pssm_dict:
            mean_value = calculate_mean(amino_acid_pssm_dict[key])
            amino_acid_pssm_dict[key] = mean_value

        return protein_name, amino_acid_pssm_dict

    # protein_name = "CBS_urn:mavedb:00000005-a"
    # protein_name = "NUDT15_urn:mavedb:00000055-0"
    protein_name = "VKOR_urn:mavedb:00000078-a"
    lopo_list = lopo_dict[protein_name]
    protein_name, amino_acid_pssm_dict_protein_one = get_PSSM_dict_for_protein(protein_name=protein_name,
                                                                    lopo_list=lopo_list,
                                                                   df=GOLD_STD_DF_HUMAN)

    # amino_acid_pssm_dict = generate_amino_pssm_dict(AMINO_ACIDS_SINGLE_LETTER)
    # for protein_name, lopo_list in lopo_dict.items():
    #     for lopo_protein in lopo_list:
    #         snp_dict = get_value_from_dataframe(df=GOLD_STD_DF_HUMAN,
    #                                  id_column=COLUMN_NAME_OF_MAVE_GOLD_STANDARD_ID,
    #                                  id_value=lopo_protein,
    #                                  value_column=COLUMN_NAME_OF_MAVE_GOLD_STANDARD_SNP_DICTIONARY)
    #         for key, value in snp_dict.items():
    #             pssm_key = remove_digits_from_key(key)
    #
    #             amino_acid_pssm_dict[pssm_key].append(value)
    #
    # for key in amino_acid_pssm_dict:
    #     mean_value = calculate_mean(amino_acid_pssm_dict[key])
    #     amino_acid_pssm_dict[key] = mean_value

    print("Debug Pause")

    # import statistics
    #
    #
    # def update_dict_with_mean_values(dictionary):
    #     for key in dictionary:
    #         dictionary[key] = statistics.mean(dictionary[key])
    #     return dictionary
    #
    #
    # amino_acid_pssm_dict = update_dict_with_mean_values(amino_acid_pssm_dict)



    # LopoBaseline.calculate_and_update_mean(df=GOLD_STD_DF_HUMAN,
    #                                        protein_name_column=COLUMN_NAME_OF_MAVE_GOLD_STANDARD_ID,
    #                                        dictionary_column=COLUMN_NAME_OF_MAVE_GOLD_STANDARD_SNP_DICTIONARY,
    #                                        mean_column="mean")
    #

    #
    # lopo_lists = LopoBaseline.leave_one_protein_out_lists(protein_name_list=protein_name_list)
    #

    #
    #
    #
    #
    # for protein_name, lopo_list in lopo_dict.items():
    #     lopo_mean = LopoBaseline.get_lopo_mean(df=GOLD_STD_DF_HUMAN,
    #                                            id_column=COLUMN_NAME_OF_MAVE_GOLD_STANDARD_ID,
    #                                            value_column="mean",
    #                                            lopo_list=lopo_list)
    #
    #     update_value_based_on_protein_name(df=GOLD_STD_DF_HUMAN,
    #                                        column_name='lopo_means',
    #                                        value=lopo_mean,
    #                                        protein_name=protein_name,
    #                                        id_column_name=COLUMN_NAME_OF_MAVE_GOLD_STANDARD_ID)
    #
    # for protein_name in protein_name_list:
    #     protein_snps_dict = get_value_from_dataframe(df=GOLD_STD_DF_HUMAN,
    #                              id_column=COLUMN_NAME_OF_MAVE_GOLD_STANDARD_ID,
    #                              id_value=protein_name,
    #                              value_column=COLUMN_NAME_OF_MAVE_GOLD_STANDARD_SNP_DICTIONARY)
    #
    #
    #     lopo_mean = get_value_from_dataframe(df=GOLD_STD_DF_HUMAN,
    #                              id_column=COLUMN_NAME_OF_MAVE_GOLD_STANDARD_ID,
    #                              id_value=protein_name,
    #                              value_column="lopo_means")
    #
    #
    #
    #     lopo_mean_dict = {key: lopo_mean for key in protein_snps_dict}
    # #
    #     GOLD_STD_DF_HUMAN = update_value_based_on_protein_name(df=GOLD_STD_DF_HUMAN,
    #                                    column_name='lopo_means_dict',
    #                                    value=lopo_mean_dict,
    #                                    protein_name=protein_name,
    #                                    id_column_name=COLUMN_NAME_OF_MAVE_GOLD_STANDARD_ID)
    #
    # pickle_dataframe(dataframe=GOLD_STD_DF_HUMAN,
    #                  file_path=PICKLED_DATAFRAMES_DIRECTORY_PATH,
    #                  file_name=MAVE_DATAFRAME_HUMAN_LOPO_MEAN_PICKLE_FILE_NAME)
    #



