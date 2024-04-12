from src.baseline_calculation import LopoBaseline
from src.constants import PICKLED_DATAFRAMES_DIRECTORY_PATH, MAVE_DATAFRAME_PICKLE_FILE_NAME, SPECIE_NAME_HUMAN, \
    COLUMN_NAME_OF_MAVE_GOLD_STANDARD_SNP_DICTIONARY, COLUMN_NAME_OF_MAVE_GOLD_STANDARD_SPECIES, \
    COLUMN_NAME_OF_MAVE_GOLD_STANDARD_ID, MAVE_DATAFRAME_HUMAN_LOPO_MEAN_PICKLE_FILE_NAME, AMINO_ACIDS_SINGLE_LETTER
from src.utils import load_dataframe, filter_dataframe_by_species, get_value_from_dataframe, \
    update_value_based_on_protein_name, pickle_dataframe, generate_amino_pssm_dict

if __name__ == '__main__':
    GOLD_STD_DF = load_dataframe(file_path=PICKLED_DATAFRAMES_DIRECTORY_PATH,
                                 file_name=MAVE_DATAFRAME_PICKLE_FILE_NAME)

    GOLD_STD_DF_HUMAN= filter_dataframe_by_species(df=GOLD_STD_DF,
                                target_species=SPECIE_NAME_HUMAN,
                                species_column=COLUMN_NAME_OF_MAVE_GOLD_STANDARD_SPECIES)

    protein_name_list = GOLD_STD_DF_HUMAN[COLUMN_NAME_OF_MAVE_GOLD_STANDARD_ID].tolist()
    lopo_dict = LopoBaseline.get_lopo_dict(protein_name_list=protein_name_list)

    protein_one = protein_name_list[0]

    amino_acid_pssm_dict = generate_amino_pssm_dict(AMINO_ACIDS_SINGLE_LETTER)




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


    """
    Sudo Code for PSSM baseline calculation
    1. Leave 1st protein and take remaining 12 proteins.
    2. Create an empty PSSM matrix.
    3. Regardless of the position of the SNP's in the 12 proteins, calculate mean for each Unique amino acid substitution.
    4. Update the mean values in the PSSM matrix.
    5. Impute values for the 1st protein using the PSSM matrix for each SNP regardless of the position.
    6. Repeat the process for all the proteins.
    """

    print("Debug Pause")