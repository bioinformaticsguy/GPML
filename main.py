from src.utils import get_dictonary_of_scores, get_list_of_mute_pred_inputs, write_file_for_mutepred, \
    get_dictionary_for_mutpred_scores, get_score_comparison_list, get_spearman_score, get_mute_pred_input

from pprint import pprint as pretty_print
from pathlib import Path


# Path Variables
mave_db_file_path = Path("Data/maveGSData/mave_db_gold_standard.fasta")
mutepred_input_file_path = Path("Data/mutpred_input_files/mutepred_36.fasta")
mutepred_score_file_dir = Path("Data/mutepred_scores/tNUDT15.out")
mutepred_scores = Path("Data/mutepred_scores/f_/")

# Protein Names
NUDT15_55_0 = "NUDT15_urn:mavedb:00000055-0"

# Data Type Variables
gs_dictionary = get_dictonary_of_scores(mave_db_file_path)
# pretty_print(gs_dictionary.keys())


# Generate file for Mutepred Input
# write_file_for_mutepred(mutepred_input_file_path, get_list_of_mute_pred_inputs(gs_dictionary))



# Calculate Spearmen
# mutepred_scores_dictionary = get_dictionary_for_mutpred_scores(mutepred_score_file_dir)
# score_comparison_list = get_score_comparison_list(tool_score_dict=mutepred_scores_dictionary,
#                                                   gs_dictionary=gs_dictionary,
#                                                   protein_name=NUDT15_55_0)
# (spearman_correlation_score, p_value) = get_spearman_score(scores_for_spearman_comparison_list = score_comparison_list)
# print(NUDT15_55_0)
# print("Spearman correlation coefficient:", spearman_correlation_score)
# print("p-value:", p_value)


# Generate 36 ProFiles
protein_names = gs_dictionary.keys()

def generate_one_file_for_each_protein(protein_names, gs_dictionary, folder_path):
    for protein_name in protein_names:
        protein_data = gs_dictionary[protein_name]
        mutepred_input = [get_mute_pred_input(protein_name, protein_data)]

        file_name = str(folder_path) + protein_name + ".fasta"

        write_file_for_mutepred(file_name, mutepred_input)


generate_one_file_for_each_protein(protein_names, gs_dictionary, folder_path=mutepred_scores)