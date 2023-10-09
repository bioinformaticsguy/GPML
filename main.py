from src.utils import get_dictonary_of_scores, get_list_of_mute_pred_inputs, write_file_for_mutepred, \
    get_dictionary_for_mutpred_scores, get_score_comparison_list, get_spearman_score, get_mute_pred_input, \
    generate_one_file_for_each_protein

from pprint import pprint as pretty_print
from pathlib import Path


# Path Variables
mave_db_file_path = Path("Data/maveGSData/mave_db_gold_standard.fasta")
mutepred_input_file_path = Path("Data/mutpred_input_files/mutepred_36.fasta")
mutepred_score_file_dir = Path("Data/mutepred_scores/tNUDT15.out")
mutepred_scores = Path("Data/mutepred_scores/f_/")
multiprocessong_folder_path = Path("Data/test_multi_processor")

# Protein Names
NUDT15_55_0 = "NUDT15_urn:mavedb:00000055-0"

# Data Type Variables
gs_dictionary = get_dictonary_of_scores(mave_db_file_path)
# pretty_print(gs_dictionary.keys())


# # Generate file for Mutepred Input
# write_file_for_mutepred(mutepred_input_file_path, get_list_of_mute_pred_inputs(gs_dictionary))



# # Calculate Spearmen
# mutepred_scores_dictionary = get_dictionary_for_mutpred_scores(mutepred_score_file_dir)
# score_comparison_list = get_score_comparison_list(tool_score_dict=mutepred_scores_dictionary,
#                                                   gs_dictionary=gs_dictionary,
#                                                   protein_name=NUDT15_55_0)
# (spearman_correlation_score, p_value) = get_spearman_score(scores_for_spearman_comparison_list = score_comparison_list)
# print(NUDT15_55_0)
# print("Spearman correlation coefficient:", spearman_correlation_score)
# print("p-value:", p_value)


# # Generate 36 ProFiles
# protein_names = gs_dictionary.keys()
# generate_one_file_for_each_protein(protein_names, gs_dictionary, folder_path=mutepred_scores, mut_count=None)



unprocessed_msa_folder_path = Path("Data/generated_msa/")
file_name = "NUDT15_uniref100_05.a2m"


output_file_name = "processed.atm"


# print(unprocessed_msa_path)

def process_msa_file_headers(unprocessed_msa_folder_path, input_file_name):
    input_file_path = unprocessed_msa_folder_path / input_file_name
    processed_lines = []
    with open(input_file_path) as file:
        for line in file:
            if line.startswith(">"):
                processed_lines.append(line.split()[0]+"\n")

            else:
                processed_lines.append(line)
    
    return processed_lines


def write_file_from_list_of_lines(list_of_lines, output_folder, output_file_name):
    output_file_path = output_folder / output_file_name
    with open(output_file_path, "w") as file:
        for line in list_of_lines:
            file.write(line)

    file.close()

list_of_lines = process_msa_file_headers(unprocessed_msa_folder_path, file_name)


write_file_from_list_of_lines(list_of_lines, 
                              output_folder=unprocessed_msa_folder_path, 
                              output_file_name=output_file_name)