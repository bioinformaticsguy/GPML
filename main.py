import csv
import pandas as pd
from pprint import pprint as pretty_print
from pathlib import Path
from src.utils import get_dictonary_of_scores, get_list_of_mute_pred_inputs, write_file_for_mutepred, \
    get_dictionary_for_mutpred_scores, get_score_comparison_list, get_spearman_score, get_mute_pred_input, \
    generate_one_file_for_each_protein, get_mutepred_dictionary_of_scores, read_fasta_iteration


# Path Variables
mave_db_gold_standard_fasta_file_path = Path("Data/maveGSData/mave_db_gold_standard.fasta")
MAVE_DB_GOLD_STANDARD_SEQUENCE_ONLY_FILE_PATH = Path("Data/maveGSData/mave_db_gold_standard_only_sequences.fasta")
nine_human_seqs_uniprot_file_path = Path("Data/maveGSData/nine_human_seqs_uniprot.fasta")
mutepred_input_file_path = Path("Data/mutpred_input_files/mutepred_36.fasta")
mutepred_score_file_dir = Path("Data/mutepred_scores/tNUDT15.out")
mutepred_scores = Path("Data/mutepred_scores/f_/")
multiprocessong_folder_path = Path("Data/test_multi_processor")
mutepred_output_file_path = Path("Data/mutepred_sore_comp/test.out")
csv_file_path = Path("Data/mutepred_sore_comp/spearman_scores.csv")


# Protein Names
NUDT15_55_0 = "NUDT15_urn:mavedb:00000055-0"

# Data Type Variables
gs_dictionary = get_dictonary_of_scores(mave_db_gold_standard_fasta_file_path)


# Get mutepred_dictionary_of_scores
mutepred_dictionary_of_scores = get_mutepred_dictionary_of_scores(mutepred_output_file_path)













# Get Spearman Scores

## generate the spearman file.
# spearman_scores_list = get_spearman_scores_for_all_mutepred_proetins(mutepred_dictionary_of_scores, gs_dictionary)
# write_csv_file_for_spearman_scores(csv_file_path, spearman_scores_list)
# print(mutepred_dictionary_of_scores["VIM-2_with_p.Met1_Phe2insGly_urn:mavedb:00000073-c"].keys())

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

# unprocessed_msa_folder_path = Path("Data/generated_msa/")
# # print(len(gs_dictionary.keys()))
# # print(gs_dictionary['VIM-2_with_p.Met1_Phe2insGly_urn:mavedb:00000073-c'])