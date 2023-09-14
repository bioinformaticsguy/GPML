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


# Generate file for Mutepred Input
write_file_for_mutepred(mutepred_input_file_path, get_list_of_mute_pred_inputs(gs_dictionary))



# Calculate Spearmen
mutepred_scores_dictionary = get_dictionary_for_mutpred_scores(mutepred_score_file_dir)
score_comparison_list = get_score_comparison_list(tool_score_dict=mutepred_scores_dictionary,
                                                  gs_dictionary=gs_dictionary,
                                                  protein_name=NUDT15_55_0)
(spearman_correlation_score, p_value) = get_spearman_score(scores_for_spearman_comparison_list = score_comparison_list)
print(NUDT15_55_0)
print("Spearman correlation coefficient:", spearman_correlation_score)
print("p-value:", p_value)


# Generate 36 ProFiles
protein_names = gs_dictionary.keys()
generate_one_file_for_each_protein(protein_names, gs_dictionary, folder_path=mutepred_scores, mut_count=None)


# from pathlib import Path

# def get_updated_path(original_path, new_subdirectory, new_extension):
#     new_subdirectory = original_path.parent / new_subdirectory
#     new_path_with_subdir = new_subdirectory / original_path.name
#     new_path_with_extension = new_path_with_subdir.with_suffix(new_extension)

#     return  new_path_with_extension

# # MultiPreprocessing
# import subprocess
# import multiprocessing

# # Function to run the command for a given input and output file
# def run_command(input_file, output_file):
#     command = f"./run_mutpred2.sh -i {input_file} -p 1 -c 1 -b 0 -t 0.1 -f 4 -o {output_file}"
#     subprocess.call(command, shell=True)

# if __name__ == "__main__":
#     input_files = [file for file in multiprocessong_folder_path.rglob('*') if file.is_file()] * 36  # Assuming you want to use the same input file for all cores
#     output_files = [get_updated_path(Path(orignal_path), "output_folder", ".out") for orignal_path in input_files] * 36  # Assuming you want to use the same output file for all cores

#     pool = multiprocessing.Pool(processes=36)  # Number of cores
#     pool.starmap(run_command, zip(input_files, output_files))
#     pool.close()
#     pool.join()



