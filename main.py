from src.utils import get_dictonary_of_scores, get_list_of_mute_pred_inputs, write_file_for_mutepred
from pprint import pprint as pretty_print
from pathlib import Path




mave_db_file_path = Path("/Users/ali/Documents/GPML/Data/maveGSData/mave_db_gold_standard.fasta")
gs_dictionary = get_dictonary_of_scores(mave_db_file_path)
mutepred_input_file_path = Path("/Data/mutpred_input_files/mutepred_36.fasta")

write_file_for_mutepred(mutepred_input_file_path, get_list_of_mute_pred_inputs(gs_dictionary))




# print(get_list_of_mute_pred_inputs(gs_dictionary)[0:1])