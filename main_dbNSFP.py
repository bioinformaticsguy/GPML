from pathlib import Path
from main import MAVE_DB_GOLD_STANDARD_SEQUENCE_ONLY_FILE_PATH, gs_dictionary
from src.dbNSFP_preprocessor import InputGenerator, dbNSFP_INPUT_FILES_DIRECTORY, OutputAnalyzer
from src.utils import get_spearman_scores_for_all_tool_proetins
import csv

MUTEPRED_INPUT_FILE_PATH = Path("Data/mutpred_input_files/mutepred_36.fasta")
META_DATA_INPUT_FILE_PATH = Path("Data/mave_gs_data/genes_meta_data.csv")

# df = pd.read_csv(Path('Data/mave_gs_data/genes_meta_data.csv'))
#
# uniprot_ids = df['Uniprot ID'].tolist()
#
#
# # Replace 'desired_value' with the value you're looking for in 'condition_column'
# desired_value = 'Human'
#
# # Create a list of elements from 'your_column_name' based on the condition
# uniprot_id_list_of_human_genes = df[df['Organism'] == desired_value]['Uniprot ID'].tolist()

# all_protein_names = list(InputGenerator.get_protein_names_and_snps_from_mutepred_input_file(MUTEPRED_INPUT_FILE_PATH).keys())


# InputGenerator.write_input_files(mutepred_input_file_path = MUTEPRED_INPUT_FILE_PATH,
#                                  genes_meta_data_file_path = META_DATA_INPUT_FILE_PATH,
#                                  dbNSFP_input_files_directory = dbNSFP_INPUT_FILES_DIRECTORY)




# OutputAnalyzer.print_protein_names_and_snp_found_not_found(OutputAnalyzer.get_count_dict_for_count_of_snps())


# csv_file_path = "/Users/ali/Documents/GPML/Data/dbNSFP_output_dir/VKOR_urn_mavedb_00000078-a_output.csv"
#
# SNP_COLUMN_NAME = "HGVSp_ANNOVAR"
# MUTEPRED_SCORE_COLUMN_NAME = "MutPred_score"


# tuples = OutputAnalyzer.get_protein_name_and_list_of_tuples(csv_file_path, SNP_COLUMN_NAME, tool_score_column_name)
# [print(i) for i in tuples]
#



# b_dict = OutputAnalyzer.get_dictionary_for_all_proteins_in_dictionary()
#
#
# spearmen_scores = get_spearman_scores_for_all_tool_proetins(b_dict ,gs_dictionary)
#

# [print(i) for i in sorted(spearmen_scores)]


tool_score_header_names_list = ["PROVEAN_score",
                                "REVEL_score",
                                "MutPred_score",
                                "Polyphen2_HDIV_score",
                                "EVE_score",
                                "AlphaMissense_score"]

def dictionary_of_spearman_corelations_of_different_tools(tool_score_header_names_list, gs_dictionary):
    dictionary_of_correlations = {}
    for tool_name in tool_score_header_names_list:
        print(tool_name)
        b_dict = OutputAnalyzer.get_dictionary_for_all_proteins_in_dictionary(tool_score_column_name=tool_name)
        spearmen_scores = get_spearman_scores_for_all_tool_proetins(b_dict, gs_dictionary)
        dictionary_of_correlations[tool_name] = spearmen_scores

    return dictionary_of_correlations


corr_dict = dictionary_of_spearman_corelations_of_different_tools(tool_score_header_names_list, gs_dictionary)
# print(corr_dict)

def get_dict_for_csv_data(corr_dict):
    data_for_csv_dict = {}
    gene_names = [entry[0] for tool_scores in corr_dict.values() for entry in tool_scores]
    # print(gene_names)
    for gene_name in gene_names:
        data_for_csv_dict[gene_name] = {}
        for tool_score_name in tool_score_header_names_list:
            a = [tup[1] for tup in corr_dict[tool_score_name] if tup[0] == gene_name]
            # print(a[0])
            data_for_csv_dict[gene_name][tool_score_name] = abs(a[0])

    return data_for_csv_dict, gene_names

data_for_csv_writing, gene_names = get_dict_for_csv_data(corr_dict)



csv_file_path = Path("Data/dbNSFP_analysis_output_dir/spearman_comparisons.csv")
rows_to_write_in_csv = []
with open(csv_file_path, mode='w', newline='') as csv_file:
    csv_writer = csv.writer(csv_file)
#
    # Write the header
    tool_names = [tool for tool in corr_dict.keys()]
    csv_writer.writerow(["Gene Name"] + [f"{tool_name}_cor" for tool_name in tool_names])
    for gene_name in gene_names:
        row_to_write = [gene_name]
        for tool_name in tool_names:
            row_to_write.append(data_for_csv_writing[gene_name][tool_name])
        # print(row_to_write)
        rows_to_write_in_csv.append(row_to_write)


    sadf = my_set_of_tuples = set(tuple(inner_list) for inner_list in rows_to_write_in_csv)
    for i in sadf:
        csv_writer.writerow(list(i))

#
#     gene_names = [entry[0] for tool_scores in corr_dict.values() for entry in tool_scores]
#
#     # Print the extracted gene names
#     # print(gene_names)
#
#     for gene_name in gene_names:
#         # Write the data
#         scores_list_row = []
#         for tool_name in tool_names:
#             # print(tool_name)
#             result = [tup[1] for tup in corr_dict[tool_name] if tup[0] == gene_name]
#             scores_list_row.append(result)
#         break
#             # row = [gene_name] + [entry[1] for entry in scores]
#             # csv_writer.writerow(row)
#
# # print(f"CSV file written successfully at: {csv_file_path}")





