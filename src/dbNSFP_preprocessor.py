import csv
from pathlib import Path
import pandas as pd
import os, re
from src.utils import get_count_of_lines_except_header
from pprint import pprint

dbNSFP_INPUT_FILES_DIRECTORY = Path("Data/dbNSFP_input_dir")
dbNSFP_OUTPUT_FILES_DIRECTORY = Path("Data/dbNSFP_output_dir")
SNP_COLUMN_NAME = "HGVSp_ANNOVAR"
MUTEPRED_SCORE_COLUMN_NAME = "MutPred_score"



class InputGenerator:
    @staticmethod
    def get_protein_names_and_snps_from_mutepred_input_file(file_path):
        dictionary_of_protein_names_and_single_char_snp_lists = {}
        with open(file_path, "r") as f:
            for line in f:
                if line[0] == ">":
                    all_protein_data_list = line.split(" ")
                    all_protein_data_list[-1] = all_protein_data_list[-1].strip()
                    dictionary_of_protein_names_and_single_char_snp_lists[
                        all_protein_data_list[0][1:]] = all_protein_data_list[1:]

        return dictionary_of_protein_names_and_single_char_snp_lists

    @staticmethod
    def get_dictionary_with_uniprot_ids(mutepred_input_file_path, genes_meta_data_file_path):
        dictionary_with_human_protein_data_only = {}
        dictionary_of_protein_names_and_single_char_snp_lists = InputGenerator.get_protein_names_and_snps_from_mutepred_input_file(
            mutepred_input_file_path)
        genes_meta_data = pd.read_csv(Path(genes_meta_data_file_path))
        protein_name_list_of_human_genes = genes_meta_data[genes_meta_data['Organism'] == "Human"][
            'Protein Name'].tolist()
        for protein_name in protein_name_list_of_human_genes:
            uniprot_id = genes_meta_data.loc[genes_meta_data['Protein Name'] == protein_name, 'Uniprot ID'].values[0]
            protein_single_letter_snp_list = dictionary_of_protein_names_and_single_char_snp_lists[protein_name]
            dictionary_with_human_protein_data_only[protein_name] = (uniprot_id, protein_single_letter_snp_list)
            # print(protein_name, uniprot_id)

        return dictionary_with_human_protein_data_only

    @staticmethod
    def write_input_files(mutepred_input_file_path, genes_meta_data_file_path, dbNSFP_input_files_directory):
        protein_data_dict = InputGenerator.get_dictionary_with_uniprot_ids(mutepred_input_file_path,
                                                                           genes_meta_data_file_path)
        for protein_name in protein_data_dict:
            uniprot_id, snps = protein_data_dict[protein_name]
            file_nmae = protein_name + ".in"
            with open(dbNSFP_input_files_directory / file_nmae, 'w') as f:

                for snp in snps:
                    f.write("UNIPROT" + ":" + uniprot_id + ":" + snp + "\n")



class OutputAnalyzer:
    @staticmethod
    def get_count_dict_for_count_of_snps(output_files_dir=dbNSFP_OUTPUT_FILES_DIRECTORY):
        dict_for_count_of_snps = {}
        snp_file_paths = [Path(file) for file in output_files_dir.iterdir() if file.is_file() and file.suffix == ".csv"]
        err_file_paths = [Path(file) for file in output_files_dir.iterdir() if file.is_file() and file.suffix == ".err"]

        for snp_file_path in snp_file_paths:
            protein_name = os.path.splitext(snp_file_path.name)[0]
            snp_count = get_count_of_lines_except_header(snp_file_path)
            dict_for_count_of_snps[protein_name] = [snp_file_path]
            dict_for_count_of_snps[protein_name].append(snp_count)

        for err_file_path in err_file_paths:
            protein_name = Path(os.path.splitext(err_file_path.name)[0]).stem
            dict_for_count_of_snps[protein_name].append(err_file_path)
            not_found_snp_count = get_count_of_lines_except_header(err_file_path)
            dict_for_count_of_snps[protein_name].append(not_found_snp_count)

        return dict_for_count_of_snps

    @staticmethod
    def print_protein_names_and_snp_found_not_found(dict_for_count_of_snps):
        list_of_print_elements = sorted([[key,
                                          dict_for_count_of_snps[key][1],
                                          dict_for_count_of_snps[key][3]] for key in dict_for_count_of_snps.keys()])
        [print(i[0], '\t', i[1], '\t', i[2]) for i in [list_elem for list_elem in list_of_print_elements]]

    @staticmethod
    def get_protein_name_and_list_of_tuples(protein_file_path, snp_column_name, tool_score_column_name):

        def _extract_snp_and_score(snp_str, score_str):
            # Extracting SNPs and their score
            snps = re.findall(r'p\.([A-Z0-9]+)', snp_str)
            # print(score_str)

            # Checking if all SNPs are the same
            # are_snps_same = snps and all(snp == snps[0] for snp in snps)
            # print(are_snps_same)

            # Checking for a numeric score before converting to float
            # score = float(score_str.split(';')[0]) if score_str.split(';')[0].replace('.', '').isdigit() else None
            scores = str(score_str).split(';')[:len(snps)]

            # Return the extracted SNP and its score or None if not valid
            return list(set(list(zip(snps, scores))))

            # return (snps[0], score) if are_snps_same else (snps, score)

        df = pd.read_csv(protein_file_path, sep='\t')
        # headers = df.columns.tolist()
        # # print(headers)

        values_list = []

        # # values_list = [tuple_snp_score for tuple_snp_score in (_extract_snp_and_score(snp_str=row[snp_column_name],
        #                                       score_str=row[tool_score_column_name])) for _, row in df.iterrows()]

        for _, row in df.iterrows():
            for tuple_snp_score in _extract_snp_and_score(snp_str=row[snp_column_name],
                                                          score_str=row[tool_score_column_name]):
                values_list.append(tuple_snp_score)

        filtered_tuples = [(snp, float(score)) for snp, score in values_list if score != '.']


        return dict(filtered_tuples)

    @staticmethod
    def get_dictionary_for_all_proteins_in_dictionary(dbNSFP_output_files_directory = dbNSFP_OUTPUT_FILES_DIRECTORY,
                                                      snp_column_name = SNP_COLUMN_NAME,
                                                      tool_score_column_name = MUTEPRED_SCORE_COLUMN_NAME):

        dictionary_of_tool_scores = {}
        protein_names = []
        snp_file_paths = [Path(file) for file in dbNSFP_output_files_directory.iterdir() if file.is_file() and file.suffix == ".csv"]
        for snp_file_path in snp_file_paths:
            protein_name = os.path.splitext(snp_file_path.name)[0].replace("_output", "")\
                                                                  .replace("_urn_", "_urn:")\
                                                                  .replace("mavedb_", "mavedb:")

            protein_names.append(protein_name)
            scores_list = OutputAnalyzer.get_protein_name_and_list_of_tuples(snp_file_path,
                                                                             snp_column_name,
                                                                             tool_score_column_name)
            dictionary_of_tool_scores[protein_name] = scores_list


        # [print(i) for i in sorted(protein_names)]


        return dictionary_of_tool_scores

