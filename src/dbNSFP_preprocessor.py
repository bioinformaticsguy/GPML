from pathlib import Path
import pandas as pd

dbNSFP_INPUT_FILES_DIRECTORY = Path("Data/dbNSFP_input_files")

class InputGenerator:
    @staticmethod
    def get_protein_names_and_snps_from_mutepred_input_file(file_path):
        dictionary_of_protein_names_and_single_char_snp_lists = {}
        with open(file_path, "r") as f:
            for line in f:
                if line[0] == ">":
                    all_protein_data_list = line.split(" ")
                    all_protein_data_list[-1] = all_protein_data_list[-1].strip()
                    dictionary_of_protein_names_and_single_char_snp_lists[all_protein_data_list[0][1:]] = all_protein_data_list[1:]

        return dictionary_of_protein_names_and_single_char_snp_lists

    @staticmethod
    def get_dictionary_with_uniprot_ids(mutepred_input_file_path, genes_meta_data_file_path):
        dictionary_with_human_protein_data_only = {}
        dictionary_of_protein_names_and_single_char_snp_lists = InputGenerator.get_protein_names_and_snps_from_mutepred_input_file(mutepred_input_file_path)
        genes_meta_data = pd.read_csv(Path(genes_meta_data_file_path))
        protein_name_list_of_human_genes = genes_meta_data[genes_meta_data['Organism'] == "Human"]['Protein Name'].tolist()
        for protein_name in protein_name_list_of_human_genes:
            uniprot_id = genes_meta_data.loc[genes_meta_data['Protein Name'] == protein_name, 'Uniprot ID'].values[0]
            protein_single_letter_snp_list = dictionary_of_protein_names_and_single_char_snp_lists[protein_name]
            dictionary_with_human_protein_data_only[protein_name] = (uniprot_id, protein_single_letter_snp_list)
            # print(protein_name, uniprot_id)

        return dictionary_with_human_protein_data_only

    @staticmethod
    def write_input_files(mutepred_input_file_path, genes_meta_data_file_path, dbNSFP_input_files_directory):
        protein_data_dict = InputGenerator.get_dictionary_with_uniprot_ids(mutepred_input_file_path, genes_meta_data_file_path)
        for protein_name in protein_data_dict:
            uniprot_id, snps = protein_data_dict[protein_name]
            file_nmae = protein_name + ".in"
            with open(dbNSFP_input_files_directory/file_nmae, 'w') as f:

                for snp in snps:
                    f.write("UNIPROT"+":"+uniprot_id+":"+snp+"\n")

                    # print("UNIPROT"+":"+uniprot_id+":"+snp)
                    # break
                # UNIPROT: Q9BQB6:H163Y


