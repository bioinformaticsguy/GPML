from pathlib import Path
import pandas as pd

from main import MAVE_DB_GOLD_STANDARD_SEQUENCE_ONLY_FILE_PATH
from src.dbNSFP_preprocessor import InputGenerator, dbNSFP_INPUT_FILES_DIRECTORY

MUTEPRED_INPUT_FILE_PATH = Path("Data/mutpred_input_files/mutepred_36.fasta")
META_DATA_INPUT_FILE_PATH = Path("Data/maveGSData/genes_meta_data.csv")

df = pd.read_csv(Path('Data/maveGSData/genes_meta_data.csv'))

uniprot_ids = df['Uniprot ID'].tolist()


# Replace 'desired_value' with the value you're looking for in 'condition_column'
desired_value = 'Human'

# Create a list of elements from 'your_column_name' based on the condition
uniprot_id_list_of_human_genes = df[df['Organism'] == desired_value]['Uniprot ID'].tolist()

all_protein_names = list(InputGenerator.get_protein_names_and_snps_from_mutepred_input_file(MUTEPRED_INPUT_FILE_PATH).keys())


dict = InputGenerator.get_protein_names_and_snps_from_mutepred_input_file(MUTEPRED_INPUT_FILE_PATH)
[print(i, len(dict[i])) for i in dict.keys()]

# InputGenerator.write_input_files(mutepred_input_file_path = MUTEPRED_INPUT_FILE_PATH,
#                                  genes_meta_data_file_path = META_DATA_INPUT_FILE_PATH,
#                                  dbNSFP_input_files_directory = dbNSFP_INPUT_FILES_DIRECTORY)



