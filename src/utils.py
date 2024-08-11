import itertools
import math
import pickle
import numpy as np
import scipy.stats as stats
import pandas as pd
from scipy.stats import spearmanr

from src.constants import COLUMN_NAME_OF_MAVE_GOLD_STANDARD_SAVS, \
    COLUMN_NAME_OF_MAVE_GOLD_STANDARD_SPECIES, SPEAR_COR_SUFFIX, USED_SAV_PERCENTAGE_SUFFIX, \
    TOOL_SCORE_COLUMN_SUFFIX, TRAINING_SAVS_COLUMN_SIFFIX, EXCLUDE_TRAINING_SAV_SUFFIX, \
    COLUMN_NAME_OF_MAVE_GOLD_STANDARD_ID, AMINO_ACIDS_SINGLE_LETTER, SAVS_STRING


def get_dictonary_of_scores_maveDB(file_path):
    """
    This function takes the file which has the mutations and the sequences, 
    then it stores the data related to each protein in one list and then, 
    adds it to the dictonary with keys as the names of the proteins.
    """
    all_data_dict = {}
    point_mutations = []
    descrip = None
    with open(file_path) as file:
        for line in file:
            if line[0] == ">":
                if descrip:
                    all_data_dict[descrip] = point_mutations
                descrip = line[1:-1]
                point_mutations = []
            else:
                point_mutations.append(line[:-1])
        all_data_dict[descrip] = point_mutations

    return all_data_dict


def convert_three_letter_to_one_letter(three_letter_code):
    """This function takes the three letter amino acids 
    in the form of a string and then returns the one letter amino acid."""
    amino_acids = {
        'Ala': 'A',
        'Arg': 'R',
        'Asn': 'N',
        'Asp': 'D',
        'Cys': 'C',
        'Gln': 'Q',
        'Glu': 'E',
        'Gly': 'G',
        'His': 'H',
        'Ile': 'I',
        'Leu': 'L',
        'Lys': 'K',
        'Met': 'M',
        'Phe': 'F',
        'Pro': 'P',
        'Ser': 'S',
        'Thr': 'T',
        'Trp': 'W',
        'Tyr': 'Y',
        'Val': 'V'
    }

    one_letter_code = amino_acids.get(three_letter_code)
    return one_letter_code


def get_single_letter_point_mutation(three_letter_point_mutation):
    """this function takes the point mutation which is in the format of 
    three letter amino acid, position and then the three letter amino acid,
    where the first amino is the orignal amino ac
    is the position where the mutation occured, last the three letter amino 
    acid is the amino acid which turned out to be in the case of mutation."""
    first = convert_three_letter_to_one_letter(three_letter_point_mutation[:3])
    second = convert_three_letter_to_one_letter(three_letter_point_mutation[-3:])
    pos = three_letter_point_mutation[3:-3]
    return first + pos + second


def get_mutations_list(mutations_data):
    point_mutation_list = []
    for line in mutations_data:
        if line[0] == "<":
            point_mutation_list.append(get_single_letter_point_mutation(line.split()[0][3:]))

    return point_mutation_list


def write_from_list(list, file_name):
    with open(file_name, 'w') as file:
        file.write("mutations\n")
        for i in list:
            file.write(i + '\n')
    file.close()


def get_dictionary_for_eve_scores(eve_score_file_path):
    scores_dictionary = {}
    with open(eve_score_file_path) as file:
        for line in file:
            if line[:12] == "protein_name":
                pass
            else:
                # print(line.split(","))
                mutation_name = line.split(",")[1]
                evol_index = line.split(",")[2]
                eve_score = line.split(",")[3]

                # print(mutation_name, evol_index, eve_score)

                scores_dictionary[mutation_name] = [evol_index, eve_score]

    return scores_dictionary


def get_score_comparison_list(tool_score_dict, gs_dictionary, protein_name):
    # print(protein_name)
    scores_list = []
    count = 0
    protein_data = gs_dictionary[protein_name]

    # pretty_print(protein_data)

    for line in protein_data:
        # print(line[0])
        if line[0] == "<":

            mutation = get_single_letter_point_mutation(line.split()[0][3:])
            scaled_effect = float(line.split(":")[1])
            # pretty_print(mutation)

            try:
                tool_score_value = tool_score_dict[mutation]
                scores_list.append([mutation, scaled_effect, tool_score_value])
                # print(eve_score_value)
            except KeyError:
                count += 1
                # print("Key not found in dictionary")

    return scores_list


def get_spearman_score(scores_for_spearman_comparison_list):
    """
    Calculate the Spearman correlation coefficient and p-value for a list of scores.

    Parameters:
    - scores_for_spearman_comparison_list (List[Tuple]): A list of tuples, where each tuple
      contains the necessary scores for Spearman comparison.

    Returns:
    - Tuple[float, float]: A tuple containing the Spearman correlation coefficient and its p-value.

    Example:
    ```
    scores_list = [(protein1, scaled_effect1, eve_score_value1), (protein2, scaled_effect2, eve_score_value2)]
    spearman_correlation, p_value = get_spearman_score(scores_list)
    ```
    """
    scaled_effects = []
    eve_score_values = []
    for score_list in scores_for_spearman_comparison_list:
        scaled_effects.append(score_list[1])
        eve_score_values.append(score_list[2])
        # print(score_list[1])

    # Calculate Spearman correlation coefficient and p-value
    spearman_correlation_score, p_value = stats.spearmanr(scaled_effects, eve_score_values)

    return (spearman_correlation_score, p_value)


def remove_stop_codon(sequence):
    if sequence[-1] == "*":
        return sequence[:-1]
    else:
        return sequence


def get_mute_pred_input(pro_name, data, mut_count=None):
    """Input
            pro_name: string --> VIM-2_with_p.Met1_Phe2insGly_urn:mavedb:00000073-c
            data:  list of strings where last string is the sequence --> see the tail 
                            ['<p.Met1Ala #scaled_effect:0.29943354034999425',
                            'MGFKPVGEVRLYQI']
            mut_count: intege to limit the number of mutations to add infront of name
                        of the protein.

        Output:
            Reutrns a touple with two elements each element is a string
            the first element contains the name of protein and all the point mutations
            second element contains the sequence of the protein itself.
            """
    seq = remove_stop_codon(data[-1])
    mut_lines = data[:-1]
    # print(len(mut_lines))

    name_and_mutations = ">" + pro_name

    if mut_count == None:
        for mut_line in mut_lines:
            mutation = get_single_letter_point_mutation(mut_line.split(".")[1].split(" ")[0])

            name_and_mutations = name_and_mutations + " " + mutation
    elif isinstance(mut_count, int):
        for i in range(mut_count):
            mutation = get_single_letter_point_mutation(mut_lines[i].split(".")[1].split(" ")[0])
            name_and_mutations = name_and_mutations + " " + mutation
            # print(mutation)

    return (name_and_mutations + "\n", seq + "\n")


def write_file_for_mutepred(name_with_ext, data_tuples):
    with open(name_with_ext, 'w') as file:
        for data_tuple in data_tuples:
            for item in data_tuple:
                # print(item)
                file.write(item)


def read_fasta_iteration(filename):
    """Given a fasta file returns a dictionary
    in which the key is the description of the fasta sequence
    and the value is the sequence itself."""

    descrip = None
    seq_dict = {}
    with open(filename) as file:
        for line in file:
            if line[0] == ">":
                if descrip:
                    seq_dict[descrip] = seq
                    # sequences.append((descrip, seq))
                descrip = line
                seq = ''
            else:
                seq = seq + line[:-1]
        seq_dict[descrip] = seq

    return seq_dict


def get_list_of_mute_pred_inputs(gs_dictionary, mut_count=None):
    """
    Input:
        gs_dictionary: this dictionary contains the names of the proteins
                        as the keys and the values are the lists containg
                        all the point mutations and the last entry in the
                        list is the sequence of that protein.

    Return:
        this function returns the list which contains the touples, each touple
        structure is described in the function mute_pred_input()
    """
    list_of_mute_pred_inputs = []
    protein_names = gs_dictionary.keys()
    for protein_name in protein_names:
        data = gs_dictionary[protein_name]
        list_of_mute_pred_inputs.append(get_mute_pred_input(protein_name, data))

    return list_of_mute_pred_inputs


def generate_one_file_for_each_protein(protein_names, gs_dictionary, folder_path, mut_count=None):
    """This functioin takes the names of the proteins and it takes the
    gold standard dictionary, also it takes the folder path. It writes all the
    files with the specific mutations in that folder, there is one small problem
    that needs to be fixed with the naming and the folder directory. I will update
    it from my mac book"""
    for protein_name in protein_names:
        protein_data = gs_dictionary[protein_name]
        mutepred_input = [get_mute_pred_input(protein_name, protein_data, mut_count=mut_count)]

        file_name = str(folder_path) + protein_name + ".fasta"

        write_file_for_mutepred(file_name, mutepred_input)


def process_msa_file_headers(unprocessed_msa_folder_path, input_file_name):
    input_file_path = unprocessed_msa_folder_path / input_file_name
    processed_lines = []
    with open(input_file_path) as file:
        for line in file:
            if line.startswith(">"):
                processed_lines.append(line.split()[0] + "\n")

            else:
                processed_lines.append(line)

    return processed_lines


def write_file_from_list_of_lines(list_of_lines, output_folder, output_file_name):
    output_file_path = output_folder / output_file_name
    with open(output_file_path, "w") as file:
        for line in list_of_lines:
            file.write(line)

    file.close()


def get_dictionary_for_mutpred_scores(score_file_path):
    """
    Input:
            This function takes the path of the output file
            with only scores using the 4th type of output
            from the mutepred Pipline.

    Return:
            It returns the dictionary, so that each key if the
            mutation name and each value is the mutepred score.
    """

    scores_dictionary = {}
    with open(score_file_path) as file:
        for line in file:
            if line[:2] == "ID":
                pass
            else:
                mutation_name = line.split(",")[1]
                mutepred_score = line.split(",")[2]
                scores_dictionary[mutation_name] = [mutepred_score]
    return scores_dictionary


def get_mutepred_dictionary_of_scores(file_path):
    """
    This function takes the output file with multiple proteins from mutepred
    and then it generates a dictionary with having keys and the names of the
    proteins and value is the list with contains the touples, where
    each touple has two elements first one for the mutation and the second
    one for the mutepred score.

    {protein_name: [('G2Y', '0.243'), ('G2W', '0.273'), ('G2V', '0.207')]}
    """
    protein_set = set()
    scores_dictionary = {}
    with open(file_path) as file:
        for line in file:
            if line[:2] == "ID":
                pass
            else:
                protein_name = line.split(",")[0]
                mutation_name = line.split(",")[1]
                mutepred_score = line.split(",")[2]
                if protein_name in protein_set:
                    scores_dictionary[protein_name][mutation_name] = float(mutepred_score)
                else:
                    scores_dictionary[protein_name] = {}
                    scores_dictionary[protein_name][mutation_name] = float(mutepred_score)
                    protein_set.add(protein_name)

    return scores_dictionary


def write_csv_file_for_spearman_scores(csv_file_path, spearman_scores_list):
    """
    This function takes the csv file path and the list of the spearman scores and then
    it generates the csv file.
    """
    with open(csv_file_path, 'w') as csv_file:
        # Writing header
        csv_file.write(','.join(map(str, spearman_scores_list[0])) + '\n')

        # Writing data rows
        for row in spearman_scores_list[1:]:
            csv_file.write(','.join(map(str, row)) + '\n')
        # print(spearman_score)


def get_spearman_scores_for_all_tool_proetins(tool_dictionary_of_scores, gs_dictionary):
    """This function takes the mutepred dictionary of scores and the gold
    standard dictionary and then it generates the spearman scores list
    for all the proteins using the data from the mutepred dictionary of scores."""
    scores_list = []
    protein_names = tool_dictionary_of_scores.keys()
    for protein_name in protein_names:
        scores_for_spearman_comparison_list = get_score_comparison_list(tool_dictionary_of_scores[protein_name],
                                                                        gs_dictionary, protein_name)
        # print(scores_for_spearman_comparison_list)
        spearman_correlation_score, p_value = get_spearman_score(scores_for_spearman_comparison_list)

        #
        scores_list.append((protein_name, spearman_correlation_score, p_value))

    return (scores_list)


def get_count_of_lines_except_header(file_path):
    """Given a file path it counts the total number
    of lines and then returns n-1 because the first
    line could be a header. This function can also be
    used for the dbNSFP .err files as the last line
    says that the how many snps were not found"""
    with open(file_path, 'r') as file:
        # Skip the header
        header = next(file)

        # Count the remaining lines
        line_count = sum(1 for line in file)

    return line_count


def get_list_to_add_in_dataframe(list_from_the_dictionary):
    """
    it takes the list entry from the dictionary value, and then
    it returns the list with three values to enter add it to the
    dataframe.
    """
    sequence = remove_stop_codon(list_from_the_dictionary[-1])

    snp_string = ""
    scaled_effect_string = ""
    snps_and_scores = list_from_the_dictionary[:-1]

    for snp_and_score in snps_and_scores:
        snp = snp_and_score.split()[0][3:]
        scaled_effect = snp_and_score.split()[1].split(":")[1]
        snp_string += snp + ";"
        scaled_effect_string += scaled_effect + ";"

    return [snp_string[:-1], scaled_effect_string[:-1], sequence]

def get_protein_names_from_db_nsfp_output_directory(directory_path):
    """
    Get protein names from CSV files in the specified directory.

    Parameters:
    - directory_path (Path): Path to the directory containing CSV files.

    Returns:
    - list: A list of protein names extracted from the CSV file names.
    """

    # List all file names in the directory
    csv_file_names = [file.name for file in directory_path.iterdir() if
                      file.is_file() and file.name.endswith('.csv')]

    # Extract protein names from CSV file names
    protein_names_from_csv_files = [filename.replace('.csv', '') for filename in csv_file_names]

    return protein_names_from_csv_files

# TODO Move these function related to the caluclations of the pearson corelation to a new module.
def get_mave_tool_scores_dataframe(mave_scores_dict: dict,
                                   tool_scores_dict: dict,
                                   mave_score_dictionary_column_name: str,
                                   tool_score_dictionary_column_name: str) -> pd.DataFrame:
    """
    Create a DataFrame with MAVE and tool scores.

    Parameters:
    - mave_scores_dict (dict): Dictionary containing MAVE scores.
    - tool_scores_dict (dict): Dictionary containing tool scores.
    - mave_score_column_name (str): Name of the column for MAVE scores.
    - tool_score_column_name (str): Name of the column for tool scores.

    Returns:
    - pd.DataFrame: DataFrame with 'SNPs', 'MAVE_Scores', and 'Tool_Scores' columns.
    """

    # Create a set of all keys from both dictionaries
    all_keys = set(mave_scores_dict.keys()).union(tool_scores_dict.keys())

    # Create a DataFrame with None as the default value
    df = pd.DataFrame({key: [mave_scores_dict.get(key), tool_scores_dict.get(key)] for key in all_keys},
                      index=[mave_score_dictionary_column_name, tool_score_dictionary_column_name]).T.reset_index()

    # Rename the columns
    df.columns = [COLUMN_NAME_OF_MAVE_GOLD_STANDARD_SAVS, mave_score_dictionary_column_name, tool_score_dictionary_column_name]

    return df


def get_correlation_and_percentage_used(df, column1_name, column2_name):
    """
    Calculate the Pearson correlation and percentage of rows used between two columns in the DataFrame.

    Parameters:
    - df (pd.DataFrame): The input DataFrame.
    - column1_name (str): Name of the first column.
    - column2_name (str): Name of the second column.

    Returns:
    - Tuple[float, float]: A tuple containing the percentage of rows used and Pearson correlation.
    """
    # Drop rows with missing values in either of the columns
    valid_data = df.dropna(subset=[column1_name, column2_name])

    # Calculate correlation
    pearson = valid_data[column1_name].corr(valid_data[column2_name])
    spearman, _ = spearmanr(valid_data[column1_name], valid_data[column2_name])
    correlation = spearman


    # Calculate the percentage of rows used
    percentage_used = (len(valid_data) / len(df)) * 100

    return percentage_used, correlation

def pickle_dataframe(dataframe, file_path, file_name):
    """
    Pickle a DataFrame to a specified file path and name.

    Parameters:
    - dataframe (pd.DataFrame): The DataFrame to pickle.
    - file_path (Path): The path where the pickle file will be saved.
    - file_name (str): The name of the pickle file.

    Returns:
    - None
    """
    # Ensure the file path exists, create it if not
    file_path.mkdir(parents=True, exist_ok=True)

    # Pickle the DataFrame
    with open(file_path / file_name, 'wb') as file:
        pickle.dump(dataframe, file)



def load_dataframe(file_path, file_name):
    """
    Load a pickled DataFrame from a specified file path and name.

    Parameters:
    - file_path (Path): The path where the pickle file is located.
    - file_name (str): The name of the pickle file.

    Returns:
    - pd.DataFrame: The loaded DataFrame.
    """
    # Check if the file exists
    file_to_load = file_path / file_name
    if not file_to_load.is_file():
        raise FileNotFoundError(f"File '{file_name}' not found in path '{file_path}'.")

    # Load the DataFrame
    with open(file_to_load, 'rb') as file:
        loaded_dataframe = pickle.load(file)

    return loaded_dataframe


def filter_dataframe_by_species(df,
                                target_species,
                                species_column=COLUMN_NAME_OF_MAVE_GOLD_STANDARD_SPECIES):
    """
    Filters a DataFrame based on a specified species.

    Parameters:
    - df (pd.DataFrame): The DataFrame to be filtered.
    - species_column (str): The name of the column containing species information.
    - target_species (str): The species to filter the DataFrame by.

    Returns:
    pd.DataFrame: Filtered DataFrame containing only rows with the specified species.
    """
    return df[df[species_column] == target_species]


def add_missing_columns(dataframe1, dataframe2):
    """
    Dynamically adds missing columns from dataframe2 to dataframe1 and copies values.

    Parameters:
    - dataframe1 (pd.DataFrame): Target DataFrame.
    - dataframe2 (pd.DataFrame): DataFrame with additional columns.

    Returns:
    pd.DataFrame: DataFrame1 with missing columns added and values copied.
    """
    missing_columns = set(dataframe2.columns) - set(dataframe1.columns)

    for column in missing_columns:
        dataframe1[column] = dataframe2[column]

    return dataframe1


def add_column_from_tool_df_to_mave_df(mave_df,
                                       tool_df,
                                       mave_df_prot_seq_col_name,
                                       tool_df_prot_seq_col_name,
                                       tool_col_to_add,
                                       name_of_new_col):
    """
    Merge a column from a tool DataFrame to a MAVE DataFrame based on a shared protein sequence column.

    Parameters:
    - mave_df (pd.DataFrame): The MAVE DataFrame to which the column will be added.
    - tool_df (pd.DataFrame): The tool DataFrame containing the additional column.
    - mave_df_prot_seq_col_name (str): The column name in the MAVE DataFrame representing protein sequences.
    - tool_df_prot_seq_col_name (str): The column name in the tool DataFrame representing protein sequences.
    - tool_col_to_add (str): The name of the column from the tool DataFrame to be added to the MAVE DataFrame.

    Returns:
    - pd.DataFrame: Merged DataFrame with the specified tool column added.

    Note:
    Rows in the resulting DataFrame will be aligned based on matching protein sequence values.
    """
    # Check if inputs are DataFrames
    if not isinstance(mave_df, pd.DataFrame) or not isinstance(tool_df, pd.DataFrame):
        raise TypeError("Both mave_df and tool_df must be pandas DataFrames.")

    # Check if specified columns exist in respective DataFrames
    if mave_df_prot_seq_col_name not in mave_df.columns or tool_df_prot_seq_col_name not in tool_df.columns:
        raise ValueError("Specified protein sequence columns do not exist in their respective DataFrames.")

    # Check if the specified column to add exists in tool_df
    if tool_col_to_add not in tool_df.columns:
        raise ValueError(f"Specified column to add ({tool_col_to_add}) does not exist in tool_df.")

    # Perform the merge
    merged_df = pd.merge(mave_df, tool_df[[tool_df_prot_seq_col_name, tool_col_to_add]],
                         left_on=mave_df_prot_seq_col_name, right_on=tool_df_prot_seq_col_name, how='left')

    # Drop duplicated protein sequence column
    merged_df = merged_df.drop(columns=[tool_df_prot_seq_col_name])

    ## Rename the column
    merged_df = merged_df.rename(columns={tool_col_to_add: name_of_new_col})

    return merged_df


def add_flag_column(df, target_column, flag_column_name):
    """
    Add a binary flag column to a DataFrame based on the presence of NaN values in a target column.

    Parameters:
    - df (pd.DataFrame): The DataFrame to which the flag column will be added.
    - target_column (str): The name of the target column to check for NaN values.
    - flag_column_name (str): The name of the new binary flag column.

    Returns:
    - pd.DataFrame: DataFrame with the new binary flag column added.
    """
    if df.empty:
        df[flag_column_name] = 0
    else:
        df[flag_column_name] = df[target_column].notnull().astype(int)
    return df


def exclude_snps(df, exclude_snps_list, savs_string=SAVS_STRING):
    """
    Exclude rows from a dataframe based on values in the 'snps' column.

    Parameters:
    - df (pd.DataFrame): The input dataframe.
    - exclude_snps_list (list): List of SNPs to exclude.

    Returns:
    - pd.DataFrame: A new dataframe with rows excluded based on the 'snps' column.
    """

    if np.any(pd.isna(exclude_snps_list)):
        return df  # Return the same DataFrame if exclude_snps_list is nan

    return df[~df[savs_string].isin(exclude_snps_list)]


def generate_tool_columns(tool_name,
                          exclude_tool_training_snps_flag,
                          pearson_suffix=SPEAR_COR_SUFFIX,
                          snps_percentage_suffix=USED_SAV_PERCENTAGE_SUFFIX,
                          score_suffix=TOOL_SCORE_COLUMN_SUFFIX,
                          exlude_training_snp_suffix=EXCLUDE_TRAINING_SAV_SUFFIX):
    """
    Generates column names for a tool based on the provided tool name and suffixes.

    Parameters:
    - tool_name (str): Name of the tool.
    - pearson_suffix (str): Suffix for the Pearson correlation column.
    - snps_percentage_suffix (str): Suffix for the SNP percentage column.
    - score_suffix (str): Suffix for the tool score column.

    Returns:
    - tuple: Tuple containing the generated column names for Pearson correlation, SNP percentage, and tool score.
    """

    if exclude_tool_training_snps_flag:
        pearson_column = tool_name + pearson_suffix +exlude_training_snp_suffix
        snps_percentage_column = tool_name + snps_percentage_suffix + exlude_training_snp_suffix
        score_column = tool_name + score_suffix

    else:
        pearson_column = tool_name + pearson_suffix
        snps_percentage_column = tool_name + snps_percentage_suffix
        score_column = tool_name + score_suffix

    return pearson_column, snps_percentage_column, score_column

def get_training_snps_column_name(tool_name, tool_training_suffix=TRAINING_SAVS_COLUMN_SIFFIX):
    """
    Generates the column name for training SNPs based on the provided tool name.

    Parameters:
    - tool_name (str): Name of the tool.
    - tool_training_suffix (str, optional): Suffix for the training SNPs column.

    Returns:
    - str: The generated column name for training SNPs.
    """
    return tool_name + tool_training_suffix

def convert_column_to_list(df, column_name):
    """
    Convert a column with comma-separated strings to lists in a DataFrame.

    Parameters:
    - df (pd.DataFrame): Input DataFrame.
    - column_name (str): Name of the column to be converted.

    Returns:
    - pd.DataFrame: Updated DataFrame with the specified column values as lists.
    """
    # Copy the DataFrame to avoid modifying the original
    df_copy = df.copy()

    # Apply the transformation to the specified column
    df_copy[column_name] = df_copy[column_name].apply(lambda x: x.split(','))

    return df_copy

def update_value_based_on_protein_name(df, column_name, value, protein_name,
                                       id_column_name=COLUMN_NAME_OF_MAVE_GOLD_STANDARD_ID):
    """
    Update a value in a specific column for a row that matches the protein_name.

    Parameters:
    - df (pd.DataFrame): The DataFrame to update.
    - column_name (str): The name of the column to update.
    - value: The new value to set.
    - protein_name (str): The name of the protein to match.

    Returns:
    - pd.DataFrame: The updated DataFrame.
    """
    # if isinstance(value, dict):
    #     value = json.dumps(value)
    df.loc[df[id_column_name] == protein_name, column_name] = value
    return df

def get_value_from_dataframe(df, id_column, id_value, value_column):
    """
    Retrieves a value from a DataFrame based on the provided id_column and id_value.

    Parameters:
    - df (pd.DataFrame): The DataFrame to retrieve the value from.
    - id_column (str): The name of the column to match the id_value against.
    - id_value: The value to match in the id_column.
    - value_column (str): The name of the column to retrieve the value from.

    Returns:
    - The value from the value_column of the first row where id_column matches id_value.
    """
    return df.loc[df[id_column] == id_value, value_column].values[0]

def generate_amino_pssm_dict(amino_acids=AMINO_ACIDS_SINGLE_LETTER):
    """
    This dictionary that mimicks the blank PSSM matrix of amino acid substitutions.
    Input: amino_acids (list): List of amino acids to generate the dictionary for.
    Output: dict: Dictionary with amino acid pairs as keys and new empty list as values.
    """
    return {aa1+aa2: [] for aa1, aa2 in itertools.product(amino_acids, repeat=2)}

def remove_digits_from_key(key):
    """
    This function removes digits from a string key.
    """
    return ''.join([char for char in key if char.isalpha()])

def calculate_mean(lst):
    """
    Calculate the mean of a list.

    Parameters:
    - lst (list): The list to calculate the mean of.

    Returns:
    - float: The mean of the list if the list is not empty.
    - None: If the list is empty.
    """
    if len(lst) == 0:
        return None
    else:
        return sum(lst) / len(lst)


def get_protein_name_list(df, column_name):
    """
    Extracts a list of protein names from a DataFrame based on the provided column name.

    Parameters:
    - df (pd.DataFrame): The DataFrame to extract the protein names from.
    - column_name (str): The name of the column containing the protein names.

    Returns:
    - list: A list of protein names.
    """
    return df[column_name].tolist()

def get_protein_list(dataframe, column_name=COLUMN_NAME_OF_MAVE_GOLD_STANDARD_ID):
    """
    This function returns a list of proteins from a given dataframe.
    Input: dataframe, column_name
    Output: list of proteins
    """
    return dataframe[column_name].tolist()

def extract_value(df, id_column, row_value, target_col_name):
    """
    Extract a value from a DataFrame given the id column name, row value, and target column name.

    Parameters:
        df (DataFrame): The DataFrame containing the data.
        id_column (str): The name of the id column.
        row_value (str): The value in the id column for the row of interest.
        target_col_name (str): The name of the column from which to extract the value.

    Returns:
        The value in column1 for the row where id_column equals row_value. If no such row exists, returns None.
    """
    row = df[df[id_column] == row_value]
    if not row.empty:
        value = row[target_col_name].values[0]
        try:
            return len(value)
        except TypeError:
            return 0
    else:
        return 0


if __name__ == '__main__':
    pass