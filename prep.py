from pprint import pprint as pretty_print
import scipy.stats as stats



def get_dictonary_of_scores(file_path):
    """
    This function takes the file which has the mutations and the sequences, 
    then it stores the date related to each protein in one list and then, 
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



# pretty_print(point_mutations)

def write_from_list(list, file_name):
    with open(file_name, 'w') as file:
        file.write("mutations\n")
        for i in list:
            file.write(i+'\n')
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



def get_dictionary_for_mutpred_scores(score_file_path):
    scores_dictionary = {}
    with open(score_file_path) as file:
        for line in file:
            # break
            if line[:2] == "ID":
                pass
            else:
                # print(line.split(","))

                mutation_name = line.split(",")[1]
                mutepred_score = line.split(",")[2]
                scores_dictionary[mutation_name] = [mutepred_score]

                # print(mutation_name, mutepred_score)
                # break
    #             # print(line.split(","))
    #             eve_score = line.split(",")[3]


    #             # print(mutation_name, evol_index, eve_score)


    return scores_dictionary


score_file_path = "/home/ali/Documents/GPML/tVIM2.out"
get_dictionary_for_mutpred_scores(score_file_path)




def get_score_comparison_list(eve_score_dict, all_data_dict, protein_name):
    print(protein_name)
    scores_list = []
    count = 0
    protein_data = all_data_dict[protein_name]


    # pretty_print(protein_data)

    for line in protein_data:
        # print(line[0])
        if line[0] == "<":

            mutation = get_single_letter_point_mutation(line.split()[0][3:])
            scaled_effect = float(line.split(":")[1])
            # pretty_print(mutation)

            try:
                eve_score_value = eve_score_dict[mutation][0]
                scores_list.append([mutation, scaled_effect, eve_score_value])
                # print(eve_score_value)
            except KeyError:
                count += 1 
                # print("Key not found in dictionary")

        

    return scores_list







def get_spearman_score(scores_for_spearman_comparison):
    scaled_effects = []
    eve_score_values = []
    for score_list in scores_for_spearman_comparison:
        scaled_effects.append(score_list[1])
        eve_score_values.append(score_list[2])
        # print(score_list[1])

    # Calculate Spearman correlation coefficient and p-value
    corr, p_value = stats.spearmanr(scaled_effects, eve_score_values)

    print("Spearman correlation coefficient:", corr)
    print("p-value:", p_value)




# protein_name = "p53_urn:mavedb:00000059-a"

data_directory = "/home/ali/Documents/GPML/Data/"

cor_folder = "spearman_coor_files/"

mave_folder_name = "maveGSData/"

file_name = "mave_db_gold_standard.fasta"


file_path = data_directory + mave_folder_name + file_name

# eve_scores_file_name = "all_EVE_scores_Jun17_VKOR1_example.csv"

# eve_scores_file_name = "all_EVE_scores_June_P53_HUMAN.csv"


# eve_scores_file_name = "all_EVE_scores_Jun8_PTEN.csv"


eve_scores_file_name = "all_EVE_scores_June17_CBS_HUMAN.csv"



# protein_name = "VKOR_urn:mavedb:00000078-a"

# protein_name = "p53_urn:mavedb:00000059-a"

# protein_name = "PTEN_urn:mavedb:00000013-a"

protein_name = "VIM-2_with_p.Met1_Phe2insGly_urn:mavedb:00000073-c\t\turn:mavedb:00000073-c,None"




# mutations_file_name = "PTEN_HUMAN_all_singles.csv"
# mutations_file_name = "all_EVE_scores_Jun17_VKOR1_example.csv"



# mutations_file_path = data_directory + mutations_file_name


all_data_dict = get_dictonary_of_scores(file_path)
# pretty_print(all_data_dict['VIM-2_with_p.Met1_Phe2insGly_urn:mavedb:00000073-c'][-1])
# pretty_print(type(all_data_dict['VIM-2_with_p.Met1_Phe2insGly_urn:mavedb:00000073-c']))
# pretty_print(all_data_dict.keys())

pro_name = "NUDT15_urn:mavedb:00000055-0		urn:mavedb:00000055-0,None"

VIM2 = all_data_dict[pro_name]


def remove_stop_codon(sequence):
    if sequence[-1] == "*":
        return sequence[:-1]
    else:
        return sequence

def mute_pred_input(pro_name, data, mut_count=None):
    """Input
            pro_name: string --> VIM-2_with_p.Met1_Phe2insGly_urn:mavedb:00000073-c
            data:  list of strings where last string is the sequence --> see the tail 
                            ['<p.Met1Ala #scaled_effect:0.29943354034999425',
                            'MGFKPVGEVRLYQI']
            mut_count: intege to limit the number of mutations to add infront of name
                        of the protein."""
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
    
    

# data = VIM2
# mut_count = None



def write_file_for_mutepred(name_with_ext, data_touple):
    with open(name_with_ext, 'w') as file:
        for item in data_touple:
            # print(item)
            file.write(item)


# name_with_ext = "NUDT15.fasta"
# data_touple = mute_pred_input(pro_name, data, mut_count)


# write_file_for_mutepred(name_with_ext, data_touple)

# # P53_data = all_data_dict['p53_urn:mavedb:00000059-a']
# # point_mutations = get_mutations_list(P53_data)

# eve_scores_directory = "/home/ali/Documents/GPML/tVIM2.out"
# eve_score_dict = get_dictionary_for_mutpred_scores(eve_scores_directory)

# scores_for_spearman_comparison = get_score_comparison_list(eve_score_dict, all_data_dict, protein_name)

# pretty_print(scores_for_spearman_comparison)

# get_spearman_score(scores_for_spearman_comparison)




# 

# # P53_data = all_data_dict['p53_urn:mavedb:00000059-a']
# # point_mutations = get_mutations_list(P53_data)

# eve_score_dict = get_dictionary_for_eve_scores(eve_scores_directory)

# scores_for_spearman_comparison = get_score_comparison_list(eve_score_dict, all_data_dict, protein_name)

# # pretty_print(scores_for_spearman_comparison)

# get_spearman_score(scores_for_spearman_comparison)



def readFastaIteration(filename):
    '''Given a fasta file returns a list of all the touples
    in which the first element of touple is the description of the fasta sequence 
    and the 2nd element is the sequence itself.'''
    sequences = []
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

seqs = list(readFastaIteration('/home/ali/Documents/GPML/Data/maveGSData/mave_db_gold_standard_only_sequences.fasta').values())

# print(len(seqs) == len(set(seqs)))
# print(len(seqs))
# print(len(set(seqs)))


def find_duplicate_values(dictionary):
    value_to_keys = {}
    duplicates = {}

    for key, value in dictionary.items():
        if value in value_to_keys:
            value_to_keys[value].append(key)
        else:
            value_to_keys[value] = [key]

    for value, keys in value_to_keys.items():
        if len(keys) > 1:
            duplicates[value] = keys

    return duplicates

# Your dictionary
# your_dictionary = {
#     "key1": "value1",
#     "key2": "value2",
#     "key3": "value1",
#     "key4": "value3",
#     "key5": "value2"
# }

your_dictionary = readFastaIteration('/home/ali/Documents/GPML/Data/maveGSData/mave_db_gold_standard_only_sequences.fasta')

duplicate_keys = find_duplicate_values(your_dictionary)

pretty_print(duplicate_keys)

if not duplicate_keys:
    print("No duplicate values found.")
else:
    print("Keys with duplicate values:")
    for value, keys in duplicate_keys.items():
        print(f"Value '{value[:5]}' appears in keys:", [key[:5] for key in keys])