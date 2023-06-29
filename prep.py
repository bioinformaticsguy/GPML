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





def get_spearman_score(scores_for_spearman_comparison):
    scaled_effects = []
    eve_score_values = []
    for score_list in scores_for_spearman_comparison:
        scaled_effects.append(score_list[1])
        eve_score_values.append(score_list[2])
        print(score_list[1])

    # Calculate Spearman correlation coefficient aid and then the number, 
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
                eve_score_value = eve_score_dict[mutation][1]
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

eve_scores_directory = data_directory + cor_folder + eve_scores_file_name


# protein_name = "VKOR_urn:mavedb:00000078-a"

# protein_name = "p53_urn:mavedb:00000059-a"

# protein_name = "PTEN_urn:mavedb:00000013-a"

protein_name = "CBS_urn:mavedb:00000005-a"




# mutations_file_name = "PTEN_HUMAN_all_singles.csv"
# mutations_file_name = "all_EVE_scores_Jun17_VKOR1_example.csv"



# mutations_file_path = data_directory + mutations_file_name


all_data_dict = get_dictonary_of_scores(file_path)


# P53_data = all_data_dict['p53_urn:mavedb:00000059-a']
# point_mutations = get_mutations_list(P53_data)

eve_score_dict = get_dictionary_for_eve_scores(eve_scores_directory)

scores_for_spearman_comparison = get_score_comparison_list(eve_score_dict, all_data_dict, protein_name)

# pretty_print(scores_for_spearman_comparison)

get_spearman_score(scores_for_spearman_comparison)




