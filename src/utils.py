import scipy.stats as stats


def get_dictonary_of_scores(file_path):
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


def get_score_comparison_list(tool_score_dict, gs_dictionary, protein_name):
    print(protein_name)
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
                tool_score_value = tool_score_dict[mutation][0]
                scores_list.append([mutation, scaled_effect, tool_score_value])
                # print(eve_score_value)
            except KeyError:
                count += 1
                # print("Key not found in dictionary")

    return scores_list


def get_spearman_score(scores_for_spearman_comparison_list):
    """
    Input:
            This function takes the scores list from the fucntion
            get_score_comparison_list()

    Return:
            This function returns the spearman score for the specific
            list.

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


def mute_pred_input(pro_name, data, mut_count=None):
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
        list_of_mute_pred_inputs.append(mute_pred_input(protein_name, data))

    return list_of_mute_pred_inputs
