# s1, s2, s3, s4 = 0,0,0,0 

# with open('Protein_mave_db_gold_standard_only_sequences_splits.tsv') as file:
#     for line in file:
#         if line[-8:] == 'Split_1\n':
#             s1+=1
#             # print(line)
#         if line[-8:] == 'Split_2\n':
#             s2+=1
#             # print(line)
#         if line[-8:] == 'Split_3\n':
#             s3+=1
#             # print(line)
#         if line[-8:] == 'Split_4\n':
#             s4+=1
#             # print(line)


# print(s1, s2, s3, s4)


# ## Testing line


data_directory = "/home/ali/Documents/GPML/Data/"

mave_folder_name = "maveGSData/"

file_name = "mave_db_gold_standard.fasta"

file_path = data_directory + mave_folder_name + file_name

list = []
with open(file_path) as file:
    for line in file:
        if line[0] == ">":
            list.append(line)
            # print(line[:-1])

print(len(list))