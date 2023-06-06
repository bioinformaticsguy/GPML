s1, s2, s3, s4 = 0,0,0,0 

with open('Protein_mave_db_gold_standard_only_sequences_splits.tsv') as file:
    for line in file:
        if line[-8:] == 'Split_1\n':
            s1+=1
            # print(line)
        if line[-8:] == 'Split_2\n':
            s2+=1
            # print(line)
        if line[-8:] == 'Split_3\n':
            s3+=1
            # print(line)
        if line[-8:] == 'Split_4\n':
            s4+=1
            # print(line)


print(s1, s2, s3, s4)