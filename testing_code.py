import pandas as pd

# Your dataframe
data = {'snps': ["K359E", "D234P", "T141M", "H513D", "T146I"],
        'mave_snps_scores_dictinary': [0.8348246019731457, 0.014415572740724913, 1.2590526645570985, 0.5783830466766263, 0.12149939270369847],
        'MutPred_score': [-0.44, None, -4.56, -2.79, -5.97]}

df = pd.DataFrame(data)

# List of SNPs to exclude
exclude_snps = ["K359E", "D234P", "H513D"]

# Exclude rows based on the "snps" column
df_filtered = df[~df['snps'].isin(exclude_snps)]

print(df_filtered)