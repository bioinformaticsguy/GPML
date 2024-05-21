from pathlib import Path

COLUMN_NAME_OF_MAVE_GOLD_STANDARD_ID = 'protein_name'
COLUMN_NAME_OF_MAVE_GOLD_STANDARD_SNP = 'snps'
COLUMN_NAME_OF_MAVE_GOLD_STANDARD_SCALED_EFFECT = 'scaled_effect'
COLUMN_NAME_OF_MAVE_GOLD_STANDARD_PROTEIN_SEQUENCE = 'Prot_sequence'

COLUMN_NAMES_OF_MAVE_GS_DATAFRAME_LIST = [COLUMN_NAME_OF_MAVE_GOLD_STANDARD_ID,
                                          COLUMN_NAME_OF_MAVE_GOLD_STANDARD_SNP,
                                          COLUMN_NAME_OF_MAVE_GOLD_STANDARD_SCALED_EFFECT,
                                          COLUMN_NAME_OF_MAVE_GOLD_STANDARD_PROTEIN_SEQUENCE]


SPECIE_NAME_HUMAN = "Human"


COLUMN_NAME_OF_MAVE_GOLD_STANDARD_SPECIES = "species"
COLUMN_NAME_OF_MAVE_GOLD_STANDARD_SNP_DICTIONARY = "mave_snps_scores_dictinary"
COLUMN_NAME_OF_BASELINE_SCORES_DICTIONARY = "baseline_scores_dictionary"
MUTEPRED_TOOL_NAME = "MutPred"

MUTEPRED_AMINO_ACID_SUBSTITUTIONS_COLUMN_NAME = "Amino_acid_substitutions"

SNP_COLUMN_NAME = "HGVSp_ANNOVAR"
AMINO_ACID_SEQUENCE_COLUMN_NAME = 'Prot_sequence'
TOOLS_LIST = ["MutPred",
              "REVEL",
              "PROVEAN",
              "Polyphen2_HDIV",
              "EVE",
              "AlphaMissense",
              ]

AMINO_ACIDS_SINGLE_LETTER = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

# TOOLS_LIST = ["MutPred"]

## Suffixes and Prefixes
TOOL_SCORE_COLUMN_SUFFIX = "_score"
PEARSON_CORELATION_SUFFIX = "_pearson_correlation"
USED_SNP_PERCENTAGE_SUFFIX = "_used_snp_percentage"
TRAINING_FLAG_SUFFIX = "_training_flag"
TRAINING_SNPS_COLUMN_SIFFIX = "_training_snps"
EXCLUDE_TRAINING_SNP_SUFFIX = "_excluded_training_snps"

## File Names and Paths
MAVE_DATAFRAME_PICKLE_FILE_NAME = "gold_std_df.pkl"
MAVE_DATAFRAME_ONLY_HUMAN_PICKLE_FILE_NAME = "gold_std_df_only_human.pkl"
MAVE_DATAFRAME_ONLY_HUMAN_WITH_BASELINE_PICKLE_FILE_NAME = "gold_std_df_only_human_with_baseline.pkl"
MAVE_DATAFRAME_HUMAN_LOPO_MEAN_PICKLE_FILE_NAME = "gold_std_df_human_lopo_mean.pkl"
MAVE_DATAFRAME_HUMAN_LOPO_MEAN_PICKLE_FILE_NAME_WITH_TOOL_SCORES = "gold_std_df_human_lopo_mean_with_tool_scores.pkl"

PICKLED_DATAFRAMES_DIRECTORY_PATH = Path("Data/pickled_dataframes")
TRAINING_DATA_FILE_PATH = Path("Data/mutepred_training_data/wo_exclusive_hgmd_mp2_training_data_MavedbData.csv")
OUTPUT_DIR_DB_NSFP = Path("Data/dbNSFP_output_dir")


PROTEIN_SPECIES_MAPPING = {
    'A0A2Z5U3Z0_9INFA_A0A2Z5U3Z0_9INFA_Doud_2016': 'Virus',
    'BLAT_ECOLX_BLAT_ECOLX_Deng_2012': 'E-Coli',
    'BLAT_ECOLX_BLAT_ECOLX_Jacquier_2013': 'E-Coli',
    'CBS_urn:mavedb:00000005-a': 'Human',
    'CCDB_ECOLI_CCDB_ECOLI_Tripathi_2016': 'E-Coli',
    'CcdB_urn:mavedb:00000084-a': 'E-Coli',
    'CCR5_urn:mavedb:00000047-c': 'Human',
    'CYP2C9_urn:mavedb:00000095-a': 'Human',
    'CYP2C9_urn:mavedb:00000095-b': 'Human',
    'IF1_ECOLI_IF1_ECOLI_Kelsic_2016': 'E-Coli',
    'KKA2_KLEPN_KKA2_KLEPN_Melnikov_2014': 'Bacteria',
    'NUDT15_urn:mavedb:00000055-0': 'Human',
    'NUDT15_urn:mavedb:00000055-a': 'Human',
    'p53_urn:mavedb:00000059-a': 'Human',
    'PTEN_urn:mavedb:00000013-a': 'Human',
    'PTEN_urn:mavedb:00000054-a': 'Human',
    'Q2N0S5_9HIV1_Q2N0S5_9HIV1_Haddox_2018': 'Virus',
    'R1AB_SARS2_R1AB_SARS2_Flynn_growth_2022': 'Virus',
    'RL401_YEAST_RL401_YEAST_Mavor_2016': 'Yeast',
    'SARS-CoV-2_receptor_binding_domain_urn:mavedb:00000044-a': 'Virus',
    'SARS-CoV-2_receptor_binding_domain_urn:mavedb:00000044-b': 'Virus',
    'SUMO1_urn:mavedb:00000001-b': 'Human',
    'TEM-1_beta-lactamase_urn:mavedb:00000070-a': 'Escherichia coli',
    'TEM-1_beta-lactamase_urn:mavedb:00000086-c': 'Escherichia coli',
    'TEM-1_beta-lactamase_urn:mavedb:00000086-d': 'Escherichia coli',
    'TEM-1_beta-lactamase_urn:mavedb:00000086-e': 'Escherichia coli',
    'TEM-17_beta-lactamase_urn:mavedb:00000085-b': 'Bacteria',
    'TP53_(P72R)_urn:mavedb:00000068-c': 'Human',
    'VIM-2_with_p.Met1_Phe2insGly_urn:mavedb:00000073-a': 'Bacteria',
    'VIM-2_with_p.Met1_Phe2insGly_urn:mavedb:00000073-c': 'Bacteria',
    'VIM-2_with_p.Met1_Phe2insGly_urn:mavedb:00000073-d': 'Bacteria',
    'VIM-2_with_p.Met1_Phe2insGly_urn:mavedb:00000073-f': 'Bacteria',
    'VIM-2_with_p.Met1_Phe2insGly_urn:mavedb:00000073-g': 'Bacteria',
    'VIM-2_with_p.Met1_Phe2insGly_urn:mavedb:00000073-h': 'Bacteria',
    'VKOR_urn:mavedb:00000078-a': 'Human',
    'VKOR_urn:mavedb:00000078-b': 'Human'
}

PROTEIN_SHORT_MAPPING = {
    'A0A2Z5U3Z0_9INFA_A0A2Z5U3Z0_9INFA_Doud_2016': '9INFA',
    'BLAT_ECOLX_BLAT_ECOLX_Deng_2012': 'BLAT_ECOL_2012',
    'BLAT_ECOLX_BLAT_ECOLX_Jacquier_2013': 'BLAT_ECOL_2013',
    'CBS_urn:mavedb:00000005-a': 'CBS_05-a',
    'CCDB_ECOLI_CCDB_ECOLI_Tripathi_2016': 'CCDB_ECOLI_2016',
    'CcdB_urn:mavedb:00000084-a': 'CcdB_84-a',
    'CCR5_urn:mavedb:00000047-c': 'CCR547-c',
    'CYP2C9_urn:mavedb:00000095-a': 'CYP2C9_95-a',
    'CYP2C9_urn:mavedb:00000095-b': 'CYP2C9_95-b',
    'IF1_ECOLI_IF1_ECOLI_Kelsic_2016': 'IF1_ECOLI_2016',
    'KKA2_KLEPN_KKA2_KLEPN_Melnikov_2014': 'KKA2_KLEPN_2014',
    'NUDT15_urn:mavedb:00000055-0': 'NUDT15_55-0',
    'NUDT15_urn:mavedb:00000055-a': 'NUDT15_55-a',
    'p53_urn:mavedb:00000059-a': 'p53_59-a',
    'PTEN_urn:mavedb:00000013-a': 'PTEN_13-a',
    'PTEN_urn:mavedb:00000054-a': 'PTEN_54-a',
    'Q2N0S5_9HIV1_Q2N0S5_9HIV1_Haddox_2018': 'Q2N0S5_9HIV1_2018',
    'R1AB_SARS2_R1AB_SARS2_Flynn_growth_2022': 'R1AB_SARS2_2022',
    'RL401_YEAST_RL401_YEAST_Mavor_2016': 'RL401_YEAST_2016',
    'SARS-CoV-2_receptor_binding_domain_urn:mavedb:00000044-a': 'SARS-CoV-2_44-a',
    'SARS-CoV-2_receptor_binding_domain_urn:mavedb:00000044-b': 'SARS-CoV-2_44-b',
    'SUMO1_urn:mavedb:00000001-b': 'SUMO1_01-b',
    'TEM-1_beta-lactamase_urn:mavedb:00000070-a': 'TEM-1_beta-70-a',
    'TEM-1_beta-lactamase_urn:mavedb:00000086-c': 'TEM-1_beta-86-c',
    'TEM-1_beta-lactamase_urn:mavedb:00000086-d': 'TEM-1_beta-86-d',
    'TEM-1_beta-lactamase_urn:mavedb:00000086-e': 'TEM-1_beta-86-e',
    'TEM-17_beta-lactamase_urn:mavedb:00000085-b': 'TEM-17_beta-85-b',
    'TP53_(P72R)_urn:mavedb:00000068-c': 'TP53_P72R_68-c',
    'VIM-2_with_p.Met1_Phe2insGly_urn:mavedb:00000073-a': 'VIM-2_73-a',
    'VIM-2_with_p.Met1_Phe2insGly_urn:mavedb:00000073-c': 'VIM-2_73-c',
    'VIM-2_with_p.Met1_Phe2insGly_urn:mavedb:00000073-d': 'VIM-2_73-d',
    'VIM-2_with_p.Met1_Phe2insGly_urn:mavedb:00000073-f': 'VIM-2_73-f',
    'VIM-2_with_p.Met1_Phe2insGly_urn:mavedb:00000073-g': 'VIM-2_73-g',
    'VIM-2_with_p.Met1_Phe2insGly_urn:mavedb:00000073-h': 'VIM-2_73-h',
    'VKOR_urn:mavedb:00000078-a': 'VKOR_78-a',
    'VKOR_urn:mavedb:00000078-b': 'VKOR_78-b'
}


if __name__ == '__main__':
    pass