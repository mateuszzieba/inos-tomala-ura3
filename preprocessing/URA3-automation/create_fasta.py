import pandas as pd

df_missense_possible = pd.read_csv("ura3_all_missense_possible.tsv",
                                   sep="\t")[["SYMBOL", "Gene", "Protein_position", "Amino_acids", "SIFT"]]


df_missense_possible = df_missense_possible[df_missense_possible.SIFT.str.contains('^tolerated')].reset_index(drop=True)


sequence = open('URA3_YEL021W.fas', 'r').readlines()
# print([sequence[0], sequence[1]])
ura3_sequence = str(sequence[1])


to_fasta = ""
for i in range(0, len(df_missense_possible)):
    position_change = df_missense_possible.iloc[i, 2] - 1
    aa_change = df_missense_possible.iloc[i, 3]
    
    name_sequence = ">" + aa_change[0] + str(df_missense_possible.iloc[i, 2]) + aa_change[2] + ".seq" + "\n"
    new_sequence = ura3_sequence[:position_change] + aa_change[2] + ura3_sequence[position_change + 1:] + "\n"
    to_fasta = to_fasta + name_sequence + new_sequence

# print(to_fasta[:-1])

file_fasta = open('ura3_all_missense_possible_sequences.fasta', 'w')
file_fasta.write(to_fasta)
file_fasta.close()

