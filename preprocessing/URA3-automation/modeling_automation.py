##########################
# import needed packages #
##########################
import math
import os
import re
import shutil
# from modeller import soap_protein_od
import sys

import pandas as pd
from modeller import *
from modeller.automodel import *


#############################################################################
# Functions needed to calculate parameters in the active site of the enzyme #
#############################################################################
def readPDB(file):
    """Function reads a pdb file and selects ATOMS"""
    atoms = {}
    f = open(file)
    for line in f:
        if line[0:6].replace(' ', '') == 'ATOM':
            chain = line[21:22]
            # print chain
            aa = line[17:20]
            pos = line[22:26].replace(' ', '')
            atom = line[12:16].replace(' ', '')
            id = chain + '*' + aa + '*' + pos + '*' + atom
            x = float(line[30:38].replace(' ', ''))
            y = float(line[38:46].replace(' ', ''))
            z = float(line[46:54].replace(' ', ''))
            atoms[id] = (x, y, z)
    return atoms


def getPOS(atoms, chain, pos):
    """Tunction gets the positions of atoms for specific amiono acids fro the pdb"""
    k = []
    coordinates = []
    for i in atoms.keys():
        if i.split('*')[0] == chain and i.split('*')[2] == pos:
            k.append(i)
            # print(k)
            for j in k:
                # print(j)
                c = atoms[j]
                # print(c)
                coordinates.append(c)
    return coordinates


def calcCENT(coordinates):
    x = 0
    y = 0
    z = 0
    for i in coordinates:
        x += i[0]
        y += i[1]
        z += i[2]
    x_cent = x / len(coordinates)
    y_cent = y / len(coordinates)
    z_cent = z / len(coordinates)
    return x_cent, y_cent, z_cent


def calcMIN(coordinates1, coordinates2):
    D = []
    for i in coordinates1:
        for j in coordinates2:
            d = math.sqrt((i[0] - j[0]) ** 2 + (i[1] - j[1]) ** 2 + (i[2] - j[2]) ** 2)
            D.append(d)
    return min(D)


def calcMINCEN(cen1, cen2):
    d = math.sqrt((cen1[0] - cen2[0]) ** 2 + (cen1[1] - cen2[1]) ** 2 + (cen1[2] - cen2[2]) ** 2)
    return d


def calculate_parameters(file, pair_aa):
    """Return list comparison parameters between two amino acids."""
    atoms = readPDB(file)
    coordinates1 = getPOS(atoms=atoms, chain=pair_aa[0][0], pos=pair_aa[0][2])
    aa1 = [key.split('*')[1] for (key, value) in atoms.items() if value == coordinates1[0]][0]
    coordinates2 = getPOS(atoms=atoms, chain=pair_aa[1][0], pos=pair_aa[1][2])
    aa2 = [key.split('*')[1] for (key, value) in atoms.items() if value == coordinates2[0]][0]
    cen1 = calcCENT(coordinates1)
    cen2 = calcCENT(coordinates2)
    calcMINCEN(cen1, cen2)
    calcMIN(coordinates1, coordinates2)

    return [file, calcMINCEN(cen1, cen2), calcMIN(coordinates1, coordinates2),
            aa1, pair_aa[0][2],
            aa2, pair_aa[1][2]]


#################
# build_profile #
#################
def build_profile(file):
    """Function build profile for sequence"""
    # Prepare the input files
    # Read in the sequence database
    file_bin = os.path.join(os.getcwd(), file[:-4], 'pdb_95.bin')
    file_prf = os.path.join(os.getcwd(), file[:-4], 'build_profile.prf')
    file_build_ali = os.path.join(os.getcwd(), file[:-4], 'build_profile.ali')

    sdb = SequenceDB(env)
    sdb.read(seq_database_file='pdb_95.pir', seq_database_format='PIR',
             chains_list='ALL', minmax_db_seq_len=(30, 4000), clean_sequences=True)

    # Write the sequence database in binary form
    sdb.write(seq_database_file=file_bin, seq_database_format='BINARY',
              chains_list='ALL')

    # Now, read in the binary database
    sdb.read(seq_database_file=file_bin, seq_database_format='BINARY',
             chains_list='ALL')

    # Read in the target sequence/alignment
    aln = Alignment(env)
    aln.append(file=file_ali, alignment_format='PIR', align_codes='ALL')

    # Convert the input sequence/alignment into
    # profile format
    prf = aln.to_profile()

    # Scan sequence database to pick up homologous sequences
    prf.build(sdb, matrix_offset=-450, rr_file='${LIB}/blosum62.sim.mat',
              gap_penalties_1d=(-500, -50), n_prof_iterations=1,
              check_profile=False, max_aln_evalue=0.01)

    # Write out the profile in text format
    prf.write(file=file_prf, profile_format='TEXT')

    # Convert the profile back to alignment format
    aln = prf.to_alignment()

    # Write out the alignment file
    aln.write(file=file_build_ali, alignment_format='PIR')


def aligne2D(file):
    file_ali_save = os.path.join(os.getcwd(), file[:-4], file[:-4] + '-3GDK.ali')
    file_pap_save = os.path.join(os.getcwd(), file[:-4], file[:-4] + '-3GDK.pap')
    env = Environ()
    aln = Alignment(env)
    mdl = Model(env, file='3GDK', model_segment=('FIRST:A', 'LAST:A'))
    aln.append_model(mdl, align_codes='3GDK', atom_files='3GDK.pdb')
    aln.append(file=file, align_codes=file[:-4])
    aln.align2d(max_gap_length=50)
    aln.write(file=file_ali_save, alignment_format='PIR')
    aln.write(file=file_pap_save, alignment_format='PAP')


def model_single(file):
    file_ali_save = os.path.join(os.getcwd(), file[:-4], file[:-4] + '-3GDK.ali')
    env = Environ()
    a = AutoModel(env, alnfile=file_ali_save,
                  knowns='3GDK', sequence=file[:-4],
                  assess_methods=(assess.DOPE,
                                  # soap_protein_od.Scorer(),
                                  assess.GA341))
    a.starting_model = 1
    a.ending_model = 5
    a.make()


def create_df(file_ali):
    # create catalog for sequence
    if os.path.exists(file_ali[:-4]):
        None
    else:
        os.makedirs(file_ali[:-4])
    # os.makedirs(file_ali[:-4])

    # Create model
    build_profile(file=file_ali)

    aligne2D(file=file_ali)

    # create file pdb and save stdout to model-single.log
    sys.stdout = open(os.path.join(os.getcwd(), file_ali[:-4], 'model-single.log'), 'w')
    model_single(file=file_ali)
    sys.stdout.close()

    # restore default stdout
    sys.stdout = sys.__stdout__

    #########################
    # choose the best model #
    #########################
    file_log = os.path.join(os.getcwd(), file_ali[:-4], 'model-single.log')
    model_single_log = open(file_log, 'rt')
    lines = model_single_log.readlines()[-6:-1]

    model_log_list = []
    for line in lines:
        model_log_list.append(line.strip().split())

    df_model_log = pd.DataFrame(data=model_log_list, columns=['filename', 'molpdf', 'DOPE', 'GA341'])
    file_pdb = df_model_log[df_model_log['DOPE'] == df_model_log['DOPE'].min()].iloc[0, 0]
    # print(file_pdb)

    ########################################
    # Calculate location aa in active side #
    ########################################
    result_list = []
    for set_aa in important_aa:
        result = calculate_parameters(file=file_pdb, pair_aa=set_aa)
        result_list.append(result)

    data_frame = pd.DataFrame(data=result_list, columns=["file", "MINCEN", "MIN",
                                                         "amino_acid1", "position1",
                                                         "amino_acid2", "position2"])

    ##################################
    # Create list of file in catalog #
    ##################################
    path = os.getcwd()
    files = os.listdir(path)
    # chose files with pattern from list
    pattern = re.compile(file_ali[:-4])
    files = list(filter(pattern.match, files))[1:]

    # move file
    for file in files:
        shutil.move(file, os.path.join(file_ali[:-4], file))

    return data_frame


# important amino acid positions of the URA3 active site
important_aa = [[["A", "GLY", "214"], ["A", "LYS", "156"]],
                [["A", "ASN", "120"], ["A", "LYS", "93"]],
                [["A", "ASN", "120"], ["A", "ASP", "91"]],
                [["A", "ASP", "91"], ["A", "LYS", "93"]],
                [["A", "ASP", "91"], ["A", "LYS", "59"]],
                [["A", "SER", "35"], ["A", "ASP", "37"]]]

log.verbose()
env = Environ()

"""
#################################################
# Calculate for 'URA3_experiment_sequences.fas' #
#################################################
sequences = open('URA3_experiment_sequences.fas', 'r')

lines = sequences.readlines()
sequences_list = []
dataframe_results = pd.DataFrame(columns=["file", "MINCEN", "MIN",
                                          "amino_acid1", "position1",
                                          "amino_acid2", "position2"])
for i in range(0, len(lines) - 1, 2):
    tmp = [">P1" + lines[i].replace(">", ";URA3_").replace(".seq", ""),
           "sequence:" + lines[i].replace(">", ";URA3_").replace(".seq\n", "") + ":: :: ::: 0.00: 0.00\n",
           lines[i + 1].replace("\n", "*\n")]
    # print(tmp)

    tmp_ali = open(tmp[0].replace(">P1;", "").replace("\n", "") + ".ali", 'w')
    tmp_ali.write(tmp[0] + tmp[1] + tmp[2])
    tmp_ali.close()

    file_ali = tmp[0].replace(">P1;", "").replace("\n", "") + ".ali"
    print(file_ali)
    df_sequence = create_df(file_ali=file_ali)
    dataframe_results = dataframe_results.append(df_sequence, ignore_index=True)

dataframe_results.to_csv("result_experimental_sequences.tsv", sep="\t")

print(dataframe_results)
"""


##############################################################
# Calculate for 'ura3_all_missense_possible_sequences.fasta' #
##############################################################
sequences = open('ura3_all_missense_possible_sequences.fasta', 'r')

lines = sequences.readlines()
sequences_list = []
dataframe_results = pd.DataFrame(columns=["file", "MINCEN", "MIN",
                                          "amino_acid1", "position1",
                                          "amino_acid2", "position2"])
for i in range(0, len(lines) - 1, 2):
    tmp = [">P1" + lines[i].replace(">", ";URA3_").replace(".seq", ""),
           "sequence:" + lines[i].replace(">", ";URA3_").replace(".seq\n", "") + ":: :: ::: 0.00: 0.00\n",
           lines[i + 1].replace("\n", "*\n")]
    # print(tmp)

    tmp_ali = open(tmp[0].replace(">P1;", "").replace("\n", "") + ".ali", 'w')
    tmp_ali.write(tmp[0] + tmp[1] + tmp[2])
    tmp_ali.close()

    file_ali = tmp[0].replace(">P1;", "").replace("\n", "") + ".ali"
    print(file_ali)
    df_sequence = create_df(file_ali=file_ali)
    dataframe_results = dataframe_results.append(df_sequence, ignore_index=True)

dataframe_results.to_csv("result_missense_possible_sequences.tsv",  sep="\t")

print(dataframe_results)

