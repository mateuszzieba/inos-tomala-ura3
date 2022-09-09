##########################
# import needed packages #
##########################
from genericpath import exists
from itertools import chain
from operator import length_hint
from tabnanny import verbose
import mathex
import os
import re
import shutil
import sys

import pandas as pd
from modeller import *
from modeller.automodel import *
from modeller import soap_protein_od
from modeller.scripts import complete_pdb

from Bio import SeqIO

# set dir of project
# os.chdir('/home/mateusz/projects/github/inos-tomala-ura3')
#########################################
# prepare needed files to create models #
#########################################
sequences = open('data/URA3_experiment_sequences.fas', 'r')

# read fas file to dictonary
fas_dict = dict()
with open('data/URA3_experiment_sequences.fas') as f:
    lines = f.readlines()
    for i in range(0, len(lines)):
        s = lines[i].strip()
        if s[0] == '>':
            key = s[1:]
        else:
            fas_dict[key] = s

# create ali files
for i in list(fas_dict.items()):
    ali_file = [">P1;URA3_" + i[0].replace(".seq", ""),
                "sequence:;URA_" + i[0].replace(".seq", "") + ":: :: ::: 0.00: 0.00",
                i[1] + "*"]
    print(ali_file)

    with open('data/ali-files/URA3_' + i[0].replace(".seq", "") +'.ali', 'w') as output_file:
        for month in ali_file:
            output_file.write(month + '\n')



# test modeller command

log.verbose()
env = Environ()
    
# Prepare the input files
pir_file='data/pdb_95.pir' # this file is the same for each model
bin_file='data/pdb_95.bin' # this file is the same for each model
template_pdb='data/3gdk.pdb' # this file is the same for each model
models_dir='data/models' # this dir is the same for each model
ali_file='data/ali-files/URA3_a1_st1.ali' # file as an input variable --> for function/script

# create variable with model name
sequence_name=ali_file.split('/')[-1].replace('.ali', '') 
model_name=sequence_name + '-' + template_pdb.split('/')[-1].replace('.pdb', '')
# create variable with model dir
model_dir=os.path.join(models_dir, model_name)
print(model_dir)
# create a directory if not exist for prepare model
if not os.path.exists(model_dir):
    os.mkdir(model_dir)

#  # created by script at the begining
# prf_file='data/URA3_a1_st1.prf' # created by script 

#################################
# prepare models using modeller #
#################################

# 1. Searching for structures --> it is ok
# Read in the sequence database
sdb = SequenceDB(env)
sdb.read(seq_database_file=pir_file, seq_database_format='PIR',
         chains_list='ALL', minmax_db_seq_len=(30, 4000), clean_sequences=True)

# Write the sequence database in binary form 
sdb.write(seq_database_file=bin_file, seq_database_format='BINARY',
          chains_list='ALL')

# Now, read in the binary database
sdb.read(seq_database_file=bin_file, seq_database_format='BINARY',
         chains_list='ALL')

# Read in the target sequence/aligment
aln = Alignment(env)
aln.append(file=ali_file, alignment_format='PIR', align_codes='ALL')

# Convert the input sequence/alignmetn into profile format
prf = aln.to_profile()

# Scan sequence database to pick up homologous sequences
prf.build(sdb, matrix_offset=-450, rr_file='${LIB}/blosum62.sim.mat',
          gap_penalties_1d=(-500, -50), n_prof_iterations=1,
          check_profile=False, max_aln_evalue=0.01)

# Write out the profile in text format
prf_file=os.path.join(model_dir, sequence_name + '.prf')
prf.write(file=prf_file, profile_format='TEXT')

# Convertthe profile back to alignment format
aln = prf.to_alignment()

# Write out the alignment file
build_profile_ali=os.path.join(model_dir, sequence_name + '.ali')
aln.write(file=build_profile_ali, alignment_format='PIR')

## 3 Aligning with the templatefrom modeller.automodel import *
from modeller import *
env = Environ()
aln = Alignment(env)
mdl = Model(env, file=template_pdb, model_segment=('FIRST:A','LAST:A'))
aln.append_model(mdl, align_codes='3gdk', atom_files=template_pdb)
aln.append(file=ali_file, align_codes=sequence_name)

aln.align2d(max_gap_length=50)

alignment_ali=os.path.join(model_dir, model_name + '.ali')
alignment_pap=os.path.join(model_dir, model_name + '.pap')

aln.write(file=alignment_ali, alignment_format='PIR')
aln.write(file=alignment_pap, alignment_format='PAP')


## 4. Model building

env = Environ()
log_file=os.path.join(model_dir, model_name + '.log')
# save stdout to file
sys.stdout = open(os.path.join(log_file), 'w')
a = AutoModel(env, alnfile=alignment_ali,
              knowns='3gdk',  sequence=sequence_name, 
              assess_methods=(assess.DOPE, 
                            #   soap_protein_od.Scores()
                              assess.GA341))

a.starting_model=1
a.ending_model=5 # this value should be input variable
a.make()
# close stdout file
sys.stdout.close()

# restore default stdout
sys.stdout = sys.__stdout__


# chose files to move to model folder
files=os.listdir(os.getcwd()) 
pattern=re.compile(sequence_name)
files_to_move = list(filter(pattern.match, files))

# move files to model folder
for file in files_to_move:
    if os.path.exists(os.path.join(model_dir, file)):
        os.remove(os.path.join(model_dir, file))
        shutil.move(file, model_dir)
    else:
        shutil.move(file, model_dir)

# chose line with info about model
f = open(log_file, "r").readlines()
pattern=re.compile('^' + sequence_name)
model_lines = list(
    map(
        lambda x: x.replace("\n", ""), list(filter(pattern.match, f))
        )
    )
model_list = list(
    map(
        lambda x: re.sub(r" +", "\t", x).split('\t'), model_lines
        )
    )
model_list

df=pd.DataFrame(
    model_list, 
    columns=['file', 'molpdf', 'DOPE_score', 'GA341_score']
    ).astype(
        {"file": str, "molpdf": float, "DOPE_score": float, "GA341_score": float}
    )

best_model = df[df.DOPE_score == df.DOPE_score.min()].iloc[0,0]


## 5. Model evaluation 


log.verbose()
env = Environ()
env.libs.topology.read(file='$(LIB)/top_heav.lib') # read topology
env.libs.parameters.read(file='$(LIB)/par.lib') # read parameters

# read model file
best_model_dir=os.path.join(model_dir, best_model)
mdl = complete_pdb(env, best_model_dir)

s = Selection(mdl)   # all atom selection
profile_file=os.path.join(model_dir, model_name + '.profile')
s.assess_dope(output='ENERGY_PROFILE NO_REPORT', file=profile_file,
              normalize_profile=True, smoothing_window=15)