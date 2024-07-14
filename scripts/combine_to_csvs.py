import os
import pandas as pd
import numpy as np
import re

output_files_dir = 'output'
genomes_dir = 'genomes'
genomes_fasta_dir = 'fasta'
sample_info_dir = 'info'
info_file_name = 'all_samples.csv'
out_dir = 'res_files'

run_path = './'
output_file_path = os.path.join(run_path, output_files_dir)
genome_path = os.path.join(run_path, genomes_dir)
genomes_fasta_path = os.path.join(genome_path , genomes_fasta_dir)
info_path = os.path.join(run_path, sample_info_dir)

genomes_fasta_files = os.listdir(genomes_fasta_path)
info_file = pd.read_csv(os.path.join(info_path,info_file_name))
info_file = info_file.astype(str)
out_path = os.path.join(run_path, out_dir)
# os.mkdir(out_path)

output_files = os.listdir(output_file_path)

for file in output_files:
    if file not in output_files:
        continue
    file_info = file.split('.')[0]
    tf, exp_prom, lib, MNase, target, run, _, _, _ = file_info.split('_')
    subset = info_file.query(
        "Strain == @tf & Exp_Promoter == @exp_prom & Library_name == @lib & Mnase_on_pep == @MNase &"
        " Target_promoter == @target & Run == @run")
    subset_files = []
    for row_i in subset.index:
        curr_row = subset.loc[row_i]
        composed_name = '_'.join(
            [tf, exp_prom, lib, MNase, target, run, curr_row['Trans_repeat'], curr_row['Time_point'],
             curr_row['Time_point_repeat']]) + '.txt'
        subset_files.append(composed_name)

    lib_genome_name = lib + '_genome.fasta'
    lib_genome = pd.read_csv(os.path.join(genomes_fasta_path, lib_genome_name), 
                             header=None)
    tmp_lib_pep_list = lib_genome.iloc[np.where(lib_genome[0].str.contains('>'))[0]][0]
    lib_pep_list = [pep.split('>')[1] for pep in tmp_lib_pep_list]

    curr_lib_data = pd.DataFrame(index=lib_pep_list)
    for file in subset_files:
        curr_file_data = pd.read_csv(os.path.join(output_file_path, file), sep=r'\s+', header=None,
                                     names=['pep', 'read_num'])
        _, _, _, _, _, _, bio_rep, tp, tp_rep = file.split('.')[0].split('_')
        curr_col_name = '_'.join([bio_rep, tp, tp_rep])

        read_counts = []
        for pep in lib_pep_list:
            if pep in curr_file_data['pep'].values:
                read_counts.append(curr_file_data.loc[curr_file_data['pep'] == pep]['read_num'].values[0])
            else:
                read_counts.append(0)
        curr_lib_data[curr_col_name] = read_counts

    curr_lib_data_ordered = curr_lib_data.sort_index(axis=1)
    csv_name = '_'.join([tf, exp_prom, lib, MNase, target, run])
    curr_lib_data_ordered.to_csv(os.path.join(out_path, csv_name + '.csv'))

    output_files = np.setdiff1d(output_files, subset_files)