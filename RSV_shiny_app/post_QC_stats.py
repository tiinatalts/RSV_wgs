#!/usr/bin/python3

'''
Script to initiate python code for sequence QC data collection, renaming of 
the fasta files and generating output csv file with collated QC data

Input: directory location of the post-assembly pipeline fasta files, csv 
  files and variant text files 
Output: renamed fasta files and collated csv files with the QC data


UKHSA, RVU, Tiina Talts 24.08.2023 V.01

'''




import subprocess
import re
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import csv
import numpy
import pandas as pd
import os
from pathlib import Path

def post_qc(run_folder, sample_name, arg1):

    run_folder = run_folder
    sample_name = sample_name
    #path = os.getcwd()
    path = os.path.dirname(arg1)
    
    curdir = os.getcwd()
    path_temp = os.path.join(curdir, 'temp')
    os.makedirs(path_temp, exist_ok=True)
    temp = f'{curdir}/temp/'
    
    
    depth_file = f'{path}\\'f'{run_folder}\\'f'QuasiBAM__ivar_trim__'f'{sample_name}.txt'
    data = pd.read_csv(depth_file, sep='\t')
    av_depth = data['Depth'].mean()

    reads_file = f'{path}\\'f'{run_folder}\\'f'QuasiBAM__ivar_trim__'f'{sample_name}.stat.csv'
    with open(reads_file, 'r') as f:
      lines = f.read().splitlines()
    for line in lines:
      line0 = line.split(',')
      mapped_reads = line0[1]
      filtered_reads = line0[2]
    
    seq_fasta = f'{path}\\'f'{run_folder}\\'f'QuasiBAM__ivar_trim__'f'{sample_name}.fas'
    with open(seq_fasta, 'r') as f:
        sequence = f.read().splitlines()
    t = ''
    d = {}
    termN1 = 0
    seq = sequence[1]
    for i in seq:
        if i == "N":
            termN1 += 1
        if i != "N":
            break
    seq0 = seq[termN1:]
    seqrev = seq0[::-1]
    termN2 = 0
    for i in seqrev:
        if i == "N":
            termN2 += 1
        if i != "N":
            break
    seqrev2 = seqrev[termN2:]
    final_seq = seqrev2[::-1]
    
    pa = re.compile(r'([HRE]{1,2}[0-9]{8,9}_OM_[ab]{1})')
    p2a = re.compile(r'(RSVAB_B[0-9]{6}_[0-9]{2}_[ab]{1}_EFAR_[0-9]{4})')
    p3a = re.compile(r'(RSVAB_B[0-9]{6}_[ab]{1}_EFAR_[0-9]{4})')
    p4a = re.compile(r'(wtr_[ab]{1}_EFAR_[0-9]{4})', re.IGNORECASE)
    
    for s in sequence:
        match = pa.search(s) or p2a.search(s) or p4a.search(s) or p3a.search(s)
        if match:
            u = t.join(match.group(1))
            d[u] = final_seq
    totalN = 0
    totaldegens = 0
    degens = 'RYSWKMBDHV'
    for k, v in d.items():
        record = SeqRecord(Seq(v), id=k, description="")
        for i in v:
            if i == "N":
                totalN += 1
            elif i in degens:
                totaldegens += 1
        try:
            percN = totalN/len(v)*100
        except:
            if len(v) == 0:
                percN = 0
        try:
            percdegens = totaldegens/len(v)*100
        except:
            if len(v) == 0:
                percdegens = 0
        with open(f'{path}\\qc_out_'f'{run_folder}\\'f'{k}.fas', "a") as f:
            SeqIO.write(record, f, "fasta")
                
        header = ['MOLIS', 
                  'sequence length', 
                  'avDepth', 'percN', 
                  'percdegens', 
                  'raw_mapped_reads', 
                  'filtered_mapped_reads']
        data = [[f'{k}', 
                len(v), 
                av_depth, 
                percN, 
                percdegens, 
                mapped_reads, 
                filtered_reads]]
        
        new_file = f'{temp}\\'f'{k}.temp.csv'
        
        with open(new_file, 'a', newline="") as file:
            csvwriter = csv.writer(file)
            csvwriter.writerow(header)
            csvwriter.writerows(data)
        return None
