#!/usr/bin/python3
'''
Script to initiate python code for sequence QC data collection, renaming of 
the fasta files and generating output csv file with collated QC data

Input: directory location of the post-assembly pipeline fasta files, csv 
  files and variant text files 
Output: renamed fasta files and collated csv files with the QC data


UKHSA, RVU, Tiina Talts 24.08.2023 V.01

'''


import sys
import os
import subprocess
import glob
import pandas as pd
import pathlib
import shutil
from pathlib import Path
import post_QC_stats as ps
from post_QC_stats import post_qc

def start_QCscript(arg1, arg2):
  
  run_folder0 = arg1.split('\\')
  run_folder1 = str(run_folder0[-1:])
  run_folder2 = run_folder1.replace("[", "")
  run_folder3 = run_folder2.replace("]", "")
  run_folder = run_folder3.replace("'", "")
  
  rack = f'post_qc_'f'{arg2}'
  curdir = os.getcwd()
  path_temp = os.path.join(curdir, 'temp')
  os.makedirs(path_temp, exist_ok=True)
  temp = f'{curdir}/temp/'
  path = os.path.dirname(arg1)
  
  df = pd.ExcelFile(f'{path}\\'f'{run_folder}\\'f'RSV_WGS_tsv_names_'f'{run_folder}.xlsx').parse('Retrieve Data')
  df = df.dropna()
  sample_list = []
  sample_list = df['For_QC:'].tolist()
  
  try:
    os.makedirs(f'{path}\\qc_out_'f'{run_folder}')
  except FileExistsError:
    pass
  
  for sample_name in sample_list:
      ps.post_qc(run_folder, sample_name, arg1)
  
  os.chdir(f'{temp}')
  extension = 'csv'
  all_filenames = [i for i in glob.glob('*.{}'.format(extension))]
  combined_csv = pd.concat([pd.read_csv(f) for f in all_filenames])
  combined_csv.to_csv(f'{rack}.csv', index=False, encoding='utf-8-sig')
  
  files = glob.glob(f'{temp}\\*.temp.csv')
  for f in files:
      try:
          os.remove(f)
      except OSError:
          print("Error while deleting file")
  
  shutil.copy(f'{temp}\\'f'{rack}.csv', f'{path}\\qc_out_'f'{run_folder}')

  file_to_delete = pathlib.Path(f'{temp}\\'f'{rack}.csv')
  file_to_delete.unlink(missing_ok = True)
  
  
  return None



