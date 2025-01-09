import os
import pandas as pd
import sys
from multiprocessing import Pool
from multiprocessing import freeze_support
import numpy as np
import argparse
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import pickle
from datetime import datetime

#### add description of script
parser = argparse.ArgumentParser(description='')

parser.add_argument('--run_path', dest='run_path', default=None,
                    help='The path to save the output of this run to. Default: Current date and time in the current directory')
parser.add_argument('--overwrite', dest='overwrite', default=False,
                    help='Whether to overwrite a previous run with the same name. Options: True/False. Default: False')
parser.add_argument('--rerun', dest='rerun', default=True,
                    help='Whether to restart a previous run. Default is True to avoid accidentally removing output of a previous run. Default: True')
parser.add_argument('--processors', dest='processors', default=1,
                    help='Number of processors to use for all multiprocessing steps. Default: 1')
parser.add_argument('--genome_folder', dest='genome_folder', default='NCBI_genomes',
                    help="Folder to download assembly summaries and genomes to. The default is to create a folder called NCBI_genomes in the current directory if it doesn't already exist. Default: NCBI_genomes")


args = parser.parse_args()

def add_to_logfile(run_path, string):
  with open(run_path+"/log_file.txt", 'a') as f:
    f.write(string+'\n')
  return

def get_checkpoint(run_path):
  checkpoint = ''
  for row in open(run_path+"/checkpoint.txt", 'r'):
    checkpoint += row.replace('\n', '')
  return checkpoint

def initiate(args):
  if args.run_path == None:
    run_path = 'coverage_checker_'+str(datetime.now()).replace(' ', '_').split('.')[0].replace(':', '.')
  else:
    run_path = args.run_path
  if os.path.exists(run_path):
    #### need to come back and check what else needs to be added here for re-starting previous run
    if not args.overwrite:
      sys.exit("This run_path already exists. If you intented to overwrite it, please set overwrite to True.")
    else:
      if not args.rerun:
        os.system('rm -r '+run_path)
      else:
        checkpoint = get_checkpoint(run_path)
        add_to_logfile(run_path, "Restarting run from checkpoint "+checkpoint)
  else:
    os.mkdir(run_path)
    options_string = 'run_path = '+run_path+'\n'
    for arg in vars(args):
      if str(arg) == 'run_path': continue
      options_string += str(arg)+' = '+str(getattr(args, arg))+'\n'
    with open(run_path+'/run_options.txt', 'w') as f:
      f.write(options_string)
    checkpoint = 0
    with open(run_path+'/checkpoint.txt', 'w') as f:
      f.write('0')
    add_to_logfile(run_path, "Run initiated with run_path "+run_path)
    add_to_logfile(run_path, "Options file and checkpoint file have been made. Current time is: "+str(datetime.now()))
  return
   

initiate(args)
