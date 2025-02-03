import os
from multiprocessing import Pool, freeze_support
import pandas as pd
import argparse
import numpy as np
import pickle
from Bio import SeqIO

parser = argparse.ArgumentParser(description="This script is to remove duplicate reads from fastq files. It is so that when coverage checker is run multiple times, we don't end up with duplicates of everything in the output files, which messes up the QUAST metrics")
parser.add_argument('--files', dest='files',
                    help="The list of fastq files to have duplicate reads removed from.")
parser.add_argument('--processors', dest='n_proc', default=1,
                    help="Number of processors to use for all multiprocessing steps")

args = parser.parse_args()
files = args.files
n_proc = int(args.n_proc)

file_list = []
for row in open(files, 'r'):
  file_list.append(row.replace('\n', ''))

def remove_duplicate_reads(fn):
  not_duplicated = []
  not_duplicated_names = []
  count = 0
  for record in SeqIO.parse(fn, "fastq"):
    if record.id not in not_duplicated_names:
      not_duplicated.append(record)
      not_duplicated_names.append(record.id)
    count += 1
  SeqIO.write(not_duplicated, fn, "fastq")
  return


def run_multiprocessing(func, i, n_processors):
    with Pool(processes=n_processors) as pool:
        return pool.map(func, i)
      
def main():
  run_multiprocessing(remove_duplicate_reads, file_list, n_proc)

if __name__ == "__main__":
    freeze_support()   # required to use multiprocessing
    main()
