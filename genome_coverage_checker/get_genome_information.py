import os
from multiprocessing import Pool, freeze_support
#import pandas as pd
import argparse
#import numpy as np
import pickle
from Bio import SeqIO

parser = argparse.ArgumentParser(description="This script is to get information about genomes.")
parser.add_argument('--files', dest='files',
                    help="The list of genome fasta files to get information about.")
parser.add_argument('--processors', dest='n_proc', default=1,
                    help="Number of processors to use for all multiprocessing steps")

args = parser.parse_args()
files = args.files
n_proc = int(args.n_proc)

file_list = []
for row in open(files, 'r'):
  file_list.append(row.replace('\n', ''))

def get_genome_information(fn):
  contig_lengths = {}
  for record in SeqIO.parse(fn, "fasta"):
    contig_lengths[record.id] = len(record.seq)
  start = 0
  contig_starts = {}
  for contig in contig_lengths:
    contig_starts[contig+'_start'] = start
    start += contig_lengths[contig]
  full_length = 0
  for contig in contig_lengths:
    full_length += contig_lengths[contig]
  for contig in contig_starts:
    contig_lengths[contig] = contig_starts[contig]
  contig_lengths['full_length'] = full_length
  with open(fn.replace('.fna', '_info.dict'), 'wb') as f:
    pickle.dump(contig_lengths, f)
  return

def run_multiprocessing(func, i, n_processors):
    with Pool(processes=n_processors) as pool:
        return pool.map(func, i)
      
def main():
  run_multiprocessing(get_genome_information, file_list, n_proc)

if __name__ == "__main__":
    freeze_support()   # required to use multiprocessing
    main()
