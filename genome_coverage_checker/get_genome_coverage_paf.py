import os
from multiprocessing import Pool, freeze_support
import pandas as pd
import argparse
import numpy as np
import pickle

parser = argparse.ArgumentParser(description='This script is to get the coverage of reads across multiple genomes at once.')
parser.add_argument('--files', dest='files',
                    help="The list of paf files and genome info to be run in a text file. This should be the file name.")
parser.add_argument('--processors', dest='n_proc', default=1,
                    help="Number of processors to use")
parser.add_argument('--mapq_threshold', dest='mapq_threshold', default=None,
                    help="MAPQ quality threshold to use")
parser.add_argument('--identity_threshold', dest='identity_threshold', default=None,
                    help="Minimum nucleotide identity % threshold to use")

args = parser.parse_args()
files = args.files
n_proc = int(args.n_proc)
mapq_threshold = args.mapq_threshold
identity_threshold = args.identity_threshold
if mapq_threshold in [None, 'None']:
  mapq_threshold = 0
else:
  mapq_threshold = float(mapq_threshold)
if identity_threshold in [None, 'None']:
  identity_threshold = 0
else:
  identity_threshold = float(identity_threshold)

file_list = []
for row in open(files, 'r'):
  file_list.append(row.replace('\n', ''))

def get_genome_coverage(file):
  files_using = file.split(' ')
  paf_file, genome_info = files_using[0], files_using[1]
  with open(genome_info, 'rb') as f:
    contig_lengths = pickle.load(f)
  all_identity, all_starts, all_ends, all_MAPQ = [], [], [], []
  for row in open(paf_file, 'r'):
    line = row.replace('\n', '').split('\t')
    ref_contig, start, end, match_bases, align_len, MAPQ = line[5], int(line[7]), int(line[8]), int(line[9]), int(line[10]),  int(line[11])
    if MAPQ <= mapq_threshold: continue
    start += contig_lengths[ref_contig+'_start']
    end += contig_lengths[ref_contig+'_start']
    identity = (match_bases/align_len)*100
    if identity < identity_threshold: continue
    all_starts.append(str(start)), all_ends.append(str(end)), all_MAPQ.append(str(MAPQ))
    all_identity.append(str(identity))
  write_string = 'all_starting_points: '+','.join(all_starts)+'\n'
  write_string += 'all_end_points: '+','.join(all_ends)+'\n'
  write_string += 'genome_identity: '+','.join(all_identity)+'\n'
  write_string += 'genome_MAPQ: '+','.join(all_MAPQ)+'\n'
  write_string += 'reads_mapped: '+str(len(all_identity))+'\n'
  positions_covered = set()
  for s, e in zip(all_starts, all_ends):
    for g in range(int(s), int(e)+1):
      positions_covered.add(g)
  genome_fraction = (len(positions_covered)/contig_lengths['full_length'])*100
  write_string += 'genome_fraction: '+str(genome_fraction)
  fn = paf_file.split('/')[-1]
  if 'bowtie2' in paf_file:
    with open(paf_file.replace('bowtie2_mapped', 'coverage').replace(fn, 'bowtie2_'+fn).replace('.paf', '.txt'), 'w') as f:
      w = f.write(write_string)
  if 'minimap2' in paf_file:
    with open(paf_file.replace('minimap2_mapped', 'coverage').replace(fn, 'minimap2_'+fn).replace('.paf', '.txt'), 'w') as f:
      w = f.write(write_string)
  return


def run_multiprocessing(func, i, n_processors):
    with Pool(processes=n_processors) as pool:
        return pool.map(func, i)
      
def main():
  run_multiprocessing(get_genome_coverage, file_list, n_proc)

if __name__ == "__main__":
    freeze_support()   # required to use multiprocessing
    main()
