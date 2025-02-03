import os
from multiprocessing import Pool, freeze_support
import pandas as pd
import argparse
import numpy as np
import pickle

parser = argparse.ArgumentParser(description='This script is to get the coverage of reads across multiple genomes at once.')
parser.add_argument('--folders', dest='folders',
                    help="The list of QUAST folders to be run in a text file. This should be the file name.")
parser.add_argument('--processors', dest='n_proc', default=1,
                    help="Number of processors to use")

args = parser.parse_args()
folders = args.folders
n_proc = int(args.n_proc)

folder_list = []
for row in open(folders, 'r'):
  folder_list.append(row.replace('\n', ''))

def get_genome_coverage(folder):
  fn = folder.split('/')[-1]
  all_starts, all_ends = [], []
  #all_coverage = {}
  try:
    ref_chromosomes, ref_chromosome_dict, chromo_before, chromosome_lines = [], {}, 0, []
    with open(folder+'/genome_stats/genome_info.txt') as f:
      for row in f:
        if 'total length' in row:
          ref_chromosomes.append(row)
    for r in range(len(ref_chromosomes)):
      chromo = ref_chromosomes[r]
      chromo_name = chromo.split(' ')[0].replace('\t', '')
      chromo_length = chromo.split(':')[1].split('bp')[0].replace(' ', '')
      ref_chromosomes[r] = [chromo_name, chromo_length]
      ref_chromosome_dict[chromo_name] = chromo_before
      chromo_before += int(chromo_length)
      chromosome_lines.append(chromo_before)
    didnt_work = False
  except:
    didnt_work = True
  alignments = pd.read_csv(folder+'/contigs_reports/all_alignments_'+fn+'.tsv', index_col=None, header=0, sep='\t')
  alignments = alignments[alignments['S1'] != 'CONTIG']
  genome_identity, genome_coverage = [], []
  for row in alignments.index.values:
    try:
      start, end = int(alignments.loc[row, 'S1']), int(alignments.loc[row, 'E1'])
    except:
      return
    genome_identity.append(str(alignments.loc[row, 'IDY']))
    adding = ref_chromosome_dict[alignments.loc[row, 'Reference']]
    s, e = min([start, end])+adding, max([start, end])+adding
    middle = np.mean([s, e])
    all_starts.append(str(s))
    all_ends.append(str(e))
    rounded = int(round(middle / 500.0) * 500.0)
    genome_coverage.append(str(rounded))
  genome_coverage_set = list(set(genome_coverage))
  #all_coverage[fn] = [genome_coverage, genome_identity, all_starts, all_ends]
  write_string = 'all_starting_points: '+','.join(all_starts)+'\n'
  write_string += 'all_end_points: '+','.join(all_ends)+'\n'
  write_string += 'genome_coverage_rounded_mid_points: '+','.join(genome_coverage)+'\n'
  write_string += 'genome_identity: '+','.join(genome_identity)
  with open(folder.replace('QUAST', 'coverage')+'.txt', 'w') as f:
    w = f.write(write_string)
  return


def run_multiprocessing(func, i, n_processors):
    with Pool(processes=n_processors) as pool:
        return pool.map(func, i)
      
def main():
  run_multiprocessing(get_genome_coverage, folder_list, n_proc)

if __name__ == "__main__":
    freeze_support()   # required to use multiprocessing
    main()
