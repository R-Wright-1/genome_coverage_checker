import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import pickle
import matplotlib as mpl
import matplotlib.cm as cm
from multiprocessing import Pool
from multiprocessing import freeze_support
import argparse

parser = argparse.ArgumentParser(description='This script is to plot the coverage of taxa after the run_coverage_checker.py script has been run')
parser.add_argument('--sample_name', dest='sample_name', default='all',
                    help='Type the sample name as it appears without any file extension. The fastq file should be unzipped and the kraken .kreport and .kraken.txt files should also be in the same directory')
parser.add_argument('--processors', dest='n_proc', default=1,
                    help="Number of processors to use for all multiprocessing steps (downloading genomes, pulling out reads, running QUAST)")
parser.add_argument('--sample_dir', dest='sample_dir', default=None,
                    help="The directory containing the fastq and kraken files, and that the output will be saved to")
parser.add_argument('--num_tax', dest='num_tax', default=10,
                    help="Look at all prokaryotic taxa with > this number of reads mapped by Kraken")
parser.add_argument('--granularity', dest='granularity', default=1000,
                    help="Look at all prokaryotic taxa with > this number of reads mapped by Kraken")

args = parser.parse_args()
sample_name = args.sample_name
wd = args.sample_dir
granularity = int(args.granularity)
n_processors = int(args.n_proc)
top_x = int(args.num_tax)

# sample_name = 'ani100_cLOW_stFalse_r0'
# wd = '/home/robyn/coverage_checker/parks_ani100_cLOW_stFalse_r0/'
# granularity = 1000
# n_processors = 24
# top_x = 30

all_tax = pd.read_csv(wd+sample_name+'_coverage.csv', index_col=1, header=0).rename(columns={'Unnamed: 0':'Species name'})
nr = 'Number of reads'
if nr not in all_tax.columns:
  nr = 'Number of reads Kraken assigned to this species'
all_tax = all_tax.sort_values(by=[nr], ascending=False)
taxids = list(all_tax.index.values)
if len(taxids) > top_x:
  taxids = taxids[:top_x]

#if '1773' not in taxids:
#  taxids = taxids+['1773']

genomes_not_working = []

def get_coverage(t):
    all_coverage = {}
    tax_name, tax_id = all_tax.loc[t, 'Species name'], str(t)
    if os.path.exists(wd+'pickle_coverage/'+str(tax_id)+'_'+str(granularity)+'_identity.dict'):
      return
    # print(tax_name, tax_id)
    folder = wd+'QUAST/'+sample_name+'_'+tax_id+'/'
    report = pd.read_csv(folder+'report.tsv', index_col=0, header=0, sep='\t')
    ref_chromosomes = []
    ref_chromosome_dict = {}
    chromo_before = 0
    chromosome_lines = []
    try:
      with open(folder+'genome_stats/genome_info.txt') as f:
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
      length = report.loc['Reference length', sample_name+'_'+tax_id]
      alignments = pd.read_csv(folder+'contigs_reports/all_alignments_'+sample_name+'_'+tax_id+'.tsv', index_col=None, header=0, sep='\t')
    except:
      genomes_not_working.append(t)
      return
    
    along_genome, along_genome_idy = {}, {}
    for a in range(1, int(length)+1):
      along_genome[a] = 0
      along_genome_idy[a] = []
    
    genome_coverage = []
    genome_identity = []
    alignments = alignments[alignments['S1'] != 'CONTIG']
    for row in alignments.index:
      if alignments.loc[row, 'IDY'] > 50:
        start, end = int(alignments.loc[row, 'S1']), int(alignments.loc[row, 'E1'])
        genome_identity.append(float(alignments.loc[row, 'IDY']))
        adding = ref_chromosome_dict[alignments.loc[row, 'Reference']]
        s, e = min([start, end])+adding, max([start, end])+adding
        for bp in range(s, e+1):
          along_genome[bp] += 1
          #along_genome_idy[bp].append(idy)
          genome_coverage.append(bp)
    
    k, vals = list(along_genome.keys()), list(along_genome.values())
    groups = [vals[x:x+granularity] for x in range(0, len(vals), granularity)]
    groups_gen = [k[x:x+granularity] for x in range(0, len(k), granularity)]
    means = [sum(group)/len(group) for group in groups]
    gen_means = [sum(group)/len(group) for group in groups_gen]
    means_df = pd.DataFrame(means, index=gen_means).transpose()
    
    all_coverage[tax_id] = [genome_coverage, means_df, chromosome_lines, genome_identity]#, means_df_idy, means_df_idy_error_low, means_df_idy_error_high]
    with open(wd+'pickle_coverage/'+str(tax_id)+'_'+str(granularity)+'_identity.dict', 'wb') as f:
      pickle.dump(all_coverage, f)
    return
  
def run_multiprocessing(func, i, n_processors):
    with Pool(processes=n_processors) as pool:
        return pool.map(func, i)

def main():
    if not os.path.exists(wd+'all_coverage_identity_'+str(granularity)+'.dict'):
      run_multiprocessing(get_coverage, taxids, int(n_processors))

if not os.path.exists(wd+'pickle_coverage'):
  os.system('mkdir '+wd+'pickle_coverage')

if __name__ == "__main__":
    freeze_support()   # required to use multiprocessing
    main()

no_coverage = []
if not os.path.exists(wd+'all_coverage_identity_'+str(granularity)+'.dict'):
  all_coverage = {}
  for t in taxids:
    try:
      with open(wd+'pickle_coverage/'+str(t)+'_'+str(granularity)+'_identity.dict', 'rb') as f:
        cov = pickle.load(f)
      for item in cov:
        all_coverage[item] = cov[item]
    except:
      no_coverage.append(t)
  with open(wd+'all_coverage_identity_'+str(granularity)+'.dict', 'wb') as f:
    pickle.dump(all_coverage, f)
else:
  with open(wd+'all_coverage_identity_'+str(granularity)+'.dict', 'rb') as f:
      all_coverage = pickle.load(f)
      
def plot_cell(ax, colormap, vmin, vmax, num):
  plt.sca(ax)
  colormap = cm.get_cmap(colormap)
  m = cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=vmin, vmax=vmax), cmap=colormap)
  a = ax.bar(0, height=1, width=1, color=m.to_rgba(num), edgecolor='k')
  xl = ax.set_xlim([-0.5, 0.5]), ax.set_ylim([0, 1])
  fc = 'k'
  if num > (vmax/2): fc = 'w'
  tx = ax.text(0, 0.5, str(round(num, 3)), ha='center', va='center', fontsize=8, color=fc)
  xt = plt.xticks([]), plt.yticks([])
  if ax == axes[0][2]: ti = plt.title('Genome fraction (%)', rotation=90, fontweight='bold')
  elif ax == axes[0][3]: ti = plt.title('Actual/expected fraction', rotation=90, fontweight='bold')
  elif ax == axes[0][4]: ti = plt.title('Kraken2 reads mapped', rotation=90, fontweight='bold')
  elif ax == axes[0][5]: ti = plt.title('Reads mapped to genome', rotation=90, fontweight='bold')
  return

fig = plt.figure(figsize=(20,134))
axes = []
for a in range(len(all_coverage)):
  ax1 = plt.subplot2grid((134,20),(a,0), colspan=7)
  ax1_idy = plt.subplot2grid((134,20),(a,8), colspan=3)
  ax2 = plt.subplot2grid((134,20),(a,12))
  ax3 = plt.subplot2grid((134,20),(a,13))
  ax4 = plt.subplot2grid((134,20),(a,14))
  ax5 = plt.subplot2grid((134,20),(a,15))
  axes.append([ax1, ax1_idy, ax2, ax3, ax4, ax5])
  
count = 0
for taxid in all_coverage:
  plt.sca(axes[count][0])
  means_df = all_coverage[taxid][1]
  gen_means = list(means_df.columns)
  means_idy = all_coverage[taxid][3]
  chromosome_lines = all_coverage[taxid][2]
  c = axes[count][0].pcolor(means_df, vmin=0, vmax=1)
  xt = list(plt.xticks()[0][:-1])
  if len(chromosome_lines) < 20:
    for l in chromosome_lines:
      li = axes[count][0].plot([l/granularity, l/granularity], [0, 1], 'w', linewidth=1, alpha=0.5)
  xt = list(plt.xticks()[0][:-1])
  t = plt.xticks(xt, [round(gen_means[int(x)]/1000000, 1) for x in xt])
  t = plt.yticks([])
  # if axes[count][0] == axes[-1][0]: 
  #   xl = plt.xlabel('Genome size (mbp)')
  xl = plt.xlabel('Genome size (mbp)')
  yl = plt.ylabel(all_tax.loc[int(taxid), 'Species name']+' ('+str(taxid)+')', fontweight='bold', rotation=0, ha='right', va='center')
  if count == 0:
    ti = plt.title('Coverage', fontweight='bold')
  plt.sca(axes[count][1])
  box = plt.boxplot(means_idy, positions=[0.5], widths=0.8, vert=False, showfliers=False)
  for item in ['boxes', 'whiskers', 'fliers', 'medians', 'caps']: bx = plt.setp(box[item], color='k')
  yl = plt.ylim([0,1])
  xl = plt.xlim([80,104])
  t = plt.yticks([])
  if axes[count][1] == axes[-1][0]:
    xl = plt.xlabel('Identity (%)')
  if count == 0:
    ti = plt.title('Identity', fontweight='bold')
  gen_frac, act_exp = all_tax.loc[int(taxid), 'Genome fraction (%)'], all_tax.loc[int(taxid), 'Actual/expected fraction']
  reads_krak, reads_quast_mapped, reads_quast_unmapped = all_tax.loc[int(taxid), nr], all_tax.loc[int(taxid), 'QUAST "contigs"'], all_tax.loc[int(taxid), 'QUAST unaligned "contigs"']
  reads_quast = reads_quast_mapped-reads_quast_unmapped
  plot_cell(axes[count][2], 'Reds', 0, 100, gen_frac)
  plot_cell(axes[count][3], 'Oranges', 0, 1, act_exp)
  plot_cell(axes[count][4], 'Greens', 0, 100000, reads_krak)
  plot_cell(axes[count][5], 'Blues', 0, 100000, reads_quast)
  count += 1

plt.subplots_adjust(hspace=2)
plt.savefig(wd+'/'+sample_name+'_'+str(top_x)+'.png', bbox_inches='tight', dpi=300)
