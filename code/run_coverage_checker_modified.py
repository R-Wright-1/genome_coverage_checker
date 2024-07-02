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
from function_file_modified import *

#Add all command line arguments
parser = argparse.ArgumentParser(description='This script is to check which taxa reads have been assigned to by Kraken, pull out these reads, download reference genomes for the taxa, and map the reads to the reference genomes.\n\
                                 It requires the bacterial and archaeal assembly summaries from NCBI, the extract_kraken_reads.py script and QUAST')
parser.add_argument('--sample_name', dest='sample_name', default='all',
                    help='Type the sample name as it appears without any file extension. The fastq file should be unzipped and the kraken .kreport and .kraken.txt files should also be in the same directory')
parser.add_argument('--processors', dest='n_proc', default=1,
                    help="Number of processors to use for all multiprocessing steps (downloading genomes, pulling out reads, running QUAST)")
parser.add_argument('--sample_dir', dest='sample_dir', default=None,
                    help="The directory containing the fastq and kraken files, and that the output will be saved to")
parser.add_argument('--fastq_dir', dest='fastq_dir', default=None,
                    help="The directory containing the fastq files")
parser.add_argument('--kraken_kreport_dir', dest='kraken_kreport_dir', default=None,
                    help="The directory containing the kraken kreport files")
parser.add_argument('--kraken_outraw_dir', dest='kraken_outraw_dir', default=None,
                    help="The directory containing the kraken outraw files")
parser.add_argument('--output_dir', dest='output_dir', default=None,
                    help="The directory to save the output files to")
parser.add_argument('--quast_loc', dest='quast_loc', default=None,
                    help="The location of the executable quast.py script (folder containing QUAST)")
parser.add_argument('--read_lim', dest='read_lim', default=None,
                    help="Look at all prokaryotic taxa with > this number of reads mapped by Kraken")
parser.add_argument('--read_mean', dest='read_mean', default=100,
                    help="Look at all prokaryotic taxa with > this number of mean reads mapped by Kraken within sample groupings")
parser.add_argument('--assembly_folder', dest='assembly_folder', default=None,
                    help="The folder containing the assembly summaries for bacteria and archaea")
parser.add_argument('--sample_metadata', dest='sample_metadata', default=None,
                    help="Location of the sample metadata file. It is expected that this will be a CSV (comma separated) file with a header and two columns. The first column contains sample names and the second contains the sample groupings")
parser.add_argument('--species', dest='species', default=None,
                    help="Location of the species file. An optional file containing a list of taxonomy ID's or species names to include")
parser.add_argument('--project_name', dest='project_name', default=None,
                    help="Name of this project")
parser.add_argument('--rerun', dest='rerun', default=False,
                    help="If this is set to True, it will re-run the extraction of reads from samples regardless of whether the files already exist")
                

#Read in the command line arguments
args = parser.parse_args()
sample_name = args.sample_name
n_proc = int(args.n_proc)
sample_dir = args.sample_dir
fastq_dir = args.fastq_dir
kraken_kreport_dir = args.kraken_kreport_dir
kraken_outraw_dir = args.kraken_outraw_dir
output_dir = args.output_dir
if sample_dir != None:
  fastq_dir = sample_dir
  kraken_kreport_dir = sample_dir
  kraken_outraw_dir = sample_dir
  output_dir = sample_dir
quast_loc = args.quast_loc
read_lim = args.read_lim
read_mean = int(args.read_mean)
assembly_folder = args.assembly_folder
if assembly_folder == None:
  assembly_folder = output_dir
sample_metadata = args.sample_metadata
species = args.species
project_name = args.project_name
rerun = args.rerun
wd = os.getcwd()
# n_proc = 24
# sample_dir = None
# sample_name = 'Miralles'
# fastq_dir = '/home/robyn/kraken_coverage/Miralles/cat_reads/'
# kraken_kreport_dir = '/home/robyn/kraken_coverage/Miralles/kraken_kreport/'
# kraken_outraw_dir = '/home/robyn/kraken_coverage/Miralles/kraken_outraw/'
# output_dir = '/home/robyn/kraken_coverage/coverage_output/'
# assembly_folder = output_dir
# quast_loc = '/home/robyn/tools/quast-5.2.0/'
# read_lim = None
# read_mean = 1000
# sample_metadata = '/home/robyn/kraken_coverage/Miralles/metadata.csv'
# species = '/home/robyn/kraken_coverage/Miralles/species_to_include.txt'
# project_name = 'Miralles'
# rerun = False
if read_lim == None: read_lim = 0
else: read_lim = int(read_lim)

wd, n_proc, sample_dir, sample_name, fastq_dir, kraken_kreport_dir, kraken_outraw_dir, output_dir, assembly_folder, quast_loc, read_lim, read_mean, sample_metadata, species, project_name, rerun, md, samples, taxid_name = run_initial_checks(wd, n_proc, sample_dir, sample_name, fastq_dir, kraken_kreport_dir, kraken_outraw_dir, output_dir, assembly_folder, quast_loc, read_lim, read_mean, sample_metadata, species, project_name, rerun) #run all of the initial checks to ensure that all of the folders and files exist before starting to try and run anything

#read in all kraken kreports and get the species with >= read_lim
kreports, kreports_to_combine, taxids_kreports = {}, [], {}
for sample_name in samples:
  krep = kraken_kreport_dir+sample_name+'.kreport'
  krep2 = kraken_kreport_dir+sample_name+'_0.0.kreport'
  if not os.path.exists(krep): krep = krep2
  kreport = pd.read_csv(krep, header=None, sep='\t')
  kreport = kreport[kreport[3] == 'S'] #keep only the species level classifications
  kreport.columns = ['Proportion of reads', 'Reads assigned to this species', 'Reads assigned to below this species', 'Rank', 'NCBI taxid', 'NCBI species name'] #rename the columns
  kreport = kreport.drop('Rank', axis=1) #drop the rank column
  kreports[sample_name] = kreport
  for row in kreport.index:
    taxids_kreports[kreport.loc[row, 'NCBI taxid']] = kreport.loc[row, 'NCBI species name'].lstrip()
  kreport_red = kreport.copy(deep=True).set_index('NCBI taxid').loc[:, ['Reads assigned to this species']]
  kreports_to_combine.append(kreport_red.rename(columns={'Reads assigned to this species':sample_name}))

#if sample_metadata != None:
kreports_to_combine = pd.concat(kreports_to_combine).fillna(value=0)
kreports_to_combine = kreports_to_combine.groupby(by=kreports_to_combine.index, axis=0).sum()
md_groups = list(set(list(md.iloc[:, 0].values)))
group_samples = {project_name:[]}
taxid_keeping = {}
for group in md_groups:
  samples_group = []
  for r in range(len(md.index.values)):
    if md.iloc[r, 0] == group:
      samples_group.append(md.index.values[r])
      group_samples[project_name].append(md.index.values[r])
  krep_group = kreports_to_combine.copy(deep=True).loc[:, samples_group]
  krep_group = krep_group[krep_group.max(axis=1) > read_lim]
  if read_mean != None:
    krep_group['Mean'] = krep_group.mean(axis=1)
    krep_group = krep_group[krep_group['Mean'] > read_mean]
    krep_group = krep_group.drop('Mean', axis=1)
  for row in krep_group.index.values:
    taxid_keeping[row] = taxids_kreports[row]
  group_samples[group] = samples_group
  
kreports_to_combine.to_csv(output_dir+project_name+'_combined_kreport.csv')

taxid_name.update(taxid_keeping)
download_list, taxids_using, genomes_got, no_genome = get_genomes(taxid_name, assembly_folder, output_dir, taxid_name, representative_only=False)

with open(output_dir+'taxid_name.dict', 'wb') as f:
    pickle.dump(taxid_name, f)

taxids_string, taxids_string_quast, taxids_list_string, all_taxid_to_combine = [], [], [], {}
for sample_name in samples:
  taxid_in_sample = ''
  for r in range(len(md.index.values)):
    if md.index.values[r] == sample_name:
      group = md.iloc[r,0]
  for t in taxids_using:
    if t not in kreports_to_combine.index.values: continue
    if kreports_to_combine.loc[t, sample_name] == 0: continue
    if project_name+'_'+str(t) not in all_taxid_to_combine:
      all_taxid_to_combine[project_name+'_'+str(t)] = []
    if group+'_'+str(t) not in all_taxid_to_combine:
      all_taxid_to_combine[group+'_'+str(t)] = []
    taxids_string.append(str(t)+','+sample_name+','+fastq_dir+','+kraken_outraw_dir+','+kraken_kreport_dir+','+output_dir+','+krep+','+str(rerun))
    taxids_string_quast.append(str(t)+','+sample_name+','+fastq_dir+','+kraken_outraw_dir+','+kraken_kreport_dir+','+output_dir+','+quast_loc+','+str(rerun))
    taxid_in_sample += str(t)+' '
    all_taxid_to_combine[project_name+'_'+str(t)].append(sample_name+'_'+str(t))
    all_taxid_to_combine[group+'_'+str(t)].append(sample_name+'_'+str(t))
  taxids_list_string.append(taxid_in_sample+','+sample_name+','+fastq_dir+','+kraken_outraw_dir+','+kraken_kreport_dir+','+output_dir+','+krep+','+str(rerun))
  
#For all of the genomes that we will download, add them to a dictionary
downloaded_genomes = {}
for g in range(len(genomes_got)):
  downloaded_genomes[str(taxids_using[g])] = genomes_got[g]
  
with open(output_dir+'downloaded_genomes.dict', 'wb') as f:
    pickle.dump(downloaded_genomes, f)
    
files_combining = []
for fn in all_taxid_to_combine:
  files_combining.append(fn+','+output_dir+','+' '.join(all_taxid_to_combine[fn]))
  taxids_string.append(fn.split('_')[-1]+','+fn.split('_')[0]+','+fastq_dir+','+kraken_outraw_dir+','+kraken_kreport_dir+','+output_dir+','+krep+','+str(rerun))
  taxids_string_quast.append(fn.split('_')[-1]+','+fn.split('_')[0]+','+fastq_dir+','+kraken_outraw_dir+','+kraken_kreport_dir+','+output_dir+','+quast_loc+','+str(rerun))

def main():
    #Start genome download with multiprocessing
    print('Starting genome download')
    run_multiprocessing(download_genome, download_list, n_proc)
    
    #Check that all genomes downloaded properly
    print('Checking that genomes were downloaded')
    genome_check(download_list, project_name)
    
    #Start getting the reads for each taxonomy ID using multiprocessing
    print('Starting get reads')
    run_multiprocessing(run_get_reads, taxids_list_string, n_proc)
    
    #combine the files for grouped taxids
    print('Starting combining files')
    run_multiprocessing(combine_files, files_combining, n_proc)
    
    #Check that all files were created
    print('Checking all files were created')
    file_check(taxids_string, project_name)
    
    #Run QUAST for all of the taxonomy ID's using multiprocessing
    print('Starting QUAST runs')
    run_multiprocessing(run_quast, taxids_string_quast, n_proc)
    
    #change back from the QUAST directory where we ran it from
    os.chdir(wd)
    
    #Check whether QUAST output files were created for each taxonomy ID
    quast_check(taxids_string, project_name)
    
    #Make bowtie2 databases and run bowtie2 for all of the taxonomy ID's using multiprocessing
    print('Making bowtie2 databases')
    run_multiprocessing(make_bowtie2_db, taxids_string, n_proc)
    print('Starting bowtie2 runs')
    run_multiprocessing(run_bowtie2, taxids_string, n_proc)
    
    #Getting the coverage and mapping of reads across the genome
    print('Starting getting genome coverage')
    run_multiprocessing(get_coverage, taxids_string, n_proc)
    return
    

if __name__ == "__main__":
    freeze_support()   # required to use multiprocessing
    main()

all_taxids_float = [t for t in taxid_name]
for tax in all_taxids_float:
  taxid_name[str(tax)] = taxid_name[tax]
  
group_samples = {project_name:[]}
group_sample_is_in = {}
for group in md_groups:
  samples_group = []
  for r in range(len(md.index.values)):
    if md.iloc[r, 0] == group:
      samples_group.append(md.index.values[r])
      group_samples[project_name].append(md.index.values[r])
      group_sample_is_in[md.index.values[r]] = group
  group_samples[group] = samples_group

columns = ['NCBI taxonomy ID', 'NCBI species name', 'Sample', 'Sample group', 'Reads Kraken assigned', 'Reads Kraken assigned below species', 'Reference length (bp)', 'Reference GC (%)', 'QUAST contigs', 'QUAST unaligned contigs', 'GC (%)', 'Genome fraction (%)', 'Duplication ratio', 'Total length (number of reads x read length; bp)', 'Number of bases/reference genome length (expected genome fraction; %)', 'Actual/expected fraction', 'Expected/actual fraction', 'Bowtie2 mapped reads', 'Genome median identity', 'Genome proportion with mapped reads', 'Proportion of reads mapped in distinct genomic areas (round 500)', 'Proportion of reads mapped in distinct genomic areas (round 1000)', 'Proportion of reads mapped in distinct genomic areas (round 5000)']
column = ''
for v in range(len(columns)):
  if v != len(columns)-1:
    column += columns[v]+'\t'
  else:
    column += columns[v]+'\n'

try:
  for sample in kreports:
    kreports[sample] = kreports[sample].set_index('NCBI taxid')
except:
  do_nothing = True

for s in range(len(taxids_string)):
  if s == 0:
    with open(output_dir+project_name+"_output.txt", "w") as f:
      w = f.write(column)
  #if s > 100: break
  string = taxids_string[s].split(',')
  #string = random_list[s].split(',')
  tax, sample = string[0], string[1]
  if sample == 'Natural': sample = 'Natural_soil'
  this_tax = [tax, taxid_name[tax], sample]
  if sample in group_sample_is_in:
    this_tax.append(group_sample_is_in[sample])
  else:
    this_tax.append(sample)
  if sample in kreports: 
    snames = [sample]
  else:
    snames = [sn for sn in group_samples[sample]]
  this_sp, below_sp = [], []
  for sn in snames:
    try:
      this_sp.append(kreports[sn].loc[int(tax), 'Reads assigned to this species'])
      below_sp.append(kreports[sn].loc[int(tax), 'Reads assigned to below this species'])
    except:
      this_sp.append(0)
      below_sp.append(0)
  this_tax.append(sum(this_sp))
  this_tax.append(sum(below_sp))
  #get quast output
  try:
    folder = output_dir+'/QUAST/'+sample+'_'+tax+'/'
    report = pd.read_csv(folder+'report.tsv', index_col=0, header=0, sep='\t')
    ref_len = report.loc['Reference length', sample+'_'+tax]
    ref_gc = report.loc['Reference GC (%)', sample+'_'+tax]
    quast_contig = nreads = report.loc['# contigs (>= 0 bp)', sample+'_'+tax]
    try:
      quast_unal_contig = float(report.loc['# unaligned contigs', sample+'_'+tax].split(' ')[0])
    except:
      quast_unal_contig = float(report.loc['# unaligned contigs', sample+'_'+tax])
    quast_gc = report.loc['GC (%)', sample+'_'+tax]
    try:
      genome_frac = report.loc['Genome fraction (%)', sample+'_'+tax]
      dup_ratio = report.loc['Duplication ratio', sample+'_'+tax]
    except:
      genome_frac = 0
      dup_ratio = 'NA'
    total_len = report.loc['Total length (>= 0 bp)', sample+'_'+tax]
    exp_frac = (int(total_len)/int(ref_len))*100
    if exp_frac == 0:
      act_exp_frac, exp_act_frac = 0, 0
    else:
      act_exp_frac = float(genome_frac)/exp_frac
      exp_act_frac = exp_frac/float(genome_frac)
  except:
    ref_len, ref_gc, quast_contig, quast_unal_contig, quast_gc, genome_frac, dup_ratio, total_len, exp_frac, act_exp_frac, exp_act_frac = 'NA', 'NA', 'NA', 'NA', 'NA', 0, 'NA', 0, 'NA', 'NA', 'NA'
  try:
    bowtie2_mapped = 0
    for record in SeqIO.parse(output_dir+'/bowtie2_mapped/'+sample+'_'+tax+'.fasta', "fasta"):
      bowtie2_mapped += 1
  except:
    bowtie2_mapped = 0
  for val in [ref_len, ref_gc, quast_contig, quast_unal_contig, quast_gc, genome_frac, dup_ratio, total_len, exp_frac, act_exp_frac, exp_act_frac, bowtie2_mapped]:
    this_tax.append(val)
  try:
    with open(output_dir+'pickle_coverage/'+str(tax)+'_'+sample+'.dict', 'rb') as f:
      coverage_dict = pickle.load(f)
    genome_coverage, genome_identity, all_starts, length = coverage_dict[sample+'_'+tax]
    this_tax.append(np.median(genome_identity))
    this_tax.append(len(genome_coverage)/(int(length)/500))
    unique_500 = list(set([round((start+50)/500) * 500 for start in all_starts]))
    unique_1000 = list(set([round((start+50)/1000) * 1000 for start in all_starts]))
    unique_5000 = list(set([round((start+50)/5000) * 5000 for start in all_starts]))
    this_tax.append(len(unique_500)/len(all_starts))
    this_tax.append(len(unique_1000)/len(all_starts))
    this_tax.append(len(unique_5000)/len(all_starts))
  except:
    for na in range(5):
      this_tax.append('NA')
  row = ''
  for v in range(len(this_tax)):
    if v != len(this_tax)-1:
      row += str(this_tax[v])+'\t'
    else:
      row += str(this_tax[v])+'\n'
  with open(output_dir+project_name+"_output.txt", "a") as f:
      w = f.write(row)
