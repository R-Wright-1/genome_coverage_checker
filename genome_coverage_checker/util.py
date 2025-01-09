#!/usr/bin/env python

import os
import pandas as pd
import sys
#from multiprocessing import Pool, freeze_support
#import argparse
#from Bio.SeqRecord import SeqRecord
#from Bio import SeqIO
#import pickle
#import numpy as np
from os.path import abspath, dirname
#from genome_coverage_checker.download_genomes_ import *

def get_checkpoint(output_dir):
  cp = ''
  for row in open(output_dir+'checkpoint.txt', 'r'):
    cp += row
  return cp

def update_checkpoint(output_dir, cp):
  with open(output_dir+'checkpoint.txt', 'w') as f:
    c = f.write(str(cp))
  return
  

def run_initial_checks(wd, n_proc, sample_dir, sample_name, fastq_dir, kraken_kreport_dir, kraken_outraw_dir, output_dir, assembly_folder, read_lim, read_mean, sample_metadata, species, project_name, rerun):
  # #check for location of QUAST
  # if quast_loc == None:
  #   sys.exit('You must supply the location of the quast.py script! Quitting the run now without doing anything.')
  
  #check all folders exist
  for folder in [fastq_dir, kraken_kreport_dir, kraken_outraw_dir, output_dir, assembly_folder]:
    if not os.path.exists(folder):
      if folder == output_dir:
        print("output_dir doesn't already exist. Making it now.")
        os.system('mkdir '+output_dir)
      else:
        sys.exit("This path doesn't exist: "+folder)

  fastq_dir, kraken_kreport_dir, kraken_outraw_dir, output_dir, assembly_folder = fastq_dir+'/', kraken_kreport_dir+'/', kraken_outraw_dir+'/', output_dir+'/', assembly_folder+'/'
  
  #directories to make if they don't exist
  all_dirs = [output_dir+'/genomes/', output_dir+'/reads_mapped/', output_dir+'/QUAST/', output_dir+'/bowtie2_db/', output_dir+'/bowtie2_mapped/', output_dir+'/pickle_coverage/']
  for direc in all_dirs:
    if not os.path.exists(direc):
      md = os.system('mkdir '+direc)
  
  #check whether at least one of sample_metadata or sample_name exist
  if sample_name == None and sample_metadata == None:
    sys.exit("You need to supply one of sample_name or sample_metadata")
  elif sample_metadata != None:
    if not os.path.exists(sample_metadata):
      sys.exit("This path doesn't exist: "+sample_metadata)
  
  if sample_metadata != None:
    md = pd.read_csv(sample_metadata, header=0, index_col=0)
    samples = list(md.index.values)
  else:
    samples = [sample_name]
    
  #check all input files exist
  for sample_name in samples:
    fq = fastq_dir+sample_name+'.fastq'
    krep = kraken_kreport_dir+sample_name+'.kreport'
    krep2 = kraken_kreport_dir+sample_name+'_0.0.kreport'
    kraw = kraken_outraw_dir+sample_name+'.kraken'
    for f in [fq, [krep, krep2], kraw]:
      if isinstance(f, str):
        if not os.path.exists(f):
          sys.exit("This file doesn't exist: "+f)
      else:
        if not os.path.exists(f[0]) and not os.path.exists(f[1]):
          sys.exit("Neither of these files exists: "+f[0]+", "+f[1])
  
  #check whether species list exists, if it does, get the taxids
  if species != None:
    if not os.path.exists(species):
      print("This file doesn't exist: "+species)
      sys.exit()
    sp_list = []
    for row in open(species, 'r'):
      sp_list.append(row.replace('\n', '').replace('\t', '').split(';s__')[-1])
    taxid_name, dont_have = get_taxid_for_names(sp_list, assembly_folder)
    with open(output_dir+project_name+'_species_couldnt_get.txt', 'w') as f:
      for sp in dont_have:
        w = f.write(sp+'\n')
  else:
    taxid_name = {}
  return wd, n_proc, sample_dir, sample_name, fastq_dir, kraken_kreport_dir, kraken_outraw_dir, output_dir, assembly_folder, read_lim, read_mean, sample_metadata, species, project_name, rerun, md, samples, taxid_name

def get_kreports(samples, kraken_kreport_dir, output_dir, md, read_lim, read_mean, project_name):
  kreports, taxids_kreports = [], {}
  for sample_name in samples:
    krep = kraken_kreport_dir+sample_name+'.kreport'
    if not os.path.exists(krep): krep = kraken_kreport_dir+sample_name+'_0.0.kreport'
    kreport = pd.read_csv(krep, header=None, sep='\t')
    kreport = kreport[kreport[3] == 'S'] #keep only the species level classifications
    kreport.columns = ['Proportion of reads', 'Reads assigned to this species', 'Reads assigned to below this species', 'Rank', 'NCBI taxid', 'NCBI species name'] #rename the columns
    kreport = kreport.drop('Rank', axis=1) #drop the rank column
    for row in kreport.index:
      taxids_kreports[kreport.loc[row, 'NCBI taxid']] = kreport.loc[row, 'NCBI species name'].lstrip()
    kreport = kreport.set_index('NCBI taxid').loc[:, ['Reads assigned to this species']]
    kreports.append(kreport.rename(columns={'Reads assigned to this species':sample_name}))
    
  kreports = pd.concat(kreports).fillna(value=0)
  kreports = kreports.groupby(by=kreports.index).sum()
  md_groups = list(set(list(md.iloc[:, 0].values)))
  group_samples = {project_name:[]}
  taxid_keeping = {}
  for group in md_groups:
    samples_group = []
    for r in range(len(md.index.values)):
      if md.iloc[r, 0] == group:
        samples_group.append(md.index.values[r])
        group_samples[project_name].append(md.index.values[r])
    krep_group = kreports.copy(deep=True).loc[:, samples_group]
    krep_group = krep_group[krep_group.max(axis=1) >= read_lim]
    if read_mean != None:
      krep_group['Mean'] = krep_group.mean(axis=1)
      krep_group = krep_group[krep_group['Mean'] >= read_mean]
      krep_group = krep_group.drop('Mean', axis=1)
    for row in krep_group.index.values:
      taxid_keeping[str(row)] = taxids_kreports[row]
    group_samples[group] = samples_group
  kreports.to_csv(output_dir+project_name+'_combined_kreport.csv')
  return group_samples, taxid_keeping, kreports

def get_assembly_summaries(assembly_folder, all_domains, representative_only):
  if not all_domains: 
    groups = ['bacteria', 'archaea']
    for group in groups:
      if not os.path.exists(assembly_folder+'assembly_summary_'+group+'.txt'):
        command_download = 'wget -q ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/'+group+'/assembly_summary.txt -O '+assembly_folder+'assembly_summary_'+group+'.txt'
        dl = os.system(command_download)
    assemblies_bacteria = pd.read_csv(assembly_folder+'assembly_summary_bacteria.txt', header=1, index_col=0, sep='\t')
    assemblies_archaea = pd.read_csv(assembly_folder+'assembly_summary_archaea.txt', header=1, index_col=0, sep='\t')
    assemblies = pd.concat([assemblies_bacteria, assemblies_archaea])
  else:
    if not os.path.exists(assembly_folder+'assembly_summary_refseq.txt'):
      command_download = 'wget -q ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt -O '+assembly_folder+'assembly_summary_refseq.txt'
      dl = os.system(command_download)
    assemblies = pd.read_csv(assembly_folder+'assembly_summary_refseq.txt', header=1, index_col=0, sep='\t')
  assemblies_ref = assemblies[assemblies['refseq_category'] == 'reference genome']
  assemblies_rep = assemblies[assemblies['refseq_category'] == 'representative genome']
  if representative_only:
    assemblies = pd.concat([assemblies_ref, assemblies_rep])
  else:
    assemblies_other = assemblies[assemblies['refseq_category'] == 'na']
    assemblies = pd.concat([assemblies_ref, assemblies_rep, assemblies_other])
  assemblies = assemblies.drop_duplicates(subset='species_taxid')
  assemblies['species_taxid'] = assemblies['species_taxid'].astype(str)
  assemblies = assemblies.set_index('species_taxid')
  return assemblies

def write_file(name, list_to_write):
  with open(name, 'w') as f:
    for i in list_to_write:
      w = f.write(i+'\n')
  return

def download_genomes(taxid, assembly_folder, output_dir, all_domains, representative_only, n_proc):
  assemblies = get_assembly_summaries(assembly_folder, all_domains, representative_only)
  download_list, unzip_genomes = [], []
  taxid_list = [tax for tax in taxid]
  assemblies = assemblies.loc[taxid_list, ['ftp_path']]
  for tax in taxid_list:
    genome_name = tax+'_'+taxid[tax].replace(' ', '_')+'.fna'
    if not os.path.exists(output_dir+'genomes/'+genome_name):
      ftp_path = assemblies.loc[tax, 'ftp_path']
      fname = ftp_path.split('/')[-1]
      ftp_path = ftp_path+'/'+fname+'_genomic.fna.gz'
      download_list.append('wget -q '+ftp_path+' -O '+output_dir+'genomes/'+genome_name+'.gz')
      unzip_genomes.append('gunzip '+output_dir+'genomes/'+genome_name+'.gz')
  write_file(output_dir+'genome_download_commands.txt', download_list)
  write_file(output_dir+'genome_unzip_commands.txt', unzip_genomes)
  os.system('python '+dirname(abspath(__file__))+'/run_commands_multiprocessing.py --commands '+output_dir+'genome_download_commands.txt --processors '+str(n_proc))
  os.system('python '+dirname(abspath(__file__))+'/run_commands_multiprocessing.py --commands '+output_dir+'genome_unzip_commands.txt --processors '+str(n_proc))
  got_genomes = {}
  no_genome = []
  for tax in taxid_list:
    genome_name = tax+'_'+taxid[tax].replace(' ', '_')+'.fna'
    if not os.path.exists(output_dir+'genomes/'+genome_name):
      no_genome.append(genome_name)
    else:
      got_genomes[tax] = taxid[tax]
  if len(no_genome) == 0:
    write_file(output_dir+'no_genome.txt', ['Got all genomes!'])
  else:
    write_file(output_dir+'no_genome.txt', no_genome)
  return got_genomes

def extract_reads(taxid, output_dir, samples, fastq_dir, kraken_kreport_dir, kraken_outraw_dir, n_proc):
  taxid_list = ' '.join([tax for tax in taxid])
  command_list = []
  for sample in samples:
    command = 'python '+dirname(abspath(__file__))+'/extract_kraken_reads_modified.py '
    command += '-k '+kraken_outraw_dir+'/'+sample+'.kraken '
    command += '-s '+fastq_dir+'/'+sample+'.fastq '
    command += '-o '+output_dir+'/reads_mapped/'+sample+'.fq '
    command += '-t '+taxid_list+' --include-children '
    if os.path.exists(kraken_kreport_dir+sample+'.kreport'):
      command += '--report '+kraken_kreport_dir+sample+'.kreport --fastq-output'
    else:
      command += '--report '+kraken_kreport_dir+sample+'.kreport --fastq-output'
    command_list.append(command)
  write_file(output_dir+'run_extract_reads_commands.txt', command_list)
  os.system('python '+dirname(abspath(__file__))+'/run_commands_multiprocessing.py --commands '+output_dir+'run_extract_reads_commands.txt --processors '+str(n_proc))
  return

def combine_convert_files(taxid, output_dir, samples, group_samples, n_proc):
  combine_commands = []
  all_fastq = []
  for tax in taxid:
    for group in group_samples:
      samples_in_group = group_samples[group]
      taxid_group = group+'_'+tax+'.fq'
      in_group = False
      for sample in samples_in_group:
        if os.path.exists(output_dir+'reads_mapped/'+sample+'_'+tax+'.fq'):
          command = 'cat '+output_dir+'reads_mapped/'+sample+'_'+tax+'.fq >> '+output_dir+'reads_mapped/'+taxid_group
          combine_commands.append(command)
          all_fastq.append(output_dir+'reads_mapped/'+sample+'_'+tax+'.fq')
          in_group = True
      if in_group:
        all_fastq.append(output_dir+'reads_mapped/'+taxid_group)
  write_file(output_dir+'run_combine_files_commands.txt', combine_commands)
  os.system('python '+dirname(abspath(__file__))+'/run_commands_multiprocessing.py --commands '+output_dir+'run_combine_files_commands.txt --processors '+str(n_proc))
  
  convert_commands = []
  all_files = []
  for fq in all_fastq:
    command =  'python '+dirname(abspath(__file__))+'/convert_fastq_to_fasta.py --fastq '+fq+' --fasta '+fq.replace('.fq', '.fa')
    convert_commands.append(command)
    all_files.append(fq.replace('.fq', ''))
  write_file(output_dir+'run_convert_fastq_commands.txt', convert_commands)
  os.system('python '+dirname(abspath(__file__))+'/run_commands_multiprocessing.py --commands '+output_dir+'run_convert_fastq_commands.txt --processors '+str(n_proc))
  return all_files

def run_quast(all_files, taxid, output_dir, n_proc):
  quast_commands = []
  for file in all_files:
    tid = file.split('_')[-1]
    genome_file = output_dir+'genomes/'+tid+'_'+taxid[tid].replace(' ', '_')+'.fna'
    command = 'quast.py '+file+'.fa -r '+genome_file+' -o '+output_dir+'/QUAST/'+file.split('/')[-1]+' --min-contig 10 --min-identity 80 --no-plots --no-html --no-icarus --silent'
    quast_commands.append(command)
  write_file(output_dir+'run_quast_commands.txt', quast_commands)
  os.system('python '+dirname(abspath(__file__))+'/run_commands_multiprocessing.py --commands '+output_dir+'run_quast_commands.txt --processors '+str(n_proc))
  return

def make_bowtie2_databases(taxid, output_dir, n_proc):
  bowtie2_commands = []
  for tid in taxid:
    genome_file = output_dir+'genomes/'+tid+'_'+taxid[tid].replace(' ', '_')+'.fna'
    bowtie2_file = output_dir+'bowtie2_db/'+tid+'_'+taxid[tid].replace(' ', '_')
    command = 'bowtie2-build '+genome_file+' '+bowtie2_file
    bowtie2_commands.append(command)
  write_file(output_dir+'run_bowtie2_database_commands.txt', bowtie2_commands)
  os.system('python '+dirname(abspath(__file__))+'/run_commands_multiprocessing.py --commands '+output_dir+'run_bowtie2_database_commands.txt --processors '+str(n_proc))
  return
