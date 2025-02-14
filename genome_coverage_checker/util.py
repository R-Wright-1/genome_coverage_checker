#!/usr/bin/env python

import os
import pandas as pd
import sys
import numpy as np
from os.path import abspath, dirname

def get_checkpoint(output_dir):
  cp = ''
  for row in open(output_dir+'/checkpoint.txt', 'r'):
    cp += row
  return cp

def update_checkpoint(output_dir, cp):
  with open(output_dir+'/checkpoint.txt', 'w') as f:
    c = f.write(str(cp))
  return cp
  

def run_initial_checks(wd, n_proc, sample_dir, sample_name, fastq_dir, kraken_kreport_dir, kraken_outraw_dir, output_dir, assembly_folder, read_lim, read_mean, sample_metadata, species, project_name, rerun, genome_dir, bowtie2_db_dir):
  #check all folders exist
  for folder in [fastq_dir, kraken_kreport_dir, kraken_outraw_dir, output_dir, assembly_folder]:
    if not os.path.exists(folder):
      if folder == output_dir:
        print(output_dir+" doesn't already exist. Making it now.")
        os.system('mkdir '+output_dir)
      else:
        sys.exit("This path doesn't exist: "+folder)

  fastq_dir, kraken_kreport_dir, kraken_outraw_dir, output_dir, assembly_folder, genome_dir, bowtie2_db_dir = fastq_dir+'/', kraken_kreport_dir+'/', kraken_outraw_dir+'/', output_dir+'/', assembly_folder+'/', genome_dir+'/', bowtie2_db_dir+'/'
  
  #directories to make if they don't exist
  all_dirs = [genome_dir, output_dir+'/reads_mapped/', output_dir+'/QUAST/', bowtie2_db_dir, output_dir+'/bowtie2_mapped/', output_dir+'/coverage/', output_dir+'/pickle_intermediates/']
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
  return wd, n_proc, sample_dir, sample_name, fastq_dir, kraken_kreport_dir, kraken_outraw_dir, output_dir, assembly_folder, read_lim, read_mean, sample_metadata, species, project_name, rerun, md, samples, taxid_name, genome_dir, bowtie2_db_dir

def get_kreports(samples, kraken_kreport_dir, output_dir, md, read_lim, read_mean, project_name, taxid_name):
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
  for tax in taxid_name:
    if tax in kreports.index.values:
      taxid_keeping[tax] = taxid_name[tax]
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
  for tax in taxid_keeping:
    sp_name = taxid_keeping[tax].replace(' ', '_')
    for s in sp_name:
      if not s.isalnum():
        if s in ['.', '_', '-']: continue
        sp_name = sp_name.replace(s, '-')
    taxid_keeping[tax] = sp_name
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
    assemblies = pd.read_csv(assembly_folder+'assembly_summary_refseq.txt', header=1, index_col=0, sep='\t', low_memory=False)
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
  assemblies.to_csv(assembly_folder+'assemblies_using.csv')
  return assemblies

def write_file(name, list_to_write):
  with open(name, 'w') as f:
    for i in list_to_write:
      w = f.write(i+'\n')
  return

def download_genomes(taxid, assembly_folder, output_dir, all_domains, representative_only, n_proc, genome_dir):
  assemblies = get_assembly_summaries(assembly_folder, all_domains, representative_only)
  download_list, unzip_genomes = [], []
  no_genome, taxid_list = [], []
  for tax in taxid:
    if tax in assemblies.index.values:
      taxid_list.append(tax)
    else:
      if representative_only:
        no_genome.append(tax+'_'+taxid[tax]+'.fna: not in NCBI assembly summary. You could try setting representative only to false incase this taxonomy ID has no representative genome.')
      else:
        no_genome.append(tax+'_'+taxid[tax]+'.fna: not in NCBI assembly summary')
  assemblies = assemblies.loc[taxid_list, ['ftp_path']]
  for tax in taxid_list:
    genome_name = tax+'_'+taxid[tax]+'.fna'
    if not os.path.exists(genome_dir+genome_name):
      if os.path.exists(genome_dir+genome_name+'.gz'):
        unzip_genomes.append('gunzip '+genome_dir+genome_name+'.gz')
        continue
      ftp_path = assemblies.loc[tax, 'ftp_path']
      fname = ftp_path.split('/')[-1]
      ftp_path = ftp_path+'/'+fname+'_genomic.fna.gz'
      download_list.append('wget -q '+ftp_path+' -O '+genome_dir+genome_name+'.gz')
      unzip_genomes.append('gunzip '+genome_dir+genome_name+'.gz')
  write_file(output_dir+'genome_download_commands.txt', download_list)
  write_file(output_dir+'genome_unzip_commands.txt', unzip_genomes)
  os.system('python '+dirname(abspath(__file__))+'/run_commands_multiprocessing.py --commands '+output_dir+'genome_download_commands.txt --processors '+str(n_proc))
  os.system('python '+dirname(abspath(__file__))+'/run_commands_multiprocessing.py --commands '+output_dir+'genome_unzip_commands.txt --processors '+str(n_proc))
  got_genomes = {}
  for tax in taxid_list:
    genome_name = tax+'_'+taxid[tax]+'.fna'
    if not os.path.exists(genome_dir+genome_name):
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

def combine_convert_files(taxid, output_dir, samples, group_samples, n_proc, skip_duplicate_check=False):
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

  if not skip_duplicate_check:
    #remove duplicate reads next
    remove_duplicates_files = [output_dir+'reads_mapped/'+f for f in os.listdir(output_dir+'reads_mapped/') if '.fq' in f]
    write_file(output_dir+'run_remove_duplicate_reads.txt', remove_duplicates_files)
    os.system('python '+dirname(abspath(__file__))+'/remove_duplicate_reads.py --files '+output_dir+'run_remove_duplicate_reads.txt --processors '+str(n_proc))

  convert_commands = []
  all_files = []
  for fq in all_fastq:
    command =  'python '+dirname(abspath(__file__))+'/convert_fastq_to_fasta.py --fastq '+fq+' --fasta '+fq.replace('.fq', '.fa')
    convert_commands.append(command)
    all_files.append(fq.replace('.fq', ''))
  write_file(output_dir+'run_convert_fastq_commands.txt', convert_commands)
  os.system('python '+dirname(abspath(__file__))+'/run_commands_multiprocessing.py --commands '+output_dir+'run_convert_fastq_commands.txt --processors '+str(n_proc))
  return all_files

def run_quast(all_files, taxid, output_dir, n_proc, genome_dir):
  quast_commands = []
  for file in all_files:
    tid = file.split('_')[-1]
    genome_file = genome_dir+tid+'_'+taxid[tid]+'.fna'
    command = 'quast.py '+file+'.fa -r '+genome_file+' -o '+output_dir+'/QUAST/'+file.split('/')[-1]+' --min-contig 10 --min-identity 80 --no-plots --no-html --no-icarus --silent'
    command += ' >> '+output_dir+'/quast_terminal_output.txt'
    quast_commands.append(command)
  write_file(output_dir+'run_quast_commands.txt', quast_commands)
  os.system('python '+dirname(abspath(__file__))+'/run_commands_multiprocessing.py --commands '+output_dir+'run_quast_commands.txt --processors '+str(n_proc))
  quast_error = []
  for file in all_files:
    if not os.path.exists(output_dir+'/QUAST/'+file.split('/')[-1]):
      quast_error.append(file.split('/')[-1])
  with open(output_dir+'quast_error.txt', 'w') as f:
    for fn in quast_error:
      w = f.write(fn+'\n')
  return

def make_bowtie2_databases(taxid, output_dir, n_proc, bowtie2_db_dir, genome_dir):
  bowtie2_commands = []
  for tid in taxid:
    genome_file = genome_dir+tid+'_'+taxid[tid]+'.fna'
    bowtie2_file = bowtie2_db_dir+tid+'_'+taxid[tid]
    bowtie2_files = [bowtie2_file+'.1.bt2', bowtie2_file+'.2.bt2', bowtie2_file+'.3.bt2', bowtie2_file+'.4.bt2', bowtie2_file+'.rev.1.bt2', bowtie2_file+'.rev.2.bt2']
    got_bt2_files = True
    for bt2f in bowtie2_files:
      if not os.path.exists(bt2f):
        got_bt2_files = False
        break
    if got_bt2_files: continue
    command = 'bowtie2-build --quiet '+genome_file+' '+bowtie2_file
    command += ' >> '+output_dir+'/bowtie2_terminal_output.txt'
    bowtie2_commands.append(command)
  write_file(output_dir+'run_bowtie2_database_commands.txt', bowtie2_commands)
  os.system('python '+dirname(abspath(__file__))+'/run_commands_multiprocessing.py --commands '+output_dir+'run_bowtie2_database_commands.txt --processors '+str(n_proc))
  return

def run_bowtie2(all_files, taxid, output_dir, n_proc, bowtie2_setting, bowtie2_db_dir):
  bowtie2_commands, view_commands, fasta_commands = [], [], []
  out_names = []
  for file in all_files:
    tid = file.split('_')[-1]
    bowtie2_file = bowtie2_db_dir+tid+'_'+taxid[tid]
    out_name = output_dir+'/bowtie2_mapped/'+file.split('/')[-1]
    out_names.append(out_name)
    command_bt2 = 'bowtie2 --quiet --threads 1 --'+bowtie2_setting+' -x '+bowtie2_file+ ' -U '+file+'.fq --no-unal -S '+out_name+'.sam >> '+output_dir+'/bowtie2_terminal_output.txt'
    #command_view = 'samtools view -b -F 4 '+out_name+'.sam > '+out_name+'.bam'
    #command_fasta = 'samtools fasta '+out_name+'.bam -0 '+out_name+'.fasta --verbosity 0'
    bowtie2_commands.append(command_bt2)
    #view_commands.append(command_view)
    #fasta_commands.append(command_fasta)
  write_file(output_dir+'run_bowtie2_commands.txt', bowtie2_commands)
  os.system('python '+dirname(abspath(__file__))+'/run_commands_multiprocessing.py --commands '+output_dir+'run_bowtie2_commands.txt --processors '+str(n_proc))
  bowtie2_error = []
  for out_name in out_names:
    if os.path.exists(out_name+'.sam'):
      command_view = 'samtools view -b -F 4 '+out_name+'.sam > '+out_name+'.bam'
      view_commands.append(command_view)
    else:
      bowtie2_error.append(out_name)
  write_file(output_dir+'run_view_commands.txt', view_commands)
  os.system('python '+dirname(abspath(__file__))+'/run_commands_multiprocessing.py --commands '+output_dir+'run_view_commands.txt --processors '+str(n_proc))
  view_error = []
  for out_name in out_names:
    if os.path.exists(out_name+'.bam'):
      command_fasta = 'samtools fasta '+out_name+'.bam -0 '+out_name+'.fasta --verbosity 0'
      fasta_commands.append(command_fasta)
    else:
      view_error.append(out_name)
  write_file(output_dir+'run_fasta_commands.txt', fasta_commands)
  os.system('python '+dirname(abspath(__file__))+'/run_commands_multiprocessing.py --commands '+output_dir+'run_fasta_commands.txt --processors '+str(n_proc))
  fasta_error = []
  for out_name in out_names:
    if not os.path.exists:
      fasta_error.append(out_name)
  with open(output_dir+'bowtie2_errors.txt', 'w') as f:
    for fn in bowtie2_error:
      w = f.write(fn+' Bowtie2 error\n')
    for fn in view_error:
      w = f.write(fn+' View error\n')
    for fn in fasta_error:
      w = f.write(fn+' Fasta error\n')
  return


def get_coverage_across_genomes(all_files, taxid, output_dir, n_proc):
  get_coverage = []
  for f in range(len(all_files)):
    fn = all_files[f].split('/')[-1]
    quast_dir = output_dir+'/QUAST/'+fn
    get_coverage.append(quast_dir)
  write_file(output_dir+'get_genome_coverage_folders.txt', get_coverage)
  os.system('python '+dirname(abspath(__file__))+'/get_genome_coverage.py --folders '+output_dir+'get_genome_coverage_folders.txt --processors '+str(n_proc))
  return


def collate_output(all_files, taxid, output_dir, kreports, samples, group_samples, skip_bowtie2, skip_coverage):
  all_files = [f.split('/')[-1] for f in all_files]
  #get quast outputs
  quast_out = {}
  for f in all_files:
    report = output_dir+'QUAST/'+f+'/report.tsv'
    if os.path.exists(report):
      report = pd.read_csv(report, index_col=0, header=0, sep='\t')
      ref_len, ref_gc = report.loc['Reference length', f], report.loc['Reference GC (%)', f]
      nreads = float(report.loc['# contigs (>= 0 bp)', f])
      quast_gc = report.loc['GC (%)', f]
      try:
        aligned_length = report.loc['Total aligned length', f]
      except:
        aligned_length = ''
      try:
        unaligned = float(report.loc['# unaligned contigs', f].split(' ')[0])
      except:
        unaligned = nreads
      try:
        genome_frac, dup_ratio = report.loc['Genome fraction (%)', f], report.loc['Duplication ratio', f]
      except:
        genome_frac, dup_ratio = 0, ''
      quast_out[f] = [ref_len, ref_gc, nreads, quast_gc, genome_frac, dup_ratio, unaligned, aligned_length]

  #get bowtie2 outputs
  if not skip_bowtie2:
    bowtie2_out = {}
    for f in all_files:
      count = 0
      if os.path.exists(output_dir+'bowtie2_mapped/'+f+'.fasta'):
        for row in open(output_dir+'bowtie2_mapped/'+f+'.fasta', 'r'):
          if row[0] == '>': count += 1
      bowtie2_out[f] = count

  #get kraken counts for each sample or each group of samples
  kraken_counts = {}
  tax_list = [t for t in taxid]
  for sample in samples:
    group_samples[sample] = [sample]
  for group in group_samples:
    for tax in tax_list:
      krak_red = kreports.loc[int(tax), group_samples[group]].values
      kraken_counts[group+'_'+tax] = sum(krak_red)

  #now compile all together
  first_row = ['Sample', 'taxid', 'Species name', 'Reference genome length (bp)', 'Kraken reads assigned', 'QUAST reads mapped', 'QUAST genome fraction (%)', 'QUAST duplication ratio', 'QUAST aligned length']
  if not skip_coverage:
    first_row.append('QUAST identity of mapped reads (%)')
  if not skip_bowtie2:
    first_row.append('Bowtie2 reads mapped')
  all_out = []
  for group in group_samples:
    for tax in taxid:
      if kraken_counts[group+'_'+tax] == 0:
        this_sample = [group, tax, taxid[tax], '', kraken_counts[group+'_'+tax], '', '', '']
        if not skip_coverage:
          this_sample.append('')
        if not skip_bowtie2:
          this_sample.append('')
        all_out.append(this_sample)
        continue
      try:
        quast_sample = quast_out[group+'_'+tax] #ref_len, ref_gc, nreads, quast_gc, genome_frac, dup_ratio, unaligned, aligned_length
        this_sample = [group, tax, taxid[tax], quast_sample[0], kraken_counts[group+'_'+tax], quast_sample[2]-quast_sample[6], quast_sample[4], quast_sample[5], quast_sample[7]]
        if not skip_coverage:
          if not quast_sample[2] == 0:
            try:
              for row in open(output_dir+'coverage/'+group+'_'+tax+'.txt', 'r'):
                if 'genome_identity' in row:
                  row = row.replace('genome_identity: ', '').replace('\n', '')
                  if row == '':
                    this_sample.append('')
                  else:
                    iden = row.split(',')
                    iden = [float(r) for r in iden]
                    iden = np.mean(iden)
                    this_sample.append(iden)
            except:
              this_sample.append('')
          else:
            this_sample.append('')
      except:
        this_sample = [group, tax, taxid[tax], '', kraken_counts[group+'_'+tax], '', '', '', '', '']
      if not skip_bowtie2:
        try:
          this_sample.append(bowtie2_out[group+'_'+tax])
        except:
          this_sample.append('')
      all_out.append(this_sample)
  out_df = pd.DataFrame(all_out, columns=first_row)
  out_df['QUAST reads mapped'], out_df['Kraken reads assigned'] = pd.to_numeric(out_df['QUAST reads mapped']), pd.to_numeric(out_df['Kraken reads assigned'])
  out_df['Proportion kraken reads mapped with QUAST'] = out_df['QUAST reads mapped']/out_df['Kraken reads assigned']
  if not skip_bowtie2:
    out_df['Bowtie2 reads mapped'] = pd.to_numeric(out_df['Bowtie2 reads mapped'])
    out_df['Proportion kraken reads mapped with Bowtie2'] = out_df['Bowtie2 reads mapped']/out_df['Kraken reads assigned']
  if not skip_bowtie2:
    out_df = out_df.loc[:, ['Sample', 'taxid', 'Species name', 'Reference genome length (bp)', 'Kraken reads assigned', 'QUAST reads mapped', 'QUAST genome fraction (%)', 'QUAST duplication ratio', 'QUAST aligned length', 'QUAST identity of mapped reads (%)', 'Bowtie2 reads mapped', 'Proportion kraken reads mapped with QUAST', 'Proportion kraken reads mapped with Bowtie2']]
  else:
    out_df = out_df.loc[:, ['Sample', 'taxid', 'Species name', 'Reference genome length (bp)', 'Kraken reads assigned', 'QUAST reads mapped', 'QUAST genome fraction (%)', 'QUAST duplication ratio', 'QUAST aligned length', 'QUAST identity of mapped reads (%)', 'Proportion kraken reads mapped with QUAST']]
  out_df.to_csv(output_dir+'coverage_checker_output.tsv', sep='\t', index=False)
  return

def clean_up(output_dir):
  files = ['genome_download_commands.txt', 'genome_unzip_commands.txt', 'get_genome_coverage_folders.txt', 'run_bowtie2_commands.txt', 'run_bowtie2_database_commands.txt', 'run_combine_files_commands.txt', 'run_convert_fastq_commands.txt', 'run_extract_reads_commands.txt', 'run_fasta_commands.txt', 'run_quast_commands.txt', 'run_view_commands.txt', 'run_remove_duplicate_reads.txt', 'quast_terminal_output.txt', 'bowtie2_terminal_output.txt']
  for f in files:
    os.system('rm '+output_dir+'/'+f)
  return

