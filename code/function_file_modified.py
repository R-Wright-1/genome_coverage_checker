import os
import pandas as pd
import sys
from multiprocessing import Pool
from multiprocessing import freeze_support
import argparse
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import pickle
import numpy as np

###FUNCTIONS
#run all of the initial checks to be sure that all of the files exist, etc.
def run_initial_checks(wd, n_proc, sample_dir, sample_name, fastq_dir, kraken_kreport_dir, kraken_outraw_dir, output_dir, assembly_folder, quast_loc, read_lim, read_mean, sample_metadata, species, project_name, rerun):
  #check for location of QUAST
  if quast_loc == None:
    sys.exit('You must supply the location of the quast.py script! Quitting the run now without doing anything.')
  
  #check all folders exist
  for folder in [fastq_dir, kraken_kreport_dir, kraken_outraw_dir, output_dir, assembly_folder, quast_loc]:
    if not os.path.exists(folder):
      sys.exit("This path doesn't exist: "+folder)
  
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
  return wd, n_proc, sample_dir, sample_name, fastq_dir, kraken_kreport_dir, kraken_outraw_dir, output_dir, assembly_folder, quast_loc, read_lim, read_mean, sample_metadata, species, project_name, rerun, md, samples, taxid_name


#function to download a genome file if it doesn't already exist
def download_genome(command):
    if os.path.exists(command.split('-O ')[-1]):
        return
    os.system(command)
    return

#multiprocessing function
def run_multiprocessing(func, i, n_processors):
    with Pool(processes=n_processors) as pool:
        return pool.map(func, i)

#get reads mapped function for taxid
def run_get_reads(in_string):
    all_taxid, sn, fq_dir, ko_dir, kk_dir, o_dir, krep, rerun = in_string.split(',') #all_taxid, sample_name, fastq_dir, kraken_outraw_dir, kraken_kreport_dir, output_dir, krep, rerun
    if rerun == "True":
      command = 'python extract_kraken_reads_modified.py -k '+ko_dir+'/'+sn+'.kraken -s '+fq_dir+'/'+sn+'.fastq -o '+o_dir+'/reads_mapped/'+sn+'.fq -t '+all_taxid+' --include-children --report '+krep+' --fastq-output'
      os.system(command)
      wf = SeqIO.convert(o_dir+'/reads_mapped/'+sn+'.fq', "fastq", o_dir+'/reads_mapped/'+sn+'.fa', "fasta")
    else:
      if os.path.exists(o_dir+'/reads_mapped/'+sn+'.fq'):
          print('Already got '+o_dir+'/reads_mapped/'+sn+'.fq')
      else:
          #command = 'python /home/robyn/coverage_checker/extract_kraken_reads.py -k '+ko_dir+'/'+sn+'.kraken -s '+fq_dir+'/'+sn+'.fastq -o '+o_dir+'/reads_mapped/'+sn+'_'+taxid+'.fq -t '+taxid+' --include-children --report '+krep+' --fastq-output'
          command = 'python extract_kraken_reads_modified.py -k '+ko_dir+'/'+sn+'.kraken -s '+fq_dir+'/'+sn+'.fastq -o '+o_dir+'/reads_mapped/'+sn+'.fq -t '+all_taxid+' --include-children --report '+krep+' --fastq-output'
          os.system(command)
      if not os.path.exists(o_dir+'/reads_mapped/'+sn+'.fa'):
          wf = SeqIO.convert(o_dir+'/reads_mapped/'+sn+'.fq', "fastq", o_dir+'/reads_mapped/'+sn+'.fa', "fasta") #convert the created fastq to fasta
    all_taxids = all_taxid.split(' ')
    fq_not_made = 0
    for taxid in all_taxids:
      try:
        int(taxid)
      except:
        continue
      if rerun == "True":
        wf = SeqIO.convert(o_dir+'/reads_mapped/'+sn+'_'+taxid+'.fq', "fastq", o_dir+'/reads_mapped/'+sn+'_'+taxid+'.fa', "fasta") #convert the created fastq to fasta
      else:
        if not os.path.exists( o_dir+'/reads_mapped/'+sn+'_'+taxid+'.fa'):
            try:
              wf = SeqIO.convert(o_dir+'/reads_mapped/'+sn+'_'+taxid+'.fq', "fastq", o_dir+'/reads_mapped/'+sn+'_'+taxid+'.fa', "fasta") #convert the created fastq to fasta
            except:
              if os.path.exists(o_dir+'/reads_mapped/'+sn+'_'+taxid+'.fq'):
                print("Couldn't convert: "+o_dir+'/reads_mapped/'+sn+'_'+taxid+'.fq')
              else:
                fq_not_made += 0
    if fq_not_made > 0: print(sn, "fastq's weren't made for:", str(fq_not_made), "files")
    return

#Run QUAST function
def run_quast(in_string):
    taxid, sn, fq_dir, ko_dir, kk_dir, o_dir, q_loc, rerun = in_string.split(',') #taxid, sample_name, fastq_dir, kraken_outraw_dir, kraken_kreport_dir, output_dir, quast_loc, rerun
    taxid = str(taxid)
    with open(o_dir+'downloaded_genomes.dict', 'rb') as f:
        downloaded_genomes = pickle.load(f)
    os.chdir(q_loc)
    command = './quast.py '+o_dir+'/reads_mapped/'+sn+'_'+taxid+'.fa -r '+o_dir+'/genomes/'+downloaded_genomes[taxid]+' -o '+o_dir+'/QUAST/'+sn+'_'+taxid+' --min-contig 10 --min-identity 80 --no-plots --no-html --no-icarus --silent'
    if rerun == "True":
      os.system(command)
    else:
      if not os.path.exists(o_dir+'/QUAST/'+sn+'_'+taxid):
          try:
              os.system(command)
          except:
              do_nothing = True
    return

#Make bowtie2 database function
def make_bowtie2_db(in_string):
    taxid, sn, fq_dir, ko_dir, kk_dir, o_dir, krep, rerun = in_string.split(',') #taxid, sample_name, fastq_dir, kraken_outraw_dir, kraken_kreport_dir, output_dir, krep, rerun
    with open(o_dir+'downloaded_genomes.dict', 'rb') as f:
        downloaded_genomes = pickle.load(f)
    if not os.path.exists(o_dir+'/bowtie2_db/'+downloaded_genomes[taxid].replace('.fna', '').replace('.gz', '')+'.1.bt2'):
        command = 'bowtie2-build '+o_dir+'/genomes/'+downloaded_genomes[taxid]+' '+o_dir+'/bowtie2_db/'+downloaded_genomes[taxid].replace('.fna', '').replace('.gz', '')
        os.system(command)
    return

#Run bowtie2 function
def run_bowtie2(in_string):
    taxid, sn, fq_dir, ko_dir, kk_dir, o_dir, krep, rerun = in_string.split(',') #taxid, sample_name, fastq_dir, kraken_outraw_dir, kraken_kreport_dir, output_dir, krep, rerun
    with open(o_dir+'downloaded_genomes.dict', 'rb') as f:
        downloaded_genomes = pickle.load(f)
    if os.path.exists(o_dir+'/bowtie2_mapped/'+sn+'_'+taxid+'.fasta'):
        return
    command1 = 'bowtie2 --threads 4 -x '+o_dir+'/bowtie2_db/'+downloaded_genomes[taxid].replace('.fna', '').replace('.gz', '')+' -U '+o_dir+'/reads_mapped/'+sn+'_'+taxid+'.fq --no-unal -S '+o_dir+'/bowtie2_mapped/'+sn+'_'+taxid+'.sam'
    command2 = 'samtools view -b -F 4 '+o_dir+'/bowtie2_mapped/'+sn+'_'+taxid+'.sam > '+o_dir+'/bowtie2_mapped/'+sn+'_'+taxid+'.bam'
    command3 = 'samtools fasta '+o_dir+'/bowtie2_mapped/'+sn+'_'+taxid+'.bam > '+o_dir+'/bowtie2_mapped/'+sn+'_'+taxid+'.fasta'
    os.system(command1)
    os.system(command2)
    os.system(command3)
    return

# #function to download assembly summaries and get genome download paths
# def get_genomes(t_name, assem_folder, o_dir, taxid_name):
#     #Get the bacteria and archaea assembly summaries - if they don't both exist already, download them
#     if not os.path.exists(assem_folder+'assembly_summary_bacteria.txt') and not os.path.exists(assem_folder+'assembly_summary_archaea.txt'):
#         print('Downloading assembly summaries')
#         command_bac = 'wget -q ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt -O '+assem_folder+'assembly_summary_bacteria.txt'
#         command_arc = 'wget -q ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/archaea/assembly_summary.txt -O '+assem_folder+'assembly_summary_archaea.txt'
#         dl = os.system(command_bac)
#         dl = os.system(command_arc)
#     #Read in the assembly summaries, combine them, and get the reference/representative genomes for each
#     assemblies_bacteria = pd.read_csv(assem_folder+'assembly_summary_bacteria.txt', header=1, index_col=0, sep='\t')
#     assemblies_archaea = pd.read_csv(assem_folder+'assembly_summary_archaea.txt', header=1, index_col=0, sep='\t')
#     assemblies = pd.concat([assemblies_bacteria, assemblies_archaea])
#     assemblies_ref = assemblies[assemblies['refseq_category'] == 'reference genome']
#     assemblies_rep = assemblies[assemblies['refseq_category'] == 'representative genome']
#     assemblies = pd.concat([assemblies_ref, assemblies_rep])
#     assemblies = assemblies.set_index('species_taxid')
#     #make a list of genomes to download - get the ftp paths for each of the taxonomy ID's in the kreport
#     download_list, taxids_using, genomes_got = [], [], []
#     for tid in taxid_name:
#         try:
#             ftp_path = assemblies.loc[tid, 'ftp_path']
#             fname = ftp_path.split('/')[-1]
#             ftp_path = ftp_path+'/'+fname+'_genomic.fna.gz'
#             out_name = str(tid)+'_'+taxid_name[tid].replace(' ', '_')+'.fna.gz'
#             if not os.path.exists(o_dir+'genomes/'+out_name):
#                 download_list.append('wget '+ftp_path+' -O '+o_dir+'genomes/'+out_name)
#             genomes_got.append(out_name)
#             taxids_using.append(tid)
#         except:
#             isnt_prokaryote = True
#     return download_list, taxids_using, genomes_got

#function to download assembly summaries if necessary and open them
def get_assembly_summary(assem_folder, representative_only = True):
    #Get the bacteria and archaea assembly summaries - if they don't both exist already, download them
  if not os.path.exists(assem_folder+'assembly_summary_bacteria.txt') and not os.path.exists(assem_folder+'assembly_summary_archaea.txt'):
    print('Downloading assembly summaries')
    command_bac = 'wget -q ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt -O '+assem_folder+'assembly_summary_bacteria.txt'
    command_arc = 'wget -q ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/archaea/assembly_summary.txt -O '+assem_folder+'assembly_summary_archaea.txt'
    dl = os.system(command_bac)
    dl = os.system(command_arc)
  #Read in the assembly summaries, combine them, and get the reference/representative genomes for each
  assemblies_bacteria = pd.read_csv(assem_folder+'assembly_summary_bacteria.txt', header=1, index_col=0, sep='\t')
  assemblies_archaea = pd.read_csv(assem_folder+'assembly_summary_archaea.txt', header=1, index_col=0, sep='\t')
  assemblies = pd.concat([assemblies_bacteria, assemblies_archaea])
  if representative_only:
    assemblies_ref = assemblies[assemblies['refseq_category'] == 'reference genome']
    assemblies_rep = assemblies[assemblies['refseq_category'] == 'representative genome']
    assemblies = pd.concat([assemblies_ref, assemblies_rep])
  assemblies = assemblies.set_index('species_taxid')
  return assemblies

# def replace_characters(string):
#   new_string = ''
#   for char in string:
#     if char not in ["(",")","/","'"]:
#       new_string += char
#   return new_string

#function to get ftp paths for all genomes that are in the assembly summary as the right type
# def get_genomes(t_name, assem_folder, o_dir, taxid_name, representative_only=True):
#   #make a list of genomes to download - get the ftp paths for each of the taxonomy ID's in the kreport
#   assemblies = get_assembly_summary(assem_folder, representative_only)
#   download_list, taxids_using, genomes_got, no_genome = [], [], [], []
#   for tid in taxid_name:
#     try:
#       ftp_path = assemblies.loc[tid, 'ftp_path']
#       fname = ftp_path.split('/')[-1]
#       ftp_path = ftp_path+'/'+fname+'_genomic.fna.gz'
#       out_name = str(tid)+'_'+taxid_name[tid].replace(' ', '_')+'.fna.gz'
#       if not os.path.exists(o_dir+'genomes/'+out_name):
#         download_list.append('wget '+ftp_path+' -O '+o_dir+'genomes/'+out_name)
#       genomes_got.append(out_name)
#       taxids_using.append(tid)
#     except:
#       no_genome.append(tid)
#       try:
#         ftp_path = assemblies.loc[tid, 'ftp_path']
#         fname = ftp_path.split('/')[-1]
#         ftp_path = ftp_path+'/'+fname+'_genomic.fna.gz'
#         out_name = str(tid)+'_'+taxid_name[tid].replace(' ', '_')+'.fna.gz'
#         print(out_name)
#       except:
#         print(tid)
#   return download_list, taxids_using, genomes_got, no_genome
def get_genomes(t_name, assem_folder, o_dir, taxid_name, representative_only=True):
  #make a list of genomes to download - get the ftp paths for each of the taxonomy ID's in the kreport
  assemblies = get_assembly_summary(assem_folder, representative_only)
  download_list, taxids_using, genomes_got, no_genome = [], [], [], []
  all_id = list(assemblies.index.values)
  for tid in taxid_name:
    #if tid != 391738: continue
    try:
      if all_id.count(tid) > 1:
        assem_list = assemblies.loc[tid, :].values
        ftp_path = ''
        for row in assem_list:
          if 'representative genome' in row:
            something = True
            ftp_path = row[17]
        if ftp_path == '': ftp_path = assem_list[0][17]
      else:
        ftp_path = assemblies.loc[tid, 'ftp_path']
      fname = ftp_path.split('/')[-1]
      ftp_path = ftp_path+'/'+fname+'_genomic.fna.gz'
      out_name = str(tid)+'_'+taxid_name[tid].replace(' ', '_')+'.fna.gz'
      if not os.path.exists(o_dir+'genomes/'+out_name):
        download_list.append('wget '+ftp_path+' -O '+o_dir+'genomes/'+out_name)
      genomes_got.append(out_name)
      taxids_using.append(tid)
    except:
      no_genome.append(tid)
  return download_list, taxids_using, genomes_got, no_genome

# function to get taxids for the list of species names
def get_taxid_for_names(species_list, assem_folder):
  assemblies = get_assembly_summary(assem_folder, representative_only=False)
  taxid_to_name, name_to_taxid = {}, {}
  for r in range(len(assemblies.index.values)):
    taxid = assemblies.index.values[r]
    sp_name = assemblies.iloc[r, assemblies.columns.get_loc('organism_name')]
    taxid_to_name[taxid] = sp_name
    name_to_taxid[sp_name] = taxid
    if sp_name.count(' ') > 1 and 'sp.' not in sp_name:
        sp_name = sp_name.split(' ')[0]+' '+sp_name.split(' ')[1]
        name_to_taxid[sp_name] = taxid
  taxid_name = {}
  taxid_name, dont_have = {}, []
  for sp in species_list:
    if sp in taxid_to_name:
      taxid_name[sp] = taxid_to_name[sp]
    elif sp in name_to_taxid:
      taxid_name[name_to_taxid[sp]] = sp
    else:
      dont_have.append(sp)
  return taxid_name, dont_have

#check whether files with mapped reads were made
def file_check(tid_using, proj_name):
    taxid, sn, fq_dir, ko_dir, kk_dir, o_dir, krep, rerun = tid_using[0].split(',') #taxid, sample_name, fastq_dir, kraken_outraw_dir, kraken_kreport_dir, output_dir, krep, rerun
    with open(o_dir+proj_name+'_unmade_files.txt', 'w') as f:
      w = f.write('\n')
    for input_string in tid_using:
        taxid, sn, fq_dir, ko_dir, kk_dir, o_dir, krep, rerun = input_string.split(',') #taxid, sample_name, fastq_dir, kraken_outraw_dir, kraken_kreport_dir, output_dir, krep, rerun
        if not os.path.exists(o_dir+'/reads_mapped/'+sn+'_'+str(taxid)+'.fq'):
          with open(o_dir+proj_name+'_unmade_files.txt', 'a') as f:
            w = f.write(o_dir+'/reads_mapped/'+sn+'_'+str(taxid)+'.fq\n')
        if not os.path.exists(o_dir+'/reads_mapped/'+sn+'_'+str(taxid)+'.fa'):
            with open(o_dir+proj_name+'_unmade_files.txt', 'a') as f:
              w = f.write(o_dir+'/reads_mapped/'+sn+'_'+str(taxid)+'.fa\n')
    return

#check whether genome files were made
def genome_check(dl_list, proj_name):
    o_dir = dl_list[0].split('-O ')[-1].split('genomes/')[0]
    with open(o_dir+proj_name+'_genomes_not_downloaded.txt', 'w') as f:
      w = f.write('\n')
    for genome in dl_list:
        gen_path = genome.split('-O ')[-1]
        if not os.path.exists(gen_path):
            with open(o_dir+proj_name+'_genomes_not_downloaded.txt', 'a') as f:
              w = f.write(gen_path+'\n')
    return

#check whether quast was run successfully
def quast_check(tid_using, proj_name):
    taxid, sn, fq_dir, ko_dir, kk_dir, o_dir, krep, rerun = tid_using[0].split(',') #taxid, sample_name, fastq_dir, kraken_outraw_dir, kraken_kreport_dir, output_dir, krep, rerun
    with open(o_dir+proj_name+'_quast_not_run.txt', 'w') as f:
      w = f.write('\n')
    for t in tid_using:
        taxid, sn, fq_dir, ko_dir, kk_dir, o_dir, krep, rerun = t.split(',') #taxid, sample_name, fastq_dir, kraken_outraw_dir, kraken_kreport_dir, output_dir, krep, rerun
        # with open(o_dir+'taxid_name.dict', 'rb') as f:
        #     taxid_name = pickle.load(f)
        # tax_name = taxid_name[int(taxid)]
        folder = o_dir+'/QUAST/'+sn+'_'+taxid+'/'
        if not os.path.exists(folder+'report.tsv'):
            with open(o_dir+proj_name+'_quast_not_run.txt', 'a') as f:
              w = f.write(o_dir+'/QUAST/'+sn+'_'+taxid+'/report.tsv'+'\n')
    return

#combine files for multifile runs
def combine_files(str_file):
  fn, o_dir, all_taxid = str_file.split(',')
  all_taxid_list = all_taxid.split(' ')
  for taxid_sn in all_taxid_list:
    if os.path.exists(o_dir+'reads_mapped/'+taxid_sn+'.fq'):
      command1 = 'cat '+o_dir+'reads_mapped/'+taxid_sn+'.fq >> '+o_dir+'reads_mapped/'+fn+'.fq'
      os.system(command1)
    if os.path.exists(o_dir+'reads_mapped/'+taxid_sn+'.fa'):
      command2 = 'cat '+o_dir+'reads_mapped/'+taxid_sn+'.fa >> '+o_dir+'reads_mapped/'+fn+'.fa'
      os.system(command2)
  return

#get coverage across genome
def get_coverage(t):
  all_coverage = {}
  string = t.split(',')
  tax, sample, o_dir = string[0], string[1], string[5]
  if os.path.exists(o_dir+'pickle_coverage/'+str(tax)+'_'+sample+'.dict'):
    return
  try:
    folder = o_dir+'QUAST/'+sample+'_'+tax+'/'
    report = pd.read_csv(folder+'report.tsv', index_col=0, header=0, sep='\t')
    ref_chromosomes, ref_chromosome_dict, chromo_before, chromosome_lines = [], {}, 0, []
  except:
    with open(o_dir+'genomes_not_working.txt', 'a') as f:
      w = f.write(sample+'_'+tax+'\n')
    return
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
    length = report.loc['Reference length', sample+'_'+tax]
    alignments = pd.read_csv(folder+'contigs_reports/all_alignments_'+sample+'_'+tax+'.tsv', index_col=None, header=0, sep='\t')
  except:
    with open(o_dir+'genomes_not_working.txt', 'a') as f:
      w = f.write(sample+'_'+tax+'\n')
    return
  genome_coverage, genome_identity, all_starts = [], [], []
  alignments = alignments[alignments['S1'] != 'CONTIG']
  for row in alignments.index:
    try:
      start, end = int(alignments.loc[row, 'S1']), int(alignments.loc[row, 'E1'])
    except:
      continue
    genome_identity.append(float(alignments.loc[row, 'IDY']))
    adding = ref_chromosome_dict[alignments.loc[row, 'Reference']]
    s, e = min([start, end])+adding, max([start, end])+adding
    middle = np.mean([s, e])
    all_starts.append(s)
    rounded = round(middle / 500.0) * 500.0
    genome_coverage.append(rounded)
  genome_coverage_set = list(set(genome_coverage))
  all_coverage[sample+'_'+tax] = [genome_coverage, genome_identity, all_starts, length]
  with open(o_dir+'pickle_coverage/'+str(tax)+'_'+sample+'.dict', 'wb') as f:
    du = pickle.dump(all_coverage, f)
  return

