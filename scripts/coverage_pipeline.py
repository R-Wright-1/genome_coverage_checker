#!/usr/bin/env python

import os
import pandas as pd
import sys
import argparse
import pickle
from genome_coverage_checker.util import *
import time

start_time = time.time()

# 0. Get command line arguments
parser = argparse.ArgumentParser(description='This script is to check which taxa reads have been assigned to by Kraken, pull out these reads, download reference genomes for the taxa, and map the reads to the reference genomes.')
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
parser.add_argument('--genome_dir', dest='genome_dir', default=None,
                    help="The directory to save the downloaded genomes to. The default is for them to be added to the output_dir")
parser.add_argument('--bowtie2_db_dir', dest='bowtie2_db_dir', default=None,
                    help="The directory to save the bowtie2 database to. The default is for them to be added to the output_dir")
# parser.add_argument('--quast_loc', dest='quast_loc', default=None,
#                     help="The location of the executable quast.py script (folder containing QUAST)")
parser.add_argument('--read_lim', dest='read_lim', default=None,
                    help="Look at all taxa with > this number of reads mapped by Kraken")
parser.add_argument('--read_mean', dest='read_mean', default=None,
                    help="Look at all taxa with > this number of mean reads mapped by Kraken within sample groupings")
parser.add_argument('--assembly_folder', dest='assembly_folder', default=None,
                    help="The folder containing the assembly summaries for bacteria and archaea")
parser.add_argument('--sample_metadata', dest='sample_metadata', default=None,
                    help="Location of the sample metadata file. It is expected that this will be a CSV (comma separated) file with a header and two columns. The first column contains sample names and the second contains the sample groupings")
parser.add_argument('--species', dest='species', default=None,
                    help="Location of the species file. An optional file containing a list of taxonomy ID's or species names to include")
parser.add_argument('--project_name', dest='project_name', default=None,
                    help="Name of this project")
parser.add_argument('--rerun', dest='rerun', default=False, action='store_true',
                    help="If this flag is added, it will re-run everything regardless of whether the folder and a checkpoint already exist")
parser.add_argument('--all_domains', dest='all_domains', default=False, action='store_true',
                    help="The default for coverage checker is for only the genomes of prokaryotes to be downloaded and checked. If you'd like to include all genomes then add this flag.")
parser.add_argument('--representative_only', dest='representative_only', default=False, action='store_true',
                    help="The default for coverage checker is to use all genomes. If you'd like to limit the checking to only representative and reference genomes then add this flag.")
parser.add_argument('--bowtie2_setting', dest='bowtie2_setting', default='sensitive',
                    help="Th default bowtie2 setting to use. Options are very-fast, fast, sensitive or very-sensitive. Default is sensitive (same as Bowtie2 default)")
parser.add_argument('--skip_bowtie2', dest='skip_bowtie2', default=False, action='store_true',
                    help="Whether to skip bowtie2 mapping. Default is to run bowtie2, but bowtie2 can be difficult to install on a mac so add this flag if you don't want to install bowtie2.")
parser.add_argument('--skip_coverage', dest='skip_coverage', default=False, action='store_true',
                    help="Whether to skip getting coverage for all reads across the genomes. This can take a little while, so you can skip it if you don't think that this output is useful to you.")
parser.add_argument('--skip_cleanup', dest='skip_cleanup', default=False, action='store_true',
                    help="If you want to skip the cleanup. This will keep the intermediate files containing some of the commands run, e.g. for making bowtie2 databases. They may be helpful if you're trying to troubleshoot issues.")
parser.add_argument('--skip_duplicate_check', dest='skip_duplicate_check', default=False, action='store_true',
                    help="If you want to skip the check for duplicates within the fastq files. Note that this step can take a while if you have a lot of samples - it was mainly added because you'll get some weird QUAST/Bowtie2 results if you have duplicate reads in your files. This can happen if you rerun coverage checker using the same output folder.")


#Read in the command line arguments
args = parser.parse_args()
sample_name, n_proc, sample_dir, fastq_dir, kraken_kreport_dir, kraken_outraw_dir, output_dir, genome_dir, bowtie2_db_dir = args.sample_name, int(args.n_proc), args.sample_dir, args.fastq_dir, args.kraken_kreport_dir, args.kraken_outraw_dir, args.output_dir, args.genome_dir, args.bowtie2_db_dir
if sample_dir != None:
  fastq_dir, kraken_kreport_dir, kraken_outraw_dir, output_dir = sample_dir, sample_dir, sample_dir, sample_dir
if genome_dir == None:
  genome_dir = output_dir+'/genomes/'
if bowtie2_db_dir == None:
  bowtie2_db_dir = output_dir+'bowtie2_db/'
read_lim, read_mean, assembly_folder = args.read_lim, args.read_mean, args.assembly_folder
if assembly_folder == None:
  assembly_folder = output_dir
sample_metadata, species, project_name, rerun, all_domains, representative_only, skip_bowtie2, skip_coverage, skip_cleanup, skip_duplicate_check, bowtie2_setting = args.sample_metadata, args.species, args.project_name, args.rerun, args.all_domains, args.representative_only, args.skip_bowtie2, args.skip_coverage, args.skip_cleanup, args.skip_duplicate_check, args.bowtie2_setting
wd = os.getcwd()
if read_lim == None: read_lim = 0
else: read_lim = int(read_lim)
if read_mean == None: read_mean = 0
else: read_mean = int(read_mean)

#check whether we've already run this and which checkpoint we're at
if rerun:
  cp = '0'
elif os.path.exists(output_dir+'/checkpoint.txt'):
  cp = get_checkpoint(output_dir)
  print("Got information from previous run: already at "+cp)
else:
  cp = '0'

#if not rerun and cp not 0, check whether the arguments given are the same?
all_args = [wd, n_proc, sample_dir, sample_name, fastq_dir, kraken_kreport_dir, kraken_outraw_dir, output_dir, assembly_folder, read_lim, read_mean, sample_metadata, species, project_name, rerun, all_domains, representative_only, skip_bowtie2, skip_coverage, skip_cleanup, skip_duplicate_check, bowtie2_setting, genome_dir, bowtie2_db_dir]
if cp != '0':
  if os.path.exists(output_dir+'/pickle_intermediates/args.pickle'):
    with open(output_dir+'/pickle_intermediates/args.pickle', 'rb') as f:
      old_args = pickle.load(f)
    if all_args[:-2] != old_args[:-2]:
      sys.exit("You've tried to run coverage checker again in a directory that already exists with different command line options. Please provide a new option for --output_dir and try again.")
    

# 1. Run the initial checks
if cp == '0':
  print("Running initial checks")
  wd, n_proc, sample_dir, sample_name, fastq_dir, kraken_kreport_dir, kraken_outraw_dir, output_dir, assembly_folder, read_lim, read_mean, sample_metadata, species, project_name, rerun, md, samples, taxid_name, genome_dir, bowtie2_db_dir = run_initial_checks(wd, n_proc, sample_dir, sample_name, fastq_dir, kraken_kreport_dir, kraken_outraw_dir, output_dir, assembly_folder, read_lim, read_mean, sample_metadata, species, project_name, rerun, genome_dir, bowtie2_db_dir) #run all of the initial checks to ensure that all of the folders and files exist before starting to try and run anything
  all_args = [wd, n_proc, sample_dir, sample_name, fastq_dir, kraken_kreport_dir, kraken_outraw_dir, output_dir, assembly_folder, read_lim, read_mean, sample_metadata, species, project_name, rerun, all_domains, representative_only, skip_bowtie2, skip_coverage, skip_cleanup, skip_duplicate_check, bowtie2_setting]
  with open(output_dir+'/pickle_intermediates/args.pickle', 'wb') as f:
    pickle.dump(all_args, f)
  cp = update_checkpoint(output_dir, "1_initial_checks_run")
  print("Completed check-point 1 initial checks")
else:
  print("Skipping initial checks because check-point wasn't 0")
  

# 2. Read in all kreports and combine them to get the genomes that we'll need
if cp == "1_initial_checks_run":
  print("Reading in kreports and combining them")
  group_samples, taxid, kreports = get_kreports(samples, kraken_kreport_dir, output_dir, md, read_lim, read_mean, project_name, taxid_name)
  names, objects = ['group_samples', 'taxid', 'kreports'], [group_samples, taxid, kreports]
  for n in range(len(names)):
    with open(output_dir+'/pickle_intermediates/'+names[n]+'.pickle', 'wb') as f:
      pickle.dump(objects[n], f)
  cp = update_checkpoint(output_dir, "2_combined_kreports")
  print("Completed check-point 2 combined kreports")
else:
  print("Skipping 2_combined_kreports because check-point wasn't 1_initial_checks_run")
  names, objects = ['group_samples', 'taxid', 'kreports'], []
  for n in range(len(names)):
    with open(output_dir+'/pickle_intermediates/'+names[n]+'.pickle', 'rb') as f:
      objects.append(pickle.load(f))
  group_samples, taxid, kreports = objects
  

# 3. Download all genomes
if cp == "2_combined_kreports":
  print("Downloading all genomes")
  taxid = download_genomes(taxid, assembly_folder, output_dir, all_domains, representative_only, n_proc, genome_dir)
  with open(output_dir+'/pickle_intermediates/taxid.pickle', 'wb') as f:
    pickle.dump(taxid, f)
  cp = update_checkpoint(output_dir, "3_downloaded_genomes")
  print("Completed check-point 3 downloaded genomes")
else:
  print("Skipping 3_downloaded_genomes because check-point wasn't 2_combined_kreports")
  with open(output_dir+'/pickle_intermediates/taxid.pickle', 'rb') as f:
    taxid = pickle.load(f)

# 4. Extract the reads for each taxonomy ID
if cp == "3_downloaded_genomes":
  print("Extracting reads for all taxonomy ID's")
  extract_reads(taxid, output_dir, samples, fastq_dir, kraken_kreport_dir, kraken_outraw_dir, n_proc)
  cp = update_checkpoint(output_dir, "4_extracted_reads")
  print("Completed check-point 4 extracted all reads")
else:
  print("Skipping 4_extracted_reads because check-point wasn't 3_downloaded_genomes")

# 5. Combine the files for each taxonomy ID (across multiple samples for the sample groupings)
if cp == "4_extracted_reads":
  print("Combining files for each taxonomy ID")
  all_files = combine_convert_files(taxid, output_dir, samples, group_samples, n_proc, skip_duplicate_check)
  with open(output_dir+'/pickle_intermediates/all_files.pickle', 'wb') as f:
    pickle.dump(all_files, f)
  cp = update_checkpoint(output_dir, "5_combined_files")
  print("Completed check-point 5 combined files")
else:
  print("Skipping 5_combined_files because check-point wasn't 4_extracted_reads")
  with open(output_dir+'/pickle_intermediates/all_files.pickle', 'rb') as f:
    all_files = pickle.load(f)

# 6. Run QUAST
if cp == "5_combined_files":
  print("Running QUAST")
  run_quast(all_files, taxid, output_dir, n_proc, genome_dir)
  cp = update_checkpoint(output_dir, "6_quast_run")
  print("Completed check-point 6 QUAST run")
else:
  print("Skipping 6_quast_run because check-point wasn't 5_combined_files")

# 7. Make bowtie2 databases & run bowtie2 for all taxonomy ID's
if cp == "6_quast_run":
  if skip_bowtie2:
    cp = update_checkpoint(output_dir, "7_bowtie2_run")
    print("Skipped check-point 7 Bowtie2 run")
  else:
    print("Running Bowtie2")
    make_bowtie2_databases(taxid, output_dir, n_proc, bowtie2_db_dir, genome_dir)
    run_bowtie2(all_files, taxid, output_dir, n_proc, bowtie2_setting, bowtie2_db_dir)
    cp = update_checkpoint(output_dir, "7_bowtie2_run")
    print("Completed check-point 7 Bowtie2 run")
else:
  print("Skipping 7_bowtie2_run because check-point wasn't 6_quast_run")

# 8. Get the coverage and mapping of reads across the genomes
if cp == "7_bowtie2_run":
  if skip_coverage:
    cp = update_checkpoint(output_dir, "8_got_coverage")
    print("Skipped check-point 8 getting coverage across genomes")
  else:
    get_coverage_across_genomes(all_files, taxid, output_dir, n_proc)
    cp = update_checkpoint(output_dir, "8_got_coverage")
    print("Completed check-point 8 got coverage across genomes")
else:
  print("Skipping 8_got_coverage because check-point wasn't 7_bowtie2_run")
  
# 9. Collate the output
if cp == "8_got_coverage":
  print("Collating all output")
  collate_output(all_files, taxid, output_dir, kreports, samples, group_samples, skip_bowtie2, skip_coverage)
  cp = update_checkpoint(output_dir, "9_collate_output")
  print("Completed check-point 9 collating output")
else:
  print("Skipping 9_collate_output because check-point wasn't 8_got_coverage")

# 10. Clean up all of the intermediate files
if cp == "9_collate_output":
  if skip_cleanup:
    cp = update_checkpoint(output_dir, "10_cleanup")
    print("Skipping check-point 10 cleanup of intermediate files")
  else:
    print("Cleaning up directory")
    clean_up(output_dir)
    cp = update_checkpoint(output_dir, "10_cleanup")
    print("Completed check-point 10 cleanup of intermediate files")
else:
  print("Skipping 10_cleanup because check-point wasn't 9_collate_output")
  
print("Finished running genome coverage checker pipeline.")
print("Running time: --- %s seconds ---" % round((time.time() - start_time), 2))
