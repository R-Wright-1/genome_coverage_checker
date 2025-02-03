#!/usr/bin/env python

import os
import pandas as pd
import sys
#from multiprocessing import Pool
#from multiprocessing import freeze_support
#import numpy as np
import argparse
#from Bio.SeqRecord import SeqRecord
#from Bio import SeqIO
import pickle
#from function_file_modified import *
from genome_coverage_checker.util import *

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
# parser.add_argument('--quast_loc', dest='quast_loc', default=None,
#                     help="The location of the executable quast.py script (folder containing QUAST)")
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
parser.add_argument('--rerun', dest='rerun', default=False, action='store_true',
                    help="If this is set to True, it will re-run the extraction of reads from samples regardless of whether the files already exist")
parser.add_argument('--all_domains', dest='all_domains', default=False, action='store_true',
                    help="The default for coverage checker is for only the genomes of prokaryotes to be downloaded and checked. If you'd like to include all genomes then add this flag.")
parser.add_argument('--representative_only', dest='representative_only', default=False, action='store_true',
                    help="The default for coverage checker is to use all genomes. If you'd like to limit the checking to only representative and reference genomes then add this flag.")
parser.add_argument('--skip_bowtie2', dest='skip_bowtie2', default=False, action='store_true',
                    help="Whether to skip bowtie2 mapping. Default is to run bowtie2, but bowtie2 can be difficult to install on a mac so add this flag if you don't want to install bowtie2.")
parser.add_argument('--skip_coverage', dest='skip_coverage', default=False, action='store_true',
                    help="Whether to skip getting coverage for all reads across the genomes. This can take a little while, so you can skip it if you don't think that this output is useful to you.")
parser.add_argument('--skip_cleanup', dest='skip_cleanup', default=False, action='store_true',
                    help="If you want to skip the cleanup. This will keep the intermediate files containing some of the commands run, e.g. for making bowtie2 databases. They may be helpful if you're trying to troubleshoot issues.")

#Read in the command line arguments
args = parser.parse_args()
sample_name, n_proc, sample_dir, fastq_dir, kraken_kreport_dir, kraken_outraw_dir, output_dir = args.sample_name, int(args.n_proc), args.sample_dir, args.fastq_dir, args.kraken_kreport_dir, args.kraken_outraw_dir, args.output_dir
if sample_dir != None:
  fastq_dir, kraken_kreport_dir, kraken_outraw_dir, output_dir = sample_dir, sample_dir, sample_dir, sample_dir
read_lim, read_mean, assembly_folder = args.read_lim, int(args.read_mean), args.assembly_folder
if assembly_folder == None:
  assembly_folder = output_dir
sample_metadata, species, project_name, rerun, all_domains, representative_only, skip_bowtie2, skip_coverage, skip_cleanup = args.sample_metadata, args.species, args.project_name, args.rerun, args.all_domains, args.representative_only, args.skip_bowtie2, args.skip_coverage, args.skip_cleanup
wd = os.getcwd()
if read_lim == None: read_lim = 0
else: read_lim = int(read_lim)

#check whether we've already run this and which checkpoint we're at
if rerun:
  cp = '0'
elif os.path.exists(output_dir+'checkpoint.txt'):
  cp = get_checkpoint(output_dir)
else:
  cp = '0'

#if not rerun and cp not 0, check whether the arguments given are the same?

# 1. Run the initial checks
wd, n_proc, sample_dir, sample_name, fastq_dir, kraken_kreport_dir, kraken_outraw_dir, output_dir, assembly_folder, read_lim, read_mean, sample_metadata, species, project_name, rerun, md, samples, taxid_name = run_initial_checks(wd, n_proc, sample_dir, sample_name, fastq_dir, kraken_kreport_dir, kraken_outraw_dir, output_dir, assembly_folder, read_lim, read_mean, sample_metadata, species, project_name, rerun) #run all of the initial checks to ensure that all of the folders and files exist before starting to try and run anything
cp = update_checkpoint(output_dir, "1_initial_checks_run")

# 2. Read in all kreports and combine them to get the genomes that we'll need
group_samples, taxid, kreports = get_kreports(samples, kraken_kreport_dir, output_dir, md, read_lim, read_mean, project_name)
cp = update_checkpoint(output_dir, "2_combined_kreports")

# 3. Download all genomes
taxid = download_genomes(taxid, assembly_folder, output_dir, all_domains, representative_only, n_proc)
cp = update_checkpoint(output_dir, "3_downloaded_genomes")

# 4. Extract the reads for each taxonomy ID
extract_reads(taxid, output_dir, samples, fastq_dir, kraken_kreport_dir, kraken_outraw_dir, n_proc)
cp = update_checkpoint(output_dir, "4_extracted_reads")

# 5. Combine the files for each taxonomy ID (across multiple samples for the sample groupings)
all_files = combine_convert_files(taxid, output_dir, samples, group_samples, n_proc)
cp = update_checkpoint(output_dir, "5_combined_files")

# 6. Run QUAST
run_quast(all_files, taxid, output_dir, n_proc)
cp = update_checkpoint(output_dir, "6_quast_run")

# 7. Make bowtie2 databases & run bowtie2 for all taxonomy ID's
if skip_bowtie2:
  cp = update_checkpoint(output_dir, "7_bowtie2_run")
else:
  make_bowtie2_databases(taxid, output_dir, n_proc)
  run_bowtie2(all_files, taxid, output_dir, n_proc)
  cp = update_checkpoint(output_dir, "7_bowtie2_run")

# 8. Get the coverage and mapping of reads across the genomes
if skip_coverage:
  cp = update_checkpoint(output_dir, "8_got_coverage")
else:
  get_coverage_across_genomes(all_files, taxid, output_dir, n_proc)
  cp = update_checkpoint(output_dir, "8_got_coverage")
  
# 9. Collate the output
collate_output(all_files, taxid, output_dir, kreports, samples, group_samples, skip_bowtie2, skip_coverage)
cp = update_checkpoint(output_dir, "9_collate_output")

# 10. Clean up all of the intermediate files
if skip_cleanup:
  cp = update_checkpoint(output_dir, "10_cleanup")
else:
  clean_up(output_dir)
  cp = update_checkpoint(output_dir, "10_cleanup")
