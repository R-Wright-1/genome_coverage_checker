# Genome Coverage Checker

This is a repository for storing the code needed to check the genome coverage for taxa identified by Kraken 2 in metagenome samples. It is currently still a work in progress and the code needs streamlining, checkpoints added, etc., but I welcome any feedback on it. 

This came from the issue that many users have mentioned regarding Kraken 2 identifying many false positive taxa within metagenome samples. While we have found [previously](https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000949) that the parameters used to run Kraken 2 can drastically reduce false positive taxa (i.e., by increasing the confidence threshold), it is not always practical to do this. In environmental or low biomass samples this typically leads to having few or no reads classified. We therefore set out to verify the taxa present in metagenome samples by doing a few things:
1. Taking the taxonomy ID's for taxa identified by Kraken 2 as being present in metagenome samples
2. Download reference genomes for these taxonomy ID's
3. Extract reads mapping to these taxonomy ID's from fastq files (with KrakenTools)
4. Map reads to the reference genome (QUAST/Bowtie2)
5. Get information on the genome fraction present, percent identity of the reads to the reference genome and the actual:expected ratio (see below)

Actual:expected ratio: Comparison of the genome fraction present (actual) with the expected fraction based on the reference genome size and the number of reads assigned by Kraken 2 x length of the reads.

## Packages required

I have installed everything needed into a conda environment on the servers that I use called ```coverage-checker```. Using these commands:
```
conda activate coverage-checker
conda list --explicit > spec-file.txt
```
I have created ```spec-file.txt``` which you can find in this repository. As I streamline the code used, I'll work on creating a better list of the packages required to run this.

Additionally, I have a folder containing the QUAST executables. 

## Example commands used

In order to run, these scripts need:
- a folder containing the kraken kreports
- a folder containing the raw kraken outputs
- a folder containing fastq files
- location of QUAST executables

Examples of these files can be found in the example_data folder (note that all files within these would need to be unzipped prior to running).

To run:
```
mkdir output_dir
python run_coverage_checker_modified.py \
  --processors 24  \
  --fastq_dir fastq_dir/ \
  --kraken_kreport_dir kraken_kreport_dir/ \
  --kraken_outraw_dir kraken_outraw_dir/ \
  --output_dir output_dir/ \
  --quast_loc /home/robyn/tools/quast-5.2.0/ \
  --sample_metadata metadata.csv \
  --project_name test-coverage \
  --read_lim 100
```

Please also note that the project and sample names should have no underscores, but hyphens are fine. Currently I've set up the naming format to use underscores for the additional files created and splits are done on these underscores, so using them in sample names will confuse the script. I plan to fix this in the future. 

Please see the help from the script for more information. This can be run with or without the ```metadata.csv``` file - this is designed to run the checker multiple times for different sample groupings as we felt it may be useful to see the genome coverage across all samples or within a metadata grouping for a taxon.

Note that this is currently designed for bacteria and archaea only. 

### What the script is doing

1. Read in command line arguments
2. Run initial checks to ensure all files and folders exist
3. Read in kraken kreports and get the taxa with >= read_lim
4. Get the sample groupings from the metadata
5. Make lists of Taxonomy IDs and files that we'll use/create
6. Add genomes that we'll download to a dictionary
7. Make lists of files that will be combined
8. Download genomes (using multiprocessing)
9. Extract reads for each taxonomy ID (using multiprocessing)
10. Check that files were created
11. Run QUAST for all taxonomy ID's (using multiprocessing)
12. Check whether QUAST output files were made for all taxonomy IDs
13. Make Bowtie2 databases and run Bowtie2 for all taxonomy IDs (using multiprocessing)
14. Get coverage and mapping of reads across genomes for each sample (using multiprocessing)
15. Summarise results

## Output

The expected output within ```output_dir``` is:
- ```assembly_summary_archaea.txt```: NCBI RefSeq assembly summary for archaea
- ```assembly_summary_bacteria.txt```: NCBI RefSeq assembly summary for bacteria
- ```project-name_combined_kreport.csv```: A table showing the number of reads assigned to taxa within all samples
- ```project-name_genomes_not_downloaded.txt```: Text file containing the names of any genomes that couldn't be downloaded. This is typically due to unexpected punctuation characters within the species names that confuse the commands
- ``project-name_output.txt```: Main output file for coverage-checker. See below for details
- ```project-name_quast_not_run.txt```: Files for which QUAST couldn't be run. This is typically because the file was empty
- ```project-name_unmade_files.txt```: Files that couldn't be made - this is likely because we were trying to create files for every taxonomy ID within every sample, but not every taxonomy ID was found in every sample
- ```bowtie2_db```: folder containing bowtie2 databases for all genomes
- ```bowtie2_mapped```: folder containing reads mapped by bowtie2 in fasta and sam format
- ```downloaded_genomes.dict```: python pickle dictionary containing all genomes that have been downloaded
- ```genomes```: folder containing all downloaded reference genomes
- ```genomes_not_working.txt```
- ```pickle_coverage```: folder containing pickle objects of coverage within each sample/sample group
- ```QUAST```: folder containing all QUAST output folders
- ```reads_mapped```: extracted reads for each taxonomy ID within each sample in fasta (fa) or fastq (fq) format
- ```taxid_name.dict```: python pickle dictionary containing all taxonomy ID's and names

### ***output.txt

This file contains most of the information that we are interested in summarised for each taxonomy ID. The columns are:
- NCBI taxonomy ID
- NCBI species name
- Sample: sample that this row shows information on (note that this could be the same as sample group if it shows the summary for the whole group)
- Sample group: metadata group that this sample belongs to
- Reads Kraken assigned: number of reads that kraken has assigned *at* the species level
- Reads Kraken assigned below species: number of reads that kraken assigned to a strain (or other lower rank) within this species. We have included all of these within this
- Reference length (bp): length of the NCBI RefSeq reference genome in bp
- Reference GC (%): GC percentage of the reference genome
- QUAST contigs: number of reads that QUAST mapped to the reference genome (note that these are only called contigs because this is the default for QUAST to run with. No assembly has been done here so they are not contigs)
- QUAST unaligned contigs: number of reads that QUAST did not map to the reference genome
- GC (%): GC percentage within the reads (calculated by QUAST)
- Genome fraction (%): genome fraction present within the reads (calculated by QUAST)
- Duplication ratio: duplication ratio within the reads (calculated by QUAST)
- Total length (number of reads x read length; bp): total length of all reads combined (note that this is based on all reads being 100bp and doesn't currently check the actual read length)
- Number of bases/reference genome length (expected genome fraction; %): fraction of the reference genome that we expect to be covered based on the total length of the reads identified by kraken as belonging to this taxonomy ID. Note that this is very simplistic as we obviously know that we wouldn't have totally equal coverage across an entire genome. 
- Actual/expected fraction: actual genome fraction present divided by expected genome fraction. The closer to 1 that this is, the better.
- Expected/actual fraction: the inverse of above - this was calculated at some point as I thought it might be interesting, but I'll probably take it out at some point
- Bowtie2 mapped reads: the number of reads that Bowtie2 mapped to the reference genome
- Genome median identity: the median identity (%) of all mapped reads to the reference genome
- Genome proportion with mapped reads: this is the proportion of the genome (in 500 bp increments) that has at least one read mapped to it. We thought that this might give an indication of whether all of the reads are mapping to the same area or not, but probably needs some refinement and testing to work out what is reasonable and whether it is actually useful
- Proportion of reads mapped in distinct genomic areas (round 500): this is looking at how many of the reads are mapped to distinct areas of the genome (in 500 bp increments)
- Proportion of reads mapped in distinct genomic areas (round 1000): this is looking at how many of the reads are mapped to distinct areas of the genome (in 1000 bp increments)
- Proportion of reads mapped in distinct genomic areas (round 5000): this is looking at how many of the reads are mapped to distinct areas of the genome (in 5000 bp increments)

## Plotting the output

This script can be run like so:
```
python plot_coverage.py --sample_name P15O --processors 24 --sample_dir /home/robyn/kraken_coverage/coverage_output/ --num_tax 30
```

It is for plotting the coverage across the genome of the reads for an individual sample. I will add more information on documenting this later. 

