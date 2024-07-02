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

## Example commands used

```
conda activate coverage-checker
cd /home/robyn/kraken_coverage
mkdir ben_coverage_test
cp -r /bigpool/benfish/GDC/coverage_test/* ben_coverage_test
rm -r ben_coverage_test/output_dir
mkdir ben_coverage_test/output_dir

python run_coverage_checker_modified.py \
  --processors 24  \
  --fastq_dir /home/robyn/kraken_coverage/ben_coverage_test/fastq_dir/ \
  --kraken_kreport_dir /home/robyn/kraken_coverage/ben_coverage_test/kraken_kreport_dir/ \
  --kraken_outraw_dir /home/robyn/kraken_coverage/ben_coverage_test/kraken_outraw_dir/ \
  --output_dir /home/robyn/kraken_coverage/ben_coverage_test/output_dir/ \
  --quast_loc /home/robyn/tools/quast-5.2.0/ \
  --sample_metadata /home/robyn/kraken_coverage/ben_coverage_test/metadata.csv \
  --project_name ben-test \
  --read_lim 100
```
