#!/usr/bin/env python

import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser(description='This script converts a fastq file to a fasta file')
parser.add_argument('--fastq', dest='fastq',
                    help='fastq file name')
parser.add_argument('--fasta', dest='fasta',
                    help="fasta file name")
                    
args = parser.parse_args()
fastq = args.fastq
fasta = args.fasta

wf = SeqIO.convert(fastq, "fastq", fasta, "fasta")
