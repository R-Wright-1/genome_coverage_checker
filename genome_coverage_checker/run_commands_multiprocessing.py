#!/usr/bin/env python

import os
from multiprocessing import Pool, freeze_support
import argparse

parser = argparse.ArgumentParser(description='This script is to run multiple commands - it is designed to work with multiple wget commands for downloading genomes, but would probably work fine with other commands, too')
parser.add_argument('--commands', dest='commands',
                    help="The list of commands to be run in a file. This should be the file name.")
parser.add_argument('--processors', dest='n_proc', default=1,
                    help="Number of processors to use for all multiprocessing steps (downloading genomes, pulling out reads, running QUAST)")

args = parser.parse_args()
commands = args.commands
n_proc = int(args.n_proc)

command_list = []
for row in open(commands, 'r'):
  command_list.append(row.replace('\n', ''))

def run_command(command):
  os.system(command)
  return

def run_multiprocessing(func, i, n_processors):
    with Pool(processes=n_processors) as pool:
        return pool.map(func, i)
      
def main():
  run_multiprocessing(run_command, command_list, n_proc)

if __name__ == "__main__":
    freeze_support()   # required to use multiprocessing
    main()
