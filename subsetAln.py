#!/usr/bin/env python

import sys
from Bio import SeqIO
from Bio.Alphabet import IUPAC,Gapped
import os
import getopt
import argparse

#####################################################
##This script subsets an aln based on a user-defined
##threshold for missing data.
####################################################

def get_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description = "Subset an alignment")
    parser.add_argument('-a', '--alignment',
        help = 'fasta alignment', required = True)
    parser.add_argument('-m', '--missing',
        help = 'missing data file gen with missingness.py',
        required = True)
    parser.add_argument('-p', '--percent',
        help = 'percent missing threshold to exclude',
        required = True, type=float)
    parser.add_argument('-o', '--outfile',
        help = 'outfile name (file.fasta)',
        required = True)
    return parser.parse_args()

def define_passing():
    """make list of samples to include"""
    keep = []
    with open(args.missing, 'r') as infile:
        for line in infile:
            line = line.strip().split('\t')
            sample = line[0]
            per_missing = float(line[2])
            if per_missing < args.percent:
                keep.append(sample)
    return keep

def subset(keep):
    """for every strain in the keep list, output to new aln"""
    keepRecords = []
    for seq_record in SeqIO.parse(args.alignment, "fasta", 
    alphabet=Gapped(IUPAC.ambiguous_dna, '-')):
        if seq_record.id in keep:
            keepRecords.append(seq_record)
    SeqIO.write(keepRecords, args.outfile, "fasta")

args = get_arguments()
keep = define_passing()
subset(keep)
