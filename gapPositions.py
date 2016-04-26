#!/usr/bin/env python

import sys
import os
import argparse
from Bio import AlignIO

# This script filters protein groups output from OrthoMCL and compares the core
# genomes of groups input by the user

def is_file(filename):
    """Checks if a file exists"""
    if not os.path.isfile(filename):
        msg = "{0} is not a file".format(filename)
        raise argparse.ArgumentTypeError(msg)
    else:
        return filename

def get_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='Count gaps at alignment sites')
    parser.add_argument("align", help="Fasta alignment", type=is_file)
    return parser.parse_args()

def find(s, ch):
    return [i for i, ltr in enumerate(s) if ltr == ch]

def find_gaps(align):
    alignment = AlignIO.read(align, "fasta")
    length = alignment.get_alignment_length()
    gap = [0]*length
    for record in alignment:
        gaps = find(record.seq, '-')
        for g in gaps:
            gap[g] += 1
    return gap

def write_gaps(align, gap):
    outfile_name = align.split(".")[0] + "_gaps.txt"
    with open(outfile_name, "w") as outfile:
        for i,j in enumerate(gap):
            outfile.write("{0}\t{1}\n".format(i,j))
        
    
args = get_args()

gap = find_gaps(args.align)
write_gaps(args.align, gap)
