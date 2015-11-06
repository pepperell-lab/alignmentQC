#!/usr/bin/env python

import argparse
import sys
from Bio import SeqIO
import os
import getopt
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

####################################################################
##This script takes an alignment (fasta)  and generates a histogram  
##and a file with the amount of missing sites (denoted by "-").
####################################################################


class FullPaths(argparse.Action):
    """Expand user- and relative-paths"""
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest,
            os.path.abspath(os.path.expanduser(values)))

def is_file(filename):
    """Checks if a file exists"""
    if not os.path.isfile(filename):
        msg = "{0} is not a file".format(filename)
        raise argparse.ArgumentTypeError(msg)
    else:
        return filename

def get_arguments(): 
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description="Quantify missing data per strain from an alignment")
    input_file = parser.add_mutually_exclusive_group(required=True)
    input_file.add_argument('-i', '--inputFile',
        help ='alignment in fasta format', 
        type=is_file)
    return parser.parse_args()

args = get_arguments()
alnIN = args.inputFile

def quantify_missingness(alnIN):
    d = {}
    """for every strain in the aln, count the number of missing data sites"""
    for seq_record in SeqIO.parse(alnIN, "fasta"):
        print("Processing sample {0}".format(seq_record.id))
        no_gaps = seq_record.seq.count("-")
        print("{0} has {1} gaps".format(seq_record.id, no_gaps))
        per_missing = float(no_gaps)/len(seq_record.seq) * 100
        d[seq_record.id] = [no_gaps, per_missing]
    return d

def make_plot(d, alnIN):
    """make a plot of the missing data"""
    plotFileName = alnIN.split(".")[0] + "missingData.png"
    x = []
    for key in d:
        x.append(d[key][1])
    n, bins, patches = plt.hist(x, 100, normed=1, alpha=0.75)
    plt.xlabel('Percent Missing')
    plt.ylabel('No. of Isolates')
    plt.title('Histogram of Missing Data')
    plt.savefig(plotFileName)
    plt.close()
            
def write_file(d, alnIN):
    """write the contents of dictionary to file"""
    outFileName = alnIN.split(".")[0] + "_missingData.txt"
    with open(outFileName, 'w') as outfile:
        for samp in d:
            outfile.write('%s\t%i\t%f\n' %
                (samp,
                d[samp][0],
                d[samp][1]))

d = quantify_missingness(alnIN)
make_plot(d, alnIN)
write_file(d, alnIN)
