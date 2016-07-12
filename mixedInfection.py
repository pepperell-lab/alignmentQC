#!/usr/bin/env python

import argparse
import sys
from Bio import SeqIO
import os
import getopt
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

####################################################################
##This script parses the pilon VCF for AF's to determine if there   
##is evidence for a mixed infection.
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
    input_file.add_argument('-v', '--vcf',
        help ='pilon VCF', 
        type=is_file)
    return parser.parse_args()

args = get_arguments()

def get_AFs():
    """extract AF of variant sites"""
    AFs = []
    with open(args.vcf, "r") as vcf:
        for line in vcf:
            if line[0] != "#":
                line = line.strip().split("\t")
                info = line[7].split(";")
                for i in info:
                    if "AF" in i:
                        freq = i.split("=")[-1]
                        if float(freq) != 0.0:
                            AFs.append(float(freq))
    return AFs

def make_plot(AFs):
    """make a plot of AF"""
    plotFileName = args.vcf + "_af.png"
    x = np.array(AFs)
    bins = np.arange(0.0, 1.0, 20) # fixed bin size
    plt.xlim([0, 1])
    plt.hist(x, bins=bins, aplpha=0.75)
    plt.xlabel("Alt Allele Frequency")
    plt.ylabel("Number of Sites")
    plt.title("Histogram of AF")
    plt.savefig(plotFileName)
    plt.close()

AFs = get_AFs()
make_plot(AFs)
