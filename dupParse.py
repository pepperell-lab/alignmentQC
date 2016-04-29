#!/usr/bin/env python

import argparse
import sys
from Bio import SeqIO
import os
import getopt
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

#######################################################################
##This script takes as input the output of 'grep DUP *vcf' and 
##outputs two files: one that outputs the pos and number of isolates
##with a duplication at that site, the other outputs a mock binary 
##genome for each isolate with a 1 indicating a duplication was
##identified.
#######################################################################


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
        description="Collate summary info from multiple bamqc runs")
    parser.add_argument('-i', '--inputFile',
        help = 'file containing paths to bamqc results that are to be collated', 
        type = is_file,
        required = True)
    parser.add_argument('-o', '--output',
        help = 'name of output file,',
        required = True)
    return parser.parse_args()

def make_genomes():
    samples = {}
    with open(args.inputFile, 'r') as infile, open(args.output + '.simplified', 'w') as outfile:
        for line in infile:
            line = line.strip().split("\t")
            if len(line) > 1:
                #extract sample name
                samp = line[0].split(":")[0].split("/")[-1].split("_")[0]
                #print(samp)
                #get 'genome' of sample
                if samp in samples:
                    genome = samples[samp]
                else:
                    genome = [0]*4411532
                #add duplication to 'genome'
                start = int(line[1])
                stop = int(line[7].split(";")[2].split("=")[1])
                for i in range((start-1),stop):
                    genome[i] = 1
                samples[samp] = genome
                outfile.write('%s\t%i\t%i\n' %
                (samp,
                start,
                stop))
            else:
                continue
    return(samples)

def sum_genomes(samples):
    master = map(sum, zip(*samples.values()))
    with open(args.output + '.summed', 'w') as outfile:
        outfile.write("pos\tNo\n")
        for i, num in enumerate(master):
            pos = str(i + 1)
            outfile.write(str(i) + '\t' + str(num) + '\n')
    return(master)

def write_file(samples):
    with open(args.output + '.all', 'w') as outfile:
        pos = range(1, 4411533)
        outfile.write('samp\t' + "\t".join(str(x) for x in pos))
        for samp in samples:
            outfile.write('%s\t%s\n' % (
            samp,
            "\t".join(str(x) for x in samples[samp])))

args = get_arguments()
samples = make_genomes()
master = sum_genomes(samples)
write_file(samples)
