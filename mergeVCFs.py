#!/usr/bin/env python

import argparse
import sys
from Bio import SeqIO
import os
import getopt
import glob
import subprocess
import shlex

####################################################################
##This script take as input a file with the paths to all VCFs you 
##would like joined, or a directory path and will join all VCFs in
##the specified directory. 
####################################################################


class FullPaths(argparse.Action):
    """Expand user- and relative-paths"""
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest,
            os.path.abspath(os.path.expanduser(values)))

def listdir_fullpath(d):
    return [os.path.join(d, f) for f in os.listdir(d)]

def is_dir(dirname):
    """Checks if a path is a directory"""
    if not os.path.isdir(dirname):
        msg = "{0} is not a directory".format(dirname)
        raise argparse.ArgumentTypeError(msg)
    else:
        return dirname
        
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
        description="Merge individual sample VCFs into a project VCF")
    input_files = parser.add_mutually_exclusive_group(required=True)
    input_files.add_argument('-i', '--inputFile',
        help = 'list of paths to VCFs to include', 
        type = is_file)
    input_files.add_argument('-d', '--directory', 
        help = 'directory containing VCFs to be merged',
        type = is_dir)
    parser.add_argument('-r', '--reference', required=True,
        help = 'fasta reference',
        type = is_file)
    parser.add_argument('-p', '--project', default='project')
    parser.add_argument('-t', '--threads', default=1)
    return parser.parse_args()

def combine_variants(inputList):
    """Merge VCFs with GATK CombineVariants"""
    print("Calling GATK CombineVariants")
    cmd = ("java -Xmx4g -jar /opt/PepPrograms/GenomeAnalysisTK.jar \
-T CombineVariants \
-R {0} ".format(args.reference) +
"-nt {0} --variant ".format(args.threads) +
" --variant ".join(inputList) +
" -o {0}.vcf ".format(args.project) +
"-genotypeMergeOptions UNIQUIFY \
--excludeNonVariants")
    subprocess.call(shlex.split(cmd))

def make_inputList_file():
    """Make inputList from file"""
    inputList = []
    with open(args.inputFile, 'r') as infile:
        for line in infile:
            line = line.strip()
            inputList.append(line)
    print("There are {0} VCFs to be merged.".format(len(inputList)))
    return inputList

def make_inputList_directory():
    """Make inputList from directory"""
    inputList = glob.glob(args.directory + "/*.vcf")
    print("There are {0} VCFs to be merged.".format(len(inputList)))
    return inputList
    
args = get_arguments()
if args.inputFile is None:
    inputList = make_inputList_directory()
else:
    inputList = make_inputList_file()

combine_variants(inputList)
