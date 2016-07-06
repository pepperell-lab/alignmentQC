# alignmentQC
QC scripts for alignments

###missingness.py
This script takes a fasta alignment as input and calculates the number
of gaps/missing sites (denoted by "-") for each strain in the alignment
and outputs this info in a file and a histogram.

Usage: missingness.py -i [inputfile]

###gapPositions.py
This script takes a fasta alignment as input and calculates the number of
gaps at each site in the alignment and outputs this to a file.

Usage: gapPositions.py [alignment]
