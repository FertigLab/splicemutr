#!/usr/bin/env python
# python python runMHCnuggets.py -t "I" -k /media/theron/My_Passport/data/references/kmers_15 -m
# /media/theron/My_Passport/splicemutr/MHC_I_alleles.txt -o /media/theron/My_Passport/data/references/scores

import sys
import sklearn
import tensorflow
from mhcnuggets.src.predict import predict
import os
import argparse

parser = argparse.ArgumentParser(description='Inputting the sequences')
parser.add_argument('-t', '--type', action="store", nargs=1, help='The MHC class: "I" or "II"', dest='type')
parser.add_argument('-k', '--kmers', action="store", nargs=1, help='The input kmer directory path', dest='kmers')
parser.add_argument("-m", "--mhc", action="store", nargs=1, help='The MHC allele',
                    dest='mhc')
parser.add_argument("-o", "--output", action="store", nargs=1, help="The output directory to write scores to",
                    dest='output')
results = parser.parse_args()

#kmers = os.listdir(results.kmers[0])  # a list of kmer_files with a fraction of the total unique number of kmers
#for file in kmers:  # iterating through each individual text file and running through mhcnuggets
#s = results.kmers[0] + "{}" + "{}"  # preparing for formatting of the path
#s = s.format("/", file)
#s_split = s.split("/")  # parsing the path
#kmer_file = s_split[len(s_split) - 1].split(".")[0]  # extracting just mhc file name
s = results.kmers[0]
s_split = s.split("/")  # parsing the path
kmer_file = s_split[len(s_split) - 1].split(".")[0]  # extracting just mhc file name
mhc=results.mhc[0]
with open(results.mhc[0], 'r') as f: # iterating through the mhc alleles chosen to run mhcnuggets with
    mhc = f.readline()
    while mhc != '':
        mhc = mhc.rstrip()
        mhc_print = mhc.replace(":", "-") # some operating systems do not recognize ":" in file names, so replace.
        out_file = results.output[0] + "{}" + "{}" + "{}" + "{}" + "{}"
        out_file = out_file.format("/", mhc_print, "_", kmer_file, ".txt")
        predict(class_=results.type[0], peptides_path=s, mhc=mhc, output=out_file)
        mhc = f.readline()
