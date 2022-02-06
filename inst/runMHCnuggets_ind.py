#!/usr/bin/env python

# python python runMHCnuggets.py -t "I" -k /media/theron/My_Passport/data/references/kmers_15 -m
# /media/theron/My_Passport/splicemutr/MHC_I_alleles.txt -o /media/theron/My_Passport/data/references/scores

import sys
import sklearn
import tensorflow
from mhcnuggets.src.predict import predict
import os

def main(options):
    s = options.kmers
    s_split = s.split("/")  # parsing the path
    kmer_file = s_split[len(s_split) - 1].split(".")[0]  # extracting just mhc file name
    mhc = options.mhc
    mhc = mhc.rstrip()
    mhc_print = mhc.replace(":", "-") # some operating systems do not recognize ":" in file names, so replace.
    out_file = options.output + "{}" + "{}" + "{}" + "{}" + "{}"
    out_file = out_file.format("/", mhc_print, "_", kmer_file, ".txt")
    predict(class_=options.type, peptides_path=s, mhc=mhc, output=out_file)

if __name__ == "__main__":

    from optparse import OptionParser

    parser = OptionParser()

    parser.add_option("-t", "--type",
                      help='The MHC class: "I" or "II"',
                      dest='type')
    parser.add_option("-k", "--kmers",
                      help='The input kmer directory path',
                      dest='kmers')
    parser.add_option("-m", "--mhc",
                      help='The MHC allele file',
                      dest='mhc')
    parser.add_option("-o", "--output",
                      help="The output directory to write scores to",
                      dest='output')
    (options, args) = parser.parse_args()

    main(options)