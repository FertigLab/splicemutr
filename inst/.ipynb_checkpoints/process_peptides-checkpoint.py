#!/usr/bin/env python

# The purpose of this script is to create a means of kmerizing and recording row locations (transcripts) that each unique kmer exists in

import sys
import os
from collections import defaultdict
import pickle
from binascii import hexlify, unhexlify
from splicemutr import BytesIntEncoder

# kmerizes a single peptide given a single kmer length
def kmerize(pep,kmer_length):
    pep_len=len(pep)
    if len(pep) < kmer_length:
        pep += "Z"*(kmer_length-pep_len)
        pep_len=len(pep)
    return([pep[i:i+kmer_length] for i in range(pep_len-kmer_length+1)])


def main(options):
    pep_dat_file = options.peptides
    out=options.out
    kmer_length=int(options.kmer_length)
    with open(pep_dat_file, "r") as peps:
        pep_dat=peps.read().replace("*","").splitlines()
    pep_dict = defaultdict(list)
    sys.stdout.write("kmer_length:%d\n"%kmer_length) # printing kmer length to standard out
    for i in range(len(pep_dat)):
        kmers=kmerize(pep_dat[i],kmer_length)
        for kmer in kmers:
            pep_dict[BytesIntEncoder.encode(kmer.encode())].append(i)
    with open(out+'/pep_dict_'+str(kmer_length)+'.pickle', 'wb') as pickle_file:
        pickle.dump(pep_dict, pickle_file)
    with open(out+'/peps_'+str(kmer_length)+'.txt','w') as peps:
        kmer_list = [BytesIntEncoder.decode(kmer).decode() for kmer in list(pep_dict.keys())]
        kmer_string = "\n".join(kmer_list)
        print(kmer_string,file=peps)

if __name__ == "__main__":

    from optparse import OptionParser

    parser = OptionParser()

    parser.add_option("-p", "--peptides", dest="peptides",
                      help="the peptide data file")
    parser.add_option("-o", "--out", dest="out",
                      help="out directory")
    parser.add_option("-k", "--kmer_length", dest="kmer_length",
                      help="the kmer length")
    (options, args) = parser.parse_args()

    main(options)