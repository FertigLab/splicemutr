#!/usr/bin/env python

# The purpose of this script is to determine the binders associated with each splicemutr data row

import sys
import os
from collections import defaultdict
import pickle
from binascii import hexlify, unhexlify



# read in the peptide data file except for the first line
# split the peptide data file by newlines
# if there is a peptide(you see a stop codon as the second to last character), kmerize the peptide,otherwise continue
# populate a dictionary with the kmers and append the row number to a list in the kmer's value
# save this dictionary as a pickle element

class BytesIntEncoder:

    @staticmethod
    def encode(b: bytes) -> int:
        return int(hexlify(b), 16) if b != b'' else 0

    @staticmethod
    def decode(i: int) -> int:
        return unhexlify('%x' % i) if i != 0 else b''

def main(options):
    binders_file=options.binders
    pickle_dir=options.pickle_dir
    out=options.out
    kmer_length=int(options.kmer_length)
    sys.stdout.write("kmer_length:%d\n"%kmer_length) # printing kmer length to standard out

    with open(binders_file, "r") as binders:
        binders_dat=binders.read().splitlines()
    tx_dict = defaultdict(list)
    with open(pickle_dir+'/pep_dict_'+str(kmer_length)+'.pickle', 'rb') as handle:
        pep_dict = pickle.load(handle)

    for i in range(len(binders_dat)):
        kmer_score=binders_dat[i].split(",")
        kmer=kmer_score[0]
        if len(kmer) != kmer_length: continue
        kmer=BytesIntEncoder.encode(kmer.encode())
        score=kmer_score[1]
        for row in pep_dict[kmer]:
            tx_dict[row].append((kmer,float(score)))

    HLA_kmers=os.path.basename(binders_file)
    HLA=HLA_kmers.split("_")[0]
    with open(out+'/'+HLA+'_tx_dict_'+str(kmer_length)+'.pickle', 'wb') as pickle_file:
        pickle.dump(tx_dict, pickle_file)

if __name__ == "__main__":

    from optparse import OptionParser

    parser = OptionParser()

    parser.add_option("-b", "--binders", dest="binders",
                      help="the binder files")
    parser.add_option("-p", "--pickle_dir", dest="pickle_dir",
                      help="the pickle directory")
    parser.add_option("-o", "--out", dest="out",
                      help="out directory")
    parser.add_option("-k", "--kmer_length", dest="kmer_length",
                      help="the kmer length")
    (options, args) = parser.parse_args()

    main(options)