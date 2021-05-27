#!/usr/bin/env python

# The purpose of this script is to create a means of kmerizing and recording row locations (transcripts) that each unique kmer exists in

import sys
import os
from collections import defaultdict
import pickle
from binascii import hexlify, unhexlify
from random import choice
from string import ascii_uppercase

class BytesIntEncoder:

    @staticmethod
    def encode(b: bytes) -> int:
        return int(hexlify(b), 16) if b != b'' else 0

    @staticmethod
    def decode(i: int) -> int:
        return unhexlify('%x' % i) if i != 0 else b''

# kmerizes a single peptide given a single kmer length
def kmerize(pep,kmer_length):
    pep_len=len(pep)
    if len(pep) < kmer_length:
        pep += "Z"*(kmer_length-pep_len)
        pep_len=len(pep)
    return([pep[i:i+kmer_length] for i in range(pep_len-kmer_length+1)])

def generate_peptides(length,num_peps):
    peps=[None]*num_peps
    for i in range(num_peps):
        peps[i]=(''.join(choice(ascii_uppercase) for i in range(length)))
    return(peps)

def main(pep_dat,kmer_length):
    out=os.getcwd()
    pep_dict = defaultdict(list)
    sys.stdout.write("kmer_length:%d\n"%kmer_length) # printing kmer length to standard out
    for i in range(len(pep_dat)):
        kmers=kmerize(pep_dat[i],kmer_length)
        for kmer in kmers:
            pep_dict[BytesIntEncoder.encode(kmer.encode())].append(i)
    with open(out+'/pep_dict_'+str(kmer_length)+'.pickle', 'wb') as pickle_file:
        pickle.dump(pep_dict, pickle_file)

if __name__ == "__main__":
    main_peptides=generate_peptides(500,1000)
    #for i in range(9,26):
    #  main(main_peptides,i)

    main_peps_log=[False]*len(main_peptides)
    for i in range(len(main_peptides)):
        if (BytesIntEncoder.decode(BytesIntEncoder.encode(main_peptides[i].encode())).decode() == main_peptides[i]):
            main_peps_log[i]=True
    if (all(main_peps_log)==True):
        print(True)