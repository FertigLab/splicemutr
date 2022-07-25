#!/usr/bin/env python

# Created: 07/30/2021
# Created by: Theron Palmer

from Bio import SeqIO
import random

# kmerizes a single peptide given a single kmer length
def kmerize(pep,kmer_length):
    pep_len=len(pep)
    if len(pep) < kmer_length:
        pep += "Z"*(kmer_length-pep_len)
        pep_len=len(pep)
    return([pep[i:i+kmer_length] for i in range(pep_len-kmer_length+1)])

def main(options,args):
    num_kmers = options.num_kmers
    pep_file = args[0]
    output_file = args[1]
    
    pep_seqs = list(SeqIO.parse(pep_file, "fasta"))
    pep_set = set()
    for elem in pep_seqs:
        pep_set.update(kmerize(str(elem.seq),9))
    pep_list = list(pep_set)
    pep_samps = random.sample(pep_list,num_kmers)
    with open(output_file,"w") as out:
        print("\n".join(pep_samps),file=out)
    
    
if __name__ == "__main__":
    
    from optparse import OptionParser
    
    parser = OptionParser()
    parser.add_option("-k", "--num_kmers", dest="num_kmers",
                  help="The number of kmers to extract",
                  type=int)
    
    (options, args) = parser.parse_args()

    main(options,args)
    