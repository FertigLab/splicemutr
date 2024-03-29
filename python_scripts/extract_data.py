#!/usr/bin/env python

# The purpose of this script is to build consensus of splicemutr rows that are immunogenic based on a patient's genotype
# although currently, only extract keys(tx row) and write to file

import pickle
from binascii import hexlify, unhexlify
from collections import defaultdict
from splicemutr import BytesIntEncoder

def main(options):
    hla_allele = options.hla_allele
    hla_allele = hla_allele.replace(":","-")
    kmer_beg = int(options.kmer_beg)
    kmer_end = int(options.kmer_end)
    pickle_dir = options.pickle_dir
    pickle_template = "%s/%s_tx_dict_%d.pickle"
    pickle_all = "%s/%s_tx_dict_summary.txt"
    tx_dict_all = defaultdict(list)
    for kmer_length in range(kmer_beg,kmer_end):
        pickle_file = pickle_template%(pickle_dir,hla_allele,kmer_length)
        with open(pickle_file,'rb') as handle:
            pep_dict = pickle.load(handle)
        rows = list(pep_dict.keys())
        for row in rows:
            data = pep_dict[row]
            tx_dict_all[row] += data
    with open(pickle_all % (pickle_dir, hla_allele), "w") as pickle_txt:
        rows = list(tx_dict_all.keys())
        for row in rows:
            data = tx_dict_all[row]
            num_peps = len(data)
            row_scores=[None]*num_peps
            kmers=[None]*num_peps
            for i in range(num_peps):
                row_scores[i]=(data[i][1])/100
                kmers[i]=str(BytesIntEncoder.decode(data[i][0]).decode())
            scores = [str(score) for score in row_scores]
            #SB_scores = [str(score) for score in row_scores if score <= 50]
            #WB_scores = [str(score) for score in row_scores if score <= 500 and score > 50]
            pickle_txt.write("%s\t%s\t%s\n"%(str(row),":".join(kmers),":".join(scores)))

if __name__ == "__main__":

    from optparse import OptionParser

    parser = OptionParser()

    parser.add_option("-a", "--allele", dest="hla_allele",
                      help="the pickle HLA allele")
    parser.add_option("-p", "--pickle_dir", dest="pickle_dir",
                      help="the pickle directory")
    parser.add_option("-b", "--kmer_beg", dest="kmer_beg",
                      help="the lowest value in kmer range")
    parser.add_option("-e", "--kmer_end", dest="kmer_end",
                      help="the highest value in kmer range")

    (options, args) = parser.parse_args()

    main(options)
