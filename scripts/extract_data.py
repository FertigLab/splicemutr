#!/usr/bin/env python

# The purpose of this script is to build consensus of splicemutr rows that are immunogenic based on a patient's genotype
# although currently, only extract keys(tx row) and write to file

import pickle
from binascii import hexlify, unhexlify
from collections import defaultdict

class BytesIntEncoder:

    @staticmethod
    def encode(b: bytes) -> int:
        return int(hexlify(b), 16) if b != b'' else 0

    @staticmethod
    def decode(i: int) -> int:
        return unhexlify('%x' % i) if i != 0 else b''

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
            for i in range(num_peps):
                row_scores[i]=(data[i][1])/100
            SB_scores = [score for score in row_scores if score <= 50]
            WB_scores = [score for score in row_scores if score <= 500 and score > 50]
            SB_NUM = len(SB_scores)
            WB_NUM = len(WB_scores)
            if SB_NUM==0:
                SB_AV=0
            else:
                SB_AV=sum(SB_scores)/SB_NUM
            if WB_NUM==0:
                WB_AV=0
            else:
                WB_AV=sum(WB_scores)/WB_NUM
            pickle_txt.write("%s\t%s\t%s\t%s\t%s\n"%(row,str(SB_AV),str(SB_NUM),str(WB_AV),str(WB_NUM)))

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