#!/usr/bin/env python

import pandas as pd
from sys import exit
import os
import splicemutr as sp

def main(options, args):
    kmers_files = options.kmers_files
    output_dir = options.output_dir
    
    if not os.path.isdir(output_dir):
        try:
            os.mkdir(output_dir)
        except:
            sys.exit("output directory does not exist and cannot be created")

    with open(kmers_files) as kmers_f:
        all_kmers_files = kmers_f.read().splitlines()
        all_kmers_counts = [[] for i in range(len(all_kmers_files))]
        sample_names = [os.path.basename(file).replace("_splicemutr_kmers.txt","") for file in all_kmers_files]
        for i in range(len(all_kmers_files)):
            kmer_file=all_kmers_files[i]
            try:
                kmer_file_df = pd.read_table(kmer_file)
                kmer_file_df.fillna('', inplace=True)
            except:
                print("failed to open file: %s"%kmer_file)
                continue
            kmers = kmer_file_df.kmers.tolist()            
            all_kmers_counts[i] = [len([k_filt for k_filt in k.split(":") if k_filt != ""]) for k in kmers]
        all_kmers_counts = pd.DataFrame(all_kmers_counts).transpose()
        all_kmers_counts.columns = sample_names
        all_kmers_counts.to_csv("%s/all_kmers.txt"%(output_dir),sep='\t',index=False)
    
if __name__ == "__main__":

    from optparse import OptionParser

    parser = OptionParser()
    
    parser.add_option("-k","--kmers_files",dest="kmers_files",
                     help="file containing one kmer file per line")
    parser.add_option("-o", "--output_dir", dest="output_dir",
                  help="the output directory")
    (options, args) = parser.parse_args()

    main(options, args)