#!/usr/bin/env python

# author: Theron Palmer
# created: 10/03/2021

# This script creates a per sample immunogenic kmers file per TCGA cancer analyzed

import pandas as pd
import os

def get_summary(tumor_kmer_file):
    tumor_kmer_frame = pd.read_table(tumor_kmer_file,sep="\t") # this part should be done in the function as well
    tumor_kmer_frame.columns = [str(i) for i in list(range(len(tumor_kmer_frame.columns)))]
    ncols = len(tumor_kmer_frame.columns)
    count = []
    for i in range(ncols-2):
        col_vals = tumor_kmer_frame[str(i)].tolist()
        col_vals = [val for val in col_vals if pd.isnull(val) == False]
        vals = [len(val.split(":")) for val in col_vals]
        count.append(sum(vals))
    return(count)

def main(options,args):
    tumor_dir=options.tumor_dir
    cancer=options.cancer
    kmer_file = "%s/%s/kmer_files.txt"%(tumor_dir,cancer)
    with open(kmer_file) as kmer:
        kmer_files = kmer.read().splitlines()
    for i in range(len(kmer_files)):
        kmer_filt_file = "%s/%s/%s_kmers_%d_filt.txt"%(tumor_dir,cancer,cancer,i+1)
        tumor_kmer_summary = get_summary(kmer_filt_file)
        tumor_kmer_summary = "\n".join([str(val) for val in tumor_kmer_summary])
        tumor_kmer_summary_file = "%s/%s/%s_kmers_summary_%d.txt"%(tumor_dir,cancer,cancer,i+1)
        
        with open(tumor_kmer_summary_file,"w") as summ:
            print(tumor_kmer_summary,file=summ)
    
if __name__ == "__main__":

    from optparse import OptionParser

    parser = OptionParser()
    
    parser.add_option("-t","--tumor_dir",dest="tumor_dir",
                     help="the tumor junction directory")
    parser.add_option("-c", "--cancer", dest="cancer",
                      help="the TCGA cancer")

    (options, args) = parser.parse_args()

    main(options, args)
