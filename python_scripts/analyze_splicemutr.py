#!/usr/bin/env python

import pandas as pd
import numpy as np
import os
import splicemutr as sp
import sys

def main(options, args):
    
    # handling command line input
    genotypes_file = pd.read_csv(options.genotypes_file,header=None,sep="\t")
    splice_dat = pd.read_csv(options.splice_dat,sep="\t") # change the save .txt files to have \t separator
    hla_dir = options.summary_dir
    output_dir = options.output_dir
    summary_type = options.summary_type
    sample_num = int(options.sample_num)-1
    fasta_file = int(options.fasta_file)

    ref_kmers = generate_ref_kmer_set(fasta_file)
    
    # processing command line arguments
    rows = range(len(splice_dat.index))
    if summary_type == "perc": # whether to load in the percentile summary or the IC50 summary
        hla_file = "%s/%s_tx_dict_summary_perc.txt"
    else:
        hla_file = "%s/%s_tx_dict_summary.txt"

    geno_length = len(genotypes_file.index)
    if sample_num > geno_length:
        print("sample_num > number of samples")
    else:
        genotypes_file_small = genotypes_file[sample_num:sample_num+1]
        try:
            splice_kmers,sample_name = sp.assign_kmers(genotypes_file_small,rows,hla_dir,hla_file)
            if "cluster" in splice_dat.columns.values: # determining whether the groups are leafcutter cluster or not
                groups=splice_dat.cluster.tolist()
            else:
                groups=splice_dat.groups.tolist()
            deltapsi=[float(i) for i in splice_dat.deltapsi.tolist()]
            groups_deltapsi_rows=pd.DataFrame({'groups':groups,'deltapsi':deltapsi,'rows':rows})
            splice_kmers_filt = sp.filter_kmers(splice_kmers,groups_deltapsi_rows,ref_kmers) # filtering by tumor and normal kmers
            splice_dat_frame = pd.DataFrame({"rows":rows,"kmers":splice_kmers_filt})
            splice_dat_frame.to_csv("%s/%s_splicemutr_kmers.txt"%(output_dir,sample_name),sep='\t',index=False)
        except:
            print("Error: no genotypes for current sample")

if __name__ == "__main__":

    from optparse import OptionParser

    parser = OptionParser()
    
    parser.add_option("-g","--genotypes_file",dest="genotypes_file",
                     help="the genotypes file containing hla alleles per cancer sample")
    parser.add_option("-s","--summary_dir",dest="summary_dir",
                 help="the hla summary directory")
    parser.add_option("-d", "--splice_dat", dest="splice_dat",
                      help="the splicemutr data peptides only")
    parser.add_option("-o", "--output_dir", dest="output_dir",
                  help="the output directory")
    parser.add_option("-t", "--summary_type", dest="summary_type",
                  help="either 'perc' (percentile) or something else")
    parser.add_option("-n", "--sample_num", dest="sample_num",
                  help="the genotype row number")
    parser.add_option("-f", "--fasta_file", dest="fasta_file",
                help="the fasta file input")
    (options, args) = parser.parse_args()

    main(options, args)