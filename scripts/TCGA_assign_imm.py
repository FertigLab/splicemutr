#!/usr/bin/env python

# author: Theron Palmer
# created: 08/13/2021

# This script takes in the genotypes file for the specific cancer type, the full splicemutr data file, the leafcutter groups file for the specific cancer, the cancer intron file, and the HLA allele summary file directory. Then it assigns kmers to the appropriate intron rows based on the genotype. 

import pandas as pd

def create_juncs(intron_dat):
    chr_dat = intron_dat.chr.tolist()
    start_dat = intron_dat.start.tolist()
    end_dat = intron_dat.end.tolist()
    juncs = [":".join([str(chr_dat[i]),str(start_dat[i]),str(end_dat[i])]) for i in range(len(chr_dat))]
    return(juncs)


def create_specific_splicemutr(splice_dat, intron_file):
    # splice_dat: the full splicemutr output pandas dataframe
    # intron file: the leafcutter output for the current cancer type
    # outputs:
        # specific_splicemutr_file: the splice dat file that only contains rows found in the leafcutter output
    
    splice_dat_juncs = create_juncs(splice_dat)
    intron_juncs = create_juncs(intron_file)
    splice_dat["juncs"] = splice_dat_juncs
    intron_file["juncs"] = intron_juncs
    
    intron_file = intron_file[intron_file.juncs.isin(splice_dat_juncs)]
    intron_juncs = intron_file.juncs.tolist()
    splice_dat = splice_dat[splice_dat.juncs.isin(intron_juncs)]
    for i in range(len(intron_juncs)):
        splice_dat_sub = splice_dat
        
    return(intron_loc)
    
def main(options, args):
        
        genotypes_file = pd.read_table(options.genotypes_file)
        splice_dat = pd.read_table(options.splice_dat,sep=" ")
        groups_file = pd.read_table(options.groups_file)
        intron_file = pd.read_table(options.intron_file)
        hla_dir = options.hla_dir
        
        groups_file.columns = ["file_name","tumor_normal"]
        num_samples = len(groups_file.index)
        
        genotypes_file = genotypes_file.values.tolist()
        intron_file = intron_file[intron_file["verdict"] != "unknown_strand"]

        print(create_specific_splicemutr(splice_dat,intron_file))

if __name__ == "__main__":

    from optparse import OptionParser

    parser = OptionParser()
    
    parser.add_option("-g","--genotypes_file",dest="genotypes_file",
                     help="the genotypes file containing hla alleles per cancer sample")
    parser.add_option("-s", "--splice_dat", dest="splice_dat",
                      help="full splicemutr data")
    parser.add_option("-r", "--groups_file", dest="groups_file",
                      help="the groups file")
    parser.add_option("-i", "--intron_file", dest="intron_file",
                       help = "the intron file")
    parser.add_option("-a", "--hla_dir", dest = "hla_dir",
                     help = "the hla summary directory")
    (options, args) = parser.parse_args()

    main(options, args)