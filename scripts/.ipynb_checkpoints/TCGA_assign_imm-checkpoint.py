#!/usr/bin/env python

def create_tcga_splicemutr(intron,splicemutr_dat):
    return()

def main(options, args):
    

if __name__ == "__main__":

    from optparse import OptionParser

    parser = OptionParser()
    
    parser.add_option("-d","--dat_file",dest="dat_file",
                     help="the splicemutr output file file")
    parser.add_option("-s", "--splice_dat", dest="splice_dat",
                      help="full splicemutr data")
    parser.add_option("-a", "--genotypes_files", dest="genotypes_files",
                      help="the genotypes file containing hla alleles")
    parser.add_option("-g", "--groups_file", dest="groups_file",
                      help="the groups file")
    parser.add_option("-m", "--junc_metadata", dest="junc_metadata",
                      help="the rse metadata file")
    (options, args) = parser.parse_args()

    main(options, args)