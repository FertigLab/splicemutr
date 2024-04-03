#!/usr/bin/env python

import pandas as pd

def main(options):
    
    neosplice_output_file = options.neosplice_output_file
    ieatlas_file = options.ieatlas_file
    out_dir = options.out_dir
    sample_id = options.sample_id

    #------------------------------------------------------------------#
    ## Loading in the neosplice output file

    neosplice_output = pd.read_csv(neosplice_output_file,sep="\t")
    neosplice_peptides = set(neosplice_output.Variant_peptide_sequence.tolist())
    
    #------------------------------------------------------------------#
    ## Loading in the non-coding orf data
    
    epitopes_in_cancer = pd.read_csv(ieatlas_file,sep="\t")
    epitopes = epitopes_in_cancer.Sequence.tolist()

    #------------------------------------------------------------------#
    ## g) preparing the per sample target data
        
    overlap = len(neosplice_peptides.intersection(epitopes))
    total = len(neosplice_peptides)
    
    out_file = "%s/%s_pep_overlap.txt"%(out_dir,sample_id)
    print(out_file)
    
    data = [[sample_id,overlap,total]]
    pd.DataFrame(data, columns=['sample_id','overlap','total']).to_csv(out_file,sep="\t",index=False,header=False)
    
if __name__ == "__main__":

    from optparse import OptionParser

    parser = OptionParser()

    parser.add_option("-i", "--ieatlas_file", dest="ieatlas_file",
                      help="the ieatlas file to be used for mass spec validation")
    parser.add_option("-n", "--neosplice_output_file", dest="neosplice_output_file",
            help="the neosplice output file")
    parser.add_option("-o", "--out_dir", dest="out_dir",
                      help="the output_directory")
    parser.add_option("-s", "--sample_id", dest="sample_id",
                help="the sample_id")
    (options, args) = parser.parse_args()

    main(options)
