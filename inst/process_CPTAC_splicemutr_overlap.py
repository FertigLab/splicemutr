#!/usr/bin/env python

# The goal of this script is to take in the CPTAC peptides file and process it  in 2 ways:
    # 1) So that there is a per sample summary of peptides found in the sample
    # 2) So that each sample is matched to its respective TCGA barcode

import pandas as pd
import os
import sys

def main(options, args):
    cptac_data = pd.read_table(options.cptac_file,sep="\t")
    cptac_row = int(options.cptac_row)
    splicemutr_files = pd.read.table(options.splicemutr_files,header=None,names=["file"])
    out_dir = options.out_dir

    external_id = cptac_data["TCGA_external_id"].tolist()[cptac_row]
    cptac_peptides = cptac_data["peptides"].split(":")
    K=9
    cptac_peptides_below_Kmer = []
    cptac_peptides_Kmer = []
    for pep in cptac_peptides:
        if len(pep)<K:
            cptac_peptides_below_Kmer.append(pep)
        else:
            cptac_peptides_Kmer+=[pep[i:(i+K)] for i in range(0,len(a)-K,1)]
    
    splicemutr_file=[i for i in splicemutr_files["file"].tolist() if external_id in i]
    splicemutr_file=splicemutr_file[0]
    splicemutr_data=pd.read_table(splicemutr_file,sep="\t")
    kmers=splicemutr_data["kmers"].tolist()
    kmers=":".join(kmers)
    kmers=kmers.split(":")

    Kmers_in_exact = [i for i in kmers if i in cptac_peptides_Kmer]
    Kmers_not_in_exact = [i for i in kmers if i not in Kmers_in_exact]
    Kmers_in_below = [i for i in Kmers_not_in_exact if i in cptac_peptides_below_Kmer]
    splicemutr_Kmers_in_cptac = Kmers_in_exact
    splicemutr_Kmers_in_cptac += Kmers_not_in_exact

    splicemutr_Kmers_in_cptac_df = pd.DataFrame({"splicemutr_kmers_count":len(kmers),
    "splicemutr_kmers_in_cptac_count":len(splicemutr_Kmers_in_cptac),
    "splicemutr_kmers":kmers,
    "splicemutr_Kmers_in_cptac":splicemutr_Kmers_in_cptac})

    splicemutr_Kmers_in_cptac_df.to_csv(out_dir+"/"+external_id+"_splicemutr_kmers_in_cptac.txt",sep="\t", index=False)
    




if __name__ == "__main__":

    from optparse import OptionParser

    parser = OptionParser()
    
    parser.add_option("-c","--cptac_file",dest="cptac_file",
                     help="the .peptides.tsv file from cptac")
    parser.add_option("-r","--cptac_row",dest="cptac_row",
                 help="the cptac file row to process")
    parser.add_option("-s","--splicemutr_files",dest="splicemutr_files",
                 help="the splicemutr files")
    parser.add_option("-o","--out_dir",dest="out_dir",
                 help="the output directory to write the summary to")
    (options, args) = parser.parse_args()

    main(options, args)