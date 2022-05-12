#!/usr/bin/env python

# The goal of this script is to take in the CPTAC peptides file and process it  in 2 ways:
    # 1) So that there is a per sample summary of peptides found in the sample
    # 2) So that each sample is matched to its respective TCGA barcode

import pandas as pd
import os
import sys
from collections import defaultdict

def main(options, args):
    cptac_data = pd.read_table(options.cptac_file,sep="\t")
    genotypes_data = pd.read_table(options.genotypes,sep="\t") # change the save .txt files to have \t separator
    out_dir = options.out_dir

    # process the genotypes data
    barcode = genotypes_data["aliquot_id"].tolist()
    external_id = genotypes_data["external_id"].tolist()
    
    # process the cptac data
    cptac_dict = defaultdict(list)

    peptides=cptac_data["Peptide"].tolist()
    samples=cptac_data["Sample"].tolist()
    samples = [i.split(";") for i in samples]
    print("Creating sample-peptide dict")
    for i in range(len(samples)):
        for sample in samples[i]:
            cptac_dict[sample].append(peptides[i])
    peptides_list=[]
    cptac_keys=list(cptac_dict.keys())
    cptac_barcode=[]
    cptac_external_id=[]
    for key in cptac_keys:
        peptides_list.append(":".join(cptac_dict[key]))
        first=key.split(":")[0]
        barcode_row=[i for i in range(len(barcode)) if first in barcode[i]]
        if (len(barcode_row)==0):
            cptac_barcode.append("NA")
            cptac_external_id.append("NA")
        else:
            cptac_barcode.append(barcode[barcode_row[0]])
            cptac_external_id.append(external_id[barcode_row[0]])
    
    cptac_df = pd.DataFrame({"cptac_sample":cptac_keys,"TCGA_barcode":cptac_barcode,"TCGA_external_id":cptac_external_id,"peptides":peptides_list})
    cptac_df.to_csv(out_dir+"/cptac_sample_peptides.txt",sep="\t")

if __name__ == "__main__":

    from optparse import OptionParser

    parser = OptionParser()
    
    parser.add_option("-c","--cptac_file",dest="cptac_file",
                     help="the .peptides.tsv file from cptac")
    parser.add_option("-g","--genotypes",dest="genotypes",
                 help="the genotypes file that contains the associated TCGA barcodes")
    parser.add_option("-o","--out_dir",dest="out_dir",
                 help="the output directory to write the summary to")
    (options, args) = parser.parse_args()

    main(options, args)