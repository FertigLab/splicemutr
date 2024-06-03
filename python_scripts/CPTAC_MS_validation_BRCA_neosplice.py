#!/usr/bin/env python

import pandas as pd
import re
from collections import defaultdict
import pickle
from os.path import basename

def main(options):
    
    psm_file = options.psm_file # the psm file for the BRCA CPTAC experiments iTRAQ4
    metadata_file = options.metadata_file # the CPTAC metdata file for the 
    out_dir = options.out_dir
    sample_id = options.sample_id
    neosplice_output_file = options.neosplice_output_file
    
    #------------------------------------------------------------------#
    # Reading in the psm file for BRCA and generating dictionary of kmers per iTRAQ sample
    #

    psm = pd.read_csv(psm_file,sep=" ",header=None)
    psm.columns = ["file","psm"]

    #------------------------------------------------------------------#
    ## Cleaning modification information from psm

    psms = [re.sub("[0-9,.,+]+",'',i) for i in psm.psm.tolist()]
    psm['psm_clean']=psms
    
    #------------------------------------------------------------------#
    ## Formatting psms as file header and list of cleaned psms dictionary
    
    psm_dict = defaultdict(list)
    psm_files = ["_".join(i.split("_")[0:9]) for i in psm.file.tolist()]
    psms = psm.psm_clean.tolist()

    for i in range(len(psm_files)):
        psm_file = psm_files[i]
        psm_dict[psm_file].append(psms[i])
    
    #------------------------------------------------------------------#
    ## Loading in metadata for psm files
    
    metadata = pd.read_excel(metadata_file,nrows=38)
    filenames = metadata.FileName.tolist()
    biospecimen_114 = metadata["114-Biospecimen"].tolist()
    biospecimen_115 = metadata["115-Biospecimen"].tolist()
    biospecimen_116 = metadata["116-Biospecimen"].tolist()
    filenames_full = filenames+filenames+filenames
    filenames_full = [i.replace("_f01.raw","") if not isinstance(i,float) else i for i in filenames_full]
    biospecimens = biospecimen_114+biospecimen_115+biospecimen_116
    biospecimens_sample_ids = ["-".join(i.split("-")[0:3]) if not isinstance(i,float) else i for i in biospecimens]

    sample_ids_filenames = defaultdict(list)

    for i in range(len(biospecimens_sample_ids)):
        sample_ids_filenames[filenames_full[i]].append(biospecimens_sample_ids[i])
    
    #------------------------------------------------------------------#
    # Obtaining psm per sample 
    
    psm_per_sample = defaultdict(set)
    psm_per_sample_list = defaultdict(list)
    psm_keys = sample_ids_filenames.keys()

    for key in psm_keys:
        for sample in sample_ids_filenames[key]:
            psm_per_sample[sample].update(psm_dict[key])

    for sample in psm_per_sample.keys():
        psm_per_sample_list[sample]=list(psm_per_sample[sample])
    
    #------------------------------------------------------------------#
    # Reading in Neosplice output file and obtaining kmer list per sample
        
    neosplice_output = pd.read_csv(neosplice_output_file,sep="\t")
    neosplice_peptides = set(neosplice_output.Variant_peptide_sequence.tolist())
    sample_peptides = set(psm_per_sample_list[sample_id])

        
    #------------------------------------------------------------------#
    ## g) preparing the per sample target data
    
    overlap = len(neosplice_peptides.intersection(sample_peptides))
    total = len(neosplice_peptides)

    out_file = "%s/%s_psm_overlap.txt"%(out_dir,sample_id)

    data = [[sample_id,overlap,total]]
    pd.DataFrame(data, columns=['sample_id','overlap','total']).to_csv(out_file,sep="\t",index=False,header=False)
    
if __name__ == "__main__":

    from optparse import OptionParser

    parser = OptionParser()
    
    parser.add_option("-p", "--psm_file", dest="psm_file",
                      help="the BRCA psm file")
    parser.add_option("-m", "--metadata_file", dest="metadata_file",
                      help="the CPTAC metadata file")
    parser.add_option("-o", "--out_dir", dest="out_dir",
                      help="the output_directory")
    parser.add_option("-s", "--sample_id", dest="sample_id",
                    help="the sample_id")
    parser.add_option("-n", "--neosplice_output_file", dest="neosplice_output_file",
                help="the neosplice output file")
    (options, args) = parser.parse_args()

    main(options)