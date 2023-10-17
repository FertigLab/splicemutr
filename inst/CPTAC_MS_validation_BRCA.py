#!/usr/bin/env python

import pandas as pd
import re
from collections import defaultdict
import pickle
from os.path import basename

def main(options):
    
    psm_file = options.psm_file # the psm file for the BRCA CPTAC experiments iTRAQ4
    metadata_file = options.metadata_file # the CPTAC metdata file for the 
    target_genes_file = options.target_genes_file
    splicemutr_data_file = options.splicemutr_data_file
    kmers_file = options.kmers_file
    genotypes_file = options.genotypes_file
    out_dir = options.out_dir
    
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
    ## a) Reading in the target genes

    target_genes = pd.read_csv(target_genes_file)
    target_genes_list = target_genes.relevant_genes.tolist()
    
    #------------------------------------------------------------------#
    ## b) Reading in the combined splicemutr file and filtering by target gene
    
    splicemutr_data = pd.read_csv(splicemutr_data_file,sep="\t",index_col=False)
    splicemutr_data_rows = range(len(splicemutr_data.gene.tolist()))
    splicemutr_data["rows_to_keep"] = list(splicemutr_data_rows)
    splicemutr_data_targets = splicemutr_data[splicemutr_data.gene.isin(target_genes_list) &
                                              splicemutr_data.deltapsi.tolist() > 0]
    rows_to_keep = splicemutr_data_targets.rows_to_keep.tolist()
    genes = splicemutr_data_targets.gene.tolist()
    juncs = splicemutr_data_targets.juncs.tolist()

    #------------------------------------------------------------------#
    ## c) Reading in a kmers file
    
    kmers = pd.read_csv(kmers_file,sep="\t")
    kmers_to_keep = kmers[kmers.rows.isin(rows_to_keep)]
    kmers_to_keep.assign(gene=genes)
    kmers_to_keep.assign(juncs=juncs)

    #------------------------------------------------------------------#
    ## d) Generating per gene kmer list
    
    genes = list(set(kmers_to_keep.gene.tolist()))
    gene_kmers = defaultdict(set)
    num_genes = len(genes)

    for gene_num in range(num_genes):
        if gene_num % 1000 == 0: print("%d/%d"%(gene_num,num_genes))
        gene = genes[gene_num]
        kmers_to_keep_small = kmers_to_keep[kmers_to_keep["gene"]==gene]
        for i in kmers_to_keep_small.kmers.tolist():
            try:
                gene_kmers_fill = [gene_kmers[gene].update(j) for j in [i.split(":")]]
            except:
                gene_kmers_fill = []
    
    #------------------------------------------------------------------#
    ## f) Reading in the genotypes file and determining the sample_ID for the sample being analyzed

    genotypes = pd.read_csv(genotypes_file,sep="\t")
    sample_ids = ["-".join(i.split("-")[0:3]) for i in genotypes.aliquot_id.tolist()]
    genotypes["sample_ids"] = sample_ids
    samples = dict(zip(genotypes.external_id.tolist(), genotypes.sample_ids.tolist()))
    external_ids= dict(zip(genotypes.sample_ids.tolist(),genotypes.external_id.tolist()))

    external_id = basename(kmers_file.replace("_splicemutr_kmers.txt",""))
    sample_id = samples[external_id]
    
    #------------------------------------------------------------------#
    ## g) preparing the per sample target data
    
    psm_per_sample_target = psm_per_sample_list[sample_id]
    
    gene_kmers_set = set()
    gene_kmers_set_fill = [gene_kmers_set.update(i) for i in list(gene_kmers.values())]
    overlap = len(gene_kmers_set.intersection(set(psm_per_sample_target)))
    total = len(gene_kmers_set)
    
    out_file = "%s/%s_psm_overlap.txt"%(out_dir,sample_id)
    
    data = [[sample_id,overlap,total]]
    pd.DataFrame(data, columns=['sample_id','overlap','total']).to_csv(out_file,sep="\t",index=False,header=False))
    
if __name__ == "__main__":

    from optparse import OptionParser

    parser = OptionParser()

    psm_file = options.psm_file
    metadata_file = options.metadata_file
    target_genes_file = options.target_genes_file
    splicemutr_data_file = options.splicemutr_data_file
    kmers_file = options.kmers_file
    genotypes_file = options.genotypes_file
    
    parser.add_option("-p", "--psm_file", dest="psm_file",
                      help="the BRCA psm file")
    parser.add_option("-m", "--metadata_file", dest="metadata_file",
                      help="the CPTAC metadata file")
    parser.add_option("-t", "--target_genes_file", dest="target_genes_file",
                      help="the BRCA target genes file output from slicemutr collate TCGA data")
    parser.add_option("-s", "--splicemutr_data_file", dest="splicemutr_data_file",
                      help="the splicemutr (BRCA) all_pep.txt file")
    parser.add_option("-g", "--genotypes_file", dest="genotypes_file",
                      help="the TCGA genotyps file containing external id and TCGA barcode mappings")
    parser.add_option("-o", "--out_dir", dest="out_dir",
                      help="the output_directory")
    (options, args) = parser.parse_args()

    main(options)