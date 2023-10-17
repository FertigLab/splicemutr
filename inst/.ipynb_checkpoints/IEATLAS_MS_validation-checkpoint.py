#!/usr/bin/env python

import pandas as pd
from collections import defaultdict
import pickle
from os.path import basename

def main(options):
    
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
    ## Loading in the non-coding orf data
    
    epitopes_in_cancer = pd.read_csv(ieatlas_file,sep="\t")
    epitopes = epitopes_in_cancer.Sequence.tolist()

    #------------------------------------------------------------------#
    ## g) preparing the per sample target data
        
    gene_kmers_set = set()
    gene_kmers_set_fill = [gene_kmers_set.update(i) for i in list(gene_kmers.values())]
    overlap = len(gene_kmers_set.intersection(set(epitopes)))
    total = len(gene_kmers_set)
    
    out_file = "%s/%s_psm_overlap.txt"%(out_dir,sample_id)
    
    data = [[sample_id,overlap,total]]
    pd.DataFrame(data, columns=['sample_id','overlap','total']).to_csv(out_file,sep="\t",index=False,header=False)
    
if __name__ == "__main__":

    from optparse import OptionParser

    parser = OptionParser()

    parser.add_option("-t", "--target_genes_file", dest="target_genes_file",
                      help="the BRCA target genes file output from slicemutr collate TCGA data")
    parser.add_option("-s", "--splicemutr_data_file", dest="splicemutr_data_file",
                      help="the splicemutr (BRCA) all_pep.txt file")
    parser.add_option("-g", "--genotypes_file", dest="genotypes_file",
                      help="the TCGA genotyps file containing external id and TCGA barcode mappings")
    parser.add_option("-o", "--out_dir", dest="out_dir",
                      help="the output_directory")
    parser.add_option("-i", "--ieatlas_file", dest="ieatlas_file",
                      help="the ieatlas epitopes file, either basic or cancer, to be used for MS validation")
    
    (options, args) = parser.parse_args()

    main(options)