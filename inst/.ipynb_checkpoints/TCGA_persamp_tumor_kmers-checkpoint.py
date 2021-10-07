#!/usr/bin/env python

# author: Theron Palmer
# created: 10/03/2021

# This script creates a per sample immunogenic kmers file per TCGA cancer analyzed

import pandas as pd
import os

def filter_tumor_kmers(kmer_row,norm_kmers_cluster):
    kmer_row=kmer_row.split("\t")
    split_kmers=[kmer_row[i].split(":") for i in range(1,len(kmer_row))]
    norm_set=set(norm_kmers_cluster)
    tum_spec_kmers=[":".join(list(set(kmers).difference(norm_set))) for kmers in split_kmers]
    return(tum_spec_kmers)

def get_norm_kmers(norm_rows,norm_kmers):
    norm_rows = [int(row) for row in norm_rows]
    norm_kmers_all = list()
    for row in norm_rows:
        norm_kmers_all += norm_kmers[row]
    norm_kmers_all = list(set(norm_kmers_all))
    return(norm_kmers_all)

def main(options,args):
    tumor_dir=options.tumor_dir
    cancer=options.cancer
    
    # reading in the necessary files
    
    genotypes_file = "%s/%s/%s_genotypes.txt"%(tumor_dir,cancer,cancer)
    with open(genotypes_file) as geno:
        genotypes = geno.read().splitlines()
    
    kmer_file = "%s/%s/kmer_files.txt"%(tumor_dir,cancer)
    with open(kmer_file) as kmer:
        kmer_files = kmer.read().splitlines()
        
    tum_norm_kmers_file = "%s/%s/%s_tumor_normal_kmers.txt"%(tumor_dir,cancer,cancer)
    with open(tum_norm_kmers_file) as tum_norm:
        tum_norm_kmers  = tum_norm.read().splitlines()
    norm_kmers = [tum_norm_kmers[i].split("\t")[2].split(":") for i in range(1,len(tum_norm_kmers))]
    
    splice_dat_file = "%s/%s/%s_splicemutr_dat.txt"%(tumor_dir,cancer,cancer)
    splice_dat = pd.read_table(splice_dat_file,sep="\t")
    splice_dat_clusters = splice_dat.cluster.tolist()
    splice_dat_rows = list(range(len(splice_dat_clusters)))
    splice_dat_deltapsi = [float(elem) for elem in splice_dat.deltapsi.tolist()]
    splice_dat_dict = {}
    for i in range(len(splice_dat_clusters)):
        cluster = splice_dat_clusters[i]
        row = splice_dat_rows[i]
        deltapsi = splice_dat_deltapsi[i]
        if cluster in splice_dat_dict:
            splice_dat_dict[cluster]["row"].append(row)
            splice_dat_dict[cluster]["deltapsi"].append(deltapsi)
        elif cluster not in splice_dat_dict:
            splice_dat_dict[cluster] = {}
            splice_dat_dict[cluster]["row"] = [row]
            splice_dat_dict[cluster]["deltapsi"] = [deltapsi]
        
    for i in range(len(kmer_files)):
        with open(kmer_files[i]) as specific_kmer:
            kmer_dat =  specific_kmer.read().splitlines()
            kmer_dat = kmer_dat[1:]
            tumor_kmers_all = []
            iter_val=1
        for cluster in splice_dat_dict.keys():
            if iter_val % 1000 == 0:
                print("%d:%d"%(iter_val,len(splice_dat_clusters)))
            rows = splice_dat_dict[cluster]["row"]
            deltapsi = splice_dat_dict[cluster]["deltapsi"]
            norm_rows = [rows[i] for i in range(len(deltapsi)) if deltapsi[i] < 0]
            tum_rows = [rows[i] for i in range(len(deltapsi)) if deltapsi[i] > 0]
            norm_kmers_cluster = get_norm_kmers(norm_rows,norm_kmers)
            for row in tum_rows:
                tumor_kmers = kmer_dat[row]
                tumor_kmers_all.append("\t".join(filter_tumor_kmers(tumor_kmers,norm_kmers_cluster)+[str(row),cluster]))
            iter_val+=1
            
        
        tumor_kmers_file_str = "\n".join(tumor_kmers_all)
        tumor_kmer_file = "%s/%s_kmers_%d_filt.txt"%(os.path.dirname(kmer_file),cancer,i+1)
        with open(tumor_kmer_file,"w") as k:
            print(tumor_kmers_file_str,file=k)
        
if __name__ == "__main__":

    from optparse import OptionParser

    parser = OptionParser()
    
    parser.add_option("-t","--tumor_dir",dest="tumor_dir",
                     help="the tumor junction directory")
    parser.add_option("-c", "--cancer", dest="cancer",
                      help="the TCGA cancer")

    (options, args) = parser.parse_args()

    main(options, args)