#!/usr/bin/env python

# author: Theron Palmer
# created: 08/13/2021

# This script takes in the genotypes file for the specific cancer type, the full splicemutr data file, the leafcutter groups file for the specific cancer, the cancer intron file, and the HLA allele summary file directory. Then it assigns kmers to the appropriate intron rows based on the genotype. 

import pandas as pd
import numpy as np
import os
import splicemutr

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
    num_splice_dat = len(splice_dat_juncs)
    splice_dat["rows"] = range(num_splice_dat)
    intron_juncs = create_juncs(intron_file)
    splice_dat["juncs"] = splice_dat_juncs
    intron_file["juncs"] = intron_juncs
    
    intron_file = intron_file[intron_file.juncs.isin(splice_dat_juncs)]
    intron_juncs = intron_file.juncs.tolist()
    splice_dat = splice_dat[splice_dat.juncs.isin(intron_juncs)]
    splice_dat_juncs = splice_dat.juncs.tolist()
    intron_loc = [intron_juncs.index(j) for j in splice_dat_juncs]
    
    all_clusters = intron_file.clusterID.tolist()
    all_verdicts = intron_file.verdict.tolist()
    all_deltapsis = intron_file.deltapsi.tolist()

    splice_dat_new = pd.DataFrame({'cluster':[all_clusters[i] for i in intron_loc],
                      'chr': splice_dat.chr.tolist(),
                      'start': splice_dat.start.tolist(),
                      'end': splice_dat.end.tolist(),
                      'gene': splice_dat.gene.tolist(),
                      'tx_id': splice_dat.tx_id.tolist(),
                      'modified': splice_dat.modified.tolist(),
                      'is_UTR': splice_dat.is_UTR.tolist(),
                      'tx_junc_loc': splice_dat.tx_junc_loc.tolist(),
                      'pep_junc_loc': splice_dat.pep_junc_loc.tolist(),
                      'verdict':[all_verdicts[i] for i in intron_loc],
                      'deltapsi':[all_deltapsis[i] for i in intron_loc],
                      'error': splice_dat.error.tolist(),
                      'peptide': splice_dat.peptide.tolist(),
                      'tx_length': splice_dat.tx_length.tolist(), 
                      'start_stop': splice_dat.start_stop.tolist(),
                      'protein_coding': splice_dat.protein_coding.tolist(),
                      'rows': splice_dat.rows.tolist(),
                      'juncs': splice_dat.juncs.tolist()})
    
    return(splice_dat_new)

def assign_kmers(genotypes_file,rows,hla_dir,cancer):
    # genotypes_file: the genotypes_file, including path
    # rows: the specific_splice_dat rows
    # outputs:
        # specific_splice_dat with imm. kmers per sample columns
        
    hla_file = "%s/%s_tx_dict_summary_perc.txt"
    genotypes_file = genotypes_file.values.tolist()
    all_kmers = {}
    hla_dict = {}
    tumor_kmers = [set() for i in range(len(rows))]
    normal_kmers = [set() for i in range(len(rows))]
    for i in range(len(genotypes_file)):
        print("%s:%d:%d"%(cancer,i,len(genotypes_file)),flush=True)
        kmers = [[] for i in range(len(rows))]
        hlas = genotypes_file[i][0:6]
        sample_type = genotypes_file[i][7]
        for hla in hlas:
            if (type(hla) is float):
                break
            else:
                if hla in hla_dict:
                    geno_dat = hla_dict[hla]
                else:
                    with open(hla_file%(hla_dir,hla)) as geno_file:
                        geno_list = geno_file.read().splitlines()
                        geno_rows = [i.split('\t')[0] for i in geno_list]
                        geno_kmers = [i.split('\t')[1] for i in geno_list]
                        hla_dict[hla] = {int(geno_rows[i]):geno_kmers[i] for i in range(len(geno_rows))}
                        geno_dat = hla_dict[hla]
                for j in range(len(rows)):
                    if rows[j] in geno_dat:
                        kmers[j].append(geno_dat[rows[j]])
        
        kmers = [":".join(k) for k in kmers]
        all_kmers[i] = kmers
        if sample_type == "T":
            for i in range(len(rows)):
                tumor_kmers[i].update(kmers[i].split(":"))
        else:
            for i in range(len(rows)):
                normal_kmers[i].update(kmers[i].split(":"))
    all_kmers = pd.DataFrame(all_kmers)
    return(all_kmers,tumor_kmers,normal_kmers)
        
def main(options, args):
        
        genotypes_file = pd.read_table(options.genotypes_file)
        splice_dat = pd.read_table(options.splice_dat,sep=" ")
        groups_file = pd.read_table(options.groups_file)
        intron_file = pd.read_table(options.intron_file)
        hla_dir = options.hla_dir
        cancer_dir = os.path.dirname(options.intron_file)
        groups_file.columns = ["file_name","tumor_normal"]
        num_samples = len(groups_file.index)
        
        intron_file = intron_file[intron_file["verdict"] != "unknown_strand"]
        
        # creating the specific splicemutr dat that corresponds to the cancer being analyzed
        specific_splice_dat = create_specific_splicemutr(splice_dat,intron_file)
        specific_splice_dat.to_csv("%s/%s_splicemutr_dat.txt"%(cancer_dir,os.path.basename(cancer_dir)),sep='\t')
        rows = specific_splice_dat.rows.tolist()
        
        # assigning the immunogenic kmers to the specific_splice_dat
        geno_length = len(genotypes_file.index)
        iter_val = 1
        #sample_num=50
        for i in range(0,geno_length,50): # original sample_num == 100
            if i+50 <= geno_length:
                end = i+50
            else:
                end = geno_length
            genotypes_file_small = genotypes_file[i:end]
            specific_splice_kmers,tumor_kmers_fill,normal_kmers_fill = assign_kmers(genotypes_file_small,
                                                                          rows,
                                                                          hla_dir,
                                                                          os.path.basename(cancer_dir))
#             if i == 0:
#                 tumor_kmers = tumor_kmers_fill
#                 normal_kmers = normal_kmers_fill
#             else:
#                 for i in range(len(rows)):
#                     tumor_kmers[i] = tumor_kmers[i].union(tumor_kmers_fill[i])
#                     normal_kmers[i] = normal_kmers[i].union(normal_kmers_fill[i])
            
#             # turning tumor and normal kmers into pandas dataframe for saving
            
            specific_splice_kmers.to_csv("%s/%s_kmers_%d.txt"%(cancer_dir,os.path.basename(cancer_dir),iter_val),sep='\t')
            iter_val+=1
            
#         tumor_normal_kmers = {"tumor":[":".join(list(i)) for i in tumor_kmers],
#                               "normal":[":".join(list(i)) for i in normal_kmers]}
#         tumor_normal_kmers = pd.DataFrame(tumor_normal_kmers)
#         tumor_normal_kmers.to_csv("%s/%s_tumor_normal_kmers.txt"%(cancer_dir,os.path.basename(cancer_dir)),sep='\t')
        
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
