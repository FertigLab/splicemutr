#!/usr/bin/env python

import pandas as pd
import numpy as np
import os
import splicemutr
import sys

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

def assign_kmers(genotypes_file,rows,hla_dir,hla_file):
    # genotypes_file: the genotypes_file, including path
    # rows: the specific_splice_dat rows
    # outputs:
        # specific_splice_dat with imm. kmers per sample columns

    genotypes_file = genotypes_file.values.tolist()
    all_kmers = {}
    tumor_kmers = [set() for i in range(len(rows))]
    normal_kmers = [set() for i in range(len(rows))]
    kmers = [[] for i in range(len(rows))]
    hlas = genotypes_file[0][0:6]
    sample_type = genotypes_file[0][7]
    sample_name = genotypes_file[0][8]

    for hla in hlas:
        if (type(hla) is float):
            break
        else:
            with open(hla_file%(hla_dir,hla)) as geno_file:
                print(hla_file%(hla_dir,hla))
                geno_list = geno_file.read().splitlines()
                geno_rows = [i.split('\t')[0] for i in geno_list]
                geno_kmers = [i.split('\t')[1] for i in geno_list]
                geno_dat = {int(geno_rows[i]):geno_kmers[i] for i in range(len(geno_rows))}
            for j in range(len(rows)):
                if rows[j] in geno_dat:
                    kmers[j].append(geno_dat[rows[j]])
        
            kmers = [":".join(k) for k in kmers]
            return(kmers,sample_name)

def filter_kmers(kmers,groups_dframe):
    groups_unique=list(set(groups_dframe.groups.tolist()))
    for group in groups_unique:
        groups_dframe_small=groups_dframe[groups_dframe["groups"]==group]
        normal_rows=groups_dframe_small[groups_dframe_small["deltapsi"]>0].rows.tolist()
        tumor_rows=groups_dframe_small[groups_dframe_small["deltapsi"]<0].rows.tolist()
        normal_kmers = set()
        [normal_kmers.update(kmers[i].split(":")) for i in normal_rows]
        tumor_kmers = set()
        [tumor_kmers.update(kmers[i].split(":")) for i in tumor_rows]
        for row in tumor_rows:
            kmers[row] = ":".join(list(set(kmers[row].split(":")).difference(normal_kmers)))
        for row in normal_rows:
            kmers[row] = ":".join(list(set(kmers[row].split(":")).difference(tumor_kmers)))
    return(kmers)

def main(options, args):
    
    # handling command line input
    genotypes_file = pd.read_table(options.genotypes_file)
    splice_dat = pd.read_table(options.splice_dat,sep=" ")
    hla_dir = options.summary_dir
    output_dir = options.output_dir
    summary_type = options.summary_type
    sample_num = int(options.sample_num)
    
    # processing command line arguments
    rows = range(len(splice_dat.index))
    if summary_type == "perc": # whether to load in the percentile summary or the IC50 summary
        hla_file = "%s/%s_tx_dict_summary_perc.txt"
    else:
        hla_file = "%s/%s_tx_dict_summary.txt"

    geno_length = len(genotypes_file.index)
    if sample_num > geno_length:
        print("sample_num > number of samples")
    else:
        genotypes_file_small = genotypes_file[sample_num:sample_num+1]
        try:
            splice_kmers,sample_name = assign_kmers(genotypes_file_small,rows,hla_dir,hla_file)
        except:
            sys.exit("Error: no genotypes for current sample")
        
    if "cluster" in splice_dat.columns.values: # determining whether the groups are leafcutter cluster or not
        groups=splice_dat.cluster.tolist()
    else:
        groups=splice_dat.groups.tolist()
    deltapsi=[float(i) for i in splice_dat.deltapsi.tolist()]
    groups_deltapsi_rows=pd.DataFrame({'groups':groups,'deltapsi':deltapsi,'rows':rows})
    splice_kmers_filt = filter_kmers(splice_kmers,groups_deltapsi_rows) # filtering by tumor and normal kmers
    splice_dat_frame = pd.DataFrame({"rows":rows,"kmers":splice_kmers_filt})
    splice_dat_frame.to_csv("%s/%s_splicemutr_kmers.txt"%(output_dir,sample_name),sep='\t',index=False)
        
if __name__ == "__main__":

    from optparse import OptionParser

    parser = OptionParser()
    
    parser.add_option("-g","--genotypes_file",dest="genotypes_file",
                     help="the genotypes file containing hla alleles per cancer sample")
    parser.add_option("-s","--summary_dir",dest="summary_dir",
                 help="the hla summary directory")
    parser.add_option("-d", "--splice_dat", dest="splice_dat",
                      help="the splicemutr data peptides only")
    parser.add_option("-o", "--output_dir", dest="output_dir",
                  help="the output directory")
    parser.add_option("-t", "--summary_type", dest="summary_type",
                  help="either 'perc' (percentile) or something else")
    parser.add_option("-n", "--sample_num", dest="sample_num",
                  help="the genotype row number")
    (options, args) = parser.parse_args()

    main(options, args)