#!/usr/bin/env python

# author: Theron Palmer
# created: 08/13/2021

# This script takes in the genotypes file for the specific cancer type, the full splicemutr data file, the leafcutter groups file for the specific cancer, the cancer intron file, and the HLA allele summary file directory. Then it assigns kmers to the appropriate intron rows based on the genotype. 

import pandas as pd
import numpy as np
import os 

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

def assign_kmers(genotypes_file,rows,hla_dir):
    # genotypes_file: the genotypes_file, including path
    # rows: the specific_splice_dat rows
    # outputs:
        # specific_splice_dat with imm. kmers per sample columns
        
    hla_file = "%s/%s_tx_dict_summary_perc.txt"
    genotypes_file = genotypes_file.values.tolist()
    all_kmers = {}
    #print(all_kmers)
    for i in range(len(genotypes_file)):
        print(i)
        kmers = [[] for i in range(len(rows))]
        hlas = genotypes_file[i][0:6]
        sample = genotypes_file[i][7]
        sample_type = genotypes_file[i][8]
        for hla in hlas:
            print(hla)
            if (type(hla) is float):
                break
            else:
                geno_file = pd.read_table(hla_file%(hla_dir,hla))
                geno_file.columns = ["row","kmer","perc"]
                geno_file = geno_file[geno_file.row.isin(rows)]
                geno_rows = geno_file.row.tolist()
                geno_kmers = geno_file.kmer.tolist()
                rows_fill = pd.DataFrame(rows)
                rows_fill.columns = ["row"]
                rows_fill = rows_fill[rows_fill.row.isin(geno_file.row.tolist())]
                rows_fill = rows_fill.row.tolist()
                geno_splice_dat_match = [geno_rows.index(j) for j in rows_fill]
                for j in range(len(rows_fill)):
                    kmers[j] = kmers[j] + geno_kmers[geno_splice_dat_match[j]].split(':')
        for j in range(len(rows)):
            kmers[j] = ":".join(kmers[j])
        all_kmers[i] = kmers
    all_kmers = pd.DataFrame(all_kmers)
    return(all_kmers)
        
def main(options, args):
        
        genotypes_file = pd.read_table(options.genotypes_file)
        splice_dat = pd.read_table(options.splice_dat,sep=" ")
        groups_file = pd.read_table(options.groups_file)
        intron_file = pd.read_table(options.intron_file)
        hla_dir = options.hla_dir
        cancer_dir = os.path.basename(os.path.dirname(options.intron_file))
        groups_file.columns = ["file_name","tumor_normal"]
        num_samples = len(groups_file.index)
        
        intron_file = intron_file[intron_file["verdict"] != "unknown_strand"]
        
        # creating the specific splicemutr dat that corresponds to the cancer being analyzed
        specific_splice_dat = create_specific_splicemutr(splice_dat,intron_file)
        specific_splice_dat.to_csv("%s/%s_splicemutr_dat.txt"%(cancer_dir,os.path.basename(cancer_dir)),sep='\t')
        rows = specific_splice_dat.rows.tolist()
        
        # assigning the immunogenic kmers to the specific_splice_dat
        specific_splice_kmers = assign_kmers(genotypes_file,rows,hla_dir)
        specific_splice_kmers.to_csv("%s/%s_kmers.txt"%(cancer_dir,os.path.basename(cancer_dir)),sep='\t')
        
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
