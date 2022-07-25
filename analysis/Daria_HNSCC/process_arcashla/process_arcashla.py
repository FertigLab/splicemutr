#!/usr/bin/env python

import sys

try:
    import pickle
except:
    print("You do not have the correct python modules installed", file=sys.stderr)
    sys.exit()

def get_abundance(log_file_contents):
        line_start = [i for i in range(len(log_file_contents)) if "Observed HLA genes" in log_file_contents[i]][0]
        line_end = [i for i in range(line_start,len(log_file_contents)) if "---" in log_file_contents[i]][0]
        gene_row = log_file_contents[line_start+1].find("gene")
        abundance_row = log_file_contents[line_start+1].find("abundance")
        
        # extracting the genes and abundance info to obtain the genotype
        gene_info = [log_file_contents[i][gene_row:(gene_row+8)].replace(' ','') 
                     for i in range(line_start+2,line_end-1)]
        abundance_info = [float(log_file_contents[i][abundance_row:(abundance_row+len("abundance")-1)].replace(' ','')) 
                          for i in range(line_start+2,line_end-1)]
        return(gene_info,abundance_info)

def cosort(genes,gene_info,abundance_info):
    locs = [i for i in range(len(gene_info)) if gene_info[i] in genes]
    abund = [abundance_info[i] for i in locs]
    genes = [gene_info[i] for i in locs]
    zipped_dat = zip(genes,abund)
    sorted_pairs = sorted(zipped_dat,reverse=True,key = lambda x: x[1])
    tuples = zip(*sorted_pairs)
    genes, abund = [list(tuple) for tuple in  tuples]
    return(genes,abund)

def get_top_alleles(gene_info,abundance_info):
    
    class_1 = ["HLA-A","HLA-B","HLA-C"]
    class_2_A = ["HLA-DPA1","HLA-DQA1","HLA-DQA2","HLA-DRA"]
    class_2_DPB = ["HLA-DPB1","HLA-DPB2"]
    class_2_DQB = ["HLA-DQB1"]
    class_2_DRB = ["HLA-DRB1","HLA-DRB5"]
    
    # finding top class 1 genes
    genes, abund = cosort(class_1,gene_info,abundance_info)
    class_1_pair = genes[0:2]
    
    # finding class 2 A gene
    genes, abund = cosort(class_2_A,gene_info,abundance_info)
    class_2_A = genes[0]
    
    if "DP" in class_2_A:
        genes,abund = cosort(class_2_DPB,gene_info,abundance_info)
        class_2_B = genes[0]
    elif "DQ" in class_2_A:
        genes,abund = cosort(class_2_DQB,gene_info,abundance_info)
        class_2_B = genes[0]
    elif "DR" in class_2_A:
        genes,abund = cosort(class_2_DRB,gene_info,abundance_info)
        class_2_B = genes[0]
    targets = class_1_pair + [class_2_A,class_2_B]
    return(targets)
    
def main(options, args):
    genotype_log_file = options.genotype_log_file
    top_alleles = options.top_alleles
    IN_ALLELES = False
    with open(genotype_log_file,"r") as log_file:
        log_file_contents = log_file.read().splitlines()
        gene_info,abundance_info = get_abundance(log_file_contents)
        targets = get_top_alleles(gene_info,abundance_info)
        
        
if __name__ == "__main__":

    from optparse import OptionParser

    parser = OptionParser()

    parser.add_option("-g", "--genotype_log_file", dest="genotype_log_file",
                      help="the genotype log file")
    parser.add_option("-t", "--top_alleles", dest="top_alleles",
                      help="The number of top alleles to extract",
                      type=int)
    (options, args) = parser.parse_args()

    main(options,args)
