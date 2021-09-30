#!/usr/bin/env python

# Created: 08/03/2021
# Created by: Theron Palmer

import numpy as np
from bisect import bisect

def main(options, args):
    data_file = args[0]
    ref_file = args[1]
    output_file = args[2]
    
    with open(data_file) as dat:
        data_list = dat.read().splitlines()
        data_list = data_list[1:]
    data_vals = [float(i.split(',')[1]) for i in data_list if i.split(',')[1] != "ic50"]
    data_kmers = [i.split(',')[0] for i in data_list if i.split(',')[0] != "peptide"]
    with open(ref_file) as ref:
        ref_list = ref.read().splitlines()
        ref_list = ref_list[1:]
    ref_vals = [float(i.split(',')[1]) for i in ref_list if i.split(',') != "ic50"]
    
    ref_vals.sort()
    data_perc = [(bisect(ref_vals,i)/len(ref_vals)) for i in data_vals]
    with open(output_file,'w') as out: 
        percentile_data = ""
        for i in range(len(data_kmers)):
            percentile_data += "%s,%.10f\n"%(data_kmers[i],data_perc[i])
        print(percentile_data,file=out)
    
if __name__ == "__main__":

    from optparse import OptionParser

    parser = OptionParser()

    (options, args) = parser.parse_args()

    main(options, args)
