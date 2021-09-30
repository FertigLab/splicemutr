#!/usr/bin/env python

import os
import statistics as st

def main(options, args):
    
    perc_file_list = args[0]
    perc_count = ""
    with open(perc_file_list) as file_list:
        perc_files = file_list.read().splitlines()
    iterr = 1
    total = len(perc_files)
    for file in perc_files:
        print("%d out of %d"%(iterr,total))
        if os.stat(file).st_size == 0: continue
        allele = os.path.basename(file).split("_")[0]
        with open(file) as percs:
            perc_vals = percs.read().splitlines()
            perc_pre = [len(i.split("\t")[2].split(":")) for i in perc_vals]
            perc_count += "%s\t%0.3f\t%0.3f\t%d\n"%(allele,st.median(perc_pre),st.mean(perc_pre),len(perc_pre))
        iterr+=1
    out_dir = os.path.dirname(file)
    with open("%s/HLA_summ.txt"%(out_dir),'w') as out:
        print(perc_count[:-1],file=out)

if __name__ == "__main__":

    from optparse import OptionParser

    parser = OptionParser()

    (options, args) = parser.parse_args()

    main(options,args)