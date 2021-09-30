#!/usr/bin/env python

# author: Theron Palmer
# created: 09/07/2021

import os 

def main(options, args):
        
        txt_file = args[0]
        with open(txt_file) as f:
            proteins = f.read().splitlines()
        out_file = os.path.dirname(txt_file)+"/"+"proteins.fa"
        with open(out_file,"w") as out:
            protein_line = [">protein%d\n%s\n"%(i,proteins[i][:len(proteins[i])-1]) for i in range(len(proteins))]
            proteins_str = "".join(protein_line)
            print(proteins_str,file=out)
            
if __name__ == "__main__":

    from optparse import OptionParser

    parser = OptionParser()
    
    (options, args) = parser.parse_args()

    main(options, args)