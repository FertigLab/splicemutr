#!/usr/bin/env python

import numpy as np

def calculate_similarity(counts_list):
    counts_np = np.array(counts_list)
    k = 4.87
    a = np.array(26)
    alpha = np.sum(np.exp(-1*k*(a - counts_list)))
    R = alpha/(1+alpha)
    return(float(R))
    
    

def main(option,args,b62):
    pep_file = args[0]
    IEDB_peps = args[1]
    output = options.output
    with open(pep_file) as file:
        peps = file.read().splitlines()
    with open(IEDB_peps) as file:
        peps_IEDB = file.read().splitlines()
        peps_IEDB = [pep.upper() for pep in peps_IEDB]
    
    similarity_IEDB =[]
    with open(output,"w") as out:
        for pep in peps:
            alignment_scores = np.array([sum([b62[(pep[i],pep_IEDB[i])] for i in range(9)]) for pep_IEDB in peps_IEDB])
            similarity_IEDB=(calculate_similarity(alignment_scores))
            out_string = "{peptide:s}\t{sim:.10e}"
            print(out_string.format(peptide=pep,sim=similarity_IEDB),file=out)
            
if __name__ == "__main__":
    
    from optparse import OptionParser
    from Bio.Align import substitution_matrices
    
    # creating the BLOSUM62 matrix for later
    b62 = substitution_matrices.load("BLOSUM62")
    
    parser = OptionParser()
    
    parser.add_option("-o", "--output", dest= "output",
                  help= "the output file")

    (options, args) = parser.parse_args()

    main(options,args,b62)
