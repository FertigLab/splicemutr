#!/usr/bin/env python

import splicemutr as sp
import os

def main(options):
    peptides_file = options.peptides
    kmer_length=int(options.kmer_length)
    
    with open(peptides_file) as peps:
        pep_dat=peps.read().replace("*","").splitlines()
    
    pep_dat_vals=[None]*len(pep_dat)
    for i in range(len(pep_dat)):
        row = pep_dat[i]
        row_vals = row.split("\t")
        if row_vals[1]=="NA":
            row_vals.append(str(0))
            row_vals.append(str(1))
            pep_dat_vals[i]="\t".join(row_vals)
        elif "-" in row_vals[0]:
            splicemutr_tx=set(sp.kmerize(row_vals[1],kmer_length))
            ref_txs = row_vals[2].split("-")
            ref_tx_1=set(sp.kmerize(ref_txs[0],kmer_length))
            ref_tx_2=set(sp.kmerize(ref_txs[1],kmer_length))
            ref_tx = set(ref_tx_1.union(ref_tx_2))
            row_vals.append(str(len(splicemutr_tx.intersection(ref_tx))/len(ref_tx)))
            row_vals.append(str(len(splicemutr_tx.difference(ref_tx))/len(ref_tx)))
            pep_dat_vals[i]="\t".join(row_vals)
        else:
            splicemutr_tx=set(sp.kmerize(row_vals[1],kmer_length))
            ref_tx=set(sp.kmerize(row_vals[2],kmer_length))
            row_vals.append(str(len(splicemutr_tx.intersection(ref_tx))/len(ref_tx)))
            row_vals.append(str(len(splicemutr_tx.difference(ref_tx))/len(ref_tx)))
            pep_dat_vals[i]="\t".join(row_vals)
            
    out_file=os.path.dirname(peptides_file)+"/"+os.path.basename(peptides_file).replace(".txt","_sim.txt")
    with open(out_file,"w") as out:
        print("\n".join(pep_dat_vals),file=out)
        
if __name__ == "__main__":
    
    from optparse import OptionParser

    parser = OptionParser()

    parser.add_option("-p", "--peptides", dest="peptides",
                      help="the transcript peptide data file")
    parser.add_option("-k","--kmer_length",dest="kmer_length",
                      help="the kmer_length for kmerizing")
    (options, args) = parser.parse_args()

    main(options)