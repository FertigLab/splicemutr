#!/usr/bin/env python
# Theron Palmer
# Creating a kmer_transcript map for human proteome then create epitope transcript map using epitope database, kmer size 9

from collections import defaultdict
from Bio import SeqIO
import pandas as pd

def kmerize(pep,kmer_length):
    # kmerize() : kmerize a peptide
    # inputs:
     # pep: the peptide to be kmerized
     # kmer_length: the kmer length to kmerize the peptide over
    # outputs: 
     # the list of kmers from the peptide
    pep_len=len(pep)
    if len(pep) < kmer_length:
        pep += "Z"*(kmer_length-pep_len)
        pep_len=len(pep)
    return([pep[i:i+kmer_length] for i in range(pep_len-kmer_length+1)])


def parse_fasta(fasta_file,kmer_size):
    fasta_dict = defaultdict(list)
    fasta_sequences = SeqIO.parse(open(fasta_file),'fasta')
    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
        transcript = name.split("|")[1]
        kmers = kmerize(sequence,kmer_size)
        for kmer in kmers:
            fasta_dict[kmer].append(transcript)
    return(fasta_dict)

def main(options,args):
    gencode_proteins_fasta = options.gencode_fasta
    epitopes_in_cancer = options.epitopes_file
    output_file = options.output_file

    kmer_transcript_dict = parse_fasta(gencode_proteins_fasta,9)

    with open(epitopes_in_cancer,"r") as epitopes:
        epitope_transcript_dict = defaultdict(str)
        for line in epitopes.readlines():
            line = line.strip()
            if line in kmer_transcript_dict.keys():
                epitope_transcript_dict[line]=";".join(kmer_transcript_dict[line])
            else:
                epitope_transcript_dict[line]="none"

    epitope_transcript_df = pd.DataFrame.from_dict(epitope_transcript_dict,orient="index")

    epitope_transcript_df.to_csv(output_file,sep="\t",header=False)

if __name__ == "__main__":

    from optparse import OptionParser

    parser = OptionParser()
    
    parser.add_option("-g","--gencode_fasta",dest="gencode_fasta",
                     help="the gencode fasta file")
    parser.add_option("-e","--epitopes_file",dest="epitopes_file",
                 help="the input epitopes file")
    parser.add_option("-o","--output_file",dest="output_file",
                help="the output file")
    (options, args) = parser.parse_args()

    main(options, args)