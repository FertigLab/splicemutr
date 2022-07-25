# !/usr/bin/env python

try:
    import argparse
    import sys
    import os
except:
    print("You do not have the correct python modules installed", file=sys.stderr)
    sys.exit()

parser = argparse.ArgumentParser(description='Inputting the sequences')
parser.add_argument("-g", '--genome', action="store", nargs=1, help="The genome fasta file", dest='genome')
parser.add_argument("-o", "--output", action="store", nargs=1, help="Output directory for the chromosome fasta files",
                    dest="output")
results = parser.parse_args()

genome = results.genome[0]
out = results.output[0]
filename = out+"/chr"+"{}"+".fa"
first_past = False
chrom_string = "paste('chr', c({}), sep="")"
chrom = ""
with(open(genome, "r")) as fasta:
    line = fasta.readline()
    while line != "":
        if ">" in line:
            chromosome = line.split(" ")[1].strip()
            print(chromosome)
            fasta_dir = filename.format(chromosome)
            if ~os.path.exists(fasta_dir):
                if (first_past == True):
                    fasta_file.close()
                chrom = chrom + chromosome + ","
                fasta_file = open(fasta_dir, "w")
                first_past = True
        fasta_file.write(line)
        line = fasta.readline()
    print(chrom_string.format(chrom[0:(len(chrom)-1)]))


