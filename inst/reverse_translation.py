#!/usr/bin/env python
2
# creator: Theron Palmer
# Reverse translating input peptide files and saving them to fasta files for use later

CODON_TABLE = {
	'A': ('GCT', 'GCC', 'GCA', 'GCG'),
	'C': ('TGT', 'TGC'),
	'D': ('GAT', 'GAC'),
	'E': ('GAA', 'GAG'),
	'F': ('TTT', 'TTC'),
	'G': ('GGT', 'GGC', 'GGA', 'GGG'),
	'I': ('ATT', 'ATC', 'ATA'),
	'H': ('CAT', 'CAC'),
	'K': ('AAA', 'AAG'),
	'L': ('TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'),
	'M': ('ATG',),
	'N': ('AAT', 'AAC'),
	'P': ('CCT', 'CCC', 'CCA', 'CCG'),
	'Q': ('CAA', 'CAG'),
	'R': ('CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'),
	'S': ('TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'),
	'T': ('ACT', 'ACC', 'ACA', 'ACG'),
	'V': ('GTT', 'GTC', 'GTA', 'GTG'),
	'W': ('TGG',),
	'Y': ('TAT', 'TAC'),
	'*': ('TAA', 'TAG', 'TGA'),
}

def reverse_translation(peptide,out_file):
	list_aa = list(peptide)
	count_to_return = 1
	to_print = ''

	for aa in list_aa:
		count_to_return *= len(CODON_TABLE[aa])

	else:
		sequences = [codon for codon in CODON_TABLE[peptide[0]]]
		count_to_return = 0
		count_to_return_aux = 0
		
		for index, amino_acid in enumerate(peptide[1:]):
			to_extend = sequences
			sequences = []
			for index_codon, codon in enumerate(CODON_TABLE[amino_acid]):
				for index_seq,sequence in enumerate(to_extend):
					sequence += codon
					sequences.append(sequence)

					if index == len(peptide[1:])-1 :
						count_to_return += 1
						count_to_return_aux += 1
						to_print+= '>'+peptide+'\n'+sequence+'\n'

						if count_to_return_aux == 10000000:
							count_to_return_aux = 0
							sequences = []

	out_file.write(to_print)
	to_print = ''
	sequences = []
	return count_to_return, peptide

def reverse_translation_peptides(peptide_file,output_file):
	peptides=list()
	with open(peptide_file,"r") as file:
		for line in file.readlines():
			peptides.append(line.strip())
	with open(output_file,"w") as out_file:
		for peptide in peptides:
			reverse_translation(peptide,out_file)
	return(True)


def main(options, args):
	peptides_file = options.peptides_file
	output_file = options.output_file
	reverse_translation_peptides(peptides_file,output_file)
    

if __name__ == "__main__":

    from optparse import OptionParser

    parser = OptionParser()
    
    parser.add_option("-p","--peptides_file",dest="peptides_file",
                     help="the peptides_file")
    parser.add_option("-o","--output_file",dest="output_file",
                 help="the output fasta file")
    (options, args) = parser.parse_args()

    main(options, args)