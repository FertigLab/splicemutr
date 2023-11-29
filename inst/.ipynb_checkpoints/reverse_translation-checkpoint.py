


def main(options, args):
    peptides_file = options.peptides_file
    output_file = options.output_file
    

if __name__ == "__main__":

    from optparse import OptionParser

    parser = OptionParser()
    
    parser.add_option("-p","--peptides_file",dest="peptides_file",
                     help="the peptides_file")
    parser.add_option("-o","--output_file",dest="output_file",
                 help="the output fasta file")
    (options, args) = parser.parse_args()

    main(options, args)