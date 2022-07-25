#!/usr/bin/env python

def main(options):
    junc_file = options.junc_file
    cov_file = options.cov_file
    samples = options.samples
    out_dir = options.out
    
    
    
if __name__ == "__main__":
    from optparse import OptionParser

    parser = OptionParser()

    parser.add_option("-b", "--junc_file", dest= "junc_file",
                      help= "the junc file")
    parser.add_option("-c","--cov_file",dest="cov_file")
    parser.add_option("-s","--samples",dest="samples")
    parser.add_option("-o", "--out", dest="out",
                      help= "the topmost output directory")

    (options, args) = parser.parse_args()

    main(options)