#!/usr/bin/env python

def main(options):
    binders_file=options.binders
    pickle_dir=options.pickle_dir

if __name__ == "__main__":

    from optparse import OptionParser

    parser = OptionParser()

    parser.add_option("-t", "--top_header", dest="top_header",
                      help="T or F whether to include header")

    (options, args) = parser.parse_args()

    main(options, args)