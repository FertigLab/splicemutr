#!/usr/bin/env python

# converts a single element per line tab file to fasta format.

def main(args):
    input_file = args[0]
    output_file = args[1]
    if not path.isfile(input_file):
        sys.exit()
    else:
        with open(output_file,'w') as out:
            line = ""
            with open(input_file) as file:
                contents = file.read().splitlines()
                for element in contents:
                    if "Z" in element:
                        continue
                    line += ">a\n%s\n"%(element)
            print(line[:-1],file=out)

if __name__ == "__main__":
    from optparse import OptionParser
    from os import path

    parser = OptionParser()
    
    (options, args) = parser.parse_args()
    
    main(args)