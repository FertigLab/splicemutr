#!/usr/bin/env python

# author: Theron Palmer
# created: 08/25/2021

# This scripts forms transcripts based on differentially-used junctions

import pyranges as py
import os
import sys
import pandas as pd

def main(options,args):
    
    # checking user input
    output_dir = options.output_dir
#     if not os.path.isdir(output_dir):
#         print("%s is not a directory"%output_dir,file=sys.stderr)
#         sys.exit()
        
    annotation_file = options.annotation_file
#     if not os.path.isfile(annotation_file):
#         print("%s is not an annotation file"%annotation_file,file=sys.stderr)
#         sys.exit()
#     if not ".gtf" in annotation_file:
#         print("%s is not a .gtf file")
#         sys.exit()
        
    junc_file = options.junc_file
#     if not os.path.isfile(junc_file):
#         print("%s is not a junc file"%junc_file)
#         sys.exit()
    
    # opening and reading files
    
    annotation_file = py.read_gtf(annotation_file)
    junc_file = pd.read_table(junc_file)
    
    all_exons = annotation_file["exon"]
    print(all_exons)
    return(True)

if __name__ == "__main__":
    
    from optparse import OptionParser

    parser = OptionParser()
    
    parser.add_option("-o","--output_dir",dest="output_dir",
                     help="The directory to write output to")
    parser.add_option("-g", "--annotation_file", dest="annotation_file",
                      help="the annotation file (gtf or gff)")
    parser.add_option("-j", "--junc_file", dest="junc_file",
                      help="the target junctions to modify transcripts with")

    (options, args) = parser.parse_args()

    main(options, args)