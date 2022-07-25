#!/usr/bin/env python

# This file creates a groups file based on the tcga file directory

import os


def main(options):
    junc_file = options.junc_file
    out = os.path.dirname(junc_file)
    normals = ""
    tumors = ""
    controls = ""
    samples_tum = ""
    samples_norm = ""
    with open(junc_file, "r") as juncs:
        junc_lines = juncs.read().splitlines()
        for sample in junc_lines:
            sample_no_junc = sample.replace(".junc", "")
            sample_split = sample_no_junc.split("-")
            sample_val = int(sample_split[3][0:2])
            if 1 <= sample_val <= 9:
                tumors += sample_no_junc + "\t" + "T" + "\n"
                samples_tum += "%s/%s\n" % (out, sample)
            elif 10 <= sample_val <= 19:
                normals += sample_no_junc + "\t" + "N" + "\n"
                samples_norm += "%s/%s\n" % (out, sample)
            else:
                controls += sample + "\t" + "N" + "\n"
    with open("%s/groups_file.txt"%out, "w") as group:
        groups = normals+tumors
        group.write(groups[:-1])
    with open("%s/junc_files.txt"%out, "w") as j:
        samples = samples_norm + samples_tum
        j.write(samples[:-1])

if __name__ == "__main__":
    from optparse import OptionParser

    parser = OptionParser()

    parser.add_option("-j", "--junc_file", dest="junc_file",
                      help="the junc file")

    (options, args) = parser.parse_args()

    main(options)