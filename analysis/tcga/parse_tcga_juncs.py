#!/usr/bin/env python

def main(options):
    junc_file = options.junc_file
    out_dir = options.out

    with open(junc_file,"r") as juncs:
        junc_lines = juncs.read().splitlines()
    samples = junc_lines[0]
    samples = samples.split()[2:len(samples)]
    bed_template = "%s\t%s\t%s\t%s\t%s\t%s\n"
    for i in range(len(samples)):
        sample_file = "%s/%s.junc"%(out_dir,samples[i])
        with open(sample_file,'w') as curr_sample:
            for j in range(2,len(junc_lines)):
                curr_line = junc_lines[j].split()
                junc = curr_line[0].split(",")
                junc_start = junc[0].split(":")
                junc_end = junc[1].split(":")
                if junc_start[2] == junc_end[2]:
                    curr_sample.write(bed_template%(junc_start[0],
                                                    junc_start[1],
                                                    junc_end[1],
                                                    "JUNC%d"%j,
                                                    curr_line[i+1],
                                                    junc_start[2]))
if __name__ == "__main__":
    from optparse import OptionParser

    parser = OptionParser()

    parser.add_option("-j", "--junc_file", dest= "junc_file",
                      help= "the junc file")
    parser.add_option("-o", "--out", dest="out",
                      help= "the output directory")

    (options, args) = parser.parse_args()

    main(options)