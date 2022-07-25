import sys

try:
    import pickle
except:
    print("You do not have the correct python modules installed", file=sys.stderr)
    sys.exit()
    

def main(options, args):
    genotype_log_file = args.genotype_log_file
    top_alleles = options.top_alleles
    IN_ALLELES = False
    with open(genotype_log_file,"r") as log_file:
        log_file_contents = log_file.read().splitlines()
        line_start = [i for i in range(len(log_file_contents)) if "Observed HLA genes" in log_file_contents[i]][0]
        line_end = [i for i in range(line_start,len(log_file_contents)) if "---" in log_file_contents[i]][0]
    print(line_start)
    print(line_end)
                
    
if __name__ == "__main__":

    from optparse import OptionParser

    parser = OptionParser()

    parser.add_argument("-g", "--genotype_log_file", dest="genotype_log_file",
                      help="the genotype log file")
    parser.add_option("-t", "--top_alleles", dest="top_alleles",
                      help="The number of top alleles to extract",
                      type=int)
    (options, args) = parser.parse_args()

    main(options,args)
