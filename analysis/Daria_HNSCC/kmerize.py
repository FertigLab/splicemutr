import sys

try:
    import pickle
    import argparse
except:
    print("You do not have the correct python modules installed", file=sys.stderr)
    sys.exit()

parser = argparse.ArgumentParser(description='Inputting the sequences')
parser.add_argument('-i', '--input', action="store", nargs=1, help='The input metadata file', dest="input")
parser.add_argument('-o', '--output', action="store", nargs=1, help='The output kmer file', dest="output")
results = parser.parse_args()


def save_dict(dictionary, file):
    """
    using cPickle to save a dictionary as a pickled binary file
    inputs:
        dictionary - a dictionary object
        File - a filename string with extension for saving the pickled dictionary
    """
    with open(file, "wb") as my_file:
        pickle.dump(dictionary, my_file)
        my_file.close()


def load_dict(file):
    """
    using cPickle to load a binary pickled dictionary from a file
    inputs:
        File - a filename string with extension for the pickled dicitonary
    outputs:
        dictionary - a dictionary object
    """
    with open(file, "rb") as my_file:
        dictionary = pickle.load(my_file)
        my_file.close()
        return dictionary


def main():
    chr = 3
    gene = 2
    junc_start = 4
    junc_end = 5
    verdict = 6
    with open(results.input[0], 'r') as f:
        line = f.readline()
        line_items = line.split('\t')
        # format of junction string is chr:gene:junc_start:junc:end:verdict
        junc_prompt = "{0}:{1}:{2}:{3}:{4}"
        junc_txt = junc_prompt.format(chr, gene, junc_start, junc_end, verdict)
