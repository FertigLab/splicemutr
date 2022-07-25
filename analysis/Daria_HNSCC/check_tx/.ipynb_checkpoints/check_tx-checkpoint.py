# !/usr/bin/env python

# This script checks whether a list of ensemble ids (transcripts) exist in a gtf

try:
    import argparse
    import sys
    import os
except:
    print("You do not have the correct python modules installed", file=sys.stderr)
    sys.exit()

parser = argparse.ArgumentParser(description='Inputting the sequences')
parser.add_argument("-g", '--gtf', action="store", nargs=1, help="The genome gtf file", dest='gtf')
parser.add_argument("-t", "--trans", action="store", nargs=1, help="tx ids path and file", dest='trans')

results = parser.parse_args()

gtf = results.gtf[0]
trans = results.trans[0]

trans_mod = ""
with(open(gtf, "r")) as reference:
    ref = reference.read().replace('\n', '')
    with (open(trans, "r")) as tx:
        tx_line = tx.readline().strip('\n')
        while tx_line != "":
            if tx_line in ref:
                trans_mod = trans_mod + tx_line + "\t+\n"
            else:
                trans_mod = trans_mod + tx_line + "\t-\n"
            tx_line = tx.readline().strip('\n')

out_file = trans[0:-4]+"_check.txt"
trans_mod_file = open(out_file, "w")
trans_mod_file.write(trans_mod)
trans_mod_file.close()
