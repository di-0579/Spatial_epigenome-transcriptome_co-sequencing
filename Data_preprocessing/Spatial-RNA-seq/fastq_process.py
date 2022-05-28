from Bio.SeqIO.QualityIO import FastqGeneralIterator
from gzip import open as gzopen

import argparse

ap = argparse.ArgumentParser()
ap.add_argument("-i", "--input", required=True, help="Path to the inputfile")
ap.add_argument("-o", "--output", required=True, help="Path to the outputfile")
args = vars(ap.parse_args())

input_filename = args["input"]
output_filename = args["output"]

with gzopen(input_filename, "rt") as in_handle:
    with open(output_filename, "w") as out_handle:
        for title, seq, qual in FastqGeneralIterator(in_handle):
            new_seq = seq[22:30] + seq[60:68] + seq[98:108]  # BC2 + BC1 + UMI
            new_qual = qual[22:30] + qual[60:78] + qual[98:108]
            out_handle.write("@%s\n%s\n+\n%s\n" % (title, new_seq, new_qual))