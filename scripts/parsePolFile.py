#!/usr/bin/python

import csv
import sys
import getopt

############################################################
# Parse the output .pol file from MAPGD
#
# Input - path to .pol file
#
# Usage:
# parsePolFile.py -i <path_to_file.pol> -o <path_to_out.csv>
#
# Author: Martin Bontrager
############################################################

COVMAX = 800
COVMIN = 40


def main():
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hi:o:",
                                   ["input=", "output="])
    except getopt.GetoptError:
        print('parsePolFile.py -i <path_to_file.pol> -o <path_to_out.csv>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('parsePolFile.py -i <path_to_file.pol> -o <path_to_out.tsv>')
            sys.exit()
        elif opt in ("-i", "--input"):
            pol = arg
        elif opt in ("-o", "--output"):
            outfile = arg

    filter_variants(pol, outfile)

# Parse the .pol file
def filter_variants(pol, outfile):

    with open(pol) as fr, open(outfile, 'w') as fw:
        next(fr)
        reader = csv.reader(fr, delimiter='\t')
        writer = csv.writer(fw, delimiter='\t')
        headers = next(reader, None)  # returns the headers or None
        if headers:
            writer.writerow(headers)
        for row in reader:
            if (row[6].startswith('1') and
                  row[7].startswith('1') and
                  row[8].startswith('1') and
                  row[9].startswith('1')):
                continue
            if (int(row[4]) > COVMIN and int(row[4]) < COVMAX):
                writer.writerow(row)

if __name__ == '__main__':
    main()
