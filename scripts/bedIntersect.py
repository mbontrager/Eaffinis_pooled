#!/usr/bin/python

import csv
import sys
import getopt

############################################################
# Give list of regions from popoolation, prep Bedtools intersect
#
# Input - path to .out file
#
# Usage:
# bedIntersect.py -i <path_to_file.out> -o <path_to_out.bed>
#
# Author: Martin Bontrager
############################################################


def main():
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hi:o:",
                                   ["input=", "output="])
    except getopt.GetoptError:
        print('bedIntersect.py -i <path_to_file.out> -o <path_to_out.bed>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('bedIntersect.py -i <path_to_file.out> -o <path_to_out.bed>')
            sys.exit()
        elif opt in ("-i", "--input"):
            infile = arg
        elif opt in ("-o", "--output"):
            outfile = arg

    create_bed(infile, outfile)


def create_bed(infile, outfile):
    # Create a bed format file from the popoolation output file. Regions are
    # 5kb so I add and subtract 2500 from each position
    
    with open(infile) as fr, open(outfile, 'w') as fw:
        reader = csv.reader(fr, delimiter=' ')
        writer = csv.writer(fw, delimiter='\t')
        for row in reader:
            newrow = [None] * 3
            newrow[0] = row[0].strip('"')
            newrow[2] = int(row[1]) + 2500
            newrow[1] = int(row[1]) - 2500
            writer.writerow(newrow)

if __name__ == '__main__':
    main()
