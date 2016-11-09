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
# bedIntersect.py -w WINDOW_SIZE -i <path_to_file.out> -o <path_to_out.bed>
#
# Author: Martin Bontrager
############################################################

def main():
    window = 2500
    
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hw:i:o:",
                                   ["input=", "output="])
    except getopt.GetoptError:
        print('bedIntersect.py -w WINDOW_SIZE -i <path_to_file.out> -o <path_to_out.bed>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('bedIntersect.py -w WINDOW_SIZE -i <path_to_file.out> -o <path_to_out.bed>')
            sys.exit()
        elif opt in ("-w", "--window"):
            window = int(arg)
        elif opt in ("-i", "--input"):
            infile = arg
        elif opt in ("-o", "--output"):
            outfile = arg

    create_bed(infile, outfile, window)


def create_bed(infile, outfile, window):
    # Create a bed format file from the popoolation output file. Regions are
    # 5kb so I add and subtract 2500 from each position
    
    with open(infile) as fr, open(outfile, 'w') as fw:
        reader = csv.reader(fr, delimiter=' ')
        writer = csv.writer(fw, delimiter='\t')
        for row in reader:
            newrow = [None] * 3
            newrow[0] = row[0].strip('"')
            newrow[2] = int(row[1]) + window
            if (int(row[1]) - window) < 0:
                newrow[1] = 0
            else:
                newrow[1] = int(row[1]) - window
            writer.writerow(newrow)

if __name__ == '__main__':
    main()
