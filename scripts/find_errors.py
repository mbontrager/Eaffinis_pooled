#!/usr/bin/python

import os, sys, getopt, glob, re

############################################################
# Find condor GATK jobs that terminated due to errors
#
# Input - path to condor .log files
#
# Usage: 
# find_errors.py -p <err_path>
#
# Author: Martin Bontrager
############################################################


def main():
    
    try:
        opts, args = getopt.getopt(sys.argv[1:],"hp:",["path="])
    except getopt.GetoptError:
        print('find_errors.py -p <err_path>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('find_errors.py -p <err_path>')
            sys.exit()
        elif opt in ("-p", "--path"):
            path = arg
                   
    find_errors(get_files(path))
    
# Get a list of all the .err log files in a directory
def get_files(path):
    
    os.chdir(path)
    full = glob.glob('*.err')
    files = []
    for f in full:
        if os.path.getsize(path + '/' + f) == 0:
            continue
        else:
            files.append(f)
    return(files)

# find GATK error messages in a list of log files
def find_errors(files):
    errors = []
    for i in files:
        if '##### ERROR' in open(i).read():
            errors.append(i.replace('GATK_HC_multi_', '').replace('.err', ''))

    print("\n".join(errors))

if __name__ == '__main__':
  main()
