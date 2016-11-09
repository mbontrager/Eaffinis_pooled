#!/usr/bin/python
# Script from https://code.google.com/p/bioman/source/browse/BlastXMLmerge.py
# Updated by Martin Bontrager for Condor output.
import sys
import subprocess
import os
import glob

path = os.path.dirname(os.path.abspath(__file__))

# f = glob.glob(path + "/*.tar.gz")
# subprocess.call(["mkdir", "unzipped"])
#subprocess.call(["rm", "cmposit.tar.gz"])

# for a in f:
    # unzip_cmd = 'tar -zxvf ' + a + ' -C unzipped'
    # print unzip_cmd
    # p = subprocess.Popen(unzip_cmd , shell=True)
    # p.communicate()

xmlfiles = glob.glob(path + "/*block*")
xmlout = path + "/BLAST_results.xml"

def merge(split_files, output_file):
        """Merging multiple XML files is non-trivial and must be done in subclasses."""
        if len(split_files) == 1:
            #For one file only, use base class method (move/copy)
            return Text.merge(split_files, output_file)
        out = open(output_file, "w")
        h = None
        for f in split_files:
            h = open(f)
            body = False
            header = h.readline()
            if not header:
                out.close()
                h.close()
                raise ValueError("BLAST XML file %s was empty" % f)
            if header.strip() != '<?xml version="1.0"?>':
                out.write(header) #for diagnosis
                out.close()
                h.close()
                raise ValueError("%s is not an XML file!" % f)
            line = h.readline()
            header += line
            if line.strip() not in ['<!DOCTYPE BlastOutput PUBLIC "-//NCBI//NCBI BlastOutput/EN" "http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd">',
                                    '<!DOCTYPE BlastOutput PUBLIC "-//NCBI//NCBI BlastOutput/EN" "NCBI_BlastOutput.dtd">']:
                out.write(header) #for diagnosis
                out.close()
                h.close()
                raise ValueError("%s is not a BLAST XML file!" % f)
            while True:
                line = h.readline()
                if not line:
                    out.write(header) #for diagnosis
                    out.close()
                    h.close()
                    raise ValueError("BLAST XML file %s ended prematurely" % f)
                header += line
                if "<Iteration>" in line:
                    break
                if len(header) > 10000:
                    #Something has gone wrong, don't load too much into memory!
                    #Write what we have to the merged file for diagnostics
                    out.write(header)
                    out.close()
                    h.close()
                    raise ValueError("BLAST XML file %s has too long a header!" % f)
            if "<BlastOutput>" not in header:
                out.close()
                h.close()
                raise ValueError("%s is not a BLAST XML file:\n%s\n..." % (f, header))
            if f == split_files[0]:
                out.write(header)
                old_header = header
            elif old_header[:300] != header[:300]:
                #Enough to check <BlastOutput_program> and <BlastOutput_version> match
                out.close()
                h.close()
                raise ValueError("BLAST XML headers don't match for %s and %s - have:\n%s\n...\n\nAnd:\n%s\n...\n" \
                                 % (split_files[0], f, old_header[:300], header[:300]))
            else:
                out.write("    <Iteration>\n")
            for line in h:
                if "</BlastOutput_iterations>" in line:
                    break
                #TODO - Increment <Iteration_iter-num> and if required automatic query names
                #like <Iteration_query-ID>Query_3</Iteration_query-ID> to be increasing?
                out.write(line)
            h.close()
        out.write("  </BlastOutput_iterations>\n")
        out.write("</BlastOutput>\n")
        out.close()

merge(xmlfiles,xmlout)
subprocess.call(["rm", "-r", "unzipped/"])
