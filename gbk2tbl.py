#!/usr/bin/python

# Script to convert from GBK to TBL. 

import sys, re
from Bio import SeqIO

metadata = str(sys.argv[1])

md = open(metadata,'r')

# Iterate over the metadata file, one line per GBK to process
for line in md:
    line = line.strip('\n') 
    md_vals = line.split('\t')

    gbk_in = "./%s.ec_numbers_mod.gbk" % (md_vals[1][:-1])
    gbk = open(gbk_in,'rU') # pull input GBK

    # Make a new file name for the outfile
    tbl_out = "./%s.tbl" % (md_vals[1][:-1])
    outfile = open(tbl_out,'w')

    records = SeqIO.parse(gbk, 'genbank')
    output_handle = open("./cor6_6.fasta", "w")
    count = SeqIO.write(records, output_handle, "fasta")
