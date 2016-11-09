#!/usr/bin/python

# Script to standardize "hypothetical" product name syntax and  tags in GBK file
# to just be "hypothetical protein". This will also remove any gene symbols that 
# precede a hypothetical protein if they are present. 

import sys, re

metadata = str(sys.argv[1])

md = open(metadata,'r')

# Iterate over the metadata file, one line per GBK to process
for line in md:
    line = line.strip('\n') 
    md_vals = line.split('\t')

    gbk_in = "./%s.common_name_mod.gbk" % (md_vals[1][:-1])
    gbk = open(gbk_in,'r') # pull input GBK

    # Make a new file name for the outfile
    gbk_out = "./%s.hypothetical_mod.gbk" % (md_vals[1][:-1])
    outfile = open(gbk_out,'w')

    for l in gbk:
        if "/product=" in l:
            l = l.replace("conserved hypothetical protein","hypothetical protein")
            l = l.replace("conserved domain protein","hypothetical protein")
            l = l.replace("conserved protein","hypothetical protein")
            outfile.write(l)
        else:
            outfile.write(l)
