#!/usr/bin/python

# Script to convert locus tags in GBK file to a different set prescribed by metadata.

import sys, re

metadata = str(sys.argv[1])
out_dir = str(sys.argv[2])

md = open(metadata,'r')

regex_for_locus = r'\s+/locus_tag="(.*_)\d+"'

# Iterate over the metadata file, one line per GBK to process
for line in md:
    line = line.strip('\n') 
    md_vals = line.split('\t')

    locus = md_vals[1] # pull locus prefix

    gbk_in = "%s/%s/cleaned.gbk" % (out_dir,md_vals[1][:-1])
    gbk = open(gbk_in,'r') # pull input GBK

    # Make a new file name for the outfile
    gbk_out = "%s/%s/locus_mod.gbk" % (out_dir,md_vals[1][:-1]) # use the locus as GBK file name
    outfile = open(gbk_out,'w')

    for l in gbk:
        if "/locus_tag=" in l:
            old_tag = re.search(regex_for_locus,l).group(1)
            l = l.replace(old_tag,locus)
            outfile.write(l)
        else:
            outfile.write(l)
