#!/usr/bin/python

# This script is the first of many which accepts a metadata file with one row for 
# each GBK to format and the following 24 columns present in each row:
# 
# gbk/file | locus tag | curate/common/names/file | delete/gene/symbols/file | 
# bioproject ID | organism name | strain name | serotype | host | date of isolation
# country | assembly method | coverage | sequencing method | contact person | 
# email of contacts | authors | title | illegal/ec/numbers/file | isolation source
# contact/list/file | biosample ID
# 
# Note that values can be null, but you must add a tab even for missing values.
#
# Script to convert locus tags in GBK file to a different set prescribed by metadata.

import sys, re

metadata = str(sys.argv[1])

md = open(metadata,'r')

regex_for_locus = r'\s+/locus_tag="(.*_)\d+"'

# Iterate over the metadata file, one line per GBK to process
for line in md:
    line = line.strip('\n') 
    md_vals = line.split('\t')

    gbk = open(md_vals[0],'r') # pull input GBK
    locus = md_vals[1] # pull locus prefix

    # Make a new file name for the outfile
    gbk_out = "./%s.locus_mod.gbk" % (md_vals[1][:-1]) # use the locus as GBK file name
    outfile = open(gbk_out,'w')

    for l in gbk:
        if "/locus_tag=" in l:
            old_tag = re.search(regex_for_locus,l).group(1)
            l = l.replace(old_tag,locus)
            outfile.write(l)
        else:
            outfile.write(l)
