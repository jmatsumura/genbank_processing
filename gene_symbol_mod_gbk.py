#!/usr/bin/python

# Script to remove (or replace) erroneous gene symbols from a GBK file. 

import sys, re

metadata = str(sys.argv[1])

md = open(metadata,'r')

regex_for_gene_sym = r'\s+/gene="(.*)"'

# Iterate over the metadata file, one line per GBK to process
for line in md:
    line = line.strip('\n') 
    md_vals = line.split('\t')

    gbk_in = "./%s.hypothetical_mod.gbk" % (md_vals[1][:-1])
    gbk = open(gbk_in,'r') # pull input GBK
    gene_symbols = open(md_vals[3],'r') # pull name map file

    gs_map = {}

    for pair in gene_symbols: # build a map from the common names for conversion
        pair = pair.strip('\n')
        names = pair.split('\t')
        gs_map[names[1]] = names[0]

    # Make a new file name for the outfile
    gbk_out = "./%s.gene_symbol_mod.gbk" % (md_vals[1][:-1])
    outfile = open(gbk_out,'w')

    for l in gbk:

        if "/gene=" in l: # found the gene symbol
            if re.search(regex_for_gene_sym,l): 
                current_symbol = re.search(regex_for_gene_sym,l).group(1) 
                if current_symbol in gs_map: # if in map, replace with mapped value
                    if not gs_map[current_symbol]: # if null, ignore and print nothing
                        pass
                    else: # if not null, use whatever mapped value is present
                        l = l.replace(current_symbol,gs_map[current_symbol])
                        outfile.write(l)
                else: # leave as is
                    outfile.write(l)  

        else:
            outfile.write(l)
