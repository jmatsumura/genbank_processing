#!/usr/bin/python

# Script to standardize "hypothetical" product name syntax and  tags in GBK file
# to just be "hypothetical protein". This will also remove any gene symbols that 
# precede a hypothetical protein if they are present. 

import sys, re

metadata = str(sys.argv[1])

md = open(metadata,'r')

hypothetical = ["conserved hypothetical protein","conserved domain protein"
    ,"conserved protein","conserved hypothetical family protein"
    ,"hypothetical protein"]
    
hypos = set(hypothetical)

regex_for_product = r'\s+/product="(.*)"'

# Iterate over the metadata file, one line per GBK to process
for line in md:
    line = line.strip('\n') 
    md_vals = line.split('\t')

    gbk_in = "./%s.common_name_mod.gbk" % (md_vals[1][:-1])
    gbk = open(gbk_in,'r') # pull input GBK

    # Make a new file name for the outfile
    gbk_out = "./%s.hypothetical_mod.gbk" % (md_vals[1][:-1])
    outfile = open(gbk_out,'w')

    between_gene_and_product = False
    hypothetical_found = False
    g_to_p_list = []

    for l in gbk:

        if between_gene_and_product == True:
            if "/product=" in l:
                if re.search(regex_for_product,l): # ignore multi-line products
                    if re.search(regex_for_product,l).group(1) in hypos:
                        hypothetical_found = True
                        l = l.replace("conserved hypothetical family protein","hypothetical protein")
                        l = l.replace("conserved hypothetical protein","hypothetical protein")
                        l = l.replace("conserved domain protein","hypothetical protein")
                        l = l.replace("conserved protein","hypothetical protein")
                g_to_p_list.append(l) # add product (potentially modified)

                for x in g_to_p_list:
                    if "/gene=" in x and hypothetical_found == True: # ignore gene symbols if hypothetical
                        pass
                    else:
                        outfile.write(x)

                g_to_p_list = [] # reinitialize
                hypothetical_found = False
                between_gene_and_product = False

            else: # still inbetween, keep building the list
                g_to_p_list.append(l)

        elif "   gene   " in l:
            g_to_p_list.append(l)
            between_gene_and_product = True

        else:
            outfile.write(l)
