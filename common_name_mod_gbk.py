#!/usr/bin/python

# Script to convert product names in GBK file to a different set prescribed by metadata.

import sys, re


# Function to print the product name to a file. Accepts the product name
# and the handle for the outfile. Will process multi-line name if the
# characters-per-line exceeds 80.
def printName(name,out): # need some handling in case of multi-line
    base = " " * 21 + "/product="
    if len(name) < 48: # good to go, can print on one line
        final = '%s"%s"\n' % (base,name)
        out.write(final)
    else: # need to do a multi-line print
        first_line = '%s"%s-\n' % (base,name[0:48]) 
        out.write(first_line)
        blank_base = " " * 21
        beg = 48
        end = len(name)+1 # end range is exclusive
        # Note that name variable should be maintained to build multiple subsets from.
        new_name = name[i:end] # first subset, might be 2 or more lines
        while len(name) > 58: # max length to add spaces and end quotes
            name = name[i:end]
        final_line = '%s%s' % (blank_base,name)

# Function to build a name from a multi-line product in a GBK file. 
# Accepts the list of lines consisting of the product names, which
# will be condensed into a single str value, and the map provided
# by the metadata file to change a name to a common name.
def resolveName(name_list,name_dict):
    last_index = len(name_list) - 1 # remove surrounding quotes for checking dict
    name_list[0] = name_list[0][1:]
    name_list[last_index] = name_list[last_index][:-1]
    name = ' '.join(name_list)
    if name in name_dict: # use mapped value if its present
        name = name_dict[name]
    return name


metadata = str(sys.argv[1])

md = open(metadata,'r')

regex_for_product = '\s+/product="(.*)"'

# Iterate over the metadata file, one line per GBK to process
for line in md:
    line = line.strip('\n') 
    md_vals = line.split('\t')

    gbk_in = "./%s.locus_mod.gbk" % (md_vals[1][:-1])
    gbk = open(gbk_in,'r') # pull input GBK
    common_names = open(md_vals[2],'r') # pull name map file

    name_map = {}

    for name_pair in common_names: # build a map from the common names for conversion
        name_pair = name_pair.strip('\n')
        names = name_pair.split('\t')
        name_map[names[1]] = names[0]

    # Make a new file name for the outfile
    gbk_out = "./%s.common_name_mod.gbk" % (md_vals[1][:-1])
    outfile = open(gbk_out,'w')

    found_product = False
    multi_product = [] # need a list for multi-line product

    for l in gbk:

        if "/product=" in l: # found the beginning of product tag
            if re.search(regex_for_product,l): # simplest case, product on one line
                single_product = re.search(regex_for_product,l).group(1) # single line product
                if single_product in name_map: # if in map, replace with mapped value
                    printName(name_map[single_product],outfile)
                else: # leave as is
                    outfile.write(l)  
            else: # multi-line, build product list
                l = l.strip('\n') # strip EOL and then leading whitespace
                l = l.lstrip()
                multi_product.append(l[9:]) 
                found_product = True

        elif found_product == True: # multi-line product
            if "/translation=" in l: # reached the end of the product, resolve it and leave
                printName(resolveName(multi_product,name_map),outfile)
                outfile.write(l)
                found_product = False # reinitialize all for next multi-line product
                multi_product = []
            else:
                l = l.strip('\n')
                multi_product.append(l.lstrip()) # strip leading white space

        else:
            outfile.write(l)