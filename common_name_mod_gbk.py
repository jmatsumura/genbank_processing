

# Script to convert product names in GBK file to a different set prescribed by metadata.

import sys, re

# Function to print the product name to a file. Accepts the product name
# and the handle for the outfile. Will process multi-line name if the
# characters-per-line exceeds 80.
def printName(name,out): # need some handling in case of multi-line
    base = " " * 21 + "/product=" # first line
    blank_base = " " * 21 # non-first lines

    # Here replace some common errors found in product names
    name = name.replace('sulphide','sulfide')
    name = name.replace('sulphur','sulfur')
    name = name.replace(' fibre ',' fiber ')

    if len(name) < 48: # good to go, can print on one line
        final = '%s"%s"\n' % (base,name)
        out.write(final)

    # If greater than that length, need to do a multi-line print. Note
    # that this approach does not absolutely maximize line length, but 
    # arbitrarily does not permit this region to extend beyond the ~70 
    # character mark. 
    else: 

        max_len = 47
        words = name.split(' ')
        multi_line = [] # store all the lines in the product qualifier
        single_line = [] # store just the current line being built

        for x in words:
            if not (len(' '.join(single_line)) + len(x)) > max_len: # if i can fit, add it
                single_line.append(x)
            else: # too long, add to multi-line and build a new one
                multi_line.append(' '.join(single_line))
                single_line = []
                single_line.append(x)

        multi_line.append(' '.join(single_line))

        for j in range(0,len(multi_line)):
            line = ""
            if j == 0: # first line
                line = '%s"%s\n' % (base,multi_line[j])
            elif j == (len(multi_line)-1): # last line
                line = '%s%s"\n' % (blank_base,multi_line[j])
            else: # middle line
                line = '%s%s\n' % (blank_base,multi_line[j])
            out.write(line)


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
out_dir = str(sys.argv[2])

md = open(metadata,'r')

regex_for_product = '\s+/product="(.*)"'

# Iterate over the metadata file, one line per GBK to process
for line in md:
    line = line.strip('\n') 
    md_vals = line.split('\t')

    locus = ""
    if md_vals[1][-1:] == "_": # trim underscore if present
        locus = md_vals[1][:-1]
    else:
        locus = md_vals[1]

    gbk_in = "%s/%s/locus_mod.gbk" % (out_dir,locus)
    gbk = open(gbk_in,'r') # pull input GBK
    common_names = open(md_vals[2],'r') # pull name map file

    name_map = {}

    for name_pair in common_names: # build a map from the common names for conversion
        name_pair = name_pair.strip('\n')
        names = name_pair.split('\t')
        name_map[names[1]] = names[0]

    # Make a new file name for the outfile
    gbk_out = "%s/%s/common_name_mod.gbk" % (out_dir,locus)
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
            if "/translation=" in l or "/transl_table=" in l: # reached the end of the product, resolve it and leave
                printName(resolveName(multi_product,name_map),outfile)
                outfile.write(l)
                found_product = False # reinitialize all for next multi-line product
                multi_product = []
            else:
                l = l.strip('\n')
                multi_product.append(l.lstrip()) # strip leading white space

        else:
            outfile.write(l)
