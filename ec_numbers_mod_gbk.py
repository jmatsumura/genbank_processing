

# Script to modify EC values in a GBK file. 
# Note that unlike the other scripts this requires a second argument
# for a enzyme.dat file. 

import sys, re
from collections import defaultdict
from shared_fxns import printName

metadata = str(sys.argv[1])
out_dir = str(sys.argv[2])
enzyme_file = str(sys.argv[3])

md = open(metadata,'r')
ed = open(enzyme_file,'r')

regex_for_ec_num = r'\s+/EC_number="(.*)"'
regex_for_ec_id = r'^ID\s+(.*)$'
regex_for_ec_transferred = r'^DE\s+Transferred entry:\s(.*)$'
regex_for_product_name = r'^DE\s+(.*)$'
regex_for_gbk_product = '/product="(.*)'

transferred_ec,old_ec_product,new_ec_product = ({} for i in range(3)) # dict for ECs that have been updated
# maintain a dictionary for output and tell users which were deleted, transferred, and multimapped
out_map = defaultdict(dict)

# First, extract all the relevant information from the enzyme.dat file
id = "" # single EC number ID
product_list = [] # use a list to handle multiline products

for line in ed:
    if line.startswith('ID'): # grab most recent ID
        id = re.search(regex_for_ec_id,line).group(1)
    elif 'Transferred entry' in line: # transferred entries are the new valid EC numbers
        all_ecs = re.search(regex_for_ec_transferred,line).group(1)
        indiv_ecs = all_ecs.split(' ')
        if len(indiv_ecs) > 1: # hit multiple, need manual curation let the discrep file catch
            out_map['multi-mapped'][id] = ' '.join(indiv_ecs)
        else:
            ec = all_ecs
            if not ec[-1:].isdigit():
                ec = ec[:-1] # trim commas and periods
            transferred_ec[id] = ec # assign to new transferred entry EC
            out_map['single-mapped'][id] = ec
    elif line.startswith('DE'): # found a product name if starting with DE and not transferred entry
        # Know that the previous ID extracted is what this product name is tied to.
        product = re.search(regex_for_product_name,line).group(1) 

        if product.endswith('.'): # got the full product    
            product_list.append(product[:-1])
            new_ec_product[id] = ''.join(product_list)
            product_list = [] # reinitialize

        else: # not the full name, just grab the first part
            product_list.append(product)


# Iterate over the metadata file, one line per GBK to process
for line in md:
    line = line.strip('\n') 
    md_vals = line.split('\t')

    locus = ""
    if md_vals[1][-1:] == "_": # trim underscore if present
        locus = md_vals[1][:-1]
    else:
        locus = md_vals[1]

    gbk_in = "%s/%s/gene_symbol_mod.gbk" % (out_dir,locus)
    gbk = open(gbk_in,'r') # pull input GBK

    # Now grab all the EC numbers that are to be deleted noted by the metadata file. 
    ec_numbers = open(md_vals[18],'r') # pull name map file
    ecs_to_delete = []
    for l in ec_numbers:
        l = l.strip('\n')
        ecs_to_delete.append(l)
    delete_us = set()
    if len(ecs_to_delete) > 0:
        delete_us = set(ecs_to_delete)

    # Make a new file name for the outfile
    gbk_out = "%s/%s/ec_numbers_mod.gbk" % (out_dir,locus)
    outfile = open(gbk_out,'w')

    ec_list = [] # gather all EC numbers and note which have changed
    product_list = [] # capture multiline products for outfile
    rename_product,within_product = (False for i in range(2))
    mod_ec,curr_ec = ("" for i in range(2))

    for l in gbk:

        if "/EC_number=" in l: # found the EC number
            curr_ec = re.search(regex_for_ec_num,l).group(1)
            ec_list.append(curr_ec)

            if curr_ec in transferred_ec: # if updated EC and single mapped, replace
                base = " " * 21 # no replacing digits easily, just build new string
                ec_str = '/EC_number="%s"\n' % (transferred_ec[curr_ec])
                l = "%s%s" % (base,ec_str)
                rename_product = True

            mod_ec = re.search(regex_for_ec_num,l).group(1) # potentially modified EC

            if mod_ec not in delete_us: # only print those that aren't to be deleted
                outfile.write(l)
        
        else:
            if rename_product == True:
                if "/product=" in l: # skip over old product, write a new one
                    printName(new_ec_product[mod_ec],outfile)
                    tmp_line = l.strip()

                    if not tmp_line.endswith('"'): # multiline product
                        within_product = True
                        product_list.append(re.search(regex_for_gbk_product,l).group(1))

                    else: # single line product
                        product = re.search(regex_for_gbk_product,l).group(1)
                        product = product[:-1]
                        old_ec_product[curr_ec] = product # put the old EC number in 

                elif within_product == True:
                    tmp_line = l.strip()

                    if tmp_line.endswith('"'): # found the end of the product
                        tmp_line = tmp_line[:-1]
                        product_list.append(tmp_line)
                        old_ec_product[curr_ec] = ''.join(product_list)

                        within_product = False
                        product_list = [] # reinitialize
                    
                    else:
                        product_list.append(tmp_line)

                elif "/translation=" in l or "/transl_table=" in l: # reached the end of the product, resolve it and leave
                    outfile.write(l)
                    rename_product = False

            else: # if don't need to rename, just keep the current structure
                outfile.write(l)

    ec_hist = "%s/%s/ec_history.tsv" % (out_dir,locus)
    histfile = open(ec_hist,'w')

    for map_type,id_key in out_map.items():
        for original_id,new_id in out_map[map_type].items():
            if original_id in ec_list and map_type == 'single-mapped':
                histfile.write("%s\t%s\t%s\t%s\t%s\n" % (map_type,original_id,old_ec_product[original_id],new_id,new_ec_product[new_id]))
            elif original_id in ec_list: # write out multimaps
                histfile.write("%s\t%s\trequires curation for which of these to change to:\t%s\n" % (map_type,original_id,out_map[map_type][original_id]))
