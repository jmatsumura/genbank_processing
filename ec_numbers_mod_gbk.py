#!/usr/bin/python

# Script to modify EC values in a GBK file. 
# Note that unlike the other scripts this requires a second argument
# for a enzyme.dat file. 

import sys, re

metadata = str(sys.argv[1])
out_dir = str(sys.argv[2])
enzyme_file = str(sys.argv[3])

md = open(metadata,'r')
ed = open(enzyme_file,'r')

regex_for_ec_num = r'\s+/EC_number="(.*)"'
regex_for_ec_id = r'^ID\s+(.*)$'
regex_for_ec_transferred = r'^DE\s+Transferred entry:\s(.*)$'

transferred_ec = {} # dict for ECs that have been updated

# First, extract all the relevant information from the enzyme.dat file
id = "" # single EC number ID

for line in ed:
    if line.startswith('ID'): # grab most recent ID
        id = re.search(regex_for_ec_id,line).group(1)
    elif 'Transferred entry' in line: # transferred entries all map to most recent ID
        all_ecs = re.search(regex_for_ec_transferred,line).group(1)
        indiv_ecs = all_ecs.split(' ')
        for ec in indiv_ecs:
            if ec[0].isdigit(): # found an EC number, not 'and' or some nonsense
                if not ec[-1:].isdigit():
                    ec = ec[:-1] # trim commas and periods
                transferred_ec[ec] = id # set all transferred to new EC

# Iterate over the metadata file, one line per GBK to process
for line in md:
    line = line.strip('\n') 
    md_vals = line.split('\t')

    gbk_in = "%s/%s/gene_symbol_mod.gbk" % (out_dir,md_vals[1][:-1])
    gbk = open(gbk_in,'r') # pull input GBK

    # Now grab all the EC numbers that are to be deleted noted by the metadata file. 
    ec_numbers = open(md_vals[18],'r') # pull name map file
    ecs_to_delete = []
    for l in ec_numbers:
        l = l.strip('\n')
        ecs_to_delete.append(l)
    delete_us = set(ecs_to_delete)

    # Make a new file name for the outfile
    gbk_out = "%s/%s/ec_numbers_mod.gbk" % (out_dir,md_vals[1][:-1])
    outfile = open(gbk_out,'w')

    for l in gbk:

        if "/EC_number=" in l: # found the EC number
            curr_ec = re.search(regex_for_ec_num,l).group(1)
            if curr_ec in transferred_ec: # if updated EC, replace
                base = " " * 21 # no replacing digits easy, just build new string
                ec_str = '/EC_number="%s"\n' % (transferred_ec[curr_ec])
                l = "%s%s" % (base,ec_str)
            mod_ec = re.search(regex_for_ec_num,l).group(1) # potentially modified EC
            if mod_ec not in delete_us: # only print those that aren't to be deleted
                outfile.write(l)
        else:
            outfile.write(l)
