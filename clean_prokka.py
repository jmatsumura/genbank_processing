#!/usr/bin/python

# Script to clean various erroneous aspects of Prokka output.
# Namely, fixing BP number in LOCUS line 

import sys, re
from Bio import SeqIO

metadata = str(sys.argv[1])
out_dir = str(sys.argv[2])

md = open(metadata,'r')

# Decide which qualifiers to output under each feature
gene_qualifiers = ["gene","locus_tag"]
gquals = set(gene_qualifiers)

# Iterate over the metadata file, one line per GBK to process
for line in md:
    line = line.strip('\n') 
    md_vals = line.split('\t')

    gbk = open(md_vals[0],'rU') # pull input GBK

    # Intermediate, kind of clean file (formatting length). This
    # also serves as an early warning for any errors determined by BioPython
    fixed_length = "%s/%s/intermediate.gbk" % (out_dir,md_vals[1][:-1])
    firstfile = open(fixed_length,'w')

    records = SeqIO.parse(gbk, 'genbank') # get all GBK entries

    # Iterate over each entry in the GBK file.
    for rec in records:
        SeqIO.write(rec, firstfile, "genbank")

    firstfile.close() # close and re-open as the GBK file
    gbk = open(fixed_length,'r')

    # Make cleaned GBK file which has been reformatted.
    cleaned_out = "%s/%s/cleaned.gbk" % (out_dir,md_vals[1][:-1])
    outfile = open(cleaned_out,'w')

    # Not Prokka until proven otherwise
    prokka = False

    # If Prokka, need to use the CDS region info to populate a gene field
    foundCDS = False

    # Now process without BioPython and relocate gene qualifiers from CDS
    # section into a separate gene section if this input is from Prokka.
    for line in gbk:
        if prokka == False:
            if line.startswith('COMMENT'):
                if 'using prokka' in line:
                    outfile.write(line)
                    prokka = True
            else:
                outfile.write(line)

        else: # Prokka file! Time to rearrange gene data embedded in CDS.
            if foundCDS == True:
                pass
            else:
                if ' CDS ' in line:
                    foundCDS = True
                else:
                    outfile.write(line)
