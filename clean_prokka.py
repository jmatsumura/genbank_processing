#!/usr/bin/python

# Script to clean various erroneous aspects of Prokka output.
# Namely, fixing BP number in LOCUS line and moving all the qualifiers
# within the CDS region that belong in a separate gene field.

import sys, re
from Bio import SeqIO

metadata = str(sys.argv[1])
out_dir = str(sys.argv[2])

md = open(metadata,'r')

regex_for_cds = r'\s+CDS\s+\d+..\d+'
# Note that capturing coords must accept complement syntax too
regex_for_coords = r'\s+[a-zA-Z]+\s+(.*\d+..\d+\)*)' 
regex_for_generic = r'\s+([a-zA-Z]+)\s+\d+..\d+'

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

    # If Prokka, need to use the rRNA/tRNA/CDS region info to populate a gene field
    extract_gene = False

    # Establish empty lists for each gene/CDS field
    gene_entry,cdsrna_entry = ([] for i in range(2))

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

        else: # Prokka file! Time to rearrange gene data embedded in CDS/RNA fields.
            if extract_gene == True:
                if line.startswith('ORIGIN'): # reached sequence region
                    for x in gene_entry:
                        outfile.write(x)
                    for x in cdsrna_entry:
                        outfile.write(x)
                    outfile.write(line) # write out ORIGIN line
                    gene_entry,cdsrna_entry = ([] for i in range(2)) # reset 
                    extract_gene = False

                # Try find some field that needs gene info extracted
                elif re.search(regex_for_generic,line):
                    if re.search(regex_for_generic,line).group(1) != 'source':
                        # Found a new CDS entry, add the gene / CDS lines
                        for x in gene_entry:
                            outfile.write(x)
                        for x in cdsrna_entry:
                            outfile.write(x)
                        gene_entry,cdsrna_entry = ([] for i in range(2)) # reset 
                        cdsrna_entry.append(line)
                        coords = re.search(regex_for_coords,line).group(1)
                        gene_str = ' '*5 + 'gene' + ' '*12 + coords + '\n'
                        gene_entry.append(gene_str)
                        extract_gene = True
                    else:
                        outfile.write(line)
                        
                else:
                    if '/gene="' in line: 
                        gene_entry.append(line)
                        cdsrna_entry.append(line)
                    elif '/locus_tag="' in line:
                        gene_entry.append(line)
                    else:
                        cdsrna_entry.append(line)
            else:
                if re.search(regex_for_generic,line):
                    if re.search(regex_for_generic,line).group(1) != 'source':
                        # Grab 'headers' for both gene and other field
                        cdsrna_entry.append(line)
                        coords = re.search(regex_for_coords,line).group(1)
                        gene_str = ' '*5 + 'gene' + ' '*12 + coords + '\n'
                        gene_entry.append(gene_str)
                        extract_gene = True
                    else:
                        outfile.write(line)
                else:
                    outfile.write(line)
