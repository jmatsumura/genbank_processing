#!/usr/bin/python

# Script to convert from GBK to TBL. 
#
# Credit to ~wanyuac (https://github.com/wanyuac/BINF_toolkit) for help 
# understanding BioPython's parser.

import sys, re
from Bio import SeqIO

metadata = str(sys.argv[1])
out_dir = str(sys.argv[2])

md = open(metadata,'r')

# Decide which qualifiers to output under each feature
gene_output_qualifiers = ["gene","locus_tag"]
cds_or_rna_output_qualifiers = ["product","prot_desc","note","go_component"
        ,"go_process","go_function","db_xref","pseudo","EC_number","protein_id"]

gquals = set(gene_output_qualifiers)
crquals = set(cds_or_rna_output_qualifiers)

cds_and_rna = ["rRNA","tRNA","CDS"] # these fields share the same quals

# Iterate over the metadata file, one line per GBK to process
for line in md:
    line = line.strip('\n') 
    md_vals = line.split('\t')

    gbk_in = "%s/%s/ec_numbers_mod.gbk" % (out_dir,md_vals[1][:-1])
    gbk = open(gbk_in,'rU') # pull input GBK

    # Make new TBL and FSA files which are needed for tbl2asn
    tbl_out = "%s/%s/out.tbl" % (out_dir,md_vals[1][:-1])
    tbl_outfile = open(tbl_out,'w')
    fsa_out = "%s/%s/out.fsa" % (out_dir,md_vals[1][:-1])
    fsa_outfile = open(fsa_out,'w')

    records = SeqIO.parse(gbk, 'genbank') # get all GBK entries

    # Iterate over each entry in the GBK file.
    # Note that a downstream script will filter out those entries <200bp
    for rec in records:

        SeqIO.write(rec, fsa_outfile, "fasta") # write the FASTA file

        # Build the .tbl file
        header = ">Feature %s\n" % (rec.name)
        tbl_outfile.write(header)

        # Build a protein_id as you go, only add if gene present.
        protein_id = ""

        # Iterate over all features in a given entry. 
        for f in rec.features:
            # First print coordinates
            if f.type == 'source':
                pass
            else: 
                coords = ""
                if f.strand == 1:
                    coords = "%d\t%d\t%s\n" % (f.location.nofuzzy_start + 1, f.location.nofuzzy_end, f.type)
                else:
                    coords = "%d\t%d\t%s\n" % (f.location.nofuzzy_end, f.location.nofuzzy_start + 1, f.type)
                tbl_outfile.write(coords)

            # Append protein_id here if it is present
            if f.type in cds_and_rna and protein_id != "":
                f.qualifiers['protein_id'] = protein_id # note this assigns one char per index due to BioPython

            # Now print any of the fields within the gene/CDS/RNA entries.
            for key,values in f.qualifiers.iteritems():

                if f.type == 'gene':
                    # First, if a 'gene' is found must introduce a protein_id
                    protein_id = "gnl|IGS|%s" % (f.qualifiers['locus_tag'][0])

                    # Now add the relevant fields
                    if key in gquals:
                        for vals in values:
                            qualifier = "\t\t\t%s\t%s\n" % (key,vals)
                            tbl_outfile.write(qualifier)

                if f.type in cds_and_rna and key in crquals:
                    if key != "protein_id":
                        for vals in values:
                            qualifier = "\t\t\t%s\t%s\n" % (key,vals)
                            tbl_outfile.write(qualifier)
                    elif key == "protein_id": # special case, need to build an aggregate string from list
                        val = ''.join(f.qualifiers['protein_id'])
                        qualifier = "\t\t\t%s\t%s\n" % (key,val)
                        tbl_outfile.write(qualifier)