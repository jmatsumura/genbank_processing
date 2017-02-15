

# Script to delete completely overlapping CDS regions. Explicitly, this means the
# following will be deleted:
# 1) A CDS entirely within a CDS
# 2) A CDS entirely within a xRNA
# The following will not be deleted:
# 1) A xRNA entirely within a CDS

import sys, re

metadata = str(sys.argv[1])
out_dir = str(sys.argv[2])

md = open(metadata,'r')

regex_for_gene = r'^\s+gene\s+[complement]*\(*\d+..\d+\)*'
regex_for_rna = r'^\s+[a-z]*RNA\s+[complement]*\(*\d+..\d+\)*'
regex_for_locus = r'^\s+/locus_tag="(.*)"'

# Iterate over the metadata file, one line per GBK to process
for line in md:
    line = line.strip('\n') 
    md_vals = line.split('\t')

    locus = ""
    if md_vals[1][-1:] == "_": # trim underscore if present
        locus = md_vals[1][:-1]
    else:
        locus = md_vals[1]

    # Pull the discrep file in generated by tbl2asn.
    disc_in = "%s/%s/%s_discrep.txt" % (out_dir,locus,locus)
    disc = open(disc_in,'r')

    # Pull the GBK file in to potentially modify
    gbk_in = "%s/%s/ec_numbers_mod.gbk" % (out_dir,locus)
    gbk = open(gbk_in,'r')

    # Write out a gbk file with overlaps deleted (if any are present)
    out = "%s/%s/delete_overlap_mod.gbk" % (out_dir,locus)
    outfile = open(out,'w')

    # Write out which genes are being deleted
    delete_out = "%s/%s/deleted_ids.txt" % (out_dir,locus)
    delete_out = open(delete_out,'w')

    overlap = False
    found_overlapping_entries = False
    delete_us = []

    # Iterate over the discrep file to find any relevant CDS/RNA within CDS/RNA.
    for line in disc:
        if line.startswith('FIND_OVERLAPPED_GENES:') and overlap == False:
            overlap = True
        # If reached detailed report and haven't found the above string, no overlap present
        elif line.startswith('Detailed Report') and overlap == False: 
            break

        if line.startswith('DiscRep_ALL:FIND_OVERLAPPED_GENES::'):
            found_overlapping_entries = True
        elif found_overlapping_entries == True:
            if line.startswith('DiscRep'): # reached next discrep set, done finding overlaps
                found_overlapping_entries = False
                break # got all that we need from discrep file
            else: # identify the loci that are entirely contained within another
                if '\t' in line:
                    line = line.strip('\n')
                    elements = line.split('\t')
                    delete_us.append(elements[-1]) # grab the locus tag to delete
                    delete_out.write("%s\t%s\n" % (elements[-1],elements[-2])

    within_gene = False
    rna_within = False
    gene_region = []

    # Iterate over the GBK file and skip any gene entries that are found entirely
    # within another. 
    for line in gbk:
        if overlap == True:

            # If we haven't seen a gene region...
            if re.search(regex_for_gene,line) and within_gene == False:
                gene_region.append(line)
                within_gene = True

            # These two cases warrant a check for whether or not we are
            # to remove this gene region due to overlap.  
            elif (re.search(regex_for_gene,line) and within_gene == True) or line.startswith('ORIGIN'):
                
                locus = ""
                for ele in gene_region:
                    if "/locus_tag" in ele:
                        locus = re.search(regex_for_locus,ele).group(1)
                    elif re.search(regex_for_rna,ele):
                        rna_within = True

                if rna_within == True:
                    for ele in gene_region:
                        outfile.write(ele)
                elif locus in delete_us: # if it is to be deleted, just skip
                    pass
                else:
                    for ele in gene_region:
                        outfile.write(ele)

                # Reinitialize
                within_gene = False
                rna_within = False
                gene_region = []

                if line.startswith('ORIGIN'):
                    outfile.write(line)

                # Need to catch the wrapper "if" check within this elif
                elif re.search(regex_for_gene,line) and within_gene == False:
                    gene_region.append(line)
                    within_gene = True

            # Anything in-between is within the gene region and should be added
            elif within_gene == True:
                gene_region.append(line)

            # Simply write out any of the lines outside these regions
            else:
                outfile.write(line)

        else:
            outfile.write(line) # simply copying over, straying away from extra dependencies (cp cmd)
