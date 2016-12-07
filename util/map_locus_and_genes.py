

# Script to create a map between the final .gbk file of the pipeline and the 
# original GBK file. Note that this script is quite greedy and probably won't 
# work in many cases, but suffices for the initial 270 cases. In its current
# state, it just aims to map loci/contigs by coordinates. The final output will be
# of the following format:
#
# modified_contig modified_locus original_contig original_locus coords
# Note that coords will be shared between the two
# 
# If you ran the GenBank pipeline, the latest GBK file that you should use to
# map back from is titled "delete_overlap_mod.gbk"
#
# Run using the following command: (inputs are the final *.gbk file and the initial *.gbk file.)
# python map_locus_and_genes.py end.gbk start.gbk /path/to/outfile

import sys, re

end = str(sys.argv[1]) # resultant/modified GBK file
beg = str(sys.argv[2]) # initial GBK file
output = str(sys.argv[3])

e = open(end,'r')
b = open(beg,'r')
o = open(output,'w')

# I know it's ugly, but the dict key will be coordinates with the value
# being locus+locus_tag. For instance: dict[1..5] = "N116P14_123+TESTING_0001"
emap,bmap = ({} for i in range(2))

regex_for_coord_string = r'^\s+CDS\s+((complement)?\(?<?\d+..>?\d+\)?)'
regex_for_locus_tag = r'^\s+/locus_tag="(.*)"'
regex_for_locus_and_bp = r'^LOCUS\s+([a-zA-Z0-9_\.]+)\s+([0-9]+\sbp)'
regex_for_translation = r'^\s+/translation="([A-Z]+)'

def isolate_locus_and_coords(gbk):
    locus,coords,ltag,translation = ("" for i in range(4))
    map = {}
    within_locus,cds_found,ltag_found,translation_found = (False for i in range(4))

    for line in gbk:
        line = line.strip('\n')

        if line.startswith('LOCUS') and within_locus == False:
            within_locus = True
            # Grab locus
            locus = re.search(regex_for_locus_and_bp,line).group(1)
            length = re.search(regex_for_locus_and_bp,line).group(2)

        elif within_locus == True:
            if line.startswith('ORIGIN'):
                within_locus = False # leave at this point, but still append to map
            if cds_found == True and ltag_found == True and translation_found == True:
                unique_id = "%s|%s|%s" % (length,coords,translation) # unique id of contig length + coords
                if unique_id not in map:
                    map[unique_id] = "%s+%s" % (locus,ltag)
                else:
                    print "ERROR, duplicates found across loci using key (contig_length|coords|translation): %s" % (unique_id)
                cds_found,ltag_found,translation_found = (False for i in range(3))
            elif "/locus_tag" in line:
                # Grab locus tag
                ltag = re.search(regex_for_locus_tag,line).group(1)
                ltag_found = True
            elif re.search(regex_for_coord_string,line):
                # Get CDS coordinate
                coords = re.search(regex_for_coord_string,line).group(1)
                cds_found = True
            elif re.search(regex_for_translation,line):
                # Get the translation sequence to make an even more unique key
                translation = re.search(regex_for_translation,line).group(1)
                translation_found = True
                
    return map

# Find the modified coordinates
emap = isolate_locus_and_coords(e)

# Find the original coordinates
bmap = isolate_locus_and_coords(b)

# Write out the mapped contigs/loci based on coordinates. Note that we want to 
# iterate over the final file that potentially has deleted genes as to only map 
# those genes that actually make it to the final GBK file. 
for k1,v1 in emap.iteritems():
    k2 = k1
    if k2 not in bmap:
        print "Uh oh, contig_length|coords|translation of %s not found in initial file." % (k2)
        v2 = "N/A+N/A"
    else:
        v2 = bmap[k2]
    e1 = v1.split('+')
    e2 = v2.split('+')
    l1 = e1[0] # grab locus
    lt1 = e1[1] # grab locus tag
    l2 = e2[0]
    lt2 = e2[1]
    # output locus(new),locus_tag(new),locus(old),locus_tag(old),coordinates
    outstr = "%s\t%s\t%s\t%s\t%s\n" % (l1,lt1,l2,lt2,k1)
    o.write(outstr)
