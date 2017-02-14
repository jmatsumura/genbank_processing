

# Script to iterate through a GenBank file and pull out the following information
# and place into a TSV. There will be eight columns:
#
# 1. Locus
# 2. Locus tag
# 3. Coordinate start
# 4. Coordinate stop
# 5. Product
# 6. Gene symbol
# 7. EC number
# 8. Notes
#
# Note that it's possible not all of these values are present, although tabs
# will be incorporated as if all are present so looking at one column at a time
# will work.
#
# python parse_gbk_info.py /path/to/in.gbk /path/to/out.tsv

import sys, re

i = str(sys.argv[1]) # in GBK file
o = str(sys.argv[2]) # out TSV file

# Note that regex in the same order as the output
regex_for_locus = r'^LOCUS\s+([a-zA-Z0-9_\.]+)\s+' # 1
regex_for_locus_tag = r'^\s+/locus_tag="(.*)"' #
regex_for_coord_string = r'^\s+[a-zA-Z]+\s+((complement)?\(?(<?\d+)..(>?\d+)\)?)'
regex_for_product = r'^\s+/product="(.*)'
regex_for_multiline_product = r'^\s+([a-zA-Z0-9\s\'\-\(\)]+)"?'
regex_for_multiline_note = r'^\s+([a-zA-Z0-9\s\'\-\(\)_;:]+)"?'
regex_for_gene = r'^\s+/gene="(.*)"'
regex_for_ec = r'^\s+/EC_number="(.*)"'
regex_for_note = r'^\s+/note="(.*)'

out_list = [] # append one entry here for every locus tag entry

with open(i,'r') as gbk:

    l,lt,c1,c2,p,g,ec,n = ("" for i in range(8))
    multi_p,multi_n = (False for i in range(2)) # loops for multiline product/notes

    for line in gbk:

        line = line.rstrip('\n')

        # Use the gene feature marker as a breaking point for when to complete
        # an entry.
        if '   gene   ' in line: 

            # If a gene feature is reached and we've already gathered data on the
            # previous entry then add this entry to the list. 
            if lt != "":
                out_list.append(("\t".join([l,lt,c1,c2,p,g,ec,n])))

                # Reinitialize all except for Locus as that one is overarching
                # across all entries in the current section.
                lt,c1,c2,p,g,ec,n = ("" for i in range(7)) 

        elif line.startswith('LOCUS'): # grab locus
            
            # Handle the last gene with a locus entry.
            if lt != "":
                out_list.append(("\t".join([l,lt,c1,c2,p,g,ec,n])))
                lt,c1,c2,p,g,ec,n = ("" for i in range(7))

            l = re.search(regex_for_locus,line).group(1)

        elif '/locus_tag=' in line: # grab locus tag
            lt = re.search(regex_for_locus_tag,line).group(1)

        # Grab the coords
        elif '   CDS   ' in line or '   tRNA   ' in line or '   rRNA   ' in line:
            results = re.search(regex_for_coord_string,line)
            if results.group(2): # check if complement is present, swap coords if so
                c1 = results.group(4)
                c1 = c1.replace(">","<")
                c2 = results.group(3)
                c2 = c2.replace("<",">")
            else: # normal 5'->3'
                c1 = results.group(3)
                c2 = results.group(4)
        
        elif multi_p == True: # loop for finding subsequent product lines
            add_p = re.search(regex_for_multiline_product,line).group(1)
            p = "{0} {1}".format(p,add_p)

            if line.endswith('"'): # found the end of the product name
                multi_p = False 

        elif '/product=' in line: # grab the annotated product
            p = re.search(regex_for_product,line).group(1)
            
            # Occasionally there will be multi-line products. Handle here. 
            if p.endswith('"'):
                p = p[:-1]
            else:
                multi_p = True
            
        elif '/gene=' in line: # grab the gene symbol
            g = re.search(regex_for_gene,line).group(1)

        elif '/EC_number=' in line: # grab the gene symbol
            ec = re.search(regex_for_ec,line).group(1)

        elif multi_n == True: # loop for finding subsequent product lines
            add_n = re.search(regex_for_multiline_note,line).group(1)
            n = "{0} {1}".format(n,add_n)

            if line.endswith('"'): # found the end of the product name
                multi_n = False 

        elif '/note=' in line: # grab the gene symbol
            n = re.search(regex_for_note,line).group(1)

            if n.endswith('"'):
                n = n[:-1]
            else:
                multi_n = True

with open(o,'w') as out:
    for entry in out_list:
        out.write("{0}\n".format(entry))
