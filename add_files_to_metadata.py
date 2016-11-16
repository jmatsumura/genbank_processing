

# Script to add the GBK/GBFF files to the first column of a 
# metadata file used in the GenBank submission pipeline. Note
# that this requires an input of paths to files that can be 
# identified by the first column of the metadata file which
# will correspond to the end of the path (either gbk/gbff ext). 
# For instance, 101ABC would map to /path/to/101ABC.gbk or 
# /path/to/101ABC.gbff.
#
# In order to generate this file, you can use the bash command
# "ls -d -1 $PWD/*.*" in the directory of your gbk/gbff files. 
# Use output redirection (">" or ">>") to build a final file For
# multiple paths, making sure that the file name can be tied 
# exactly to column 1 of the metadata file. 
#
# Run this script like so:
# python add_files_to_metadata.py input_metadata.tsv file_path_list.txt /path/to/outputfile.txt

import sys, re

metadata = str(sys.argv[1])
file_list = str(sys.argv[2])
output = str(sys.argv[3])

md = open(metadata,'r')
fl = open(file_list,'r')
out = open(output,'w')

regex_for_id = r'.*/(.*)\.gb.*'

map = {}

for file_path in fl:
    file_path = file_path.strip('\n')
    id = re.search(regex_for_id,file_path).group(1)
    map[id] = file_path

# Have the map, now build the new output file.
for row in md:
    row = row.strip('\n')
    elements = row.split('\t')
    if elements[0] not in map:
        print "ERROR. First column value not found in file path list."
    else:
        elements[0] = map[elements[0]]
    out.write('\t'.join(elements))
