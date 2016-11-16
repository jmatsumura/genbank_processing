

# Script to arbitrarily join all contigs into one GBK file. Note that this
# does not use a linker sequence, simply concatenates all files. It is 
# recommended that a linker is incorporated inbetween the contigs.
#
# Run using the following command: (input file is a list of all contig file paths)
# python concat_contigs_to_gbk.py input_file output_file

import sys

contig_list_file = str(sys.argv[1])
out = str(sys.argv[2])

infile = open(contig_list_file,'r')
outfile = open(out,'a')

# Iterate over each line of each file and build the final file. Again, no linker
# being incorporated in this version. 
for single_file in infile:
    single_file = single_file.strip('\n')
    infile = open(single_file,'r')
    for line in infile:
        outfile.write(line)
