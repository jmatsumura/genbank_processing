#!/bin/bash

# genbank_submission.sh
# Purpose: 
#	Runs a series of scripts to curate genbank files, 
#	 and produce relevant files for submission to NCBI
# Authors:
#	Shaun Adkins (sadkins@som.umaryland.edu)
#	James Matsumura (jmatsumura@som.umaryland.edu)

while [[ $# -ge 1 ]]
do
i="$1"
arg=$(echo $i | cut -f1 -d "=")
val=$(echo $i | cut -f2 -d "=")

case $arg in
    --metadata_list)
    metadata_list="$val"
    ;;
    --output_dir)
    output_dir="$val"
    ;;
esac
shift
done

# Read through metadata file, grab the locus, and mkdir with locus name
for locus in `awk '{print $2}' $metadata_list`; do
	new_locus=$(echo $locus | sed 's/.$//')
	mkdir -p 777 $output_dir/$new_locus
done

exit 0


cmd=/usr/bin/perl /path/to/script -arg1 arg1 -arg2 yes
echo "$cmd"
$cmd

#TODO: More Perl scripts in place

exit 0
