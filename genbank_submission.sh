#!/bin/bash

# genbank_submission.sh
# Purpose: 
#	Runs a series of scripts to curate genbank files, 
#	 and produce relevant files for submission to NCBI
# Authors:
#	Shaun Adkins (sadkins@som.umaryland.edu)
#	James Matsumura (jmatsumura@som.umaryland.edu)

# Python executable
PY_EXE=/usr/local/packages/Python-2.7.8/bin/python
# Directory of shell script (and other scripts)
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
# Path to used Expasy enzyme database file
EC_DAT=/local/projects/aengine/bin/enzyme.dat

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
    --enzyme_dat)
    enzyme_dat="$val"
    ;;
esac
shift
done

# Read through metadata file, grab the locus, and mkdir with locus name
for locus in `awk '{print $2}' $metadata_list`; do
	new_locus=$(echo $locus | sed 's/.$//')
	mkdir -p 777 $output_dir/$new_locus
done

cmd=$PY_EXE $DIR/locus_mod_gbk.py $metadata_list $output_dir
echo "$cmd"
$cmd

cmd=$PY_EXE $DIR/common_name_mod_gbk.py $metadata_list $output_dir
echo "$cmd"
$cmd

cmd=$PY_EXE $DIR/hypothetical_mod_gbk.py $metadata_list $output_dir
echo "$cmd"
$cmd

cmd=$PY_EXE $DIR/gene_symbol_mod_gbk.py $metadata_list $output_dir
echo "$cmd"
$cmd

cmd=$PY_EXE $DIR/ec_numbers_mod_gbk.py $metadata_list $output_dir $enzyme_dat
echo "$cmd"
$cmd

cmd=$PY_EXE $DIR/gbk2tbl.py $metadata_list $output_dir
echo "$cmd"
$cmd

exit 0
