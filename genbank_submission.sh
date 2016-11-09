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
    --input_list)
    input_list="$val"
    ;;
    --output_dir)
    output_dir="$val"
    ;;
esac
shift
done

# If output directory does not exist, then create it
if [ ! -d $output_dir ]; then 
	mkdir -p 777 $output_dir
fi

cmd=/usr/bin/perl /path/to/script -arg1 arg1 -arg2 yes
echo "$cmd"
$cmd

#TODO: More Perl scripts in place

exit 0
