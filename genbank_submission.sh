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
	mkdir -p -m 777 $output_dir/$new_locus 
done

# modify locus tag names in genbank file
cmd=$PY_EXE $DIR/locus_mod_gbk.py $metadata_list $output_dir
echo "$cmd"
$cmd || { echo 'locus_mod_gbk.py failed!' ; exit 1; }

# curate common names in genbank file
cmd=$PY_EXE $DIR/common_name_mod_gbk.py $metadata_list $output_dir
echo "$cmd"
$cmd || { echo 'common_name_mod_gbk.py failed!' ; exit 1; }

# remove hypotheticals in genbank file
cmd=$PY_EXE $DIR/hypothetical_mod_gbk.py $metadata_list $output_dir
echo "$cmd"
$cmd || { echo 'hypothetical_mod_gbk.py failed!' ; exit 1; }

# modify gene symbols in genbank file
cmd=$PY_EXE $DIR/gene_symbol_mod_gbk.py $metadata_list $output_dir
echo "$cmd"
$cmd || { echo 'gene_symbol_mod_gbk.py failed!' ; exit 1; }

# modify EC numbers in genbank file
cmd=$PY_EXE $DIR/ec_numbers_mod_gbk.py $metadata_list $output_dir $EC_DAT
echo "$cmd"
$cmd || { echo 'ec_numbers_mod_gbk.py failed!' ; exit 1; }

# gbk2tbl - round 1
cmd=$PY_EXE $DIR/gbk2tbl.py $metadata_list $output_dir
echo "$cmd"
$cmd || { echo 'gbk2tbl.py failed!' ; exit 1; }

# create sbt and cmt files
cmd=perl $DIR/prepare_sbt_cmt.pl --input_file=$metadata_list --output_dir=$output_dir
echo "$cmd"
$cmd || { echo 'prepare_sbt_cmt.pl failed!' ; exit 1; }

# run tbl2asn - round 1
# delete overlapping genes
# gbk2tbl - round 2
# fix gene symbols in tbl file
for locus in `awk '{print $2}' $metadata_list`; do
	new_locus=$(echo $locus | sed 's/.$//')

	# create agp file
	cmd=perl $DIR/format_velvet_assembled_AGP.pl --tbl_file=$output_dir/$new_locus/${new_locus}_gs_corrected.tbl --fsa_file=$output_dir/$new_locus/${new_locus}.fsa --split_param=10 --min_contig_len=200 --output_dir=$output_dir/$new_locus
	echo "$cmd"
	$cmd || { echo 'prepare_sbt_cmt.pl failed!' ; exit 1; }

	# rm 2nd tbl file
	/bin/rm $output_dir/$new_locus/${new_locus}_gs_corrected.tbl

	# rename corrected tbl file
	/bin/mv $output_dir/$new_locus/${new_locus}_gs_corrected_corrected.tbl $output_dir/$new_locus/${new_locus}.tbl

	# rename corrected fsa file
	/bin/mv $output_dir/$new_locus/${new_locus}_corrected.fsa $output_dir/$new_locus/${new_locus}.fsa
	
	# rename corrected agp file
	/bin/mv $output_dir/$new_locus/${new_locus}_gs_corrected.agp $output_dir/$new_locus/${new_locus}.agp
done
# run tbl2asn - round 2
exit 0
