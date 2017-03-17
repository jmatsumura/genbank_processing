#!/usr/bin/python
#
# Script to parse out the location of contaminants from a file
# with a format like:
#
#------------------
#abcdefghi...
#
#Trim:
#Sequence name, length, span(s), apparent source
#N103P20H1_10_1  34983   1..54,34959..34983      adaptor:NGB00735.1
#N103P20H1_23_1  31351   1..63   adaptor:NGB00735.1
#
# abcdefghi...
#------------------
#
# HOW TO RUN: python extract_contaminants.py -contaminant_file contaminants.txt -out_dir /path/to/out_dir_prefix -tbl total
#
# This will produce a number of files...
# *_fsa = IDs to be used for trimming contaminants from FSA files
# *_tbl = IDs to be used for trimming contaminants from TBL files
# contigs_to_remove.txt = IDs to be used as input for remove_contigs.py
#
# Author: James Matsumura

import sys, argparse

def main():

	parser = argparse.ArgumentParser(description='Script to prep files for trim_contamination.pl.')
	parser.add_argument('-contaminant_file', type=str, required=True, help='Path to a contamination file from Path to map.tsv output from format_for_assembly.py or final_verdict.py.')
	parser.add_argument('-out_dir', type=str, required=True, help='Prefix to directory where all this output will be written.')
	parser.add_argument('-tbl', type=str, required=True, help='Path to a TBL file that contains ALL entries that are going to be considered for contaminant purging.')
	args = parser.parse_args()
	
	tbl_ids = set()

	with open(args.tbl,'r') as tbl:
		for line in tbl:
			line = line.rstrip()
			tbl_ids.add(line.split(' ')[1])

	with open(args.contaminant_file,'r') as contamination:

		contaminant_section,exclusion_section = (False for i in range(2))
		contaminant_list,exclusion_list = ([] for i in range(2))

		for line in contamination:

			# Find those we want to trim.
			if line.startswith('Trim:'):
				contaminant_section = True
				contaminant_list = []

			elif contaminant_section == True and line.startswith('Sequence name'):
				pass

			# If we've found all the contaminants, leave
			elif contaminant_section == True and line.strip() =='':
				genome = contaminant_list[0].split('\t')[0].split('_')[0]

				first_round,second_round = ([] for i in range(2))

				for contam in contaminant_list:

					# One contig has multiple areas to trim
					if ',' in contam.split('\t')[2]:
						print(contam.rstrip())
						coords = contam.split('\t')[2]
						str_1 = contam.split('\t')
						str_2 = contam.split('\t')
						str_1[2] = coords.split(',')[0] # reassign coordinates to individual entries
						# The second coords need a little bit of work as they are now
						# shifted due to trimming of the earlier coords. 
						nums = coords.split(',')[0].split('..')
						shift_val = (int(nums[1]) - int(nums[0]) + 1)
						new_nums = coords.split(',')[1].split('..')
						new_nums[0] = (int(new_nums[0]) - shift_val)
						new_nums[1] = (int(new_nums[1]) - shift_val)
						new_coords = "{0}..{1}".format(new_nums[0],new_nums[1])
						str_2[2] = new_coords
						first_round.append(('\t').join(str_1))
						second_round.append(('\t').join(str_2))
					else:
						first_round.append(contam)

				with open("{0}/{1}_fsa_ids".format(args.out_dir,genome),'w') as out:
					for x in first_round:
						out.write(x)
				with open("{0}/{1}_tbl_ids".format(args.out_dir,genome),'w') as out:
					for x in first_round:
						if x.split('\t')[0] in tbl_ids:
							out.write(x)

				if len(second_round) > 0:
					with open("{0}/{1}_fsa_ids_2".format(args.out_dir,genome),'w') as out:
						for x in second_round:
							out.write(x)
					with open("{0}/{1}_tbl_ids_2".format(args.out_dir,genome),'w') as out:
						for x in second_round:
							if x.split('\t')[0] in tbl_ids:
								out.write(x)
				
				contaminant_section = False

			elif contaminant_section == True:
				contaminant_list.append(line)

			# Find those we want to exclude.
			if line.startswith('Exclude:'):
				exclusion_section = True
				exclusion_list = []

			elif exclusion_section == True and line.startswith('Sequence name'):
				pass

			# If we've found all the contaminants, leave
			elif exclusion_section == True and line.strip() =='':
				
				with open("{0}/contigs_to_remove.txt".format(args.out_dir),'a') as out:
					for x in exclusion_list:
						out.write(x)

				exclusion_section = False

			elif exclusion_section == True:
				exclusion_list.append(line)

if __name__ == '__main__':
	main()
