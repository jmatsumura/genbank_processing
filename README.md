# genbank_processing
Scripts for reformatting GenBank files given a metadata file. Will
do the following:

- Remove overlapping genes (not *RNAs)
- Correct gene symbol syntax for prokaryotic genes
- Remove contigs that are <200BP and reorder the others accordingly
- Replace locus tags
- Update EC numbers
- Remove hypothetical gene symbols

## Order of the scripts:
1. locus_mod_gbk.py
2. common_name_mod_gbk.py
3. hypothetical_mod_gbk.py
4. gene_symbol_mod_gbk.py
5. ec_numbers_mod_gbk.py
6. gbk2tbl.py