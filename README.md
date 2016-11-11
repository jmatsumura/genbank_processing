# genbank_processing
Scripts for reformatting GenBank files given a metadata file. Will
do the following:

- Remove overlapping genes (not *RNAs)
- Correct gene symbol syntax for prokaryotic genes
- Remove contigs that are <200BP and reorder the others accordingly
- Replace locus tags
- Update EC numbers
- Remove hypothetical gene symbols

## Steps (found explicitly in genbank_submission.sh):
1. clean_prokka.py 
2. locus_mod_gbk.py
3. common_name_mod_gbk.py
4. hypothetical_mod_gbk.py
5. gene_symbol_mod_gbk.py
6. ec_numbers_mod_gbk.py
7. gbk2tbl.py               - Round 1
8. prepare_sbt_cmt.pl
9. wrap_tbl2asn.pl          - Round 1
10. delete_overlap_mod_gbk.py
11. gbk2tbl.py              - Round 2
12. gene_symbol_mod_tbl.pl
13. format_velvet_assembled_AGP.pl
14. File rearrangement
15. wrap_tbl2asn.pl         - Round 2


Many of these scripts expect a metadata file with one row for each GBK entry to
reformatted and the following 24 columns present in each row:

gbk/file | locus tag | curate/common/names/file | delete/gene/symbols/file | 
bioproject ID | organism name | strain name | serotype | host | date of isolation
country | assembly method | coverage | sequencing method | contact person | 
email of contacts | authors | title | illegal/ec/numbers/file | isolation source
contact/list/file | biosample ID

Note that values can be null, but you must add a tab even for missing values.
