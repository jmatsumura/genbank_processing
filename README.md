# genbank_processing
Scripts for reformatting GenBank files given a metadata file. Will
do the following:

- Remove overlapping genes (not *RNAs)
- Correct gene symbol syntax for prokaryotic genes
- Remove contigs that are <200BP and reorder the others accordingly
- Replace locus tags
- Update EC numbers
- Remove hypothetical gene symbols
