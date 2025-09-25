# OLE R1285 analysis

## motif analysis

To get motif information the following was run
```
rna_motif -s OLE.pdb > OLE.pdb.rna_motif.out

tar -xf R1285o.tar.gz

cd ../R1285o
for f in R1285TS???_?o; do mv $f ${f}.pdb; done
for f in R1285TS???_?o.pdb; do echo $f; rna_motif -s ${f} > ${f}.rna_motif.out; done
```
Three PDBs errored: 
```
R1285TS167_1o
R1285TS208_2o
R1285TS450_1o
```

(process_data.py)[process_data.py] compares the native interactions in (native)[native] to predicted motifs in (R1285o)[R1285o] to produce the summary (interaction_results.csv)[interaction_results.csv].

(all_intermolecular_interaction.py)[all_intermolecular_interaction.py] reads all predicted motifs in (R1285o)[R1285o] and outputs a summary, (all_intermolecular_interactions.csv)[all_intermolecular_interactions.csv].


Visuals are found in (images)[images] and the pymol session, (R1285.pse)[R1285_pymol.zip].