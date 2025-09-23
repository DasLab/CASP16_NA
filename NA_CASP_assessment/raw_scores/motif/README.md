# Motif analysis

Motif analysis was run using Rosetta rna_motif. Rosetta can be installed following [instructions here](https://docs.rosettacommons.org/demos/latest/tutorials/install_build/install_build). It can in theory be run via the [Rosie server](https://rosie.rosettacommons.org/rna_info).

rna_motif will output the ordered list of residues in each motif it identified in the stdout, so the relevant information can be saved via the command: `rna_motif -s $PDB > rna_motif_outputs/RNA_Data_rna_motif/${PDB}.stdout`.

This can than be read using code in (extract_rna_motif.py)[sets/extract_rna_motif.py]

Scoring and plots can be found in (NA_motif.ipynb)[../analysis/NA_motif.ipynb].