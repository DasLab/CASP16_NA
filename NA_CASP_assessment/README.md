# CASP16 NA assessment

    Kretsch, R. C.; Hummer, A. M.; He, S.; Yuan, R.; Zhang, J.; Karagianes, T.; Cong, Q.; Kryshtafovych, A.; Das, R. Assessment of Nucleic Acid Structure Prediction in CASP16. bioRxiv, 2025, 2025.05.06.652459. https://doi.org/10.1101/2025.05.06.652459.


## Ligand

CASP caculated scores, based on OST, can be found in (raw_scores/ligand)[raw_scores/ligand]. Figure code can be found at (analysis/NA_ligand.ipynb)[analysis/NA_ligand.ipynb]. Pymol sessions can be found in (visuals/D1273_ligand.pse)[visuals/D1273_ligand.pse], (visuals/R1288_ligand.pse)[visuals/R1288_ligand.pse], and (visuals/R1261_ligand_images.pse)[visuals/R1261_ligand_images.pse]. Overview of scores can be found in (analysis/summary_tables/ligand.csv)[analysis/summary_tables/ligand.csv] with final Z sum scores in (analysis/summary_tables/final_Z_NA_ligand.csv)[analysis/summary_tables/final_Z_NA_ligand.csv].

## Hybrid

CASP caculated scores can be found in (raw_scores/RNA_hybrid_results_table.csv)[raw_scores/RNA_hybrid_results_table.csv]. Figure code can be found at (analysis/NA_hybrid.ipynb)[analysis/NA_hybrid.ipynb]. The code to seperate these scores by interaction time (eg protein-protein v RNA-protein) is in (raw_scores/RNP)[raw_scores/RNP]. A summary of scores can be found at (analysis/summary_tables/hybrid_scores.csv)[analysis/summary_tables/hybrid_scores.csv] with final Z sum scores in (analysis/summary_tables/final_Z_NA_hybrid.csv)[analysis/summary_tables/final_Z_NA_hybrid.csv]. ChimeraX sessions to visualize prediction of M1221, M1224, M1276, M1212 can be found in [visuals](visuals).

## Motif analysis

Scripts to obtain motif annotation can be found at (raw_scores/motif)[raw_scores/motif]. Figure code can be found at (analysis/NA_motif.ipynb)[analysis/NA_motif.ipynb]. F1 scores for each group, target, and motif type can be found in (analysis/summary_tables/motif_monomer.csv)[analysis/summary_tables/motif_monomer.csv] and (analysis/summary_tables/motif_multimer.csv)[analysis/summary_tables/motif_multimer.csv]. Summary by group-target and by group can be found at (analysis/summary_tables/motif_performance.csv)[analysis/summary_tables/motif_performance.csv] and (analysis/summary_tables/motif_summ.csv)[analysis/summary_tables/motif_summ.csv] respectively.

## Other information

- Group annotation can be found at (raw_scores/CASP16_groups.csv)[raw_scores/CASP16_groups.csv].


### TODO REMOVE

### Secondary structure - metric calculation

This is analysis for only canonical base pairs, including the subset of singlets and cross pairs. The raw data from which this is pulled is found in this repository at `NA_assessment/raw_scores/secondary_structure/base_pair_tables/*canonical.csv`.

- The results of the CASP16 predictors is found at [summary_tables/CASP16_3Dpreds_SS.parquet](summary_tables/CASP16_3Dpreds_SS.parquet)
- The results of the secondary structure predictors is found at [summary_tables/CASP16_2Dpreds_SS.parquet](summary_tables/CASP16_2Dpreds_SS.parquet)
  - Secondary structure predictions that were the incorrect length were listed in [incorrect_length_db.csv](incorrect_length_db.csv). This will be addressed and this file removed. TODO

### Secondary structure - ranking

Based on the metrics calculated above, the ranking and figures were created as seen in [NA_motifs.ipynb](NA_motifs.ipynb).

This folder currently contains:

1. [raw_scores](raw_scores). Tables with the metric scores for each prediction including scripts for how to obtain these scores when applicable.
2. [analysis](analysis). Scripts for summarizing and analyzing the metrics including ranking
3. [chimerax_sessions](chimerax_sessions). Molecular visualizations
4. [NA_assessment_figures](NA_assessment_figures). The multi panel figures for publication.
5. [NA_assessment_tables](NA_assessment_tables). The tables produced for publication after some beautification.