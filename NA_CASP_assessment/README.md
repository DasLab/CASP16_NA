# CASP16 NA assessment

> Kretsch, R. C.; Hummer, A. M.; He, S.; Yuan, R.; Zhang, J.; Karagianes, T.; Cong, Q.; Kryshtafovych, A.; Das, R. Assessment of Nucleic Acid Structure Prediction in CASP16. bioRxiv, 2025, 2025.05.06.652459. https://doi.org/10.1101/2025.05.06.652459.


## Monomer

CASP monomer analysis can be found at [analysis/NA_monomer.ipynb](analysis/NA_monomer.ipynb). CASP caculated scores can be found in [raw_scores/RNA_monomer_results_table.csv](raw_scores/RNA_monomer_results_table.csv) with summary scores found at [analysis/summary_tables/final_Z_monomer.csv](analysis/summary_tables/final_Z_monomer.csv). Chimerax sessions can be found in [visuals/cite.cxc](visuals/cite.cxc), [visuals/dna.cxc](visuals/dna.cxc), [visuals/ole.cxc](visuals/ole.cxc), [visuals/R1205.cxc](visuals/R1205.cxc), [visuals/R1288.cxc](visuals/R1288.cxc).


## Ligand

CASP caculated scores, based on OST, can be found in [raw_scores/ligand](raw_scores/ligand). Figure code can be found at [analysis/NA_ligand.ipynb](analysis/NA_ligand.ipynb). Pymol sessions can be found in [visuals/D1273_ligand.pse](visuals/D1273_ligand.pse], [visuals/R1288_ligand.pse](visuals/R1288_ligand.pse], and [visuals/R1261_ligand_images.pse](visuals/R1261_ligand_images.pse). Overview of scores can be found in [analysis/summary_tables/ligand.csv](analysis/summary_tables/ligand.csv) with final Z sum scores in [analysis/summary_tables/final_Z_NA_ligand.csv](analysis/summary_tables/final_Z_NA_ligand.csv).

## RNA multimer

Calculation of stiochemtry and symmetry can be found at [raw_scores/stochiometry_and_symmetry](raw_scores/stochiometry_and_symmetry). CASP caculated scores can be found in [raw_scores/RNA_multimer_results_table.csv](raw_scores/RNA_multimer_results_table.csv). Figure code can be found at [analysis/NA_multimer.ipynb](analysis/NA_multimer.ipynb). Chimerax sessions can be found in [visuals/R1251.cxc](visuals/R1251.cxc). Final Z sum scores in [analysis/summary_tables/final_Z_NA_multimer.csv](analysis/summary_tables/final_Z_NA_multimer.csv).

## Hybrid

Calculation of stiochemtry can be found at [raw_scores/stochiometry_and_symmetry](raw_scores/stochiometry_and_symmetry).CASP caculated scores can be found in [raw_scores/RNA_hybrid_results_table.csv](raw_scores/RNA_hybrid_results_table.csv). Figure code can be found at [analysis/NA_hybrid.ipynb](analysis/NA_hybrid.ipynb). The code to seperate these scores by interaction time [eg protein-protein v RNA-protein) is in [raw_scores/RNP](raw_scores/RNP). A summary of scores can be found at [analysis/summary_tables/hybrid_scores.csv](analysis/summary_tables/hybrid_scores.csv) with final Z sum scores in [analysis/summary_tables/final_Z_NA_hybrid.csv](analysis/summary_tables/final_Z_NA_hybrid.csv). ChimeraX sessions to visualize prediction of M1221, M1224, M1276, M1212 can be found in [visuals](visuals).

## Base pairing

Base pairing analysis can be found in [analysis/NA_secstruct.ipynb](analysis/NA_secstruct.ipynb). Theidentification of base pairs using can be found in [raw_scores/secondary_structure](raw_scores/secondary_structure). The results of the CASP16 predictors is found at [analysis/summary_tables/CASP16_3Dpreds_SS.parquet](analysis/summary_tables/CASP16_3Dpreds_SS.parquet), [analysis/summary_tables/CASP16_3Dpreds_SS_o.parquet](analysis/summary_tables/CASP16_3Dpreds_SS_o.parquet), and [analysis/summary_tables/CASP16_3Dpreds_SS_NC.parquet](analysis/summary_tables/CASP16_3Dpreds_SS_NC.parquet). The results of the secondary structure predictors is found at [analysis/summary_tables/CASP16_2Dpreds_SS.parquet](analysis/summary_tables/CASP16_2Dpreds_SS.parquet)

## Motif analysis

Scripts to obtain motif annotation can be found at [raw_scores/motif](raw_scores/motif). Figure code can be found at [analysis/NA_motif.ipynb](analysis/NA_motif.ipynb). F1 scores for each group, target, and motif type can be found in [analysis/summary_tables/motif_monomer.csv](analysis/summary_tables/motif_monomer.csv) and [analysis/summary_tables/motif_multimer.csv](analysis/summary_tables/motif_multimer.csv). Summary by group-target and by group can be found at [analysis/summary_tables/motif_performance.csv](analysis/summary_tables/motif_performance.csv) and [analysis/summary_tables/motif_summ.csv](analysis/summary_tables/motif_summ.csv) respectively.

## Overall performance

The perfromance over categories is summarized and the perfromance over time is dispalyed in [analysis/overall_performance.ipynb](analysis/overall_performance.ipynb). RNA puzzle and CASP15 performance is found in [raw_scores/all_rnapz_model_eval.csv](raw_scores/all_rnapz_model_eval.csv) and [raw_scores/casp15_rna_all_score_per_group_withZ_230403.csv](raw_scores/casp15_rna_all_score_per_group_withZ_230403.csv) respectively. Labling of servers for RNA puzzles can be found at [raw_scores/rnapuzzles_human_webserver_stastic.csv](raw_scores/rnapuzzles_human_webserver_stastic.csv). The histotical overview and pairs of indepdent expiremental structures can be found at [raw_scores/historical_overview.csv](raw_scores/historical_overview.csv) and [raw_scores/IndepdentStructureComparisons.csv](raw_scores/IndepdentStructureComparisons.csv). 

## Other information

- Supplemental table can be found at [metric_scores/Supplemental_Tables_CASP16_NA_assessment.xlsx](metric_scores/Supplemental_Tables_CASP16_NA_assessment.xlsx). Note in the all_scores tab, the entries without a casp group number, eg secondary structure predictors and AF3 with and without MSA were run after CASP and are not participants.
- Group annotation can be found at [raw_scores/CASP16_groups.csv](raw_scores/CASP16_groups.csv) and targets annotations can be found in [raw_scores/casp16_targets.csv](raw_scores/casp16_targets.csv) and [raw_scores/target_names.csv](raw_scores/target_names.csv).
- The list of residues unresolved in the expiremental structures can be found at [raw_scores/unresolved_residues.csv](raw_scores/unresolved_residues.csv)
