

# initially cleaning of ensembles


# get predictions and untar them
wget "https://predictioncenter.org/casp16/assessors/TARBALLS/predictions/R1260" --recursive --no-parent --level=2 --user assessor --password Casp16MexicO -nH --cut-dirs=4
rm robots.txt
rm -r index.html* hybrid/ oligo/ QA/ regular/ RNA/ selected/ ligands/
cd R1260
rm index*
for tar in *tgz;
do
  tar -xf $tar
  rm $tar
done
cd ..


# rename these inproperly named pdbs
mv R1260/R1260TS417/R1260TS417_1_checked_reordered.pdb R1260/R1260TS417/R1260TS417_1.pdb
mv R1260/R1260TS183/R1260TS183_1_checked_reordered.pdb R1260/R1260TS183/R1260TS183_1.pdb

# convert pdbs to csvs, cleaning up the predictions to a standard format
mkdir R1260_csvs
python scripts/0-casp16_water-pdb_to_csv.py \
    --input_folder R1260 \
    --output_folder R1260_serial \
    --min_dist 5 \
    --csv_output summary_atom_counts.csv

# similiar for the baselines
python scripts/0-casp16_water-pdb_to_csv.py \
    --input_folder base_line \
    --output_folder R1260_serial \
    --min_dist 5 \
    --csv_output summary_atom_counts_baseline.csv

# check all these contents
python scripts/1-casp16_water-check-atom-contents.py \
    --native_pdb 9CBU_cleaned.pdb \
    --input_folder R1260_serial \
    --native_output scripts/9CBU-rna.csv

