# Base pair ananlysis

## Base pair annotation

DSSR outputs were generated from pdb files using the command

```x3dna-dssr -i={pdb} -o > output.txt```

The following code was run to parse DSSR output:

```bash
cp /oak/stanford/groups/rhiju/sherlock/home/shujun/casp16_dssr_outputs.zip .
unzip casp16_dssr_outputs.zip
rm casp16_dssr_outputs.zip

python process_dssr_outputs.py 
```


Then this was run to separate canonical and non-conical analysis:

```bash
python filter_cononical.py
```


The secondary structure algorithm predictions were in dot bracket notation and can be found in [combined_secondary_structure_prediction_algorithmns_rck_edited](combined_secondary_structure_prediction_algorithmns_rck_edited).

