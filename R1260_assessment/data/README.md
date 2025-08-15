# Data for R1260

## Data locations

__Predictions:__ predictions of atomic ensembles from 28 prediction groups can be downloaded from https://predictioncenter.org/download_area/CASP16/predictions/R1260/.

__Atomic model derived maps:__ The predicted ensembles were converted to maps for direct comparison against the cryo-EM maps. These predicted maps can be found XXXX.

__Base line ensembles__: Base line ensembles are derived from selection 1,000 models from the various simulations in [Kretsch et al 2025](https://doi.org/10.1038/s41586-025-08855-w). They can be found at https://doi.org/10.25740/vz022mq8177. The full simulations are found at https://doi.org/10.25740/sw275qs6749.

__Cryo-EM maps:__ The cryo-em maps used as targets in R1260 can be found in [EMD-42499](https://www.ebi.ac.uk/emdb/EMD-42499) and [EMD-42498](https://www.ebi.ac.uk/emdb/EMD-42498).

__Cryo-EM derived models:__ The above cryo-EM maps were modeled in the traditional sense of placing water and magnesium ion in well-defined peaks: PDB [9CBU](https://doi.org/10.2210/pdb9CBU/pdb), [9CBX](https://doi.org/10.2210/pdb9CBX/pdb), [9CBW](https://doi.org/10.2210/pdb9CBW/pdb), and [9CBY](https://doi.org/10.2210/pdb9CBY/pdb). Due to the various biases in modeling, predicted ensembles were compared against the cryo-EM maps, not these models.

__Locally aligned neighborhoods:__ Intermediates in the process of converting atomic models to densities is the conversion to locally aligned neighborhoods which can be found XXXX.



## Conversion of atomic models to densities



### 1. Cleaning of the models

This was run using the script [1-casp16_water-clean_ensembles.sh](1-casp16_water-clean_ensembles.sh), which included the use of these python scripts [scripts/0-casp16_water-pdb_to_csv.py](scripts/0-casp16_water-pdb_to_csv.py) and [scripts/1-casp16_water-check-atom-contents.py](scripts/1-casp16_water-check-atom-contents.py). All models were untarred, cleaned including remaining residues to common names, and reorganized for downstream processing. At this time all atoms > 5 Å from a RNA heavy-atom were removed. The number of models submitted and the number of water and ions submitted were recorded. These cleaned ensemble files can be found in the locally aligned neighborhoods files _R1260_serial.tar.xz_.

### 2. Creation of locally aligned neighborhoods

Neighborhoods were created and ensembles were split into these neighborhoods using script [2-casp16_water-get_neighborhoods.sh](2-casp16_water-get_neighborhoods.sh) which uses the python file [scripts/2-casp16_water-get_areas_of_interest.py](scripts/2-casp16_water-get_areas_of_interest.py).

Neighborhoods around each RNA residue were decided based on 9CBU. Any residues that had a heavy-atom within X Å in 9CBU were included in the neighborhood. X = 6, 10, 12, 20 were all tested. At this time predicted ensembles were also split into neighborhoods. For each RNA residues, the previously defined neighborhood of RNA residues was extracted, in addition to any water and ions within 5 Å. Hence, each atom is likely in multiple neighborhoods. These neighborhoods have been uploaded in the link above.

_Note, if a model was missing any RNA atoms (eg residue 408 was missing in a handful of the groups predictions) in the neighborhood, that model was not included in that neighborhood._

### 3. Conversion of atoms the density and stitching of neighborhoods

From these neighborhoods, density were calculated and stitched together using the following script [3-casp16_water-get_densities.sh](3-casp16_water-get_densities.sh) which uses the python script in [scripts/3-casp16_water-get_density.py](scripts/3-casp16_water-get_density.py).