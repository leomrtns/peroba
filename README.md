<img src="recipe/peroba.png" height="100">

__Leonardo de Oliveira Martins<sup>1</sup>__
<br>
<sub>1. Quadram Institute Bioscience, Norwich Research Park, NR4 7UQ, UK</sub>

## Introduction
**peroba** is a tool for continuous updating of outbreak information as new sequences are
incorporated. 
It is been developed as a phylogenetic tracking tool to aggregate samples sequenced at the QIB 
with global information from [COG-UK](https://www.cogconsortium.uk/) and [GISAID](https://www.gisaid.org/). 
Therefore you will not find any real data here, although [all COG-UK data are available
online](https://www.cogconsortium.uk/data/).
If you find any report/results/data here, it will be random/rubish (due to privacy reasons) and cannot be used or
interpreted. 

This tool is not useful (yet) to the scientific community at large, you may need to be familiar and have access to the 
COGUK consortium to make sense of some variables.
Peroba is under active testing and development, being employed at the QIB but with hope others may find it useful.
If you are looking for more stable COGUK-related tools, please see the ones available at 
[https://github.com/COG-UK](https://github.com/COG-UK) (in particular [civet](https://github.com/COG-UK/civet) or 
[phylo reports](https://github.com/COG-UK/phylo-reports)) and [https://github.com/cov-lineages](https://github.com/cov-lineages).


**peroba** is the name of an [endangered Brazilian timber tree](https://en.wikipedia.org/wiki/Aspidosperma_polyneuron).
But if you like acronyms it stands for Phylogenetic Epidemiology with ROBust Assignment. 

## Modules
**peroba** is composed of three modules, that should be run in order:
1. **`peroba_database`**: This script collects information from several sources and generates a set of `perobaDB` files.
2. **`peroba_backbone`**: This script selects a set of "global" sequences from `perobaDB` to be analysed together with the local ones 
(`NORW`). It finds local sequences within the database, but the user should also include other local sequences. 
3. **`peroba_report`**: once the user finishes the analysis (i.e. has a phylogenetic tree using suggestions from
   `peroba_backbone`), this script will estimate ancestral states and generate a PDF report.

## Installation

Before installing peroba, you will need to download and copy the shapefiles for plotting the maps, which we cannot
distribute here due to copyright issues.

### Requirements

* `conda`
* `texlive`
* linux 
* python > 3.6 
* internet access to download shapefiles

### Download shapefiles
Shapefiles can be downloaded however, and the postcode shapefiles are kindly provided by [OpenDoorLogistics](https://www.opendoorlogistics.com) (please check
[their license terms](https://www.opendoorlogistics.com/data)):
```bash
wget https://www.opendoorlogistics.com/wp-content/uploads/Data/UK-postcode-boundaries-Jan-2015.zip
unzip  UK-postcode-boundaries-Jan-2015.zip -d postcodes
cp postcodes/Distribution/Districts.* ${perobadir}/peroba/data/
```
Where `${perobadir}` is the root directory of your `peroba` installation. 
The directory `${perobadir}/peroba/data` should already exist when you cloned this repo.

Likewise, the `adm2` location correspond to NUTS 2 regions, and can be downloaded from
[GADM](https://gadm.org/download_country_v3.html):
```bash
wget https://biogeo.ucdavis.edu/data/gadm3.6/shp/gadm36_GBR_shp.zip
unzip gadm36_GBR_shp.zip -d adm2
cp adm2/gadm36_GBR_2.* ${perobadir}/peroba/data/
```

You will also need to download by hand the Sars-cov2 sequence and metadata files, which is not covered here.
I have never tried, but it might work to download the publicly available data on [COG-UK](https://github.com/COG-UK/data).
### Generate a conda environment

This software depends on several other packages, installable through conda or pip.
The suggested installation procedure is to create a conda environment (to take care of dependencies) and then installing
the python package:
```bash
conda update -n base -c defaults conda # probably not needed, but some machines complained about it
conda env create -f environment.yml  
conda activate peroba
python setup.py install
```

Since this software is still under development, these two commands are quite useful:
```bash
conda env update -f environment.yml # updat conda evironment after changing dependencies
pip install -e . # install  # installs in development mode (modifications to python files are live)
```

There are two system packages that might need to be installed outside conda, `libGL` and `texlive`:
```
apt-get install libgl1-mesa-glx texlive-full
```
`ete3` complains about missing `lbGL` [which is safer to install system-wise](https://github.com/conda-forge/pygridgen-feedstock/issues/10);
To test it, type `from PyQt5 import QtGui` in a python console and see if you are good to go.
And `texlive` is for the PDF report generation (`texlive-full` is a monster, but you won't need to worry about missing
fonts again :D )


The report generation relies on the [Eisvogel latex template for pandoc](https://github.com/Wandmalfarbe/pandoc-latex-template), 
which is included here (it's released under a BSD 3-clause).
The complete list of dependencies is described in the file [environment.yml](./environment.yml).
Please let me know if there are missing dependencies, although `peroba` is under active development and its behaviour
may change without notice. 


## Instructions
You can find [a tutorial on using the software here](docs/023.peroba_pipeline.ipynb).

## Caveats

When interpreting any result, please remember that not all sequences pass the sequencing quality control. 
Those that do may be excluded from COGUK phylogenetic analysis,
which means we won't have metadata (in particular `sequence_name`, which allows mapping between tree, sequence, and epi
data) information from them. 
We minimise this by using local information whenever possible, but still the reasons for exclusion remain.

This is not a general-purpose software. 
It is being released publicly in the hope that other researchers can build upon it. 

## License 
SPDX-License-Identifier: GPL-3.0-or-later

Copyright (C) 2020-today  [Leonardo de Oliveira Martins](https://github.com/leomrtns)

peroba is free software; you can redistribute it and/or modify it under the terms of the GNU General Public
License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later
version (http://www.gnu.org/copyleft/gpl.html).
