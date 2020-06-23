<img src="recipe/peroba.png" height="100">

__Leonardo de Oliveira Martins<sup>1</sup>__
<br>
<sub>1. Quadram Institute Bioscience, Norwich Research Park, NR4 7UQ, UK</sub>

## Introduction
**peroba** is a tool for continuous updating of outbreak information as new sequences are
incorporated. 
It is been developed as a phylogenetic tracking tool to aggregate samples sequenced at the QIB 
with global information from [COG-UK](https://www.cogconsortium.uk/) and [GISAID](https://www.gisaid.org/). 
Therefore you will not find any data here, although [all COG-UK data are available
online](https://www.cogconsortium.uk/data/).

**peroba** is the name of a Brazilian timber tree, but if you like acronyms it stands 
for Phylogenetic Epidemiology with ROBust Assignment. 

## Modules
**peroba** is composed of three modules, that should be run in order:
1. **`peroba_database`**: This script collects information from several sources and generates a set of `perobaDB` files.
2. **`peroba_backbone`**: This script selects a set of "global" sequences from `perobaDB` to be analysed together with the local ones 
(`NORW`). It finds local sequences within the database, but the user should also include other local sequences. 
3. **`peroba_report`**: once the user finishes the analysis (i.e. has a phylogenetic tree using suggestions from
   `peroba_backbone`), this script will estimate ancestral states and generate a PDF report.

## Instructions
You can find [a tutorial on the software here](docs/023.peroba_pipeline.html).

## Caveats
Not all sequences pass the sequencing quality control. Those that do may be excluded from COGUK phylogenetic analysis,
which means we won't have metadata (importantly, `sequence_name` which allows mapping between tree, sequence, and epi
data) information from them. 
We minimise this by using local information whenever possible, but still the reasons for exclusion remain. 

## License 
SPDX-License-Identifier: GPL-3.0-or-later

Copyright (C) 2020-today  [Leonardo de Oliveira Martins](https://github.com/leomrtns)

peroba is free software; you can redistribute it and/or modify it under the terms of the GNU General Public
License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later
version (http://www.gnu.org/copyleft/gpl.html).
