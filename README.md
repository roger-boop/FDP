# FDP

Welcome to the documentation for the **Identification of the recurring domain interaction motifs**! This README provides an overview of the various scripts included in the project and instructions on how to use them.

## Table of Contents
- [Introduction](#introduction)
- [Scripts](#scripts)
- [Dependencies](#dependencies)

## Introduction

These are the scripts used in my FDP project.

## Scripts

Here is a list of scripts included in this project and a brief description of their functionalities:

1. `steps.sh`: Shell script with general steps. It executes all the python files in the correct order.
2. `steps.PDB.sh`: The `steps.sh` version when using the PDB dataset.
3. `pdb2fasta.py`: Takes a directory containing the structures and generates a fasta format file (_multifasta.fa_). Used by `steps.sh` and `steps.PDB.sh`.
4. `domtblParser.py`: Parses the domtbl file (the output from hmmscan). Used by `steps.sh` and `steps.PDB.sh`.
5. `AlphaFoldParser.py`: Parses AlphaFold structure prediction files (pdb format). Used by `steps.sh`.
6. `AlphaFoldParser.PDB.py`: Adaptation from `AlphaFoldParser.py` that parses PDB structure files (pdb format). Used by `steps.PDB.sh`.
7. `domainInfo.py`: Middle step that parses data from a shelve file to a new one. Used by `steps.sh` and `steps.PDB.sh`.
8. `peptides_interaction.py`: Analyzes interactions involving peptides and generates a _csv_ file containing the results. Used by `steps.sh` and `steps.PDB.sh`.
7. `domain_separator.py`: Creates a folder for each domain and saves a structure for each domain prediction (only the section that should be the domain) with data from the hmmscan.
9. `ce_roger.py`: It finds a reference structure for each domain.
10. `calculate_rmsd.py`: Calculates RMSD (Root Mean Square Deviation) between protein structures. Used by `ce_roger.py`.

12. `pdbview.ipynb`: Jupyter Notebook for viewing PDB structures.
13. `sprot_check.py`: Checks protein information in SwissProt database.
14. `sprot_check.R`: R script for protein information checking.
15. `sprotParser.py`: Parses SwissProt files.
18. `heatmap.R`: R script for generating heatmaps.
19. `analysis.R`: R script for analysis (please provide description).

6. `domain_interpretation.py`: Interprets protein domain information from the hmmscan results. Analysis file.

## Dependencies
This scripts are used in machine with ubuntu 22.10. You should have python and R installed with the following packages:
- Python:
  - Biopython
  - NetworkX
  - Pyvis.network
  - shelve
  - gzip
  - matplotlib
  - shutil
  - subprocess
  - csv
  - webcolors
  - nglview
  - piywidgets
- R:
  - dplyr
  - ggplot2
  - plotly
  - stringr
