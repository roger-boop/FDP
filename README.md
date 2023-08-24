# FDP

Welcome to the documentation for the **Identification of the recurring domain interaction motifs**! This README provides an overview of the various scripts included in the project and instructions on how to use them.

## Table of Contents
- [Introduction](#introduction)
- [Scripts](#scripts)
- [Usage](#usage)
- [Dependencies](#dependencies)

## Introduction

These are the scripts used in my FDP project.

## Scripts

Here is a list of scripts included in this project and a brief description of their functionalities:

1. `steps.sh`: Shell script with general steps. It executes all the python files in the correct order.
2. `steps.PDB.sh`: The `steps.sh` version when using the PDB dataset.
3. `pdb2fasta.py`: Takes a directory containing the structures and generates a fasta format file (_multifasta.fa_).
1. `AlphaFoldParser.py`: Parses AlphaFold structure prediction files (pdb format).
2. `AlphaFoldParser.PDB.py`: Adaptation from `AlphaFoldParser.py` that parses PDB structure files (pdb format).
3. `ce_roger.py`: CE calculations using the ROGER method.
4. `calculate_rmsd.py`: Calculates RMSD (Root Mean Square Deviation) between protein structures.
5. `domainInfo.py`: Provides information about protein domains.
6. `domain_interpretation.py`: Interprets protein domain information.
7. `domain_separator.py`: Separates domains from protein structures.
8. `domtblParser.py`: Parses domtbl files.
9. `motif3DCheck.py`: Checks 3D motifs in protein structures.
10. `peptides_interaction.py`: Analyzes interactions involving peptides.
12. `pdbview.ipynb`: Jupyter Notebook for viewing PDB structures.
13. `sprot_check.py`: Checks protein information in SwissProt database.
14. `sprot_check.R`: R script for protein information checking.
15. `sprotParser.py`: Parses SwissProt files.
18. `heatmap.R`: R script for generating heatmaps.
19. `analysis.R`: R script for analysis (please provide description).

## Usage

Provide instructions on how to use each script. You can break this down into steps, including any command-line arguments or input files that need to be specified.

### `AlphaFoldParser.PDB.py`

Explain the purpose of the script and how to use it:

1. **Installation**: Mention any setup or installation required, such as dependencies, libraries, or packages that need to be installed.
2. **Running the Script**: Provide the command to run the script, including any command-line arguments.

Example:
```sh
python AlphaFoldParser.PDB.py --input input_file.pdb --output output_file.txt
```

## Dependencies
