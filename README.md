# Giotto for NEP comparison
Giotto from Dries et al. 2021 for NEP analysis. https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02286-2

Giotto is a suite for many spatial analysis directions, while we are interested in the cell-cell neighbor preference (NEP) part. Giotto identifies distinct cell type(ct)/ct interaction patterns by evaluating the enrichment of the frequency that each pair of cell-types is proximal to each other. For NBH definition, you can either use a kNN or Delaunay neighborhood network. The ratio of observed over expected frequencies between two ct is calculated by permutation testing. P-values are calculated by observing how often the observed values are higher or lower than expected. A wrapper of the function is written in cellProximityEnrichment() of the package.

# Usage

## Installation

The Rpackage "Giotto" (Version: Giotto_1.1.2) was installed from the RubD/Giotto github repo.
Installation is described here: https://rubd.github.io/Giotto_site/ and copied below

`library(devtools)  # if not installed: install.packages('devtools')
library(remotes)  # if not installed: install.packages('remotes')
remotes::install_github("RubD/Giotto")`

The installed packages and versions can be found in `session_info.txt`

## miniconda env

For running, you need any miniconda env that is working for Giotto. The authors propose two ways here: https://rubd.github.io/Giotto_site/articles/tut0_giotto_environment.html (28th Feb 2023).
The second options creating a Giotto env did not work for me, so I created another miniconda env (python version 3.8). You don't need to install any libraries, just set the initial lines of the script to the python version. 

`my_instructions = createGiottoInstructions(python_path = 'your_env/bin/python')
Sys.setenv(RETICULATE_PYTHON = "your_env/bin/python") 
py_config()`

## Data

### In silico tissue (IST) data
Simulated .csv data with x, y, and ct annotation columns were used. The asymmetric and symmetric in silico tissue (IST) datasets were generated as described here: https://github.com/SchapiroLabor/NEP_IST_generation. 

# Scripts

`/src`:
- `NEP_Giotto.R`: This script runs Giotto on the simulated data (symmetric or asymmetric dataset) with a neighborhood definition of KNN=5. We output the CPScores. 

