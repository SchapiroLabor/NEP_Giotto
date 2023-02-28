# Giotto_SCNA
Giotto from Dries et al. 2021 for SCNA

# Introduction
Giotto is a suite for many spatial analysis directions, while we are interested in the cell-cell interaction part. Giotto identifies distinct cell type(ct)/ct interaction patterns ny evaluating the enrichment of the frequency that each pair of cell-types is proximal to each other. For NBH definition, you can either use a kNN or Delaunay neighborhood network. For the Delaunay network, they cite i-niche. Further, the ratio of observed over expected frequencies between two ct is calculated by permutating (Monte Carlo sampling?) expected frequencies. P-values are calculated by observing how often the observed values are higher or lower than expected. A wrapper of the fucntion is written in cellProximityEnrichment() of the package.

# Installation

Please clone this github repository. Input data are. csv files in a folder named "data" two levels above the folder in which the cloned repo is. The file structure are 3 columns, first x and y coordinates, then the "ct" annotation column.

## Dependencies
The Rpackage "Giotto" (Version: Giotto_1.1.2) was installed from the RubD/Giotto github repo.
Installation is described here: https://rubd.github.io/Giotto_site/ and copied below

`library(devtools)  # if not installed: install.packages('devtools')
library(remotes)  # if not installed: install.packages('remotes')
remotes::install_github("RubD/Giotto")`

Additionally, the packages reticulate (1.28), readr (2.1.3) and tidyverse (1.3.2) are required to be installed.

## miniconda env

For running, you need any miniconda env that is working for Giotto. The authors propose two ways here: https://rubd.github.io/Giotto_site/articles/tut0_giotto_environment.html (28th Feb 2023).
The second options creating a Giotto env did not work for me, so I created another miniconda env (python version 3.8). You don't need to install any libraries, just set the initial lines of the script to the python version. 

`my_instructions = createGiottoInstructions(python_path = 'your_env/bin/python')
Sys.setenv(RETICULATE_PYTHON = "your_env/bin/python") 
py_config()`

Replace these lines in the script so far, will create a better solution at another point.

# Usage

For analysing your data, make sure the data structure is as described above. An output table per file will be saved in an "output" folder on the same level as your cloned github repo.

For running Giotto, you can just run this line in the terminal:
`Rscript src/Giotto_SCNA.R`





