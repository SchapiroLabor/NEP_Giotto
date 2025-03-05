# Giotto for NEP comparison
Giotto from Dries et al. 2021 for SCNA. https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02286-2

# Introduction
Giotto is a suite for many spatial analysis directions, while we are interested in the cell-cell interaction part. Giotto identifies distinct cell type(ct)/ct interaction patterns by evaluating the enrichment of the frequency that each pair of cell-types is proximal to each other. For NBH definition, you can either use a kNN or Delaunay neighborhood network. The ratio of observed over expected frequencies between two ct is calculated by permutation testing. P-values are calculated by observing how often the observed values are higher or lower than expected. A wrapper of the function is written in cellProximityEnrichment() of the package.

# Installation

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


# Usage

Set input and output paths at the beginning of the scripts. 
The script scr/Giotto_SCNA.R was used to create the results on simulated data (asymmetric and symmetric) in the manuscript Schiller et al. (2025) on NEP analysis comparison. 







