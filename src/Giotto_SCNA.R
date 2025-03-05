#####################################################
### Script for running Giotto for simulated data  ###
### Schiller et at. 2025.                         ###
### Author: Chiara Schiller                       ###
#####################################################

# Load packages
library(Giotto)
library(reticulate)
library(readr)
library(tidyverse)
library(foreach)
library(doParallel)
library(doRNG)
library(here)
library(pheatmap)


### Data paths
# Define paths
data_path <- here::here("../../../../data/20250217_sym00_nbh2_1000dim_grid200_300iter_50swaps")
output_path <- here::here("../../../Comparison/20250218_results_sym/Giotto_delaunay_4ct_self00.csv")
# Get list of CSV files
files <- list.files(data_path, pattern = ".csv", full.names = TRUE)

### Configuration
# Set seeds
set.seed(42)
reticulate::py_run_string("import numpy as np; np.random.seed(42)")
# Giotto and Python configurations
python_path <- "/Users/chiaraschiller/miniconda3/bin/python"
my_instructions <- createGiottoInstructions(python_path = python_path)
Sys.setenv(RETICULATE_PYTHON = python_path)
reticulate::py_config()
# Set up parallel processing
num_cores <- detectCores() - 1  # Use one less than available cores
cl <- makeCluster(num_cores)
registerDoParallel(cl)
registerDoRNG(42)

### Parallelized NEP analysis with Giotto and Dleaunay nbh def
process_file <- function(file_path, index) {
  set.seed(42 + index)  # Ensure reproducibility
  
  data_list <- read_csv(file_path)
  df <- tibble(M1 = 1, M2 = 1, cell_ID = seq_len(nrow(data_list)))
  
  metadata <- tibble(ct = as.character(data_list$ct), cell_ID = seq_len(nrow(data_list)))
  
  giotto_obj <- createGiottoObject(
    raw_exprs = as.data.frame(t(df)),
    spatial_locs = select(data_list, -ct),
    cell_metadata = metadata
  )
  
  giotto_obj <- normalizeGiotto(gobject = giotto_obj)
  
  giotto_network <- createSpatialNetwork(
    gobject = giotto_obj, 
    minimum_k = 2, 
    method = "Delaunay"
  )
  
  cell_proximities <- cellProximityEnrichment(
    gobject = giotto_network,
    cluster_column = "ct",
    spatial_network_name = "Delaunay_network",
    adjust_method = "fdr",
    number_of_simulations = 100
  )
  
  cell_PI <- separate(cell_proximities[[2]], unified_int, into = c("label1", "label2"), sep = "--")
  cell_PI$sample <- basename(file_path)  # Store sample name
  return(cell_PI)
}

# Run the function in parallel
cell_PI_list <- foreach(i = seq_along(files), .packages = c("Giotto", "reticulate", "tidyverse")) %dopar% {
  process_file(files[i], i)
}

# Stop parallel cluster
stopCluster(cl)

### create standard matrix for output comparison
# Combine results into a single dataframe
data <- bind_rows(cell_PI_list) %>%
  mutate(sample = sub(".csv$", "", sample)) %>%
  select(label1, label2, PI_value, sample) %>%
  mutate(key = paste(label1, label2, sep = "_")) %>%
  select(-label1, -label2) %>%
  spread(key = key, value = PI_value)

# Set row names and replace NA values with 0
rownames(data) <- data$sample
data <- data %>% select(-sample) %>% mutate_all(~ replace_na(.x, 0))

# Save results
write.csv(data, output_path, row.names = TRUE)

### Heatmap
pheatmap(data)

sessionInfo()


