# load packages
library(Giotto)
library(reticulate)
library(readr)
library(tidyverse)
library(foreach)
library(doParallel)
library(doRNG)

#### This script uses simulated data to run Giotto and save results as interaction by sample matrix ####

# Reproducibility
set.seed(42) 
reticulate::py_run_string("import numpy as np; np.random.seed(42)")

# Giotto and python configurations
my_instructions = createGiottoInstructions(python_path = '/Users/chiaraschiller/miniconda3/bin/python')
Sys.setenv(RETICULATE_PYTHON = "/Users/chiaraschiller/miniconda3/bin/python") 
reticulate::py_config()

# paths and files
files = list.files("./../../../../data/Sim_nbh2_asym01_1000_grid0.2_300iter_50swap/", pattern = ".csv")
data_path = "./../../../../data/Sim_nbh2_asym01_1000_grid0.2_300iter_50swap/"
output_path = "./../../../Comparison/20250218_results_asym/Giotto_delaunay_4ct_cross01.csv"

data_list = list()

## Paralellize
# Set the number of cores to use
num_cores <- 4
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Use doRNG to ensure reproducibility in parallel execution
registerDoRNG(42)  

# create dataframe with random values to imitate expression data to generate a giotto object

foreach(i = 1:length(files), .packages = c("Giotto","reticulate", "tidyverse")) %dopar%{
  # Set seed for reproducibility inside the loop
  set.seed(42 + i)
  data_list <- readr::read_csv(paste0(data_path,files[i]))
  df = data.frame(cbind(rep(1,nrow(data_list)), rep(1,nrow(data_list))))
  df$cell_ID = seq(1, nrow(df), by = 1)
  colnames(df) = c("M1", "M2", "cell_ID")

  metadata = cbind(as.character(data_list$ct), seq(1, nrow(df), by = 1))
  colnames(metadata) = c("ct", "cell_ID")

  data = createGiottoObject(raw_exprs = as.data.frame(t(df)),
                          spatial_locs = as.data.frame(data_list[, !(names(data_list) %in% c("ct"))]),
                          cell_metadata = metadata)

  data <- normalizeGiotto(gobject = data)
  
  # create network (required for binSpect methods)
  giotto_list = createSpatialNetwork(gobject = data, minimum_k = 2,
                                                  method = "Delaunay",
                                                  #method = "kNN",
                                                  #k = 8,
                                                  )
  
  cell_proximities = cellProximityEnrichment(gobject = giotto_list,
                                             cluster_column = 'ct',
                                             spatial_network_name = 'Delaunay_network',
                                             #spatial_network_name = 'kNN_network',
                                             adjust_method = 'fdr',
                                             number_of_simulations = 100)
  
  cell_PI= separate(cell_proximities[[2]], unified_int, into = c("label1", "label2"), sep = "--")
} -> cell_PI_list

for(i in 1:length(cell_PI_list)){
  cell_PI_list[[i]]$sample = rep(files[i], nrow(cell_PI_list[[i]]))
}

## create matrix for comparison
data = do.call("rbind", cell_PI_list)
data$sample <- sub(".csv$", "", data$sample)
data = data %>% select(label1, label2, PI_value, sample)
data$key = paste(data$label1, data$label2, sep = "_")
data = data %>% select(-c(label1, label2))


data = spread(data, key = key, value = PI_value)
rownames(data) = data$sample
data = data[,-1]

data <- data.frame(lapply(data, function(x) ifelse(is.na(x) | is.nan(x), 0, x)))
data = as.data.frame(data)
rownames(data) = sub(".csv", "", files )

write.csv(data, output_path, row.names = TRUE)

#df_heatmap = data %>% select(`Endocardial.cells_Mono...Macros.Ccr2.`, Endocardial.cells_Neutrophils)

library(pheatmap)
pheatmap(data)

sessionInfo()



# Load necessary packages
library(Giotto)
library(reticulate)
library(readr)
library(tidyverse)
library(foreach)
library(doParallel)
library(doRNG)
library(here)
library(pheatmap)

# Set up reproducibility
set.seed(42)
reticulate::py_run_string("import numpy as np; np.random.seed(42)")

# Giotto and Python configurations
python_path <- "/Users/chiaraschiller/miniconda3/bin/python"
my_instructions <- createGiottoInstructions(python_path = python_path)
Sys.setenv(RETICULATE_PYTHON = python_path)
reticulate::py_config()

# Define paths
data_path <- here::here("../../../../data/20250217_sym00_nbh2_1000dim_grid200_300iter_50swaps")
output_path <- here::here("../../../Comparison/20250218_results_sym/Giotto_delaunay_4ct_cross01.csv")

# Get list of CSV files
files <- list.files(data_path, pattern = ".csv", full.names = TRUE)

# Set up parallel processing
num_cores <- detectCores() - 1  # Use one less than available cores
cl <- makeCluster(num_cores)
registerDoParallel(cl)
registerDoRNG(42)

# Function to process each file
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

# Plot heatmap
pheatmap(data)

# Print session info
sessionInfo()


