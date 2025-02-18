#Giotto usage for simulated data as input

library(Giotto)
library(reticulate)
library(readr)
library(tidyverse)
library(foreach)
library(doParallel)


#### This script uses simulated data to run Giotto and save results as interaction by sample matrix ####

# Giotto and python configurations
# Create Giotto path
#my_instructions = createGiottoInstructions(python_path = '/Users/chiaraschiller/Library/Library/r-miniconda-arm64/bin/python3')
#Sys.setenv(RETICULATE_PYTHON = "/Users/chiaraschiller/Library/r-miniconda-arm64/bin/python3") 

my_instructions = createGiottoInstructions(python_path = '/Users/chiaraschiller/miniconda3/bin/python')
Sys.setenv(RETICULATE_PYTHON = "/Users/chiaraschiller/miniconda3/bin/python") 
#py_config()

#installGiottoEnvironment()

# load simulated data

files = list.files("./../../../../data/Sim_nbh2_asym01_1000_grid0.2_300iter_50swap/", pattern = ".csv")
data_path = "./../../../../data/Sim_nbh2_asym01_1000_grid0.2_300iter_50swap/"

#files = list.files("./../../../../../data/Risom_breast/processed_files_celllineage/", pattern = ".csv")
#data_path = "./../../../../../data/Risom_breast/processed_files_celllineage/"

#files = list.files("./../../../../../MI_heart_paper/data/MI_heart_lunaphore_split_csv/", pattern = ".csv")
#data_path = "./../../../../../MI_heart_paper/data/MI_heart_lunaphore_split_csv/"

data_list = list()

## Paralellize
# Set the number of cores to use
num_cores <- 4
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# create dataframe with random values to pretend to have expression data to generate a giotto object
#for (i in data_list){
foreach(i = 1:length(files), .packages = c("Giotto","reticulate", "tidyverse")) %dopar%{
  data_list <- readr::read_csv(paste0(data_path,files[i]))
  df = data.frame(cbind(rep(1,nrow(data_list)), rep(1,nrow(data_list))))
  #df <- data.frame(matrix(runif(nrow(data_list[[i]]) * 2), nrow = nrow(data_list[[i]])))
  df$cell_ID = seq(1, nrow(df), by = 1)
  colnames(df) = c("M1", "M2", "cell_ID")

  metadata = cbind(as.character(data_list$ct), seq(1, nrow(df), by = 1))
  colnames(metadata) = c("ct", "cell_ID")
  # colbind the two datasets
  #data = cbind(df, as.data.frame(data))
  data = createGiottoObject(raw_exprs = as.data.frame(t(df)),
                          spatial_locs = as.data.frame(data_list[, !(names(data_list) %in% c("ct"))]),
                          #instructions = my_instructions,
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
#data$sample <- sub(".csv.[0-9]*$", "", data$sample)
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

#write.csv(data,"./../../../Comparison/results_4ct_asym_0.2grid_self/Giotto_delaunay_4ct_cross01.csv")
#write.csv(data,file=paste0("./../../../Comparison/results_Risom/Giotto_delaunay_celllineage.csv"),row.names = TRUE)
#write.csv(data,file=paste0("./../../../../../MI_heart_paper/results_coloc/Giotto_knn8_MIdata.csv"),row.names = TRUE)
write.csv(data,"./../../../Comparison/202411_results_simfixed_asym//Giotto_delaunay_4ct_cross01.csv", row.names = TRUE)

#df_heatmap = data %>% select(`Endocardial.cells_Mono...Macros.Ccr2.`, Endocardial.cells_Neutrophils)

library(pheatmap)
pheatmap(data)

sessionInfo()

