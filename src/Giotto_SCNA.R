#Giotto usage for simulated data as input

library(Giotto)
library(reticulate)
library(readr)
library(tidyverse)
library(stats)

# Giotto and python configurations
# Create Giotto path
my_instructions = createGiottoInstructions(python_path = '/Users/chiaraschiller/Library/r-miniconda-arm64/bin/python3')
Sys.setenv(RETICULATE_PYTHON = "/Users/chiaraschiller/Library/r-miniconda-arm64/bin/python3") 
#py_config()

#installGiottoEnvironment()

# load simulated data

files = list.files("./../../../../data/", pattern = ".csv")
data_path = "./../../../../data/"
data_list = list()

for (i in files){
  print(i)
data_list[[i]] <- readr::read_csv(paste0(data_path,i))
}

giotto_list = list()
x=0
# create dataframe with random values to pretend to have expression data to generate a giotto object
for (i in data_list){
  x = x+1
  df = data.frame(cbind(rep(1,nrow(i)), rep(1,nrow(i))))
  #df <- data.frame(matrix(runif(nrow(data_list[[i]]) * 2), nrow = nrow(data_list[[i]])))
  df$cell_ID = seq(1, nrow(df), by = 1)
  colnames(df) = c("M1", "M2", "cell_ID")

  metadata = cbind(as.character(i$ct), seq(1, nrow(df), by = 1))
  colnames(metadata) = c("ct", "cell_ID")
  # colbind the two datasets
  #data = cbind(df, as.data.frame(data))
  data = createGiottoObject(raw_exprs = as.data.frame(t(df)),
                          spatial_locs = as.data.frame(i[, !(names(i) %in% c("ct"))]),
                          #instructions = my_instructions,
                          cell_metadata = metadata)

  data <- normalizeGiotto(gobject = data)
  # create network (required for binSpect methods)
  giotto_list[[files[x]]] = createSpatialNetwork(gobject = data, minimum_k = 2)
}

# calculate cell proximities
# PI value is the proximity index of the feature in the cluster
cell_PI=list()
x=0

for (i in giotto_list){
  x=x+1
  cell_proximities = cellProximityEnrichment(gobject = i,
                        cluster_column = 'ct',
                        spatial_network_name = 'Delaunay_network',
                        adjust_method = 'fdr',
                        number_of_simulations = 1000)

  cell_PI[[files[x]]] = separate(cell_proximities[[2]], unified_int, into = c("label1", "label2"), sep = "--")
  #names(cell_PI)[i]=gsub("\\..*", "", files[i])
  write.csv(cell_PI[[files[x]]],file=paste0("./../../output/",gsub("\\..*", "", files[x]) ,"_proximities.csv"),row.names = TRUE)
}


