#Giotto usage for multiple samples
# script created by Margot Chazotte for analyzing multiple samples with Giotto

library(Giotto)
library(reticulate)
library(readr)
library(tidyverse)
library(ggplot2)
library(pheatmap)

getwd()
# Giotto and python configurations
# Create Giotto path
my_instructions = createGiottoInstructions(python_path = '/Users/margotchazotte/Library/r-miniconda/bin/python3')
Sys.setenv(RETICULATE_PYTHON = "/Users/margotchazotte/Library/r-miniconda/bin/python3") 
#py_config()


#dataframe of markers (M1 to Mn), celltypes (celltype), Cell ID (CellID), X and Y Positions (X , Y), Sample (ID)
data <- read_csv("Path/to/CSV.csv")

#use sample ID to identify different samples
sample <- setNames(as.character(unique(data$ID)), as.character(unique(data$ID)))
list <- lapply(sample , function(i) filter(data, ID == i))

#create giotto object and fill with necessary computations
giotto_list <- lapply(sample , function(i) createGiottoObject(raw_exprs = as.data.frame(t(select(list[[i]], M1:Mn))),
                                                             cell_metadata = select(list[[i]], ct = celltype, CellID),
                                                             spatial_locs = select(list[[i]], x = X, y = Y)))
giotto_list <- lapply(sample, function(i) normalizeGiotto(gobject = giotto_list[[i]]))
giotto_list <- lapply(sample, function(i) addStatistics(gobject = giotto_list[[i]]))
giotto_list <- lapply(sample, function(i) createSpatialNetwork(gobject = giotto_list[[i]], minimum_k = 2))

#calculate cell proximities
cell_proximities <- lapply(sample, function(i) cellProximityEnrichment(gobject = giotto_list[[i]],
                                                                       cluster_column = 'ct',
                                                                       spatial_network_name = 'Delaunay_network',
                                                                       adjust_method = 'fdr',
                                                                       number_of_simulations = 1000))

cell_PI <- lapply(sample, function(i) separate(cell_proximities[[i]][[2]], unified_int, into = c("label1", "label2"), sep = "--"))

#save cell proximity dataframes in created folder
dir.create("Path/To/Outputfolder")
for (i in 1:length(cell_PI)) {
  write.csv(cell_PI[[i]],file=paste0("Path/To/Outputfolder/", cores[i] ,"_proximities.csv"),row.names = TRUE)
}

#heatmap of all interactions for every sample
cell_PI <- lapply(sample, function(i) cbind(cell_PI[[i]], ID = i))
out <- rbindlist(cell_PI)
df = as.data.frame(out) %>% select(label1, label2, PI_value, ID)
df$key = paste(df$label1, df$label2, sep = "_")
df = df %>% select(-c(label1, label2))
df = spread(df, key = key, value = PI_value)

rownames(df) = df$ID
pheatmap(df[,-1], treeheight_row = 1, treeheight_col = 1)


#enrichment and interaction number plot for all samples
sampleID <- setNames(names(cell_proximities), names(cell_proximities))
enrichment <- lapply(sampleID, function(i) bind_cols(cell_proximities[[i]]$enrichm_res, sample = i))
enrichment <- rbindlist( enrichment ,use.names = TRUE, fill = TRUE)
ordered <- enrichment%>%
  group_by(unified_int)%>%
  summarise( avg = mean(enrichm))
enrichment$unified_int <- factor(enrichment$unified_int, levels=as.data.frame(ordered[order(ordered$avg), 1])[,1])

cowplot::plot_grid( 
  enrichment %>%
    ggplot( aes(unified_int, enrichm, fill = type_int)) +
    geom_boxplot() +
    coord_flip() +
    theme(legend.position = "none") +
    ylab("enrichment/depletion")+
    xlab("interaction"),
  enrichment %>%
    ggplot(aes(unified_int, log10(original), fill = type_int)) +
    geom_boxplot() +
    coord_flip() +
    xlab("") +
    ylab("log 10 # of interactions") + 
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
    scale_fill_discrete(name = "type"),
  rel_widths = c(2,1))

