
# Authors: Ludwig Baldaszti - lbaldaszti@rbge.org.uk
# Date: 24/02/2026


# Libraries --------------------------------------------------------------------
library(tidyverse)
library(data.table)
library(parquetize)
library(arrow)
library(duckdb)

# Defining functions -----------------------------------------------------------



options(digits = 5)  
# Import data ------------------------------------------------------------------

path <- "D:/samples/all_diversity_parquet_gbif/phy/" # change path here 

sample_data_con <- open_dataset(
  sources = path, 
  format = "parquet"
)


samples_all <- sample_data_con %>% 
  map_batches(function(batch) {
    
    batch %>%
      as.data.frame() %>%
      filter(corsp_median >= 0.95) %>%
      group_by(continent, datasetkey) %>% 
      filter(sp == min(sp)) %>% 
      as_record_batch()
    
    }) %>%
  
  collect()


threshold <- samples_all %>% 
  group_by(continent, datasetkey) %>% 
  slice_min(sp,  n = 1)


threshold_n <- threshold %>%
  group_by(datasetkey) %>%
  summarise(n = n())



write.table(threshold,paste0("Average_sp_threshold_median_", unique(threshold$index),"_", unique(threshold$source), ".txt"))
