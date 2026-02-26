
# Authors: Ludwig Baldaszti - lbaldaszti@rbge.org.uk
# Date: 24/02/2026


# Libraries --------------------------------------------------------------------
library(tidyr)
library(dplyr)
library(purrr)
library(stringr)
library(data.table)
library(sf)
library(rWCVP)
library(plotrix)
library(parquetize)
library(arrow)


# Defining functions -----------------------------------------------------------

error <- function(x) {
  qnorm(0.975)*(sd(x)/sqrt(length(x)))
} 


options(digits = 5)  
# Import data ------------------------------------------------------------------

path <- "E:/samples/all_diversity_parquet_gbif/phy/"

sample_data_con <- open_dataset(
  sources = path, 
  format = "parquet"
)



# Reframe ----------------------------------------------------------------------



sample_data <- sample_data_con %>% 
  dplyr::select(continent, corsp_mean, corsp_median, 
                corsp_sd, corsp_iqr, corsp_cov, n_prop, index, source, sp, datasetkey) %>% 
  filter(continent == "GLOBAL") %>%  
  collect() 


setDT(sample_data)

index_source <- c(unique(sample_data$index), unique(sample_data$source))





samples_cumulative_rel <- rollup(sample_data, 
                                 list(
                                                "mean_corgw" = mean(corsp_median, na.rm =T), 
                                                "median_corgw" = median(corsp_median, na.rm =T), 
                                                
                                                
                                                "sd_corgw_mean" = mean(corsp_sd, na.rm =T), 
                                                "se_corgw" = std.error(corsp_median, na.rm =T),
                                                "cov_corgw" = std.error(corsp_cov, na.rm =T),
                                                
                                                "lower_ci_corgw" = mean(corsp_median, na.rm =T)  - error(corsp_mean), 
                                                "upper_ci_corgw" = mean(corsp_median, na.rm =T)  + error(corsp_mean),
                                              
                                                
                                                "mean_n_prop" = mean(n_prop, na.rm =T),
                                                "mean_iqr" = mean(corsp_iqr, na.rm =T), 
                                                "samples" = length(corsp_mean)
                                                
                                                ),
                                 by = c("continent", "sp")) |> 
  drop_na() %>% 
  mutate(index = index_source[1], 
         source = index_source[2])


rm(sample_data)
gc()

write_parquet(samples_cumulative_rel, paste0("output/samples_aggregated/Global_diversity_parquetfullsamples_", index_source[1], "_", index_source[2], ".parquet"))
#Output analysis----------------------------------------------------------------



