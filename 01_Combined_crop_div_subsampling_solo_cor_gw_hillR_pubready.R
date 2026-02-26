
# Authors: Ludwig Baldaszti - lbaldaszti@rbge.org.uk
# Date: 24/02/2026



# Libraries --------------------------------------------------------------------
library(tidyverse)
library(data.table)
library(sf)
library(spmodel)
library(ape)
library(phyloregion)
library(GWmodel)
library(arrow)
library(vegan)
library(hillR)
library(moments)





# Defining functions -----------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
diversity_dimension <- args[1]
print(diversity_dimension)



cov <- function(x) {
  sd(x) / mean(x)     
}

cor.gw <- function(rich_rel_shp) {
  bw <- try(bw.gwr(diversity ~ diversity_sample, data=rich_rel_shp,
                   adaptive = T, 
                   longlat = T))
  
  if (!is.numeric(bw)) {
    bw <- bw_global
  }
  
  stats <- try(as.data.frame(gwss(rich_rel_shp, rich_rel_shp, 
                                  vars = c("diversity", "diversity_sample"),
                                  adaptive = T, bw = bw, longlat = T)$SDF))
  
  if (class(stats)[1]=="try-error") {
    error_data <- output_model %>% 
      replace(!is.numeric(.) ==T, 0) %>% 
      mutate(id = output_model$id, 
             sp = nrow(species_sample))
    return(error_data)
  }
  
  else {
    
    stats_total <-  cbind(stats, rich_rel_shp) %>% 
      dplyr::select(cor.sp_gwr =Spearman_rho_diversity.diversity_sample, LEVEL3_COD, 
                    continent = LEVEL1_NAM) %>% 
      replace(is.na(.), 0) %>% 
      setDT()
    

    
    funs <- list(mean, median, sd, skewness, kurtosis, min, max, cov, IQR)
    names(funs) <- c("mean", "median",  "sd", "skewness", "kurtosis", "min", "max", "cov", "iqr")
    
    global <- stats_total %>%  
      summarise(across(
        .cols = where(is.numeric), 
        .fns = funs,  
        .names = "corsp_{fn}"
      )) %>% 
      mutate(continent = "GLOBAL") 
    
    continents <- stats_total %>% 
      group_by(continent) %>% 
      summarise(across(
        .cols = where(is.numeric), 
        .fns = funs, 
        .names = "corsp_{fn}"
      )) 
    
    output_file <- rbind(continents, global) 
    return(output_file)
  }
}

subsampling.plants <- function(spec_n, dimension) {
  
  if (nrow(plantlist_names) > spec_n) {
    species_sample <- sample_n(plantlist_names, spec_n)
    # too include last sample
  } else {
    species_sample <- plantlist_names
  }
  
  dist <- plantlist_dist_phylo_growth %>% 
    filter(plant_name_id %in% species_sample$plant_name_id)
  
  
  sp_num <- dist %>%
    left_join(continent_names, by = "LEVEL3_COD") %>% 
    group_by(LEVEL1_NAM) %>%  
    summarise(n = n_distinct(plant_name_id))
  
  if(dimension == "diversity_phy") {
    subset_matrix <- dist.mat[ ,colnames(dist.mat) %in% dist$label]
    prune_tree <- try(match_phylo_comm(output_tree,subset_matrix))
    
    if(class(prune_tree) == "try-error") {
      error_data <- output_model %>% 
        replace(!is.numeric(.) ==T, 0) %>% 
        mutate(
          n = nrow(species_sample), 
          index = diversity_dimension, 
          n_prop = round(n / sp_num_abs$n, 4), 
          source = "wcvp")
      return(error_data)
    }  
    
    # doesnt work right now with old tree
    rich_rel <- as.data.frame(PD(subset_matrix, prune_tree$phy)) %>% 
      mutate(LEVEL3_COD = rownames(.)) %>% 
      rename(diversity_sample = `PD(subset_matrix, prune_tree$phy)`)
  }
  
  if(dimension == "diversity_fun") {
    table_data_sample <- table(dist$LEVEL3_COD, dist$growth_form)
    rich_rel <- as.data.frame(hill_taxa(table_data_sample, q = 1 )) 
    rich_rel$C <- rownames(rich_rel)
    names(rich_rel) <- c("diversity_sample", "LEVEL3_COD")
    
  }
  
  
  if(dimension == "diversity_tax") {
    table_data_sample <- table(dist$LEVEL3_COD, dist$family)
    rich_rel <- as.data.frame(hill_taxa(table_data_sample, q = 1 )) 
    rich_rel$C <- rownames(rich_rel)
    names(rich_rel) <- c("diversity_sample", "LEVEL3_COD")
    
  }
  
  if(dimension == "diversity_sp") {
    rich_rel <- dist %>% 
      group_by(LEVEL3_COD) %>% 
      summarise(diversity_sample = length(unique(plant_name_id)))
  }

  rich_rel_shp <- rich_rel  %>% 
    right_join(rich_overall_bru_mid_used, by = c("LEVEL3_COD")) %>%
    
    #right_join(rich_overall_bru, by = c("LEVEL3_COD")) %>%
    #left_join(midpoints_red, by = c("LEVEL3_COD")) %>% 
    replace(is.na(.), 0) %>% 
    mutate(sp = nrow(species_sample),
           diversity_sample = diversity_sample,
           diversity = diversity) %>% 
    st_sf()  %>% 
    as("Spatial")
  #filter(!LEVEL1_NAM == "ANTARCTICA")
  

  correlation <- cor.gw(rich_rel_shp = rich_rel_shp)
  
  model_comp_out <- correlation %>% 
    left_join(sp_num, by = c("continent" = "LEVEL1_NAM")) %>% 
    mutate(n = case_when(continent == "GLOBAL" ~ nrow(species_sample), 
                         is.na(n) ~ 0, 
                         TRUE ~n ), 
           n_prop = round(n / sp_num_abs$n, 4), 
           index = diversity_dimension, 
           source = "wcvp")
  
  #cumulative pattern
  if (nrow(species_sample) %in% seq(2,nrow(plantlist_names),50000)){
    tf1 <- paste0("data/",dimension, "/checkpoints/checkpoint_",dimension,"_", nrow(species_sample),".parquet")
    write_parquet(model_comp_out, tf1) 
  }
  
  return(model_comp_out)
}


# Import data ------------------------------------------------------------------
dist_native <- dist_native <- fread("data/dist_native.txt") 
plants_full_raw <- fread("data/wcvp_accepted_merged.txt")
output_tree <- read.tree("data/output_tree_20012024.tre")

tdwg_3 <- st_read(dsn ="data/wgsrpd-master/level3") %>% 
  filter(!LEVEL3_COD == "BOU")

midpoints_red <- st_read(dsn ="data/wgsrpd-master/level3_midpoints")

continent_names <- tdwg_3 %>% 
  dplyr::select(LEVEL1_NAM, LEVEL3_COD) %>% 
  st_drop_geometry()

growth <- fread("data/growth_forms_nnet.txt") %>% 
  dplyr::select(growth_form, plant_name_id)

plants_full <- plants_full_raw %>% 
  left_join(growth, by = "plant_name_id")


continent_names_vec <- unique(continent_names$LEVEL1_NAM)

# manipulate data --------------------------------------------------------------

plantlist_names <- plants_full %>%  
  dplyr::select(plant_name_id, taxon_rank, family, taxon_name, growth_form)

plantlist_dist_phylo_growth <- dist_native %>% 
  left_join(plantlist_names, by = "plant_name_id") %>% 
  dplyr::select(plant_name_id, taxon_name, LEVEL3_COD = area_code_l3, growth_form, family) %>% 
  mutate(label = gsub(" ", "_", taxon_name))

dist.mat <- long2sparse(plantlist_dist_phylo_growth, grids = "LEVEL3_COD", species = "label")


table_data <- table(plantlist_dist_phylo_growth$LEVEL3_COD, plantlist_dist_phylo_growth$growth_form)
diversity_patterns_fun <- as.data.frame(hill_taxa(table_data, q = 1 )) 
diversity_patterns_fun$C <- rownames(diversity_patterns_fun)
names(diversity_patterns_fun) <- c("diversity_fun", "LEVEL3_COD")


table_data <- table(plantlist_dist_phylo_growth$LEVEL3_COD, plantlist_dist_phylo_growth$family)
diversity_patterns_tax <- as.data.frame(hill_taxa(table_data, q = 1 )) 
diversity_patterns_tax$C <- rownames(diversity_patterns_tax)
names(diversity_patterns_tax) <- c("diversity_tax", "LEVEL3_COD")



prune_tree <- match_phylo_comm(output_tree,dist.mat)
diversity_patterns_phy <- as.data.frame(PD(dist.mat, prune_tree$phy), 
                                       col.names =c("diversity")) %>% 
  mutate(LEVEL3_COD = rownames(.)) %>% 
  rename(diversity_phy = `PD(dist.mat, prune_tree$phy)`) 


diversity_patterns_rich <- plantlist_dist_phylo_growth %>% 
  group_by(LEVEL3_COD) %>% 
  summarise(diversity_sp = length(unique(plant_name_id))) 

diversity_patterns <- diversity_patterns_rich %>% 
  left_join(diversity_patterns_tax, by = "LEVEL3_COD") %>% 
  left_join(diversity_patterns_phy, by = "LEVEL3_COD") %>% 
  left_join(diversity_patterns_fun, by = "LEVEL3_COD") %>% 
  left_join(continent_names, by = "LEVEL3_COD") 

rich_patterns_long <- diversity_patterns %>% 
  pivot_longer(cols = starts_with("diversity"),
               names_to = "dim",
               values_to = "diversity") %>% 
  mutate(diversity_scaled = scale(diversity)[,1])


rich_overall_bru <- rich_patterns_long


rich_overall_bru_mid <- rich_patterns_long %>% 
  left_join(midpoints_red, by = "LEVEL3_COD")  %>% 
  mutate(diversity2 = diversity) %>% 
  st_as_sf()

rich_overall_bru_shp <- rich_patterns_long %>% 
  left_join(tdwg_3, by = "LEVEL3_COD") %>% 
  st_as_sf()


rich_overall_bru_mid_used <- rich_overall_bru_mid %>% 
  filter(dim == diversity_dimension) 

rich_overall_bru_mid_shp_used <- rich_overall_bru_mid_used %>% 
  as("Spatial")

bw_global <- bw.gwr(diversity ~ diversity2, 
                   data = rich_overall_bru_mid_shp_used, 
                   approach = "CV",
                   adaptive = T)

#rich_overall_bru_mid_test <- split(rich_overall_bru_mid_shp,rich_overall_bru_mid_shp$dim)

sp_num_abs <- plantlist_dist_phylo_growth %>%
  left_join(continent_names, by = "LEVEL3_COD") %>% 
  group_by(LEVEL1_NAM) %>%  
  summarise(n = n_distinct(plant_name_id))

sp_num_abs[10,1] <- c("GLOBAL")

sp_num_abs[10,2] <- nrow(plantlist_names)


setDT(plantlist_names)
setDT(dist_native)
setDT(rich_overall_bru)
setDT(midpoints_red)

output_model <- subsampling.plants(349113, dimension = diversity_dimension)


# compare models ---------------------------------------------------------------

sample <- sample(1:10000000, 1)
#set.seed(sample)

start <- Sys.time()
subsampling.plants(349113, dimension = diversity_dimension)
Sys.time() - start

result_list <- lapply(seq(1,nrow(plantlist_names),1), 
                      subsampling.plants,  dimension = diversity_dimension)

results <- rbindlist(result_list)

file <- paste0("data/",diversity_dimension,"/Sample_hillr_",sample, "_", 
               diversity_dimension,"_wcvp.parquet")
write_parquet(results, file) 

