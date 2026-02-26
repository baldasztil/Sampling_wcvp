
# Authors: Ludwig Baldaszti - lbaldaszti@rbge.org.uk
# Date: 24/02/2026


# Libraries --------------------------------------------------------------------
library(tidyverse)
library(data.table)
library(ggplot2)
library(sf)
library(GWmodel)
library(feather)
library(datawizard)
library(plotrix)

# Defining functions -----------------------------------------------------------
quiet <- function (x, print_cat = TRUE, message = TRUE, warning = TRUE) 
{
  stopifnot(is.logical(print_cat) && length(print_cat) == 1)
  stopifnot(is.logical(message) && length(message) == 1)
  stopifnot(is.logical(warning) && length(warning) == 1)
  if (print_cat) 
    sink(tempfile(), type = "out")
  on.exit(if (print_cat) sink())
  if (warning && message) 
    invisible(force(suppressMessages(suppressWarnings(x))))
  else if (warning && !message) 
    invisible(suppressWarnings(force(x)))
  else if (!warning && message) 
    invisible(suppressMessages(force(x)))
  else invisible(force(x))
}

cov <- function(x) {
  sd(x) / mean(x)     
}
cor.gw <- function(rich_rel_shp) {
  
  continent <- unique(rich_rel_shp$continent_name)
  if (length(continent) > 1) {
    continent <- "GLOBAL"
  }
    
  
  bw <- try(quiet(bw.gwr(richness ~ richness_sample, data=rich_rel_shp,
                         adaptive = T, 
                         longlat = T)))
  
  if (!is.numeric(bw)) {
    bw <- nrow(rich_rel_shp)
  }
  
  stats <- try(as.data.frame(gwss(rich_rel_shp, rich_rel_shp, 
                                  vars = c("richness", "richness_sample"),
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
      dplyr::select(cor.sp_gwr =Spearman_rho_richness.richness_sample, LEVEL3_COD, 
                    continent = LEVEL1_NAM) %>% 
      replace(is.na(.), 0) %>% 
      setDT()
    
 
    
    funs <- list(mean, median, sd, min, max, cov, IQR)
    names(funs) <- c("mean", "median",  "sd", "min", "max", "cov", "iqr")
    
    
    
    output_file <- stats_total %>%  
      reframe(across(
        .cols = where(is.numeric), 
        .fns = funs,  
        .names = "corsp_{fn}"
      )) %>% 
      mutate(continent = continent)
    

    
    return(output_file)
  }
}

subsampling.plants.tax <- function(tax, source_names, source_dist, fam = T) {
  
  if(fam == T) {
    number <- which(families %in% fam)
    
    print(paste0("This is the ", tax, " family"," [", number, "/", length(families), "]"))
    
    species_sample <- source_names %>% 
      filter(family %in% tax)
    
  } else {
    number <- which(genera %in% tax)
    
    print(paste0("This is ", tax, " [", number, "/", length(genera), "]"))
    
    species_sample <- source_names %>% 
      filter(genus %in% tax)
  }


  
  sp_origin <- dist_native %>% 
    filter(plant_name_id %in% species_sample$plant_name_id) %>% 
    setDT() %>% 
    cube(., j= length(unique(plant_name_id)), by = "continent") %>% 
    mutate(continent = ifelse(is.na(continent), "GLOBAL", continent)) %>% 
    dplyr::select(continent, sp_origin = V1)

  
  dist <- source_dist %>% 
    filter(plant_name_id %in% species_sample$plant_name_id) %>% 
    left_join(continent_names, by = "LEVEL3_COD") %>% 
    rename(continent = LEVEL1_NAM)
  
  sp_source <- dist %>% 
    filter(plant_name_id %in% species_sample$plant_name_id) %>% 
    setDT() %>% 
    cube(., j= length(unique(plant_name_id)), by = "continent") %>% 
    mutate(continent = ifelse(is.na(continent), "GLOBAL", continent)) %>% 
    dplyr::select(continent, sp_source = V1)
  
  # richness patterns across brus
  rich_rel <- dist %>% 
    group_by(continent) %>% 
    mutate(sp = n_distinct(plant_name_id)) %>% 
    ungroup() %>%  
    group_by(LEVEL3_COD) %>% 
    summarise(richness_sample = n_distinct(plant_name_id), 
              sp = unique(sp)) %>% 
    right_join(midpoints_red, by = c("LEVEL3_COD")) %>% 
    right_join(rich_overall_bru, by = "LEVEL3_COD") %>%
    right_join(continent_names, by = "LEVEL3_COD") %>%
    mutate(continent_name = LEVEL1_NAM) %>% 
    replace(is.na(.), 0) %>% 
    mutate(
           richness_sample = richness_sample,
           richness = richness, 
           ) %>% 
    st_sf() 
  
  
  rich_rel_shp <-  rich_rel %>%  
  as("Spatial")

  correlation <- cor.gw(rich_rel_shp = rich_rel_shp)
  
  rich_rel_shp_continent <- lapply(split(rich_rel,rich_rel$LEVEL1_NAM), FUN = function(x)
    as(x,"Spatial"))
  
  
  model_comp_out <- lapply(rich_rel_shp_continent, cor.gw) %>% 
    rbindlist() %>% 
    rbind(correlation) %>% 
    left_join(sp_origin, by = c("continent")) %>% 
    left_join(sp_source, by = c("continent")) %>% 
    mutate(
      sp_miss = sp_origin - sp_source, 
      tax =   paste(tax, collapse = " x "))
  

  
  
  #cumulative pattern 
  
  return(model_comp_out)
}

# Import data ------------------------------------------------------------------
dist_native <- dist_native <- fread("data/dist_native.txt") 


plants_full <- fread("data/wcvp_accepted_merged.txt")

tdwg_3 <- st_read(dsn ="data/wgsrpd-master/level3") %>% 
  filter(!LEVEL3_COD == "BOU")

midpoints_red <- st_read(dsn ="data/wgsrpd-master/level3_midpoints")

continent_names <- tdwg_3 %>% 
  dplyr::select(LEVEL1_NAM, LEVEL3_COD) %>% 
  st_drop_geometry()

continent_names <- tdwg_3 %>% 
  dplyr::select(LEVEL1_NAM, LEVEL3_COD) %>% 
  st_drop_geometry()

continent_names_vec <- unique(continent_names$LEVEL1_NAM)


# manipulate data --------------------------------------------------------------

plantlist_dist <- dist_native %>% 
  left_join(plants_full, by = c("plant_name_id")) %>% 
  dplyr::select(plant_name_id, LEVEL3_COD = area_code_l3)


plantlist_names <- plants_full %>%  
  dplyr::select(plant_name_id, taxon_rank, family, taxon_name)

nrow(plantlist_names)

rich_overall_bru <- plantlist_dist %>% 
  group_by(LEVEL3_COD) %>% 
  summarise(richness = length(unique(plant_name_id))) 


rich_overall_bru_mid <- rich_overall_bru %>% 
  mutate(richness2 = richness) %>% 
  left_join(midpoints_red, by = "LEVEL3_COD") %>% 
  st_as_sf()

rich_overall_bru_shp <- rich_overall_bru %>% 
  left_join(tdwg_3, by = "LEVEL3_COD") %>%
  mutate(richness2 = richness) %>% 
  st_as_sf() 

rich_overall_bru_mid_shp <- rich_overall_bru_mid %>% 
  as("Spatial") 

bw_global <- 368



setDT(plantlist_names)
setDT(plantlist_dist)
#setDT(gbif_dist)
setDT(rich_overall_bru)
setDT(midpoints_red)

families <- unique(plants_full$family)
names(families) <- "fam"

genera <- unique(plants_full$genus)
names(genera) <- "gen"
# run function ==---------------------------------------------------------------

test <- "Poaceae" 


result_list_fam <- lapply(families, 
                      subsampling.plants.tax, source_dist = plantlist_dist, 
                      source_names = plants_full, fam = T)

results <- rbindlist(result_list_fam) %>% 
  mutate(source = "wcvp")

res_file <- paste0("output/Corgw_rich_family.txt")
fwrite(results, file=res_file) 


result_list_gen <- lapply(genera, 
                          subsampling.plants.tax, source_dist = plantlist_dist, 
                          source_names = plants_full, fam = F)

results <- rbindlist(result_list_gen) %>% 
  mutate(source = "wcvp")

res_file <- paste0("output/Corgw_rich_genus.txt")
fwrite(results, file=res_file) 


