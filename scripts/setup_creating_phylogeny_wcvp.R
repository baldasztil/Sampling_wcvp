
# Authors: Ludwig Baldaszti - lbaldaszti@rbge.org.uk
# Date: 24/02/2026



# Libraries --------------------------------------------------------------------
library(tidyverse)
library(data.table)
library(sf)
library(ape)
library(U.PhyloMaker)

# Import data ------------------------------------------------------------------


# downloaded from https://github.com/megatrees
megatree <- read.tree("data/phylomaker/plant_megatree.tre")
gen_list <- read.csv("data/phylomaker/plant_genus_list.csv")

# downloaded from https://sftp.kew.org/pub/data-repositories/WCVP/Archive/  version 12
plants_full <- fread("data/wcvp/wcvp_names_full.csv")
plants_full$taxon_name <- str_replace(plants_full$taxon_name, " ", "_")


plants_full_accepted <- fread("data/wcvp_accepted_merged.txt")
plants_full_accepted$taxon_name <- str_replace(plants_full_accepted$taxon_name, " ", "_")


# name matching ----------------------------------------------------------------


search <- plants_full_accepted %>%  
  filter (taxon_name %in% megatree$tip.label) %>% 
  reframe(genus_n = n_distinct(genus))

gen_number <- length(unique(plants_full_accepted$genus))



plants_formated_sample <- plants_full_accepted


sp_list <- plants_formated_sample %>% 
  dplyr::select(species = taxon_name, genus, family) %>% 
  mutate(
         species.relative = NA, 
         genus.relative = NA)  


gen_list_wcvp <- plants_formated_sample %>% 
  filter(taxon_rank == "Genus")


gen_list <- plants_full %>% 
  filter(taxon_rank == "Genus") %>% 
  group_by(genus, family, accepted_plant_name_id) %>% 
  summarise()


tip_labels <- as.data.frame(megatree$tip.label) 

tip_labels <- tip_labels %>% 
  rename(taxon_name = `megatree$tip.label`) %>% 
  mutate(taxon_name2 = taxon_name, 
         test_original = taxon_name == taxon_name2)


name_list_tips <- plants_full %>% 
  filter(taxon_name %in% tip_labels$taxon_name) %>%  
  filter(!taxon_status == "Accepted") %>% 
  dplyr::select(taxon_name, accepted_plant_name_id) %>% 
  filter(!is.na(accepted_plant_name_id))


accepted_names_list <- plants_full %>% 
  filter(accepted_plant_name_id %in% name_list_tips$accepted_plant_name_id) %>% 
  filter(taxon_status == "Accepted") %>% 
  mutate(accepted_name = taxon_name) %>% 
  dplyr::select(accepted_name, accepted_plant_name_id)


names_reformed <- name_list_tips %>% 
  left_join(accepted_names_list, by = "accepted_plant_name_id") %>% 
  mutate(accepted_name = word(accepted_name, 1)) %>% 
  dplyr::select(accepted_name, taxon_name)

names_good <- plants_full_accepted %>% 
  filter(taxon_name %in% tip_labels$taxon_name) %>%  
  filter(taxon_status == "Accepted")  %>% 
  mutate(accepted_name = taxon_name) %>% 
  dplyr::select(accepted_name, taxon_name) %>% 
  rbind(names_reformed) %>% 
  filter(!duplicated(.)) %>% 
  mutate(test = taxon_name == accepted_name)


names_comb <- names_good %>% 
  right_join(tip_labels, by = "taxon_name") %>% 
  filter(!duplicated(taxon_name)) %>% 
  mutate(accepted_name = ifelse(is.na(accepted_name), taxon_name, accepted_name)) %>% 
  dplyr::select(accepted_name, taxon_name, test)


# adjust tip lables
tip_labels_new <- tip_labels %>% 
  left_join(names_comb, by = "taxon_name") %>% 
  dplyr::select(accepted_name)

megatree$tip.label <- tip_labels_new$accepted_name


genus_megatree <- unique(str_split(megatree$tip.label, "_", simplify = T)[,1])

gen_list_wcvp_accepted <- plants_formated_sample %>% 
  group_by(genus, family) %>% 
  summarise()


# tree building ----------------------------------------------------------------

result <- phylo.maker(sp_list, megatree, gen_list_wcvp_accepted, nodes.type = 2, scenario = 3, output.tree	= T)
#save.image("output_tree_wcvp")
output_tree <- result$phylo

check <- result$sp.list

unique(result$sp.list$output.note)


check2 <- check %>% 
  filter(output.note == "no insertion")

check2$species <- gsub(" ", "_", check2$species)

library(phytools)
final <- add.random(result$phylo, tips = check2$species)


write.tree(final$phylo, "data/output_tree.tre")
write.table(final$sp.list, "data/output_splist.txt")

