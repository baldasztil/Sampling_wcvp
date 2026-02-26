# Subsampling Global Plant Biodiversity 13/02/2022
# data available @ ??

# Authors: Ludwig Baldaszti - lbaldaszti@rbge.org.uk
# Date: 13/02/2022

# Intro ------------------------------------------------------------------------


# Libraries --------------------------------------------------------------------
library(tidyverse)
library(data.table)
library(sf)
library(spmodel)
library(ape)
library(phyloregion)
library(GWmodel)
library(hillR)
library(feather)
library(gmodels)
library(datawizard)
library(moments)
library(paletteer)
library(nationalparkcolors)
library(caret)
library(rstatix)
library(patchwork)
library(ggpubr)
library(GGally)
library("ggnewscale")
library(egg)
library(plotrix)





# Defining functions -----------------------------------------------------------
tag_facet2 <- function(p, open = "(", close = ")", tag_pool = letters, x = -Inf, y = Inf, 
                       hjust = -0.5, vjust = 1.5, fontface = 2, family = "", ...) {
  
  gb <- ggplot_build(p)
  lay <- gb$layout$layout
  tags <- cbind(lay, label = paste0(open, tag_pool[lay$PANEL], close), x = x, y = y)
  p + geom_text(data = tags, aes_string(x = "x", y = "y", label = "label"), ..., hjust = hjust, 
                vjust = vjust, fontface = fontface, family = family, inherit.aes = FALSE) +
    labs(x = NULL, y = NULL)
}
cov <- function(x) {
  sd(x) / mean(x)     
}
min_max_norm <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}
mad_iqr <-  function(x) {
  (x - median(x)) / IQR(x)
}



equal_breaks <- function(n = 3, s = 0.05, ...){
  function(x){
    # rescaling
    #d <- s * diff(range(x)) / (1+2*s)
    plyr::round_any(seq(min(x), max(x), length=n), 5)
  }
}


symmetric_limits <- function (x) 
{
  max <- plyr::round_any(max(abs(x)) * 1.1, 100,  f = ceiling)
  if (min(abs(x)) > max)
    max <- plyr::round_any(min(abs(x)) * 1.1, 100, f = ceiling)
  
  if (max - 370 > 0 )
    max <- 350
  
  c(-max, max)
}

mapping.plants.random <- function(spec_n) {
  
  
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
  
  subset_matrix <- dist.mat[ ,colnames(dist.mat) %in% dist$label]
  prune_tree <- try(match_phylo_comm(output_tree,subset_matrix))
  
  # doesnt work right now with old tree
  rich_rel_pd <- as.data.frame(PD(subset_matrix, prune_tree$phy)) %>%   
    mutate(LEVEL3_COD = rownames(.)) %>%  
    rename(diversity_phy = `PD(subset_matrix, prune_tree$phy)`)
  
  table_data_sample <- table(dist$LEVEL3_COD, dist$growth_form)
  rich_rel_fun <- as.data.frame(hill_taxa(table_data_sample, q = 1 )) 
  rich_rel_fun$C <- rownames(rich_rel_fun)
  names(rich_rel_fun) <- c("diversity_fun", "LEVEL3_COD")
  
  
  table_data_sample <- table(dist$LEVEL3_COD, dist$family)
  rich_rel_tax <- as.data.frame(hill_taxa(table_data_sample, q = 1 )) 
  rich_rel_tax$C <- rownames(rich_rel_tax)
  names(rich_rel_tax) <- c("diversity_tax", "LEVEL3_COD")
  
  
  rich_rel <- dist %>%   
    group_by(LEVEL3_COD) %>%   
    summarise(diversity_sp = length(unique(plant_name_id)))
  
  richness_rel <- rich_rel %>%   
    left_join(rich_rel_pd, by = "LEVEL3_COD") %>%   
    left_join(rich_rel_tax, by = "LEVEL3_COD") %>%   
    left_join(rich_rel_fun, by = "LEVEL3_COD")
  
  rich_rel_long <- richness_rel %>%   
    pivot_longer(cols = starts_with("diversity"),
                 names_to = "dim",
                 values_to = "richness_sample") %>%   
    mutate(richness_sample_scaled = scale(richness_sample)[,1]) %>%  
    right_join(rich_overall_bru, by = c("LEVEL3_COD", "dim")) %>%    
    replace(is.na(.), 0) 
  
  return(rich_rel_long)
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


rich_overall_bru <- rich_patterns_long %>%   
  mutate(dataset = "wcvp")


rich_overall_bru_mid <- rich_patterns_long %>%   
  left_join(midpoints_red, by = "LEVEL3_COD") %>%  
  mutate(diversity2 = diversity) %>%   
  st_as_sf()

rich_overall_bru_shp <- rich_patterns_long %>%   
  left_join(tdwg_3, by = "LEVEL3_COD") %>%   
  st_as_sf()

bw_global <- 368

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



# extract samples --------------------------------------------------------------
set.seed(11)

spec_n <- 1048

wcvp_mapping_random_1 <- lapply(rep(spec_n,100), 
                                mapping.plants.random) %>%   
  rbindlist(idcol = "sample") %>%  
  mutate(dataset = paste0("sample_size_", spec_n))




spec_n <- 67772

wcvp_mapping_random_2 <- lapply(rep(spec_n,100), 
                                mapping.plants.random) %>%   
  rbindlist(idcol = "sample") %>%  
  mutate(dataset = paste0("sample_size_", spec_n))


sample <- rbind(wcvp_mapping_random_1, wcvp_mapping_random_2)

# assemble diversity stats -----------------------------------------------------

overall_joined <- sample %>%   
  group_by(dim) %>%   
  left_join(tdwg_3, by = "LEVEL3_COD") %>%   
  st_drop_geometry() %>%   
  dplyr::select(LEVEL3_COD, dim, diversity = richness_sample, dataset, continent = LEVEL1_NAM.x)

combined_rich <- rich_overall_bru %>%   
  mutate(sample = 1) %>%   
  dplyr::select(LEVEL3_COD, dim, diversity , dataset, continent = LEVEL1_NAM) %>%   
  rbind(overall_joined) 


xx <- sample_stats %>% 
  st_drop_geometry()

xxx <- combined_rich %>% 
  filter(LEVEL3_COD == "MCS") %>% 
  filter(dim == "diversity_sp") %>% 
  group_by(dataset) %>% 
  summarise(mean =  ci(diversity)[1], 
            median = median(diversity), 
            sd = sd(diversity), 
            coef = sd / mean * 100)

mean(xxx$diversity, na.rm = T)

summary_stats  <- combined_rich %>%   
  group_by(LEVEL3_COD, dim, dataset) %>%  
  summarise(mean = ci(diversity)[1], 
            sd = sd(diversity),
            lower_ci =  ci(diversity)[2], 
            upper_ci =  ci(diversity)[3], 
            se =  ci(diversity)[4], 
            coef = sd / mean * 100) %>%   
  ungroup() %>%   
  group_by(dim, dataset) %>%   
  mutate(ranked = 369 - ranktransform(mean), 
         scaled_richness = scale(mean), 
         scaled_diversity = min_max_norm(mean), 
         scaled_mad = mad_iqr(mean)) %>%   
  left_join(tdwg_3, by = "LEVEL3_COD") %>%   
  st_as_sf() %>%   
  mutate(dim = factor(dim, levels = c("diversity_sp", "diversity_phy", "diversity_fun",  "diversity_tax"), 
                      labels = c("Species richness", "Phylogenetic diversity", 
                                 "E. growth form diversity ", "E. family diversity")), 
         dataset = factor(dataset, levels = c("wcvp",  "sample_size_1048",  "sample_size_67772"), 
                          labels = c("Reference: 349,113 spp.", "Sample: 1,048 spp.",  "Sample: 67,772 spp.")))



reference <- summary_stats %>%   
  filter(dataset == "Reference: 349,113 spp.") %>%   
  dplyr::select(rank_reference = ranked, scaled_richness_org = scaled_richness, 
                scaled_diversity_org = scaled_diversity, mean_org = mean, 
                LEVEL3_COD, dim) %>%  
  st_drop_geometry()


sample_stats <- summary_stats %>%   
  filter(!dataset == "Reference: 349,113 spp.") %>%   
  group_by(dataset) %>%    
  right_join(reference, by = c("LEVEL3_COD", "dim")) %>%  
  mutate(rank_diff = (ranked - rank_reference) * -1, 
         scaled_richness_diff = (scaled_richness - scaled_richness_org), 
         scaled_diversity_diff = (scaled_diversity - scaled_diversity_org), 
         log_coef = log10(coef), 
         CV = coef / 100, 
         mad = mad(mean)) %>%   
  dplyr::select(LEVEL3_NAM, LEVEL3_COD, LEVEL1_NAM, dim, ranked, rank_reference, rank_diff, 
                scaled_richness, scaled_richness_org, scaled_richness_diff, 
                scaled_diversity, scaled_diversity_org, scaled_diversity_diff, 
                coef, log_coef, CV, se, mad, dataset) %>%   
  mutate(CV = ifelse(is.na(CV), 0, CV), 
         coef = ifelse(is.na(coef), 0, coef))


sample_stats_all <-  sample_stats %>% 
  st_drop_geometry()


sample_stats_all[] <- lapply(sample_stats_all, unname)

stats <- sample_stats_all %>% 
  group_by(dim, dataset) %>% 
  slice_max(abs(rank_diff), n = 1)

bias <- sample_stats_all %>% 
  st_drop_geometry() %>% 
  group_by(dataset, dim) %>% 
  #summarise(mean_rank = mean(abs(rank_diff), na.rm = T))
  mutate(across(where(is.numeric), abs)) %>% 
  summarise(across(where(is.numeric), mean, na.rm = T)) %>% 
  mutate(across(where(is.numeric), round,3))

write.csv(bias, "bias_table.csv")



coef_plot <- ggplot(sample_stats_all, aes(y = coef, x = ranked, fill = dataset)) +
  geom_point(aes(col = dataset), alpha = 0.3, 
              show.legend = T) +
  scale_color_manual(values = palette) +
  scale_x_continuous( name = "Average rank") +  
  scale_y_continuous( name = "Coefficient of variation (%)") +
  theme_bw(base_size = 14)+
  stat_cor(aes(col = dataset), method = "spearman", cor.coef.name = "rho")


ggsave("scatter_cv_rank_corr.png",coef_plot, dpi = 600,width = 7, height = 5)

# Mapping ----------------------------------------------------------------------


eqearth_crs <- "+proj=eqearth"


summary_stats_eq <- st_transform(summary_stats, crs = eqearth_crs)


border <- st_graticule() %>%   
  st_bbox() %>%   
  st_as_sfc() %>%   
  st_transform(3857) %>%   
  st_segmentize(500000) %>%   
  st_transform(st_crs(eqearth_crs)) %>%   
  st_cast("POLYGON")


diversity_plot <- ggplot(summary_stats_eq) +
  geom_sf(data = border, fill = "azure", color = "black", size = 0.2) +
  geom_sf(aes(fill = scaled_diversity), color = "grey40", size = 0.2) +
  
  scale_fill_distiller(
    palette = "YlOrRd",
    direction = 1,
    name = "Diversity",
    breaks = pretty(summary_stats_eq$scaled_diversity, n = 10)
  ) +
  facet_grid(dim ~ dataset, switch = "y") +
  coord_sf(crs = eqearth_crs, expand = T) +
  theme_minimal(base_size = 10) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    strip.background = element_rect(fill = "black"),
    strip.text = element_text(color = "white", face = "bold", size = 6.5),
    panel.grid.major = element_line(color = "grey85", size = 0.2),
    legend.position = "bottom",
    legend.box = "horizontal",
    axis.text = element_blank(),
    axis.ticks = element_blank()
  ) +
  guides(fill = guide_colorbar(
    barwidth = 15,
    barheight = 0.6,
    title.position = "top"
  ))


diversity_plot2 <- tag_facet2(diversity_plot) 
ggsave("div_sample_maps_scaled_minmax_geomsf_sr.png", diversity_plot2, dpi = 600)




brks <- quantile(sample_stats$coef,
                 probs = seq(0, 1, 0.1),
                 na.rm = TRUE)


brks_fmt <- format(round(brks, 2),
                   nsmall = 2,
                   scientific = FALSE)

xxx <-  x_clean <- gsub(" ", "", head(brks_fmt, -1))
labels <- paste0(
  "[",
  gsub(" ", "",head(brks_fmt, -1)),
  " - ",
  gsub(" ", "",tail(brks_fmt, -1)),
  "]"
)


sample_stats_cv <- sample_stats %>%  
  ungroup() %>%   
  mutate(cv_decile = cv_decile <- cut(
    coef,
    breaks = brks,
    include.lowest = TRUE,
    labels = labels
  )) 



greys <- gray(seq(0.1, 0.9, length.out = 10))


bias_mapping <- ggplot(sample_stats_cv) +
  geom_sf(data = border, fill = "azure", color = "black", size = 0.2) +
  geom_sf(aes(fill = cv_decile), color = "grey40", size = 0.2) +
  scale_fill_manual(values = rev(greys), na.value = "#ff6961") +
  facet_grid(dim ~ dataset, switch = "y") +
  coord_sf(crs = eqearth_crs, expand = T) +
  theme_minimal(base_size = 10) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    strip.background = element_rect(fill = "black"),
    strip.text = element_text(color = "white", face = "bold", size = 12),
    panel.grid.major = element_line(color = "grey85", size = 0.2),
    legend.position = "bottom",
    legend.box = "horizontal",
    axis.text = element_blank(),
    axis.ticks = element_blank()
  ) +
  guides(fill = guide_legend(
    title = "Coefficient of variation (in %)",
    title.position = "top"
  ))


bias_mapping2 <- tag_facet2(bias_mapping) 

ggsave("bias_maps_CV_decile_geomsf.png", bias_mapping2, dpi = 600)







# Difference -------------------------------------------------------------------

list <- as_labeller(c("Species richness" =  "Species richness diff.", 
                      "functional diversity" = "Effective growth form diversity  diff.", 
                      "phylogenetic diversity" = "Phylogenetic diversity diff.", 
                      "Effective family diversity" ="Effective family diversity diff."))




palette <- c("darkorange", "#B576CC")  


box <- ggplot(sample_stats, aes(y = dim, x = rank_diff, fill = dataset)) +
  geom_jitter(aes(col = dataset), alpha = 0.3, position = position_jitterdodge(dodge.width = 0.7, 
                                                                               jitter.width = 0.1), 
              show.legend = F) +
  geom_boxplot( alpha = 0.01, show.legend = F,position = position_dodge(width = 0.7), 
                outlier.alpha = 0.1, outlier.colour = NA, width = 0.5) +
  geom_hline(yintercept = 0, col = "grey60", linetype = 2, linewidth = 0.7) +
  scale_color_manual(values = palette) +
  scale_x_continuous(limits = ggpmisc::symmetric_limits, name = "Average rank difference") +  
  scale_y_discrete( name = NULL,  guide = guide_axis(angle = 0),
                    position = "left", expand = expansion(mult = c(0,0))) +
  theme_bw(base_size = 14) +
  theme(
    axis.ticks.length = unit(0.25, "cm"),  
    axis.text.y = element_text(margin = margin(r = 0)),  
    axis.text.x = element_text(margin = margin(t = 0))   
  )

box


box_continent <- ggplot(sample_stats, aes(y = dim, x = rank_diff, group = interaction(dataset,dim))) +
  geom_jitter(aes(color = dataset), alpha = 0.3,
              position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.1)) +
  geom_boxplot(      
               color = "black",     
               alpha = 0.0001, 
               position = position_dodge(width = 0.7),
               outlier.alpha = 0.1,
               outlier.colour = NA,
               width = 0.5) +
  geom_vline(xintercept = 0, col = "grey60", linetype = 2, linewidth = 0.7) +
  scale_color_manual(values = palette)+
  facet_wrap(~LEVEL1_NAM, scales = "fixed", nrow = 3) +
  scale_x_continuous(limits = ggpmisc::symmetric_limits, name = "Average rank difference") +
  scale_y_discrete(name = NULL,
                   guide = guide_axis(angle = 0),
                   position = "left",
                   expand = expansion(mult = c(0,0))) +
  theme_bw(base_size = 14) +
  ggplot2::guides(col = ggplot2::guide_legend(override.aes = base::list(size = 5,
                                                                   alpha = 1))) +
  theme(strip.background = element_blank(), 
        strip.placement = "outside",
        strip.text = element_text(face = "bold"), 
        axis.ticks.length = unit(0.25, "cm"),
        legend.title = element_blank(), 
        legend.text = element_text(size = 14)) 

box_continent_tag <- tag_facet2(box_continent, tag_pool = letters[-1])


lay <- c(
  "#AAA#
  BBBBB
  BBBBB
  "
)

all_boxes <- tag_facet(box) / free(box_continent_tag) + plot_layout(axis_titles = "collect", guides = "collect",  design = lay) &
  theme(legend.direction = "horizontal", legend.byrow = T, legend.position = "top") 



ggsave("bias_boxes_seed_left.png", all_boxes, width = 10, height = 13, dpi = 600)




