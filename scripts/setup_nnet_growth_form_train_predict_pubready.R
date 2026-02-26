# Libraries --------------------------------------------------------------------
library(tidyverse)
library(data.table)
library(sf)
library(ranger)
library(missRanger)
library(missForest)
library(mice)
library(miceRanger)
library(rcompanion)
library(vcd)
library(mltools)
library(regrrr)
library(datawizard)
library(regrrr)
library(brms)
library(tidymodels)


which.max.adj <- function(x) {
  
  max <- which.max(table(x))
  if(is_empty(max)) {
    name <- NA
    return(name)
  } else 
    name <- names(max)
  return(name)
}
# Import data ------------------------------------------------------------------
plants_full <- fread("data/wcvp_accepted_merged.txt") 
dist_native <- fread("data/dist_native.txt") 


# dataset from Taylor, A., Weigelt, P., Denelle, P., Cai, L. and Kreft, H. (2023), The contribution of plant life and growth forms to global gradients of vascular plant diversity. New Phytol, 240: 1548-1560. https://doi.org/10.1111/nph.19011

growthform <- read.csv("data/traits_all/traits_all.csv")



# manipulate data --------------------------------------------------------------



wcvp_accepted_merger_filt <- plants_full %>% 
  filter(taxon_name %in% growthform$work_species)

growth_gift <- growthform %>% 
  dplyr::select(work_species, growth_form) %>% 
  right_join(wcvp_accepted_merger_filt, by = c("work_species" = "taxon_name")) %>% 
  filter(!duplicated(.))

dups_growth_gift_names <- data_duplicated(growth_gift, select = "plant_name_id") %>% 
  group_by(plant_name_id) %>% 
  mutate(combo=paste(growth_form, collapse='-'), 
         growth_form = case_when(grepl("climber", combo) ~ "climber",
                                 grepl("herb", combo) ~ "herb",
                                 grepl("subshrub", combo) ~ "subshrub",
                                 grepl("shrub", combo) ~ "shrub",
                                 
                                 TRUE ~ "check" )) %>% 
  ungroup() %>%
  dplyr::select(-one_of("V1", "Row", "count_na")) %>% 
  distinct() 



gift_growth_merge <- growth_gift %>% 
  filter(!plant_name_id %in% dups_growth_gift_names$plant_name_id) %>% 
  mutate(combo = growth_form) %>% 
  dplyr::select(-one_of("V1")) %>% 
  rbind(dups_growth_gift_names)




gift_growth_join <- gift_growth_merge %>% 
  dplyr::select(growth_gift = growth_form , plant_name_id)

data <- plants_full %>% 
  left_join(gift_growth_join, by = "plant_name_id") 



data_ref <- data %>% 
  filter(growth_gift == "check" | is.na(growth_gift)) %>% 
  mutate(growth_gift = NA) 


data_merger <- data %>%  
  filter(!plant_name_id %in% data_ref$plant_name_id) %>% 
  rbind(data_ref) 

round(nrow(data_ref) / nrow(data_merger),3) # prop na for gift growth



additional_parameters_dist <- dist_native %>%  
  group_by(plant_name_id)  %>% 
  summarize (major_continent = which.max.adj(continent), 
             major_region = which.max.adj(region))

data_model <-  data_merger %>% 
  group_by(genus) %>% 
  summarize (major_growth = which.max.adj(growth_gift)) %>% 
  left_join(data_merger, by = "genus") %>% 
  left_join(additional_parameters_dist, by = "plant_name_id") %>% 
  dplyr::select(plant_name_id, genus, growth_gift, family, major_growth, 
                major_continent, major_region) %>% 
  mutate(across(c(2:6), as.factor))


# modelling and setup data -----------------------------------------------------


set.seed(123)


sample <- data_model %>%
  filter(!plant_name_id %in% data_ref$plant_name_id) %>% 
  slice_sample(prop = 0.222) # corresponding to gift growth form na proportion

data_sample <- data_model %>%  
  filter(!plant_name_id %in% data_ref$plant_name_id) %>% 
  mutate(growth_gift = ifelse(plant_name_id %in% sample$plant_name_id, NA, growth_gift)) 


data_to_predict <- data_model %>%  
  filter(is.na(growth_gift))


#  modelling -------------------------------------------------------------------

# use data_sample for training and data_model for prediction
model_fit <- nnet::multinom(growth_gift ~ family + major_region + major_growth , data = data_model, 
                            MaxNWts = 5000, maxit =1000)
glance(model_fit)

predict_growth <- as.data.frame(predict(model_fit, newdata = data_model, "probs")) %>% 
  mutate_if(is.numeric,
            round,
            digits = 3)

data_model$pred <- colnames(predict_growth)[max.col(predict_growth,ties.method="random")]


test <- data_model %>% 
  filter(is.na(pred))

model_fit2 <- nnet::multinom(growth_gift ~ family + major_region , data = data_model, 
                             MaxNWts = 5000, maxit =1000)


predict_growth2 <- as.data.frame(predict(model_fit2, newdata = data_model, "probs")) %>% 
  mutate_if(is.numeric,
            round,
            digits = 3)

data_model$pred2 <- colnames(predict_growth2)[max.col(predict_growth2,ties.method="random")]


test <- data_model %>% 
  filter(is.na(pred2))

table(data_model$pred2)
table(data_model$pred)


growth_forms_imputed <- data_model %>% 
  mutate(across(c(2:8), as.character)) %>% 
  mutate(
    filled1 = ifelse(is.na(pred), "yes", "no"), 
    filled2 = ifelse(is.na(growth_gift), "yes", "no"), 
    filled = paste0(filled1, filled2), 
    pred = ifelse(is.na(pred), pred2, pred)
  ) %>% 
  mutate(growth_gift = ifelse(is.na(growth_gift), pred, growth_gift)) %>%
  dplyr::select(growth_form = growth_gift, predicted_growth_form = pred,  pred2, plant_name_id, genus, family, filled)

table(growth_forms_imputed$filled)

growth_forms_imputed_as <- growth_forms_imputed %>% 
  dplyr::select(growth_gift_full = growth_form, pred = predicted_growth_form, pred2, plant_name_id) %>% 
  left_join(data, by = "plant_name_id") %>% 
  mutate(growth_gift = as.factor(growth_gift), 
         pred = as.factor(pred), 
         pred2 = as.factor(pred2))

# extract data -----------------------------------------------------------------


x <- conf_mat(growth_forms_imputed_as, truth = "growth_gift", estimate = "pred")[[1]]
round(prop.table(x, margin = 2),2)

x <- conf_mat(growth_forms_imputed_as, truth = "growth_gift", estimate = "pred2")[[1]]
round(prop.table(x, margin = 2),2)

model1 <- summary(conf_mat(growth_forms_imputed_as, truth = "growth_gift", estimate = "pred")) %>% 
  mutate(model = "model_family + major_region + major_growth ")
model2 <-summary(conf_mat(growth_forms_imputed_as, truth = "growth_gift", estimate = "pred2"))%>% 
  mutate(model = "model_family + major_region")

model_perf <- rbind(model1, model2)


write.csv(model_perf, "performance_growth_forms_nnet_2026.csv")



autoplot(conf_mat(growth_forms_imputed_as, truth = "growth_gift", estimate = "pred"), type = "heatmap", )
autoplot(conf_mat(growth_forms_imputed_as, truth = "growth_gift", estimate = "pred2"), type = "heatmap")

growth_forms_imputed_csv <- growth_forms_imputed %>% 
  dplyr::select(growth_form, plant_name_id)


write.table(growth_forms_imputed_csv, "growth_forms_nnet_2026.txt")



