
# Authors: Ludwig Baldaszti - lbaldaszti@rbge.org.uk
# Date: 24/02/2026

# Libraries --------------------------------------------------------------------
library(tidyverse)
library(data.table)
library(plotrix)
library(paletteer)
library(multcompView)
library(FSA)
library(ggplot2)
library(rstatix)
library(ggpubr)
library(nationalparkcolors)
library(patchwork) 

options(scipen = 999999)


error <- function(x) {
  qnorm(0.975)*(sd(x)/sqrt(length(x)))
}




# data -------------------------------------------------------------------------

a <- nationalparkcolors::park_palette("Acadia")
palette <- a[c(1,2,3,5)]

threshold_files <-  list.files(path = "output/threshold/median/", full.names = T)
threshold <- lapply(threshold_files, fread) %>% 
  rbindlist() %>% 
  filter(continent == "GLOBAL") %>%
  mutate(continent = as.factor(continent), 
         #index = as.factor(index), 
         index = factor(x = as.character(index), levels = c("diversity_sp", "diversity_phy",  "diversity_fun", "diversity_tax"), 
                        labels = c("Species richness", "Phylogenetic diversity", "Effective growth form diversity",  "Effective family diversity")))  %>% 
  filter(source == "wcvp")
levels(threshold$index)




global_threshold_line <- threshold %>% 
  group_by(index) %>% 
  summarise( mean = mean(sp))



# stats ------------------------------------------------------------------------


dunn_res <- dunnTest(global_threshold$sp,  global_threshold$index,  method = "bonferroni")


res <- dunn_res$res

pv <- res$P.adj
names(pv) <- gsub(" - ", "-", res$Comparison)  

cld <- multcompLetters(pv, threshold = 0.05)

letters_df <- data.frame(group = levels(global_threshold$index),
                         letters = c("a", "b","c", "d")) %>% 
  mutate(index = factor(x = as.character(group), levels = c("Species richness", "Phylogenetic diversity", "Effective growth form diversity",  "Effective family diversity")))


# thresholds -------------------------------------------------------------------


global_plot <- ggplot(global_threshold, aes(x = reorder(index,sp), y = log10(sp))) +
  geom_violin(aes(fill = index, col = "black"), alpha = 1, col = "black", show.legend = F, lwd = 0.4,  
              position = position_dodge(width = 0.5), lty = 1) +
  geom_text(data = letters_df, aes(x = index, y = log10(max(global_threshold$sp)) +0.5, label = letters),
            size = 5) +
  scale_fill_manual(values = palette) +
  scale_y_continuous(name = "Number of species included (log10)", breaks = c(0,1,2,3,4,5, 6), limits = c(0,6)) +  
  scale_x_discrete(name = NULL) +
  theme_bw(base_size = 12)  +
  labs(fill = "Diversity metric") + 
  coord_flip()

global_plot



global_summary <- threshold %>% 
  group_by(index,  source) %>% 
  summarise(mean = mean(sp), 
            min = min(sp), 
            max = max(sp), 
            median = median(sp),
            mad = mad(sp), 
            se = std.error(sp), 
            upper_CI = mean + error(sp),
            lower_CI = mean - error(sp), 
            se_upper = mean + se, 
            se_lower = mean - se) %>% 
  mutate(across(where(is.numeric), ceiling)) 



# curves -----------------------------------------------------------------------


sample_data <- open_dataset(
  sources = "output/samples_aggregated/hillr/", 
  format = "parquet"
)


samples_cumulative_rel <- sample_data %>% 
  collect() %>% 
  setDT() %>% 
  mutate(index = factor(index, levels = c("diversity_sp", "diversity_phy", "diversity_fun",  "diversity_tax"), 
                        labels = c("Species richness", "Phylogenetic diversity", 
                                   "Effective growth form diversity", "Effective family diversity"))) 



curves_global_x <-  ggplot(data = samples_cumulative_rel, aes(x=log10(sp), y=mean_corgw)) +
  geom_hline(yintercept=0.95, linetype="longdash", color = alpha("red", 0.6), lwd = 0.5) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.95, ymax = Inf,
           alpha = 0.2, fill = "black") +
  ylab("Global correlation coefficient") +
  geom_ribbon(aes(ymin = lower_ci_corgw, ymax = upper_ci_corgw, fill = index), 
              outline.type = "both", 
              linetype = 0, alpha = 0.5, show.legend = F) +
  geom_line(aes(col = index), alpha = 1, lwd = 0.5,  show.legend = F) +
  geom_point(aes(col = index, fill = NULL), alpha = 1, size = 0.5, stroke=0, show.legend = F) +
  scale_x_continuous(name = "Number of species included (log10)", breaks = c(0,1,2,3,4,5, 6), limits = c(0,6)) +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  theme_bw(base_size = 12) 

curves_global_x
combined_global <- curves_global_x / global_plot +  
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") & theme(legend.position = 'right',  
                                                                                legend.text = element_text(size=12), 
                                                                                plot.tag.position = c(0.3,0.96), 
                                                                                plot.tag = element_text(face = "bold", size = 14), 
                                                                                plot.tag.location = "plot")

ggsave(paste0("Global_average_sp_threshold_median_violin_curves_sr", unique(threshold$source), "nodens_sameaxes.png"),
       combined_global,
       dpi = 1200,
       height = 11,
       width = 8)

         