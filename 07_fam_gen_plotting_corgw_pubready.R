
# Authors: Ludwig Baldaszti - lbaldaszti@rbge.org.uk
# Date: 24/02/2026

# Libraries --------------------------------------------------------------------
library(tidyverse)
library(data.table)
library(RColorBrewer)
library(ggdist)
library(arrow)
library(ggrepel)
library(paletteer)
library(patchwork)

options(scipen = 999)

paletteer_d("wesanderson::AsteroidCity3")
two_col <-   c("#FBA72AFF", "#5785C1FF") 

#  needs output from 02_samples_output_assembly_eddie_parquet_hillr
sample_data <- open_dataset(
  sources = "output/samples_aggregated/hillr/", 
  format = "parquet"
)



samples_cumulative_rel <- sample_data %>% 
  dplyr::select(mean_corgw, sp, lower_ci_corgw, upper_ci_corgw, continent, index, source, cov_corgw) %>% 
  filter(source == "wcvp") %>% 
  collect() %>% 
  setDT() %>% 
  mutate(index = factor(index, levels = c("diversity_sp", "diversity_phy", "diversity_fun",  "diversity_tax"), 
                        labels = c("species diversity", "phylogenetic diversity", 
                                   "Effective growth form diversity", "Effective family diversity"))) %>% 
  filter(index == "species diversity")



results <- fread("output/Corgw_rich_genera.txt")


results_global <- results %>%  
  filter(continent == "GLOBAL")


data_only_global2 <- results_global %>% 
  left_join(samples_cumulative_rel, by = c("sp_origin" = "sp")) %>% 
  mutate(
    difference_baseline = corsp_median - mean_corgw, 
    diff = ifelse(difference_baseline < 0, "lower", "higher"))

top5_better <- data_only_global2 %>% 
  filter(difference_baseline > 0) %>% 
  slice_max(sp_origin, n = 5)

top5 <- data_only_global2 %>% 
  slice_max(corsp_median, n = 5)

top5_dev_max <- data_only_global2 %>% 
  slice_max(difference_baseline, n = 10)

top5_dev_min <- data_only_global2 %>% 
  slice_min(difference_baseline, n = 10)

gen_label <- c(top5$gen, top5_dev_max$gen, top5_dev_min$gen, top5_better$gen)


data_only_global3 <- data_only_global2 %>% 
  mutate(
    lab = ifelse(gen %in% gen_label, gen, NA))


xx <- data_only_global3 %>% 
  filter(difference_baseline > 0) %>% 
  slice_max(difference_baseline, n = 10)

plot_gen <- ggplot() +
  geom_point(data = data_only_global3, aes(x = sp_source, y = corsp_median, col = diff), alpha = 0.8) +
  #geom_smooth() +
  geom_text_repel(data = data_only_global3, aes(x = sp_source, y = corsp_median, label = lab, fontface = "italic"), size = 3, 
                  min.segment.length = 0, 
                  seed = 42, 
                  box.padding = 0.25) +
   geom_hline(yintercept=0.95, linetype="longdash", color = alpha("red", 0.6), lwd = 0.5) +
   annotate("rect", xmin = 0, xmax = Inf, ymin = 0.95, ymax = Inf,
            alpha = 0.2, fill = "black") +
  #ylab("Correlation coefficient") +
  geom_ribbon(data = samples_cumulative_rel,
              aes(x=sp, y=mean_corgw, ymin = lower_ci_corgw, ymax = upper_ci_corgw), 
              outline.type = "both", 
              linetype = 0, alpha = 0.9, show.legend = F, fill = "black") +
  theme_bw(base_size = 14) +
  scale_color_discrete(palette =  two_col) +
  scale_x_continuous(trans = "log10", breaks = c(1,10, 100, 1000, 10000, 100000)) +

  labs(col = "Diff. to random", 
       y = "Global correlation coefficient", 
       x = "Number of species") 





results <- fread("output/Corgw_rich_family.txt")


results_global <- results %>%  
  filter(continent == "GLOBAL")


data_only_global2 <- results_global %>% 
  left_join(samples_cumulative_rel, by = c("sp_origin" = "sp")) %>% 
  mutate(
    difference_baseline = corsp_median - mean_corgw, 
    diff = ifelse(difference_baseline < 0, "lower", "higher"))

top5_better <- data_only_global2 %>% 
  filter(difference_baseline > 0) %>% 
  slice_max(sp_origin, n = 5)

top5 <- data_only_global2 %>% 
  slice_max(corsp_median, n = 5)

top5_dev_max <- data_only_global2 %>% 
  slice_max(difference_baseline, n = 10)

top5_dev_min <- data_only_global2 %>% 
  slice_min(difference_baseline, n = 10)

fam_label <- c(top5$fam, top5_dev_max$fam, top5_dev_min$fam, top5_better$fam)


data_only_global3 <- data_only_global2 %>% 
  mutate(
    lab = ifelse(fam %in% fam_label, fam, NA))



plot_fam <- ggplot() +
  geom_point(data = data_only_global3, aes(x = sp_source, y = corsp_median, col = diff), alpha = 0.8) +
  geom_text_repel(data = data_only_global3, aes(x = sp_source, y = corsp_median, label = lab, fontface = "plain"), size = 3, 
                  min.segment.length = 0, 
                  seed = 42, 
                  box.padding = 0.25,
                 ) +
  geom_hline(yintercept=0.95, linetype="longdash", color = alpha("red", 0.6), lwd = 0.5) +
  annotate("rect", xmin = 0, xmax = Inf, ymin = 0.95, ymax = Inf,
           alpha = 0.2, fill = "black") +

  geom_ribbon(data = samples_cumulative_rel,
              aes(x=sp, y=mean_corgw, ymin = lower_ci_corgw, ymax = upper_ci_corgw), 
              outline.type = "both", 
              linetype = 0, alpha = 0.9, show.legend = F, fill = "black") +
  theme_bw(base_size = 14) +
  scale_color_discrete(palette =  two_col) +
  scale_x_continuous(trans = "log10", breaks = c(1,10, 100, 1000, 10000, 100000)) +
  
  labs(col = "Diff. to random", 
       y = "Global correlation coefficient", 
       x = "Number of species") 





comb <- plot_fam / plot_gen  +
  patchwork::plot_layout(guides = "collect") +  # axes = "collect"
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") & theme(legend.position = 'right',  
                                            legend.text = element_text(size=12), 
                                            plot.tag.position = c(0.1,0.96), 
                                            plot.tag = element_text(face = "bold", size = 14), 
                                            plot.tag.location = "plot") & guides(color = guide_legend(override.aes = base::list(size=5))) 



ggsave("gen_fam_curves_comp_med.svg", comb, height = 12, width = 9)
ggsave("gen_fam_curves_comp_med.png", comb, height = 12, width = 9, dpi = 600)


