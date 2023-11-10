library(tidyverse)

for(i in 2:100){
  print(paste0("iter",i))
  if(i ==2 ){
    results_saemodel <- readRDS("results/sae_scores/scores_sae_approx2.rds")
  }else{
    if(file.exists(paste0("results/sae_scores/scores_sae_approx",i,".rds"))){
      read_models <- readRDS(paste0("results/sae_scores/scores_sae_approx",i,".rds"))
      results_saemodel <- rbind(results_saemodel, read_models) 
    }else{
      print(i)
    }
  }
}

true_score <- results_saemodel %>%
  filter(method == "EXACT POPULATION" & type_of_score == "TRUE MRP")%>%
  group_by(model, score, iter, level, variable)%>%
  summarise(mean_truth = mean(value))%>% #small variations across model runs
  ungroup()

comparison_score <- results_saemodel %>%
  filter(method == "APPROX LOCO" & type_of_score == "MRP CELLWISE")%>%
  left_join(true_score)

comparison_score %>%
  filter(score == "SQUARED ERROR")%>%
ggplot(aes(x = value, y = mean_truth, colour = level, shape = model))+
  geom_point(size = 1, alpha = .7)+
  facet_grid(variable~score, scales = "free")+
  theme_bw()+
  ggthemes::scale_color_colorblind()+  guides(color = guide_legend(nrow = 4))+
  theme(legend.position = "bottom", 
        legend.title = element_blank())

ggsave("figures/fullreference_vs_psisloco.png", width = 20, height = 12, units = "cm")
