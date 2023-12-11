library(tidyverse)

#summarise sample coverage
for(i in 1:100){
  print(paste0("iter",i))
  if(i ==1 ){
    sample_summary <- readRDS("results/combined_reference_model/sample_summary1.rds")
  }else{
    if(file.exists(paste0("results/combined_reference_model/sample_summary",i,".rds"))){
      read_models <- readRDS(paste0("results/combined_reference_model/sample_summary",i,".rds"))
      sample_summary <- rbind(sample_summary, read_models) 
    }else{
      print(i)
    }
  }
}
summary(sample_summary$prop_cells_covered)

#summarise results
for(i in 1:100){
  print(paste0("iter",i))
  if(i ==1 ){
    results_comb_ref_model <- readRDS("results/combined_reference_model/scores_reference_results_approx1.rds")
  }else{
    if(file.exists(paste0("results/combined_reference_model/scores_reference_results_approx",i,".rds"))){
      read_models <- readRDS(paste0("results/combined_reference_model/scores_reference_results_approx",i,".rds"))
      results_comb_ref_model <- rbind(results_comb_ref_model, read_models) 
    }else{
      print(i)
    }
  }
}

true_score <- results_comb_ref_model %>%
  filter(type_of_score == "TRUE MRP") %>%
  rename(true_score = "value") %>%
  select(-method, -type_of_score, -reference_model)

results_comb_ref_model %>%
  filter(type_of_score %in% c("PARTIAL MRP CELLWISE")) %>%
  pivot_wider(names_from = "method", values_from = "value") %>%
  ggplot(., aes(x = `APPROX LOCO REFERENCE`, y = `APPROX LOCO`, colour = model))+
  geom_point(size = 1, alpha = .7)+
  facet_wrap(score~reference_model, scales = "free")+
  xlab("LOCO Reference Model Approach")+
  ylab("LOCO Error Approach")+
  theme_bw()+
  ggthemes::scale_color_colorblind()+
  guides(color = guide_legend(nrow = 4))+
  theme(legend.position = "bottom",
        legend.title = element_blank())

ggsave("figures/partialreference_vs_partialloco.png", width = 15, height = 15, units = "cm")

results_comb_ref_model %>%
  filter(method %in% c("APPROX COMBINED LOCO-REFERENCE"), reference_model == "FULL") %>%
  left_join(true_score)%>%
  ggplot(., aes(x = value, y = true_score, colour = model))+
  geom_abline(intercept = 0, slope = 1 )+
  geom_point(size = 1, alpha = .7)+
  facet_wrap(.~score, scales = "free")+
  xlab("LOCO Combined Reference Model Approach")+
  ylab("True Error")+
  ggthemes::scale_color_colorblind()+
  guides(color = guide_legend(nrow = 4))+
  theme_bw()+
  theme(legend.position = "bottom",
        legend.title = element_blank())

ggsave("figures/combined_referenceloco_vs_truth.png", width = 20, height = 12, units = "cm")
