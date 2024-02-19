library(tidyverse)

for(i in 1:100){
  print(paste0("iter",i))
  if(i ==1 ){
    results_refmodel <- readRDS("results/reference_model/scores_reference_results_approx1.rds")
  }else{
    if(file.exists(paste0("results/reference_model/scores_reference_results_approx",i,".rds"))){
      read_models <- readRDS(paste0("results/reference_model/scores_reference_results_approx",i,".rds"))
      results_refmodel <- rbind(results_refmodel, read_models) 
    }else{
      print(i)
    }
  }
}

results_refmodel <- results_refmodel %>%
  mutate(value = ifelse(score == "CRPS", -value, value),
         score = ifelse(score == "CRPS", "-CRPS",score))

for(j in 1:4){
  for(i in 1:100){
    print(paste0("model no",j, "iter",i))
    if(j==1 & i ==1 ){
      results_df_approx <- readRDS("results/model1/scores_validation_results_approx1.rds")
    }else{
      if(file.exists(paste0("results/model",j,"/scores_validation_results_approx",i,".rds"))){
        read_models_approx <- readRDS(paste0("results/model",j,"/scores_validation_results_approx",i,".rds"))
        results_df_approx <- rbind(results_df_approx, read_models_approx) 
      }
    }
  }
}

results_df_approx <- results_df_approx %>%
  mutate(value = ifelse(score == "CRPS", -value, value),
         score = ifelse(score == "CRPS", "-CRPS",score))

psis_score <- results_df_approx %>%
  filter(type_of_score == "MRP CELLWISE") %>%
  rename(psis_approx = "value") %>%
  select(-method, -type_of_score)

results_refmodel %>%
  mutate(model = factor(model))%>%
  mutate(model = fct_recode(model, "y_count | trials(n_j) ~ (1 | X1) + (1 | X3) + (1 | X2)" = "y_count | trials(n_j) ~ (1 | X1) + (1 | X2) + (1 | X3)")) %>%
  left_join(psis_score)%>%
  ggplot(aes(x = value, y = psis_approx, colour = model))+
  geom_abline() +
  geom_point(size = 1, alpha = .7)+
  facet_wrap(type_of_score~score, scales = "free")+
  theme_bw()+
  ggthemes::scale_color_colorblind()+
  xlab("Reference model score (Mk)") + ylab("PSIS Score(Mk)")+
  guides(color = guide_legend(nrow = 2))+
  theme(legend.position = "bottom", 
        legend.title = element_blank())

ggsave("figures/fullreference_vs_psisloco.png", width = 20, height = 12, units = "cm")
