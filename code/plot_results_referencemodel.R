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

diff_scores <- results_df_approx %>%
  filter(type_of_score == "MRP CELLWISE") %>%
  pivot_wider(names_from = "model", values_from = "value") %>%
  mutate(`y_count | trials(n_j) ~ (1 | X1) + (1 | X3)` = 
           `y_count | trials(n_j) ~ (1 | X1) + (1 | X3)` - `y_count | trials(n_j) ~ (1 | X1) + (1 | X2) + (1 | X3) + (1 | X4)`,
         `y_count | trials(n_j) ~ (1 | X1) + (1 | X2) + (1 | X3)` = 
           `y_count | trials(n_j) ~ (1 | X1) + (1 | X3) + (1 | X2)` - `y_count | trials(n_j) ~ (1 | X1) + (1 | X2) + (1 | X3) + (1 | X4)`,
         `y_count | trials(n_j) ~ (1 | X1) + (1 | X3) + (1 | X4)` = 
          `y_count | trials(n_j) ~ (1 | X1) + (1 | X3) + (1 | X4)` -  `y_count | trials(n_j) ~ (1 | X1) + (1 | X2) + (1 | X3) + (1 | X4)`)%>%
  select(-c(`y_count | trials(n_j) ~ (1 | X1) + (1 | X2) + (1 | X3) + (1 | X4)`,`y_count | trials(n_j) ~ (1 | X1) + (1 | X3) + (1 | X2)`))%>%
  pivot_longer(`y_count | trials(n_j) ~ (1 | X1) + (1 | X3) + (1 | X4)`:`y_count | trials(n_j) ~ (1 | X1) + (1 | X2) + (1 | X3)`, names_to = "model", values_to = "approx_score_diff")%>%
  select(score, iter, model, approx_score_diff)

results_refmodel %>%
  left_join(diff_scores)%>%
  ggplot(aes(x = value, y = approx_score_diff, colour = model))+
  geom_abline() +
  geom_point(size = 1, alpha = .7)+
  facet_wrap(type_of_score~score, scales = "free")+
  theme_bw()+
  ggthemes::scale_color_colorblind()+
  xlab("Reference model score") + ylab("PSIS Score(Mk) - PSIS Score(M*)")+
  guides(color = guide_legend(nrow = 4))+
  theme(legend.position = "bottom", 
        legend.title = element_blank())

ggsave("figures/fullreference_vs_psisloco.png", width = 20, height = 12, units = "cm")
