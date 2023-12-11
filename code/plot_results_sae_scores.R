library(tidyverse)

for(i in 1:100){
  if(i ==1 ){
    results_saemodel <- readRDS("results/sae_scores/scores_sae_approx1.rds")
  }else{
    if(file.exists(paste0("results/sae_scores/scores_sae_approx",i,".rds"))){
      read_models <- readRDS(paste0("results/sae_scores/scores_sae_approx",i,".rds"))
      results_saemodel <- rbind(results_saemodel, read_models) 
    }else{
      print(paste0("there's an issue with iteration", i))
    }
  }
}

true_score <- results_saemodel %>%
  filter(method == "EXACT POPULATION" & type_of_score == "TRUE MRP")%>%
  rename(true_score = value)%>%
  select(-type_of_score, - method)

comparison_score <- results_saemodel %>%
  filter(method == "APPROX LOCO" & type_of_score == "MRP CELLWISE")%>%
  left_join(true_score)

comparison_score %>%
  mutate(variable = forcats::fct_recode(variable,`SAE: X1` = "X1",
                                      `SAE: X2` = "X2",
                                      `SAE: X3` = "X3",
                                      `SAE: X4` = "X4"))%>%
  rename(Model = model, `SAE Level` = level)%>%
ggplot(aes(x = value, y = true_score, shape = `SAE Level`, colour = Model))+
  geom_abline(slope = 1, intercept = 0)+
  geom_point(size = 1, alpha = .7)+
  facet_wrap(variable~score, scales = "free", ncol = 4)+
  theme_bw()+
  ggthemes::scale_color_colorblind()+  
  guides(color = guide_legend(nrow = 6),shape = guide_legend(nrow = 5))+
  theme(legend.position = "bottom")+
  xlab("LOCO estimated score")+
  ylab("True score")

ggsave("figures/saelevels_truth_vs_psisloco.png", width = 20, height = 15, units = "cm")


comparison_score %>%
  group_by(iter, model, variable,score, type_of_score,method) %>%
  summarise(mean_of_score = mean(value),
            mean_of_truescore = mean(true_score))%>%
  ungroup()  %>%
  mutate(variable = forcats::fct_recode(variable,`SAE: X1` = "X1",
                                        `SAE: X2` = "X2",
                                        `SAE: X3` = "X3",
                                        `SAE: X4` = "X4"))%>%
  ggplot(aes(x = mean_of_score, y = mean_of_truescore,  colour = model))+
  geom_abline(slope = 1, intercept = 0)+
  geom_point(size = 1, alpha = .7)+
  facet_wrap(variable~score, scales = "free", ncol = 4)+
  theme_bw()+
  ggthemes::scale_color_colorblind()+  
  guides(color = guide_legend(nrow = 6),shape = guide_legend(nrow = 5))+
  theme(legend.position = "bottom", legend.title = element_blank())+
  xlab("Mean LOCO estimated score")+
  ylab("Mean true score")


ggsave("figures/sae_sumk_truth_vs_psisloco.png", width = 20, height = 15, units = "cm")

comparison_score %>%
  filter(score == "SQUARED ERROR")%>%
  group_by(iter, model, variable,score, type_of_score,method) %>%
  summarise(mean_of_score = mean(value),
            mean_of_truescore = mean(true_score))%>%
  ungroup()  %>%
  mutate(variable = forcats::fct_recode(variable,`SAE: X1` = "X1",
                                        `SAE: X2` = "X2",
                                        `SAE: X3` = "X3",
                                        `SAE: X4` = "X4"))%>%
  ggplot(aes(x = mean_of_score, y = mean_of_truescore,  colour = model))+
  geom_abline(slope = 1, intercept = 0)+
  geom_point(size = 1, alpha = .7)+
  facet_wrap(variable~., scales = "free", ncol = 2)+
  theme_bw()+
  ggthemes::scale_color_colorblind()+  
  guides(color = guide_legend(nrow = 6),shape = guide_legend(nrow = 5))+
  theme(legend.position = "right", legend.title = element_blank())+
  xlab("Mean LOCO estimated score")+
  ylab("Mean true score")


ggsave("figures/slideversion_squarederror_sae_sumk_truth_vs_psisloco.png", width = 20, height = 10, units = "cm")

