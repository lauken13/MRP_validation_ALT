library(tidyverse)

for(j in 1:4){
  for(i in 1:100){
    print(paste0("model no",j, "iter",i))
    if(j==1 & i ==1 ){
      results_df <- readRDS("results/model1/scores_validation_results1.rds")
    }else{
      if(file.exists(paste0("results/model",j,"/scores_validation_results",i,".rds"))){
        read_models <- readRDS(paste0("results/model",j,"/scores_validation_results",i,".rds"))
        results_df <- rbind(results_df, read_models) 
      }
    }
  }
}

results_df <- results_df %>%
  filter(method != "APPROX LOCO")

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

results_df <- rbind(results_df, results_df_approx)

results_df <- results_df %>%
  mutate(value = ifelse(score == "CRPS", -value, value),
         score = ifelse(score == "CRPS", "-CRPS",score))

#Compare exact population methods

true_score <- results_df %>%
  filter(method == "EXACT POPULATION" & type_of_score == "TRUE MRP")%>%
  group_by(model, score, iter)%>%
  summarise(mean_truth = mean(value))%>% #small variations across model runs
  ungroup()

comparison_score <- results_df %>%
  filter(method == "EXACT POPULATION" & type_of_score == "MRP CELLWISE")%>%
  left_join(true_score)

ggplot(comparison_score, aes(x = value, y = mean_truth, colour = model)) +
  geom_abline()+
  geom_point(size = 1, alpha = .7)+
  facet_wrap(.~score, scales = "free")+
  theme_bw()+
  ggthemes::scale_color_colorblind()+
  xlab("Cellwise formula") + ylab("Population score")+
  guides(color = guide_legend(nrow = 2))+
  theme(legend.position = "bottom", 
        legend.title = element_blank())

ggsave("figures/proof_of_score.png", width = 20, height = 12, units = "cm")


#Compare with mean of score method

true_score <- results_df %>%
  filter(method == "EXACT POPULATION" & type_of_score == "TRUE MRP")%>%
  group_by(model, score, iter)%>%
  summarise(mean_truth = mean(value))%>% #small variations across model runs
  ungroup()

comparison_score <- results_df %>%
  filter(method == "EXACT POPULATION" & type_of_score == "MEAN CELLWISE")%>%
  left_join(true_score)

ggplot(comparison_score, aes(x = value, y = mean_truth, colour = model)) +
  geom_abline()+
  geom_point(size = 1, alpha = .7)+
  facet_wrap(.~score, scales = "free")+
  theme_bw()+
  ggthemes::scale_color_colorblind()+
  xlab("Mean of cellwise score") + ylab("Population score")+
  guides(color = guide_legend(nrow = 2))+
  theme(legend.position = "bottom", 
        legend.title = element_blank())

ggsave("figures/meancellwise_vs_population.png", width = 20, height = 12, units = "cm")


#Compare with mean of score to true mean of score

true_score_cellwise <- results_df %>%
  filter(method == "EXACT POPULATION" & type_of_score == "MEAN CELLWISE")%>%
  group_by(model, score, iter)%>%
  summarise(mean_truth = mean(value))%>% #small variations across model runs
  ungroup()

comparison_score_cellwise <- results_df %>%
  filter(method == "APPROX LOCO" & type_of_score == "MEAN CELLWISE")%>%
  left_join(true_score_cellwise)

ggplot(comparison_score_cellwise, aes(x = value, y = mean_truth, colour = model)) +
  geom_abline()+
  geom_point(size = 1, alpha = .7)+
  facet_wrap(.~score, scales = "free")+
  theme_bw()+
  ggthemes::scale_color_colorblind()+
  xlab("Mean of approx cellwise score") + ylab("Sum of cellwise scores in population")+
  guides(color = guide_legend(nrow = 2))+
  theme(legend.position = "bottom", 
        legend.title = element_blank())

ggsave("figures/meancellwise_vs_populationcellwise.png", width = 20, height = 12, units = "cm")


#Compare scores using sample as proxxy for population 
sample_score <- results_df %>%
  filter(method %in% c("SAMPLE ESTIMATE")  & type_of_score == "MRP CELLWISE")%>%
  left_join(true_score)

ggplot(sample_score, aes(x = value, y = mean_truth, colour = model)) +
  geom_abline()+
  geom_point(size = 1, alpha = .7)+
  facet_wrap(.~score, scales = "free")+
  theme_bw()+
  ggthemes::scale_color_colorblind()+
  xlab("Sample estimate") + ylab("Truth")+
  guides(color = guide_legend(nrow = 2))+
  theme(legend.position = "bottom", 
        legend.title = element_blank())

ggsave("figures/sample_score.png", width = 20, height = 12, units = "cm")

#Brute force LOCO
bruteforce_loco <- results_df %>%
  filter(method == "BRUTE FORCE LOCO" & type_of_score == "MRP CELLWISE")%>%
  left_join(true_score)

ggplot(bruteforce_loco, aes(x = value, y = mean_truth, colour = model)) +
  geom_abline() +
  geom_point(size = 1, alpha = .7)+
  facet_wrap(.~score, scales = "free")+
  theme_bw()+
  ggthemes::scale_color_colorblind()+
  xlab("Brute-force LOCO") + ylab("True Score")+
  guides(color = guide_legend(nrow = 2))+
  theme(legend.position = "bottom", 
        legend.title = element_blank())

ggsave("figures/bruteforce_vs_truth.png", width = 20, height = 12, units = "cm")

#Does brute force loco deviate from just using the sample score
sample_score <- results_df %>%
  filter(method %in% c("SAMPLE ESTIMATE", "BRUTE FORCE LOCO") & type_of_score == "MRP CELLWISE")%>%
  pivot_wider(names_from = method, values_from = value)

ggplot(sample_score, aes(x = `SAMPLE ESTIMATE`, y = `BRUTE FORCE LOCO`, colour = model)) +
  geom_abline() +
  geom_point(size = 1, alpha = .7)+
  facet_wrap(.~score, scales = "free")+
  theme_bw()+
  ggthemes::scale_color_colorblind()+
  xlab("Sample Estimate") + ylab("Brute force LOCO")+
  guides(color = guide_legend(nrow = 2))+
  theme(legend.position = "bottom", 
        legend.title = element_blank())

ggsave("figures/bruteforce_vs_sample_score.png", width = 20, height = 12, units = "cm")

#PSIS LOCO
approx_loco <- results_df %>%
  filter(method == "APPROX LOCO" & type_of_score == "MRP CELLWISE")%>%
  left_join(true_score)

ggplot(approx_loco, aes(x = value, y = mean_truth, colour = model)) +
  geom_point(size = 1, alpha = .7)+
  facet_wrap(type_of_score~score, scales = "free")+
  theme_bw()+
  ggthemes::scale_color_colorblind()+
  xlab("Estimate") + ylab("Truth")+
  guides(color = guide_legend(nrow = 2))+
  theme(legend.position = "bottom", 
        legend.title = element_blank())


#Does psis deviate from just using the brute force
sample_score <- results_df %>%
  filter(method %in% c("APPROX LOCO", "BRUTE FORCE LOCO"),
         type_of_score %in% c("MRP CELLWISE"))%>%
  pivot_wider(names_from = method, values_from = value)

sample_score %>%
  filter(type_of_score == "MRP CELLWISE")%>%
  ggplot(., aes(x = `BRUTE FORCE LOCO`, y = `APPROX LOCO`, colour = model)) +
  geom_abline() +
  geom_point(size = 1, alpha = .7)+
  facet_wrap(type_of_score~score, scales = "free")+
  theme_bw()+
  ggthemes::scale_color_colorblind()+
  xlab("Brute-force LOCO") + ylab("PSIS-LOCO")+
  guides(color = guide_legend(nrow = 2))+
  theme(legend.position = "bottom", 
        legend.title = element_blank())

ggsave("figures/bruteforce_vs_psisloco.png", width = 20, height = 12, units = "cm")
