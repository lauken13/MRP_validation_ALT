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

#Compare exact population methods

true_score <- results_df %>%
  filter(method == "EXACT POPULATION" & type_of_score == "TRUE MRP")%>%
  group_by(model, score, iter)%>%
  summarise(mean_truth = mean(value))%>% #small variations across model runs
  ungroup()

comparison_score <- results_df %>%
  filter(method == "EXACT POPULATION" & type_of_score != "TRUE MRP")%>%
  left_join(true_score)

ggplot(comparison_score, aes(x = value, y = mean_truth, colour = model)) +
  geom_abline()+
  geom_point(size = 1, alpha = .7)+
  facet_wrap(type_of_score~score, scales = "free")+
  theme_bw()+
  ggthemes::scale_color_colorblind()+
  xlab("Estimate") + ylab("Truth")+
  guides(color = guide_legend(nrow = 4))+
  theme(legend.position = "bottom", 
        legend.title = element_blank())

#Compare scores using sample as proxxy for population 
sample_score <- results_df %>%
  filter(method == "SAMPLE ESTIMATE")%>%
  left_join(true_score)

ggplot(sample_score, aes(x = value, y = mean_truth, colour = model)) +
  geom_abline()+
  geom_point(size = 1, alpha = .7)+
  facet_wrap(type_of_score~score, scales = "free")+
  theme_bw()+
  ggthemes::scale_color_colorblind()+
  xlab("Estimate") + ylab("Truth")+
  guides(color = guide_legend(nrow = 4))+
  theme(legend.position = "bottom", 
        legend.title = element_blank())

#Brute force LOCO
bruteforce_loco <- results_df %>%
  filter(method == "BRUTE FORCE LOCO")%>%
  left_join(true_score)

ggplot(bruteforce_loco, aes(x = value, y = mean_truth, colour = model)) +
  geom_point(size = 1, alpha = .7)+
  facet_wrap(type_of_score~score, scales = "free")+
  theme_bw()+
  ggthemes::scale_color_colorblind()+
  xlab("Estimate") + ylab("Truth")+
  guides(color = guide_legend(nrow = 4))+
  theme(legend.position = "bottom", 
        legend.title = element_blank())

#Does brute force loco deviate from just using the sample score
sample_score <- results_df %>%
  filter(method %in% c("SAMPLE ESTIMATE", "BRUTE FORCE LOCO"))%>%
  pivot_wider(names_from = method, values_from = value)

ggplot(sample_score, aes(x = `SAMPLE ESTIMATE`, y = `BRUTE FORCE LOCO`, colour = model)) +
  geom_abline() +
  geom_point(size = 1, alpha = .7)+
  facet_wrap(type_of_score~score, scales = "free")+
  theme_bw()+
  ggthemes::scale_color_colorblind()+
  xlab("SAMPLE ESTIMATE") + ylab("BRUTE FORCE LOCO")+
  guides(color = guide_legend(nrow = 4))+
  theme(legend.position = "bottom", 
        legend.title = element_blank())

sample_score %>%
  filter(score == "SQUARED ERROR" & type_of_score == "MRP CELLWISE")%>%
  ggplot(., aes(x = `SAMPLE ESTIMATE`, y = `BRUTE FORCE LOCO`, colour = model)) +
  geom_abline() +
  geom_point(size = 1, alpha = .7)+
  coord_fixed() +
  facet_wrap(type_of_score~score)+
  theme_bw()+
  ggthemes::scale_color_colorblind()+
  xlab("SAMPLE ESTIMATE") + ylab("BRUTE FORCE LOCO")+
  guides(color = guide_legend(nrow = 4))+
  theme(legend.position = "bottom", 
        legend.title = element_blank())


sample_score %>%
  filter(score == "CRPS" & type_of_score == "MRP CELLWISE")%>%
  ggplot(., aes(x = `SAMPLE ESTIMATE`, y = `BRUTE FORCE LOCO`, colour = model)) +
  geom_abline() +
  geom_point(size = 1, alpha = .7)+
  coord_fixed() +
  facet_wrap(type_of_score~score)+
  theme_bw()+
  ggthemes::scale_color_colorblind()+
  xlab("SAMPLE ESTIMATE") + ylab("BRUTE FORCE LOCO")+
  guides(color = guide_legend(nrow = 4))+
  theme(legend.position = "bottom", 
        legend.title = element_blank())

#PSIS LOCO
approx_loco <- results_df %>%
  filter(method == "APPROX LOCO")%>%
  left_join(true_score)

ggplot(approx_loco, aes(x = value, y = mean_truth, colour = model)) +
  geom_point(size = 1, alpha = .7)+
  facet_wrap(type_of_score~score, scales = "free")+
  theme_bw()+
  ggthemes::scale_color_colorblind()+
  xlab("Estimate") + ylab("Truth")+
  guides(color = guide_legend(nrow = 4))+
  theme(legend.position = "bottom", 
        legend.title = element_blank())


#Does psis deviate from just using the brute force
sample_score <- results_df %>%
  filter(method %in% c("APPROX LOCO", "BRUTE FORCE LOCO"))%>%
  pivot_wider(names_from = method, values_from = value)

sample_score %>%
  filter(score == "CRPS" & type_of_score == "MRP CELLWISE")%>%
  ggplot(., aes(x = `BRUTE FORCE LOCO`, y = `APPROX LOCO`, colour = model)) +
  geom_abline() +
  geom_point(size = 1, alpha = .7)+
  coord_fixed() +
  facet_wrap(type_of_score~score)+
  theme_bw()+
  ggthemes::scale_color_colorblind()+
  xlab("BRUTE FORCE LOCO") + ylab("PSIS LOCO")+
  guides(color = guide_legend(nrow = 4))+
  theme(legend.position = "bottom", 
        legend.title = element_blank())


sample_score %>%
  filter(score == "SQUARED ERROR" & type_of_score == "MRP CELLWISE")%>%
  ggplot(., aes(x = `BRUTE FORCE LOCO`, y = `APPROX LOCO`, colour = model)) +
  geom_abline() +
  geom_point(size = 1, alpha = .7)+
  coord_fixed() +
  facet_wrap(type_of_score~score)+
  theme_bw()+
  ggthemes::scale_color_colorblind()+
  xlab("BRUTE FORCE LOCO") + ylab("PSIS LOCO")+
  guides(color = guide_legend(nrow = 4))+
  theme(legend.position = "bottom", 
        legend.title = element_blank())


