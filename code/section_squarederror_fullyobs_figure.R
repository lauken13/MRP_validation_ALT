library(tidyverse)
library(ggplot2)


results_bias <- readRDS(paste0("results/section3_1/simulation_3_1_biasmodel_iteration",1,".rds"))
for(ITE in c(2:100)){
  if(file.exists(paste0("results/section3_1/simulation_3_1_biasmodel_iteration",ITE,".rds"))){
    tmp <-readRDS(paste0("results/section3_1/simulation_3_1_biasmodel_iteration",ITE,".rds"))
    results_bias <- rbind(results_bias,tmp)
  }
}
results_precision <- readRDS(paste0("results/section3_1/simulation_3_1_precisionmodel_iteration",1,".rds"))
for(ITE in c(2:100)){
  if(file.exists(paste0("results/section3_1/simulation_3_1_precisionmodel_iteration",ITE,".rds"))){
    tmp <-readRDS(paste0("results/section3_1/simulation_3_1_precisionmodel_iteration",ITE,".rds"))
    results_precision <- rbind(results_precision,tmp)
  }
}
results_nuisance <- readRDS(paste0("results/section3_1/simulation_3_1_nuisancemodel_iteration",1,".rds"))
for(ITE in c(2:100)){
  if(file.exists(paste0("results/section3_1/simulation_3_1_nuisancemodel_iteration",ITE,".rds"))){
    tmp <-readRDS(paste0("results/section3_1/simulation_3_1_nuisancemodel_iteration",ITE,".rds"))
    results_nuisance <- rbind(results_nuisance,tmp)
  }
}
results_full <- readRDS(paste0("results/section3_1/simulation_3_1_fullmodel_iteration",1,".rds"))
for(ITE in c(2:100)){
  if(file.exists(paste0("results/section3_1/simulation_3_1_fullmodel_iteration",ITE,".rds"))){
    tmp <-readRDS(paste0("results/section3_1/simulation_3_1_fullmodel_iteration",ITE,".rds"))
    results_full <- rbind(results_full,tmp)
  }
}
results_df<-rbind(results_full, results_nuisance, results_bias, results_precision) 

results_df$estimate[results_df$type == "no_estimate_mrp_loo"] <- results_df$estimate[results_df$type == "no_estimate_mrp_loo"]^2

extract_truth <- results_df %>%
  filter(type == "mrp_error")%>%
  rename(truth = estimate)%>%
  select(-type)

results_df_noloo <- results_df %>%
  filter(type %in% c("pointwise_mrp_error","pointwise_mrp_error_estimate","no_estimate_mrp_loo"))%>%
  left_join(extract_truth)
  
ggplot(results_df_noloo, aes(x = estimate, y = truth, colour = model)) +
  geom_point() + 
  facet_wrap(.~type)+
  ggthemes::scale_colour_colorblind()+
  theme_bw()+
  xlab("Calculated Score")+ ylab("True MSE of MRP Estimate")+
  theme(legend.position = "bottom",
        legend.title = element_blank())

results_df_loo <- results_df %>%
  filter(type %in% c("no_estimate_mrp_loo","pointwise_mrp_error","loo_pointwise_mrp_error_estimate","loo_pointwise_mrp_error_estimate_alt"))%>%
  left_join(extract_truth)

ggplot(results_df_loo, aes(x = estimate, y = truth, colour = model)) +
  geom_point() + 
  coord_fixed(ratio = 1)+
  facet_wrap(.~type, ncol = 2)+
  ggthemes::scale_colour_colorblind()+
  theme_bw()+
  xlab("Calculated Score")+ ylab("True MSE of MRP Estimate")+
  theme(legend.position = "bottom",
        legend.title = element_blank())


ggsave("figure/section_squarederror_fullyobs.png", width =10, height = 6)

