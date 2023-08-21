library(tidyverse)
library(ggplot2)


results_df <- readRDS(paste0("results/section_3_2/simulation_3_2_iteration",1,".rds"))

for(ITE in c(2:100)){
  if(file.exists(paste0("results/section_3_2/simulation_3_2_iteration",ITE,".rds"))){
    tmp <-readRDS(paste0("results/section_3_2/simulation_3_2_iteration",ITE,".rds"))
    results_df <- rbind(results_df,tmp)
  }
}

results_df_crps_est <- results_df %>%
  filter(type %in% c("mrp_crps","true_crps_mrp"))%>%
  pivot_wider(names_from = c(type), values_from = estimate)
  
ggplot(results_df_crps_est, aes(x = mrp_crps, y = true_crps_mrp, colour = model)) +
  geom_point() +
  coord_fixed(ratio = 1) +
  ggthemes::scale_colour_colorblind()+
  theme_bw()+ 
  xlab("Calculated CRPS using PSIS-LOO")+ ylab("True CRPS MRP Estimate")+
  theme(legend.position = "bottom",
        legend.title = element_blank())


results_df_crps_est <- results_df %>%
  filter(type %in% c("popn_crps","true_crps_mrp"))%>%
  pivot_wider(names_from = c(type), values_from = estimate)

ggplot(results_df_crps_est, aes(x = popn_crps, y = true_crps_mrp, colour = model)) +
  geom_point() +
  coord_fixed(ratio = 1) +
  ggthemes::scale_colour_colorblind()+
  theme_bw()+ 
  xlab("Calculated CRPS using PSIS-LOO")+ ylab("True CRPS MRP Estimate")+
  theme(legend.position = "bottom",
        legend.title = element_blank())

results_df_crps_est <- results_df %>%
  filter(type %in% c("mrp_crps","true_crps_mrp"))%>%
  pivot_wider(names_from = c(type), values_from = estimate)

ggplot(results_df_crps_est, aes(x = mrp_crps, y = true_crps_mrp, colour = model)) +
  geom_point() +
  coord_fixed(ratio = 1) +
  ggthemes::scale_colour_colorblind()+
  theme_bw()+ 
  xlab("Calculated CRPS using PSIS-LOO")+ ylab("True CRPS MRP Estimate")+
  theme(legend.position = "bottom",
        legend.title = element_blank())

results_df_crps <- results_df %>%
  filter(type %in% c("true_crps_mrp","mrp_crps_truth"))%>%
  pivot_wider(names_from = c(type), values_from = estimate)

ggplot(results_df_crps, aes(x = mrp_crps_truth, y = true_crps_mrp, colour = model)) +
  geom_point() +
  coord_fixed(ratio =1)+
  ggthemes::scale_colour_colorblind()+
  theme_bw()+ 
  xlab("Calculated CRPS using PSIS-LOO")+ ylab("True CRPS MRP Estimate")+
  theme(legend.position = "bottom",
        legend.title = element_blank())

results_df_test <- results_df %>%
  filter(type %in% c("sample_crps",
                     "loo_pckg_sample_crps"))%>%
  pivot_wider(names_from = c(type), values_from = estimate)


ggplot(results_df_test, aes(x = sample_crps, y = loo_pckg_sample_crps, colour = model)) +
  geom_point() +
  coord_fixed(ratio =1)+
  ggthemes::scale_colour_colorblind()+
  theme_bw()+ 
  xlab("Calculated CRPS using PSIS-LOO")+ ylab("True CRPS MRP Estimate")+
  theme(legend.position = "bottom",
        legend.title = element_blank())
