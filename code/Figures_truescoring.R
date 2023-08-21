library(tidyverse)
library(ggplot2)

### Read in CRPS scores
crps_results_df <- readRDS(paste0("results/section_3_2/simulation_3_2_iteration",1,".rds"))

for(ITE in c(2:100)){
  if(file.exists(paste0("results/section_3_2/simulation_3_2_iteration",ITE,".rds"))){
    tmp <-readRDS(paste0("results/section_3_2/simulation_3_2_iteration",ITE,".rds"))
    crps_results_df <- rbind(crps_results_df,tmp)
  }
}

results_df_crps_est <- crps_results_df %>%
  filter(type %in% c("mrp_crps_truth","true_crps_mrp","loo_pckg_popn_crps"))


ggplot(results_df_crps_est, aes(x = mrp_crps, y = true_crps_mrp, colour = model)) +
  geom_point() +
  coord_fixed(ratio = 1) +
  ggthemes::scale_colour_colorblind()+
  theme_bw()+ 
  xlab("Calculated CRPS using PSIS-LOO")+ ylab("True CRPS MRP Estimate")+
  theme(legend.position = "bottom",
        legend.title = element_blank())
