on_cluster = TRUE

if(on_cluster){
  wd <- '/gpfs/users/a1193023/MRP_validation_ALT'
  setwd(wd)
  slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
  ITE = as.numeric(slurm_arrayid)
  print(ITE)
  .libPaths(c("/gpfs/users/a1193023/local/RLibs",.libPaths()))
}else{
  #For testing
  ITE = 1
}

library(tidyverse)
library(brms)
library(survey)
library(loo)
library(here)

source("code/create_data.R")
source("code/error_helper_functions.R")

n_samp = 1000

data_use <- gen_dat(N= 20000, samp_size = n_samp, ITE = ITE)

sample <- data_use$samp_data%>%
  select(-X1_cont, -X2_cont, -X3_cont, -X4_cont)%>%
  mutate(y_obs = as.numeric(y_obs))%>%
  mutate(unique_categories = paste0(X2,X4))

summary(sample)

set.seed(65438)

pn = 100 # number of different population
seed = round(runif(pn, min=10, max=100000),0) # fixed seed number

# setting seed using array ID
set.seed(seed[ITE])

sample_ps <- sample %>%
  group_by(X1,X2,X3,X4)%>%
  summarise(n_j = n(), y_count = sum(y_obs), y_prob = mean(y_prob))%>%
  ungroup()

population <- data_use$popn_data

popn_ps <- population %>%
  group_by(X1,X2,X3,X4)%>%
  summarise(Nj = n(), y_count = sum(y_obs), y_prob = round(mean(y_prob),2))%>%
  ungroup()


full_model_fit <- brm(y_count|trials(n_j) ~ (1|X1) + (1|X2) +(1|X3) + (1|X4), 
                      data = sample_ps, 
                      family = binomial(link = "logit"), 
                      backend = "rstan", 
                      cores = 1,
                      save_pars = save_pars(all = TRUE))

validate_score(model = full_model_fit,
           popn_counts = popn_ps$Nj,
           popn_obs = popn_ps$y_count)

model_score_full <- validate_model(model = full_model_fit,
                                   popn_counts = popn_ps$Nj,
                                   sample_counts =sample_ps$n_j,
                                   sample_obs = sample_ps$y_count,
                                   popn_obs = popn_ps$y_count)

model_score_full_manual <- validate_manually(model = full_model_fit,
                                             popn_counts = popn_ps$Nj,
                                             sample_counts =sample_ps$n_j,
                                             sample_obs = sample_ps$y_count,
                                             popn_obs = popn_ps$y_count)
model_score_full <- rbind(model_score_full, model_score_full_manual)
gc()

final_df <- data.frame(model_score_full, model = "full_only")

final_df$iter = ITE

saveRDS(final_df, paste0("results/section3_1/simulation_3_1_fullmodel_iteration",ITE,".rds"))
