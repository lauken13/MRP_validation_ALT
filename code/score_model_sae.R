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
library(posterior)

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

#full model
full_model_fit <- brm(y_count|trials(n_j) ~ (1|X1) + (1|X2) +(1|X3) + (1|X4), 
                      data = sample_ps, 
                      family = binomial(link = "logit"), 
                      backend = "rstan", 
                      cores = 1,
                      save_pars = save_pars(all = TRUE))

approx_loco_sae_score_x1 <- approx_loco_sae_score(model = full_model_fit,
                                                small_area_var = "X1",
                                                popn_ps = popn_ps,
                                                sample_ps = sample_ps)

approx_loco_sae_score_x2 <- approx_loco_sae_score(model = full_model_fit,
                                                  small_area_var = "X2",
                                                  popn_ps = popn_ps,
                                                  sample_ps = sample_ps)

approx_loco_sae_score_x3 <- approx_loco_sae_score(model = full_model_fit,
                                                  small_area_var = "X3",
                                                  popn_ps = popn_ps,
                                                  sample_ps = sample_ps)

approx_loco_sae_score_x4 <- approx_loco_sae_score(model = full_model_fit,
                                                  small_area_var = "X4",
                                                  popn_ps = popn_ps,
                                                  sample_ps = sample_ps)

final_df_fullmodel <- rbind(approx_loco_sae_score_x1,approx_loco_sae_score_x2,approx_loco_sae_score_x3,approx_loco_sae_score_x4)

precision_model_fit <- brm(y_count|trials(n_j) ~ (1|X1) + (1|X2) +(1|X3), 
                      data = sample_ps, 
                      family = binomial(link = "logit"), 
                      backend = "rstan", 
                      cores = 1,
                      save_pars = save_pars(all = TRUE))

precision_approx_loco_sae_score_x1 <- approx_loco_sae_score(model = precision_model_fit,
                                                  small_area_var = "X1",
                                                  popn_ps = popn_ps,
                                                  sample_ps = sample_ps)

precision_approx_loco_sae_score_x2 <- approx_loco_sae_score(model = precision_model_fit,
                                                  small_area_var = "X2",
                                                  popn_ps = popn_ps,
                                                  sample_ps = sample_ps)

precision_approx_loco_sae_score_x3 <- approx_loco_sae_score(model = precision_model_fit,
                                                  small_area_var = "X3",
                                                  popn_ps = popn_ps,
                                                  sample_ps = sample_ps)

precision_approx_loco_sae_score_x4 <- approx_loco_sae_score(model = precision_model_fit,
                                                  small_area_var = "X4",
                                                  popn_ps = popn_ps,
                                                  sample_ps = sample_ps)

final_df_precisionmodel <- rbind(precision_approx_loco_sae_score_x1,precision_approx_loco_sae_score_x2,precision_approx_loco_sae_score_x3,precision_approx_loco_sae_score_x4)

bias_model_fit <- brm(y_count|trials(n_j) ~ (1|X1) + (1|X4) +(1|X3), 
                      data = sample_ps, 
                      family = binomial(link = "logit"), 
                      backend = "rstan", 
                      cores = 1,
                      save_pars = save_pars(all = TRUE))

bias_approx_loco_sae_score_x1 <- approx_loco_sae_score(model = bias_model_fit,
                                                       small_area_var = "X1",
                                                       popn_ps = popn_ps,
                                                       sample_ps = sample_ps)

bias_approx_loco_sae_score_x2 <- approx_loco_sae_score(model = bias_model_fit,
                                                       small_area_var = "X2",
                                                       popn_ps = popn_ps,
                                                       sample_ps = sample_ps)

bias_approx_loco_sae_score_x3 <- approx_loco_sae_score(model = bias_model_fit,
                                                       small_area_var = "X3",
                                                       popn_ps = popn_ps,
                                                       sample_ps = sample_ps)

bias_approx_loco_sae_score_x4 <- approx_loco_sae_score(model = bias_model_fit,
                                                       small_area_var = "X4",
                                                       popn_ps = popn_ps,
                                                       sample_ps = sample_ps)

final_df_biasmodel <- rbind(bias_approx_loco_sae_score_x1,bias_approx_loco_sae_score_x2,bias_approx_loco_sae_score_x3,bias_approx_loco_sae_score_x4)

nuisance_model_fit <- brm(y_count|trials(n_j) ~ (1|X1) + (1|X3), 
                           data = sample_ps, 
                           family = binomial(link = "logit"), 
                           backend = "rstan", 
                           cores = 1,
                           save_pars = save_pars(all = TRUE))

nuisance_approx_loco_sae_score_x1 <- approx_loco_sae_score(model = nuisance_model_fit,
                                                            small_area_var = "X1",
                                                            popn_ps = popn_ps,
                                                            sample_ps = sample_ps)

nuisance_approx_loco_sae_score_x2 <- approx_loco_sae_score(model = nuisance_model_fit,
                                                            small_area_var = "X2",
                                                            popn_ps = popn_ps,
                                                            sample_ps = sample_ps)

nuisance_approx_loco_sae_score_x3 <- approx_loco_sae_score(model = nuisance_model_fit,
                                                            small_area_var = "X3",
                                                            popn_ps = popn_ps,
                                                            sample_ps = sample_ps)

nuisance_approx_loco_sae_score_x4 <- approx_loco_sae_score(model = nuisance_model_fit,
                                                            small_area_var = "X4",
                                                            popn_ps = popn_ps,
                                                            sample_ps = sample_ps)


final_df_nuisancemodel <- rbind(nuisance_approx_loco_sae_score_x1,nuisance_approx_loco_sae_score_x2,nuisance_approx_loco_sae_score_x3,nuisance_approx_loco_sae_score_x4)


nuisance_x1_model_fit <- brm(y_count|trials(n_j) ~ (1|X1), 
                          data = sample_ps, 
                          family = binomial(link = "logit"), 
                          backend = "rstan", 
                          cores = 1,
                          save_pars = save_pars(all = TRUE))

nuisance_x1_approx_loco_sae_score_x1 <- approx_loco_sae_score(model = nuisance_x1_model_fit,
                                                           small_area_var = "X1",
                                                           popn_ps = popn_ps,
                                                           sample_ps = sample_ps)

nuisance_x1_approx_loco_sae_score_x2 <- approx_loco_sae_score(model = nuisance_x1_model_fit,
                                                           small_area_var = "X2",
                                                           popn_ps = popn_ps,
                                                           sample_ps = sample_ps)

nuisance_x1_approx_loco_sae_score_x3 <- approx_loco_sae_score(model = nuisance_x1_model_fit,
                                                           small_area_var = "X3",
                                                           popn_ps = popn_ps,
                                                           sample_ps = sample_ps)

nuisance_x1_approx_loco_sae_score_x4 <- approx_loco_sae_score(model = nuisance_x1_model_fit,
                                                           small_area_var = "X4",
                                                           popn_ps = popn_ps,
                                                           sample_ps = sample_ps)


final_df_nuisance_x1model <- rbind(nuisance_x1_approx_loco_sae_score_x1,nuisance_x1_approx_loco_sae_score_x2,nuisance_x1_approx_loco_sae_score_x3,nuisance_x1_approx_loco_sae_score_x4)

nuisance_x3_model_fit <- brm(y_count|trials(n_j) ~ (1|X3), 
                             data = sample_ps, 
                             family = binomial(link = "logit"), 
                             backend = "rstan", 
                             cores = 1,
                             save_pars = save_pars(all = TRUE))

nuisance_x3_approx_loco_sae_score_x1 <- approx_loco_sae_score(model = nuisance_x3_model_fit,
                                                              small_area_var = "X1",
                                                              popn_ps = popn_ps,
                                                              sample_ps = sample_ps)

nuisance_x3_approx_loco_sae_score_x2 <- approx_loco_sae_score(model = nuisance_x3_model_fit,
                                                              small_area_var = "X2",
                                                              popn_ps = popn_ps,
                                                              sample_ps = sample_ps)

nuisance_x3_approx_loco_sae_score_x3 <- approx_loco_sae_score(model = nuisance_x3_model_fit,
                                                              small_area_var = "X3",
                                                              popn_ps = popn_ps,
                                                              sample_ps = sample_ps)

nuisance_x3_approx_loco_sae_score_x4 <- approx_loco_sae_score(model = nuisance_x3_model_fit,
                                                              small_area_var = "X4",
                                                              popn_ps = popn_ps,
                                                              sample_ps = sample_ps)


final_df_nuisance_x3model <- rbind(nuisance_x3_approx_loco_sae_score_x1,nuisance_x3_approx_loco_sae_score_x2,nuisance_x3_approx_loco_sae_score_x3,nuisance_x3_approx_loco_sae_score_x4)



final_df <- rbind(final_df_fullmodel,final_df_biasmodel, final_df_precisionmodel, final_df_nuisancemodel,final_df_nuisance_x1model,final_df_nuisance_x3model)

final_df$iter = ITE

saveRDS(final_df, paste0("results/sae_scores/scores_sae_approx",ITE,".rds"))
