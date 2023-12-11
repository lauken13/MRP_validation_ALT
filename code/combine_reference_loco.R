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

data_use <- gen_dat(N= 20000, samp_size = n_samp, ITE = ITE, fully_filled = FALSE)

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

sample_summary <- data.frame(prop_cells_covered = nrow(sample_ps)/nrow(popn_ps), iter = ITE)

saveRDS(sample_summary, paste0("results/combined_reference_model/sample_summary",ITE,".rds"))

reference_model_fit <- brm(y_count|trials(n_j) ~ (1|X1) + (1|X2) +(1|X3) + (1|X4), 
                           data = sample_ps, 
                           family = binomial(link = "logit"), 
                           backend = "rstan", 
                           cores = 1,
                           save_pars = save_pars(all = TRUE))


candidate_model2_fit <- brm(y_count|trials(n_j) ~ (1|X1)  +(1|X3) + (1|X4), 
                            data = sample_ps, 
                            family = binomial(link = "logit"), 
                            backend = "rstan", 
                            cores = 1,
                            save_pars = save_pars(all = TRUE))


candidate_model3_fit <- brm(y_count|trials(n_j) ~ (1|X1)  +(1|X2) + (1|X3), 
                            data = sample_ps, 
                            family = binomial(link = "logit"), 
                            backend = "rstan", 
                            cores = 1,
                            save_pars = save_pars(all = TRUE))



candidate_model4_fit <- brm(y_count|trials(n_j) ~ (1|X1)  + (1|X3), 
                            data = sample_ps, 
                            family = binomial(link = "logit"), 
                            backend = "rstan", 
                            cores = 1,
                            save_pars = save_pars(all = TRUE))

combined_reference_model2 <- partially_obs_approx_loco_score(reference_model = reference_model_fit,
                                                             candidate_model = candidate_model2_fit,
                                                             popn_ps = popn_ps,
                                                             sample_ps = sample_ps)

combined_reference_model3 <- partially_obs_approx_loco_score(reference_model = reference_model_fit,
                                                             candidate_model = candidate_model3_fit,
                                                             popn_ps = popn_ps,
                                                             sample_ps = sample_ps)

combined_reference_model4 <- partially_obs_approx_loco_score(reference_model = reference_model_fit,
                                                             candidate_model = candidate_model4_fit,
                                                             popn_ps = popn_ps,
                                                             sample_ps = sample_ps)

combined_incorrectreference_model1 <- partially_obs_approx_loco_score(reference_model = candidate_model3_fit,
                                                             candidate_model = reference_model_fit,
                                                             popn_ps = popn_ps,
                                                             sample_ps = sample_ps)

combined_incorrectreference_model2 <- partially_obs_approx_loco_score(reference_model = candidate_model3_fit,
                                                             candidate_model = candidate_model2_fit,
                                                             popn_ps = popn_ps,
                                                             sample_ps = sample_ps)

combined_incorrectreference_model4 <- partially_obs_approx_loco_score(reference_model = candidate_model3_fit,
                                                             candidate_model = candidate_model4_fit,
                                                             popn_ps = popn_ps,
                                                             sample_ps = sample_ps)

true_score_model2 <- population_score(candidate_model2_fit,
           popn_counts = popn_ps$Nj,
           popn_obs = popn_ps$y_count,
           popn_ps = popn_ps %>%
             mutate(n_j = Nj)) 

true_score_model3 <- population_score(candidate_model3_fit,
                                      popn_counts = popn_ps$Nj,
                                      popn_obs = popn_ps$y_count,
                                      popn_ps = popn_ps %>%
                                        mutate(n_j = Nj)) 

true_score_model4 <- population_score(candidate_model4_fit,
                                      popn_counts = popn_ps$Nj,
                                      popn_obs = popn_ps$y_count,
                                      popn_ps = popn_ps %>%
                                        mutate(n_j = Nj)) 

final_df <- rbind(rbind(combined_reference_model2, combined_reference_model3,combined_reference_model4) %>%
                    mutate(reference_model = "FULL"),
                  rbind(combined_incorrectreference_model1,combined_incorrectreference_model2,combined_incorrectreference_model4) %>%
                    mutate(reference_model = "PRECISION ONLY"),
                  rbind(true_score_model2,true_score_model3,true_score_model4) %>%
                    mutate(reference_model = NA))


final_df$iter = ITE

saveRDS(final_df, paste0("results/combined_reference_model/scores_reference_results_approx",ITE,".rds"))

