library(tidyverse)
library(brms)
library(survey)
library(loo)
library(here)

wd <- '/mnt/lustre/projects/pMona0070/lkennedy/MRP_validation_ALT'
setwd(wd)

source("code/create_data.R")

slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
ITE = as.numeric(slurm_arrayid)

n_samp = 1000

data_use <- gen_dat(N= 20000, samp_size = n_samp, ITE = ITE)

sample <- data_use$samp_data%>%
  select(-X1_cont, -X2_cont, -X3_cont, -X4_cont)%>%
  mutate(y_obs = as.numeric(y_obs))%>%
  mutate(unique_categories = paste0(X2,X4))

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
                      backend = "cmdstanr", 
                      cores = 4)


validate_model <- function(model, population_counts,sample_truth, popn_truth){
  model_preds <-  posterior_linpred(model, transform = TRUE)
  sample_truth <- t(matrix(rep(sample_truth,4000),nrow =length(sample_truth)))
  raw_error = model_preds - sample_truth
  squared_error = (model_preds - sample_truth)^2
  
  normal_loo <- loo(model, save_psis = TRUE, reloo = TRUE)
  
  log_lik <- log_lik(model)
  psis_obj <- normal_loo$psis_object
  
  psis_var_raw_error <- E_loo(raw_error, psis_obj, type = "variance", log_ratios = -log_lik)$value
  psis_var_squared_error <- E_loo(squared_error, psis_obj, type = "variance", log_ratios = -log_lik)$value
  psis_mean_raw_error <- E_loo(raw_error, psis_obj, type = "mean", log_ratios = -log_lik)$value
  psis_mean_squared_error <- E_loo(squared_error, psis_obj, type = "mean", log_ratios = -log_lik)$value
  
  estimated_mse_mrp <- (sum(population_counts*psis_mean_raw_error)/sum(population_counts))^2
  estimated_mse_sample_sum <- sum(psis_mean_squared_error)
  estimated_mse_popn_sum <- (sum(population_counts*psis_mean_squared_error)/sum(population_counts))
  var_mse_mrp <- (sum(population_counts*psis_var_raw_error)/sum(population_counts))^2
  var_mse_sample_sum <- sum(psis_var_squared_error)
  var_mse_popn_sum <- (sum(population_counts*psis_var_squared_error)/sum(population_counts))
  
  true_error = (sum(popn_truth*population_counts)/sum(population_counts) - sum(colMeans(model_preds)*population_counts)/sum(population_counts))^2
  
  return(data.frame(estimate = c(normal_loo$estimates[1,1], estimated_mse_mrp,estimated_mse_sample_sum,estimated_mse_popn_sum, true_error), 
                    var = c(normal_loo$estimates[1,2], var_mse_mrp, var_mse_sample_sum, var_mse_popn_sum, NA),
                    type = c("elppd","mse_mrp","mse_sample","mse_population","truth")))
}

model_score_full <- validate_model(full_model_fit, popn_ps$Nj, sample_ps$y_count/sample_ps$n_j, popn_truth = popn_ps$y_count/popn_ps$Nj)

model_bias_only <- brm(y_count|trials(n_j)  ~ (1|X1) + (1|X3) + (1|X4), 
                       data = sample_ps, 
                       family = binomial(link = "logit"), 
                       backend = "cmdstanr", 
                       cores = 4)

score_bias_only <- validate_model(model_bias_only, popn_ps$Nj, sample_ps$y_count/sample_ps$n_j, popn_truth = popn_ps$y_count/popn_ps$Nj)


model_precision_only <- brm(y_count|trials(n_j)  ~ (1|X1) + (1|X2) + (1|X3), 
                            data = sample_ps, 
                            family = binomial(link = "logit"), 
                            backend = "cmdstanr", 
                            cores = 4,
                            control = list(adapt_delta = .90))


score_precision_only <- validate_model(model_precision_only, popn_ps$Nj, sample_ps$y_count/sample_ps$n_j, popn_truth = popn_ps$y_count/popn_ps$Nj)

model_nuisance_only <- brm(y_count|trials(n_j)  ~ (1|X1) + (1|X3), 
                            data = sample_ps, 
                            family = binomial(link = "logit"), 
                            backend = "cmdstanr", 
                            cores = 4)


score_nuisance_only <- validate_model(model_nuisance_only, popn_ps$Nj, sample_ps$y_count/sample_ps$n_j, popn_truth = popn_ps$y_count/popn_ps$Nj)


final_df <- rbind(data.frame(score_nuisance_only, model = "nuisance_only"),
      data.frame(score_bias_only, model = "bias_only"),
      data.frame(score_precision_only, model = "precision_only"),
      data.frame(model_score_full, model = "full_only"))

final_df$iter = ITE

saveRDS(final_df, paste0("results/section_3_1/simulation_3_1_iteration",ITE,".rds"))
