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
#implementation with CRPS & CRPS-like difference between distributions
#LOO implementation modified from the CRPS and SCRPS PR https://github.com/stan-dev/loo/pull/203/files

source("code/create_data.R")

n_samp = 1000
  print(ITE)
  data_use <- gen_dat(N= 20000, samp_size = n_samp, ITE = ITE, fully_filled =TRUE)
  
  sample <- data_use$samp_data%>%
    select(-X1_cont, -X2_cont, -X3_cont, -X4_cont)%>%
    mutate(y_obs = as.numeric(y_obs))%>%
    mutate(unique_categories = paste0(X2,X4))
  
  sample_ps <- sample %>%
    group_by(X1,X2,X3,X4)%>%
    summarise(n_j = n(), y_count = sum(y_obs), y_prob = mean(y_prob))%>%
    ungroup()%>%
    mutate(observed_samp = TRUE,
           sample_loc = 1:n())
  
  population <- data_use$popn_data
  
  popn_ps <- population %>%
    group_by(X1,X2,X3,X4)%>%
    summarise(Nj = n(),n_j = n(), y_count = sum(y_obs), y_prob = round(mean(y_prob),2))%>%
    ungroup()
  
  observed_cells <- popn_ps %>% select(X1,X2,X3,X4)%>%
    left_join(sample_ps[c("X1","X2","X3","X4", "observed_samp","sample_loc")])%>%
    mutate(observed_samp = ifelse(is.na(observed_samp),FALSE, TRUE))
  
  observed_cells_popn_loc <- observed_cells %>%
    select(observed_samp)%>%
    deframe()
  
  full_model_fit <- brm(y_count|trials(n_j) ~ (1|X1) + (1|X2) +(1|X3) + (1|X4), 
                        data = sample_ps, 
                        family = binomial(link = "logit"), 
                        backend = "cmdstanr", 
                        cores = 4)
  
  popn_predict <- posterior_linpred(full_model_fit, transform = TRUE, newdata = popn_ps) 
  popn_predict_e <- colMeans(popn_predict)
  
  validate_model <- function(model, population_ps, population_counts,sample_truth, reference_model, popn_truth,sample_counts){
    model_preds <-  posterior_linpred(model, transform = TRUE, newdata = population_ps)
    S <- nrow(model_preds)
    shuffle <- sample(1:S, size = S, replace = FALSE)
    model_preds_prime <- model_preds[shuffle,]
  
    focus_model_loo <- loo(model, save_psis = TRUE, reloo = TRUE)
    
    log_lik <- log_lik(model)
    log_lik_prime <- log_lik[shuffle,]
    psis_obj_focus_model <- focus_model_loo$psis_object
    psis_obj_focus_prime_joint <- psis(-log_lik-log_lik_prime)
    
    #regular crps
    EXX_abs = E_loo(abs(model_preds - model_preds_prime), psis_obj_focus_prime_joint, log_ratios = -log_lik - log_lik_prime)
    EXY_abs = E_loo(abs(sweep(model_preds,2,sample_truth)), psis_obj_focus_model, log_ratios = -log_lik)
  
    #alternate crps
    EXX_alt = E_loo(model_preds - model_preds_prime, psis_obj_focus_prime_joint, log_ratios = -log_lik - log_lik_prime)
    EXY_alt = E_loo(sweep(model_preds,2,sample_truth), psis_obj_focus_model, log_ratios = -log_lik)
    EXY_alt_truth = E_loo(sweep(model_preds,2,popn_truth), psis_obj_focus_model, log_ratios = -log_lik)
    
    sample_crps = sum(.5*sample_counts*EXX_abs$value-sample_counts*EXY_abs$value)/sum(sample_counts)
    popn_crps = sum(.5*population_counts*EXX_abs$value-population_counts*EXY_abs$value)/sum(population_counts)
    mrp_crps = (abs(sum(.5*population_counts*EXX_alt$value))-abs(sum(population_counts*EXY_alt$value)))/sum(population_counts)
    mrp_crps_truth = (abs(sum(.5*population_counts*EXX_alt$value))-abs(sum(population_counts*EXY_alt_truth$value)))/sum(population_counts)
    
    
    loo_pckg_sample_crps  = sum(sample_counts*loo_crps(model_preds, model_preds_prime, y = sample_truth, log_lik  = log_lik)$pointwise)/sum(sample_counts)
    loo_pckg_popn_crps  = sum(population_counts*loo_crps(model_preds, model_preds_prime, y = sample_truth, log_lik  = log_lik)$pointwise)/sum(population_counts)
    
    true_crps_mrp = crps(colSums(apply(model_preds,1,function(x)x*population_counts))/sum(population_counts), 
                         colSums(apply(model_preds_prime,1,function(x)x*population_counts))/sum(population_counts),
                         sum(popn_truth*population_counts)/sum(population_counts))
    
    return(data.frame(estimate = c(sample_crps, popn_crps, mrp_crps,mrp_crps_truth, loo_pckg_sample_crps,loo_pckg_popn_crps, true_crps_mrp$estimates[[1]]), 
                      type = c("sample_crps","popn_crps","mrp_crps","mrp_crps_truth","loo_pckg_sample_crps","loo_pckg_popn_crps","true_crps_mrp")))
  }
  
  model_score_full <- validate_model(full_model_fit, popn_ps, popn_ps$Nj, sample_truth = sample_ps$y_count/sample_ps$n_j, popn_truth = popn_ps$y_count/popn_ps$Nj,observed_cells_popn_loc,sample_counts = sample_ps$n_j)
  
  model_bias_only <- brm(y_count|trials(n_j)  ~ (1|X1) + (1|X3) + (1|X4), 
                         data = sample_ps, 
                         family = binomial(link = "logit"), 
                         backend = "cmdstanr", 
                         cores = 4)
  
  score_bias_only <- validate_model(model_bias_only, popn_ps, popn_ps$Nj, sample_truth = sample_ps$y_count/sample_ps$n_j, popn_truth = popn_ps$y_count/popn_ps$Nj,observed_cells_popn_loc,sample_counts = sample_ps$n_j)
  
  
  model_precision_only <- brm(y_count|trials(n_j)  ~ (1|X1) + (1|X2) + (1|X3), 
                              data = sample_ps, 
                              family = binomial(link = "logit"), 
                              backend = "cmdstanr", 
                              cores = 4,
                              control = list(adapt_delta = .90))
  
  
  score_precision_only <- validate_model(model_precision_only, popn_ps, popn_ps$Nj, sample_truth = sample_ps$y_count/sample_ps$n_j, popn_truth = popn_ps$y_count/popn_ps$Nj,observed_cells_popn_loc,sample_counts = sample_ps$n_j)
  
  model_nuisance_only <- brm(y_count|trials(n_j)  ~ (1|X1) + (1|X3), 
                             data = sample_ps, 
                             family = binomial(link = "logit"), 
                             backend = "cmdstanr", 
                             cores = 4)
  
  
  score_nuisance_only <- validate_model(model_nuisance_only, popn_ps, popn_ps$Nj, sample_truth = sample_ps$y_count/sample_ps$n_j, popn_truth = popn_ps$y_count/popn_ps$Nj,observed_cells_popn_loc,sample_counts = sample_ps$n_j)
  
  final_df <- rbind(data.frame(score_nuisance_only, model = "nuisance_only"),
                    data.frame(score_bias_only, model = "bias_only"),
                    data.frame(score_precision_only, model = "precision_only"),
                    data.frame(model_score_full, model = "full_only"))
  
  final_df$iter = ITE
  
  saveRDS(final_df, paste0("results/section_3_2/simulation_3_2_iteration",ITE,".rds"))

