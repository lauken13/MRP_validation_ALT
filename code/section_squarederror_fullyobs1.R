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


validate_model <-
  function(model,
           popn_counts,
           sample_counts,
           sample_obs,
           popn_obs) {
    sample_truth = sample_obs/sample_counts
    popn_truth = popn_obs/popn_counts
    
    model_preds <-  posterior_linpred(model, transform = TRUE)
    print(summary(colMeans(model_preds)))
    sample_truth_matrix <-
      t(matrix(rep(sample_truth, 4000), nrow = length(sample_truth)))
    popn_truth_matrix <-
      t(matrix(rep(popn_truth, 4000), nrow = length(popn_truth)))
    raw_error = model_preds - sample_truth_matrix
    squared_error = (model_preds - sample_truth_matrix) ^ 2
    true_error = model_preds - popn_truth_matrix
    squared_true_error = model_preds - popn_truth_matrix
  
    normal_loo <- loo(model, save_psis = TRUE, reloo = TRUE, cores = 1)
    print(normal_loo)
    
    log_lik <- log_lik(model, moment_match = TRUE)
    print(summary(colMeans(log_lik)))
    psis_obj <- normal_loo$psis_object
    print(summary(psis_obj$diagnostics$pareto_k))
    print(summary(psis_obj$log_weights[,1]))
    print(summary(colMeans(psis_obj$log_weights)))
    
    psis_mean_raw_error <-
      E_loo(raw_error,
            psis_obj,
            type = "mean",
            log_ratios = -log_lik)$value
    psis_mean_squared_raw_error <-
      E_loo(squared_error,
            psis_obj,
            type = "mean",
            log_ratios = -log_lik)$value
    psis_mean_true_error <-
      E_loo(true_error,
            psis_obj,
            type = "mean",
            log_ratios = -log_lik)$value
    psis_mean_squared_true_error <-
      E_loo(squared_true_error,
            psis_obj,
            type = "mean",
            log_ratios = -log_lik)$value
    
    N = sum(popn_counts)
    Nj = popn_counts
    J = ncol(sample_truth_matrix)
    pointwise_mrp_error_estimate <-
      Pointwise_MRP_error_calc(Nj, N, sample_truth, colMeans(model_preds))
    pointwise_mrp_error <-
      Pointwise_MRP_error_calc(Nj, N, popn_truth, colMeans(model_preds))
    mrp_error <-
      True_MRP_error_calc(Nj, N, popn_truth, colMeans(model_preds))
    
    print(mrp_error)
    
    loo_pointwise_mrp_error <- ((sum(Nj * psis_mean_true_error) / N))^2
    print(loo_pointwise_mrp_error)
    
    loo_pointwise_mrp_error_estimate <-
      (sum(Nj * psis_mean_raw_error) / N)^2 
    print(loo_pointwise_mrp_error_estimate)
    
    
    ###Alternate score focussing on misclassification
    tmp_error_var = rep(0,max(sample_counts))
    for(k in 1:max(sample_counts)){
      sample_truth_alt <- sample_truth
      more_than_one_obs <- which(sample_counts>1)
      for (i in more_than_one_obs){
        sample_truth_alt[i] <- sample(c(rep(1,sample_obs[i]),rep(0,sample_counts[i]-sample_obs[i])),1)
      }
      sample_raw_alt_matrix <-
        t(matrix(rep(sample_truth_alt, 4000), nrow = length(popn_truth)))
      raw_error_alt = model_preds - sample_raw_alt_matrix
      psis_mean_raw_error_alt <-
        E_loo(raw_error_alt,
              psis_obj,
              type = "mean",
              log_ratios = -log_lik)$value
      
      tmp_error_var[k] <-
        (sum(Nj * psis_mean_raw_error_alt) / N) 
    }
    loo_pointwise_mrp_error_estimate_alt <- mean((tmp_error_var)^2)
    print(loo_pointwise_mrp_error_estimate_alt)
    return(data.frame(
      estimate = c(
        normal_loo$estimates[1, 1],
        pointwise_mrp_error_estimate,
        pointwise_mrp_error,
        mrp_error,
        loo_pointwise_mrp_error_estimate,
        loo_pointwise_mrp_error_estimate_alt,
        loo_pointwise_mrp_error
      ),
      type = c(
        "elppd",
        "pointwise_mrp_error_estimate",
        "pointwise_mrp_error",
        "mrp_error",
        "loo_pointwise_mrp_error_estimate",
        "loo_pointwise_mrp_error_estimate_alt",
        "loo_pointwise_mrp_error"
      )
    ))
  }

validate_manually <- function(model, popn_counts, sample_counts, sample_obs, popn_obs){
  sample_truth = sample_obs/sample_counts
  popn_truth = popn_obs/popn_counts
  n_cells = length(sample_counts)
  cell_median_error <- rep(NA,n_cells)
  for(i in 1:n_cells){
    print(i/n_cells)
    test_data <- model$data[i,]
    train_data <- model$data[(1:n_cells)[1:n_cells!=i],]
    refit_model <- brm(model$formula,
                       data = train_data,
                       family = binomial(link = "logit"), 
                       backend = "rstan", 
                       cores = 1,
                       save_pars = save_pars(all = TRUE))
    loo_pred_one <- posterior_linpred(refit_model, newdata = test_data, transform = TRUE)
    cell_median_error[i] <- median(loo_pred_one) - sample_truth[i]
    rm(refit_model)
    gc()
  }
  mrp_error_est_loo <- (sum(popn_counts*cell_median_error)/sum(popn_counts))^2
  return(data.frame(estimate = mrp_error_est_loo, type = "no_estimate_mrp_loo"))
}

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

# model_bias_only <- brm(y_count|trials(n_j)  ~ (1|X1) + (1|X3) + (1|X4), 
#                        data = sample_ps, 
#                        family = binomial(link = "logit"), 
#                        backend = "rstan", 
#                        cores = 1,
#                        save_pars = save_pars(all = TRUE))
# 
# score_bias_only <- validate_model(model_bias_only, 
#                                   popn_counts = popn_ps$Nj, 
#                                   sample_counts =sample_ps$n_j, 
#                                   sample_obs = sample_ps$y_count, 
#                                   popn_obs = popn_ps$y_count)
# 
# score_bias_only_manual <- validate_manually(model = model_bias_only, 
#                                              popn_counts = popn_ps$Nj, 
#                                              sample_counts =sample_ps$n_j, 
#                                              sample_obs = sample_ps$y_count, 
#                                              popn_obs = popn_ps$y_count)
# score_bias_only <- rbind(score_bias_only, score_bias_only_manual)
# 
# 
# final_df <- data.frame(score_bias_only, model = "bias_only")
# 
# final_df$iter = ITE
# 
# saveRDS(final_df, paste0("results/section3_1/simulation_3_1_biasmodel_iteration",ITE,".rds"))

# 
# model_precision_only <- brm(y_count|trials(n_j)  ~ (1|X1) + (1|X2) + (1|X3), 
#                             data = sample_ps, 
#                             family = binomial(link = "logit"), 
#                             backend = "rstan", 
#                             cores = 1,
#                             control = list(adapt_delta = .90),
#                             save_pars = save_pars(all = TRUE))
# 
# 
# score_precision_only <- validate_model(model_precision_only, 
#                                        popn_counts = popn_ps$Nj, 
#                                        sample_counts =sample_ps$n_j, 
#                                        sample_obs = sample_ps$y_count, 
#                                        popn_obs = popn_ps$y_count)
# 
# score_precision_only_manual <- validate_manually(model = model_precision_only, 
#                                                  popn_counts = popn_ps$Nj, 
#                                                  sample_counts =sample_ps$n_j, 
#                                                  sample_obs = sample_ps$y_count, 
#                                                  popn_obs = popn_ps$y_count)
# score_precision_only <- rbind(score_precision_only, score_precision_only_manual)
# 
# 
# final_df <- data.frame(score_precision_only, model = "precision_only")
# 
# final_df$iter = ITE
# 
# saveRDS(final_df, paste0("results/section3_1/simulation_3_1_precisionmodel_iteration",ITE,".rds"))
# 
# model_nuisance_only <- brm(y_count|trials(n_j)  ~ (1|X1) + (1|X3), 
#                             data = sample_ps, 
#                             family = binomial(link = "logit"), 
#                             backend = "rstan", 
#                             cores = 1,
#                            save_pars = save_pars(all = TRUE))
# 
# score_nuisance_only <- validate_model(model_nuisance_only, 
#                                       popn_counts = popn_ps$Nj, 
#                                       sample_counts =sample_ps$n_j, 
#                                       sample_obs = sample_ps$y_count, 
#                                       popn_obs = popn_ps$y_count)
# 
# 
# score_nuisance_only_manual <- validate_manually(model = model_nuisance_only, 
#                                                 popn_counts = popn_ps$Nj, 
#                                                 sample_counts =sample_ps$n_j, 
#                                                 sample_obs = sample_ps$y_count, 
#                                                 popn_obs = popn_ps$y_count)
# score_nuisance_only <- rbind(score_nuisance_only, score_nuisance_only_manual)
# 
# 
# 
# final_df <- data.frame(score_nuisance_only, model = "nuisance_only")
# 
# final_df$iter = ITE
# 
# saveRDS(final_df, paste0("results/section3_1/simulation_3_1_nuisancemodel_iteration",ITE,".rds"))
# 
# final_df <- rbind(data.frame(score_nuisance_only, model = "nuisance_only"),
#       data.frame(score_bias_only, model = "bias_only"),
#       data.frame(score_precision_only, model = "precision_only"),
#       data.frame(model_score_full, model = "full_only"))
# 
# final_df$iter = ITE
# 
# saveRDS(final_df, paste0("results/section3_1/simulation_3_1_iteration",ITE,".rds"))
