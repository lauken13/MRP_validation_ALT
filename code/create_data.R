## Create data for bias/precision variables
## Taken from https://github.com/lauken13/LOO_MRP 26/1/2023

library(tidyverse)
library(brms) # need cmdstanr to use with brms
library(loo) # calculating loo and elpd
library(survey) # creating raked weights 
library(tidybayes)


## generating data - 5 continuous predictors/covariates and a binary outcome 
gen_dat <- function(N, samp_size, ITE){
  set.seed(65438)
  
  pn = 100 # number of different population
  seed = round(runif(pn, min=10, max=100000),0) # fixed seed number
  
  # setting seed using array ID
  set.seed(seed[ITE])
  popn_data <- data.frame(X1_cont = rnorm(N, 0, 2), 
                          X2_cont = rnorm(N, 0, 2),
                          X3_cont = rnorm(N, 0, 2), 
                          X4_cont = rnorm(N, 0, 2))
  
  
  wkly1 = 0.1
  strg1 = 1
  
  ## generating continuous and binary outcome
  popn_data$y_prob <- inv_logit_scaled(wkly1*popn_data$X1_cont +
                                         strg1*popn_data$X2_cont +
                                         wkly1*popn_data$X3_cont +
                                         strg1*popn_data$X4_cont)
  popn_data$y_obs <- rbinom(N,1,popn_data$y_prob)
  
  
  ## generate inclusion prob. for each individual
  # weakly predictive - 0.1 (sd), strongly predictive - 1 (sd)
  wkly2 = 0.1
  strg2 = 1
  popn_data$incl_prob <- inv_logit_scaled(wkly2*popn_data$X1_cont + 
                                            wkly2*popn_data$X2_cont + 
                                            strg2*popn_data$X3_cont +
                                            strg2*popn_data$X4_cont)
  
  
  ## categorising the continuous covariates 
  J = 5
  
  popn_data <- popn_data %>% 
    mutate(X1_fct = cut_interval(X1_cont,J),
           X2_fct = cut_interval(X2_cont,J),
           X3_fct = cut_interval(X3_cont,J),
           X4_fct = cut_interval(X4_cont,J)) %>% 
    mutate(across(X1_fct:X4_fct, ~ as.numeric(.x))) 
  
  popn_data <- popn_data %>%
    rename(X1 = X1_fct,
           X2 = X2_fct,
           X3 = X3_fct,
           X4 = X4_fct)
  
  ## generating samples
  samp_loc = sample(1:nrow(popn_data), size = samp_size-(J*4), replace=F, prob = popn_data$incl_prob)
  
  ## making sure at least each level of the covariates are sampled
  for(j in 1:J){
    samp_loc[length(samp_loc)+1] = sample(which(popn_data$X1 == j), size=1)
    samp_loc[length(samp_loc)+1] = sample(which(popn_data$X2 == j), size=1)
    samp_loc[length(samp_loc)+1] = sample(which(popn_data$X3 == j), size=1)
    samp_loc[length(samp_loc)+1] = sample(which(popn_data$X4 == j), size=1)
  }
  
  samp_data = popn_data[samp_loc,]
  
  ## creating survey design
  svy1 = svydesign(ids=~1, # cluster id, ~1 for no clusters
                   weights=~rep(1,nrow(samp_data)), # equal weights for each unit
                   data=samp_data)
  
  ## calculating population totals for each level
  X1_margin = xtabs(~X1, data=popn_data)
  X2_margin = xtabs(~X2, data=popn_data)
  X3_margin = xtabs(~X3, data=popn_data)
  X4_margin = xtabs(~X4, data=popn_data)
  
  ## raked to the population
  rake1 = rake(design = svy1, sample.margins = list(~X1,~X2,~X3,~X4), 
               population.margins = list(X1_margin, X2_margin, X3_margin, X4_margin))
  
  ## raked weights ####
  samp_data$wts = weights(rake1)
  
  ## factorise the relevant variables
  str(samp_data)
  ind = c(1:4,6)
  samp_data[,ind] = apply(samp_data[,ind], 2, function(x)as.factor(x))
  
  ## make poststrat table for popn
  popn_ps = popn_data %>% 
    group_by(X1, X2, X3, X4) %>% 
    summarise(Nj = n(), prob_out = mean(y_prob), .groups='drop') %>% 
    ungroup()
  
  all_list <- list(samp_data, popn_ps, popn_data, N, J)
  names(all_list) = c('samp_data', 'popn_ps', 'popn_data', 'N', 'J')
  all_list
}
