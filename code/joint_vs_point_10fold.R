source("code/create_data.R")

n_samp = 2000
for(ITE in 1:10){
  print(ITE)
  data_use <- gen_dat(N= 20000, samp_size = n_samp, ITE = ITE)
  
  sample <- data_use$samp_data%>%
    select(-X1_cont, -X2_cont, -X3_cont, -X4_cont)%>%
    mutate(y_obs = as.numeric(y_obs))%>%
    mutate(unique_categories = paste0(X2,X4))
  
  population <- data_use$popn_data
  
  full_model_fit <- brm(y_obs ~ (1|X1) + (1|X2) +(1|X3) + (1|X4), 
                        data = sample, 
                        family = binomial(link = "logit"), 
                        backend = "cmdstanr", 
                        cores = 4)
  
  model_bias_only <- brm(y_obs  ~ (1|X1) + (1|X3) + (1|X4), 
                         data = sample, 
                         family = binomial(link = "logit"), 
                         backend = "cmdstanr", 
                         cores = 4)
  
  
  model_precision_only <- brm(y_obs  ~ (1|X1) + (1|X2) + (1|X3), 
                         data = sample, 
                         family = binomial(link = "logit"), 
                         backend = "cmdstanr", 
                         cores = 4)
  
  
  calculate_model_validation <- function(K = 10, model, model_name, sample, sample_size){
    # Use random fold as some combinations have 1 observations
    ids <- kfold_split_random(K = 10, N = sample_size)
    
    
    loo_model <- loo(model)
    ll<-log_lik(model)
    
    cvii <- rep(1:K,each=(ceiling(sample_size/K)))
    ll2 <- matrix(nrow=nrow(ll),ncol=K)
    
    #joint log likelihood
    for (i in 1:K) { ll2[,i] <- apply(ll[,cvii==i],1,sum) }
    loo_foldK_joint<-loo(ll2, save_psis=TRUE)
    
    #pointwise log likelihood for folds
    ll_weights <- loo_foldK_joint$psis_object$log_weights
    w <- matrix(nrow=nrow(ll),ncol=n_samp)
    for (i in 1:K) { w[,cvii==i] <- rep(exp(ll_weights[,i])/sum(exp(ll_weights[,i])),times=sum(cvii==i)) }
    
    sum(log(colSums(exp(ll)*w)))
    
    model_validation_full <- data.frame(type = c("loo","10k_joint","10k_point"),
                                        validation_score = c(loo_model$estimates['elpd_loo','Estimate'],
                                                  loo_foldK_joint$estimates['elpd_loo','Estimate'],
                                                  sum(log(colSums(exp(ll)*w)))),
                                        model = model_name)
    
  }
  
  mrp_goodness <- function(model, population, model_name){
    
    popn_ps <- population %>%
      group_by(X1,X2,X3,X4)%>%
      summarise(Nj = n())%>%
      ungroup()
    
    predictions <- posterior_linpred(model, transform = TRUE, newdata = popn_ps)
    population_estimates  <- predictions %>% 
      t() %>%
      cbind(popn_ps) %>%
      pivot_longer(cols = c(-X1,-X2,-X3,-X4, -Nj),
                   names_to = "posterior_iter", 
                   values_to = "posterior_val")%>%
      group_by(posterior_iter)%>%
      summarise(mrp_popn_est = sum(posterior_val*Nj)/sum(Nj))%>%
      ungroup() %>%
      summarise(lower = as.numeric(quantile(mrp_popn_est,.05)),
                med = as.numeric(quantile(mrp_popn_est,.50)),
                upper = as.numeric(quantile(mrp_popn_est,.95)))
    
    truth = mean(population$y_prob)
    
    popnestX95 = population_estimates$upper
    popnestX5 = population_estimates$lower
    
    interval_score =   (popnestX95 - popnestX5) + 
      ((2 / .1 * (popnestX5 - truth)) * ifelse(truth < popnestX5, 1, 0)) + 
      ((2 / .1 * (truth- popnestX95)) * ifelse(truth > popnestX95, 1, 0)) 
    
    scored_truth <- data.frame(interval_score = interval_score, model = model_name)
  }
  
  
  #score without truth
  full_model_scores <- calculate_model_validation(10, full_model_fit, "full",sample, sample_size = nrow(sample))
  
  precision_model_scores <- calculate_model_validation(10, model_precision_only, "precision",sample, sample_size = nrow(sample))
  
  bias_model_scores <- calculate_model_validation(10, model_bias_only, "bias",sample, sample_size = nrow(sample))
  
  all_scores <- rbind(full_model_scores,
                      bias_model_scores,
                      precision_model_scores)
  
  #score with truth
  full_model_truth <- mrp_goodness(full_model_fit, population, "full")
  precision_model_truth <- mrp_goodness(model_precision_only, population, "precision")
  bias_model_truth <- mrp_goodness(model_bias_only, population, "bias")
  
  known_goodness <- rbind(full_model_truth, precision_model_truth, bias_model_truth)
  
  compare_scores <- all_scores %>%
    left_join(known_goodness)
  
  saveRDS(compare_scores,paste0("results/joint_vs_point10fold/iter_",ITE,".rds"))
}
complete_scores <- data.frame(readRDS(paste0("results/joint_vs_point10fold/iter_",1,".rds")), iter = 1)

for(i in 2: 10){
  complete_scores <- rbind(complete_scores, data.frame(readRDS(paste0("results/joint_vs_point10fold/iter_",i,".rds")), iter = i))
}

ggplot(complete_scores, aes(x = validation_score, y = interval_score, colour = model))+
  geom_point()+
  facet_grid(.~type)
