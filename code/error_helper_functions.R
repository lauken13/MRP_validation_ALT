population_score <-
  function(model,
           popn_counts,
           popn_obs) {
    popn_truth = popn_obs / popn_counts
    N = sum(popn_counts)
    Nj = popn_counts
    J = length(popn_counts)
    S = 4000
    
    #Predict probability for each cell
    model_preds <-  posterior_linpred(model, transform = TRUE)
    
    #calculate mrp estimate and distribution
    mrp_estimate <- apply(model_preds,1,function(x) sum(Nj*x)/N)
    
    #calculate truth
    true_value <- sum(popn_truth*Nj)/N
    
    #squared error
    true_squarederror <- median((mrp_estimate - true_value)^2) #true 
    cellwise_error <- (apply(model_preds,2,median)-popn_truth) 
    mean_cellwise_squarederror <- sum(Nj*(cellwise_error)^2)/N #summing squared cellwise using popn totals
    mrp_cellwise_squarederror <- (sum(cellwise_error*Nj)/N)^2 #estimating mrp squared error using pointwise
    
    #crps
    shuffle <- sample(1:S, size = S, replace = FALSE)
    model_preds_prime <- model_preds[shuffle,]
    mrp_estimate_prime <- apply(model_preds_prime,1,function(x) sum(Nj*x)/N)
    EXY = apply(sweep(model_preds,2,popn_truth),1,function(x) sum(Nj*x)/N)
    EXX = apply(model_preds - model_preds_prime,1,function(x) sum(Nj*x)/N)
    true_crps <- crps(mrp_estimate, mrp_estimate_prime, true_value)$estimates[1] #true crps for mrp estimate
    mean_cellwise_crps <- sum(Nj*crps(model_preds, model_preds_prime, popn_truth)$pointwise)/N #sum of pointwise crps
    mrp_cellwise_crps = .5*mean(abs(EXX))-mean(abs(EXY))
    
    results_df <- tribble(
      ~model,                   ~method,            ~score,          ~type_of_score, ~value,
      paste(formula(model))[1], "EXACT POPULATION", "CRPS",           "TRUE MRP",     true_crps,
      paste(formula(model))[1], "EXACT POPULATION", "CRPS",           "MEAN CELLWISE",mean_cellwise_crps,
      paste(formula(model))[1], "EXACT POPULATION", "CRPS",           "MRP CELLWISE",  mrp_cellwise_crps,
      paste(formula(model))[1], "EXACT POPULATION", "SQUARED ERROR",  "TRUE MRP",     true_squarederror,
      paste(formula(model))[1], "EXACT POPULATION", "SQUARED ERROR",  "MEAN CELLWISE",mean_cellwise_squarederror,
      paste(formula(model))[1], "EXACT POPULATION", "SQUARED ERROR",  "MRP CELLWISE",  mrp_cellwise_squarederror,
    )
    return(results_df)
  }

sample_score <-
  function(model,
           popn_counts,
           popn_obs,
           sample_counts,
           sample_obs) {
    popn_truth = popn_obs / popn_counts
    sample_est = sample_obs / sample_counts
    N = sum(popn_counts)
    Nj = popn_counts
    J = length(popn_counts)
    S = 4000
    
    #Predict probability for each cell
    model_preds <-  posterior_linpred(model, transform = TRUE)
    
    #calculate mrp estimate and distribution
    mrp_estimate <- apply(model_preds,1,function(x) sum(Nj*x)/N)
    
    #calculate truth
    true_value <- sum(popn_truth*Nj)/N
    
    #squared error
    true_squarederror <- median((mrp_estimate - true_value)^2) #true 
    cellwise_error <- (apply(model_preds,2,median)-sample_est) 
    mean_cellwise_squarederror <- sum(Nj*(cellwise_error)^2)/N #summing squared cellwise using popn totals
    mrp_cellwise_squarederror <- (sum(cellwise_error*Nj)/N)^2 #estimating mrp squared error using pointwise
    
    #crps
    shuffle <- sample(1:S, size = S, replace = FALSE)
    model_preds_prime <- model_preds[shuffle,]
    mrp_estimate_prime <- apply(model_preds_prime,1,function(x) sum(Nj*x)/N)
    EXY = apply(sweep(model_preds,2,sample_est),1,function(x) sum(Nj*x)/N)
    EXX = apply(model_preds - model_preds_prime,1,function(x) sum(Nj*x)/N)
    true_crps <- crps(mrp_estimate, mrp_estimate_prime, true_value)$estimates[1] #true crps for mrp estimate
    mean_cellwise_crps <- sum(Nj*crps(model_preds, model_preds_prime, sample_est)$pointwise)/N #sum of pointwise crps
    mrp_cellwise_crps = .5*mean(abs(EXX))-mean(abs(EXY))
    
    results_df <- tribble(
      ~model,                   ~method,            ~score,          ~type_of_score, ~value,
      paste(formula(model))[1], "EXACT POPULATION", "CRPS",           "TRUE MRP",     true_crps,
      paste(formula(model))[1], "SAMPLE ESTIMATE" , "CRPS",           "MEAN CELLWISE",mean_cellwise_crps,
      paste(formula(model))[1], "SAMPLE ESTIMATE" , "CRPS",           "MRP CELLWISE",  mrp_cellwise_crps,
      paste(formula(model))[1], "EXACT POPULATION", "SQUARED ERROR",  "TRUE MRP",     true_squarederror,
      paste(formula(model))[1], "SAMPLE ESTIMATE" , "SQUARED ERROR",  "MEAN CELLWISE",mean_cellwise_squarederror,
      paste(formula(model))[1], "SAMPLE ESTIMATE" , "SQUARED ERROR",  "MRP CELLWISE",  mrp_cellwise_squarederror,
    )
    return(results_df)
  }

bruteforce_loco_score <-
  function(model,
           popn_counts,
           popn_obs,
           sample_counts,
           sample_obs) {
    sample_truth = sample_obs / sample_counts
    popn_truth = popn_obs / popn_counts
    N = sum(popn_counts)
    Nj = popn_counts
    J = length(popn_counts)
    S = 4000
    
    #Predict probability for each cell
    model_preds <-  posterior_linpred(model, transform = TRUE)
    
    #calculate mrp estimate and distribution
    mrp_estimate <- apply(model_preds,1,function(x) sum(Nj*x)/N)
    
    #calculate truth
    true_value <- sum(popn_truth*Nj)/N
    
    #Leave out each cell and predict the error
    loco_prediction <- matrix(nrow = S, ncol = J)
    model_diagnostics <- data.frame(simulation_iter = ITE, cell_iter = rep(NA,J), n_divergent = rep(NA,J), rhat_larger_1_1 = rep(NA,J))
    for(i in 1:J){
      print(i/J) #progress tracker
      test_data <- model$data[i,]
      train_data <- model$data[(1:J)[1:J!=i],]
      refit_model <- brm(model$formula,
                         data = train_data,
                         family = binomial(link = "logit"), 
                         backend = "rstan", 
                         cores = 1,
                         iter = 1000,
                         save_pars = save_pars(all = TRUE))
      #save model diagnostics
      model_diagnostics$cell_iter <- i
      model_diagnostics$rhat_larger_1_1[i] <- sum(rhat(refit_model)>1.1)
      model_diagnostics$n_divergent[i] <- sum(subset(nuts_params(refit_model), Parameter == "divergent__")$Value)
      model_diagnostics$n_maxtreedepth[i] <- sum(subset(nuts_params(refit_model), Parameter == "treedepth__")$Value ==10)
      saveRDS(model_diagnostics, paste0("results/model_diagnostics/model_",paste(formula(model))[1],"_simulation_iter",ITE,"cell_iter",i,".rds"))
      #cellwise prediction
      loco_prediction[,i] <- posterior_linpred(refit_model, newdata = test_data, transform = TRUE)
      rm(refit_model)
      gc()
    }
    
    #squared error
    true_squarederror <- median((mrp_estimate - true_value)^2) #true 
    cellwise_error <- (apply(loco_prediction,2,median)-sample_truth) 
    mean_cellwise_squarederror <- sum(Nj*(cellwise_error)^2)/N #summing squared cellwise using popn totals
    mrp_cellwise_squarederror <- (sum(cellwise_error*Nj)/N)^2 #estimating mrp squared error using pointwise
    
    #crps
    shuffle <- sample(1:S, size = S, replace = FALSE)
    #needed for true crps
    model_preds_prime <- model_preds[shuffle,]
    mrp_estimate_prime <- apply(model_preds_prime,1,function(x) sum(Nj*x)/N)
    EXY = sweep(model_preds,2,popn_truth)
    EXX = model_preds - model_preds_prime
    #loco crps
    loco_prediction_prime <- loco_prediction[shuffle,]
    loco_EXY = apply(sweep(loco_prediction,2,sample_truth),1,function(x) sum(Nj*x)/N)
    loco_EXX = apply(loco_prediction - loco_prediction_prime,1,function(x) sum(Nj*x)/N)
    
    true_crps <- crps(mrp_estimate, mrp_estimate_prime, true_value)$estimates[1] #true crps for mrp estimate
    mean_cellwise_crps <- mean(.5*colMeans(abs(EXX)) + colMeans(abs(EXY)))
    mrp_cellwise_crps = .5*mean(abs(loco_EXX))-mean(abs(loco_EXY))
    
    results_df <- tribble(
      ~model,                   ~method,            ~score,          ~type_of_score, ~value,
      paste(formula(model))[1], "EXACT POPULATION", "CRPS",           "TRUE MRP",     true_crps,
      paste(formula(model))[1], "BRUTE FORCE LOCO", "CRPS",           "MEAN CELLWISE",mean_cellwise_crps,
      paste(formula(model))[1], "BRUTE FORCE LOCO", "CRPS",           "MRP CELLWISE", mrp_cellwise_crps,
      paste(formula(model))[1], "EXACT POPULATION", "SQUARED ERROR",  "TRUE MRP",     true_squarederror,
      paste(formula(model))[1], "BRUTE FORCE LOCO", "SQUARED ERROR",  "MEAN CELLWISE",mean_cellwise_squarederror,
      paste(formula(model))[1], "BRUTE FORCE LOCO", "SQUARED ERROR",  "MRP CELLWISE", mrp_cellwise_squarederror,
    )
    return(results_df)
}

approx_loco_score <-
  function(model,
           popn_counts,
           popn_obs,
           sample_counts,
           sample_obs) {
    sample_truth = sample_obs / sample_counts
    popn_truth = popn_obs / popn_counts
    N = sum(popn_counts)
    Nj = popn_counts
    J = length(popn_counts)
    S = 4000
    
    #Predict probability for each cell
    model_preds <-  posterior_linpred(model, transform = TRUE)
    
    #calculate mrp estimate and distribution
    mrp_estimate <- apply(model_preds,1,function(x) sum(Nj*x)/N)
    
    #calculate truth
    true_value <- sum(popn_truth*Nj)/N
    
    #Set up matrices for approximate LOCO using PSIS LOO
    sample_truth_matrix <- t(matrix(rep(sample_truth, S), nrow = length(sample_truth)))
    error = model_preds - sample_truth_matrix
    squared_error = (model_preds - sample_truth_matrix) ^ 2
    
    #Normal loco
    psis_loco <- loo(model, save_psis = TRUE, reloo = TRUE, cores = 1)
    psis_obj <- psis_loco$psis_object
    psis_diagnostics <- data.frame(iter = ITE, model = paste(formula(model))[1], n_paretok0_7 = sum(psis_obj$diagnostics$pareto_k>.7), percent_paretok0_7 = mean(psis_obj$diagnostics$pareto_k>.7))
    saveRDS(psis_diagnostics,paste0("results/psis_diagnostics/model_",paste(formula(model))[1],"_simulation_iter",ITE,".rds"))
    
    #Custom score
    log_lik_loco <- log_lik(model)
    
    #squared error
    psis_error <-
      E_loo(error,
            psis_obj,
            type = "mean",
            log_ratios = -log_lik_loco)$value
    
    psis_squared_error <-
      E_loo(squared_error,
            psis_obj,
            type = "mean",
            log_ratios = -log_lik_loco)$value
    
    true_squarederror <- median((mrp_estimate - true_value)^2) #true 
    mean_cellwise_squarederror <- sum(Nj*(psis_squared_error))/N #summing squared cellwise using popn totals
    mrp_cellwise_squarederror <- (sum(psis_error*Nj)/N)^2 #estimating mrp squared error using pointwise
    
    #crps
    shuffle <- sample(1:S, size = S, replace = FALSE)
    log_lik_loco_prime <- log_lik_loco[shuffle,]
    psis_obj_prime <- psis(-log_lik_loco-log_lik_loco_prime, r_eff = relative_eff(exp(-log_lik_loco-log_lik_loco_prime), chain_id = rep(1:4, each = 1000)))
   
    #needed for true crps
    model_preds_prime <- model_preds[shuffle,]
    mrp_estimate_prime <- apply(model_preds_prime,1,function(x) sum(Nj*x)/N)
    EXY = apply(sweep(model_preds,2,popn_truth),1,function(x) sum(Nj*x)/N)
    EXX = apply(model_preds - model_preds_prime,1,function(x) sum(Nj*x)/N)
    
    #loco crps
    XX = model_preds - model_preds_prime
    XY = sweep(model_preds,2,sample_truth)
    prime_weights = weights(psis_obj_prime, log = FALSE)
    weights = weights(psis_obj, log = FALSE)
    
    psis_loco_EXX = E_loo(model_preds - model_preds_prime, psis_obj_prime, log_ratios = -log_lik_loco - log_lik_loco_prime)$value
    psis_loco_EXY = E_loo(sweep(model_preds,2,sample_truth), psis_obj, log_ratios = -log_lik_loco)$value
    
    loo_crps <- loo_crps(model_preds,model_preds_prime, sample_truth, log_lik = log_lik_loco, r_eff = relative_eff(exp(-log_lik_loco), chain_id = rep(1:4, each = 1000)))
    
    Nj_mat = matrix(rep(Nj,S), nrow = S, byrow = TRUE)
    true_crps <- crps(mrp_estimate, mrp_estimate_prime, true_value)$estimates[1] #true crps for mrp estimate
    mean_cellwise_crps <- sum(Nj*loo_crps(model_preds,model_preds_prime, sample_truth, log_lik = log_lik(model), r_eff = relative_eff(exp(log_lik(model)), chain_id = rep(1:4, 1000)))$pointwise)/N #sum of pointwise loco crps
    mrp_cellwise_crps = .5*sum(abs(rowSums(prime_weights*XX*Nj_mat)/N)) - sum(abs(rowSums(weights*XY*Nj_mat)/N))
    
    results_df <- tribble(
      ~model,                   ~method,            ~score,          ~type_of_score, ~value,
      paste(formula(model))[1], "EXACT POPULATION", "CRPS",           "TRUE MRP",     true_crps,
      paste(formula(model))[1], "APPROX LOCO",      "CRPS",           "MEAN CELLWISE",mean_cellwise_crps,
      paste(formula(model))[1], "APPROX LOCO",      "CRPS",           "MRP CELLWISE", mrp_cellwise_crps,
      paste(formula(model))[1], "EXACT POPULATION", "SQUARED ERROR",  "TRUE MRP",     true_squarederror,
      paste(formula(model))[1], "APPROX LOCO",      "SQUARED ERROR",  "MEAN CELLWISE",mean_cellwise_squarederror,
      paste(formula(model))[1], "APPROX LOCO",      "SQUARED ERROR",  "MRP CELLWISE", mrp_cellwise_squarederror,
    )
    return(results_df)
}
