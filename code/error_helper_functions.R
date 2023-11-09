population_score <-
  function(model,
           popn_counts,
           popn_obs,
           popn_ps = NULL) {
    popn_truth = popn_obs / popn_counts
    N = sum(popn_counts)
    Nj = popn_counts
    J = length(popn_counts)
    S = 4000
    
    #Predict probability for each cell
    if(is.null(popn_ps)){
      model_preds <-  posterior_linpred(model,transform = TRUE)
    }else{
      model_preds <-  posterior_linpred(model,newdata = popn_ps, transform = TRUE)
    }

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
           popn_ps = NULL,
           sample_counts,
           sample_obs) {
    popn_truth = popn_obs / popn_counts
    sample_est = sample_obs / sample_counts
    N = sum(popn_counts)
    Nj = popn_counts
    J = length(popn_counts)
    S = 4000
    
    #Predict probability for each cell
    if(is.null(popn_ps)){
      model_preds <-  posterior_linpred(model,transform = TRUE)
    }else{
      model_preds <-  posterior_linpred(model,newdata = popn_ps, transform = TRUE)
    }
    
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
           popn_ps = NULL,
           sample_counts,
           sample_obs) {
    sample_truth = sample_obs / sample_counts
    popn_truth = popn_obs / popn_counts
    N = sum(popn_counts)
    Nj = popn_counts
    J = length(popn_counts)
    S = 4000
    
    #Predict probability for each cell
    if(is.null(popn_ps)){
      model_preds <-  posterior_linpred(model,transform = TRUE)
    }else{
      model_preds <-  posterior_linpred(model,newdata = popn_ps, transform = TRUE)
    }
    
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
      #cellwise prediction
      loco_prediction[,i] <- posterior_linpred(refit_model, newdata = test_data, transform = TRUE)
      rm(refit_model)
      gc()
    }
    saveRDS(model_diagnostics, paste0("results/model_diagnostics/model_",paste(formula(model))[1],"_simulation_iter",ITE,".rds"))
    
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
           popn_ps = NULL,
           observed_cells = NULL,
           sample_counts,
           sample_obs) {
    sample_truth = sample_obs / sample_counts
    popn_truth = popn_obs / popn_counts
    N = sum(popn_counts)
    Nj = popn_counts
    J = length(popn_counts)
    S = 4000
    
    #Predict probability for each cell
    if(is_empty(popn_ps)){
      model_preds <-  posterior_linpred(model,transform = TRUE)
    }else{
      model_preds <-  posterior_linpred(model,newdata = popn_ps, transform = TRUE)
    }
    
    #if working with full ps
    if(is_empty(observed_cells)){
      observed_cells = rep(TRUE,length(popn_truth))
    }
    
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
    error_full <- matrix(0,nrow = S, ncol = dim(psis_obj)[2])
    error_full[,observed_cells] <- error
    psis_error <-
      E_loo(error_full,
            psis_obj,
            type = "mean")$value[observed_cells]
    squared_error_full <- matrix(0,nrow = S, ncol = dim(psis_obj)[2])
    squared_error_full[,observed_cells] <- squared_error  
    psis_squared_error <-
      E_loo(squared_error_full,
            psis_obj,
            type = "mean")$value[observed_cells]
    
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
    prime_weights = weights(psis_obj_prime, log = FALSE)[,observed_cells]
    weights = weights(psis_obj, log = FALSE)[,observed_cells]
    
    model_preds_resample = matrix(nrow = S, ncol = J)
    model_preds_resample_prime = matrix(nrow = S, ncol = J)
    for(j in 1:J){
      model_preds_resample[,j] <- resample_draws(as_draws_matrix(model_preds)[,j],
                                             weights = weights[,j])
      model_preds_resample_prime[,j] <- resample_draws(as_draws_matrix(model_preds_prime)[,j],
                                                       weights = prime_weights[,j]) 
    }
    XX_resample = model_preds_resample - model_preds_resample_prime
    XY_resample = sweep(model_preds_resample,2,sample_truth)
    
    loo_crps <- loo_crps(model_preds,model_preds_prime, sample_truth, log_lik = log_lik_loco[,observed_cells], r_eff = relative_eff(exp(-log_lik_loco[,observed_cells]), chain_id = rep(1:4, each = 1000)))
    
    Nj_mat = matrix(rep(Nj,S), nrow = S, byrow = TRUE)
    true_crps <- crps(mrp_estimate, mrp_estimate_prime, true_value)$estimates[1] #true crps for mrp estimate
    mean_cellwise_crps <- sum(Nj*loo_crps(model_preds,model_preds_prime, sample_truth, log_lik = log_lik(model)[,observed_cells], r_eff = relative_eff(exp(log_lik(model)[,observed_cells]), chain_id = rep(1:4, 1000)))$pointwise)/N #sum of pointwise loco crps
    mrp_cellwise_crps = (1/S)*(.5*sum(abs(rowSums(XX_resample*Nj_mat)/N)) - sum(abs(rowSums(XY_resample*Nj_mat)/N)))
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

approx_loco_referencemodel <-
  function(reference_model, candidate_model,
           popn_counts,
           popn_obs,
           popn_ps = NULL,
           sample_counts,
           sample_obs) {
    sample_truth = sample_obs / sample_counts
    popn_truth = popn_obs / popn_counts
    N = sum(popn_counts)
    Nj = popn_counts
    J = length(popn_counts)
    S = 4000
    
    #Predict probability for each cell
    if(is_null(popn_ps)){
      model_preds_ref <-  posterior_linpred(reference_model, transform = TRUE)
      model_preds_candidate <-  posterior_linpred(candidate_model, transform = TRUE)      
    }else{
      model_preds_ref <-  posterior_linpred(reference_model, newdata = popn_ps, transform = TRUE)
      model_preds_candidate <-  posterior_linpred(candidate_model, newdata = popn_ps, transform = TRUE)   
    }

    
    #calculate mrp estimate and distribution
    mrp_estimate_ref <- apply(model_preds_ref,1,function(x) sum(Nj*x)/N)
    mrp_estimate_candidate <- apply(model_preds_candidate,1,function(x) sum(Nj*x)/N)
    
    #calculate truth
    true_value <- sum(popn_truth*Nj)/N
    
    #PSIS
    psis_loco_reference<- loo(reference_model, save_psis = TRUE, reloo = TRUE, cores = 1)
    psis_loco_candidate <- loo(candidate_model, save_psis = TRUE, reloo = TRUE, cores = 1)
    psis_obj_reference <- psis_loco_reference$psis_object
    psis_obj_candidate <- psis_loco_candidate$psis_object
    psis_diagnostics_candidate <- data.frame(iter = ITE, model = paste(formula(candidate_model))[1], n_paretok0_7 = sum(psis_obj_candidate$diagnostics$pareto_k>.7), percent_paretok0_7 = mean(psis_obj_candidate$diagnostics$pareto_k>.7))
    psis_diagnostics_reference <- data.frame(iter = ITE, model = paste(formula(reference_model))[1], n_paretok0_7 = sum(psis_obj_reference$diagnostics$pareto_k>.7), percent_paretok0_7 = mean(psis_obj_reference$diagnostics$pareto_k>.7))
    psis_diagnostics <- rbind(psis_diagnostics_candidate,psis_diagnostics_reference)
    saveRDS(psis_diagnostics,paste0("results/psis_diagnostics/model_",paste(formula(candidate_model))[1],"_simulation_iter",ITE,".rds"))
    
    #Obtain weights
    weights_reference = weights(psis_obj_reference, log = FALSE)
    weights_candidate = weights(psis_obj_candidate, log = FALSE)
    
    #Shuffle
    log_lik_reference = log_lik(reference_model)
    log_lik_candidate = log_lik(candidate_model)
    
    shuffle_reference <- sample(1:S, size = S, replace = FALSE)
    shuffle_candidate <- sample(1:S, size = S, replace = FALSE)
    
    log_lik_reference_prime <- log_lik_reference[shuffle_reference,]
    log_lik_candidate_prime <- log_lik_candidate[shuffle_candidate,]
    
    psis_obj_reference_prime <- psis(-log_lik_reference-log_lik_reference_prime, r_eff = relative_eff(exp(-log_lik_reference-log_lik_reference_prime), chain_id = rep(1:4, each = 1000)))
    psis_obj_candidate_prime <- psis(-log_lik_candidate-log_lik_candidate_prime, r_eff = relative_eff(exp(-log_lik_candidate-log_lik_candidate_prime), chain_id = rep(1:4, each = 1000)))
    
    weights_reference_prime = weights(psis_obj_reference_prime, log = FALSE)
    weights_candidate_prime = weights(psis_obj_candidate_prime, log = FALSE)
    
    model_preds_ref_prime <- model_preds_ref[shuffle_reference,]
    model_preds_candidate_prime <- model_preds_candidate[shuffle_candidate,]
    
    #resample
    reference_preds_resample = matrix(nrow = S, ncol = J)
    reference_preds_resample_prime = matrix(nrow = S, ncol = J)
    candidate_preds_resample = matrix(nrow = S, ncol = J)
    candidate_preds_resample_prime = matrix(nrow = S, ncol = J)
    for(j in 1:J){
      reference_preds_resample[,j] <- resample_draws(as_draws_matrix(model_preds_ref)[,j],
                                                 weights = weights_reference[,j])
      reference_preds_resample_prime[,j] <- resample_draws(as_draws_matrix(model_preds_ref_prime)[,j],
                                                       weights = weights_reference_prime[,j]) 
      candidate_preds_resample[,j] <- resample_draws(as_draws_matrix(model_preds_candidate)[,j],
                                                     weights = weights_candidate[,j])
      candidate_preds_resample_prime[,j] <- resample_draws(as_draws_matrix(model_preds_candidate_prime)[,j],
                                                           weights = weights_candidate_prime[,j]) 
    }
    #Squared error
    sample_truth_matrix <- t(matrix(rep(sample_truth, S), nrow = length(sample_truth)))
    error_reference = model_preds_ref - sample_truth_matrix
    error_candidate = model_preds_candidate - sample_truth_matrix
    psis_error_reference <-
      E_loo(error_reference,
            psis_obj_reference,
            type = "mean",
            log_ratios = -log_lik_reference)$value
    psis_error_candidate <-
      E_loo(error_candidate,
            psis_obj_candidate,
            type = "mean",
            log_ratios = -log_lik_candidate)$value
  
    
    mrp_cellwise_squared_candidate_ref <- (sum((psis_error_candidate - psis_error_reference)*Nj)/N)^2 #reference vs candidate
    mrp_true_squarederror <- median((mrp_estimate_candidate - mrp_estimate_ref)^2) #true 
    
    #CRPS
    
    #loco crps
    XX_resample_candidate = candidate_preds_resample - candidate_preds_resample_prime
    YY_resample_reference = reference_preds_resample - reference_preds_resample_prime
    XY_resample = reference_preds_resample - candidate_preds_resample
    
    Nj_mat = matrix(rep(Nj,S), nrow = S, byrow = TRUE)
    mrp_reference_crps = (1/S)*(.5*sum(abs(rowSums(XX_resample_candidate*Nj_mat)/N)) +
                                  .5*sum(abs(rowSums(YY_resample_reference*Nj_mat)/N)) -
                                 sum(abs(rowSums(XY_resample*Nj_mat)/N)))
    
    mrp_cellwise_crps_resample = 
      results_df <- tribble(
        ~model,                   ~method,            ~score,          ~type_of_score, ~value,
        paste(formula(candidate_model))[1], "APPROX LOCO REFERENCE",      "CRPS",           "MRP CELLWISE", mrp_reference_crps,
        paste(formula(candidate_model))[1], "APPROX LOCO REFERENCE",      "SQUARED ERROR",  "MRP CELLWISE", mrp_cellwise_squared_candidate_ref
      )
    return(results_df)
  }

approx_loco_sae_score <-
  function(model,
           small_area_var,
           popn_ps,
           sample_ps) {
    
    levels_of_small_area_var = levels(as.factor(popn_ps[[small_area_var]]))
    
    for(k in 1:length(levels_of_small_area_var)){
      target_area <- levels_of_small_area_var[[k]]
      popn_ps_sae <- popn_ps %>%
        filter(.data[[small_area_var]] == target_area) %>%
        mutate(n_j = Nj)
      sample_ps_sae <- sample_ps %>%
        filter(.data[[small_area_var]] == target_area)
      popn_counts_sae <-popn_ps_sae$Nj
      popn_obs_sae <- popn_ps_sae$y_count
      sample_counts_sae <- sample_ps_sae$n_j
      sample_obs_sae <- sample_ps_sae$y_count
      loco_approx_score_sae = approx_loco_score(model,
                 popn_counts = popn_counts_sae,
                 popn_obs = popn_obs_sae,
                 popn_ps = popn_ps_sae,
                 observed_cells = popn_ps[[small_area_var]] == target_area, 
                 sample_counts = sample_counts_sae,
                 sample_obs = sample_obs_sae) 
      if(target_area == levels_of_small_area_var[1]){
        results_df = data.frame(loco_approx_score_sae, level = target_area, variable = small_area_var)
      }else{
        results_df = rbind(results_df,
                           data.frame(loco_approx_score_sae, level = target_area, variable = small_area_var))
      }
    }
    return(results_df)
  }

