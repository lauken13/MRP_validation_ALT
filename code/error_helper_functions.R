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
           observed_cells_sample = NULL, #allows you to estimate errors for saes
           sample_counts,
           sample_obs,
           method = "probability" # allows you to choose between using probability or
           #nj bernoulli draws when estimating cell probability.
  ) {
    sample_truth = sample_obs / sample_counts
    popn_truth = popn_obs / popn_counts
    N = sum(popn_counts)
    Nj = popn_counts
    J = length(popn_counts)
    S = 4000
    
    #Predict probability for each cell
    if(is_empty(popn_ps)){
      if(method == "bernoulli"){
        model_preds <-  posterior_predict(model)
        for(i in 1:nrow(model_preds)){
          model_preds[i,] <- model_preds[i,]/sample_counts
        }
      }else if(method == "probability"){
        model_preds <-  posterior_linpred(model,transform = TRUE)
      } else{
        print("method input should be probability or bernoulli")
        stop()
      }
    }else{
      if(method == "bernoulli"){
        model_preds <-  posterior_predict(model, newdata = cbind(popn_ps, nj = sample_counts))
        for(i in 1:nrow(model_preds)){
          model_preds[i,] <- model_preds[i,]/sample_counts
        }
      }else if(method == "probability"){
        model_preds <-  posterior_linpred(model,newdata = popn_ps,transform = TRUE)
      } else{
        print("method input should be probability or bernoulli")
        stop()
      }
    }
    
    #if working with full ps
    if(is_empty(observed_cells_sample)){
      observed_cells_sample = rep(TRUE,length(popn_truth))
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
    error_full[,observed_cells_sample] <- error
    psis_error <-
      E_loo(error_full,
            psis_obj,
            type = "mean")$value[observed_cells_sample]
    squared_error_full <- matrix(0,nrow = S, ncol = dim(psis_obj)[2])
    squared_error_full[,observed_cells_sample] <- squared_error  
    psis_squared_error <-
      E_loo(squared_error_full,
            psis_obj,
            type = "mean")$value[observed_cells_sample]
    
    true_squarederror <- median((mrp_estimate - true_value)^2) #true 
    mean_cellwise_squarederror <- sum(Nj*(psis_squared_error))/N #summing squared cellwise using popn totals
    mrp_cellwise_squarederror <- (sum(psis_error*Nj)/N)^2 #estimating mrp squared error using pointwise
    true_mean_cellwise_squarederror = sum(Nj*(colMeans(model_preds) - popn_truth)^2)/N #true pointwise sum of cellwise squared error
    
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
    prime_weights = weights(psis_obj_prime, log = FALSE)[,observed_cells_sample]
    weights = weights(psis_obj, log = FALSE)[,observed_cells_sample]
    
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
    
    loo_crps <- loo_crps(model_preds,model_preds_prime, sample_truth, log_lik = log_lik_loco[,observed_cells_sample], r_eff = relative_eff(exp(-log_lik_loco[,observed_cells_sample]), chain_id = rep(1:4, each = 1000)))
    
    #cellwise true crps
    cellwise_true_crps <- crps(model_preds, model_preds_prime, popn_truth)$pointwise
    true_mean_cellwise_crps <- sum(Nj*cellwise_true_crps)/N
    
    Nj_mat = matrix(rep(Nj,S), nrow = S, byrow = TRUE)
    true_crps <- crps(mrp_estimate, mrp_estimate_prime, true_value)$estimates[1] #true crps for mrp estimate
    mean_cellwise_crps <- sum(Nj*loo_crps(model_preds,model_preds_prime, sample_truth, log_lik = log_lik(model)[,observed_cells_sample], r_eff = relative_eff(exp(log_lik(model)[,observed_cells_sample]), chain_id = rep(1:4, 1000)))$pointwise)/N #sum of pointwise loco crps
    mrp_cellwise_crps = (1/S)*(.5*sum(abs(rowSums(XX_resample*Nj_mat)/N)) - sum(abs(rowSums(XY_resample*Nj_mat)/N)))
    results_df <- tribble(
      ~model,                   ~method,            ~score,          ~type_of_score, ~value,
      paste(formula(model))[1], "EXACT POPULATION", "CRPS",           "TRUE MRP",     true_crps,
      paste(formula(model))[1], "APPROX LOCO",      "CRPS",           "MEAN CELLWISE",mean_cellwise_crps,
      paste(formula(model))[1], "EXACT POPULATION", "CRPS",           "MEAN CELLWISE",true_mean_cellwise_crps,
      paste(formula(model))[1], "APPROX LOCO",      "CRPS",           "MRP CELLWISE", mrp_cellwise_crps,
      paste(formula(model))[1], "EXACT POPULATION", "SQUARED ERROR",  "TRUE MRP",     true_squarederror,
      paste(formula(model))[1], "EXACT POPULATION", "SQUARED ERROR",  "MEAN CELLWISE",true_mean_cellwise_squarederror,
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
           observed_cells_sample = NULL, #allows you to estimate errors for saes
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

    #if working with full ps
    if(is_empty(observed_cells_sample)){
      observed_cells_sample = rep(TRUE,length(popn_truth))
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
    weights_reference = weights(psis_obj_reference, log = FALSE)[,observed_cells_sample]
    weights_candidate = weights(psis_obj_candidate, log = FALSE)[,observed_cells_sample]
    
    #Shuffle
    log_lik_reference = log_lik(reference_model)
    log_lik_candidate = log_lik(candidate_model)
    
    shuffle_reference <- sample(1:S, size = S, replace = FALSE)
    shuffle_candidate <- sample(1:S, size = S, replace = FALSE)
    
    log_lik_reference_prime <- log_lik_reference[shuffle_reference,observed_cells_sample]
    log_lik_candidate_prime <- log_lik_candidate[shuffle_candidate,observed_cells_sample]
    
    psis_obj_reference_prime <- psis(-log_lik_reference-log_lik_reference_prime, r_eff = relative_eff(exp(-log_lik_reference-log_lik_reference_prime), chain_id = rep(1:4, each = 1000)))
    psis_obj_candidate_prime <- psis(-log_lik_candidate-log_lik_candidate_prime, r_eff = relative_eff(exp(-log_lik_candidate-log_lik_candidate_prime), chain_id = rep(1:4, each = 1000)))
    
    weights_reference_prime = weights(psis_obj_reference_prime, log = FALSE)[,observed_cells_sample]
    weights_candidate_prime = weights(psis_obj_candidate_prime, log = FALSE)[,observed_cells_sample]
    
    model_preds_ref_prime <- model_preds_ref[shuffle_reference,observed_cells_sample]
    model_preds_candidate_prime <- model_preds_candidate[shuffle_candidate,observed_cells_sample]
    
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
    psis_error_reference <-
      E_loo(model_preds_ref,
            psis_obj_reference,
            type = "mean")$value[observed_cells_sample]

    psis_error_candidate <-
      E_loo(model_preds_candidate,
            psis_obj_candidate,
            type = "mean")$value[observed_cells_sample]
    
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

approx_combined_loco_referencemodel <-
  function(reference_model, candidate_model,
           popn_ps,
           sample_ps) {
    
    N = sum(popn_ps$Nj)
    Nj =  popn_ps$Nj
    J = length(popn_ps$Nj)
    
    ps_obs_indicator <- popn_ps %>%
      left_join(.,
                data.frame(sample_ps[c("X1","X2","X3","X4", "n_j","y_count")], observed = TRUE) %>%
                  rename(y_count_samp = "y_count")) %>%
      mutate(observed = ifelse(is.na(observed),FALSE, observed))
    
    #calculate truth
    popn_truth = popn_ps$y_count / popn_ps$Nj
    true_value <- sum(popn_truth*popn_ps$Nj)/sum(popn_ps$Nj)
    
    ###### first focus on observed portion. #############
    popn_portion_obs = ps_obs_indicator %>%
      filter(observed == TRUE)
    
    sample_truth_obs = popn_portion_obs$y_count_samp / popn_portion_obs$n_j
    popn_truth_obs = popn_portion_obs$y_count / popn_portion_obs$Nj
    N_obs = sum( popn_portion_obs$Nj)
    Nj_obs =  popn_portion_obs$Nj
    J_obs = length(popn_portion_obs$Nj)
    S = 4000
    
    #predict
    model_preds_candidate_obs <-  posterior_linpred(candidate_model,transform = TRUE)
    psis_loco_candidate_obs <- loo(model_preds_candidate_obs, save_psis = TRUE, reloo = TRUE, cores = 1)
    psis_obj_candidate_obs <- psis_loco_candidate_obs$psis_object
    weights_candidate_obs = weights(psis_obj_candidate_obs, log = FALSE)
    log_lik_candidate_obs = log_lik(candidate_model)

    psis_diagnostics <- data.frame(iter = ITE, model = paste(formula(candidate_model))[1], n_paretok0_7 = sum(psis_obj_candidate_obs$diagnostics$pareto_k>.7), percent_paretok0_7 = mean(psis_obj_candidate_obs$diagnostics$pareto_k>.7))
    saveRDS(psis_diagnostics,paste0("results/psis_diagnostics/combined_reference/model_",paste(formula(candidate_model))[1],"_simulation_iter",ITE,".rds"))
    
    #Shuffle
    shuffle_candidate_obs <- sample(1:S, size = S, replace = FALSE)
    log_lik_candidate_prime_obs <- log_lik_candidate_obs[shuffle_candidate_obs,]
    psis_obj_candidate_prime_obs <- psis(-log_lik_candidate_obs-log_lik_candidate_prime_obs, r_eff = relative_eff(exp(-log_lik_candidate_obs-log_lik_candidate_prime_obs), chain_id = rep(1:4, each = 1000)))
    weights_candidate_prime_obs = weights(psis_obj_candidate_prime_obs, log = FALSE)[,]
    model_preds_candidate_prime_obs <- model_preds_candidate_obs[shuffle_candidate_obs,]
    
    # squared error
    sample_truth_matrix_obs <- t(matrix(rep(sample_truth_obs, S), nrow = length(sample_truth_obs)))
    candidate_error_obs = model_preds_candidate_obs - sample_truth_matrix_obs
    psis_candidate_error_obs <- E_loo(candidate_error_obs,
            psis_obj_candidate_obs,
            type = "mean")$value
    
    #loco crps
    XX_obs = model_preds_candidate_obs - model_preds_candidate_prime_obs
    XY_obs = sweep(model_preds_candidate_obs,2,sample_truth_obs)
    
    candidate_model_preds_resample = matrix(nrow = S, ncol = J_obs)
    candidate_model_preds_resample_prime = matrix(nrow = S, ncol = J_obs)
    for(j in 1:J_obs){
      candidate_model_preds_resample[,j] <- resample_draws(as_draws_matrix(model_preds_candidate_obs)[,j],
                                                 weights = weights_candidate_obs[,j])
      candidate_model_preds_resample_prime[,j] <- resample_draws(as_draws_matrix(model_preds_candidate_prime_obs)[,j],
                                                       weights = weights_candidate_prime_obs[,j]) 
    }
    XX_resample_obs = candidate_model_preds_resample - candidate_model_preds_resample_prime
    XY_resample_obs = sweep(candidate_model_preds_resample,2,sample_truth_obs)
    
    #Now focus on unobserved portion
    popn_portion_unobs = ps_obs_indicator %>%
      filter(observed == FALSE)
    
    sample_truth_unobs = popn_portion_unobs$y_count_samp / popn_portion_unobs$n_j
    popn_truth_unobs = popn_portion_unobs$y_count / popn_portion_unobs$Nj
    N_unobs = sum( popn_portion_unobs$Nj)
    Nj_unobs =  popn_portion_unobs$Nj
    J_unobs = length(popn_portion_unobs$Nj)
    
    #Predict probability for each unobserved cell
    model_preds_ref_unobs <-  posterior_linpred(reference_model, 
                                                newdata = popn_portion_unobs %>%
                                                  mutate(n_j = Nj), 
                                                transform = TRUE,
                                                allow_new_levels = TRUE,
                                                sample_new_levels = "gaussian")
    model_preds_candidate_unobs <-  posterior_linpred(candidate_model, 
                                                      newdata = popn_portion_unobs %>%
                                                        mutate(n_j = Nj), 
                                                      transform = TRUE,
                                                      allow_new_levels = TRUE,
                                                      sample_new_levels = "gaussian") 
    
    #squared error
    reference_model_error_unobs <- colMeans(model_preds_candidate_unobs)-colMeans(model_preds_ref_unobs)
    
    #crps
    shuffle_reference_unobs <- sample(1:S, size = S, replace = FALSE)
    shuffle_candidate_unobs <- sample(1:S, size = S, replace = FALSE)
    
    model_preds_reference_unobs_prime <- model_preds_ref_unobs[shuffle_reference_unobs,]
    model_preds_candidate_unobs_prime <- model_preds_candidate_unobs[shuffle_candidate_unobs,]
    
    XX_candidate_unobs <- model_preds_candidate_unobs - model_preds_candidate_unobs_prime
    YY_reference_unobs <- model_preds_ref_unobs - model_preds_reference_unobs_prime
    XY_unobs <- model_preds_candidate_unobs - model_preds_ref_unobs
  
    mrp_combined_reference_squared_error <- (sum(c(Nj_obs*psis_candidate_error_obs,Nj_unobs*reference_model_error_unobs))/N)^2
    
    Nj_mat = matrix(rep(Nj,S), nrow = S, byrow = TRUE)  
    Nj_mat_obs = matrix(rep(Nj_obs,S), nrow = S, byrow = TRUE)  
    Nj_mat_unobs = matrix(rep(Nj_unobs,S), nrow = S, byrow = TRUE)  
    mrp_combined_reference_crps = (1/S)*(.5*sum(abs((rowSums(XX_resample_obs*Nj_mat_obs)+ #loco score for obs
                                                      rowSums(XX_candidate_unobs*Nj_mat_unobs))/N)) + #ref for unobs
                                  .5*sum(abs((rowSums(YY_reference_unobs*Nj_mat_unobs))/N_unobs)) - #ref for unobs as Y-Y' = 0 when truth known
                                  sum(abs((rowSums(XY_resample_obs*Nj_mat_obs)+ #loco score for obs
                                             rowSums(XY_unobs*Nj_mat_unobs))/N))) #ref for unobs
    
    mrp_cellwise_crps_resample = 
      results_df <- tribble(
        ~model,                   ~method,            ~score,          ~type_of_score, ~value,
        paste(formula(candidate_model))[1], "APPROX COMBINED LOCO-REFERENCE",      "CRPS",           "MRP CELLWISE", mrp_combined_reference_crps,
        paste(formula(candidate_model))[1], "APPROX COMBINED LOCO-REFERENCE",      "SQUARED ERROR",  "MRP CELLWISE", mrp_combined_reference_squared_error
      )
    return(results_df)
  }


partially_obs_approx_loco_score <-
  function(reference_model,
           candidate_model,
           popn_ps,
           sample_ps) {
    
    stopifnot("The population is not partially observed" = nrow(sample_ps) != nrow(popn_ps))
    
    popn_ps_obs <- popn_ps %>%
      left_join(.,
                data.frame(sample_ps[c("X1","X2","X3","X4", "n_j","y_count")], observed = TRUE) %>%
                  rename(y_count_samp = "y_count")) %>%
      mutate(observed = ifelse(is.na(observed),FALSE, observed))
    
    #calculate scores for observed portion
    
    observed_cells_scores <- approx_loco_score(candidate_model,
                                               popn_counts = popn_ps_obs$Nj[popn_ps_obs$observed==TRUE],
                                               popn_obs = popn_ps_obs$y_count[popn_ps_obs$observed==TRUE],
                                               sample_counts = popn_ps_obs$n_j[popn_ps_obs$observed==TRUE],
                                               sample_obs = popn_ps_obs$y_count_samp[popn_ps_obs$observed==TRUE])%>%
      mutate(type_of_score = paste0("PARTIAL ", type_of_score))
    
    observed_cells_referencemodel<- approx_loco_referencemodel(reference_model = reference_model,
                                                               candidate_model = candidate_model,
                                                               popn_counts = popn_ps_obs$Nj[popn_ps_obs$observed==TRUE],
                                                               popn_obs = popn_ps_obs$y_count[popn_ps_obs$observed==TRUE],
                                                               sample_counts = popn_ps_obs$n_j[popn_ps_obs$observed==TRUE],
                                                               sample_obs = popn_ps_obs$y_count_samp[popn_ps_obs$observed==TRUE]) %>%
      mutate(type_of_score = paste0("PARTIAL ", type_of_score))
    
    
    combined_reference_model_approach <-approx_combined_loco_referencemodel(reference_model = reference_model, 
                                                                            candidate_model = candidate_model,
                                                                            popn_ps = popn_ps,
                                                                            sample_ps = sample_ps) 
    results_df = rbind(observed_cells_scores, observed_cells_referencemodel, combined_reference_model_approach)
    return(results_df)
  }


