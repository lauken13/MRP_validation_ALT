validate_score <-
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
    mean_cellwise_crps <- crps(model_preds, model_preds_prime, popn_truth)$estimates[1] #sum of pointwise crps
    mrp_cellwise_crps = .5*mean(abs(EXX))-mean(abs(EXY))
    
    results_df <- tribble(
      ~model,                   ~score,          ~type_of_score, ~value,
      paste(formula(model))[1], "CRPS",           "TRUE MRP",     true_crps,
      paste(formula(model))[1], "CRPS",           "MEAN CELLWISE",mean_cellwise_crps,
      paste(formula(model))[1], "CRPS",           "MRP CELLWISE",  mrp_cellwise_crps,
      paste(formula(model))[1], "SQUARED ERROR",  "TRUE MRP",     true_squarederror,
      paste(formula(model))[1], "SQUARED ERROR",  "MEAN CELLWISE",mean_cellwise_squarederror,
      paste(formula(model))[1], "SQUARED ERROR",  "MRP CELLWISE",  mrp_cellwise_squarederror,
    )
    return(results_df)
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




MSE_calc <- function(Nj,N,truth, est){
  sum(Nj*(est-truth)^2)/N 
}

True_MRP_error_calc <- function(Nj,N,truth, est){
  (sum(est*Nj)/N-sum(truth*Nj)/N )^2
}

Pointwise_MRP_error_calc <- function(Nj, N, truth, est, error){
  if(missing(error)){
    error = est - truth
  }
  EE = matrix(Nj*error, nrow = length(Nj)) %*% matrix(Nj*error, ncol = length(Nj))
  EE[lower.tri(EE,diag = TRUE)]<-0
  ones = matrix(rep(1,length(Nj)), ncol = 1)
  sum((Nj/N*error)^2) +2/N^2*t(ones)%*%EE%*%ones
}


