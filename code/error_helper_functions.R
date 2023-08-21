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
