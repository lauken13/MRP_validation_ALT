set.seed(27472)
Pointwise_MRP_error_calc <- function(Nj, N, truth, est, error){
  if(missing(error)){
    error = truth - est
  }
  EE = matrix(Nj*error, nrow = length(Nj)) %*% matrix(Nj*error, ncol = length(Nj))
  EE[lower.tri(EE,diag = TRUE)]<-0
  ones = matrix(rep(1,length(Nj)), ncol = 1)
  sum((Nj/N*error)^2) +2/N^2*t(ones)%*%EE%*%ones
}

Pointwise_MRP_error_calc(c(1,2,3),6,c(.1,.5,.9),c(.3,.5,.7))