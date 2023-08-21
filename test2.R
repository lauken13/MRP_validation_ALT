set.seed(27472)
#.libPaths(c("/gpfs/users/a1193023/local/RLibs",.libPaths()))
library(brms)
library(loo)
library(tidyverse)

rnorm(1,0,1)
test_df_sample <- expand.grid(x1 = 1:5, x2 = 1:5)%>%
  mutate(nj = sample(c(25:40), 25,replace = TRUE),
        y_prob = inv_logit_scaled(rnorm(5,0,1)[x1]+rnorm(5,0,1)[x2]),
        y_count = round(y_prob*nj))

Nj = seq(100,5000,197)
summary(test_df_sample)

model <- brm(y_count|trials(nj) ~ (1|x1) + (1|x2) ,
             data = test_df_sample, 
             family = binomial(link = "logit"), 
             backend = "cmdstanr", 
             cores = 4)
summary(model)
ranef(model)

loo_model <- loo(model, save_psis = TRUE, reloo = FALSE)
loo_model

model_preds <-  posterior_linpred(model, transform = TRUE)

colMeans(model_preds)



Pointwise_MRP_error_calc <- function(Nj, N, truth, est, error){
  if(missing(error)){
    error = truth - est
  }
  EE = matrix(Nj*error, nrow = length(Nj)) %*% matrix(Nj*error, ncol = length(Nj))
  EE[lower.tri(EE,diag = TRUE)]<-0
  ones = matrix(rep(1,length(Nj)), ncol = 1)
  sum((Nj/N*error)^2) +2/N^2*t(ones)%*%EE%*%ones
}

Pointwise_MRP_error_calc(Nj,sum(Nj),colMeans(model_preds),test_df_sample$y_prob)

(sum(colMeans(model_preds)*Nj)/sum(Nj) -sum(test_df_sample$y_prob*Nj)/sum(Nj))^2

