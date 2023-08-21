library(loo)
library(rstan)
set.seed(48583)

test_model <- readRDS("test_model.rds")
loo(test_model, cores=4)

t <- loo(test_model, save_psis = TRUE, cores=1,reloo = TRUE)
t
summary(t$psis_object$log_weights[,1])
summary(colMeans(t$psis_object$log_weights))

