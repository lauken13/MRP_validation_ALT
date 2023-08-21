library(loo)

set.seed(48583)

test_model <- readRDS("test_model.rds")
loo(test_model)
