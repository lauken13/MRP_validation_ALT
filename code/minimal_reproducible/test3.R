library(loo)
test_model <- readRDS("test_model.rds")
loo(test_model)
