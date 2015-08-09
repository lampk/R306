## packages need to be import
pkg <- c("survival", "glmnet", "penalized",
         "ncvreg", "mgcv", "rpart", "randomForest", "gbm", "glmulti", "caret", "ROCR", "c060")
sapply(pkg, function(x) devtools::use_package(x))

## build documentation
devtools::document()
