#' Methods to get response from fitted models
#'
#' @param fitModel Fitted model to get response
#' @param data Data to get response
#' @param ... Other arguments
#'
#' @return Vector of response
#' @export
fit.response <- function(fitModel, data, type = c('lm',
                                                  'glm', 'glm.step.AIC', 'glm.step.BIC',
                                                  'glm.best.AIC', 'glm.best.BIC',
                                                  'gam',
                                                  'glmnet1', 'glmnet2',
                                                  "adaptiveLASSO", "SCAD",
                                                  'penalized',
                                                  'stability',
                                                  'cox', 'cox.step.AIC', 'cox.step.BIC',
                                                  'cox.best.AIC', 'cox.best.BIC',
                                                  'glmnet.cox1', 'glmnet.cox2',
                                                  'rpart1', 'rpart2',
                                                  'randomForest',
                                                  'gbm', 'gbm.cox',
                                                  'coxtv', 'coxtv.step.AIC', 'coxtv.step.BIC'), ...){
  response <- switch(type,
                     lm               = fit.response.lm(fitModel, data, ...),
                     glm              = fit.response.glm(fitModel, data, ...),
                     glm.step.AIC     = fit.response.glm(fitModel, data, ...),
                     glm.step.BIC     = fit.response.glm(fitModel, data, ...),
                     glm.best.AIC     = fit.response.glm(fitModel, data, ...),
                     glm.best.BIC     = fit.response.glm(fitModel, data, ...),
                     gam              = fit.response.gam(fitModel, data, ...),
                     glmnet1          = fit.response.glmnet1(fitModel, data, ...),
                     glmnet2          = fit.response.glmnet2(fitModel, data, ...),
                     adaptiveLASSO    = fit.response.glmnet1(fitModel, data, ...),
                     SCAD             = fit.response.scad(fitModel, data, ...),
                     penalized        = fit.response.penalized(fitModel, data, ...),
                     glmnet.stability = fit.response.glm(fitModel, data, ...),
                     cox              = fit.response.cox(fitModel, data, ...),
                     cox.step.AIC     = fit.response.cox(fitModel, data, ...),
                     cox.step.BIC     = fit.response.cox(fitModel, data, ...),
                     cox.best.AIC     = fit.response.cox(fitModel, data, ...),
                     cox.best.BIC     = fit.response.cox(fitModel, data, ...),
                     glmnet.cox1      = fit.response.glmnet1(fitModel, data, ...),
                     glmnet.cox2      = fit.response.glmnet2(fitModel, data, ...),
                     rpart1           = fit.response.rpart(fitModel, data, ...),
                     rpart2           = fit.response.rpart(fitModel, data, ...),
                     randomForest     = fit.response.randomForest(fitModel, data, ...),
                     gbm              = fit.response.gbm(fitModel, data, ...),
                     gbm.cox          = fit.response.gbm(fitModel, data, ...),
                     coxtv            = fit.response.coxtv(fitModel, data, ...),
                     coxtv.step.AIC   = fit.response.coxtv(fitModel, data, ...),
                     coxtv.step.BIC   = fit.response.coxtv(fitModel, data, ...))
  response
}

#' @describeIn fit.response Linear regression
fit.response.lm <- function(fitModel, data, ...) {
  predict.lm(fitModel, newdata = data, type = "response")
}

#' @describeIn fit.response Logistic regression
fit.response.glm <- function(fitModel, data, ...) {
  predict.glm(fitModel, newdata = data, type = "response")
}

#' @describeIn fit.response Cox regression
fit.response.cox <- function(fitModel, data, ...) {
  survival::predict.coxph(fitModel, newdata = data, type = "risk")
}

#' @describeIn fit.response Cox regression with time-dependent variable
fit.response.coxtv <- function(fitModel, data, restype = c("risk", "lp"), id, save.data = TRUE, t.hor, ...) {
  if (missing(restype)) restype <- "risk"
  # change model formula for response
  model <- formula(fitModel)
  model[[2]][2] <- NULL
  # get response
  Y <- model.extract(model.frame(formula = model, data = data), "response")
  # get names of response variables (time and status)
  resvar <- as.character(model[[2]])[-1]

  if (restype == "lp"){
    newdat <- survival::survSplit(data = data, cut = sort(unique(Y[Y[, "status"] == 1, "time"])),
                        end = resvar[1], event = resvar[2], start = "Tstart")
    response <- survival::predict(fitModel, newdata = newdat, type = "lp")
  } else {
    if (restype == "risk"){
      data[, resvar[1]] <- max(Y[, "time"])
      newdat <- survival::survSplit(data = data, cut = 0:max(Y[, "time"]),
                          end = resvar[1], event = resvar[2], start = "Tstart")
      response <- getProb(fit = fitModel, newdata = newdat, t.hor = ifelse(missing(t.hor), max(Y[, "time"]), t.hor), t.lm = 0, id = id)
    } else {stop("Please check restype ! This type of response is not implemented yet.")}
  }

  # return
  y <- model.response(model.frame(model, newdat))
  if (save.data) {response <- list(response = response, data = newdat, y = as.matrix(y))}
  return(response)
}

#' @describeIn fit.response GLM LASSO (using glmnet, choose the minimum lambda)
fit.response.glmnet1 <- function(fitModel, data, ...){
  glmnet::predict.cv.glmnet(fitModel, newx = model.matrix(fitModel$model, data)[, -1], type = 'response', s = "lambda.min", ...)
}

#' @describeIn fit.response GLM LASSO (using glmnet, choose the minimum lambda within 1 se)
fit.response.glmnet2 <- function(fitModel, data, ...){
  glmnet::predict.cv.glmnet(fitModel, newx = model.matrix(fitModel$model, data)[, -1], type = 'response', s = "lambda.1se", ...)
}

#' @describeIn fit.response GLM LASSO (using penalized)
fit.response.penalized <- function(fitModel, data, ...) {
  penalized::predict(fitModel, fitModel@formula$penalized, data = data)
}

#' @describeIn fit.response GLM with SCAD (choose the minimum lambda)
fit.response.scad <- function(fitModel, data, ...) {
  ncvreg:::predict.ncvreg(fitModel, X = model.matrix(fitModel$model, data)[, -1], type = "response", lambda = fitModel$lambda.min)
}

#' @describeIn fit.response GAM
fit.response.gam <- function(fitModel, data, ...) {
  as.vector(mgcv::predict.gam(fitModel, data, type = "response", ...))
}

#' @describeIn fit.response CART
fit.response.rpart <- function(fitModel, data, ...) {
  rpart:::predict.rpart(fitModel, newdata = data, ...)
}

#' @describeIn fit.response Random forest
fit.response.randomForest <- function(fitModel, data, ...) {
  randomForest:::predict.randomForest(fitModel, newdata = data, type = "response", ...)
}

#' @describeIn fit.response Boosting
fit.response.gbm <- function(fitModel, data, ...){
  method <- ifelse(fitModel["cv.folds"] == 0, "OOB", "cv")
  best.iter <- gbm::gbm.perf(fitModel, method = method, plot.it = FALSE)
  gbm::predict.gbm(fitModel, data, n.trees = best.iter, type = "response")
}
