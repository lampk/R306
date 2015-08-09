#' Functions to assess model performance
#'
#' @param predicted Predicted values from fitted models
#' @param y Response
#'
#' @return numeric value of model's performance
#' @export
myaccuracy <- function(predicted, y, criteria = mse, ...){

  output <- vector("list", length(criteria))
  for (i in (1:length(criteria))){
    output[[i]] <- criteria[[i]](predicted, y, ...)
  }
  return(unlist(output))
}

#' @describeIn myaccuracy Mean squared errors
#' @export
mse <- function(predicted, y){
  output <- mean((y - predicted)^2)
  names(output) <- "MSE"
  return(output)
}

#' @describeIn myaccuracy log-likelihood (binary case)
#' @export
loglik_binary <- function(predicted, y){
  #- replace 0 and 1 predictions (occur e.g. in rpart) by "almost" 0 and 1
  predicted[predicted < .Machine$double.eps]     <- .Machine$double.eps
  predicted[predicted > 1 - .Machine$double.eps] <- 1 - .Machine$double.eps
  #- get output
  output <- sum(y * log(predicted) + (1 - y) * log(1 - predicted))
  names(output) <- 'loglik'
  #- out
  return(output)
}

#' @describeIn myaccuracy Brier score
#' @export
Brier <- function(predicted, y, confint = FALSE){
  #- replace 0 and 1 predictions (occur e.g. in rpart) by "almost" 0 and 1
  predicted[predicted < .Machine$double.eps]     <- .Machine$double.eps
  predicted[predicted > 1 - .Machine$double.eps] <- 1 - .Machine$double.eps
  #- get output
  output <- mean((predicted - y)^2)
  names(output) <- 'Brier'
  #- out
  return(output)
}

#' @describeIn myaccuracy Discrimination slope
#' @export
dis_slope <- function(predicted, y){
  #- replace 0 and 1 predictions (occur e.g. in rpart) by "almost" 0 and 1
  predicted[predicted < .Machine$double.eps]     <- .Machine$double.eps
  predicted[predicted > 1 - .Machine$double.eps] <- 1 - .Machine$double.eps
  #- get output
  output <- mean(predicted[which(y == 1)]) - mean(predicted[which(y == 0)])
  names(output) <- 'disSlope'
  #- out
  return(output)
}

#' @describeIn myaccuracy AUC (with confidence interval)
#' @export
c_index <- function(predicted, y, confint = FALSE, ci_collapse = "-", ci_parentheses = FALSE){
  #- replace 0 and 1 predictions (occur e.g. in rpart) by "almost" 0 and 1
  predicted[predicted < .Machine$double.eps] <- .Machine$double.eps
  predicted[predicted > 1 - .Machine$double.eps] <- 1 - .Machine$double.eps
  #- get output
  w <- Hmisc::rcorr.cens(predicted, y)
  C <- w['C Index']
  if (confint){
    se <- w['S.D.']/2; CI <- paste(formatC(C + 1.96 * c(-se, se), format = 'f', digits = 2), collapse = ci_collapse)
    output <- c(formatC(C, format = 'f', digits = 2), ifelse(ci_parentheses, paste("(", CI, ")", sep = ""), CI)); names(output) <- c('AUC', 'AUC.95CI')
  } else {
    output <- C; names(output) <- 'AUC'
  }
  #- out
  return(output)
}

#' @describeIn myaccuracy Concordance
#' @export
concordance <- function(x, y){
  # y: Surv object in newdata
  # x: linear predictor in newdata
  temp <- survival::survConcordance.fit(y, x)
  out <- c((temp[1] + temp[3]/2)/sum(temp[1:3]), temp[5]/(2 * sum(temp[1:3])))
  names(out) <- c("concordance", "se")
  out
}

#' @describeIn myaccuracy Calibration slope and intercept (with confidence intervals)
#' @export
calib <- function(predicted, y, confint = FALSE){
  #- replace 0 and 1 predictions (occur e.g. in rpart) by "almost" 0 and 1
  predicted[predicted < .Machine$double.eps] <- .Machine$double.eps
  predicted[predicted > 1 - .Machine$double.eps] <- 1 - .Machine$double.eps
  #- get output
  logitPredicted <- qlogis(predicted)
  cal_fit1 <- glm(y ~ offset(logitPredicted), family = 'binomial')
  cal_fit2 <- glm(y ~ logitPredicted, family = 'binomial')
  cal_intheLarge <- coef(cal_fit1)
  cal_other <- coef(glm(y ~ logitPredicted, family = 'binomial'))
  name_out <- paste('Calibration', c('in-the-large', 'intercept', 'slope'))
  if (confint){
    confint_large <- formatC(confint(cal_fit1), format = 'f', digits = 2); confint_other <- formatC(confint(cal_fit2), format = 'f', digits = 2)
    intheLarge <- c(formatC(cal_intheLarge, format = 'f', digits = 2), paste(confint_large, collapse = '-'))
    other <- c(rbind(formatC(cal_other, format = 'f', digits = 2), paste(confint_other[,1], confint_other[,2], sep = '-')))
    output <- c(intheLarge, other)
    name_out2 <- expand.grid(c('', '.95CI'), name_out); names(output) <- paste(name_out2[, 2], name_out2[, 1], sep = '')
  } else {
    output <- c(cal_intheLarge, cal_other)
    names(output) <- name_out
  }
  #- out
  return(output)
}

#' @describeIn myaccuracy Other measures of binary test
#' @export
binperform <- function(predicted, y, cutoff = 0.5){
  #- replace 0 and 1 predictions (occur e.g. in rpart) by "almost" 0 and 1
  predicted[predicted < .Machine$double.eps] <- .Machine$double.eps
  predicted[predicted > 1 - .Machine$double.eps] <- 1 - .Machine$double.eps
  #- get output
  tmp <- caret::confusionMatrix(factor(predicted>cutoff,levels=c(T,F)),factor(as.logical(y),levels=c(T,F)))
  output <- c(tmp$table, tmp$overall[1],tmp$byClass[1:4])
  names(output)[1:4] <- c("TP", "FN", "FP", "TN")
  names(output) <- paste(names(output)," (cutoff ",cutoff,")",sep="")
  #- out
  return(output)
}

#' @describeIn myaccuracy Sensitivity
#' @export
sen <- function(predicted, y, cutoff) {
  sum(predicted[y == 1] > cutoff, na.rm = TRUE)/sum(y == 1 & !is.na(predicted), na.rm = TRUE)
}

#' @describeIn myaccuracy Specificity
#' @export
spe <- function(predicted, y, cutoff) {
  sum(predicted[y == 0] <= cutoff, na.rm = TRUE)/sum(y == 0 & !is.na(predicted), na.rm = TRUE)
}

#' @describeIn myaccuracy Positive predictive value
#' @export
ppv <- function(predicted, y, cutoff) {
  sum(y[predicted > cutoff] == 1, na.rm = TRUE)/sum(predicted > cutoff, na.rm = TRUE)
}

#' @describeIn myaccuracy Negative predictive value
#' @export
npv <- function(predicted, y, cutoff) {
  sum(y[predicted <= cutoff] == 0, na.rm = TRUE)/sum(predicted <= cutoff, na.rm = TRUE)
}

#' @describeIn myaccuracy Accuracy rate
#' @export
acc <- function(predicted, y, cutoff) {
  (sum(y[predicted > cutoff] == 1, na.rm = TRUE) + sum(y[predicted <= cutoff] == 0, na.rm = TRUE))/length(y)
}

#' @describeIn myaccuracy AUC
#' @export
auc <- function(predicted, y, alpha = 0.05, digits = 2) {
  w <- Hmisc::rcorr.cens(predicted, y)
  C <- w['C Index']; se <- w['S.D.']/2
  output <- c(C, C + qnorm(1 - 0.5 * alpha) * c(-se, se))
  names(output) <- c("AUC", "AUC.lower", "AUC.upper")
  return(round(output, digits = digits))
}

#' @describeIn myaccuracy Best Youden's index
#' @export
best.youden <- function(predicted, y, by = 0.01, digits = 2, name = "") {
  # to find "best" cut-off
  cutoff <- seq(from = 0, to = 1, by = by)
  y <- y[!is.na(predicted)]; predicted <- predicted[!is.na(predicted)]
  tmp <- sapply(cutoff, function(x, .y = y, .p = predicted){
    sen(predicted = .p, y = .y, cutoff = x) + spe(predicted = .p, y = .y, cutoff = x)
  })
  best <- cutoff[tmp == max(tmp)]
  data.frame(model = name,
             cbind(cutoff = best,
                   sensitivity = as.numeric(sapply(best, sen, predicted = predicted, y = y)),
                   specificity = as.numeric(sapply(best, spe, predicted = predicted, y = y)),
                   PPV = as.numeric(sapply(best, ppv, predicted = predicted, y = y)),
                   NPV = as.numeric(sapply(best, npv, predicted = predicted, y = y))))
}

#' @describeIn myaccuracy Best accuracy rate
#' @export
best.acc <- function(predicted, y, by = 0.01, digits = 2, name = "") {
  # to find cut-off with best accuracy [(TP + TN) / (P + N)]
  cutoff <- seq(from = 0, to = 1, by = by)
  y <- y[!is.na(predicted)]; predicted <- predicted[!is.na(predicted)]
  tmp <- sapply(cutoff, function(x){acc(predicted = predicted, y = y, cutoff = x)})
  best <- max(cutoff[tmp == max(tmp)])
  data.frame(model = name,
             cbind(cutoff = best,
                   accuracy = as.numeric(sapply(best, acc, predicted = predicted, y = y)),
                   sensitivity = as.numeric(sapply(best, sen, predicted = predicted, y = y)),
                   specificity = as.numeric(sapply(best, spe, predicted = predicted, y = y)),
                   PPV = as.numeric(sapply(best, ppv, predicted = predicted, y = y)),
                   NPV = as.numeric(sapply(best, npv, predicted = predicted, y = y))))
}

#' @describeIn myaccuracy Accuracy of binary test
#' @export
new.accuracy <- function(predicted, y, cutoff=0.5){
  require(caret); require(Hmisc)
  #- replace 0 and 1 predictions (occur e.g. in rpart) by "almost" 0 and 1
  predicted[predicted<.Machine$double.eps] <- .Machine$double.eps
  predicted[predicted>1-.Machine$double.eps] <- 1-.Machine$double.eps
  # Likelihood (Lam added)
  loglik <- sum(y*log(predicted) + (1-y)*log(1-predicted))
  names(loglik) <- 'log-likelihood'
  # Brier's score (Lam added)
  brier <- mean((predicted - y)^2)
  names(brier) <- 'Brier'
  # Discrimination slope
  disSlope <- mean(predicted[which(y==1)]) - mean(predicted[which(y==0)])
  names(disSlope) <- 'Discrimination slope'
  #- AUC and calibrarion intercept and slope
  auc <- somers2(predicted,y)["C"]
  names(auc) <- "AUC"
  #browser()
  logitPredicted <- binomial()$linkfun(predicted)
  intSlope <- coef(glm(y~logitPredicted,family=binomial()))
  names(intSlope) <- paste("calibration",c("intercept","slope"))
  #- Measures for the binary test
  tmp <- confusionMatrix(factor(predicted>cutoff,levels=c(T,F)),factor(as.logical(y),levels=c(T,F)))
  binperformance <- c(tmp$overall[1],tmp$byClass[1:4])
  names(binperformance) <- paste(names(binperformance)," (cutoff ",cutoff,")",sep="")
  #- all in one vector
  c(loglik,brier,disSlope,auc,intSlope,binperformance)
}