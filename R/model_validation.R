#' Validation of binary classification rule
#'
#' @param model Model's formula
#' @param data Data to fit and assess model performance
#' @param type Type of statistical model
#' @param scope Scope in stepwise variable selection
#' @param validation_method "cv" (cross-validation) or "bootstrap"
#' @param B number of bootstrap samples / B-fold cross-validation
#' @param B_mat Data matrix for cross-validation
#' @param test_set_method How is the test set defined (either "all" -> use original dataset as test data (default for bootstrap), or "excluded" -> use patients excluded from cv/bootstrap sample as test data)
#' @param show_selected_coef If TRUE, show the selected coefficients (only for glm's with model selection)
#' @param criteria
#' @param roc If TRUE, true positive rates corresponding to false positive rates given in argument fpr are also evaluated and cv/bootstrap corrected
#' @param fpr
#' @param ... Other parameters for fit_method and fit_response
#'
#' @return List of validation output
#' @export
classValidation <- function(model, data, type, scope=NULL,
                             validation_method = "bootstrap", B = NULL, B_mat = NULL,
                             test_set_method = ifelse(validation_method == "bootstrap", "all", "excluded"),
                             show_selected_coef = F, criteria = c(Brier, c_index, calib),
                             roc = F, fpr = seq(0, 1, length = 101), ...){

  n <- nrow(data)
  #-- checks
  if (is.null(B) && is.null(B_mat)) stop("Either argument B or B_mat must be provided")
  if ((!is.null(B)) && (!is.null(B_mat)) && (ncol(B_mat) != B)) stop("Arguments B and B_mat incompatible")
  #-- Create B and B_mat (containing selected observations in each bootstrap/cv sample)
  if (is.null(B)) B <- ncol(B_mat)
  if (is.null(B_mat)){
    if (validation_method == "bootstrap") B_mat <- rmultinom(B, size = n, prob = rep(1,n)/n)
    if (validation_method == "cv"){
      which.group <- sample(rep(1:B, ceiling(n/B))[1:n], n, replace = F)
      B_mat <- matrix(NA, nrow = n, ncol = B)
      for (i in 1:B) B_mat[,i] <- as.numeric(which.group != i)
    }
  }

  #-- function to calculate true positive rate at points fpr (corresponding false positive rate)
  true_positive_rate <- function(predicted, y, fpr){
    perf <- ROCR::performance(ROCR::prediction(predicted,y), "tpr", "fpr")
    tpr <- rep(NA, length(fpr))
    for (i in 1:length(tpr)) tpr[i] <- (perf@y.values[[1]])[sum(perf@x.values[[1]] <= fpr[i])]
    tpr
  }
  #-- set up dummy data structures to store accuracy indices, selected coefficients, etc.
  # selected variables on resampled data
  if (show_selected_coef == T){
    names.allcov <- colnames(model.matrix(model, data = data)) #includes intercept (is removed at the end)
    num.allcov   <- length(names.allcov)
    selectedVars <- matrix(0, ncol = num.allcov, nrow = B)
    colnames(selectedVars) <- names.allcov
  }
  # true positive rates on resampled data
  tpr.test <- tpr.train <- matrix(NA, nrow = B, ncol = length(fpr))
  #-- Indices and tpr on original data
  #browser()
  fitted.orig <- fit.response(fit.method(model, data, type, ...), data, type, ...)
  y.orig <- as.numeric(model.response(model.frame(update(model,.~1), data)))
  # (note: "update(model,.~1)" instead of "model" is necessary, because model.frame gives an error otherwise for gam's with model terms s(.)
  perf.orig <- myaccuracy(fitted.orig, y.orig, criteria)

  # aggrageted indices
  indices <- matrix(nrow = length(perf.orig), ncol = 5)
  rownames(indices) <- names(perf.orig)
  colnames(indices) <- c("index.orig", "training", "test", "optimism", "index.corrected")
  indices[,1] <- perf.orig
  if (roc) tpr.orig <- true_positive_rate(as.numeric(fitted.orig), y.orig, fpr)
  # indices for resampled data
  #browser()
  indices.resample <- matrix(ncol = length(perf.orig)*2, nrow = B)
  colnames(indices.resample) <- c(paste(names(perf.orig), "train"),
                                  paste(names(perf.orig), "test"))

  #-- Start the actual valdation
  cat("Start: \n")
  for (i in 1:B){
    cat(i, "\r")
    # datasets and response
    trainData <- data[rep(1:n, B_mat[,i]), ]
    if (test_set_method == "excluded") testData  <- data[B_mat[,i] == 0, ]
    if (test_set_method == "all") testData  <- data
    y.train <- as.numeric(model.response(model.frame(update(model, .~1), trainData)))
    y.test <- as.numeric(model.response(model.frame(update(model, .~1), testData)))
    # Predicted responses
    fitModel <- fit.method(model, data= trainData, type, ...)
    fitted.risk.train <- fit.response(fitModel, data = trainData, type, ...)
    fitted.risk.test <- fit.response(fitModel, data = testData, type, ...)
    # accuracy measures
    indices.resample[i, 1:length(perf.orig)] <- myaccuracy(fitted.risk.train, y.train, criteria)
    indices.resample[i, (length(perf.orig)+1):(2*length(perf.orig))] <- myaccuracy(fitted.risk.test, y.test, criteria)
    if (show_selected_coef == T){
      selectedVars[i, colnames(model.matrix(fitModel$formula, data = data))] <- 1
    }
    if (roc){
      tpr.train[i,] <- true_positive_rate(as.numeric(fitted.risk.train), y.train, fpr)
      tpr.test[i,] <-  true_positive_rate(as.numeric(fitted.risk.test), y.test, fpr)
    }
  }
  #browser()
  #-- aggregate bootstrap/cv indices
  summary.result <- apply(indices.resample, 2, function(x) mean(x, na.rm = TRUE))
  indices[, "training"] <- summary.result[1:length(perf.orig)]
  indices[, "test"] <- summary.result[(length(perf.orig)+1):(2*length(perf.orig))]
  indices[, "optimism"] <- indices[, "training"] - indices[, "test"]
  if (validation_method == "cv") indices[, "optimism"] <- indices[, "index.orig"] - indices[, "test"]
  indices[, "index.corrected"] <- indices[, "index.orig"] - indices[, "optimism"]


  #-- get estimate of "corrected" tpr
  if (roc){
    tpr.corrected <- tpr.orig - (apply(tpr.train, 2, mean) - apply(tpr.test, 2, mean))
    # "correct" tpr.corrected to be between 0 and 1 and monotonically increasing
    tpr.corrected <- pmin(pmax(tpr.corrected, 0), 1)
    for (i in (2:length(tpr.corrected))){
      if (tpr.corrected[i] < tpr.corrected[i-1]) tpr.corrected[i] <- tpr.corrected[i-1]
    }
  }
  #-- final result
  result <- list(indices = indices, indices.resample = indices.resample)
  #add selected coefficients
  if (show_selected_coef == T) {
    if (colnames(selectedVars)[1] == "(Intercept)") selectedVars <- selectedVars[, -1]
    numvar <- rowSums(selectedVars)
    selectedVars_propChosen <- c(apply(selectedVars,2,mean), "mean.numberOfVariables" = mean(numvar), "sd.numberOfVariables" = sd(numvar))
    selectedVars <- cbind(selectedVars, "numberOfVariables" = numvar)
    result$selectedVars_propChosen <- selectedVars_propChosen
    result$selectedVars.resample <- selectedVars
  }
  # add ROC information
  if (roc){
    result$ROC <- list(fpr = fpr, tpr.orig = tpr.orig, tpr.corrected = tpr.corrected, tpr.train = tpr.train, tpr.test = tpr.test)
  }
  # return result
  result
}