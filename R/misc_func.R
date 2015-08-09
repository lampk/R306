#' @export
# Matrix for cross-validation
cv.mat <- function(n, B = 10, times = 10) {
  B_mat <- NULL
  for (i in (1:times)){
    folds <- sample(rep(1:B, length = n))
    tmp <- sapply(1:B, function(x){as.numeric(folds != x)})
    B_mat <- cbind(B_mat, tmp)
  }
  return(B_mat)
}

gety <- function(model, data, type){
  if (type %in% c("coxtv", "coxtv.step.AIC", "coxtv.step.BIC")){
    y <- model.response(model.frame(model, data))
  } else {
    y <- as.numeric(model.response(model.frame(update(model,.~1), data)))
  }
  return(y)
}