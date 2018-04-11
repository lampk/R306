#' Methods to fit models
#'
#' @param model Model formula to fit
#' @param data Data to fit model
#' @param ... Other arguments
#'
#' @return fitted model
#' @export
fit.method <- function(model, data, type = c('lm',
                                             'glm', 'glm.step.AIC', 'glm.step.BIC',
                                             'glm.best.AIC', 'glm.best.BIC',
                                             'gam',
                                             'glmnet', 'glmnet.1se', 'glmnet.cv', 'glmnet.boot',
                                             "adaptiveLASSO",
                                             "SCAD",
                                             'penalized',
                                             'stability',
                                             'stability.step.AIC',
                                             'cox', 'cox.step.AIC', 'cox.step.BIC',
                                             'cox.best.AIC', 'cox.best.BIC',
                                             'glmnet.cox1', 'glmnet.cox2',
                                             'rpart1', 'rpart2',
                                             'randomForest',
                                             'gbm', 'gbm.cox',
                                             'coxtv', 'coxtv.step.AIC', 'coxtv.step.BIC'), 
                       pi = 0.8, rep = 100, size = 0.5, seed = NULL, ...){
  fit <- switch(type,
                lm               = fit.method.lm(model, data, ...),
                glm              = fit.method.glm(model, data, ...),
                glm.step.AIC     = fit.method.glm.step(model, data, ...),
                glm.step.BIC     = fit.method.glm.step(model, data, k = log(nrow(data)), ...),
                glm.best.AIC     = fit.method.glm.glmulti(model, data, crit = "aic", ...),
                glm.best.BIC     = fit.method.glm.glmulti(model, data, crit = "bic", ...),
                gam              = fit.method.gam(model, data, ...),
                glmnet           = fit.method.glmnet(model, data, ...),
                glmnet.1se       = fit.method.glmnet(model, data, ...),
                glmnet.cv        = fit.method.glmnet.cv(model, data, ...),
                glmnet.boot      = fit.method.glmnet.boot(model, data, ...),
                adaptiveLASSO    = fit.method.alasso(model, data, ...),
                SCAD             = fit.method.scad(model, data, ...),
                penalized        = fit.method.penalized(model, data, ...),
                stability        = fit.method.stability(model, data, ...),
                stability.step.AIC = fit.method.stab.step.AIC(model, data, pi = pi, B = rep, size = size, seed = seed, parallel = parallel, ...),
                cox              = fit.method.cox(model, data, ...),
                cox.step.AIC     = fit.method.cox.step(model, data, ...),
                cox.step.BIC     = fit.method.cox.step(model, data, k = log(nrow(data)), ...),
                cox.best.AIC     = fit.method.cox.glmulti(model, data, crit = "aic", ...),
                cox.best.BIC     = fit.method.cox.glmulti(model, data, crit = "bic", ...),
                glmnet.cox1      = fit.method.glmnet(model, data, family = "cox", ...),
                glmnet.cox2      = fit.method.glmnet(model, data, family = "cox", ...),
                rpart1           = fit.method.rpart1(model, data, ...),
                rpart2           = fit.method.rpart2(model, data, ...),
                randomForest     = fit.method.randomForest(model, data, ...),
                gbm              = fit.method.gbm(model, data, ...),
                gbm.cox          = fit.method.gbm(model, data, distribution = "coxph", ...),
                coxtv            = fit.method.coxtv(model, data, ...),
                coxtv.step.AIC   = fit.method.coxtv.step(model, data, ...),
                coxtv.step.BIC   = fit.method.coxtv.step(model, data, k = log(nrow(data)), ...))
  fit
}

#' @describeIn fit.method Linear regression
fit.method.lm <- function(model, data, ...) {
  lm(model, data)
}

#' @describeIn fit.method Linear regression with stepwise variable selection using AIC as stopping rule
fit.method.lm.step.AIC <- function(model, data, ...) {
  step(lm(model, data), trace = 0)
}

#' @describeIn fit.method Linear regression with stepwise variable selection using BIC as stopping rule
fit.method.lm.step.BIC <- function(model, data, ...) {
  step(lm(model, data), trace = 0, k = log(nrow(data)))
}


#' @describeIn fit.method Logistic regression
fit.method.glm <- function(model, data, ...) {
  glm(model, data, family = "binomial")
}

#' @describeIn fit.method Logistic regression with stepwise variable selection using AIC as stopping rule
fit.method.glm.step.AIC <- function(model, data, ...) {
  step(glm(model, data, family = "binomial"), trace = 0, ...)
}

#' @describeIn fit.method Logistic regression with stepwise variable selection using AIC as stopping rule + stability selection
fit.method.stab.step.AIC <- function(model, data, pi = 0.8, B = 100, size = 0.5, seed = NULL, parallel = FALSE, ...) {
  # stability selection
  n <- nrow(data)
  if (!is.null(seed)) set.seed(seed)
  id <- lapply(1:B, function(x) sample.int(n = n, size = round(size * n)))
  #browser()
  if (parallel) {
    require(foreach)
    require(doParallel)
    registerDoParallel(detectCores() - 1)
    res <- foreach(i = 1:B, .combine = list, .multicombine = TRUE)  %dopar%  {
      datai <- data[id[[i]], ]
      fiti <- fit.method.glm.step.AIC(model = model, data = datai, ...)
      all.vars(formula(fiti))[-1]
    }
    stopImplicitCluster()
  } else {
    res <- vector("list", B)
    ## get selected variables
    for (i in (1:B)) {
      cat("\r Step :", i)
      datai <- data[id[[i]], ]
      fiti <- fit.method.glm.step.AIC(model = model, data = datai, ...)
      res[[i]] <- all.vars(formula(fiti))[-1]
    }
  }
  
  #browser()
  ## derive output
  phi <- table(unlist(res))/B
  if (sum(phi >= pi) == 0) {
    selected <- NULL
  } else {selected <- names(phi[phi >= pi])}
  if (is.null(selected)) cat("No variable was selected.\n")
  
  ## re-fit with glm
  if (is.null(selected)) {
    tmp <- 1
  } else {tmp <- paste(selected, collapse = "+")}
  return(fit.method.glm(model = update(model, paste(". ~", tmp)), data, ...))
}

#' @describeIn fit.method Logistic regression with stepwise variable selection using BIC as stopping rule
fit.method.glm.step.BIC <- function(model, data, ...) {
  step(glm(model, data, family = "binomial"), trace = 0, k = log(nrow(data)), ...)
}

#' @describeIn fit.method Logistic regression with stepwise variable selection
fit.method.glm.step <- function(model, data, ...) {
  step(glm(model, data, family = "binomial"), trace = 0, ...)
}

#' @describeIn fit.method Cox regression
fit.method.cox <- function(model, data, ...) {
  survival::coxph(model, data, ...)
}

#' @describeIn fit.method Cox regression with stepwise variable selection using AIC as stopping rule
fit.method.cox.step.AIC <- function(model, data, ...) {
  step(survival::coxph(model, data), trace = 0)
}

#' @describeIn fit.method Cox regression with stepwise variable selection using BIC as stopping rule
fit.method.cox.step.BIC <- function(model, data, ...) {
  step(survival::coxph(model, data), trace = 0, k = log(n))
}

#' @describeIn fit.method Cox regression with stepwise variable selection
fit.method.cox.step <- function(model, data, ...) {
  fit <- do.call("survival::coxph", list(model, data))
  step(fit, trace = 0, ...)
}

#' @describeIn fit.method Cox regression with time-dependent variable
fit.method.coxtv <- function(model, data, ...) {
  # get response
  mf <- model.frame(formula = model, data = data)
  Y <- model.extract(mf, "response")
  # get names of response variables (time and status)
  resvar <- as.character(model[[2]])[-1]
  # create new dataset with many rows for each observation
  tmp <- survival::survSplit(data = data, cut = sort(unique(Y[Y[, "status"] == 1, "time"])),
                   end = resvar[1], event = resvar[2], start = "Tstart")
  # update call of outcome
  res <- model[[2]]
  model[[2]][2] <- expression(Tstart); model[[2]][3] <- res[2]; model[[2]][4] <- res[3]
  # fit
  survival::coxph(model, tmp, model = TRUE)
}

#' @describeIn fit.method Cox regression with time-dependent variable using stepwise variable selection
fit.method.coxtv.step <- function(model, data, ...) {
  # get response
  mf <- model.frame(formula = model, data = data)
  Y <- model.extract(mf, "response")
  # get names of response variables (time and status)
  resvar <- as.character(model[[2]])[-1]
  # create new dataset with many rows for each observation
  tmp <- survival::survSplit(data = data, cut = sort(unique(Y[Y[, "status"] == 1, "time"])),
                   end = resvar[1], event = resvar[2], start = "Tstart")
  # update call of outcome
  res <- model[[2]]
  model[[2]][2] <- expression(Tstart); model[[2]][3] <- res[2]; model[[2]][4] <- res[3]

  # fit
  fit <- do.call("survival::coxph", list(model, data = tmp, model = TRUE), envir = environment())
  step(fit, trace = FALSE, ...)
}

#' @describeIn fit.method Logistic regression with LASSO (using glmnet)
fit.method.glmnet <- function(model, data,
                              family = "binomial", lambda = NULL, nfolds = 10,
                              type.measure = "deviance", standardize = TRUE, alpha = 1, ...){
  #!!! be careful: not work with categorical variables, for that consider 'grouped lasso'

  # define response and covariate
  x <- model.matrix(model, data)[, -1]
  if (family == "cox"){
    y <- eval(parse(text = paste("with(data, as.matrix(", model[2], "))")))
  } else {y <- unlist(model.frame(model, data)[1])}

  # fit glmnet model using 10-fold cross-validation to choose the best tuning parameter
  out <- glmnet::cv.glmnet(x, y, family = family, lambda = lambda, nfolds = nfolds, type.measure = type.measure, standardize = standardize, alpha = alpha, ...)
  out$model <- model
  out
}

#' @describeIn fit.method Logistic regression with LASSO (using glmnet) and refit (only for logistic regression)
fit.method.glmnet.cv <- function (model, data, family = "binomial", lambda = NULL, refit = TRUE, nfolds = 10, 
                                  type.measure = "deviance", standardize = TRUE, alpha = 1, 
                                  krepeat = 10, nse = 1, nlambda = 100, lambda.min.ratio = 0.0001, parallel = FALSE, ...) {
  #!!! be careful: not work with categorical variables, for that consider 'grouped lasso'
  
  ## function to match varnames
  match.varname <- function(varname, model) {
    model.vars <- all.vars(model)
    model.x <- attr(terms(model), "term.labels")
    tmp1 <- sapply(varname, function(x) {paste(as.numeric(sapply(model.vars, function(y) grepl(y, x))), collapse = "")})
    tmp2 <- sapply(model.x, function(x) {paste(as.numeric(sapply(model.vars, function(y) grepl(y, x))), collapse = "")})
    out <- names(tmp2)[tmp2 %in% tmp1]
    return(out)
  }
  
  # define response and covariate
  x <- model.matrix(model, data)[, -1]
  y <- unlist(model.frame(model, data)[1])
  
  ## pre-specify lambda
  lambda.max <- glmnet::glmnet(x, y, family = family, standardize = standardize, alpha = alpha, nlambda = nlambda, ...)$lambda[1]
  lambda <- seq(from = lambda.max, to = lambda.max * lambda.min.ratio, length = nlambda)
  
  # fit glmnet model using k-fold cross-validation to choose the best tuning parameter
  ## to reduce randomness: https://stats.stackexchange.com/questions/97777/variablity-in-cv-glmnet-results
  glmnet.fit <- glmnet::cv.glmnet(x, y, family = family, 
                                  lambda = lambda, nfolds = nfolds, type.measure = type.measure, 
                                  standardize = standardize, alpha = alpha, ...)
  if (parallel) {
    require(doParallel)
    cl <- makeCluster(detectCores())
    registerDoParallel(cl)
    err <- foreach(i = 1:krepeat,.combine = cbind) %dopar% {
      glmnet.fit <- glmnet::cv.glmnet(x, y, family = family, 
                                      lambda = lambda, nfolds = nfolds, type.measure = type.measure, 
                                      standardize = standardize, alpha = alpha, ...)
      cvmi <- data.frame(cvm = glmnet.fit$cvm)
      names(cvmi) <- paste0("cvm", i)
      return(cvmi)
    }
    stopCluster(cl)
    cvms <- cbind(lambda = glmnet.fit$lambda, err)
  } else {
    cvms <- NULL
    for (i in 1:krepeat) {
      cat("\r", i)
      glmnet.fiti <- glmnet::cv.glmnet(x, y, family = family, 
                                      lambda = lambda, nfolds = nfolds, type.measure = type.measure, 
                                      standardize = standardize, alpha = alpha, ...)
      cvmi <- data.frame(lambda = glmnet.fiti$lambda, cvm = glmnet.fit$cvm)
      names(cvmi)[2] <- paste0("cvm", i)
      if (is.null(cvms)) {
        cvms <- cvmi
      }
      else {
        cvms <- merge(cvms, cvmi, by = "lambda", all = TRUE)
      }
    }
  }
  
  cvms <- merge(data.frame(lambda = glmnet.fit$lambda, cvsd = glmnet.fit$cvsd), 
                cvms[cvms$lambda %in% glmnet.fit$lambda, ],
                by = "lambda", all = TRUE)
  cvmmeans <- rowMeans(cvms[, -c(1:2), drop = FALSE], na.rm = TRUE)
  cvmin <- min(cvmmeans, na.rm = TRUE)
  idmin <- cvmmeans <= cvmin
  lambda.min <- max(cvms$lambda[idmin], na.rm = TRUE)
  semin <- (cvmmeans + nse * cvms$cvsd)[idmin]
  idmin <- cvmmeans <= semin
  lambda.final <- max(cvms$lambda[idmin], na.rm = TRUE)
  
  if (refit) {
    coefs <- coef(glmnet.fit, s = lambda.final)
    
    # get selected variables
    selected.vars <- match.varname(varname = coefs@Dimnames[[1]][coefs@i + 1], model = model)
    selected.vars <- gsub(pattern = ":", replacement = "*", x = selected.vars)
    new.model <- update(model, reformulate(selected.vars, response = "."))
    # refit model
    out <- glm(formula = new.model, data = data, family = "binomial")
  } else {
    out <- glmnet.fit
    out$model <- model
    out$s <- lambda.final
  }
  return(out)
}

#' @describeIn fit.method Logistic regression with LASSO (using glmnet) and refit (only for logistic regression) ## use bootstrap
fit.method.glmnet.boot <- function (model, data, family = "binomial", lambda = NULL, refit = TRUE, krepeat = 100, 
                                    type.measure = "deviance", standardize = TRUE, alpha = 1, 
                                    nse = 1, nlambda = 100, lambda.min.ratio = 0.0001, parallel = FALSE, ...) {
  #!!! be careful: not work with categorical variables, for that consider 'grouped lasso'
  
  ## function to match varnames
  match.varname <- function(varname, model) {
    model.vars <- all.vars(model)
    model.x <- attr(terms(model), "term.labels")
    tmp1 <- sapply(varname, function(x) {paste(as.numeric(sapply(model.vars, function(y) grepl(y, x))), collapse = "")})
    tmp2 <- sapply(model.x, function(x) {paste(as.numeric(sapply(model.vars, function(y) grepl(y, x))), collapse = "")})
    out <- names(tmp2)[tmp2 %in% tmp1]
    return(out)
  }
  
  get.deviance <- function(y, p) {
    ## modify y
    nc <- dim(y)
    if (is.null(nc)) {
      y <- as.factor(y)
      ntab <- table(y)
      nc <- as.integer(length(ntab))
      y <- diag(nc)[as.numeric(y), ]
    }
    ## modify p
    p <- pmin(pmax(p, 1e-05), 1 - 1e-05)
    
    ## get deviance
    lp <- y[, 1] * log(1 - p) + y[, 2] * log(p)
    ly <- log(y)
    ly[y == 0] <- 0
    ly <- drop((y * ly) %*% c(1, 1))
    return(2 * (ly - lp))
  }
  
  # define response and covariate
  x <- model.matrix(model, data)[, -1]
  y <- unlist(model.frame(model, data)[1])
  
  # fit glmnet model using bootstrap to choose the best tuning parameter
  m <- nrow(data)
  boot_mat <- rmultinom(krepeat, size = m, prob = rep(1, m)/m)
  
  ## pre-specify lambda
  lambda.max <- glmnet::glmnet(x, y, family = family, standardize = standardize, alpha = alpha, nlambda = nlambda, ...)$lambda[1]
  lambda <- seq(from = lambda.max, to = lambda.max * lambda.min.ratio, length = nlambda)
  
  if (parallel) {
    require(doParallel)
    cl <- makeCluster(detectCores())
    registerDoParallel(cl)
    err <- foreach(i = 1:krepeat, .combine = cbind) %dopar% {
      x.train <- x[rep(1:m, boot_mat[,i]), , drop = FALSE]
      y.train <- y[rep(1:m, boot_mat[,i])]
      glmnet.fit <- glmnet::glmnet(x.train, y.train, family = family, standardize = standardize, alpha = alpha, lambda = lambda, ...)
      sapply(lambda, function(l) {
        p <- plogis(glmnet::predict.glmnet(glmnet.fit, newx = x, s = l, type = "response"))
        sum(get.deviance(y = y, p = p))
      })
    }
    stopCluster(cl)
    err <- cbind(lambda = lambda, err)
  } else {
    err <- cbind(lambda = lambda, matrix(ncol = krepeat, nrow = length(lambda)))
    for (i in 1:krepeat) {
      cat("\r", i)
      ## fit
      x.train <- x[rep(1:m, boot_mat[,i]), , drop = FALSE]
      y.train <- y[rep(1:m, boot_mat[,i])]
      glmnet.fit <- glmnet::glmnet(x.train, y.train, family = family, standardize = standardize, alpha = alpha, lambda = lambda, ...)
      err[, i + 1] <- sapply(lambda, function(l) {
        p <- plogis(glmnet::predict.glmnet(glmnet.fit, newx = x, s = l, type = "response"))
        sum(get.deviance(y = y, p = p))
      })
    }
  }
  
  #browser()
  errmeans <- rowMeans(err[, -1, drop = FALSE], na.rm = TRUE)
  errse <- apply(err[, -1, drop = FALSE], 1, function(x){
    sd(x)/sqrt(krepeat)
  })
  
  errmin <- min(errmeans, na.rm = TRUE)
  idmin <- errmeans <= errmin
  lambda.min <- max(lambda[idmin], na.rm = TRUE)
  semin <- (errmeans + nse * errse)[idmin]
  idmin <- errmeans <= semin
  lambda.final <- max(lambda[idmin], na.rm = TRUE)
  
  glmnet.fit <- glmnet::glmnet(x, y, family = family, standardize = standardize, alpha = alpha, lambda = lambda, ...)
  
  if (refit) {
    coefs <- coef(glmnet.fit, s = lambda.final)
    
    # get selected variables
    selected.vars <- match.varname(varname = coefs@Dimnames[[1]][coefs@i + 1], model = model)
    selected.vars <- gsub(pattern = ":", replacement = "*", x = selected.vars)
    new.model <- update(model, reformulate(selected.vars, response = "."))
    # refit model
    out <- glm(formula = new.model, data = data, family = "binomial")
  } else {
    out <- glmnet.fit
    out$model <- model
    out$s <- lambda.final
  }
  return(out)
}

#' @describeIn fit.method Logistic regression with LASSO (using penalized)
fit.method.penalized <- function(model, data, ...) {
  opt <-  penalized::optL1(model, data = data, model = 'logistic', fold = 10, standardize = TRUE, trace = FALSE)
  penalized::penalized(model, data = data, model = 'logistic', lambda1 = opt$lambda, standardize = TRUE, trace = FALSE)
}

#' @describeIn fit.method Logistic regression with adaptive LASSO
fit.method.alasso <- function(model, data,
                              family = "binomial", lambda = NULL, nfolds = 10, type.measure = "deviance", standardize = TRUE, alpha = 1,
                              ...) {
  # get coefficient from logistic regression
  betas <- glm(model, family = family, data = data)$coef[-1]
  # fit lasso
  ## define response and covariate
  x <- model.matrix(model, data)[, -1]
  if (family == "cox"){
    y <- eval(parse(text = paste("with(data, as.matrix(", model[2], "))")))
  } else {y <- unlist(model.frame(model, data)[1])}

  ## fit glmnet model using 10-fold cross-validation to choose the best tuning parameter
  out <- glmnet::cv.glmnet(x, y, family = family, lambda = lambda, nfolds = nfolds, type.measure = type.measure, standardize = standardize, alpha = alpha, penalty.factor = 1/abs(betas), ...)
  out$model <- model
  out
}

#' @describeIn fit.method Logistic regression with SCAD
fit.method.scad <- function(model, data,
                            family = "binomial", penalty = "SCAD", gamma = 3.7, alpha = 1, lambda.min = .001, nlambda = 100) {
  ## get data in
  x <- model.matrix(model, data)[, -1]
  if (family == "cox"){
    y <- eval(parse(text = paste("with(data, as.matrix(", model[2], "))")))
  } else {y <- unlist(model.frame(model, data)[1])}

  ## fit SCAD
  out <- ncvreg::cv.ncvreg(x, y, family = family, penalty = penalty, gamma = gamma, alpha = alpha, lambda.min = lambda.min, nlambda = nlambda, nfolds = 10)
  out$model <- model
  out
}

#' @describeIn fit.method GLM LASSO with stability selection
fit.method.stability <- function(model, data,
                                        size = 0.632, steps = 100, weakness = 1, error = 0.05, pi_thr = 0.6, error.type = "pfer",
                                        family = "binomial", standardize = TRUE, intercept = TRUE, alpha = 1,
                                        ...){
  # define response and covariate
  x <- model.matrix(model, data)[, -1]
  y <- unlist(model.frame(model, data)[1])

  # stability selection
  stability_path <- c060::stabpath(y, x, size = size, steps = steps, weakness = weakness, family = family, standardize = standardize, intercept = intercept, alpha = alpha, ...)
  idx <- c060::stabsel(stability_path, error = error, type = error.type, pi_thr = pi_thr)$stable

  # fit a glm with variable selected by stability selection
  if (is.character(family))
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family))
    family <- family()

  if (length(idx) == 0){
    fit <- glm(y ~ 1, family = family)
  } else {
    fit <- glm.fit(x = model.matrix(model, data)[, c(1, idx + 1)], y = y, family = family)
    class(fit) <- c(fit$class, c("glm", "lm"))
  }
  fit
}

#' @describeIn fit.method GAM Logistic regression
fit.method.gam <- function(model, data, family = "binomial", ...) {
  mgcv::gam(model, data, family = family, ...)
}

#' @describeIn fit.method CART: build and prune tree at the same time
fit.method.rpart1 <- function(model, data, ...) {
  rpart::rpart(model, data, method = "anova", ...)
  }

#' @describeIn fit.method CART: build, then prune tree
fit.method.rpart2 <- function(model, data,
                              method = "anova", control = rpart.control(minsplit = 20, xval = 10, cp = 0),
                              ...) {
  # fit large tree
  tree <- rpart::rpart(model, data, method = method, control = control, ...)
  # find optimal complexity parameter (according to 1-std rule)
  index <- which.min(tree$cptable[, "xerror"])
  minrange <- tree$cptable[index, "xerror"] + c(-1, 1) * tree$cptable[index, "xstd"]
  cp.opt <- with(subset(as.data.frame(tree$cptable), xerror <= max(minrange) & xerror >= min(minrange)), max(CP))
  # prune tree with optimal complexity parameter
  rpart::prune.rpart(tree, cp = cp.opt)
}

#' @describeIn fit.method Random forest
fit.method.randomForest <- function(model, data, ntree = 500, ...) {
  randomForest::randomForest(formula = model, data = data, ntree = ntree, ...)
}

#' @describeIn fit.method GLM Boosting
fit.method.gbm <- function(model, data,
                           distribution = "bernoulli", cv.folds = 10, n.trees = 3000, interaction.depth = 2, verbose = FALSE, shrinkage = 0.001, n.cores = 1, ...) {
  gbm::gbm(model, data = data, distribution = distribution, verbose = verbose, cv.folds = cv.folds, n.trees = n.trees, interaction.depth = interaction.depth, shrinkage = shrinkage, n.cores = n.cores, ...)
}

#' @describeIn fit.method Logistic regression with best subset using AIC as criterion
fit.method.glm.glmulti.AIC <- function(model, data, ...) {
  bestmodel <- as.formula(paste(summary(glmulti::glmulti(model, data = data,
                                                         fitfunc = glm, crit = aic,
                                                         marginality = T,
                                                         plotty = F, report = F, ...))$bestmodel,
                                collapse = ""))
  glm(bestmodel, data, family = "binomial")
}

#' @describeIn fit.method Logistic regression with best subset using BIC as criterion
fit.method.glm.glmulti.BIC <- function(model, data, ...) {
  bestmodel <- as.formula(paste(summary(glmulti::glmulti(model, data = data,
                                                         fitfunc = glm, crit = bic,
                                                         marginality = T,
                                                         plotty = F, report = F, ...))$bestmodel,
                                collapse = ""))
  glm(bestmodel, data, family = "binomial")
}

#' @describeIn fit.method Cox regression with best subset using AIC as criterion
fit.method.cox.glmulti.AIC <- function(model, data, ...) {
  bestmodel <- as.formula(paste(summary(glmulti::glmulti(model, data = data,
                                                         fitfunc = coxph, crit = aic,
                                                         marginality = T,
                                                         plotty = F, report = F, ...))$bestmodel,
                                collapse = ""))
  survival::coxph(formula = bestmodel, data = data)
}

#' @describeIn fit.method Cox regression with best subset using BIC as criterion
fit.method.cox.glmulti.BIC <- function(model, data, ...) {
  bestmodel <- as.formula(paste(summary(glmulti::glmulti(model, data = data,
                                                         fitfunc = coxph, crit = bic,
                                                         marginality = T,
                                                         plotty = F, report = F, ...))$bestmodel,
                                collapse = ""))
  survival::coxph(formula = bestmodel, data = data)
}

#' @describeIn fit.method Logistic regression with best subset
fit.method.glm.glmulti <- function(model, data, ...) {
  fit <- glmulti::glmulti(model, data = data, fitfunc = glm, marginality = T, plotty = F, report = F, ...)
  bestmodel <- as.formula(paste(summary(fit)$bestmodel, collapse = ""))
  glm(bestmodel, data, family = "binomial")
}

#' @describeIn fit.method Cox regression with best subset
fit.method.cox.glmulti <- function(model, data, ...) {
  fit <- glmulti::glmulti(model, data = data, fitfunc = coxph, marginality = T, plotty = F, report = F, ...)
  bestmodel <- as.formula(paste(summary(fit)$bestmodel, collapse = ""))
  survival::coxph(formula = bestmodel, data = data)
}

