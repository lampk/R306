#' Summarise dataset
#'
#' @param formula Specification of dependent and independent variables
#' @param data Data to summarise
#' @param pooledGroup Pooling all subgroups of dependent variable?
#' @param contSummary Summary statistics for continuous variables?
#' @param caption Caption of final kable table
#' @param kable Create kable table
#' @param test Display p values from statistical test comparing independent variables between sub-group of dependent variable?
#' @param continuous Independent variables are continuous variables?
#' @param digits Number of decimal digit number
#'
#' @return a matrix

mySummary.onevar <- function(varname, variable, group = NULL, continuous = NA, contSummary = "med.IQR", 
                             test = FALSE, digits = 1, pcutoff = 0.0001, 
                             chisq.test = FALSE, correct = FALSE, workspace = 1000000,
                             hybrid = FALSE, simulate.p.value = FALSE, B = 2000){
  if (is.na(continuous)) continuous <- ifelse(is.factor(variable) | length(unique(na.omit(variable))) <= 5, FALSE, TRUE)

  mycont.summary <- function(variable,group,test,digits) {
    if (is.null(group)) {
      ngroup <- 1
      summarystat.nice <- cont.Summary(unclass(variable), contSummary = contSummary, digits = digits, n = FALSE)
      n <- length(na.omit(variable))
    } else {
      ngroup <- length(levels(group))
      summarystat.nice <- by(unclass(variable), group, cont.Summary, contSummary = contSummary, digits = digits, n = FALSE)
      n <- c(by(variable, group, function(x) length(na.omit(x))))
    }

    result <- matrix("", ncol = ngroup * 2 + 1, nrow = 1)
    result[1, seq(2, ncol(result), by = 2)] <- n
    result[1, seq(3, ncol(result), by = 2)] <- unlist(summarystat.nice)
    if (test == TRUE & !is.null(group)) {
      # overall Kruskal-Wallis test for group differences
      pval <- myformat.pval(kruskal.test(variable ~ group)$p.value, cutoff = pcutoff)
      result <- cbind(result, pval)
    }
    result
  }

  mycat.summary <- function(variable, group, test, digits, chisq.test, correct, workspace, hybrid, simulate.p.value, B) {
    if (is.null(group)) {
      ngroup <- 1
      ta <- table(variable)
      ta.prop <- ta/sum(ta)
      dim(ta) <- c(ngroup, length(ta))
      colnames(ta) <- names(table(variable))
    } else {
      ngroup <- length(levels(group))
      ta <- table(group, variable)
      ta.prop <- unclass(ta/apply(ta, 1, sum))
    }

    ta.nice <- matrix(paste(ta," (", formatC(100*unclass(ta.prop), digits, format = "f"), "%", ")", sep = ""),
                      nrow = nrow(ta), ncol = ncol(ta))
    result <- matrix("", ncol = ngroup * 2 + 1, nrow = ncol(ta) + 1)
    result[2:nrow(result), 1] <- paste("- ", colnames(ta))
    result[2:nrow(result), seq(3, ncol(result), by = 2)] <- t(ta.nice)
    result[1, seq(2, ncol(result), by = 2)] <- apply(ta, 1, sum) # n's
    if (test) {
      # Fisher's exact test for group differences
      pval <- switch(chisq.test,
                     TRUE = myformat.pval(chisq.test(ta, correct = correct, 
                                                     simulate.p.value = simulate.p.value, B = B)$p.value, cutoff = pcutoff),
                     FALSE = myformat.pval(fisher.test(ta, workspace = workspace, hybrid = hybrid, 
                                        simulate.p.value = simulate.p.value, B = B)$p.value, cutoff = pcutoff))
      result <- cbind(result, "")
      result[1,ngroup * 2 + 2] <- pval
    }
    result
  }

  if (continuous) r <- mycont.summary(variable, group, test, digits)
  else r <- mycat.summary(variable, group, test, digits, chisq.test, correct, workspace, hybrid, simulate.p.value, B)
  r[1, 1] <- varname
  r
}

#' @export
mySummary.allvar <- function(formula, data, pooledGroup = FALSE, contSummary = "med.IQR",
                             caption = NULL, kable = FALSE, test = FALSE, continuous = NA,
                             digits = 1, pcutoff = 0.0001, chisq.test = FALSE, correct = FALSE, workspace = 1000000,
                             hybrid = FALSE, simulate.p.value = FALSE, B = 2000){
  # contSummary can be median (90% range) "med.90" or median (IQR) "med.IQR" or median (range) "med.range" or "mean.sd"

  if (pooledGroup&test){
    warning("Display of both pooled groups and test currently not implemented. Tests not displayed.")
    test <- FALSE
  }

  dat <- model.frame(formula, data = data, na.action = NULL)
  if (length(formula)==2) {blvars <- dat; group <- factor(rep("All patients",nrow(data)),levels="All patients"); gr.lev <- levels(group)}
  else {
    blvars <- dat[,-1]
    if (is.null(ncol(blvars))) {
      dim(blvars) <- c(length(blvars), 1)
      colnames(blvars) <- as.character(formula[[3]])
    }
    group <- droplevels(factor(dat[,1]))
    if (is.logical(dat[,1])) gr.lev <- as.character(unique(dat[,1])) else gr.lev <- levels(dat[,1])
    if (pooledGroup) {
      mylabels <- getlabel(blvars)
      blvars <- rbind(blvars,blvars)
      for (i in 1:ncol(blvars)) attr(blvars[,i], "label") <- mylabels[i]
      group <- c(as.character(group),rep("All patients",nrow(data)))
      group <- factor(group,levels=c("All patients",gr.lev))
      gr.lev <- levels(group)
    }
  }

  gr.lev <- levels(group)
  header1 <- c("", c(rbind(rep("", length(gr.lev)), paste(gr.lev, " (N=", table(group), ")", sep = ""))))
  header2 <- c("Characteristic", rep(c("n", "Summary statistic"), length(gr.lev)))
  if (test) {
    header1 <- c(header1, "Comparison")
    header2 <- c(header2, "(p-value)")
  }
  result <-  rbind(header1, header2)
  if (length(continuous) == 1) continuous <- rep(continuous, ncol(blvars))
  if (length(digits) == 1) digits <- rep(digits, ncol(blvars))

  for (i in 1:ncol(blvars)){
    result.i <- mySummary.onevar(varname = ifelse(getlabel(blvars[, i]) != "", getlabel(blvars[, i]), getlabel(blvars)[i]),
                                 blvars[, i], group, contSummary = contSummary, test = test,
                                 continuous = continuous[i], digits = digits[i], pcutoff = pcutoff)
    result <- rbind(result, result.i)
  }
  rownames(result) <- rep("", nrow(result))

  if (kable){
    if (test){
      align <- c("l", "r", "c", "r", "c", "c")
    } else {
      align <- c("l", "r", "c", "r", "c")
    }
    result[2,] <- paste("_", result[2,], "_", sep = "_")
    tab <- knitr::kable(result[-1,],
                 row.names = FALSE,
                 col.names = header1,
                 align = align,
                 caption = caption)
    footnote <- paste("_Summary statistic is absolute count (%) for categorical variables and",
                      switch(contSummary,
                             med.IQR = "median (IQR)",
                             med.90  = "median (90% range)",
                             med.range = "median (range)",
                             mean.sd = "mean (sd)"),
                      "for continuous data._")
    if (test) footnote <- c(footnote,
                            "",
                            paste0("_p-values based on Kruskal-Wallis/Mann-Whitney U-test (continuous data) and ",
                                   ifelse(chisq.test == FALSE, "Fisher's exact test", "Chi-squared test"), " (categorical data)._"))
    ## output
    structure(c(tab, "", footnote), format = "markdown", class = "knitr_kable")
  } else {
    result
  }
}

## Special function
#' @export
mySummary.simple <- function(formula, data, pooledGroup = TRUE, continuous = NA) {
  
  mySummary.onevar2 <- function(variable, group = NULL, continuous = NA, pooledGroup){
    if (is.na(continuous)) continuous <- ifelse(is.factor(variable) | length(unique(na.omit(variable))) <= 5, FALSE, TRUE)
    
    tmpfunc <- function(x, range) {
      return((median(x, na.rm = TRUE) - range[1])/(range[2] - range[1]))
    }
    
    mycont.summary <- function(variable,group) {
      tmp <- quantile(variable, probs = c(0.1, 0.9), na.rm = TRUE)
      if (is.null(group)) {
        ngroup <- 1
        n <- length(na.omit(variable))
        xscale <- tmpfunc(variable, range = tmp)
        xmean <- mean(variable, na.rm = TRUE)
        xvar  <- var(variable, na.rm = TRUE)
        
      } else {
        ngroup <- length(levels(group))
        n <- c(by(variable, group, function(x) length(na.omit(x))))
        xscale <- by(variable, group, tmpfunc, range = tmp)
        xmean <- by(variable, group, mean, na.rm = TRUE)
        xvar <- by(variable, group, var, na.rm = TRUE)
        
      }
      
      result <- matrix(NA, ncol = ngroup * 4, nrow = 1)
      result[1, seq(1, ncol(result), by = 4)] <- n
      result[1, seq(2, ncol(result), by = 4)] <- unlist(xscale)
      result[1, seq(3, ncol(result), by = 4)] <- unlist(xmean)
      result[1, seq(4, ncol(result), by = 4)] <- unlist(xvar)
      return(result)
    }
    
    mycat.summary <- function(variable, group) {
      
      if (is.null(group)) {
        ngroup <- 1
        ta <- ta2 <- table(factor(variable, levels = c(FALSE, TRUE)))
        ta2[ta2 == 0] <- 0.5
        ta.prop <- ta/sum(ta)
        ta2.prop <- ta2/sum(ta2)
        dim(ta) <- c(ngroup, length(ta))
        colnames(ta) <- names(table(variable))
      } else {
        ngroup <- length(levels(group))
        ta <- ta2 <- table(group, factor(variable, levels = c(FALSE, TRUE)))
        ta2[ta2 == 0] <- 0.5
        ta.prop <- unclass(ta/apply(ta, 1, sum))
        ta2.prop <- unclass(ta2/apply(ta2, 1, sum))
      }
      
      n <- apply(ta, 1, sum)
      xmean <- ta2.prop[, "TRUE"]
      xvar <- xmean * (1 - xmean)/n
      
      result <- matrix(NA, ncol = ngroup * 4, nrow = 1)
      result[1, seq(1, ncol(result), by = 4)] <- n
      result[1, seq(2, ncol(result), by = 4)] <- ta.prop[, "TRUE"]
      result[1, seq(3, ncol(result), by = 4)] <- xmean
      result[1, seq(4, ncol(result), by = 4)] <- xvar
      return(result)
    }
    
    if (continuous) r <- mycont.summary(variable, group)
    else r <- mycat.summary(variable, group)
    
    ## calculate heterogeneity index
    #browser()
    if (pooledGroup) {r2 <- t(r[, -c(1:4)])}
    M <- r2[, seq(3, ncol(r2), by = 4)]
    W <- 1/r2[, seq(4, ncol(r2), by = 4)]
    Q <- sum(W * (M - sum(W * M)/sum(W))^2)
    df <- ncol(r2)/4 - 1
    I2 <- pmax(100 * (Q - df)/Q, 0)
    H <- 1/sqrt(1 - I2/100)
    
    return(cbind(r, Q, df, I2, H))
  }
  
  dat <- model.frame(formula, data = data, na.action = NULL)
  if (length(formula) == 2) {blvars <- dat; group <- factor(rep("All patients",nrow(data)),levels = "All patients"); gr.lev <- levels(group)}
  else {
    blvars <- dat[,-1]
    if (is.null(ncol(blvars))) {
      dim(blvars) <- c(length(blvars), 1)
      colnames(blvars) <- as.character(formula[[3]])
    }
    group <- droplevels(factor(dat[,1]))
    if (is.logical(dat[,1])) gr.lev <- as.character(unique(dat[,1])) else gr.lev <- levels(dat[,1])
    if (pooledGroup) {
      mylabels <- getlabel(blvars)
      blvars <- rbind(blvars,blvars)
      for (i in 1:ncol(blvars)) attr(blvars[,i], "label") <- mylabels[i]
      group <- c(as.character(group),rep("All patients",nrow(data)))
      group <- factor(group,levels = c("All patients",gr.lev))
      gr.lev <- levels(group)
    }
  }
  
  gr.lev <- levels(group)
  result <-  NULL
  if (length(continuous) == 1) continuous <- rep(continuous, ncol(blvars))
  
  for (i in 1:ncol(blvars)) {
    result.i <- mySummary.onevar2(blvars[, i], group, continuous = continuous[i], pooledGroup = pooledGroup)
    result <- rbind(result, result.i)
  }
  colnames(result) <- c(rbind(rep("n", length(gr.lev)), 
                              paste(gr.lev, " (N=", table(group), ")", sep = ""),
                              paste(gr.lev, "mean", sep = "."),
                              paste(gr.lev, "var", sep = ".")),
                        "Q", "df", "I2", "H")
  rownames(result) <- getlabel(blvars)
  return(result)
}

#' @export
mySummary.simple2 <- function(formula, data, pooledGroup = TRUE, continuous = NA) {
  require(metafor)
  mySummary.onevar2 <- function(variable, group = NULL, continuous = NA, pooledGroup){
    if (is.na(continuous)) continuous <- ifelse(is.factor(variable) | length(unique(na.omit(variable))) <= 5, FALSE, TRUE)
    
    tmpfunc <- function(x, range) {
      return((median(x, na.rm = TRUE) - range[1])/(range[2] - range[1]))
    }
    
    mycont.summary <- function(variable, group, pooledGroup) {
      tmp <- quantile(variable, probs = c(0.1, 0.9), na.rm = TRUE)
      if (is.null(group)) {
        ngroup <- 1
        n <- length(na.omit(variable))
        xscale <- tmpfunc(variable, range = tmp)
        xmean <- mean(variable, na.rm = TRUE)
        xvar  <- var(variable, na.rm = TRUE)
        
      } else {
        ngroup <- length(levels(group))
        n <- c(by(variable, group, function(x) length(na.omit(x))))
        xscale <- by(variable, group, tmpfunc, range = tmp)
        xmean <- by(variable, group, mean, na.rm = TRUE)
        xvar <- by(variable, group, var, na.rm = TRUE)
        
      }
      
      result <- matrix(NA, ncol = ngroup * 4, nrow = 1)
      result[1, seq(1, ncol(result), by = 4)] <- n
      result[1, seq(2, ncol(result), by = 4)] <- unlist(xscale)
      result[1, seq(3, ncol(result), by = 4)] <- unlist(xmean)
      result[1, seq(4, ncol(result), by = 4)] <- unlist(xvar)
      
      if (pooledGroup) {
        rmai <- metafor::rma(measure = "MN", mi = as.numeric(unlist(xmean)[-1]), sdi = as.numeric(sqrt(unlist(xvar))[-1]), ni = n[-1], method = "DL")
      } else {
        rmai <- metafor::rma(measure = "MN", mi = as.numeric(unlist(xmean)), sdi = as.numeric(sqrt(unlist(xvar))), ni = n, method = "DL")
      }
      return(cbind(result, I2 = rmai$I2))
    }
    
    mycat.summary <- function(variable, group, pooledGroup) {
      
      if (is.null(group)) {
        ngroup <- 1
        ta <- ta2 <- table(factor(variable, levels = c(FALSE, TRUE)))
        ta2[ta2 == 0] <- 0.5
        ta.prop <- ta/sum(ta)
        ta2.prop <- ta2/sum(ta2)
        dim(ta) <- c(ngroup, length(ta))
        colnames(ta) <- names(table(variable))
      } else {
        ngroup <- length(levels(group))
        ta <- ta2 <- table(group, factor(variable, levels = c(FALSE, TRUE)))
        ta2[ta2 == 0] <- 0.5
        ta.prop <- unclass(ta/apply(ta, 1, sum))
        ta2.prop <- unclass(ta2/apply(ta2, 1, sum))
      }
      
      n <- apply(ta, 1, sum)
      xmean <- ta2.prop[, "TRUE"]
      xvar <- xmean * (1 - xmean)/n
      
      result <- matrix(NA, ncol = ngroup * 4, nrow = 1)
      result[1, seq(1, ncol(result), by = 4)] <- n
      result[1, seq(2, ncol(result), by = 4)] <- ta.prop[, "TRUE"]
      result[1, seq(3, ncol(result), by = 4)] <- xmean
      result[1, seq(4, ncol(result), by = 4)] <- xvar
      
      if (pooledGroup) {
        rmai <- metafor::rma(measure = "PR", xi = ta[, "TRUE"][-1], ni = n[-1], method = "DL")
      } else {
        rmai <- metafor::rma(measure = "PR", xi = ta[, "TRUE"], ni = n, method = "DL")
      }
      return(cbind(result, I2 = rmai$I2))
    }
    
    if (continuous) return(mycont.summary(variable, group, pooledGroup))
    else return(mycat.summary(variable, group, pooledGroup))
  }
  
  dat <- model.frame(formula, data = data, na.action = NULL)
  if (length(formula) == 2) {blvars <- dat; group <- factor(rep("All patients",nrow(data)),levels = "All patients"); gr.lev <- levels(group)}
  else {
    blvars <- dat[,-1]
    if (is.null(ncol(blvars))) {
      dim(blvars) <- c(length(blvars), 1)
      colnames(blvars) <- as.character(formula[[3]])
    }
    group <- droplevels(factor(dat[,1]))
    if (is.logical(dat[,1])) gr.lev <- as.character(unique(dat[,1])) else gr.lev <- levels(dat[,1])
    if (pooledGroup) {
      mylabels <- getlabel(blvars)
      blvars <- rbind(blvars,blvars)
      for (i in 1:ncol(blvars)) attr(blvars[,i], "label") <- mylabels[i]
      group <- c(as.character(group),rep("All patients",nrow(data)))
      group <- factor(group,levels = c("All patients",gr.lev))
      gr.lev <- levels(group)
    }
  }
  
  gr.lev <- levels(group)
  result <-  NULL
  if (length(continuous) == 1) continuous <- rep(continuous, ncol(blvars))
  
  for (i in 1:ncol(blvars)) {
    result.i <- mySummary.onevar2(blvars[, i], group, continuous = continuous[i], pooledGroup = pooledGroup)
    result <- rbind(result, result.i)
  }
  colnames(result) <- c(rbind(rep("n", length(gr.lev)), 
                              paste(gr.lev, " (N=", table(group), ")", sep = ""),
                              paste(gr.lev, "mean", sep = "."),
                              paste(gr.lev, "var", sep = ".")),
                        "I2")
  rownames(result) <- getlabel(blvars)
  return(result)
}