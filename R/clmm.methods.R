print.clmm <- 
  function(x, digits = max(3, getOption("digits") - 3), ...)
{
  if(x$nAGQ >= 2)
    cat(paste("Cumulative Link Mixed Model fitted with the adaptive",
              "Gauss-Hermite \nquadrature approximation with",
              x$nAGQ ,"quadrature points"), "\n\n")
  else if(x$nAGQ <= -1)
    cat(paste("Cumulative Link Mixed Model fitted with the",
              "Gauss-Hermite \nquadrature approximation with",
              abs(x$nAGQ) ,"quadrature points"), "\n\n")
  else
    cat("Cumulative Link Mixed Model fitted with the Laplace approximation\n",
      fill=TRUE)
  cat("formula:", deparse(x$call$formula), fill=TRUE)
  if(!is.null(data.name <- x$call$data))
    cat("data:   ", deparse(data.name), fill=TRUE)
  if(!is.null(x$call$subset))
    cat("subset: ", deparse(x$call$subset), fill=TRUE)
  cat("\n")

  print(x$info, row.names=FALSE, right=FALSE)

  cat("\nRandom effects:\n")
  print(x$varMat, digits=digits, ...)
  nlev.char <- paste(names(x$nlev), " ", x$nlev, sep="", collapse=",  ")
  cat("Number of groups: ", nlev.char, "\n")

  if(length(x$beta)) {
    cat("\nCoefficients:\n")
    print(x$beta, digits=digits, ...)
  } else {
    cat("\nNo Coefficients\n")
  }
  if(length(x$alpha) > 0) {
    cat("\nThresholds:\n")
    print(x$alpha, digits=digits, ...)
  }

  if(nzchar(mess <- naprint(x$na.action))) cat("(", mess, ")\n", sep="")
  return(invisible(x))
}

vcov.clmm <- function(object, ...) vcov.clm(object, ...)

summary.clmm <- function(object, correlation = FALSE, ...)
{
  if(is.null(object$Hessian))
    stop("Model needs to be fitted with Hess = TRUE")
  
  npar <- length(object$alpha) + length(object$beta)
  coef <- matrix(0, npar, 4,
                 dimnames = list(names(object$coefficients[1:npar]),
                   c("Estimate", "Std. Error", "z value", "Pr(>|z|)")))
  coef[, 1] <- object$coefficients[1:npar]
  vc <- try(vcov(object), silent = TRUE)
  if(class(vc) == "try-error") {
    warning("Variance-covariance matrix of the parameters is not defined")
    coef[, 2:4] <- NaN
    if(correlation) warning("Correlation matrix is unavailable")
    object$condHess <- NaN
  }
  else {
    coef[, 2] <- sd <- sqrt(diag(vc)[1:npar])
    ## Cond is Inf if Hessian contains NaNs:
    object$condHess <-
      if(any(is.na(object$Hessian))) Inf
      else with(eigen(object$Hessian, only.values = TRUE),
                abs(max(values) / min(values)))
    coef[, 3] <- coef[, 1]/coef[, 2]
    coef[, 4] <- 2 * pnorm(abs(coef[, 3]), lower.tail=FALSE)
    if(correlation) ## {
      ## sd <- sqrt(diag(vc))
      object$correlation <- cov2cor(vc)
    ## (vc / sd) / rep(sd, rep(object$edf, object$edf))
  }
  object$info$cond.H <- formatC(object$condHess, digits=1, format="e")
  object$coefficients <- coef
  class(object) <- "summary.clmm"
  return(object)
}

print.summary.clmm <-
  function(x, digits = max(3, getOption("digits") - 3),
           signif.stars = getOption("show.signif.stars"), ...)
{
  if(x$nAGQ >= 2)
    cat(paste("Cumulative Link Mixed Model fitted with the adaptive",
              "Gauss-Hermite \nquadrature approximation with",
              x$nAGQ ,"quadrature points"), "\n\n")
  else if(x$nAGQ <= -1)
    cat(paste("Cumulative Link Mixed Model fitted with the",
              "Gauss-Hermite \nquadrature approximation with",
              abs(x$nAGQ) ,"quadrature points"), "\n\n")
  else
    cat("Cumulative Link Mixed Model fitted with the Laplace approximation\n",
      fill=TRUE)
  cat("formula:", deparse(x$call$formula), fill=TRUE)
  if(!is.null(data.name <- x$call$data))
    cat("data:   ", deparse(data.name), fill=TRUE)
  if(!is.null(x$call$subset))
    cat("subset: ", deparse(x$call$subset), fill=TRUE)
  cat("\n")

  print(x$info, row.names=FALSE, right=FALSE)

  cat("\nRandom effects:\n")
  print(x$varMat, digits=digits, ...)
  nlev.char <- paste(names(x$nlev), " ", x$nlev, sep="", collapse=",  ")
  cat("Number of groups: ", nlev.char, "\n")
  
  nbeta <- length(x$beta)
  nalpha <- length(x$alpha)
  if(nbeta > 0) {
    cat("\nCoefficients:\n")
    printCoefmat(x$coefficients[nalpha + 1:nbeta, , drop=FALSE],
                 digits=digits, signif.stars=signif.stars,
                 has.Pvalue=TRUE, ...) 
  } else {
    cat("\nNo Coefficients\n")
  }
  if(nalpha > 0) { ## always true
    cat("\nThreshold coefficients:\n")
    printCoefmat(x$coefficients[seq_len(nalpha), -4, drop=FALSE],
                 digits=digits, has.Pvalue=FALSE, signif.stars=FALSE,
                 ...) 
  }
  
  if(nzchar(mess <- naprint(x$na.action))) cat("(", mess, ")\n", sep="")
  if(!is.null(correl <- x$correlation)) {
    cat("\nCorrelation of Coefficients:\n")
    ll <- lower.tri(correl)
    correl[ll] <- format(round(correl[ll], digits))
    correl[!ll] <- ""
    print(correl[-1, -ncol(correl)], quote = FALSE, ...)
  }
  return(invisible(x))
}

## anova.clmm <- function(object, ...)
##   anova.clm(object, ...)

anova.clmm <- function(object, ...) {
### This essentially calls anova.clm(object, ...), but the names of
### the models were not displayed correctly in the printed output
### unless the following dodge is enforced.
  mc <- match.call()
  arg.list <- as.list(mc)
  arg.list[[1]] <- NULL
  return(do.call(anova.clm, arg.list))
}

logLik.clmm <- function(object, ...)
  structure(object$logLik, df = object$edf, class = "logLik")

extractAIC.clmm <- function(fit, scale = 0, k = 2, ...) {
  edf <- fit$edf
  c(edf, -2*fit$logLik + k * edf)
}
