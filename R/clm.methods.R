print.clm <- 
  function(x, digits = max(3, getOption("digits") - 3), ...)
{
  cat("formula:", deparse(x$call$formula), fill=TRUE)
  if(!is.null(x$call$scale))
    cat("scale:  ", deparse(x$call$scale), fill=TRUE)
  if(!is.null(x$call$nominal))
    cat("nominal:", deparse(x$call$nominal), fill=TRUE)
  if(!is.null(data.name <- x$call$data))
    cat("data:   ", deparse(data.name), fill=TRUE)
  if(!is.null(x$call$subset))
    cat("subset: ", deparse(x$call$subset), fill=TRUE)
  cat("\n")

  print(x$info, row.names=FALSE, right=FALSE)

  if(length(x$beta)) {
    if(sum(x$aliased$beta) > 0) {
      cat("\nCoefficients: (", sum(x$aliased$beta),
          " not defined because of singularities)\n", sep = "")
    }
    else cat("\nCoefficients:\n")
    print.default(format(x$beta, digits = digits), 
                  quote = FALSE)
  }
    if(length(x$zeta)) {
    if(sum(x$aliased$zeta) > 0) 
      cat("\nlog-scale coefficients: (", sum(x$aliased$zeta),
          " not defined because of singularities)\n", sep = "") 
    else cat("\nlog-scale coefficients:\n")
    print.default(format(x$zeta, digits = digits), 
                  quote = FALSE)
  }
  if(length(x$alpha) > 0) {
    if(sum(x$aliased$alpha) > 0) 
      cat("\nThreshold coefficients: (", sum(x$aliased$alpha),
          " not defined because of singularities)\n", sep = "") 
    else cat("\nThreshold coefficients:\n")
    print.default(format(x$alpha, digits = digits), 
                  quote = FALSE)
  }

  if(nzchar(mess <- naprint(x$na.action))) cat("(", mess, ")\n", sep="")
  return(invisible(x))
}

vcov.clm <- function(object, ...)
{
  if(is.null(object$Hessian)) 
    stop("Model needs to be fitted with Hess = TRUE")
  dn <- dimnames(object$Hessian)
  H <- object$Hessian
  ## To handle NaNs in the Hessian resulting from parameter
  ## unidentifiability:  
  if(any(His.na <- !is.finite(H))) {
    H[His.na] <- 0
    VCOV <- MASS::ginv(H)
    VCOV[His.na] <- NaN
  }
  else
    VCOV <- solve(H) ## MASS::ginv(H)
  return(structure(VCOV, dimnames = dn))
}

summary.clm <- function(object, correlation = FALSE, ...)
{
  if(is.null(object$Hessian))
    stop("Model needs to be fitted with Hess = TRUE")
  coefs <- matrix(NA, length(object$coefficients), 4,
                  dimnames = list(names(object$coefficients),
                    c("Estimate", "Std. Error", "z value", "Pr(>|z|)")))
  coefs[, 1] <- object$coefficients
  cov <- try(vcov(object), silent = TRUE)
  if(class(cov) == "try-error") {
    warning("Variance-covariance matrix of the parameters is not defined")
    coefs[, 2:4] <- NaN
    if(correlation) warning("Correlation matrix is unavailable")
    object$condHess <- NaN
  }
  else {
    alias <- unlist(object$aliased)
    coefs[!alias, 2] <- sd <- sqrt(diag(cov))
    ## Cond is Inf if Hessian contains NaNs:
    object$condHess <-
      if(any(is.na(object$Hessian))) Inf
      else with(eigen(object$Hessian, symmetric=TRUE, only.values = TRUE),
                abs(max(values) / min(values)))
    coefs[!alias, 3] <- coefs[!alias, 1]/coefs[!alias, 2]
    coefs[!alias, 4] <- 2 * pnorm(abs(coefs[!alias, 3]),
                                  lower.tail=FALSE) 
    if(correlation)
      object$correlation <-
        (cov / sd) / rep(sd, rep(object$edf, object$edf))
    object$cov <- cov
  }
  object$info$cond.H <- formatC(object$condHess, digits=1, format="e")
  object$coefficients <- coefs
  class(object) <- "summary.clm"
  return(object)
}

print.summary.clm <- 
  function(x, digits = max(3, getOption("digits") - 3),
           signif.stars = getOption("show.signif.stars"), ...)
{
  cat("formula:", deparse(x$call$formula), fill=TRUE)
  if(!is.null(x$call$scale))
    cat("scale:  ", deparse(x$call$scale), fill=TRUE)
  if(!is.null(x$call$nominal))
    cat("nominal:", deparse(x$call$nominal), fill=TRUE)
  if(!is.null(data.name <- x$call$data))
    cat("data:   ", deparse(data.name), fill=TRUE)
  if(!is.null(x$call$subset))
    cat("subset: ", deparse(x$call$subset), fill=TRUE)
  cat("\n")

  print(x$info, row.names=FALSE, right=FALSE)
  
  nalpha <- length(x$alpha)
  nbeta <- length(x$beta)
  nzeta <- length(x$zeta)
  if(nbeta > 0) {
    if(sum(x$aliased$beta) > 0) 
      cat("\nCoefficients: (", sum(x$aliased$beta),
          " not defined because of singularities)\n", sep = "") 
    else cat("\nCoefficients:\n")
    printCoefmat(x$coefficients[nalpha + 1:nbeta, , drop=FALSE],
                 digits=digits, signif.stars=signif.stars,
                 has.Pvalue=TRUE, ...) 
  } ## else  cat("\nNo Coefficients\n")
  if(nzeta > 0) {
    if(sum(x$aliased$zeta) > 0) 
      cat("\nlog-scale coefficients: (", sum(x$aliased$zeta),
          " not defined because of singularities)\n", sep = "") 
    else cat("\nlog-scale coefficients:\n")
    printCoefmat(x$coefficients[nalpha + nbeta + 1:nzeta, , drop=FALSE],
                 digits=digits, signif.stars=signif.stars,
                 has.Pvalue=TRUE, ...) 
  }
  if(nalpha > 0) { ## always true
    if(sum(x$aliased$alpha) > 0) 
      cat("\nThreshold coefficients: (", sum(x$aliased$alpha),
          " not defined because of singularities)\n", sep = "") 
    else cat("\nThreshold coefficients:\n")
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
  
logLik.clm <- function(object, ...)
  structure(object$logLik, df = object$edf, class = "logLik")

extractAIC.clm <- function(fit, scale = 0, k = 2, ...) {
  edf <- fit$edf
  c(edf, -2*fit$logLik + k * edf)
}

### NOTE: AIC.clm implicitly defined via logLik.clm

anova.clm <- function(object, ...)
### requires that clm objects have components:
###  edf: no. parameters used
###  call$formula
###  link (character)
###  threshold (character)
###  logLik
###  
{
  mc <- match.call()
  dots <- list(...)
  if (length(dots) == 0)
    stop('anova is not implemented for a single "clm" object')
  mlist <- list(object, ...)
  if(!all(sapply(mlist, function(model)
                 inherits(model, c("clm", "clmm")))))
    stop("only 'clm' and 'clmm' objects are allowed")
  nfitted <- sapply(mlist, function(x) length(x$fitted.values))
  if(any(nfitted != nfitted[1L])) 
    stop("models were not all fitted to the same dataset")
### FIXME: consider comparing y returned by the models for a better
### check? 
  no.par <- sapply(mlist, function(x) x$edf)
  ## order list with increasing no. par:
  ord <- order(no.par, decreasing=FALSE)
  mlist <- mlist[ord]
  no.par <- no.par[ord]
  no.tests <- length(mlist)
  ## extract formulas, links, thresholds, scale formulas, nominal
  ## formulas:
  forms <- sapply(mlist, function(x) deparse(x$call$formula))
  links <- sapply(mlist, function(x) x$link)
  thres <- sapply(mlist, function(x) x$threshold)
  nominal <- sapply(mlist, function(x) deparse(x$call$nominal))
  scale <- sapply(mlist, function(x) deparse(x$call$scale))
  models <- data.frame(forms)
  models.names <- 'formula:'
  if(any(!nominal %in% c("~1", "NULL"))) {
    nominal[nominal == "NULL"] <- "~1"
    models$nominal <- nominal
    models.names <- c(models.names, "nominal:")
  }
  if(any(!scale %in% c("~1", "NULL"))) {
    scale[scale == "NULL"] <- "~1"
    models$scale <- scale
    models.names <- c(models.names, "scale:")
  }
  models.names <- c(models.names, "link:", "threshold:")
  models <- cbind(models, data.frame(links, thres))
  ## extract AIC, logLik, statistics, df, p-values:
  AIC <- sapply(mlist, function(x) -2*x$logLik + 2*x$edf)
  logLiks <- sapply(mlist, function(x) x$logLik)
  statistic <- c(NA, 2*diff(sapply(mlist, function(x) x$logLik)))
  df <- c(NA, diff(no.par))
  pval <- c(NA, pchisq(statistic[-1], df[-1], lower.tail=FALSE))
  pval[!is.na(df) & df==0] <- NA
  ## collect results in data.frames:
  tab <- data.frame(no.par, AIC, logLiks, statistic, df, pval) 
  tab.names <- c("no.par", "AIC", "logLik", "LR.stat", "df",
                 "Pr(>Chisq)")
  mnames <- sapply(as.list(mc), deparse)[-1]
  colnames(tab) <- tab.names
  rownames(tab) <- rownames(models) <- mnames[ord]
  colnames(models) <- models.names
  attr(tab, "models") <- models
  attr(tab, "heading") <-
    "Likelihood ratio tests of cumulative link models:\n"
  class(tab) <- c("anova.clm", "data.frame")
  tab
}

print.anova.clm <-
  function(x, digits=max(getOption("digits") - 2, 3),
           signif.stars=getOption("show.signif.stars"), ...) 
{
  if (!is.null(heading <- attr(x, "heading"))) 
    cat(heading, "\n")
  models <- attr(x, "models")
  print(models, right=FALSE)
  cat("\n")
  printCoefmat(x, digits=digits, signif.stars=signif.stars,
               tst.ind=4, cs.ind=NULL, # zap.ind=2, #c(1,5),
               P.values=TRUE, has.Pvalue=TRUE, na.print="", ...)
  return(invisible(x))
}

model.matrix.clm <- function(object, type = c("design", "B"), ...)
### returns a list of model matrices for the X, NOM and S formulas or
### the B matrices (including S if present) actually used for the
### fitting
### Aliased columns are retained in the former but dropped in the
### latter. 
{
  type <- match.arg(type)
  if(type == "design") {
    mf <- update(object, method="model.frame")
    keep <- c("X", "NOM", "S")
    select <- match(keep, names(mf), nomatch=0)
    return(mf[select])
  } else {
    env <- update(object, doFit=FALSE)
    ans <- list(B1 = env$B1, B2 = env$B2)
    ans$S <- env$S ## may not exist
    return(ans)
  }
}

model.frame.clm <- function(formula, ...) {
### returns a model frame with *all* variables used for fitting. 
  if(is.null(mod <- formula$model))
    update(formula, method="model.frame")$mf
  else
    mod
}

coef.clm <- function(object, na.rm = FALSE, ...) {
  if(na.rm) {
    coefs <- object$coefficients
    coefs[!is.na(coefs)]
  }
  else 
    object$coefficients
}

coef.summary.clm <- function(object, na.rm = FALSE, ...) {
  if(na.rm) {
    coefs <- object$coefficients
    coefs[!is.na(coefs[,1]), , drop=FALSE]
  }
  else
    object$coefficients
}
  
nobs.clm <- function(object, ...) object$nobs

