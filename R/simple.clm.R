simple_clm <-
  function(formula, data, weights, start, subset, offset,
           doFit = TRUE, na.action, contrasts, model = TRUE,
           control = list(),
           link = c("logit", "probit", "cloglog", "loglog"), 
           threshold = c("flexible", "symmetric", "equidistant"), ...)
{
  ## Initial argument matching and testing:
  mc <- match.call(expand.dots = FALSE)
  link <- match.arg(link)
  threshold <- match.arg(threshold)
  ## check for presence of formula:
  if(missing(formula)) stop("Model needs a formula")
  if(missing(contrasts)) contrasts <- NULL
  ## set control parameters:
  control <- do.call(clm.control, c(control, list(...)))

  ## Compute: y, X, wts, off, mf:
  if (missing(data)) 
    data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "na.action", 
               "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  ## Return model.frame?
  if(control$method == "model.frame") return(mf)
  y <- model.response(mf, "any") ## any storage mode
  if(!is.factor(y)) stop("response needs to be a factor", call.=FALSE)
  ## design matrix:
  mt <- attr(mf, "terms")
  X <- if (!is.empty.model(mt)) 
    model.matrix(mt, mf, contrasts)
  else cbind("(Intercept)" = rep(1, NROW(y)))
  ## Test for intercept in X:
  Xint <- match("(Intercept)", colnames(X), nomatch = 0)
  if(Xint <= 0) {
    X <- cbind("(Intercept)" = rep(1, NROW(y)), X)
    warning("an intercept is needed and assumed in 'formula'",
            call.=FALSE)
  } ## intercept in X is guaranteed.
  frames <- list(y=y, X=X)
  wts <- ordinal:::getWeights(mf)
  off <- ordinal:::getOffset(mf)

  ## Compute the transpose of the Jacobian for the threshold function,
  ## tJac and the names of the threshold parameters, alpha.names:
  frames$ths <- ordinal:::makeThresholds(y, threshold)
  ## test for column rank deficiency in design matrices:
  frames <- ordinal:::drop.cols(frames, silent=TRUE)

  ## Set envir rho with variables: B1, B2, o1, o2, wts, fitted:
  rho <- ordinal:::eclm.newRho(parent.frame(), y=frames$y, X=frames$X,
                     NOM=NULL, S=NULL, 
                     weights=wts, offset=off, S.offset=NULL,
                     tJac=frames$ths$tJac)

  ## Set appropriate logLik and deriv functions in rho:
  rho$clm.nll <- ordinal:::clm.nll
  rho$clm.grad <- ordinal:::clm.grad
  rho$clm.hess <- ordinal:::clm.hess

  ## Set starting values for the parameters:
  start <- ordinal:::set.start(rho, start=start, get.start=missing(start),
                     threshold=threshold, link=link, frames=frames)
  rho$par <- as.vector(start) ## remove attributes
  
  ## Set pfun, dfun and gfun in rho:
  ordinal:::setLinks(rho, link)

  ## Possibly return the environment rho without fitting:
  if(!doFit) return(rho)
  
  ## Fit the clm:
  if(control$method == "Newton") 
    fit <- ordinal:::clm.fit.env(rho, control)
  else
    fit <- ordinal:::clm.fit.optim(rho, control$method, control$ctrl)
### NOTE: we could add arg non.conv = c("error", "warn", "message") to
### allow non-converged fits to be returned.

  ## Modify and return results:
  res <- ordinal:::eclm.finalize(fit, weights=wts,
                       coef.names=frames$coef.names,
                       aliased=frames$aliased) 
  res$link <- link
  res$start <- start
  if(control$method == "Newton" &&
     !is.null(start.iter <- attr(start, "start.iter")))
    res$niter <- res$niter + start.iter
  res$threshold <- threshold
  res$call <- match.call()
  res$contrasts <- attr(frames$X, "contrasts")
  res$na.action <- attr(mf, "na.action")
  res$terms <- mt
  res$xlevels <- .getXlevels(mt, mf)
  res$tJac <- frames$ths$tJac
  res$y.levels <- levels(frames$y)
  res$info <- with(res, {
    data.frame("link" = link,
               "threshold" = threshold,
               "nobs" = nobs, 
               "logLik" = formatC(logLik, digits=2, format="f"),
               "AIC" = formatC(-2*logLik + 2*edf, digits=2,
                 format="f"),
               "niter" = paste(niter[1], "(", niter[2], ")", sep=""),
### NOTE: iterations to get starting values for scale models are not
### included here.
               "max.grad" = formatC(maxGradient, digits=2,
                 format="e") 
               ## BIC is not part of output since it is not clear what
               ## the no. observations are. 
               )
  })
  class(res) <- c("sclm", "clm")
  ## add model.frame to results list?
  if(model) res$model <- mf
  
  return(res)
}
