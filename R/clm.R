## clm <-
##   function(formula, data, weights, start, subset, doFit = TRUE,
##            na.action, contrasts, model = TRUE, control,
##            link = c("logit", "probit", "cloglog", "loglog", "cauchit"), 
##            threshold = c("flexible", "symmetric", "equidistant"), ...)
## {
##   ## Initial argument matching and testing:
##   mc <- match.call(expand.dots = FALSE)
##   link <- match.arg(link)
##   threshold <- match.arg(threshold)
##   ## check for presence of formula:
##   if(missing(formula)) stop("Model needs a formula")
##   if(missing(contrasts)) contrasts <- NULL
##   
##   ## set control parameters:
##   if(missing(control)) control <- clm.control(...)
##   if(!setequal(names(control), c("method", "trace", "maxIter",
##                                  "gradTol", "maxLineIter")))
##     stop("specify 'control' via clm.control()")
## 
##   ## make data a data.frame:
##   ##    if(missing(data)) mc$data <- environment(formula)
##   ##    if(is.matrix(eval.parent(mc$data)))
##   ##      mc$data <- as.data.frame(mc$data)
##   
##   ## Compute: y, X, wts, off, mf:
##   frames <- clm.model.frame(mc, contrasts)
##   if(control$method == "model.frame") return(frames)
## 
##   ## Test rank deficiency and possibly drop some parameters:
##   ## X is guarantied to have an intercept at this point.
##   frames$X <- drop.coef(frames$X)
##   
##   ## Compute the transpose of the Jacobian for the threshold function,
##   ## tJac and the names of the threshold parameters, alpha.names:
##   ths <- makeThresholds(frames$y, threshold)
## 
##   ## set envir rho with variables: B1, B2, o1, o2, wts, fitted:
##   rho <- with(frames,
##               clm.newRho(parent.frame(), y=y, X=X, weights=wts,
##                          offset=off, tJac=ths$tJac) )
## 
##   ## set starting values for the parameters:
##   if(missing(start)) 
##     start <- clm.start(frames$y, frames$X, has.intercept = TRUE,
##                        threshold = threshold)
##   stopifnot(is.numeric(start) && 
##             length(start) == (ths$nalpha + ncol(frames$X) - 1) )
##   rho$par <- start
##   ## start cauchit models at the probit estimates:
##   if(link == "cauchit" && is.null(match.call()$start)) {
##     setLinks(rho, link="probit")
##     fit <- try(clm.fit.env(rho), silent=TRUE) ## standard control values
##     if(class(fit) == "try-error") 
##       stop("Failed to find suitable starting values: please supply some")
##     rho$par <- fit$par
##   }
##   
##   ## Set pfun, dfun and gfun in rho:
##   setLinks(rho, link)
## 
##   ## possibly return the environment rho without fitting:
##   if(!doFit) return(rho)
## 
##   ## fit the clm:
##   control$method <- NULL
##   fit <- clm.fit.env(rho, control)
## ### FIXME: add arg non.conv = c("error", "warn", "message") to allow
## ### non-converged fits to be returned?
## 
##   ## Modify and return results:
##   res <- clm.finalize(fit, weights = frames$wts,
##                       alpha.names = ths$alpha.names,
##                       beta.names = colnames(frames$X)[-1])
##   res$link <- link
##   res$start <- start
##   res$threshold <- threshold
##   res$call <- match.call()
##   res$contrasts <- contrasts
##   res$na.action <- attr(frames$mf, "na.action")
##   res$terms <- frames$terms
##   res$tJac <- ths$tJac
##   res$info <- with(res, {
##     data.frame("link" = link,
##                "threshold" = threshold,
##                "nobs" = nobs, 
##                "logLik" = formatC(logLik, digits=2, format="f"),
##                "AIC" = formatC(-2*logLik + 2*edf, digits=2,
##                  format="f"),
##                "niter" = paste(niter[1], "(", niter[2], ")", sep=""),
##                "max.grad" = formatC(maxGradient, digits=2,
##                  format="e") 
##                ## BIC is not part of output since it is not clear what
##                ## the no. observations are. 
##                )
##   })
##   ## add model.frame to results list?
##   if(model) res$model <- frames$mf
##   
##   return(res)
## }
## 
## clm.finalize <- function(fit, weights, alpha.names, beta.names,
##                          zeta.names=NULL)
## {   
##   nalpha <- length(alpha.names)
##   nbeta <- length(beta.names)
##   nzeta <- if(!is.null(zeta.names)) length(zeta.names) else 0L
##   stopifnot(length(fit$par) == nalpha + nbeta + nzeta)
##   
##   fit <- within(fit, {
##     alpha <- par[1:nalpha]
##     names(alpha) <- alpha.names
##     beta <- if(nbeta > 0) par[nalpha + 1:nbeta] else numeric(0)
##     zeta <- if(nzeta > 0) par[nalpha + nbeta + 1:nzeta]
##     else numeric(0) 
##     names(beta) <- beta.names
##     names(zeta) <- zeta.names
##     coefficients <- c(alpha, beta, zeta)
##     
##     names(gradient) <- names(coefficients)
##     dimnames(Hessian) <- list(names(coefficients),
##                               names(coefficients))
##     edf <- length(coefficients) ## estimated degrees of freedom
##     nobs <- sum(weights)
##     n <- length(weights)
##     fitted.values <- fitted
##     df.residual = nobs - edf
##     rm(list = c("par"))
##   })
##   class(fit) <- "clm"
##   return(fit)
## }
## 
## clm.model.frame <- function(mc, contrasts) {
## ### mc - the matched call
## ### contrasts - contrasts for the fixed model terms
## 
##   ## Extract model.frame(mf):
##   m <- match(c("formula", "data", "subset", "weights", "na.action",
##                "offset"), names(mc), 0)
##   mf <- mc[c(1, m)]
##   mf$drop.unused.levels <- TRUE
##   mf[[1]] <- as.name("model.frame")
##   mf <- eval(mf, envir = parent.frame(2))
## 
##   ## Extract model response:
##   y <- model.response(mf, "any") ## any storage mode
##   if(!is.factor(y))
##     stop("response needs to be a factor")
## 
##   ## Extract X:
##   terms <- attr(mf, "terms")
##   X <- model.matrix(terms, mf, contrasts)
##   n <- nrow(X)
##   ## Test for intercept in X:
##   Xint <- match("(Intercept)", colnames(X), nomatch = 0)
##   if(Xint <= 0) {
##     X <- cbind("(Intercept)" = rep(1, n), X)
##     warning("an intercept is needed and assumed")
##   } ## intercept in X is guarantied.
## 
##   ## Extract the weights and offset:
##   if(is.null(wts <- model.weights(mf))) wts <- rep(1, n)
##   if(is.null(off <- model.offset(mf))) off <- rep(0, n)
##   
##   ## check weights and offset:
##   if (any(wts <= 0))
##     stop(gettextf("non-positive weights are not allowed"))
## ### FIXME: Wouldn't it be usefull to just remove any observations with
## ### zero weights?
##   if(length(wts) && length(wts) != NROW(y))
##     stop(gettextf("number of weights is %d should equal %d
## (number of observations)", length(wts), NROW(y)))
##   if(length(off) && length(off) != NROW(y))
##     stop(gettextf("number of offsets is %d should equal %d
## (number of observations)", length(off), NROW(y)))
##   
##   ## return list of results:
##   res <- list(y = y, X = X, wts = as.double(wts),
##               off = as.double(off), mf = mf, terms = terms)
##   ## Note: X is with dimnames and an intercept is guarantied.
##   return(res)
## }
## ### Parameters in rho needed to optimize a clm:
## ### Initial:
## ### B1, B2, nxi, p, par(initial), o1, o2, wts, pfun, dfun, gfun,
## ### control 
## 
## ### Generated during fitting:
## ### eta1, eta2, pr, wtpr, p1, p2, dS.psi, dpi.psi
## 
## ### Variables needed to set starting values for a clm:
## ### Thresholds: ntheta/nlev.y/y, threshold
## ### Regression par: y, X, wts, off, link, OR regParStart <- rep(0, p)
## 
clm.newRho <-
  function(parent, y, X, weights, offset, tJac)
### Setting variables in rho: B1, B2, o1, o2, wts.
{
  rho <- new.env(parent = parent)

  ## Make B1, B2, o1, o2 based on y, X and tJac:
  ntheta <- nlevels(y) - 1
  n <- nrow(X)
  B2 <- 1 * (col(matrix(0, n, ntheta + 1)) == c(unclass(y)))
  rho$o1 <- c(1e5 * B2[, ntheta + 1]) - offset
  rho$o2 <- c(-1e5 * B2[,1]) - offset
  B1 <- B2[, -(ntheta + 1), drop = FALSE]
  B2 <- B2[, -1, drop = FALSE]
  ## adjust B1 and B2 for structured thresholds:
  rho$B1 <- B1 %*% tJac
  rho$B2 <- B2 %*% tJac
  ## update B1 and B2 with location effects (X):
  nbeta <- NCOL(X) - 1
  if(nbeta > 0) {
    rho$B1 <- cbind(rho$B1, -X[, -1, drop = FALSE])
    rho$B2 <- cbind(rho$B2, -X[, -1, drop = FALSE])
  }
  dimnames(rho$B1) <- NULL
  dimnames(rho$B2) <- NULL

  rho$fitted <- numeric(length = n)
  rho$wts <- weights

  return(rho)
}

clm.fit <-
  function(y, X, weights = rep(1, nrow(X)), offset = rep(0, nrow(X)),
           control = list(), start, 
           link = c("logit", "probit", "cloglog", "loglog", "cauchit"), 
           threshold = c("flexible", "symmetric", "equidistant"))
{
  ## Initial argument matching and testing:
  threshold <- match.arg(threshold)
  link <- match.arg(link)
  control <- do.call(clm.control, control)
  stopifnot(is.factor(y) &&
            is.matrix(X))
  stopifnot(length(y) == nrow(X) &&
            length(y) == length(weights) &&
            length(y) == length(offset))
  ## Test for intercept in X:
  Xint <- match("(Intercept)", colnames(X), nomatch = 0)
  if(Xint <= 0) {
    X <- cbind("(Intercept)" = rep(1, nrow(X)), X)
    warning("an intercept is needed and assumed")
  } ## intercept in X is guarantied.

  ## ensure X has full rank, generate threshold structure and set the
  ## rho environment:
  X <- drop.coef(X)
  ths <- makeThresholds(y, threshold)
  rho <- clm.newRho(parent.frame(), y=y, X=X, weights=weights,
                    offset=offset, tJac=ths$tJac)
  rho$clm.nll <- clm.nll
  rho$clm.grad <- clm.grad
  rho$clm.hess <- clm.hess

  ## set starting values for the clm:
  if(missing(start))
    start <- clm.start(y, X, has.intercept = TRUE, threshold = threshold)
  stopifnot(is.numeric(start) && 
            length(start) == (ths$nalpha + ncol(X) - 1) )
  rho$par <- start
  ## start cauchit models at the probit estimates:
  if(link == "cauchit" && is.null(match.call()$start)) {
    setLinks(rho, link="probit")
    fit <- try(clm.fit.env(rho), silent=TRUE) ## standard control values
    if(class(fit) == "try-error") 
      stop("Failed to find suitable starting values: please supply some")
    rho$par <- fit$par
  }

  ## Set pfun, dfun and gfun in rho:
  setLinks(rho, link)
  
  ## fit the clm:
  fit <- clm.fit.env(rho, control = control)
  return(fit)
### FIXME: should 'par' be 'coefficients' to allow coef(fit) etc.?
### FIXME: should fit contain 'vcov'?
}
  
clm.fit.env <-
  function(rho, control = list())
### The main work horse: Where the actual fitting of the clm goes on.
### Fitting the clm via Newton-Raphson with step halving. 

### -------- Assumes the existence of the following functions:
### clm.nll - negative log-likelihood
### clm.grad - gradient of nll wrt. par
### clm.hess - hessian of nll wrt. par
### Trace - for trace information

### Assumes that the Hessian is positive definite (such that step will
### be in direction of smaller gradient).
{
  control <- do.call(clm.control, control)
  stepFactor <- 1
  innerIter <- 0
  conv <- 1  ## Convergence flag
  message <- "iteration limit reached"
  nll <- rho$clm.nll(rho)
  if(!is.finite(nll))
    stop("Non-finite log-likelihood at starting value")
  gradient <- rho$clm.grad(rho)
  maxGrad <- max(abs(gradient))
  hessian <- rho$clm.hess(rho)
  if(control$trace > 0)
    Trace(iter=0, stepFactor, nll, maxGrad, rho$par, first=TRUE)
  
  ## Newton-Raphson algorithm:
  for(i in 1:control$maxIter) {
    if(maxGrad < control$gradTol) {
      message <- "max|gradient| < tol, so current iterate is probably solution"
      if(control$trace > 0)
        cat("\nOptimizer converged! ", "max|grad|:",
            maxGrad, message, fill = TRUE)
      conv <- 0
      break
    } ## end convergence test

    ## Compute Newton step and update parameters
    ## print(min(diag(hessian)))
    step <- .Call("La_dgesv", hessian, gradient, .Machine$double.eps,
                  PACKAGE = "base") ## solve H*step = g for 'step'
### FIXME: possibly use ginv instead for unidentifiable model?
    step <- as.vector(step)
    rho$par <- rho$par - stepFactor * step
    nllTry <- rho$clm.nll(rho)
    lineIter <- 0

    ## Step halfing if nll increases:
    while(nllTry > nll) {
      ## Step-half if nll increases or if nll is flat and maxGrad is
      ## not decreased. Thus if nll is flat, but maxGrad is decreased,
      ## the full step is justified:
      ## while(nllTry > nll || nllTry == nll && maxGradTry >= maxGrad) {
### Can nll increase while maxGrad decreases?

### It seems peculiar to ponder about a change in parameters that is
### so small that the likelihood is entirely flat. Still, the
### parameters can be determined with higher precision than what is
### visible from changes in nll.

### FIXME: Takes the step even if nllTry == nll - is this reasonable?
### This does happen with the rhyme data set. Apparently this is due
### to discreteness in the log-likelhood function, so even though a
### step can be taken that reduces max|grad|, the nll is flat and not
### reduced by the step.
### Could and should this be avoided by terminating if the step is not
### large enough? e.g., max(abs(step)) > 1e-8
### In a very well determined problem, the step can be very small -
### not if the gradient is very small, but if the Hessian is very
### large. Effectively par are at the MLEs.
      ## if(nllTry == nll) {
      ##   gradientTry <- clm.grad(rho)
      ##   if(maxGrad <= max(abs(clm.grad(rho)))) {
      ##     conv <- 3
      ##     message <- ""
      ## }
      stepFactor <- stepFactor/2
      rho$par <- rho$par + stepFactor * step
      nllTry <- rho$clm.nll(rho)
      lineIter <- lineIter + 1
      if(control$trace > 0) {
        cat("step halving:\n")
        Trace(i+innerIter, stepFactor, nll, maxGrad,
              rho$par, first = FALSE)
      }
      if(lineIter > control$maxLineIter){
        message <- "step factor reduced below minimum"
        conv <- 2
        break
      }
      innerIter <- innerIter + 1
    } ## end step halfing
    if(conv == 2) break

    ## Update nll, gradient, maxGrad and hessian:
    gradient <- rho$clm.grad(rho)
    maxGrad <- max(abs(gradient))
    hessian <- rho$clm.hess(rho)
    nll <- nllTry
    if(control$trace > 0) {
      if(control$trace > 1) {
        cat("\tgrad: ")
        cat(paste(formatC(gradient, digits=3, format="e")))
        cat("\n\tstep: ")
        cat(paste(formatC(step, digits=3, format="e")))
        cat("\n\teigen: ")
        cat(paste(formatC(eigen(hessian, symmetric=TRUE,
                                only.values=TRUE)$values, digits=3,
                          format="e")))
        cat("\n")
      }
      Trace(iter=i, stepFactor, nll,
            maxGrad, rho$par, first = FALSE)
    }
    ## Double stepFactor if needed:
    stepFactor <- min(1, 2 * stepFactor)
  } ## end Newton iterations

### FIXME: warn if fitted values are close to 0 or 1?
  
  if(conv > 0) { ## optimization failed
    if(control$trace > 0) cat(message, fill = TRUE)
    ## stop(gettextf("optimization failed: %s", message), call. = FALSE)
    warning(gettextf("optimization failed: %s", message), call. = FALSE)
  }
  
  ## return results:
  res <- list(par = rho$par,
              gradient = as.vector(gradient),
              Hessian = rho$clm.hess(rho), ## ensure hessian is evaluated
              ## at optimum
              logLik = -nll, 
              convergence = conv,
              ## 0: successful convergence
              ## 1: iteration limit reached
              ## 2: step factor reduced below minimum
              message = message,
              maxGradient = maxGrad,
              niter = c(outer = i-1, inner = innerIter),
              fitted = rho$fitted)
  return(res)
}

clm.nll <- function(rho) { ## negative log-likelihood
### For linear models
  with(rho, {
    eta1 <- drop(B1 %*% par) + o1
    eta2 <- drop(B2 %*% par) + o2
    fitted <- pfun(eta1) - pfun(eta2)
    if(all(fitted > 0))
### NOTE: Need test here because some fitted <= 0 if thresholds are
### not ordered increasingly.
### It is assumed that 'all(is.finite(pr)) == TRUE' 
      -sum(wts * log(fitted))
    else Inf
  })
}

clm.grad <- function(rho) { ## gradient of the negative log-likelihood
### return: vector of gradients
### For linear models
  with(rho, {
    p1 <- dfun(eta1)
    p2 <- dfun(eta2)
    wtpr <- wts/fitted
    dpi.psi <- B1 * p1 - B2 * p2
    -crossprod(dpi.psi, wtpr)
### NOTE: It is assumed that all(fitted > 0) == TRUE and that
### all(is.finite(c(p1, p2))) == TRUE
  })
}

clm.hess <- function(rho) { ## hessian of the negative log-likelihood
### return Hessian matrix
### For linear models
  with(rho, {
    dg.psi <- crossprod(B1 * gfun(eta1) * wtpr, B1) -
      crossprod(B2 * gfun(eta2) * wtpr, B2)
    -dg.psi + crossprod(dpi.psi, (dpi.psi * wtpr / fitted))
### NOTE: It is assumed that all(fitted > 0) == TRUE and that
### all(is.finite(c(g1, g2))) == TRUE
  })
}

clm.control <-
  function(method = c("Newton", "model.frame", "ucminf", "nlminb",
             "optim"), ...,  trace = 0, maxIter = 100, gradTol = 1e-6,
           maxLineIter = 15) 
{
### FIXME: change Newton to clm.fit?
  method <- match.arg(method)
  
  if(!all(is.numeric(c(maxIter, gradTol, maxLineIter))))
    stop("maxIter, gradTol, and maxLineIter should all be numeric")
  
  ## method <- match.arg(method)
  ctrl <- list(method = method,
               trace = as.integer(trace),
               maxIter = as.integer(maxIter),
               gradTol = as.numeric(gradTol),
               maxLineIter = as.integer(maxLineIter))
  if(method %in% c("ucminf", "nlminb", "optim"))
    ctrl$ctrl <- list(trace = as.integer(abs(trace)), ...)

  return(ctrl)
}

