clm.newRho <-
  function(parent, y, X, weights, offset, tJac)
### Set variables in rho: B1, B2, o1, o2 and wts.
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
  function(y, X, S, N, weights = rep(1, nrow(X)),
           offset = rep(0, nrow(X)), S.offset = rep(0, nrow(X)),
           control = list(), start, 
           link = c("logit", "probit", "cloglog", "loglog", "cauchit"), 
           threshold = c("flexible", "symmetric", "equidistant"))
### This function basically does the same as clm, but without setting
### up the model matrices from formulae, and with minimal post
### processing after parameter estimation. 
{
  ## Initial argument matching and testing:
  threshold <- match.arg(threshold)
  link <- match.arg(link)
  control <- do.call(clm.control, control)
  if(missing(y)) stop("please specify y")
  if(missing(X)) X <- cbind("(Intercept)" = rep(1, length(y)))
  stopifnot(is.factor(y) &&
            is.matrix(X))
  stopifnot(length(y) == nrow(X) &&
            length(y) == length(weights) &&
            length(y) == length(offset) && 
            length(y) == length(S.offset))
  frames <- list(y=y, X=X)
  ## S and N are optional:
  if(!missing(S)) {
    frames$S <- S
    stopifnot(is.matrix(S) &&
              length(y) == nrow(S))
  }
  if(!missing(N)) {
    frames$NOM <- N
    stopifnot(is.matrix(N) &&
              length(y) == nrow(N))
  }

  ## Identify model as 'simple clm' or 'extended clm':
  Class <- if(any(S.offset != 0) || !missing(S) || link == "cauchit")
    c("eclm", "clm") else c("sclm", "clm")
  ## get  threshold structure:
  frames$ths <- makeThresholds(y, threshold)
  ## test for column rank deficiency in design matrices:
  frames <- drop.cols(frames, silent=TRUE)
  
  ## Set envir rho with variables: B1, B2, o1, o2, wts, fitted...:
  rho <- eclm.newRho(parent.frame(), y=frames$y, X=frames$X,
                     NOM=frames$NOM, S=frames$S, weights=weights,
                     offset=offset, S.offset=S.offset,
                     tJac=frames$ths$tJac)

  ## Set appropriate logLik and deriv functions in rho:
  if("eclm" %in% Class) {
    rho$clm.nll <- eclm.nll
    rho$clm.grad <- eclm.grad
    rho$clm.hess <- eclm.hess
  } else {
    rho$clm.nll <- clm.nll
    rho$clm.grad <- clm.grad
    rho$clm.hess <- clm.hess
  }

  ## Set starting values for the parameters:
  start <- set.start(rho, start=start, get.start=missing(start),
                     threshold=threshold, link=link, frames=frames)
  rho$par <- as.vector(start) ## remove attributes
  
  ## Set inverse link function and its derivatives (pfun, dfun and
  ## gfun) in rho: 
  setLinks(rho, link)

  ## Fit the model:
  if(control$method == "Newton") {
    fit <- if("eclm" %in% Class) clm.fit.NR(rho, control) else
    clm.fit.env(rho, control)
  } else
  fit <- clm.fit.optim(rho, control$method, control$ctrl)

  ## Format and return the fit:
  fit$coef.names <- frames$coef.names
  fit$aliased <- lapply(frames$aliased, as.logical)
  if(control$method == "Newton" &&
     !is.null(start.iter <- attr(start, "start.iter")))
    fit$niter <- fit$niter + start.iter

  return(fit)
### FIXME: should 'par' be 'coefficients' to allow coef(fit) etc.?
### FIXME: should fit contain 'vcov'?
}
  
clm.fit.env <-
  function(rho, control = list())
### The main work horse: Where the actual fitting of the clm goes on.
### Fitting the clm via Newton-Raphson with step halving. 

### -------- Assumes the existence of the following functions in rho: 
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
    ## step <- .Call("La_dgesv", hessian, gradient, .Machine$double.eps,
    ##               PACKAGE = "base") ## solve H*step = g for 'step'
##     step <- try(.Call("La_dgesv", hessian, gradient, .Machine$double.eps,
##                       PACKAGE = "base"), silent=TRUE)
##     if(class(step) == "try-error") {
##       warning("Unable to compute Newton step:
## possibly some parameter estimates are infinite or the model is not
## identifiable")
##       break
##     }
    ## step <- as.vector(step)
    step <- as.vector(solve(hessian, gradient))
    rho$par <- rho$par - stepFactor * step
    nllTry <- rho$clm.nll(rho)
    lineIter <- 0

    ## Step halfing if nll increases:
    while(nllTry > nll) {
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
      Trace(iter=i+innerIter, stepFactor, nll,
            maxGrad, rho$par, first = FALSE)
    }
    ## Double stepFactor if needed:
    stepFactor <- min(1, 2 * stepFactor)
  } ## end Newton iterations

  if(conv > 0) { ## optimization failed
    if(control$trace > 0) cat(message, fill = TRUE)
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
  })
### NOTE: getFitted is not found from within rho, so we have to
### evalueate it outside of rho
  rho$fitted <- getFittedC(rho$eta1, rho$eta2, rho$link)
  if(all(rho$fitted > 0))
### NOTE: Need test here because some fitted <= 0 if thresholds are
### not ordered increasingly.
### It is assumed that 'all(is.finite(pr)) == TRUE' 
    -sum(rho$wts * log(rho$fitted))
  else Inf
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
  method <- match.arg(method)
  
  if(!all(is.numeric(c(maxIter, gradTol, maxLineIter))))
    stop("maxIter, gradTol, and maxLineIter should all be numeric")
  
  ctrl <- list(method = method,
               trace = as.integer(trace),
               maxIter = as.integer(maxIter),
               gradTol = as.numeric(gradTol),
               maxLineIter = as.integer(maxLineIter))
  if(method %in% c("ucminf", "nlminb", "optim"))
    ctrl$ctrl <- list(trace = as.integer(abs(trace)), ...)

  return(ctrl)
}

