## This file contains functions for fitting CLMMs with a single simple
## RE term (ssr).



rho.clm2clmm.ssr <- function(rho, grList, ctrl)
### Version of rho.clm2clmm that is set up to use the C
### implementations of Laplace, AGQ and GHQ for a single random
### effect. 
{
  rho$ctrl <- ctrl
  rho$grFac <- grList[[1]]
  rho$sigma <- rep(1, nrow(rho$B1))
  rho$lambda <- 0
  rho$nlev <- as.vector(sapply(grList, nlevels))
  rho$random.names <- sapply(grList, levels)
  rho$tau.names <- names(grList)
  rho$nrandom <- sum(rho$nlev) ## no. random effects
  rho$Niter <- 0L
  rho$u <- rho$uStart <- rep(0, rho$nrandom)
  rho$linkInt <- switch(rho$link,
                        logit = 1L,
                        probit = 2L,
                        cloglog = 3L,
                        loglog = 4L,
                        cauchit = 5L)
}

set.AGQ <- function(rho, nAGQ) {
  rho$nAGQ <- nAGQ
  if(nAGQ %in% c(0L, 1L)) return(invisible())
  ghq <- gauss.hermite(abs(nAGQ))
  rho$ghqns <- ghq$nodes
  rho$ghqws <-
    if(nAGQ > 0) ghq$weights ## AGQ
    else log(ghq$weights) + (ghq$nodes^2)/2 ## GHQ
}

clmm.fit.ssr <- 
  function(rho, control = list(), Hess = FALSE)
### Fit a clmm with a single simple random effects term using AGQ, GHQ
### or Laplace.
{
  ## Set appropriate objective function:
  obj.fun <-
    if(rho$nAGQ < 0) getNGHQ.ssr
    else if(rho$nAGQ > 1) getNAGQ.ssr
    else getNLA.ssr ## nAGQ %in% c(0, 1)
  
  ## Fit the model with Laplace:
  fit <- ucminf(rho$par, function(par) obj.fun(rho, par),
                control = control)
  ## rho$par <- fit$par
  
  ## Ensure random mode estimation at optimum:
  nllBase.uC(rho)
  update.uC(rho)

  ## Format ranef modes and condVar:
  ranef <- rho$u * rho$tau
  condVar <- 1/rho$D * rho$tau^2
  names(ranef) <- names(condVar) <- rho$random.names
  ranef <- list(ranef)
  condVar <- list(condVar)
  names(ranef) <- names(condVar) <- rho$tau.names

  ## Prepare list of results:
  res <- list(coefficients = fit$par,
              optRes = fit,
              logLik = -fit$value,
              fitted.values = rho$fitted,
              ranef = ranef,
              condVar = condVar,
              Niter = rho$Niter)
  
  ## Compute hessian?
  if(Hess)
    res$Hessian <-
      hessian(function(par) getNLA.ssr(rho, par), x = fit$par,
              method.args = list(r = 2, show.details = TRUE))
  
  return(res)
}

getNGHQ.ssr <- function(rho, par) {
### negative log-likelihood by standard Gauss-Hermite quadrature
### implemented in C:
  if(!missing(par))
    rho$par <- par
  nllBase.uC(rho) ## Update tau, eta1Fix and eta2Fix
  with(rho, {
    .C("getNGHQ",
       nll = double(1),
       as.integer(grFac),
       as.double(tau),
       as.double(eta1Fix),
       as.double(eta2Fix),
       as.double(o1),
       as.double(o2),
       as.double(sigma),
       as.double(wts),
       length(sigma),
       length(uStart),
       as.double(ghqns),
       as.double(ghqws),
       as.integer(abs(nAGQ)),
       as.integer(linkInt),
       as.double(ghqns * tau),
       as.double(lambda))$nll
  })
}

getNAGQ.ssr <- function(rho, par) {
### negative log-likelihood by adaptive Gauss-Hermite quadrature
### implemented in C:
  if(!missing(par))
    rho$par <- par
  if(!update.uC(rho)) return(Inf)
  if(any(rho$D < 0)) return(Inf)
  with(rho, {
    .C("getNAGQ",
       nll = double(1),
       as.integer(grFac),
       as.double(tau),
       as.double(eta1Fix),
       as.double(eta2Fix),
       as.double(o1),
       as.double(o2),
       as.double(sigma),
       as.double(wts),
       length(sigma),
       length(uStart),
       as.double(ghqns),
       as.double(log(ghqws)),
       as.double(ghqns^2),
       as.double(u),
       as.double(D),
       as.integer(abs(nAGQ)),
       as.integer(linkInt),
       as.double(lambda))$nll
  })
}

getNLA.ssr <- function(rho, par) {
### negative log-likelihood by the Laplace approximation
### (with update.u2 in C or R):
  if(!missing(par)) rho$par <- par
  if(!update.uC(rho)) return(Inf)
  if(any(rho$D < 0)) return(Inf)
  logDetD <- sum(log(rho$D)) - rho$nrandom * log(2*pi)
  rho$negLogLik + logDetD / 2
}

nllBase.uC <- function(rho) {
### updates tau, eta1Fix and eta2Fix given new parameter values
  with(rho, {
    tau <- exp(par[nalpha + nbeta + 1:ntau])
    eta1Fix <- drop(B1 %*% par[1:(nalpha + nbeta)])
    eta2Fix <- drop(B2 %*% par[1:(nalpha + nbeta)])
  })
  return(invisible())
}

update.uC <- function(rho) {
### C-implementation of NR-algorithm.
  nllBase.uC(rho) ## update: tau, eta1Fix, eta2Fix
  fit <- with(rho, {
    .C("NRalgv3",
       as.integer(ctrl$trace),
       as.integer(ctrl$maxIter),
       as.double(ctrl$gradTol),
       as.integer(ctrl$maxLineIter),
       as.integer(grFac), ## OBS
       as.double(tau), # stDev
       as.double(o1),
       as.double(o2),
       as.double(eta1Fix),
       as.double(eta2Fix),
       as.double(sigma), ## rep(1, n)
       as.integer(linkInt), ## 
       as.double(wts), ## pre. weights
       u = as.double(uStart),
       fitted = as.double(fitted), ## pre. pr
       funValue = double(1),
       gradValues = as.double(uStart),
       hessValues = as.double(uStart),
       length(fitted),
       length(uStart),
       maxGrad = double(1),
       conv = 0L,
       as.double(lambda), ## 
       Niter = as.integer(Niter) ## OBS
       )[c("u", "fitted", "funValue", "gradValues",
           "hessValues", "maxGrad", "conv", "Niter")] })
    ## Get message:
  message <- switch(as.character(fit$conv),
                    "1" = "max|gradient| < tol, so current iterate is probably solution",
                    "0" = "Non finite negative log-likelihood",
                    "-1" = "iteration limit reached when updating the random effects",
                    "-2" = "step factor reduced below minimum when updating the random effects")
  ## Check for convergence and report warning/error:
  if(rho$ctrl$trace > 0 && fit$conv == 1)
    cat("\nOptimizer converged! ", "max|grad|:",
        fit$maxGrad, message, fill = TRUE)
  if(fit$conv != 1 && rho$ctrl$innerCtrl == "warnOnly")
    warning(message, "\n  at iteration ", rho$Niter)
  else if(fit$conv != 1 && rho$ctrl$innerCtrl == "giveError")
    stop(message, "\n  at iteration ", rho$Niter)
  ## Store values and return:
  rho$Niter <- fit$Niter
  rho$fitted <- fit$fitted
  rho$u <- fit$u
  rho$D <- fit$hessValue
  rho$gradient <- fit$gradValue
  if(!is.finite(rho$negLogLik <- fit$funValue))
    return(FALSE)
  return(TRUE)
}

