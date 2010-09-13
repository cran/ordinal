clm <-
  function(formula, scale, nominal, data, weights, start, subset,
           doFit = TRUE, na.action, contrasts, model = TRUE,
           control = list(),
           link = c("logit", "probit", "cloglog", "loglog", "cauchit"), 
           threshold = c("flexible", "symmetric", "equidistant"), ...)
### deliberately no offset argument - include offset in the relevant
### formula/scale instead. 
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
  
  ## identify model as 'simple clm' or 'extended clm':
  Class <- if(!missing(scale) || link == "cauchit")
    c("eclm", "clm") else c("sclm", "clm")
  
  ## Compute: y, X, wts, off, mf:
  frames <- eclm.model.frame(mc, contrasts)

  ## Compute the transpose of the Jacobian for the threshold function,
  ## tJac and the names of the threshold parameters, alpha.names:
  frames$ths <- makeThresholds(frames$y, threshold)

  ## Return model.frame?
  if(control$method == "model.frame") return(frames)
  
  ## Test column rank deficiency and possibly drop some
  ## parameters. Also set the lists aliased and coef.names:
  ## (X is with intercept at this point.)
  frames <- drop.cols(frames, silent=TRUE)
### Note: intercept could be dropped from X and S in drop.cols?
### Currently they are dropped in eclm.newRho instead.
  
  ## Set envir rho with variables: B1, B2, o1, o2, wts, fitted:
  rho <- with(frames, {
    eclm.newRho(parent.frame(), y=y, X=X, NOM=frames$NOM, S=frames$S, 
                weights=wts, offset=off, S.offset=frames$S.off,
                tJac=ths$tJac) })

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
  
  ## Set pfun, dfun and gfun in rho:
  setLinks(rho, link)

  ## Possibly return the environment rho without fitting:
  if(!doFit) return(rho)
  
  ## Fit the clm:
  if(control$method == "Newton") {
    fit <- if("eclm" %in% Class) clm.fit.NR(rho, control) else
    clm.fit.env(rho, control)
  } else
  fit <- clm.fit.optim(rho, control$method, control$ctrl)
### NOTE: we could add arg non.conv = c("error", "warn", "message") to
### allow non-converged fits to be returned.

  ## Modify and return results:
  res <- eclm.finalize(fit, weights=frames$wts,
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
  res$na.action <- attr(frames$mf, "na.action")
  res$terms <- frames$terms
  res$xlevels <- .getXlevels(res$terms, frames$mf)
  if(!is.null(frames$NOM)) {
    res$nom.contrasts <- attr(frames$NOM, "contrasts")
    res$nom.terms <- frames$nom.terms
    res$nom.xlevels <- .getXlevels(res$nom.terms, frames$mf)
  }
  if(!is.null(frames$S)) {
    res$S.contrasts <- attr(frames$S, "contrasts")
    res$S.terms <- frames$S.terms
    res$S.xlevels <- .getXlevels(res$S.terms, frames$mf)
  }
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
  class(res) <- Class
  ## add model.frame to results list?
  if(model) res$model <- frames$mf
  
  return(res)
}

eclm.model.frame <- function(mc, contrasts) {
### mc - the matched call
### contrasts - contrasts for the model terms

  ## Collect all variables in a full formula:
  forms <- list(mc$formula)
  if(!is.null(mc$scale)) forms$scale <- mc$scale
  if(!is.null(mc$nominal)) forms$nominal <- mc$nominal
  ## ensure formula, scale and nominal are formulas:
  forms <- lapply(forms, function(x) {
    try(formula(deparse(x), env = parent.frame(2)), silent=TRUE) })
  if(any(sapply(forms, function(f) class(f) == "try-error")))
    stop("unable to interpret 'formula', 'scale' or 'nominal'")

  fullForm <- do.call(getFullForm, forms)
  ## set environment of 'fullForm' to the environment of 'formula': 
  form.envir <- environment(eval(mc$formula)) 
  environment(fullForm) <- form.envir

  ## Extract the full model.frame(mf):
  m <- match(c("data", "subset", "weights", "na.action"),
             names(mc), 0) 
  mf <- mc[c(1, m)]
  mf$formula <- fullForm
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  fullmf <- eval(mf, envir = parent.frame(2))
  mf$na.action <- "na.pass" ## filter NAs by hand below

  ## Extract X:
  ## get X from fullmf to handle e.g., NAs and subset correctly
  mf$formula <- mc$formula
  X.mf <- eval(mf, envir = parent.frame(2))
  ## X.mf <- eval(mf)
  X.terms <- attr(X.mf, "terms")
  X <- model.matrix(X.terms, fullmf, contrasts)
  n <- nrow(X)
  ## Test for intercept in X:
  Xint <- match("(Intercept)", colnames(X), nomatch = 0)
  if(Xint <= 0) {
    X <- cbind("(Intercept)" = rep(1, n), X)
    warning("an intercept is needed and assumed in 'formula'",
            call.=FALSE)
  } ## intercept in X is guaranteed.
  off <- getOffset(X.mf)
  ## Filter NAs if any:
  if(!is.null(naa <- na.action(fullmf)))
    off <- off[-naa]

  ## Extract model response:
  y <- model.response(fullmf, "any") ## any storage mode
  if(!is.factor(y)) stop("response needs to be a factor", call.=FALSE)
  
  ## list of results:
  res <- list(y=y, X=X, wts=getWeights(fullmf), off=off,
              mf=fullmf, terms=X.terms)

  ## Extract S (design matrix for the scale effects):
  if(!is.null(mc$scale)) {
    mf$formula <- forms$scale
    S.mf <- eval(mf, envir = parent.frame(2))
    if(!is.null(model.response(S.mf)))
      stop("response not allowed in 'scale'", call.=FALSE)
    res$S.terms <- attr(S.mf, "terms")
    S <- model.matrix(res$S.terms, fullmf, contrasts)
    ## Test for intercept in S:
    Sint <- match("(Intercept)", colnames(S), nomatch = 0)
    if(Sint <= 0) {
      S <- cbind("(Intercept)" = rep(1, n), S)
      warning("an intercept is needed and assumed in 'scale'",
              call.=FALSE)
    } ## intercept in S is guaranteed.
    res$S <- S
    res$S.off <- getOffset(S.mf)
    ## Filter NAs if any:
    if(!is.null(naa <- na.action(fullmf)))
      res$S.off <- res$S.off[-naa]
  }
  
  ## Extract NOM (design matrix for the nominal effects):
  if(!is.null(mc$nominal)) {
    mf$formula <- forms$nominal
    nom.mf <- eval(mf, envir = parent.frame(2))
    ## nom.mf <- eval(mf)
    if(!is.null(model.response(nom.mf)))
      stop("response not allowed in 'nominal'", call.=FALSE)
    if(!is.null(model.offset(nom.mf)))
      stop("offset not allowed in 'nominal'", call.=FALSE)
    res$nom.terms <- attr(nom.mf, "terms")
    NOM <- model.matrix(res$nom.terms, fullmf, contrasts)
    NOMint <- match("(Intercept)", colnames(NOM), nomatch = 0)
    if(NOMint <= 0) {
      NOM <- cbind("(Intercept)" = rep(1, n), NOM)
      warning("an intercept is needed and assumed in 'nominal'",
              call.=FALSE)
    } ## intercept in NOM is guarantied.
    res$NOM <- NOM
  }

  ## return results:
  return(res)
  ## Note: X, S and NOM are with dimnames and intercepts are
  ## guaranteed. They may be column rank defecient. 
}

eclm.newRho <-
  function(parent=parent.frame(), y, X, NOM=NULL, S=NULL, weights,
           offset, S.offset=NULL, tJac) 
### Setting variables in rho: B1, B2, o1, o2, wts.
{
  rho <- new.env(parent = parent)
  ## Make B1, B2, o1, o2 based on y, X and tJac:
  ## rho <- list2env(getB(y=y, NOM=NOM, X=X, offset=offset,
  ## tJac=tJac), parent=parent) 
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
  ## update B1 and B2 with nominal effects:
  if(NCOL(NOM) > 1) { ## !is.null(NOM) && ncol(NOM) > 1) {
    ## if !is.null(NOM) and NOM is more than an intercept:
    LL1 <- lapply(1:ncol(NOM), function(x) rho$B1 * NOM[,x])
    rho$B1 <- do.call(cbind, LL1)
    LL2 <- lapply(1:ncol(NOM), function(x) rho$B2 * NOM[,x])
    rho$B2 <- do.call(cbind, LL2)
  }
  ## update B1 and B2 with location effects (X):
  nbeta <- NCOL(X) - 1
  if(nbeta > 0) {
    rho$B1 <- cbind(rho$B1, -X[, -1, drop = FALSE])
    rho$B2 <- cbind(rho$B2, -X[, -1, drop = FALSE])
  }
  dimnames(rho$B1) <- NULL
  dimnames(rho$B2) <- NULL
  rho$n.psi <- ncol(rho$B1) ## no. linear model parameters
  rho$k <- 0
  ## there may be scale offset without scale predictors:
  rho$sigma <- rho$Soff <-
    if(is.null(S.offset)) rep(1, n) else exp(S.offset) 
  ## save scale model:
  if(!is.null(S)) {
    rho$S <- S[, -1, drop=FALSE]
    dimnames(rho$S) <- NULL
    rho$k <- ncol(rho$S) ## no. scale parameters
  }
  rho$has.scale <- ## TRUE if scale has to be considered.
    if(!is.null(S) || any(S.offset != 0)) TRUE else FALSE
  ## initialize fitted values and weights:
  rho$fitted <- numeric(length = n)
  rho$wts <- weights
  ## return:
  return(rho)
}

clm.fit.optim <-
  function(rho, method = c("ucminf", "nlminb", "optim"), control=list()) 
{
  method <- match.arg(method)
  ## optimize the likelihood:
  optRes <-
    switch(method,
           "nlminb" = nlminb(rho$par,
             function(par) eclm.nll(rho, par),
             function(par) eclm.grad2(rho, par),
             control=control),        
           "ucminf" = ucminf(rho$par, 
             function(par) eclm.nll(rho, par),
             function(par) eclm.grad2(rho, par),
             control=control),             
           "optim" = optim(rho$par,
             function(par) eclm.nll(rho, par),
             function(par) eclm.grad2(rho, par),
             method="BFGS",
             control=control),
           )
  ## save results:
  rho$par <- optRes[[1]]
  res <- list(par = rho$par,
              logLik = -eclm.nll(rho),
              gradient = eclm.grad(rho),
              Hessian = eclm.hess(rho),
              fitted = rho$fitted)
  res$maxGradient = max(abs(res$gradient))
  res$optRes <- optRes
  res$niter <- switch(method, "nlminb" = optRes$evaluations,
                      "ucminf" = c(optRes$info["neval"], 0),
                      "optim" = optRes$counts)
  res$convergence <-
    switch(method, "nlminb" = optRes$convergence,
           "ucminf" = optRes$convergence,
           "optim" = optRes$convergence)
  
  return(res)
}

eclm.nll <- function(rho, par) {
  if(!missing(par)) rho$par <- par
  with(rho, {
    if(k > 0) 
      sigma <- Soff * exp(drop(S %*% par[n.psi + 1:k]))
### NOTE: we have to divide by sigma even if k=0 since there may be an
### offset but no predictors in the scale model:
    eta1 <- (drop(B1 %*% par[1:n.psi]) + o1)/sigma
    eta2 <- (drop(B2 %*% par[1:n.psi]) + o2)/sigma
    fitted <- pfun(eta1) - pfun(eta2)
    if(all(is.finite(fitted)) && all(fitted > 0))
### NOTE: Need test here because some fitted <= 0 if thresholds are
### not ordered increasingly.
      -sum(wts * log(fitted))
    else Inf
  })
}

eclm.grad <- function(rho) {
### requires that eclm.nll has been called prior to
### eclm.grad. 
  with(rho, {
    p1 <- dfun(eta1)
    p2 <- dfun(eta2)
    wtpr <- wts/fitted
    C2 <- B1*p1/sigma - B2*p2/sigma
    if(k <= 0) return(-crossprod(C2, wtpr))
    C3 <- -(eta1 * p1 - eta2 * p2) * S
    return(-crossprod(cbind(C2, C3), wtpr))
### NOTE: C2 and C3 are used by eclm.hess
  })
}

eclm.grad2 <- function(rho, par) {
### does not require that eclm.nll has been called prior to
### eclm.grad.
  eclm.nll(rho, par)
  eclm.grad(rho)
}

eclm.hess <- function(rho) {
### requires that eclm.grad has been called prior to this.
  with(rho, {
    g1 <- gfun(eta1)
    g2 <- gfun(eta2)
    wtprpr <- wtpr/fitted ## Phi3
    dg.psi <- crossprod(B1 * gfun(eta1) * wtpr / sigma^2, B1) -
      crossprod(B2 * gfun(eta2) * wtpr / sigma^2, B2)
    ## upper left:
    D <- dg.psi - crossprod(C2, (C2 * wtprpr))
    if(k <= 0) return(-D) ## no scale predictors
    ## upper right (lower left transpose):
    wtprsig <- wtpr/sigma
    epg1 <- p1 + g1*eta1
    epg2 <- p2 + g2*eta2
    Et <- crossprod(B1, -wtprsig * epg1 * S) -
      crossprod(B2, -wtprsig * epg2 * S) -
        crossprod(C2, wtprpr * C3)  
    ## lower right:
    F <- -crossprod(S, wtpr * ((eta1*p1 - eta2*p2)^2 / fitted -
                               (eta1*epg1 - eta2*epg2)) * S)
    ## combine and return hessian:
    H <- rbind(cbind(D    , Et),
               cbind(t(Et), F))
    return(-H)
  })
}

clm.fit.NR <-
  function(rho, control = list())
### The main work horse: Where the actual fitting of the clm goes on.
### Fitting the clm via Newton-Raphson with step halving. 

### -------- Assumes the existence of the following functions:
### eclm.nll - negative log-likelihood
### eclm.grad - gradient of nll wrt. par
### eclm.hess - hessian of nll wrt. par
### Trace - for trace information
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
    ## e.val <- eigen(hessian, symmetric=TRUE, only.values=TRUE)$values
    min.ev <- min(eigen(hessian, symmetric=TRUE, only.values=TRUE)$values)
    tol <- 1e-5
    if(min.ev < 0) {
### if(min.ev <= tol) {
      inflate <- abs(min.ev) + 1 # + tol
      ## inflate <- ceiling(abs(min.ev + tol))
      hessian <- hessian + diag(inflate, nrow(hessian))
      if(control$trace > 0)
        cat(paste("X is not positive definite, inflating diagonal with",
                  formatC(inflate, digits=5, format="fg"), "\n"))
    }
    step <- .Call("La_dgesv", hessian, gradient, .Machine$double.eps,
                  PACKAGE = "base") ## solve H*step = g for 'step'
    step <- as.vector(step)
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
      Trace(iter=i, stepFactor, nll,
            maxGrad, rho$par, first = FALSE)
    }
    ## Double stepFactor if needed:
    stepFactor <- min(1, 2 * stepFactor)
  } ## end Newton iterations
  if(conv > 0)  ## optimization failed
    ## if(control$trace > 0) cat(message, fill = TRUE)
    ## stop(gettextf("optimization failed: %s", message), call. = FALSE)
    warning(gettextf("optimization failed: %s", message), call. = FALSE)
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
