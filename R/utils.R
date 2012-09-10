getFitted <- function(eta1, eta2, pfun, ...) {
  ## eta1, eta2: linear predictors
  ## pfun: cumulative distribution function
  ##
  ## Compute fitted values while maintaining high precision in the
  ## result - if eta1 and eta2 are both large, fitted is the
  ## difference between two numbers very close to 1, which leads to
  ## imprecision and potentially errors.
  ##
  ## Note that (eta1 > eta2) always holds, hence (eta2 > 0) happens
  ## relatively rarely.
  k2 <- eta2 > 0
  fitted <- pfun(eta1) - pfun(eta2)
  fitted[k2] <- pfun(eta2[k2], lower.tail=FALSE) -
    pfun(eta1[k2], lower.tail=FALSE)
  fitted
}

getFittedC <-
  function(eta1, eta2,
           link = c("logit", "probit", "cloglog", "loglog", "cauchit",
             "Aranda-Ordaz", "log-gamma"), lambda=1)
### Same as getFitted only this is implemented in C and handles all
### link functions including the flexible ones.
{
  link <- match.arg(link)
  .Call("get_fitted", eta1, eta2, link, lambda)
}

getWeights <- function(mf) {
### mf - model.frame
  n <- nrow(mf)
  if(is.null(wts <- model.weights(mf))) wts <- rep(1, n)
  if (any(wts <= 0))
    stop(gettextf("non-positive weights are not allowed"),
         call.=FALSE)
### NOTE: We do not remove observations where weights == 0, because
### that could be a somewhat surprising behaviour. It would also
### require that the model.frame be evaluated all over again to get
### the right response vector with the right number of levels.
  if(length(wts) && length(wts) != n)
    stop(gettextf("number of weights is %d should equal %d (number of observations)",
                  length(wts), n), call.=FALSE)
  return(as.double(wts))
}

getOffset <- function(mf) {
### mf - model.frame
  n <- nrow(mf)
  if(is.null(off <- model.offset(mf))) off <- rep(0, n)
  if(length(off) && length(off) != n)
    stop(gettextf("number of offsets is %d should equal %d (number of observations)",
                  length(off), n), call.=FALSE)
  return(as.double(off))
}

getFullForm <- function(form, ..., envir=parent.frame()) {
### collect terms in several formulas in a single formula
### sets the environment of the resulting formula to envir. 
  forms <- list(...)
  if(lf <- length(forms)) {
    rhs <- character(0)
    ## Collect rhs terms in a single vector of rh-sides:
    for(i in 1:lf) {
      rhs <- c(rhs, Deparse(forms[[i]][[2]]))
      if(length(forms[[i]]) >= 3)
        rhs <- c(rhs, Deparse(forms[[i]][[3]]))
    }
    ## add '+' inbetween terms:
    rhs <- paste(rhs, collapse=" + ")
    ## combine if 'deparse(form)' is a (long) vector:
    form2 <- paste(deparse(form, width.cutoff=500L), collapse=" ")
    ## combine form2 and rhs into a single string:
    form <- paste(form2, rhs, sep=" + ")
  }
  return(as.formula(form, env=envir))
}

## getFullForm <- function(form, ..., envir=parent.frame()) {
## ### collect terms in several formulas in a single formula (on the rhs)
## ### sets the environment of the resulting formula to envir. 
##   forms <- list(form, ...)
##   allVars <- unlist(sapply(forms, all.vars))
##   rhs <- paste(allVars, collapse=" + ")
##   form <- paste("~", rhs)
##   return(as.formula(form, env=envir))
## }

drop.coef2 <- function(X, tol = 1e-7, silent = FALSE, test.ans = FALSE)
### works if ncol(X) >= 0 and nrow(X) >= 0 
{
  ## test and match arguments:
  stopifnot(is.matrix(X))
  silent <- as.logical(silent)[1]
  aliased <- rep.int(0, ncol(X))
  ## perform the qr-decomposition of X using LINPACK methods:
  qr.X <- qr(X, tol = tol, LAPACK = FALSE)
  if(qr.X$rank == ncol(X)) {
    ## return X if X has full column rank
    attr(X, "aliased") <- aliased
    attr(X, "orig.colnames") <- colnames(X)
    return(X)
  }
  if(!silent) ## message the no. dropped columns:
    message(gettextf("design is column rank deficient so dropping %d coef",
                     ncol(X) - qr.X$rank))
  ## return the columns correponding to the first qr.x$rank pivot
  ## elements of X:
  newX <- X[, qr.X$pivot[1:qr.X$rank], drop = FALSE]
  sel <- qr.X$pivot[-(1:qr.X$rank)]
  aliased[sel] <- 1
  attr(newX, "aliased") <- aliased
  attr(newX, "orig.colnames") <- colnames(X)
  ## did we succeed? stop-if-not:
  if(test.ans && qr.X$rank != qr(newX)$rank)
    stop(gettextf("determination of full column rank design matrix failed"),
         call. = FALSE)
  return(newX)
}


drop.cols <- function(mf, silent = FALSE) 
### drop columns from X and possibly NOM and S to ensure full column
### rank.
### mf - list with X and possibly NOM and S design matrices. Includes
### a ths object containing ths$alpha.names
### 
### returns: updated version of mf.
{
  nalpha <- length(mf$ths$alpha.names)
  ## X is assumed to contain an intercept at this point:
  Xint <- match("(Intercept)", colnames(mf$X), nomatch = 0)
    if(Xint <= 0) {
    mf$X <- cbind("(Intercept)" = rep(1, nrow(mf$X)), mf$X)
    warning("an intercept is needed and assumed")
  } ## intercept in X is guaranteed.
  if(!is.null(mf$NOM)){
    ## store coef names:
    mf$coef.names <- list()
    mf$coef.names$alpha <-
      paste(rep(mf$ths$alpha.names, ncol(mf$NOM)), ".",
            rep(colnames(mf$NOM), each=nalpha), sep="")
    mf$coef.names$beta <- colnames(mf$X)[-1]
    ## drop columns from NOM:
    mf$NOM <- drop.coef2(mf$NOM, silent=silent)
    ## drop columns from X:
    NOMX <- drop.coef2(cbind(mf$NOM, mf$X[,-1, drop=FALSE]),
                      silent=silent) 
    ## extract and store X:
    mf$X <- cbind("(Intercept)" = rep(1, nrow(mf$X)),
                  NOMX[,-seq_len(ncol(mf$NOM)), drop=FALSE])
    ## store alias information:
    mf$aliased <- list(alpha = rep(attr(mf$NOM, "aliased"),
                         each=nalpha)) 
    mf$aliased$beta <- attr(NOMX, "aliased")[-seq_len(ncol(mf$NOM))]
    if(!is.null(mf$S)) {
      mf$coef.names$zeta <- colnames(mf$S)[-1]
      ## drop columns from S:
      NOMS <- drop.coef2(cbind(mf$NOM, mf$S[,-1, drop=FALSE]),
                         silent=silent) 
      ## extract and store X:
      mf$S <- cbind("(Intercept)" = rep(1, nrow(mf$S)),
                    NOMS[,-seq_len(ncol(mf$NOM)), drop=FALSE])
      mf$aliased$zeta <- attr(NOMS, "aliased")[-seq_len(ncol(mf$NOM))] 
    }
    return(mf)
  }
  ## drop columns from X assuming an intercept:
  mf$coef.names <- list(alpha = mf$ths$alpha.names,
                        beta = colnames(mf$X)[-1])
  mf$X <- drop.coef2(mf$X, silent=silent)
  mf$aliased <- list(alpha = rep(0, nalpha),
                     beta = attr(mf$X, "aliased")[-1])
  ## drop columns from S if relevant:
  if(!is.null(mf$S)) { 
    Sint <- match("(Intercept)", colnames(mf$S), nomatch = 0)
    if(Sint <= 0) {
      mf$S <- cbind("(Intercept)" = rep(1, nrow(mf$S)), mf$S)
      warning("an intercept is needed and assumed in 'scale'",
              call.=FALSE)
    } ## intercept in S is guaranteed.
    mf$coef.names$zeta <- colnames(mf$S)[-1]
    mf$S <- drop.coef2(mf$S, silent=silent)
    mf$aliased$zeta <- attr(mf$S, "aliased")[-1]
  }
  return(mf)
}

eclm.finalize <- function(fit, weights, coef.names, aliased)
### destinguishing between par and coef where the former does not
### contain aliased coefficients.
{   
  nalpha <- length(aliased$alpha) 
  nbeta <- length(aliased$beta)
  nzeta <- length(aliased$zeta)
  ## nalias <- sum(unlist(alised))
  ncoef <- nalpha + nbeta + nzeta ## including aliased coef
  
  npar <- sum(!unlist(aliased)) ## excluding aliased coef
  stopifnot(length(fit$par) == npar)

  fit <- within(fit, {
    coefficients <- rep(NA, nalpha + nbeta + nzeta)
    ## insure correct order of alpha, beta and zeta:
    keep <- match(c("alpha", "beta", "zeta"), names(aliased),
                  nomatch=0) 
    aliased <- lapply(aliased[keep], as.logical)
    for(i in names(aliased))
      names(aliased[[i]]) <- coef.names[keep][[i]] 
    names(coefficients) <- unlist(coef.names[keep])
    par.names <- names(coefficients)[!unlist(aliased)]
    coefficients[!unlist(aliased)] <- fit$par
    alpha <- coefficients[1:nalpha]
    if(nbeta) beta <- coefficients[nalpha + 1:nbeta]
    if(nzeta) zeta <- coefficients[nalpha + nbeta + 1:nzeta]
    names(gradient) <- par.names ## names(coefficients)
    dimnames(Hessian) <- list(par.names, par.names)
    edf <- npar ## estimated degrees of freedom
    nobs <- sum(weights)
    n <- length(weights)
    fitted.values <- fitted
    df.residual = nobs - edf
    rm(list = c("par", "par.names", "keep", "i", "fitted"))
  })
  class(fit) <- "clm"
  return(fit)
}

setLinks <- function(rho, link) {
### The Aranda-Ordaz and log-gamma links are not supported in this
### version of clm.
  rho$pfun <- switch(link,
                     logit = plogis,
                     probit = pnorm,
                     cloglog = function(x, lower.tail=TRUE) pgumbel(x,
                                             lower.tail=lower.tail, max=FALSE),
                     cauchit = pcauchy,
                     loglog = pgumbel,
                     "Aranda-Ordaz" = function(x, lambda) pAO(x, lambda),
                     "log-gamma" = function(x, lambda) plgamma(x, lambda))
  rho$dfun <- switch(link,
                     logit = dlogis,
                     probit = dnorm,
                     cloglog = function(x) dgumbel(x, max=FALSE),
                     cauchit = dcauchy,
                     loglog = dgumbel,
                     "Aranda-Ordaz" = function(x, lambda) dAO(x, lambda),
                     "log-gamma" = function(x, lambda) dlgamma(x, lambda))
  rho$gfun <- switch(link,
                     logit = glogis,
                     probit = gnorm, 
                     cloglog = function(x) ggumbel(x, max=FALSE),
                     loglog = ggumbel,
                     cauchit = gcauchy,
                     "Aranda-Ordaz" = function(x, lambda) gAO(x, lambda), ## shouldn't happen
                     "log-gamma" = function(x, lambda) glgamma(x, lambda)
                     )
  rho$link <- link
}

makeThresholds <- function(y, threshold) {
### Generate the threshold structure summarized in the transpose of
### the Jacobian matrix, tJac. Also generating nalpha and alpha.names.

### args:
### y - response variable, a factor
### threshold - one of "flexible", "symmetric" or "equidistant"
  stopifnot(is.factor(y))
  lev <- levels(y)
  ntheta <- nlevels(y) - 1
  
  if(threshold == "flexible") {
    tJac <- diag(ntheta)
    nalpha <- ntheta
    alpha.names <- paste(lev[-length(lev)], lev[-1], sep="|")
  }
  
  if(threshold == "symmetric") {
    if(!ntheta >=2)
      stop("symmetric thresholds are only meaningful for responses with 3 or more levels", call.=FALSE)
    if(ntheta %% 2) { ## ntheta is odd
      nalpha <- (ntheta + 1)/2 ## No. threshold parameters
      tJac <- t(cbind(diag(-1, nalpha)[nalpha:1, 1:(nalpha-1)],
                      diag(nalpha)))
      tJac[,1] <- 1
      alpha.names <-
        c("central", paste("spacing.", 1:(nalpha-1), sep=""))
    }
    else { ## ntheta is even
      nalpha <- (ntheta + 2)/2
      tJac <- cbind(rep(1:0, each = ntheta / 2),
                    rbind(diag(-1, ntheta / 2)[(ntheta / 2):1,],
                          diag(ntheta / 2)))
      tJac[,2] <- rep(0:1, each = ntheta / 2)
      alpha.names <- c("central.1", "central.2",
                      paste("spacing.", 1:(nalpha-2), sep=""))
    }
  }
  
  if(threshold == "equidistant") {
    if(!ntheta >=2)
      stop("equidistant thresholds are only meaningful for responses with 3 or more levels", call.=FALSE)
    tJac <- cbind(1, 0:(ntheta-1))
    nalpha <- 2
    alpha.names <- c("threshold.1", "spacing")
  }
  return(list(tJac = tJac, nalpha = nalpha, alpha.names = alpha.names))
}

Trace <- function(iter, stepFactor, val, maxGrad, par, first=FALSE) {
    t1 <- sprintf(" %3d:  %-5e:    %.3f:   %1.3e:  ",
                  iter, stepFactor, val, maxGrad)
    t2 <- formatC(par)
    if(first)
        cat("iter:  step factor:     Value:     max|grad|:   Parameters:\n")
    cat(t1, t2, "\n")
}

##################################################################
## Functions for starting values:

start.threshold <-
  function(y, threshold = c("flexible", "symmetric", "equidistant")) 
### args:
### y - model response, a factor with at least two levels
### threshold - threshold structure, character.
{
  ## match and test arguments:
  threshold <- match.arg(threshold)
  stopifnot(is.factor(y) && nlevels(y) >= 2)
  ntheta <- nlevels(y) - 1L
  if(threshold %in% c("symmetric", "equidistant") && nlevels(y) < 3)
    stop(gettextf("symmetric and equidistant thresholds are only
meaningful for responses with 3 or more levels"))
  
  ## default starting values:
  start <- qlogis((1:ntheta) / (ntheta + 1) ) # just a guess
  
  ## adjusting for threshold functions:
  if(threshold == "symmetric" && ntheta %% 2) { ## ntheta odd >= 3
    nalpha <- (ntheta + 1) / 2
    start <- c(start[nalpha], diff(start[nalpha:ntheta])) ## works for
    ## ntheta >= 1
  }
  if(threshold == "symmetric" && !ntheta %% 2) {## ntheta even >= 4
    nalpha <- (ntheta + 2) / 2
    start <- c(start[c(nalpha - 1, nalpha)],
               diff(start[nalpha:ntheta])) ## works for ntheta >= 2
  }
  if(threshold == "equidistant")
    start <- c(start[1], mean(diff(start)))

  ## return starting values for the threshold parameters:
  return(as.vector(start))
}

start.beta <- function(X, has.intercept = TRUE)
  return(rep(0, NCOL(X) - has.intercept))

clm.start <- function(y, threshold, X, has.intercept = TRUE)
### could use eclm.start instead
  return(c(start.threshold(y, threshold),
           start.beta(X, has.intercept)))  

eclm.start <- function(y, threshold, X, NOM=NULL, S=NULL,
                       has.intercept=TRUE)
{
  st <- start.threshold(y, threshold)
  if(NCOL(NOM) > 1)
    st <- c(st, rep(rep(0, length(st)), ncol(NOM)-1))
  start <- c(st, start.beta(X, has.intercept))
  if(NCOL(S) > 1)
    start <- c(start, rep(0, ncol(S) - 1))
  start
}


clmm.start <- function(frames, link, threshold) {
  ## get starting values from clm:
  fit <- with(frames,
              clm.fit(y=y, X=X, weights=wts, offset=off, link=link,
                      threshold=threshold)) 
  
  ## initialize variance parameters to zero:
  start <- c(fit$par, rep(0, length(frames$grList)))
  return(start)
}

## set.start <-
##   function(rho, start=NULL, get.start=TRUE, threshold, link, frames) 
## {
##   ## set starting values for the parameters:
##   if(get.start) 
##     start <-
##       eclm.start(y=frames$y, threshold=threshold, X=frames$X,
##                  NOM=frames$NOM, S=frames$S, has.intercept=TRUE)
##   ## test start:
##   stopifnot(is.numeric(start))
##   length.start <- ncol(rho$B1) + NCOL(frames$S) - 1 #- length(rho$alised)
##   if(length(start) != length.start)
##     stop(gettextf("length of start is %d should equal %d",
##                   length(start), length.start), call.=FALSE)
##   ## start cauchit models at the probit estimates if start is not
##   ## supplied: 
##   if(link == "cauchit" && get.start) {
##     rho$par <- start
##     setLinks(rho, link="probit")
## ### Update this fit if class is eclm:
##     fit <- try(clm.fit.env(rho), silent=TRUE) ## standard control values
##     if(class(fit) == "try-error") 
##       stop("Failed to find suitable starting values: please supply some",
##            call.=FALSE)
##     start <- fit$par
##   }
##   return(start)
## }

set.start <-
  function(rho, start=NULL, get.start=TRUE, threshold, link, frames) 
{
  ## set starting values for the parameters:
  if(get.start) {
    start <- ## not 'starting' scale effects:
      eclm.start(y=frames$y, threshold=threshold, X=frames$X,
                 NOM=frames$NOM, has.intercept=TRUE)
    if(NCOL(frames$S) > 1 || link == "cauchit") {
### NOTE: only special start if NCOL(frames$S) > 1 (no reason for
### special start if scale is only offset and no predictors).
### NOTE: start cauchit models at the probit estimates if start is not
### supplied: 
      rho$par <- start
      if(link == "cauchit") setLinks(rho, link="probit")
      else setLinks(rho, link)
      tempk <- rho$k
      rho$k <- 0
      ## increased gradTol:
      fit <- try(clm.fit.env(rho, control=list(gradTol=1e-3)),
                 silent=TRUE) 
      if(class(fit) == "try-error") 
        stop("Failed to find suitable starting values: please supply some",
             call.=FALSE)
      start <- c(fit$par, rep(0, NCOL(frames$S) - 1))
      attr(start, "start.iter") <- fit$niter
      rho$k <- tempk
    }
  }
  ## test start:
  stopifnot(is.numeric(start))
  length.start <- ncol(rho$B1) + NCOL(frames$S) - 1 #- length(rho$alised)
  if(length(start) != length.start)
    stop(gettextf("length of start is %d should equal %d",
                  length(start), length.start), call.=FALSE)

  return(start)
}

##################################################################

response.name <- function(terms) {
  vars <- as.character(attr(terms, "variables"))
  vars[1 + attr(terms, "response")]
}

getB <- function(y, NOM=NULL, X=NULL, offset=NULL, tJac=NULL) {
### NOTE: no tests that arguments conform.
  nlev <- nlevels(y)
  n <- length(y)
  B2 <- 1 * (col(matrix(0, n, nlev)) == c(unclass(y)))
  o1 <- c(1e5 * B2[, nlev]) - offset
  o2 <- c(-1e5 * B2[,1]) - offset
  B1 <- B2[, -(nlev), drop = FALSE]
  B2 <- B2[, -1, drop = FALSE]
  ## adjust B1 and B2 for structured thresholds:
  if(!is.null(tJac)) {
    B1 <- B1 %*% tJac
    B2 <- B2 %*% tJac
  }
  ## update B1 and B2 with nominal effects:
  if(NCOL(NOM) > 1) { ## !is.null(NOM) && ncol(NOM) > 1) {
    ## if !is.null(NOM) and NOM is more than an intercept:
    LL1 <- lapply(1:ncol(NOM), function(x) B1 * NOM[,x])
    B1 <- do.call(cbind, LL1)
    LL2 <- lapply(1:ncol(NOM), function(x) B2 * NOM[,x])
    B2 <- do.call(cbind, LL2)
  }
  ## update B1 and B2 with location effects (X):
  nbeta <- NCOL(X) - 1
  if(NCOL(X) > 1) {
    B1 <- cbind(B1, -X[, -1, drop = FALSE])
    B2 <- cbind(B2, -X[, -1, drop = FALSE])
  }
  dimnames(B1) <- NULL
  dimnames(B2) <- NULL
  list(B1=B1, B2=B2, o1=o1, o2=o2) 
}

Deparse <-
  function(expr, width.cutoff = 500L, backtick = mode(expr) %in%  
           c("call", "expression", "(", "function"),
           control = c("keepInteger", "showAttributes", "keepNA"),
           nlines = -1L)
  deparse(expr=expr, width.cutoff= width.cutoff, backtick=backtick,
          control=control, nlines=nlines) 
