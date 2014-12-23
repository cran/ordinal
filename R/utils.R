## This file contains:
## Various utility functions.

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

makeThresholds <- function(ylevels, threshold) { ## , tJac) {
### Generate the threshold structure summarized in the transpose of
### the Jacobian matrix, tJac. Also generating nalpha and alpha.names.

### args:
### y - response variable, a factor
### threshold - one of "flexible", "symmetric" or "equidistant"
  ## stopifnot(is.factor(y))
  lev <- ylevels
  ntheta <- length(lev) - 1

  ## if(!is.null(tJac)) {
  ##   stopifnot(nrow(tJac) == ntheta)
  ##   nalpha <- ncol(tJac)
  ##   alpha.names <- colnames(tJac)
  ##   if(is.null(alpha.names) || anyDuplicated(alpha.names))
  ##     alpha.names <- as.character(1:nalpha)
  ##   dimnames(tJac) <- NULL
  ## }
  ## else { ## threshold structure identified by threshold argument:
    if(threshold == "flexible") {
      tJac <- diag(ntheta)
      nalpha <- ntheta
      alpha.names <- paste(lev[-length(lev)], lev[-1], sep="|")
    }

    if(threshold == "symmetric") {
      if(!ntheta >=2)
        stop("symmetric thresholds are only meaningful for responses with 3 or more levels",
             call.=FALSE)
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
    ## Assumes latent mean is zero:
    if(threshold == "symmetric2") {
      if(!ntheta >=2)
        stop("symmetric thresholds are only meaningful for responses with 3 or more levels",
             call.=FALSE)
      if(ntheta %% 2) { ## ntheta is odd
        nalpha <- (ntheta - 1)/2 ## No. threshold parameters
        tJac <- rbind(apply(-diag(nalpha), 1, rev),
                      rep(0, nalpha),
                      diag(nalpha))
      }
      else { ## ntheta is even
        nalpha <- ntheta/2
        tJac <- rbind(apply(-diag(nalpha), 1, rev),
                      diag(nalpha))
      }
      alpha.names <- paste("spacing.", 1:nalpha, sep="")
    }

    if(threshold == "equidistant") {
      if(!ntheta >=2)
        stop("equidistant thresholds are only meaningful for responses with 3 or more levels",
             call.=FALSE)
      tJac <- cbind(1, 0:(ntheta-1))
      nalpha <- 2
      alpha.names <- c("threshold.1", "spacing")
    }
  ## }
  return(list(tJac = tJac, nalpha = nalpha, alpha.names = alpha.names))
}

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
    ## if (any(wts <= 0))
    ##   stop(gettextf("non-positive weights are not allowed"),
    ##        call.=FALSE)
### NOTE: We do not remove observations where weights == 0, because
### that could be a somewhat surprising behaviour. It would also
### require that the model.frame be evaluated all over again to get
### the right response vector with the right number of levels.
    if(length(wts) && length(wts) != n)
        stop(gettextf("number of weights is %d should equal %d (number of observations)",
                      length(wts), n), call.=FALSE)
    if(any(wts < 0))
        stop(gettextf("negative weights are not allowed"),
             call.=FALSE)
    ## if(any(wts == 0)) {
    ##     y <- model.response(mf, "any")
    ##     if(any(table(y[wts > 0]) == 0))
    ##         stop(gettextf("zero positive weights for one or more response categories"),
    ##              call.=FALSE)
    ## }
    return(as.double(wts))
}

getOffset <- function(mf, terms) {
### mf - model.frame
    n <- nrow(mf)
    off <- rep(0, n)
    if(!is.null(o <- attr(terms, "offset"))) {
        if(length(o) > 1)
            stop("only one offset term allowed in each formula", call.=FALSE)
        varnm <- attr(terms, "variables")
        ## deparse all variable names - character vector:
        varnm <- unlist(lapply(as.list(varnm), deparse)[-1])
        off <- mf[, varnm[o]]
    }
    ## off <- as.vector(mf[, o])
    if(length(off) && length(off) != n)
        stop(gettextf("number of offsets is %d should equal %d (number of observations)",
                      length(off), n), call.=FALSE)
    return(as.double(off))
}

getOffsetStd <- function(mf) {
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

## getCtrlArgs <- function(control, extras) {
## ### Recover control arguments from clmm.control and extras (...):
## ###
##   ## Collect control arguments in list:
##   ctrl.args <- c(extras, control$method, control$useMatrix,
##                  control$ctrl, control$optCtrl)
##   ## Identify the two occurences "trace", delete them, and add trace=1
##   ## or trace=-1 to the list of arguments:
##   which.trace <- which(names(ctrl.args) == "trace")
##   trace.sum <- sum(unlist(ctrl.args[which.trace]))
##   ctrl.args <- ctrl.args[-which.trace]
##   ## remove duplicated arguments:
##   ctrl.args <- ctrl.args[!duplicated(names(ctrl.args))]
##   if(trace.sum >= 1) ctrl.args$trace <- 1
##   if(trace.sum >= 2 || trace.sum <= -1) ctrl.args$trace <- -1
##   ## return the updated list of control parameters:
##   do.call("clmm.control", ctrl.args)
## }

getCtrlArgs <- function(control, extras) {
### Recover control arguments from clmm.control and extras (...):
###
    if(!is.list(control))
        stop("'control' should be a list")
    ## Collect control arguments in list:
    ## 1) assuming 'control' is a call to clmm.control:
        ctrl.args <-
        if(setequal(names(control), names(clmm.control())))
            c(extras, control["method"], control["useMatrix"],
              control$ctrl, control$optCtrl)
    ## assuming 'control' is specified with control=list( 'args'):
        else
            c(extras, control)
### NOTE: having c(extras, control) rather than c(control, extras)
### means that extras have precedence over control.
    ## Identify the two occurences "trace", delete them, and add trace=1
    ## or trace=-1 to the list of arguments:
    which.trace <- which(names(ctrl.args) == "trace")
    trace.sum <- sum(unlist(ctrl.args[which.trace]))
    if(trace.sum)
        ctrl.args <- ctrl.args[-which.trace]
    ## remove duplicated arguments:
    ctrl.args <- ctrl.args[!duplicated(names(ctrl.args))]
    if(trace.sum >= 1) ctrl.args$trace <- 1
    if(trace.sum >= 2 || trace.sum <= -1) ctrl.args$trace <- -1
    ## return the updated list of control parameters:
    do.call("clmm.control", ctrl.args)
}

Trace <- function(iter, stepFactor, val, maxGrad, par, first=FALSE) {
    t1 <- sprintf(" %3d:  %-5e:    %.3f:   %1.3e:  ",
                  iter, stepFactor, val, maxGrad)
    t2 <- formatC(par)
    if(first)
        cat("iter:  step factor:     Value:     max|grad|:   Parameters:\n")
    cat(t1, t2, "\n")
}

response.name <- function(terms) {
  vars <- as.character(attr(terms, "variables"))
  vars[1 + attr(terms, "response")]
}

getB <- function(y, NOM=NULL, X=NULL, offset=NULL, tJac=NULL) {
### FIXME: Is this function ever used?
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
### FIXME: test if formals(Deparse) == formals(deparse)??
  deparse(expr=expr, width.cutoff= width.cutoff, backtick=backtick,
          control=control, nlines=nlines)

getContrasts <- function(terms, contrasts) {
    if(is.null(contrasts)) return(NULL)
    term.labels <- attr(terms, "term.labels")
    contrasts[names(contrasts) %in% term.labels]
}

checkContrasts <- function(terms, contrasts) {
### Check that contrasts are not specified for absent factors and warn
### if about them
    term.labels <- attr(terms, "term.labels")
    nm.contr <- names(contrasts)
    notkeep <- nm.contr[!nm.contr %in% term.labels]
    msg <-
        if(length(notkeep) > 2)
            "variables '%s' are absent: their contrasts will be ignored"
        else "variable '%s' is absent: its contrasts will be ignored"
    if(length(notkeep))
        warning(gettextf(msg, paste(notkeep, collapse=", ")),
                call.=FALSE)
    invisible()
}

get_clmInfoTab <- function(object, ...) {
    names <- c("link", "threshold", "nobs", "logLik", "edf", "niter",
               "maxGradient", "condHess")
    stopifnot(all(names %in% names(object)))
    info <- with(object, {
        data.frame("link" = link,
                   "threshold" = threshold,
                   "nobs" = nobs,
                   "logLik" = formatC(logLik, digits=2, format="f"),
                   "AIC" = formatC(-2*logLik + 2*edf, digits=2,
                   format="f"),
                   "niter" = paste(niter[1], "(", niter[2], ")", sep=""),
### NOTE: iterations to get starting values for scale models *are*
### included here.
                   "max.grad" = formatC(maxGradient, digits=2,
                   format="e"),
                   "cond.H" = formatC(condHess, digits=1, format="e")
                   ## BIC is not part of output since it is not clear what
                   ## the no. observations are.
                   )
    })
    info
}

format_tJac <- function(frames) {
    lev <- frames$ylevels
    tJac <- frames$ths$tJac
    rownames(tJac) <- paste(lev[-length(lev)], lev[-1], sep="|")
    colnames(tJac) <- frames$ths$alpha.names
    tJac
}

extractFromFrames <- function(frames, fullmf) {
    lst <- list(y.levels=frames$ylevels,
                na.action=attr(fullmf, "na.action"),
                tJac=format_tJac(frames))

    lstX <- list(contrasts=attr(frames$X, "contrasts"), terms=frames$terms,
                 xlevels=.getXlevels(frames$terms, fullmf))
    lst <- c(lst, lstX)

    if(!is.null(frames$S))
        lst <- c(lst, list(S.contrasts=attr(frames$S, "contrasts"),
                           S.terms=frames$S.terms,
                           xlevels=.getXlevels(frames$S.terms, fullmf)))
    if(!is.null(frames$NOM))
        lst <- c(lst, list(nom.contrasts=attr(frames$NOM, "contrasts"),
                           nom.terms=frames$nom.terms,
                           nom.xlevels=.getXlevels(frames$nom.terms, fullmf)))
    lst
}

formatTheta <- function(object, NOM, fullmf) {
    res <- object
    Theta.ok <- TRUE

    ## Get matrix of thresholds; Theta:
    Theta.list <-
        getThetamat(terms=res$nom.terms,
                    alpha=res$alpha,
                    assign=attr(NOM, "assign"),
                    contrasts=res$nom.contrasts,
                    tJac=res$tJac,
                    xlevels=res$nom.xlevels,
                    dataClasses=get_dataClasses(fullmf))
    ## Test that thresholds are increasing:
        if(all(is.finite(unlist(Theta.list$Theta)))) {
        th.increasing <- apply(Theta.list$Theta, 1, function(th)
                               all(diff(th) >= 0))
        if(!all(th.increasing))
            Theta.ok <- FALSE
    }
    Theta <- if(length(Theta.list) == 2)
        with(Theta.list, cbind(mf.basic, Theta)) else Theta.list$Theta
    alpha.mat <- matrix(res$alpha, ncol=ncol(res$tJac), byrow=TRUE)
    colnames(alpha.mat) <- colnames(res$tJac)
    rownames(alpha.mat) <- attr(NOM, "orig.colnames")
    ## Return
    list(Theta=Theta, alpha.mat=alpha.mat, Theta.ok=Theta.ok)
}

get_dataClasses <- function(mf) {
    if(!is.null(Terms <- attr(mf, "terms")) &&
       !is.null(dataCl <- attr(Terms, "dataClasses")))
        return(dataCl)
    sapply(mf, .MFclass)
}
