## This file contains:
## The function clm.fit() - an lm.fit or glm.fit equivalent for CLMs.

clm.fit <-
  function(y, X, S, N, weights = rep(1, nrow(X)),
           offset = rep(0, nrow(X)), S.offset = rep(0, nrow(X)),
           control = list(), start,
           link = c("logit", "probit", "cloglog", "loglog", "cauchit"),
           threshold = c("flexible", "symmetric", "symmetric2", "equidistant"))
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
  stopifnot(is.factor(y), is.matrix(X))
  if(missing(weights) || is.null(weights))
      weights <- rep(1, length(y))
  if(missing(offset) || is.null(offset))
      offset <- rep(0, length(y))
  if(missing(S.offset) || is.null(S.offset))
      S.offset <- rep(0, length(y))
  stopifnot(length(y) == nrow(X) &&
            length(y) == length(weights) &&
            length(y) == length(offset) &&
            length(y) == length(S.offset))
  frames <- list(y=y, X=X)
  y[weights <= 0] <- NA
  frames$ylevels <- levels(droplevels(y))
  ## S and N are optional:
  if(!missing(S) && !is.null(S)) {
    frames$S <- S
    stopifnot(is.matrix(S),
              length(y) == nrow(S))
  }
  if(!missing(N) && !is.null(N)) {
    frames$NOM <- N
    stopifnot(is.matrix(N),
              length(y) == nrow(N))
  }

  ## Get  threshold structure:
  frames$ths <- makeThresholds(frames$ylevels, threshold)
  ## Test for column rank deficiency in design matrices:
  frames <- drop.cols(frames, silent=TRUE)

  ## Set envir rho with variables: B1, B2, o1, o2, wts, fitted...:
  rho <- clm.newRho(parent.frame(), y=frames$y, X=frames$X,
                    NOM=frames$NOM, S=frames$S, weights=weights,
                    offset=offset, S.offset=S.offset,
                    tJac=frames$ths$tJac)

  ## Set starting values for the parameters:
  start <- set.start(rho, start=start, get.start=missing(start),
                     threshold=threshold, link=link, frames=frames)
  rho$par <- as.vector(start) ## remove attributes

  ## Set inverse link function and its derivatives (pfun, dfun and
  ## gfun) in rho:
  setLinks(rho, link)

  ## Fit the model:
  fit <- if(control$method == "Newton") {
      clm.fit.NR(rho, control) } else {
          clm.fit.optim(rho, control$method, control$ctrl) }

  ## Format and return the fit:
  fit$control <- control
  fit$coef.names <- frames$coef.names
  fit$aliased <- lapply(frames$aliased, as.logical)
  if(control$method == "Newton" &&
     !is.null(start.iter <- attr(start, "start.iter")))
    fit$niter <- fit$niter + start.iter

  ## Check convergence:
  conv <- conv.check(fit, Theta.ok=TRUE, tol=control$tol)
### NOTE: we are not checking if the thresholds are increasing (if
### there are nominal effects) even though potentially we could.
  print.conv.check(conv, action=control$convergence) ## print convergence message
  fit$vcov <- conv$vcov
  fit$condHess <- conv$cond.H
  fit$convergence <- conv[!names(conv) %in% c("vcov", "cond.H")]

  return(fit)
### FIXME: should 'par' be 'coefficients' to allow coef(fit) etc.?
### FIXME: should fit contain 'vcov'?
}
