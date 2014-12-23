## This file contains:
## The main clm function and some auxiliary functions to generate the
## model frame and handle the model environment.

clm <-
  function(formula, scale, nominal, data, weights, start, subset,
           doFit = TRUE, na.action, contrasts, model = TRUE,
           control = list(),
           link = c("logit", "probit", "cloglog", "loglog", "cauchit"),
           threshold = c("flexible", "symmetric", "symmetric2", "equidistant"), ...)
### deliberately no offset argument - include offset in the relevant
### formula/scale instead.
###
### FIXME: drop the tJac argument and let threshold accept it as a
### numeric matrix. Also test that ncol(tJac) <= nrow(tJac) or,
### perhaps better: rank(tJac) <= nlevels(y).
###
### FIXME: allow "threshold="fixed", theta=NULL" arguments to make it
### possible to fit certain kinds of models.
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

  ## Extract and process formulas:
  formulas <- get_clmFormulas(mc)
  ## Get full model.frame:
  fullmf <- get_clmFullmf(mc, fullForm=formulas$fullForm,
                          form.envir=formulas$form.envir)
  ## Compute: y, X, wts, off, fullmf:
  frames <- get_clmDesign(fullmf, formulas, contrasts)
  ## frames <- clm.model.frame(mc, contrasts)

  ## Compute the transpose of the Jacobian for the threshold function,
  ## tJac and the names of the threshold parameters, alpha.names:
  frames$ths <- makeThresholds(frames$ylevels, threshold)

  ## Return model.frame?
  if(control$method == "model.frame") {
      frames$mf <- fullmf
      return(frames)
  }

  ## Test column rank deficiency and possibly drop some
  ## parameters. Also set the lists 'aliased' and 'coef.names':
  ## (X is with intercept at this point.)
  frames <- drop.cols(frames, drop.scale=FALSE, silent=TRUE)
### Note: intercept could be dropped from X and S in drop.cols?
### Currently they are dropped in clm.newRho instead.

  ## Set envir rho with variables: B1, B2, o1, o2, wts, fitted:
  rho <- with(frames, {
      clm.newRho(parent.frame(), y=y, X=X, NOM=frames$NOM, S=frames$S,
                  weights=wts, offset=off, S.offset=frames$S.off,
                  tJac=ths$tJac)
  })

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
      fit <- clm.fit.NR(rho, control)
  } else
      fit <- clm.fit.optim(rho, control$method, control$ctrl)

  ## Format fit information:
  res <- clm.finalize(fit, weights=frames$wts,
                      coef.names=frames$coef.names,
                      aliased=frames$aliased)
  ## Save stuff in model object:
  lst <- list(formulas=formulas, formula=formulas$formula, link=link,
              start=start, control=control, threshold=threshold,
              call=match.call())
  res <- c(res, lst, extractFromFrames(frames, fullmf))
### NOTE: need formula in addition to formulas to aviod formulas being
### picked by accident (e.g. by update, drop1, etc.) with partial
### matching of formula.

  ## Add iterations to get starting values to niter count:
  if(control$method == "Newton" &&
     !is.null(start.iter <- attr(start, "start.iter")))
      res$niter <- res$niter + start.iter

  ## Add Theta table to results list:
  Theta.ok <- TRUE ## convergence condition
  if(!is.null(frames$NOM)) {
      tmp <- formatTheta(res, frames$NOM, fullmf)
      Theta.ok <- tmp$Theta.ok
      res <- c(res, tmp[!names(tmp) %in% "Theta.ok"])
  } else {
      res$Theta <- res$alpha %*% t(res$tJac)
      ## colnames(res$Theta) <- names(res$alpha)
  }
### NOTE: we need to check convergence *after* generation of the Theta
### matrix, since here we check if thresholds are increasing.

  ## Check convergence:
  conv <- conv.check(res, Theta.ok=Theta.ok, tol=control$tol)
  print.conv.check(conv, action=control$convergence) ## print convergence message

  ## Add more stuff to the results list:
  res <- c(res, list(vcov=conv$vcov, condHess=conv$cond.H))
  res$convergence=conv[!names(conv) %in% c("vcov", "cond.H")]
  res$info <- get_clmInfoTab(res)
  if(model) res$model <- fullmf
  res$par <- res$fitted <- res$niter <- NULL

  ## order elements of result alphabetically:
  res <- res[order(tolower(names(res)))]
  class(res) <- "clm"
  res
}

## clm.model.frame <- function(mc, contrasts) {
## ### mc - the matched call
## ### contrasts - contrasts for the model terms
##
##     ## Collect all variables in a full formula:
##   ## evaluate the formulae in the enviroment in which clm was called
##   ## (parent.frame(2)) to get them evaluated properly:
##   forms <- list(eval.parent(mc$formula, 2))
##   if(!is.null(mc$scale)) forms$scale <- eval.parent(mc$scale, 2)
##   if(!is.null(mc$nominal)) forms$nominal <- eval.parent(mc$nominal, 2)
##   ## get the environment of the formula. If this does not have an
##   ## enviroment (it could be a character), then use the parent frame.
##   form.envir <-
##     if(!is.null(env <- environment(forms[[1]]))) env
##     else parent.frame(2)
##   ## ensure formula, scale and nominal are formulas:
##   ## forms <- lapply(forms, function(x) {
##   ##   try(formula(deparse(x), env = form.envir), silent=TRUE) })
##   for(i in 1:length(forms)) {
##     forms[[i]] <- try(formula(deparse(forms[[i]]),
##                               env = form.envir), silent=TRUE)
##   }
##   if(any(sapply(forms, function(f) class(f) == "try-error")))
##     stop("unable to interpret 'formula', 'scale' or 'nominal'")
##   ## collect all variables in a full formula:
##   fullForm <- do.call("getFullForm", forms)
##   ## set environment of 'fullForm' to the environment of 'formula':
##   environment(fullForm) <- form.envir
##
##   ## Extract the full model.frame(mf):
##   m <- match(c("data", "subset", "weights", "na.action"),
##              names(mc), 0)
##   mf <- mc[c(1, m)]
##   mf$formula <- fullForm
##   mf$drop.unused.levels <- TRUE
##   mf[[1]] <- as.name("model.frame")
##   if(is.null(mf$data)) mf$data <- form.envir
##   fullmf <- eval(mf, envir = parent.frame(2))
##   mf$na.action <- "na.pass" ## filter NAs by hand below
##   ## Check specification of contrasts:
##   checkContrasts(terms=attr(fullmf, "terms"), contrasts=contrasts)
##
##   ## Extract X:
##   ## get X from fullmf to handle e.g., NAs and subset correctly
##   mf$formula <- forms[[1]]
##   X.mf <- eval(mf, envir = parent.frame(2))
##   X.terms <- attr(X.mf, "terms")
##   X <- model.matrix(X.terms, data=fullmf,
##                     contrasts.arg=getContrasts(X.terms, contrasts))
##   n <- nrow(X)
##   ## Test for intercept in X:
##   Xint <- match("(Intercept)", colnames(X), nomatch = 0)
##   if(Xint <= 0) {
##     X <- cbind("(Intercept)" = rep(1, n), X)
##     warning("an intercept is needed and assumed in 'formula'",
##             call.=FALSE)
##   } ## intercept in X is guaranteed.
##   off <- getOffset(X.mf)
##   ## Filter NAs if any:
##   if(!is.null(naa <- na.action(fullmf)))
##     off <- off[-naa]
##
##   wts <- getWeights(fullmf)
##   ## Extract model response:
##   y <- model.response(fullmf, "any") ## any storage mode
##   if(!is.factor(y)) stop("response needs to be a factor", call.=FALSE)
##   ## ylevels are the levels of y with positive weights
##   ylevels <- levels(droplevels(y[wts > 0]))
##   ## check that y has at least two levels:
##   if(length(ylevels) == 1L)
##       stop(gettextf("response has only 1 level ('%s'); expecting at least two levels",
##                     ylevels), call.=FALSE)
##   if(!length(ylevels))
##       stop("response should be a factor with at least two levels")
##   ## stop("response factor should have at least two levels")
##
##   ## list of results:
##   res <- list(y=y, ylevels=ylevels, X=X, wts=wts, off=off,
##               mf=fullmf, terms=X.terms)
##
##   ## Extract S (design matrix for the scale effects):
##   if(!is.null(mc$scale)) {
##     mf$formula <- forms$scale
##     S.mf <- eval(mf, envir = parent.frame(2))
##     if(!is.null(model.response(S.mf)))
##       stop("response not allowed in 'scale'", call.=FALSE)
##     res$S.terms <- attr(S.mf, "terms")
##     S <- model.matrix(res$S.terms, data=fullmf,
##                       contrasts.arg=getContrasts(res$S.terms, contrasts))
##     ## Test for intercept in S:
##     Sint <- match("(Intercept)", colnames(S), nomatch = 0)
##     if(Sint <= 0) {
##       S <- cbind("(Intercept)" = rep(1, n), S)
##       warning("an intercept is needed and assumed in 'scale'",
##               call.=FALSE)
##     } ## intercept in S is guaranteed.
##     res$S <- S
##     res$S.off <- getOffset(S.mf)
##     ## Filter NAs if any:
##     if(!is.null(naa <- na.action(fullmf)))
##       res$S.off <- res$S.off[-naa]
##   }
##
##   ## Extract NOM (design matrix for the nominal effects):
##   if(!is.null(mc$nominal)) {
##     mf$formula <- forms$nominal
##     nom.mf <- eval(mf, envir = parent.frame(2))
##     if(!is.null(model.response(nom.mf)))
##       stop("response not allowed in 'nominal'", call.=FALSE)
##     if(!is.null(model.offset(nom.mf)))
##       stop("offset not allowed in 'nominal'", call.=FALSE)
##     res$nom.terms <- attr(nom.mf, "terms")
##     NOM <- model.matrix(res$nom.terms, data=fullmf,
##                         contrasts.arg=getContrasts(res$nom.terms, contrasts))
##     NOMint <- match("(Intercept)", colnames(NOM), nomatch = 0)
##     if(NOMint <= 0) {
##       NOM <- cbind("(Intercept)" = rep(1, n), NOM)
##       warning("an intercept is needed and assumed in 'nominal'",
##               call.=FALSE)
##     } ## intercept in NOM is guarantied.
##     res$NOM <- NOM
##   }
##
##   ## return results:
##   return(res)
##   ## Note: X, S and NOM are with dimnames and intercepts are
##   ## guaranteed. They may be column rank deficient.
## }

## clm.newRho <-
##   function(parent, y, X, weights, offset, tJac)
## ### Set variables in rho: B1, B2, o1, o2 and wts.
## {
##   rho <- new.env(parent = parent)
##
##   ## Make B1, B2, o1, o2 based on y, X and tJac:
##   ntheta <- nlevels(y) - 1
##   n <- nrow(X)
##   B2 <- 1 * (col(matrix(0, n, ntheta + 1)) == c(unclass(y)))
##   rho$o1 <- c(1e5 * B2[, ntheta + 1]) - offset
##   rho$o2 <- c(-1e5 * B2[,1]) - offset
##   B1 <- B2[, -(ntheta + 1), drop = FALSE]
##   B2 <- B2[, -1, drop = FALSE]
##   ## adjust B1 and B2 for structured thresholds:
##   rho$B1 <- B1 %*% tJac
##   rho$B2 <- B2 %*% tJac
##   ## update B1 and B2 with location effects (X):
##   nbeta <- NCOL(X) - 1
##   if(nbeta > 0) {
##     rho$B1 <- cbind(rho$B1, -X[, -1, drop = FALSE])
##     rho$B2 <- cbind(rho$B2, -X[, -1, drop = FALSE])
##   }
##   dimnames(rho$B1) <- NULL
##   dimnames(rho$B2) <- NULL
##
##   rho$fitted <- numeric(length = n)
##   rho$wts <- weights
##
##   return(rho)
## }

clm.newRho <-
  function(parent=parent.frame(), y, X, NOM=NULL, S=NULL, weights,
           offset, S.offset=NULL, tJac)
### Setting variables in rho: B1, B2, o1, o2, wts.
{
  rho <- new.env(parent = parent)
  ## Make B1, B2, o1, o2 based on y, X and tJac:
  ## rho <- list2env(getB(y=y, NOM=NOM, X=X, offset=offset,
  ## tJac=tJac), parent=parent)
  keep <- weights > 0
  y[!keep] <- NA
  y <- droplevels(y)
  ntheta <- nlevels(y) - 1
  y <- c(unclass(y))
  y[is.na(y)] <- 0
  n <- sum(keep)
  B2 <- 1 * (col(matrix(0, nrow(X), ntheta + 1)) == y)
  rho$o1 <- c(1e5 * B2[keep, ntheta + 1]) - offset[keep]
  rho$o2 <- c(-1e5 * B2[keep, 1]) - offset[keep]
  B1 <- B2[keep, -(ntheta + 1), drop = FALSE]
  B2 <- B2[keep, -1, drop = FALSE]
  ## adjust B1 and B2 for structured thresholds:
  rho$B1 <- B1 %*% tJac
  rho$B2 <- B2 %*% tJac
  ## update B1 and B2 with nominal effects:
  if(NCOL(NOM) > 1) { ## !is.null(NOM) && ncol(NOM) > 1) {
    ## if !is.null(NOM) and NOM is more than an intercept:
    LL1 <- lapply(1:ncol(NOM), function(x) rho$B1 * NOM[keep, x])
    rho$B1 <- do.call(cbind, LL1)
    LL2 <- lapply(1:ncol(NOM), function(x) rho$B2 * NOM[keep, x])
    rho$B2 <- do.call(cbind, LL2)
  }
  ## update B1 and B2 with location effects (X):
  nbeta <- NCOL(X) - 1
  if(nbeta > 0) {
    rho$B1 <- cbind(rho$B1, -X[keep, -1, drop = FALSE])
    rho$B2 <- cbind(rho$B2, -X[keep, -1, drop = FALSE])
  }
  dimnames(rho$B1) <- NULL
  dimnames(rho$B2) <- NULL
  rho$n.psi <- ncol(rho$B1) ## no. linear model parameters
  rho$k <- 0
  ## there may be scale offset without scale predictors:
  rho$sigma <- rho$Soff <-
    if(is.null(S.offset)) rep(1, n) else exp(S.offset[keep])
  ## save scale model:
  if(!is.null(S)) {
    rho$S <- S[keep, -1, drop=FALSE]
    dimnames(rho$S) <- NULL
    rho$k <- ncol(rho$S) ## no. scale parameters
  }
  rho$has.scale <- ## TRUE if scale has to be considered.
    (!is.null(S) || any(S.offset != 0))
  ## initialize fitted values and weights:
  rho$fitted <- numeric(length = n)
  rho$wts <- weights[keep]
  ## Setting likelihood, gradient and Hessian functions:
  rho$clm.nll <- clm.nll
  rho$clm.grad <- clm.grad
  rho$clm.hess <- clm.hess
  ## return:
  return(rho)
}

clm.finalize <- function(fit, weights, coef.names, aliased)
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
    rm(list = c("par.names", "keep", "i"))
  })
  class(fit) <- "clm"
  return(fit)
}
