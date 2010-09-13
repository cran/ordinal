predict.clm <-
  function(object, newdata, se.fit = FALSE, interval = FALSE,
           level = 0.95, 
           type = c("prob", "class", "cum.prob", "linear.predictor"),
           na.action = na.pass, ...) 
### result - a list of predictions (fit)
### FIXME: restore names of the fitted values
###
### Assumes object has terms, xlevels, contrasts, tJac
{
  ## match and test arguments:
  type <- match.arg(type)
  se.fit <- as.logical(se.fit)[1]
  interval <- as.logical(interval)[1]
  stopifnot(length(level) == 1 && is.numeric(level) && level < 1 &&
            level > 0)
  if(type == "class" && (se.fit || interval)) {
    warning("se.fit and interval set to FALSE for type = 'class'")
    se.fit <- interval <- FALSE
  }
  cov <- if(se.fit || interval) vcov(object) else NULL
### Get newdata object; fill in response if missing and always for
### type=="class":
  has.response <- TRUE
  if(type == "class" && missing(newdata)) 
    ## newdata <- update(object, method="model.frame")$mf
    newdata <- model.frame(object)
  ## newdata supplied or type=="class":
  if(!(missing(newdata) || is.null(newdata)) || type=="class") {
    if(!(missing(newdata) || is.null(newdata)) &&
       sum(unlist(object$aliased)) > 0)
      warning("predictions from column rank-deficient fit may be misleading") 
    newdata <- as.data.frame(newdata)
    ## Test if response is in newdata:
    resp <- response.name(object$terms)
    ## remove response from newdata if type == "class"
    if(type == "class") newdata <- newdata[!names(newdata) %in% resp]
    has.response <- resp %in% names(newdata) ##  FALSE for type == "class"
    if(!has.response) { 
      ## fill in response variable in newdata if missing:
      ylev <- object$y.levels
      nlev <- length(ylev)
      nnd <- nrow(newdata)
      newdata <-
        cbind(newdata[rep(1:nnd, each=nlev) , , drop=FALSE],
              factor(rep(ylev, nnd), levels=ylev, ordered=TRUE))
      names(newdata)[ncol(newdata)] <- resp
    }
### Set model matrices:
    mf <- model.frame(object$terms, newdata, na.action=na.action,
                      xlev=object$xlevels)
    ## model.frame will warn, but here we also throw an error:
    if(nrow(mf) != nrow(newdata))
      stop("length of variable(s) found do not match nrow(newdata)")
    ## check that variables are of the right type:
    if (!is.null(cl <- attr(object$terms, "dataClasses")))
      .checkMFClasses(cl, mf)
    ## make model.matrix:
    X <- model.matrix(object$terms, mf, contrasts = object$contrasts)
    Xint <- match("(Intercept)", colnames(X), nomatch = 0L)
    n <- nrow(X)
    if(Xint <= 0) X <- cbind("(Intercept)" = rep(1, n), X)
    ## drop aliased columns:
    if(sum(object$aliased$beta) > 0)
      X <- X[, !c(FALSE, object$aliased$beta), drop=FALSE]
    ## handle offset (from predict.lm):
    offset <- rep(0, nrow(X))
    if(!is.null(off.num <- attr(object$terms, "offset"))) 
      for(i in off.num) offset <- offset +
        eval(attr(object$terms, "variables")[[i + 1]], newdata)
    y <- model.response(mf)
### make NOMINAL model.matrix:
    if(is.nom <- !is.null(object$nom.terms)) {
      ## allows NAs to pass through to fit, se.fit, lwr and upr:
      nom.mf <- model.frame(object$nom.terms, newdata,
                            na.action=na.pass,
                            xlev=object$nom.xlevels)  
      ## model.frame will warn, but here we also throw an error:
      if(nrow(nom.mf) != nrow(newdata))
        stop("length of variable(s) found do not match nrow(newdata)")
      if (!is.null(cl <- attr(object$nom.terms, "dataClasses")))
        .checkMFClasses(cl, nom.mf)
      NOM <- model.matrix(object$nom.terms, nom.mf,
                          contrasts=object$nom.contrasts)
      NOMint <- match("(Intercept)", colnames(NOM), nomatch = 0L)
      if(NOMint <= 0) NOM <- cbind("(Intercept)" = rep(1, n), NOM)
      alias <- t(matrix(object$aliased$alpha,
                        nrow=length(object$y.levels) - 1))[,1] 
      if(sum(alias) > 0)
        NOM <- NOM[, !c(FALSE, alias), drop=FALSE]
    }
### make SCALE model.matrix:
    if(is.scale <- !is.null(object$S.terms)) {
      ## allows NAs to pass through to fit, se.fit, lwr and upr:
      S.mf <- model.frame(object$S.terms, newdata,
                          na.action=na.pass,
                          xlev=object$S.xlevels)
      ## model.frame will warn, but here we also throw an error:
      if(nrow(S.mf) != nrow(newdata))
        stop("length of variable(s) found do not match nrow(newdata)")
      if (!is.null(cl <- attr(object$S.terms, "dataClasses")))
        .checkMFClasses(cl, S.mf)
      S <- model.matrix(object$S.terms, S.mf,
                        contrasts=object$S.contrasts)
      Sint <- match("(Intercept)", colnames(S), nomatch = 0L)
      if(Sint <= 0) S <- cbind("(Intercept)" = rep(1, n), S)
      if(sum(object$aliased$zeta) > 0)
        S <- S[, !c(FALSE, object$aliased$zeta), drop=FALSE]
      Soff <- rep(0, nrow(S))
      if(!is.null(off.num <- attr(object$S.terms, "offset"))) 
        for(i in off.num) Soff <- Soff +
          eval(attr(object$S.terms, "variables")[[i + 1]], newdata)
    }
### Construct model environment:
    env <- eclm.newRho(parent.frame(), y=y, X=X,
                       NOM=if(is.nom) NOM else NULL,
                       S=if(is.scale) S else NULL,
                       weights=rep(1, n), offset=offset,
                       S.offset=if(is.scale) Soff else rep(0, n),
                       tJac=object$tJac)
    setLinks(env, link=object$link)
  } ## end !missing(newdata) or type == "class"
  else  env <- update(object, doFit=FALSE)
  env$par <- as.vector(coef(object))
  env$par <- env$par[!is.na(env$par)]
### FIXME: better way to handle NAs in coef?
  ## if(length(env$par) != ncol(env$B1))
  ##   stop(gettextf("design matrix has %d columns, but expecting %d (number of parameters)",
  ##                 ncol(env$B1), length(env$par)))
## Get predictions:
  pred <- 
    switch(type,
           "prob" = prob.predict.clm(env=env, cov=cov, se.fit=se.fit,
             interval=interval, level=level), 
           "class" = prob.predict.clm(env=env, cov=cov, se.fit=se.fit,
             interval=interval, level=level), 
           "cum.prob" = cum.prob.predict.clm(env=env, cov=cov,
             se.fit=se.fit, interval=interval, level=level),
           "linear.predictor" = lin.pred.predict.clm(env=env, cov=cov,
             se.fit=se.fit, interval=interval, level=level)
           )
### Arrange predictions in matrices if response is missing from
### newdata arg or type=="class":
  if(!has.response || type == "class") {
    pred <- lapply(pred, function(x) {
      x <- matrix(unlist(x), ncol=nlev, byrow=TRUE)
      dimnames(x) <- list(1:nrow(x), ylev)
      x
    })
    if(type == "class")
      pred <- lapply(pred, function(x) {
        factor(max.col(x), levels=seq_along(ylev), labels=ylev) })
  }
### Filter missing values (if relevant):
  if(missing(newdata) && !is.null(object$na.action))
    pred <- lapply(pred, function(x) napredict(object$na.action, x))
  return(pred)
}

prob.predict.clm <-
  function(env, cov, se.fit=FALSE, interval=FALSE, level=0.95)
### Works for linear and scale models:
### env - model environment with par set.
### cov - vcov for the parameters
{
  ## evaluate nll and grad to set dpi.psi in env:
  eclm.nll(env)
  pred <- list(fit = as.vector(env$fitted))
  if(se.fit || interval) {
    se.pr <- get.se(env, cov, type="prob")
    if(se.fit)
      pred$se.fit <- se.pr
    if(interval) {
      pred.logit <- qlogis(pred$fit)
      ## se.logit <- dlogis(pred$fit) * se.pr
      se.logit <- se.pr / (pred$fit * (1 - pred$fit))
      a <- (1 - level)/2
      pred$lwr <- plogis(pred.logit + qnorm(a) * se.logit)
      pred$upr <- plogis(pred.logit - qnorm(a) * se.logit)
    }
  }
  return(pred)
}

lin.pred.predict.clm <-
  function(env, cov, se.fit=FALSE, interval=FALSE, level=0.95)
### get predictions on the scale of the linear predictor
{
  ## evaluate nll and grad to set dpi.psi in env:
  eclm.nll(env)
  pred <- list(eta1=env$eta1, eta2=env$eta2)
  if(se.fit || interval) {
    se <- get.se(env, cov, type="eta")
    if(se.fit) {
      pred$se.eta1 <- se[[1]]
      pred$se.eta2 <- se[[2]]
    }
    if(interval) {
      a <- (1 - level)/2
      pred$lwr1 <- env$eta1 + qnorm(a) * se[[1]]
      pred$lwr2 <- env$eta2 + qnorm(a) * se[[2]]
      pred$upr1 <- env$eta1 - qnorm(a) * se[[1]]
      pred$upr2 <- env$eta2 - qnorm(a) * se[[2]]
    }
  }
  return(pred) ## list with predictions.
}

cum.prob.predict.clm <-
  function(env, cov, se.fit=FALSE, interval=FALSE, level=0.95)
{
  ## evaluate nll and grad to set dpi.psi in env:
  eclm.nll(env)
  pred <- list(cprob1=env$pfun(env$eta1), cprob2=env$pfun(env$eta2)) 
  if(se.fit || interval) {
    se <- get.se(env, cov, type="gamma")
    if(se.fit) {
      pred$se.cprob1 <- se[[1]]
      pred$se.cprob2 <- se[[2]]
    }
    if(interval) {
      a <- (1 - level)/2
      pred$lwr1 <- pred$cprob1 + qnorm(a) * se[[1]]
      pred$lwr2 <- pred$cprob2 + qnorm(a) * se[[2]]
      pred$upr1 <- pred$cprob1 - qnorm(a) * se[[1]]
      pred$upr2 <- pred$cprob2 - qnorm(a) * se[[2]]
    }
  }
  return(pred)
}

get.se <- function(rho, cov, type=c("eta", "gamma", "prob")) {
### Computes standard errors of predicted probabilities (prob),
### cumulative probabilities (gamma) or values of the linear
### predictor (eta) for linear (k<=0) or location-scale models
### (k>0). 
  rho$type <- match.arg(type)
  rho$cov <- cov
  eclm.nll(rho) ## just to be safe
  with(rho, {
### First compute d[eta, gamma, prob] / d par; then compute variance
### covariance matrix of the observations and extract SEs as the
### square root of the diagonal elements:
    if(type %in% c("eta", "gamma")) {
      D1 <- B1
      D2 <- B2
      if(k > 0) {
        D1 <- cbind(D1/sigma, -S*eta1)
        D2 <- cbind(D2/sigma, -S*eta2)
      }
      if(type == "gamma") {
        D1 <- D1*dfun(eta1)
        D2 <- D2*dfun(eta2)
      }
      se <- list(se1=sqrt(diag(D1 %*% cov %*% D1)),
                 se2=sqrt(diag(D2 %*% cov %*% D2)))
    }
    if(type == "prob") {
      p1 <- dfun(eta1)
      p2 <- dfun(eta2)
      C2 <- if(k <= 0) B1*p1 - B2*p2 else
      cbind(B1*p1/sigma - B2*p2/sigma,
            -(eta1 * p1 - eta2 * p2) * S)
      se <- sqrt(diag(C2 %*% cov %*% t(C2)))
    }
    return(se)
  })
}

## clm.nll.pred <- function(rho) { 
##   with(rho, {
##     eta1 <- drop(B1 %*% par) + o1
##     eta2 <- drop(B2 %*% par) + o2
##     fitted <- pfun(eta1) - pfun(eta2)
##     ## -sum(wts * log(fitted))
##   })
## }

