convergence.clmm <-
  function(object, digits = max(3, getOption("digits") - 3), ...)
### Results: data.frame with columns:
### Estimate
### Std. Error
### Gradient - gradient of the coefficients at optimizer termination 
### Error - the signed error in the coefficients at termination
### Rel. Error - the relative error in the coefficeints at termination
###
### The (signed) Error is determined as the Newton step, so this is
### only valid close to the optimum where the likelihood function is
### quadratic.
###
### The relative error equals step/Estimate.
{
  ## get gradient:
  if(is.null(object$Hessian))
    stop("Model has to be fitted with 'Hess = TRUE'")
  summ <- summary(object)
  rho <- update(object, doFit=FALSE)
  rho$par <- as.vector(coef(object))
  if(rho$ntau == 1L) {
    obj.fun <-
      if(object$nAGQ < 0) getNGHQ.ssr
      else if(object$nAGQ > 1) getNAGQ.ssr
      else getNLA.ssr ## nAGQ %in% c(0, 1)
  }
  else
    obj.fun <- getNLA
  g <- grad(function(p) obj.fun(rho, p), x=as.vector(coef(object)))
  h <- object$Hessian
  info <- summ$info[c("nobs", "logLik", "niter", "max.grad",
                      "cond.H")] 
  ## Compute approximate error in the coefficients:
  step <- solve(h, g)
  if(max(abs(step)) > 1e-2)
    warning("convergence assessment may be unreliable due to large numerical error")
  ## Compute approximate error in the log-likelihood function:
  rho$par <- coef(object) - step
  new.logLik <- -obj.fun(rho)
  logLik.err <- object$logLik - new.logLik
  if(new.logLik < object$logLik)
    stop("Cannot assess convergence: ",
         "please assess the likelihood with slice()")
  info$logLik.Error <- formatC(logLik.err, digits=2, format="e")
  se <- sqrt(diag(vcov(object)))
  tab <- cbind(coef(object), se, g, step, cor.dec(step),
               signif.digits(coef(object), step)) 
  dimnames(tab) <-
    list(names(coef(object)),
         c("Estimate", "Std.Err", "Gradient",
           "Error", "Cor.Dec", "Sig.Dig"))
  tab.print <- tab
  for(i in 1:2) 
    tab.print[,i] <- format(c(tab[,i]), digits=digits)
  for(i in 3:4) tab.print[,i] <-
    format(c(tab[,i]), digits=max(1, digits - 1))
  print(info, row.names=FALSE, right=FALSE)
  cat("\n")
  print(tab.print, quote=FALSE, right=TRUE, ...)
  e.val <- eigen(object$Hessian, symmetric=TRUE,
                 only.values=TRUE)$values
  cat("\nEigen values of Hessian:\n")
  cat(format(e.val, digits=digits), "\n", fill=TRUE)
  if(any(e.val <=0))
    cat("\nNegative eigen values occured so model did not converge\n")
  return(invisible(tab))  
}

nll.u.ssr <- function(rho) {
  with(rho, {
    tau <- exp(par[nalpha + nbeta + 1:ntau])
    eta1 <- drop(B1 %*% par[1:(nalpha + nbeta)]) + o1 - u[grFac] * tau 
    eta2 <- drop(B2 %*% par[1:(nalpha + nbeta)]) + o2 - u[grFac] * tau
  })
  rho$pr <- getFittedC(rho$eta1, rho$eta2, rho$link)
  if(all(is.finite(rho$pr)) && all(rho$pr > 0))
    rho$nll <- -sum(rho$wts * log(rho$pr)) -
      sum(dnorm(x=rho$u, mean=0, sd=1, log=TRUE))
  else
    rho$nll <- Inf
  rho$nll
}  

jnll.u.ssr <- function(rho)
### Compute the contributions to the joint nll for each level of grFac
### in a clmm with a single RE term.
### result: a vector of jnll contributions for each level of grFac. 
{
  ## evaluate the multinomial contribution (y|u) to the joint log 
  ## likelihood: 
  with(rho, {
    tau <- exp(par[nalpha + nbeta + 1:ntau])
    eta1 <- drop(B1 %*% par[1:(nalpha + nbeta)]) + o1 - u[grFac] * tau 
    eta2 <- drop(B2 %*% par[1:(nalpha + nbeta)]) + o2 - u[grFac] * tau
  })
  rho$pr <- getFittedC(rho$eta1, rho$eta2, rho$link)
  ## split the contributions according to grFac:
  y.u <- split(rho$wts * log(rho$pr), rho$grFac)
  ## compute the contribution to the joint logLik for each level of
  ## grFac: 
  nll.contrib <- sapply(seq_along(y.u), function(i) {
    -sum(y.u[[i]]) - dnorm(x=rho$u[i], mean=0, sd=1, log=TRUE) })
  return(nll.contrib)
}

slice.u <-
  function(object, alpha=.001, grid=100, quad.approx=TRUE, ...)
### Compute the joint nll for a clmm object for each of the random
### effects. Also compute the corresponding quadratic approximation
### and the Laplace contributions to the (marginal) log-likelihood for
### each level of grFac.
###
### result: a list of data.frames with components: "u", "jnll",
### "quad" and some attributes.
### 
### FIXME: make it possible to select one or more random effects. 
{
  ## argument matching and testing:
  stopifnot(is.numeric(alpha) && alpha > 0)
  stopifnot(is.numeric(grid) && grid >= 1)
  grid <- as.integer(round(grid))
  ## get model environment:
  rho <- update(object, doFit=FALSE)
  if(rho$ntau != 1L)
    stop("only models with a single RE term are allowed") 
  rho$par <- as.vector(coef(object))
  ## update the conditional modes of the random effects:
  nllBase.uC(rho)
  update.uC(rho)
  ranef <- rho$u ## save mode of u
  cond.hess <- rho$D

  ## compute contributions to the Laplace likelihood:
  jnll <- jnll.u.ssr(rho)
  la.contrib <- -jnll - log(cond.hess/(2*pi)) / 2
  
  ## compute range of u-values at which to compute the joint nll: 
  lim <- c(1, -1) * qnorm(alpha/2)
  lims <- outer(lim, rho$u, FUN="+")
  ## compute the u-values at which to evaluate the log.lik
  u.vals <- lapply(seq_len(ncol(lims)), function(i) {
    seq(from=lims[1,i], to=lims[2,i], length.out=grid) })
  ## compute [n, nu] matrix of joint logLik:
  jll.mat <- t(apply(do.call(cbind, u.vals), 1, function(uu) {
    rho$u <- uu
    -jnll.u.ssr(rho)
  }))
  
  ## collect parameter sequences and relative logLik in a list of
  ## data.frames: 
  res <- lapply(seq_along(u.vals), function(i) {
    structure(data.frame(u.vals[[ i ]], jll.mat[,i]),
              names = c("u", "jll"))
  })

  names(res) <- levels(rho$grFac)
  attr(res, "original.fit") <- object
  attr(res, "la.contrib") <- la.contrib
  ## attr(res, "mode") <- modes
  class(res) <- "slice.clmm"

  if(!quad.approx) return(res)
  ## compute quadratic approximation to the joint nll
  for(i in seq_along(u.vals)) 
    res[[ i ]]$quad <-
      jll.mat[,i] - (u.vals[[ i ]] - ranef[i])^2 / cond.hess[i] / 2

  return(res)
}

plot.slice.clmm <-
  function(x, Log = FALSE,
           ask = prod(par("mfcol")) < length(x) && dev.interactive(),
           ...)
### Plot the joint log-likelihood or joint likelihood functions with
### the quadratic/gaussian approximations (if quad elements are
### present). 
{
  
  if(ask) {
    oask <- devAskNewPage(TRUE)
    on.exit(devAskNewPage(oask))
  }
 
  for(i in seq_along(x)) {
    z <- x[[ i ]]
    if(!Log) z[,2:3] <- exp(z[,2:3])
    plot(z[[1]], z[[2]], type="l", ylab="", xlab="", ...)
    if(!is.null(z$quad))
      lines(z[[1]], z[[3]], lty=2)
  }
  
  return(invisible())
}
