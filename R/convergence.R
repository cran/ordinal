convergence <- function(object, ...) {
  UseMethod("convergence")
}

convergence.clm <-
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
  summ <- summary(object)
  info <- summ$info[c("nobs", "logLik", "niter", "max.grad",
                      "cond.H")] 
  ## Compute approximate error in the coefficients:
  step <- with(object, solve(Hessian, gradient))
  if(max(abs(step)) > 1e-2)
    warning("convergence assessment may be unreliable due to large numerical error")
  ## Compute approximate error in the log-likelihood function:
  env <- update(object, doFit=FALSE)
  env$par <- coef(object, na.rm=TRUE) - step
  new.logLik <- -ordinal:::eclm.nll(env)
  logLik.err <- object$logLik - new.logLik
  if(new.logLik < object$logLik)
    stop("Cannot assess convergence: ",
         "please assess the likelihood with slice()")
  info$logLik.Error <- formatC(logLik.err, digits=2, format="e")
  tab <- coef(summ, na.rm=TRUE)[,1:2]
  se <- tab[,2]
    tab <- cbind(tab, object$gradient, step, cor.dec(step),
                 signif.digits(tab[,1], step)) 
  dimnames(tab) <-
    list(names(coef(object, na.rm=TRUE)),
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
  e.val <- eigen(object$Hessian, symmetric=TRUE, only.values=)$values
  cat("\nEigen values of Hessian:\n")
  cat(format(e.val, digits=digits), "\n", fill=TRUE)
  if(any(e.val <=0))
    cat("\nNegative eigen values occured so model did not converge\n")
  return(invisible(tab))  
}


cor.dec <- function(error) {
### computes the no. correct decimals in a number if 'error' is the
### error in the number. 
### The function is vectorized.
  xx <- -log(abs(error), 10)
  lead <- floor(xx)
  res <- ifelse(xx < lead - log(.5, 10), lead-1, lead)
  res[abs(error) >= .05] <- 0
  res
}

signif.digits <- function(value, error) {
### Determines the number of significant digits in 'value' if the
### absolute error in 'value' is 'error'.
### The function is vectorized.
  res <- cor.dec(error) + ceiling(log(abs(value), 10))
  res[res < 0] <- 0
  res
}
