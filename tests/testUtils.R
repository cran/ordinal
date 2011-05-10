
## Not all R-functions give the right answer, so these tests are disabled.
library(ordinal)
## 
## stopifnot(isTRUE(all.equal(ordinal:::pgumbelR(-10:10), ordinal:::pgumbel(-10:10))))
## stopifnot(isTRUE(all.equal(ordinal:::pgumbel2R(-10:10), ordinal:::pgumbel2(-10:10))))
## stopifnot(isTRUE(all.equal(ordinal:::pAOR(-10:10, 1), ordinal:::pAO(-10:10, 1))))
## stopifnot(isTRUE(all.equal(ordinal:::plgammaR(-10:10, 1), ordinal:::plgamma(-10:10, 1))))
## 
## stopifnot(isTRUE(all.equal(ordinal:::dgumbelR(-10:10), ordinal:::dgumbel(-10:10))))
## stopifnot(isTRUE(all.equal(ordinal:::dgumbel2R(-10:10), ordinal:::dgumbel2(-10:10))))
## stopifnot(isTRUE(all.equal(ordinal:::dAOR(-10:10, 1), ordinal:::dAO(-10:10, 1))))
## stopifnot(isTRUE(all.equal(ordinal:::dlgammaR(-10:10, 1), ordinal:::dlgamma(-10:10, 1))))
## 
## stopifnot(isTRUE(all.equal(ordinal:::ggumbelR(-10:10), ordinal:::ggumbel(-10:10))))
## stopifnot(isTRUE(all.equal(ordinal:::ggumbel2R(-10:10), ordinal:::ggumbel2(-10:10))))
## stopifnot(isTRUE(all.equal(ordinal:::glogisR(-10:10), ordinal:::glogis(-10:10))))
## stopifnot(isTRUE(all.equal(ordinal:::gnormR(-10:10), ordinal:::gnorm(-10:10))))
## stopifnot(isTRUE(all.equal(ordinal:::gcauchyR(-10:10), ordinal:::gcauchy(-10:10))))
## stopifnot(isTRUE(all.equal(ordinal:::gAOR(-10:10, 1), ordinal:::gAO(-10:10, 1))))
## stopifnot(isTRUE(all.equal(ordinal:::glgammaR(-10:10, 1), ordinal:::glgamma(-10:10, 1))))

## Do we need a test suite for PFUN, DFUN and GFUN as well?

##################################################################


## None of these functions are exported from ordinal, so we would need
## to explicitly call e.g. ordinal:::pgumbel.


##   ## Testing link functions:
##   ## Test values:
##   nans <- c(NA, NaN)
##   infs <- c(-Inf, Inf)
##   usuals <- -10:10
##   limits <- c(t(c(-1,1) %o% 10^(1:10)))
##   limits <- limits[order(limits)]
##   all <- c(nans, infs, limits, usuals)
##   set.seed(12)
##   mixed <- sample(all)
##   all.list <- list(nans, infs, limits, usuals, mixed)
##   
##   ## Usual links:
##   ## pfun:
##   pfuns <- c("pnorm", "plogis", "pgumbel", "pgumbel2", "pcauchy")
##   ## dfun:
##   dfuns <- c("dnorm", "dlogis", "dgumbel", "dgumbel2", "dcauchy")
##   ## gfun:
##   gfuns <- c("gnorm", "glogis", "ggumbel", "ggumbel", "gcauchy")
##   
##   pres <- sapply(pfuns, function(pfun) {
##     lapply(all, function(vals) {
##       y <- list(as.name(pfun), vals)
##       y <- as.call(y)
##       eval(y)
##     } ) } )
##   rownames(pres) <- all
##   pres
##   
##   dres <- sapply(dfuns, function(fun) {
##     lapply(all, function(vals) {
##       y <- list(as.name(fun), vals)
##       y <- as.call(y)
##       eval(y)
##     } ) } )
##   rownames(dres) <- all
##   dres
##   
##   gres <- sapply(gfuns, function(fun) {
##     lapply(all, function(vals) {
##       y <- list(as.name(fun), vals)
##       y <- as.call(y)
##       eval(y)
##     } ) } )
##   rownames(gres) <- all
##   gres
##   
##   ## Flexible link functions:
##   
##   res <- sapply(usuals, function(lambda) plgamma(all, lambda))
##   res <- sapply(usuals, function(lambda) dlgamma(all, lambda))
##   res <- sapply(usuals, function(lambda) glgamma(all, lambda))
##   rownames(res) <- all
##   colnames(res) <- usuals
##   round(res, 4)
##   ## What if lambda is Inf or -Inf?
##   
##   
##   res <- sapply(usuals[usuals > 0], function(lambda) pAO(all, lambda))
##   res <- sapply(usuals[usuals > 0], function(lambda) dAO(all, lambda))
##   res <- sapply(usuals[usuals > 0], function(lambda) gAO(all, lambda))
##   rownames(res) <- all
##   colnames(res) <- usuals[usuals > 0]
##   round(res, 4)
##   
##   ### Not exactly sure how to make hard-coded tests for all these
##   ### cases... 
##   
##   ## These CDFs, PDFs and their gradients are designed to work for NA,
##   ## NaN, Inf, -Inf and a give reasonable results for extreme
##   ## arguments.
##   
##   #################################
##   ## Plots:
##   z <- seq(-10, 10, len = 1e3)
##   plot(z, gnorm(z), type = "l", ylim = c(-.4, .4))
##   lines(z, glogis(z))
##   lines(z, ggumbel(z), lty = 2)
##   lines(z, ggumbel2(z), lty = 3)
##   lines(z, gcauchy(z), lty = 3, col = "red")


