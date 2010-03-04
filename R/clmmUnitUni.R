.negLogLikM <- function(rho) { ## negative log-likelihood
    with(rho, {
        locat <- Loffset
        locat <- locat + u[grFac] * exp(par[p+nxi+k+estimLambda+ 1:s])
        sigma <- exp(Soffset)
        if(estimLambda > 0)
            lambda <- par[nxi + p + k + 1:estimLambda]
        if(k > 0)
            sigma <- sigma * exp(drop(Z %*% par[nxi+p + 1:k]))
        eta1 <- (drop(B1 %*% par[1:(nxi + p)]) + o1 - locat)/sigma
        eta2 <- (drop(B2 %*% par[1:(nxi + p)]) + o2 - locat)/sigma
        if(nlambda)
            pr <- pfun(eta1, lambda) - pfun(eta2, lambda)
        else
            pr <- pfun(eta1) - pfun(eta2)
        if (all(pr > 0))
            -sum(weights * log(pr)) -
                sum(dnorm(x=u, mean=0, sd=1, log=TRUE))
        else Inf
    })
}

.gradM <- function(rho) { ## gradient of the negative log-likelihood
    with(rho, {
        if(nlambda) {
            p1 <- dfun(eta1, lambda)
            p2 <- dfun(eta2, lambda)
        }
        else {
            p1 <- dfun(eta1)
            p2 <- dfun(eta2)
        }
        wtprSig <- weights/pr/sigma
        tapply(exp(par[nxi+p+k+estimLambda+ 1:s]) * wtprSig * (p1 - p2),
               grFac, sum) + u
### Potentially use rowsum instead of tapply - it might be faster
    })
}

.hessianM <- function(rho)  ## hessian of the negative log-likelihood
    with(rho,{
        if(nlambda) {
            g1 <- gfun(eta1, lambda)
            g2 <- gfun(eta2, lambda)
        }
        else {
            g1 <- gfun(eta1)
            g2 <- gfun(eta2)
        }
        tapply(((p1 - p2)^2 / pr - g1 + g2) * wtprSig, grFac, sum) *
            exp(2 * par[nxi+p+k+estimLambda+ 1:s]) +1
    })

### Potentially use rowsum instead of tapply - it might be faster
##    with(rho,
##          tapply(((p1 - p2)^2 / pr) * wtprSig,
##                 grFac, sum) +  exp(-2*par[p+q+k+ 1:s]) )

update.u <- function(rho)
{
    ctrl <- rho$ctrl
    stepFactor <- 1
    innerIter <- 0
    rho$u <- rho$uStart
    rho$negLogLik <- .negLogLikM(rho)
    if(rho$negLogLik == Inf)
        stop("Non-finite log-likelihood at starting value")
### FIXME: Is this the most reasonable action to take? Could this case
### be handled by letting getNLA return Inf in this case?
    rho$gradient <- .gradM(rho)
    maxGrad <- max(abs(rho$gradient))
    if(ctrl$trace > 0) {
        Trace(iter=0, stepFactor, rho$negLogLik, maxGrad, rho$u, first=TRUE)
        conv <- -1  ## Convergence flag
        message <- "iteration limit reached"
    }

    ## Newton-Raphson algorithm:
    for(i in 1:rho$ctrl$maxIter) {
        if(maxGrad < rho$ctrl$gradTol) {
            if(rho$ctrl$trace > 0) {
                message <- "max|gradient| < tol, so current iterate is probably solution"
                cat("\nOptimizer converged! ", "max|grad|:",
                    maxGrad, message, fill = TRUE)
                conv <- 0
            }
            break
        }
        rho$D <- .hessianM(rho)
        step <- rho$gradient / rho$D
        rho$u <- rho$u - stepFactor * step
        negLogLikTry <- .negLogLikM(rho)
        lineIter <- 0
        ## simple line search, i.e. step halfing:
        while(negLogLikTry > rho$negLogLik) {
            stepFactor <- stepFactor/2
            rho$u <- rho$u + stepFactor * step
            negLogLikTry <- .negLogLikM(rho)
            lineIter <- lineIter + 1
            if(rho$ctrl$trace > 0)
                Trace(i+innerIter, stepFactor, rho$negLogLik, maxGrad,
                      rho$u, first=FALSE)
            if(lineIter > rho$ctrl$maxLineIter){
                message <- "step factor reduced below minimum"
                conv <- 1
                break
            }
            innerIter <- innerIter + 1
        }
        rho$negLogLik <- negLogLikTry
        rho$gradient <- .gradM(rho)
        maxGrad <- max(abs(rho$gradient))
        if(ctrl$trace > 0)
            Trace(i+innerIter, stepFactor, rho$negLogLik, maxGrad, rho$u, first=FALSE)
        stepFactor <- min(1, 2 * stepFactor)
    }
    if(ctrl$trace > 0)
        if(conv != 0) cat(message, fill = TRUE)
    rho$Niter <- rho$Niter + i
    rho$D <- .hessianM(rho)
}

clmm <-
  function(location, scale, nominal, random, data, weights, start, subset,
           na.action, contrasts = NULL, Hess = FALSE, model = TRUE,
           method = c("ucminf", "nlminb"),
           link = c("logistic", "probit", "cloglog", "loglog",
           "cauchit", "Aranda-Ordaz", "log-gamma"), lambda = NULL,
           doFit = TRUE,  control = list(), nAGQ = 1,
           threshold = c("flexible", "symmetric", "equidistant"), ...)
    ## Handle if model = FALSE
{
    R <- Call <- match.call(expand.dots = FALSE)
    if(missing(random)) {
        Call[[1]] <- as.name("clm")
        return(eval.parent(Call))
    }
    Call$random <- Call$control <- Call$start <- NULL
    Call$doFit <- Call$Hess <- FALSE
    Call[[1]] <- as.name("clm")
    rhoM <- eval.parent(Call)
    rhoM$randomName <- as.character(R$random)
    rhoM$call <- match.call()
### Set grouping factor:
    if (missing(data)) data <- environment(location)
    if (is.matrix(eval.parent(R$data)))
        R$data <- as.data.frame(data)
    m <- match(c("location", "random", "data", "subset",
                 "weights", "na.action"), names(R), 0)
    R <- R[c(1, m)]
    R$drop.unused.levels <- TRUE
    names(R)[names(R) == "location"] <- "formula"
    R[[1]] <- as.name("model.frame")
    R <- eval.parent(R)
    rhoM$grFac <- R[,"(random)"]
    with(rhoM, {
        r <- nlevels(grFac) ## no. random effects
        if(r <= 2)
            stop("Grouping factor must have 3 or more levels")
        s <- 1 ## no. variance parameters
        Niter <- 0L
    })
### set starting values:
    if(missing(start)) {
        fitCLM(rhoM)
        rhoM$start <- rhoM$par <- c(rhoM$par, log(1))
    }
    else
        rhoM$start <- rhoM$par <- start
    rhoM$uStart <- rhoM$u <- rep(0, rhoM$r)
### Test starting values:
    if(length(rhoM$start) != with(rhoM, nxi + p + k + estimLambda + s))
        stop("'start' is not of the correct length")
    if(rhoM$ncolXX == 0) {
        if(!all(diff(c(rhoM$tJac %*% rhoM$start[1:rhoM$nalpha])) > 0))
            stop("Threshold starting values are not of increasing size")
    }
    if(!.negLogLikM(rhoM) < Inf)
        stop("Non-finite log-likelihood at starting values")
### Set AGQ parameters:
    ObjFun <- getNLA
    rhoM$nAGQ <- as.integer(nAGQ)
    if(rhoM$nAGQ >= 2) {
        ObjFun <- getNAGQ
        rhoM$PRnn <- array(0, dim=c(rhoM$nobs, rhoM$nAGQ))
        rhoM$PRrn <- array(0, dim=c(rhoM$r, rhoM$nAGQ))
        ghq <- gauss.hermite(rhoM$nAGQ)
        rhoM$ghqns <- ghq$nodes
        rhoM$ghqws <- ghq$weights * exp(rhoM$ghqns^2)
    }
### Set control parameters:
    rhoM$ctrl <- list(trace = 0, maxIter = 10,
                      gradTol = 1e-3, maxLineIter = 20)
    if("trace" %in% names(control) && control["trace"] < 0) {
        rhoM$ctrl$trace <- 1
        control$trace <- 1
    }
    if(rhoM$method == "ucminf") {
        ctrl <- list(grad="central", grtol=1e-4)
        control <- c(control, ctrl[!names(ctrl) %in% names(control)])
        rhoM$ctrl$grtol <- control$grtol
    }
### Fit the model:
    if(!doFit)
        return(rhoM)
    update.u(rhoM) # Try updating the random effects
    rhoM$optRes <- switch(rhoM$method,
                       "ucminf" = ucminf(rhoM$start, function(x)
                       ObjFun(rhoM, x), control=control),
                       "nlminb" = nlminb(rhoM$start, function(x)
                       ObjFun(rhoM, x), control=control))
    rhoM$par <- rhoM$optRes[[1]]
    if(Hess) {
        rhoM$Hessian <- hessian(function(x) getNLA(rhoM, x),
                                method.args = list(r=2), rhoM$par)
        rhoM$par <- rhoM$optRes[[1]]
    }
    update.u(rhoM) # Makes sure ranef's are evaluated at the optimum
### Post processing:
    res <- finalizeRhoM(rhoM)
    res$call <- match.call()
    res$na.action <- attr(R, "na.action")
    res$contrasts <- contrasts
    class(res) <- c("clmm", "clm")
    res
}

getNLA <- function(rho, par) {
    if(!missing(par))
        rho$par <- par
    update.u(rho)
    logDetD <- sum(log(rho$D/(2*pi)))
    rho$negLogLik + logDetD / 2
}

getNAGQ <- function(rho, par) {
    if(!missing(par))
        rho$par <- par
    update.u(rho)
    with(rho, {
        K <- sqrt(2/D)
        agqws <- K %*% t(ghqws)
        agqns <- apply(K %*% t(ghqns), 2, function(x) x + u)
        locat2 <- locat - u[grFac] * exp(par[nxi+p+k+estimLambda+ 1:s])
        locatTmp <- apply(agqns, 2, function(x) locat2
                          + x[grFac] * exp(par[nxi+p+k+estimLambda+ 1:s]))

        eta1Tmp <- (drop(B1 %*% par[1:(nxi + p)]) + o1 - locatTmp)/sigma
        eta2Tmp <- (drop(B2 %*% par[1:(nxi + p)]) + o2 - locatTmp)/sigma
        if(nlambda)
            PRnn <- (pfun(eta1Tmp, lambda) - pfun(eta2Tmp, lambda))^weights
        else
            PRnn <- (pfun(eta1Tmp) - pfun(eta2Tmp))^weights
### Potentially a more numerically safe solution:
        ## PRnn <- exp(weights * log(pfun(eta1Tmp) - pfun(eta2Tmp)))
        for(i in 1:r)
            PRrn[i,] <- apply(PRnn[unclass(grFac) == i, ], 2, prod)
        PRrn <- PRrn * agqws * dnorm(x=agqns, mean=0, sd=1)
    })
    -sum(log(rowSums(rho$PRrn)))
}

finalizeRhoM <- function(rhoM) {
    if(rhoM$method == "ucminf") {
### FIXME: Test of convergence would be nice with nlminb as well.
        if(rhoM$optRes$info[1] > rhoM$ctrl[["grtol"]])
            warning(sprintf("Optimizer 'ucminf' terminated with max|gradient|: %e",
                            rhoM$optRes$info[1]), call.=FALSE)
        rhoM$convergence <-
            ifelse(rhoM$optRes$info[1] > rhoM$ctrl[["grtol"]], FALSE, TRUE)
    }

    with(rhoM, {
        xi <- par[1:nxi]
        names(xi) <- xiNames
        thetaNames <- paste(lev[-length(lev)], lev[-1], sep="|")
        Alpha <- Theta <- matrix(par[1:nxi], nrow=ncolXX, byrow=TRUE)
        Theta <- t(apply(Theta, 1, function(x) c(tJac %*% x)))
        if(ncolXX > 1){
            dimnames(Theta) <- list(dnXX[[2]], thetaNames)
            dimnames(Alpha) <- list(dnXX[[2]], alphaNames)
        }
        else {
            Theta <- c(Theta)
            Alpha <- c(Alpha)
            names(Theta) <- thetaNames
            names(Alpha) <- alphaNames
        }
        coefficients <- xi
        if(p > 0) {
            beta <- par[nxi + 1:p]
            names(beta) <- dnX[[2]]
            coefficients <- c(coefficients, beta)
        }
        if(k > 0) {
            zeta <- par[nxi+p + 1:k]
            names(zeta) <- dnZ[[2]]
            coefficients <- c(coefficients, zeta)
        }
        if(estimLambda > 0) {
            names(lambda) <- "lambda"
            coefficients <- c(coefficients, lambda)
        }
        if(s > 0) {
            stDev <- exp(par[nxi+p+k + estimLambda + 1:s])
            names(stDev) <- randomName
            coefficients <- c(coefficients, stDev)
        }
        if(exists("Hessian", inherits=FALSE))
            dimnames(Hessian) <- list(names(coefficients),
                                      names(coefficients))
        edf <- p + nxi + k + estimLambda + s
        nobs <- sum(weights)
        fitted.values <- pr
        df.residual = nobs - edf
        ranef <- u * exp(par[nxi+p+k + estimLambda + 1:s])
        condVar <- 1/D * exp(2*par[nxi+p+k + estimLambda + 1:s])
        logLik <- -optRes[[2]]
    })
    res <- as.list(rhoM)
###     notKeep <- c("wtpr", "g", "g1", "g2", "dpi.psi", "X", "Z",
###                  "dS.psi", "p2", "p1", "pr", "eta1", "eta2", "gfun",
###                  "dfun", "pfun", "dnX", "B1", "B2", "locat", "ctheta",
###                  "prSig", "gradSigma", "gradThetaBeta", "sigma",
###                  "Soffset", "Loffset", "dnZ", "weights", "par",
###                  "start", "bStart", "u", "Hess",
###                  "D", "wtprSig", "negLogLik", "grLev", "randomName",
###                  "convTol")
###     keep <- names(res)[!(names(res) %in% notKeep)]
    keepNames <-
        c("ranef", "df.residual", "fitted.values", "edf", "start",
          "stDev", "beta", "coefficients", "zeta", "Alpha", "Theta",
          "xi", "lambda", "convergence", "Hessian", # "convTol",
          "gradient", "optRes", "logLik", "Niter", "grFac", "call",
          "scale", "location", "nominal", "method", "y", "lev",
          "nobs", "threshold", "estimLambda", "link", "nAGQ",
          "condVar", "contrasts", "na.action")
    m <- match(keepNames, names(res), 0)
    names(res)
    res <- res[m]
    ## class(res) <- "clm.fit"
    res
}

anova.clmm <- function (object, ..., test = c("Chisq", "none"))
{
    anova.clm(object, ..., test = c("Chisq", "none"))
}

print.clmm <- function(x, ...)
{
    if(x$nAGQ >=2)
        cat(paste("Cumulative Link Mixed Model fitted with the adaptive",
                  "Gauss-Hermite \nquadrature approximation with",
                  x$nAGQ ,"quadrature points"), "\n\n")
    else
        cat("Cumulative Link Mixed Model fitted with the Laplace approximation\n",
            fill=TRUE)
    if(!is.null(cl <- x$call)) {
        cat("Call:\n")
        dput(cl, control=NULL)
    }
    if(length(x$stDev)) {
        cat("\nRandom effects:\n")
        varMat <- matrix(c(x$stDev^2, x$stDev), nrow =
                         length(x$stDev), ncol=2)
        rownames(varMat) <- names(x$stDev)
        colnames(varMat) <- c("Var", "Std.Dev")
        print(varMat, ...)
    } else {
        cat("\nNo random effects\n")
    }
    if(length(x$beta)) {
        cat("\nLocation coefficients:\n")
        print(x$beta, ...)
    } else {
        cat("\nNo location coefficients\n")
    }
    if(length(x$zeta)) {
        cat("\nScale coefficients:\n")
        print(x$zeta, ...)
    } else {
        cat("\nNo Scale coefficients\n")
    }
    if(x$estimLambda > 0) {
        cat("\nLink coefficient:\n")
        print(x$lambda)
    }
    cat("\nThreshold coefficients:\n")
    print(x$Alpha, ...)
    if(x$threshold != "flexible") {
        cat("\nThresholds:\n")
        print(x$Theta, ...)
    }
    cat("\nlog-likelihood:", format(x$logLik, nsmall=2), "\n")
    cat("AIC:", format(-2*x$logLik + 2*x$edf, nsmall=2), "\n")
    if(nzchar(mess <- naprint(x$na.action))) cat("(", mess, ")\n", sep="")
    invisible(x)
}

vcov.clmm <- function(object, ...)
{
    if(is.null(object$Hessian)) {
        stop("Model needs to be fitted with Hess = TRUE")
    }
    dn <- names(object$coefficients)
    structure(solve(object$Hessian), dimnames = list(dn, dn))
}

summary.clmm <- function(object, digits = max(3, .Options$digits - 3),
                         correlation = FALSE, ...)
{
    coef <- with(object,
                 matrix(0, edf-length(stDev), 4,
                        dimnames = list(names(coefficients[-edf]),
                        c("Estimate", "Std. Error", "z value", "Pr(>|z|)"))))
    coef[, 1] <- object$coefficients[-object$edf]
    vc <- vcov(object)
    sd <- sqrt(diag(vc))
    coef[, 2] <- sd[-object$edf]
    coef[, 3] <- coef[, 1]/coef[, 2]
    coef[, 4] <- 2*pnorm(abs(coef[, 3]), lower.tail=FALSE)
    object$coefficients <- coef
    object$digits <- digits
    varMat <- matrix(c(object$stDev^2, object$stDev),
                     nrow = length(object$stDev), ncol=2)
    rownames(varMat) <- names(object$stDev)
    colnames(varMat) <- c("Var", "Std.Dev")
    object$varMat <- varMat

    if(correlation)
        object$correlation <- (vc/sd)/rep(sd, rep(object$edf, object$edf))
    class(object) <- "summary.clmm"
    object
}

print.summary.clmm <- function(x, digits = x$digits, signif.stars =
                              getOption("show.signif.stars"), ...)
{
    cat("Cumulative Link Mixed Model fitted with the Laplace approximation",
        fill=TRUE)
    if(!is.null(cl <- x$call)) {
        cat("Call:\n")
        dput(cl, control=NULL)
    }
    if(length(x$stDev)) {
        cat("\nRandom effects:\n")
        print(x$varMat, ...)
    } else {
        cat("\nNo random effects\n")
    }
### FIXME: Should the number of obs. and the number of groups be
### displayed as in lmer?
    coef <- format(round(x$coefficients, digits=digits))
    coef[,4] <- format.pval(x$coefficients[, 4])
    p <- length(x$beta); nxi <- length(x$xi)
    k <- length(x$zeta); u <- x$estimLambda
    if(p > 0) {
        cat("\nLocation coefficients:\n")
        print(coef[nxi + 1:p, , drop=FALSE],
              quote = FALSE, ...)
    } else {
        cat("\nNo location coefficients\n")
    }
    if(k > 0) {
      cat("\nScale coefficients:\n")
      print(coef[(nxi+p+1):(nxi+p+k), , drop=FALSE],
            quote = FALSE, ...)
    } else {
      cat("\nNo scale coefficients\n")
    }
    if(x$estimLambda > 0) {
        cat("\nLink coefficients:\n")
        print(coef[(nxi+p+k+1):(nxi+p+k+u), , drop=FALSE],
              quote = FALSE, ...)
    }
    cat("\nThreshold coefficients:\n")
    print(coef[1:nxi, -4, drop=FALSE], quote = FALSE, ...)

    cat("\nlog-likelihood:", format(x$logLik, nsmall=2), "\n")
    cat("AIC:", format(-2*x$logLik + 2*x$edf, nsmall=2), "\n")
    if(nzchar(mess <- naprint(x$na.action))) cat("(", mess, ")\n", sep="")
    if(!is.null(correl <- x$correlation)) {
        cat("\nCorrelation of Coefficients:\n")
        ll <- lower.tri(correl)
        correl[ll] <- format(round(correl[ll], digits))
        correl[!ll] <- ""
        print(correl[-1, -ncol(correl)], quote = FALSE, ...)
    }
    invisible(x)
}

## ranef.clmm <- function(x) {
##     x$ranef
## }

Trace <- function(iter, stepFactor, val, maxGrad, par, first=FALSE) {
    t1 <- sprintf(" %3d:     %-5e:   %.3f:   %1.3e:  ",
                  iter, stepFactor, val, maxGrad)
    t2 <- formatC(par)
    if(first)
        cat("iter:  step factor:     Value:   max|grad|:   Parameters:\n")
    cat(t1, t2, "\n")
}

gauss.hermite <- function (n)
{
    n <- as.integer(n)
    if (n < 0)
        stop("need non-negative number of nodes")
    if (n == 0)
        return(list(nodes = numeric(0), weights = numeric(0)))
    i <- 1:n
    i1 <- i[-n]
    muzero <- sqrt(pi)
    a <- rep(0, n)
    b <- sqrt(i1/2)

    A <- rep(0, n * n)
    A[(n + 1) * (i1 - 1) + 2] <- b
    A[(n + 1) * i1] <- b
    dim(A) <- c(n, n)
    vd <- eigen(A, symmetric = TRUE)
    w <- rev(as.vector(vd$vectors[1, ]))
    w <- muzero * w^2
    x <- rev(vd$values)
    list(nodes = x, weights = w)
}

