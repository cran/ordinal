\name{clmm}
\alias{clmm}
\title{Cumulative link mixed models}
\description{
  Fits cumulative link mixed models, i.e. cumulative link models with
  random effects via the Laplace approximation or the standard and the
  adaptive Gauss-Hermite quadrature approximation. The functionality in
  \code{\link{clm}} is also implemented here. Currently only a single
  random term is allowed in the location-part of the model.
}
\usage{
clmm(location, scale, nominal, random, data, weights, start, subset,
     na.action, contrasts, Hess = FALSE, model = TRUE, sdFixed,
     link = c("logistic", "probit", "cloglog", "loglog",
     "cauchit", "Aranda-Ordaz", "log-gamma"), lambda,
     doFit = TRUE, control, nAGQ = 1,
     threshold = c("flexible", "symmetric", "equidistant"), ...)
}
\arguments{
  \item{location}{
    as in \code{\link{clm}}.
  }
  \item{scale}{
    as in \code{\link{clm}}.
  }
  \item{nominal}{
    as in \code{\link{clm}}.
  }
  \item{random}{
    a factor for the random effects in the location-part of the model.
  }
  \item{data}{
    as in \code{\link{clm}}.
}
\item{weights}{
    as in \code{\link{clm}}.
}
\item{start}{
  initial values for the parameters in the format
  \code{c(alpha, beta, log(zeta), lambda, log(stDev))} where
  \code{stDev} is the standard deviation of the random effects.
}
\item{subset}{
    as in \code{\link{clm}}.
}
\item{na.action}{
    as in \code{\link{clm}}.
}
\item{contrasts}{
    as in \code{\link{clm}}.
}
\item{Hess}{
  logical for whether the Hessian (the inverse of the observed
  information matrix) should be computed.
  Use \code{Hess = TRUE} if you intend to call \code{summary} or
  \code{vcov} on the fit and \code{Hess = FALSE} in all other instances
  to save computing time.
}
\item{model}{
  as in \code{\link{clm}}.
}
\item{sdFixed}{
  If \code{sdFixed} is specified (a positive scalar), a model is fitted
  where the standard deviation for the random term is fixed at the value
  of \code{sdFixed}. If \code{sdFixed} is left unspecified, the standard
  deviation of the random term is estimated from data.
}
\item{link}{
  as in \code{\link{clm}}.
}
\item{lambda}{
    as in \code{\link{clm}}.
}
\item{doFit}{
  as in \code{\link{clm}} although it can also be one of \code{c("no",
    "R" "C")}, where \code{"R"} use the R-implementation for fitting,
  \code{"C"} (default) use C-implementation for fitting and \code{"no"}
  behaves as \code{FALSE} and returns the environment. 
}
\item{control}{
  a call to \code{\link{clmm.control}}.
}
\item{threshold}{
    as in \code{\link{clm}}.
}
\item{nAGQ}{
    the number of quadrature points to be used in the adaptive
    Gauss-Hermite quadrature approximation to the marginal
    likelihood. Defaults to \code{1} which leads to the Laplace
    approximation. An odd number of quadrature points is encouraged and
    3, 5 or 7 are usually enough to achive high precision. Negative
    values give the standard, i.e. non-adaptive Gauss-Hermite
    quadrature. 
}
\item{\dots}{
  additional arguments are passed on to \code{\link{clm.control}} and
  possibly further on to the optimizer, which can lead to surprising
  error or warning messages when mistyping arguments etc.
}

}
\details{
  There are methods for the standard model-fitting functions, including
  \code{\link{summary}}, \code{\link{vcov}},
  \code{\link[=profile.clmm]{profile}},
  \code{\link[=profile.clmm]{plot.profile}},
  \code{\link[=confint.profile.clmm]{confint}},
  \code{\link[ordinal]{anova}}, \code{\link{logLik}},
  and an \code{extractAIC} method.

  A Newton scheme is used to obtain the conditional modes of the random
  effects for Laplace and AGQ approximations, and a non-linear
  optimization is performed over the fixed parameter set to get the
  maximum likelihood estimates.
  The Newton
  scheme uses the observed Hessian rather than the expected as is done
  in e.g. \code{\link[lme4]{glmer}}, so results from the Laplace
  approximation for binomial fits should in general be more precise -
  particularly for other links than the \code{"logistic"}.

  Core parts of the function are implemented in C-code for speed. For
  binomial fits with a single random factor \code{clmm} is now faster
  than \code{\link[lme4]{glmer}} for at least some models.

  The function calls \code{\link{clm}} to up an
  environment and to get starting values.
}
\value{
  If \code{doFit = FALSE} the result is an environment
  representing the model ready to be optimized.
  If \code{doFit = TRUE} the result is an
  object of class \code{"clmm"} with the following components:

  \item{stDev}{
    the standard deviation of the random effects.
  }
  \item{Niter}{
    the total number of iterations in the Newton updates of the
    conditional modes of the random effects.
  }
  \item{grFac}{
    the grouping factor defining the random effects.
  }
  \item{nAGQ}{
    the number of quadrature points used in the adaptive Gauss-Hermite
    Quadrature approximation to the marginal likelihood.
  }
  \item{condVar}{
    the conditional variances of the random effects at their conditional
    modes.
  }
  \item{beta}{the parameter estimates of the location part.
  }
  \item{zeta}{the parameter estimates of the scale part on the log
    scale; the scale parameter estimates on the original scale are given
    by \code{exp(zeta)}.
  }
  \item{Alpha}{vector or matrix of the threshold parameters.
  }
  \item{Theta}{vector or matrix of the thresholds.
  }
  \item{xi}{vector of threshold parameters, which, given a
    threshold function (e.g. \code{"equidistant"}), and possible nominal
    effects define the class boundaries, \code{Theta}.
  }
  \item{lambda}{the value of lambda if lambda is supplied or estimated,
    otherwise missing.
  }
  \item{coefficients}{the coefficients of the intercepts
    (\code{theta}), the location  (\code{beta}), the scale
    (\code{zeta}), and the link function parameter (\code{lambda}).
  }
  \item{df.residual}{the number of residual degrees of freedoms,
    calculated using the weights.
  }
  \item{fitted.values}{vector of fitted values for each
    observation. An observation here is each of the scalar elements of
    the multinomial table and not a multinomial vector.
  }
  \item{convergence}{\code{TRUE} if the optimizer terminates wihtout
    error and \code{FALSE} otherwise.
  }
  \item{gradient}{vector of gradients for the unit-variance random
    effects at their conditional modes.
  }
  \item{optRes}{list with results from the optimizer. The contents of the
    list depends on the choice of optimizer.
  }
  \item{logLik}{the log likelihood of the model at
    optimizer termination.
  }
  \item{Hessian}{if the model was fitted with \code{Hess = TRUE}, this
    is the Hessian matrix of the parameters at the optimum.
  }
  \item{scale}{\code{model.frame} for the scale model.
  }
  \item{location}{\code{model.frame} for the location model.
  }
  \item{nominal}{\code{model.frame} for the nominal model.
  }
  \item{edf}{the (effective) number of degrees of freedom used by the
    model.
  }
  \item{start}{the starting values.
  }
  \item{method}{character, the optimizer.
  }
  \item{y}{the response variable.
  }
  \item{lev}{the names of the levels of the response variable.
  }
  \item{nobs}{the (effective) number of observations, calculated
    as the sum of the weights.
  }
  \item{threshold}{character, the threshold function used in the model.
  }
  \item{estimLambda}{\code{1} if lambda is estimated in one of the
    flexible link functions and \code{0} otherwise.
  }
  \item{link}{character, the link function used in the model.
  }
  \item{call}{the matched call.
  }
  \item{contrasts}{contrasts applied to terms in location and scale
    models.
  }
  \item{na.action}{the function used to filter missing data.
  }
}
\author{Rune Haubo B Christensen}
\references{
  Agresti, A. (2002) \emph{Categorical Data.} Second edition.  Wiley.

}
\examples{
options(contrasts = c("contr.treatment", "contr.poly"))
data(soup)

## More manageable data set:
dat <- subset(soup, as.numeric(as.character(RESP)) <=  24)
dat$RESP <- dat$RESP[drop=TRUE]

m1 <- clmm(SURENESS ~ PROD, random = RESP, data = dat, link="probit",
           Hess = TRUE, method="ucminf", threshold = "symmetric")

m1
summary(m1)
logLik(m1)
vcov(m1)
extractAIC(m1)
anova(m1, update(m1, location = SURENESS ~ 1, Hess = FALSE))
anova(m1, update(m1, random = NULL))

## Use adaptive Gauss-Hermite quadrature rather than the Laplace
## approximation:
update(m1, Hess = FALSE, nAGQ = 3)

## Use standard Gauss-Hermite quadrature:
update(m1, Hess = FALSE, nAGQ = -7)

##################################################################
## Binomial example with data from the lme4-package:
data(cbpp, package = "lme4")
cbpp2 <- rbind(cbpp[,-(2:3)], cbpp[,-(2:3)])
cbpp2 <- within(cbpp2, {
    incidence <- as.factor(rep(0:1, each=nrow(cbpp)))
    freq <- with(cbpp, c(incidence, size - incidence))
})

## Fit with Laplace approximation:
fm1 <- clmm(incidence ~ period, random = herd, weights = freq,
            data = cbpp2, Hess = 1)
summary(fm1)

## Fit with the adaptive Gauss-Hermite quadrature approximation:
fm2 <- clmm(incidence ~ period, random = herd, weights = freq,
            data = cbpp2, Hess = 1, nAGQ = 7)
summary(fm2)

}
\keyword{models}
