\name{clm.fit}
\alias{clm.fit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Fit Cumulative Link Models
%%  ~~function to do ... ~~
}
\description{
  A direct fitter of cumulative link models.
}
\usage{

clm.fit(y, X, S, N, weights = rep(1, nrow(X)),
  offset = rep(0, nrow(X)), S.offset = rep(0, nrow(X)),
  control = list(), start,
  link = c("logit", "probit", "cloglog", "loglog", "cauchit"),
  threshold = c("flexible", "symmetric", "symmetric2", "equidistant"))

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{the response variable; a factor, preferably and ordered
    factor.
  }
  \item{X, S, N}{optional design matrices for the regression parameters,
    scale parameters and nominal parameters respectively.
  }
  \item{weights}{optional case weights.
  }
  \item{offset}{an optional offset.
  }
  \item{S.offset}{an optional offset for the scale part of the model.
  }
  \item{control}{a list of control parameters, optionally a call to
    \code{\link{clm.control}}.
  }
  \item{start}{an optional list of starting values of the form
    \code{c(alpha, beta, zeta)} for the thresholds and nominal effects
    (\code{alpha}), regression parameters (\code{beta}) and scale
    parameters (\code{zeta}).
  }
  \item{link}{the link function.
  }
  \item{threshold}{the threshold structure, see further at
    \code{\link{clm}}.
  }
}
\details{
  This function does almost the same thing that \code{\link{clm}} does:
  it fits a cumulative link model. The main differences are that
  \code{clm.fit} does not setup design matrices from formulae and only
  does minimal post processing after parameter estimation.

  Compared to \code{\link{clm}}, \code{clm.fit} does little to warn the
  user of any problems with data or model. However, \code{clm.fit} will
  attempt to identify column rank defecient designs. Any unidentified
  parameters are indicated in the \code{aliased} component of the fit.

  In contrast to \code{\link{clm}}, \code{clm.fit} allows non-positive
  weights.
}

\value{
  A list with the following components:
  \item{par}{the maximum likelihood estimate of the model parameters on
    the form \code{c(alpha, beta, zeta)} for the thresholds and nominal
    effects (\code{alpha}), regression parameters (\code{beta}) and
    scale parameters (\code{zeta}).
  }
  \item{gradient}{gradient of the log-likelihood function for the
    parameters at the optimum.
  }
  \item{Hessian}{the Hessian of the parameters at the optimum.
  }
  \item{logLik}{the value of the log-likelihood function at the
    optimum.
  }
  \item{convergence}{a convergence code; 1: successful convergence, 2:
    iteration limit reached, 3: step factor reduced below minimum.
  }
  \item{message}{convergence message.
  }
  \item{maxGradient}{the maximum absolute gradient:
    \code{max(abs(gradient))}.
  }
  \item{niter}{depends on the optimizer. For the default Newton
    optimizer it is a vector of length 2: the number of iteration and
    the number of step halvings respectively.
  }
  \item{fitted}{the fitted values (probabilities).
  }
  \item{coef.names}{a list of vectors of names of the coefficients as
    inferred from the column names of the design matrices.
  }
  \item{aliased}{a list of logical vectors; \code{TRUE} is a parameter
    is not identified due to column rank deficiency in the design
    matrices and \code{FALSE} otherwise.
  }
}
%% \references{ bla
%%   %% ~put references to the literature/web site here ~
%% }
\author{
  Rune Haubo B Christensen
}
%% \note{ bla
%% %%  ~~further notes~~
%% }
%%
%% %% ~Make other sections like Warning with \section{Warning }{....} ~
%%
\seealso{
  \code{\link{clm}}
}
\examples{

## A simple example:
fm1 <- clm(rating ~ contact + temp, data=wine)
summary(fm1)
## get the model frame containing y and X:
mf1 <- update(fm1, method="model.frame")
clm.fit(mf1$y, mf1$X)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
