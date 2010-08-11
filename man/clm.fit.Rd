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
clm.fit(y, X, weights = rep(1, nrow(X)), offset = rep(0, nrow(X)),
control = list(), start,
link = c("logit", "probit", "cloglog", "loglog", "cauchit"),
threshold = c("flexible", "symmetric", "equidistant"))

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{the response variable; a factor, preferably and ordered
    factor. 
    %%     ~~Describe \code{y} here~~
}
\item{X}{the design matrix for the regression parameters.
%%     ~~Describe \code{X} here~~
}
  \item{weights}{optional case weights.
%%     ~~Describe \code{weights} here~~
}
  \item{offset}{an optional offset.
%%     ~~Describe \code{offset} here~~
}
  \item{control}{a list of control parameters, optionally a call to
    \code{\link{clm.control}}. 
%%     ~~Describe \code{control} here~~
}
  \item{start}{an optional list of starting values of the form
    \code{c(alpha, beta)} for the thresholds (\code{alpha}) and
    regression parameters (\code{beta}). 
%%     ~~Describe \code{start} here~~
}
  \item{link}{the link function.
%%     ~~Describe \code{link} here~~
}
  \item{threshold}{the threshold structure, see further at
    \code{\link{clm}}.
%%     ~~Describe \code{threshold} here~~
}
}
%% \details{ bla
%% %%  ~~ If necessary, more details than the description above ~~
%% }
\value{
  A list with the following components:
  \item{par}{the maximum likelihood estimate of the model parameters on 
    the form \code{c(alpha, beta)} for the thresholds (\code{alpha}) and  
    regression parameters (\code{beta}). 
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
  \item{niter}{a vector of length 2: the number of iteration and the
    number of step halvings respectively.
  }
  \item{fitted}{the fitted values (probabilities).
  }
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
%% \references{ bla
%%   %% ~put references to the literature/web site here ~
%% }
\author{
  Rune Haubo B Christensen
%%  ~~who you are~~
}
%% \note{ bla
%% %%  ~~further notes~~
%% }
%% 
%% %% ~Make other sections like Warning with \section{Warning }{....} ~
%% 
%% \seealso{ bla
%% %% ~~objects to See Also as \code{\link{help}}, ~~~
%% }
\examples{

## A simple example:
data(wine)
fm1 <- clm(rating ~ contact + temp, data=wine)
summary(fm1)
## get the model frame containing y and X:
mf1 <- update(fm1, method="model.frame")
clm.fit(mf1$y, mf1$X)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
