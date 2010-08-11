\name{clm.control}
\alias{clm.control}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Set control parameters for cumulative link models}
\description{
  Set control parameters for cumulative link models
}
\usage{
clm.control(method = c("Newton", "model.frame", "ucminf", "nlminb",
   "optim"), ..., trace = 0, 
   maxIter = 100, gradTol = 1e-06, maxLineIter = 15)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{method}{\code{"Newton"} fits the model by maximum likelihood and
    \code{"model.frame"} cause \code{clm} to return the \code{model.frame}.
  }
  \item{trace}{numerical, if > 0 information is printed about and during
    the optimization process. Defaults to \code{0}.
  }
  \item{maxIter}{the maximum number of Newton-Raphson iterations.
    Defaults to \code{100}.
  }
  \item{gradTol}{the maximum absolute gradient. This is the termination
    criterion and defaults to \code{1e-6}.
  }
  \item{maxLineIter}{the maximum number of step halfings allowed if
    a Newton(-Raphson) step over shoots. Defaults to \code{10}.
  }
  \item{\dots}{control arguments parsed on to \code{\link{ucminf}},
    \code{\link{nlminb}} or \code{\link{optim}}.
  } 
}
\value{
  a list of control parameters.
}
\author{Rune Haubo B Christensen}
\seealso{
  \code{\link{clm}}
}
\keyword{models}
