\name{clmm.control}
\alias{clmm.control}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Set control parameters for cumulative link mixed models  
}
\description{
  Set control parameters for cumulative link mixed models
}
\usage{
clmm.control(method = c("ucminf", "model.frame"), ..., trace = 0,
maxIter = 50, gradTol = 1e-4, maxLineIter = 50, innerCtrl =
c("warnOnly", "noWarn", "giveError")) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{method}{
    the optimizer used to maximize the marginal likelihood function.
  }
  \item{\dots}{control arguments passed on to the optimizer; see
    \code{\link[ucminf]{ucminf}} for details.
    \code{ucminf} for details.
  }
  \item{trace}{numerical, if > 0 information is printed about and during
    the outer optimization process, if < 0 information is also printed
    about the inner optimization process. Defaults to \code{0}.
  }
  \item{maxIter}{the maximum number of Newton updates of the inner
    optimization. \code{50}.
  }
  \item{gradTol}{the maximum absolute gradient of the inner
    optimization.
  }
  \item{maxLineIter}{the maximum number of step halfings allowed if
    a Newton(-Raphson) step over shoots during the inner optimization.
  }
  \item{innerCtrl}{the use of warnings/errors if the inner optimization
    fails to converge.
  }
}
\value{
  a list of control parameters
}
\author{
  Rune Haubo B Christensen
}
\seealso{
  \code{\link{clmm}}
  }
\keyword{models}

