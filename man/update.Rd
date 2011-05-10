\name{update.clm}
\alias{update.clm}
\alias{update.clmm}
\title{Update method for cumulative link models}
\description{
  Update method for cumulative link models fitted with \code{clm}.
  This makes it possible to use e.g.
  \code{update(obj, location = ~ . - var1, scale = ~ . + var2)}
}
\usage{
\method{update}{clm}(object, formula., location, scale, nominal,...,
        evaluate = TRUE)
\method{update}{clmm}(object, formula., location, scale, nominal,...,
        evaluate = TRUE)
}
\arguments{
  \item{object}{a \code{\link{clm}} object.
  }
  \item{formula.}{not used---unfortunately this argument is part of the
    default method.
  }
  \item{location}{an optional new formula for the location; see
    \code{\link{update.formula}} for details.
  }
  \item{scale}{an optional new formula for the scale; see
    \code{\link{update.formula}} for details.
  }
  \item{nominal}{an optional new formula for nominal effects; see
    \code{\link{update.formula}} for details.
  }
  \item{\dots}{additional arguments to the call, or arguments with
    changed values.
  }
  \item{evaluate}{if true evaluate the new call else return the call.
  }
}
\value{
  If \code{evaluate = TRUE} the fitted object is returned,
  otherwise the updated call.
}
\author{Rune Haubo B Christensen}
\examples{
options(contrasts = c("contr.treatment", "contr.poly"))
data(soup)
m1 <-  clm(SURENESS ~ PROD, scale = ~PROD, data = soup,
           link = "logistic")

m2 <- update(m1, link = "probit")
m3 <- update(m1, link = "cloglog")
m4 <- update(m1, link = "loglog")
anova(m1, update(m1, scale = ~.-PROD))
mT1 <- update(m1, threshold = "symmetric")

## Fit model from polr example:
data(housing, package = "MASS")
fm1 <- clm(Sat ~ Infl + Type + Cont, weights = Freq, data = housing)
update(fm1, location = ~ . - Type, scale =~ Cont)

}
\keyword{models}
