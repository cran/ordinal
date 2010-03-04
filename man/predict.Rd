\name{predict.clm}
\alias{predict.clm}
\title{Predict Method for CLM fits}
\description{
  Obtains predictions from a cumulative link model.
}
\usage{
\method{predict}{clm}(object, newdata, ...)
}
\arguments{
  \item{object}{a fitted object of class inheriting from \code{"clm"}.}
  \item{newdata}{optionally, a data frame in which to look for variables
    with which to predict.}
  \item{\dots}{further arguments passed to or from other methods.}
}
\value{
  A vector of predicted probabilities.
}
\seealso{
    \code{\link[ordinal]{clm}}.
}
\examples{
options(contrasts = c("contr.treatment", "contr.poly"))
data(soup)

## More manageable data set for less voluminous printing:
(tab26 <- with(soup, table("Product" = PROD, "Response" = SURENESS)))
dimnames(tab26)[[2]] <- c("Sure", "Not Sure", "Guess", "Guess", "Not Sure", "Sure")
dat26 <- expand.grid(sureness = as.factor(1:6), prod = c("Ref", "Test"))
dat26$wghts <- c(t(tab26))
dat26

m1 <- clm(sureness ~ prod, scale = ~prod, data = dat26,
          weights = wghts, link = "logistic")
predict(m1)

mN1 <-  clm(sureness ~ 1, nominal = ~prod, data = dat26,
            weights = wghts)
predict(mN1)

predict(update(m1, scale = ~.-prod))

## Fit model from polr example:
data(housing, package = "MASS")
fm1 <- clm(Sat ~ Infl + Type + Cont, weights = Freq, data = housing)
predict(fm1)
}
\keyword{models}
