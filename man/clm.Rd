\name{clm}
\alias{clm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Cumulative Link Models
%%  ~~function to do ... ~~
}
\description{
  Fits cumulative link models (CLMs) such as the propotional odds
  model. The model allows for various link functions and structured
  thresholds that restricts the thresholds or cut-points to be e.g.,
  equidistant or symmetrically arranged around the central
  threshold(s). Nominal effects (partial proportional odds with the
  logit link) are also allowed.
  A modified Newton algorithm is used to optimize the likelihood function.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
clm(formula, scale, nominal, data, weights, start, subset, doFit = TRUE,
  na.action, contrasts, model = TRUE, control=list(),
  link = c("logit", "probit", "cloglog", "loglog", "cauchit"),
  threshold = c("flexible", "symmetric", "equidistant"), ...) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
    a formula expression as for regression models, of the form
    \code{response ~ predictors}. The response should be a factor
    (preferably an ordered factor), which will be interpreted as an
    ordinal response with levels ordered as in the factor.
    The model must have an intercept: attempts to remove one will
    lead to a warning and will be ignored. An offset may be used. See the
    documentation of \code{\link{formula}} for other details.
  }
  \item{scale}{
    an optional formula expression, of the form
    \code{ ~ predictors}, i.e. with an empty left hand side.
    An offset may be used. Variables included here will have
    multiplicative effects and can be interpreted as effects on the
    scale (or dispersion) of a latent distribution.
  }
  \item{nominal}{
    an optional formula of the form \code{ ~ predictors}, i.e. with an
    empty left hand side. The effects of the predictors in this formula
    are assumed to nominal - this corresponds to the so-called partial
    proportional odds (with the logit link).
}
  \item{data}{
    an optional data frame in which to interpret the variables occurring
    in the formulas.
  }
  \item{weights}{
    optional case weights in fitting. Defaults to 1.
  }
  \item{start}{
    initial values for the parameters in the format
    \code{c(alpha, beta, zeta)}, where \code{alpha} are the threshold
    parameters (adjusted for nominal effects), \code{beta} are the
    regression parameters and \code{zeta} are the scale parameters. 
  }
  \item{subset}{
    expression saying which subset of the rows of the data should  be used
    in the fit. All observations are included by default.
  }
  \item{doFit}{
    logical for whether the model should be fit or the model
    environment should be returned.
  }
  \item{na.action}{
    a function to filter missing data. Applies to terms in all three
    formulae. 
  }
  \item{contrasts}{
    a list of contrasts to be used for some or all of
    the factors appearing as variables in the model formula.
  }
  \item{model}{
    logical for whether the model frame should be part of the returned
    object.
  }
  \item{control}{
    a list of control parameters passed on to
    \code{\link{clm.control}}. 
  }
  \item{link}{
    link function, i.e., the type of location-scale distribution
    assumed for the latent distribution. The default \code{"logit"} link
    gives the proportional odds model.
  }
  \item{threshold}{
    specifies a potential structure for the thresholds
    (cut-points). \code{"flexible"} provides the standard unstructured
    thresholds, \code{"symmetric"} restricts the distance between the
    thresholds to be symmetric around the central one or two thresholds
    for odd or equal numbers or thresholds respectively, and
    \code{"equidistant"} restricts the distance between consecutive
    thresholds to be of the same size.
  }
  \item{\dots}{
    additional arguments are passed on to \code{\link{clm.control}}.
  }
}
\details{
  This is a new (as of August 2011) improved implementation of CLMs. The
  old implementation is available in \code{\link{clm2}}, but will
  probably be removed at some point.
  
  There are methods for the standard model-fitting functions, including
  \code{\link{summary}},
  \code{\link{anova}},
  \code{\link{model.frame}},
  \code{\link{model.matrix}},
  \code{\link{drop1}},
  \code{\link{dropterm}},
  \code{\link{step}},
  \code{\link{stepAIC}},
  \code{\link{extractAIC}},
  \code{\link{AIC}},
  \code{\link{coef}},
  \code{\link{nobs}},
  \code{\link{profile}},
  \code{\link{confint}},
  \code{\link{vcov}} and
  \code{\link[=slice.clm]{slice}}.
  %% \code{slice}.
  
}
\value{
  If \code{doFit = FALSE} the result is an environment
  representing the model ready to be optimized.
  If \code{doFit = TRUE} the result is an
  object of class \code{"clm2"} with the components listed below.
  
  Note that some components are only present if \code{scale} and
  \code{nominal} are used.

  \item{alpha}{a vector of threshold parameters.
  }
  \item{beta}{a vector of regression parameters.
  }
  \item{zeta}{a vector of scale regression parameters.
  }
  \item{coefficients}{a vector of coefficients of the form
    \code{c(alpha, beta, zeta)}
  }
  \item{aliased}{list of length 3 or less with components \code{alpha},
    \code{beta} and \code{zeta} each being logical vectors containing
    alias information for the parameters of the same names. 
  }
  \item{start}{the parameter values at which optimizations has started. 
  }
  \item{fitted.values}{the fitted probabilities.
  }
  \item{logLik}{the value of the log-likelihood at the estimated
    optimum. 
  }
  \item{edf}{the estimated degrees of freedom, i.e., the number of
    parameters in the model fit.
  }
  \item{nobs}{the number of observations counted as \code{sum(weights)}.
  }
  \item{n}{the number of observations counted as \code{nrow(X)}, where
    \code{X} is the design matrix.
  }
  \item{gradient}{a vector of gradients for the coefficients at the
    estimated optimum.
  }
  \item{maxGradient}{the maximum absolute gradient, i.e.,
    \code{max(abs(gradient))}. 
  }
  \item{Hessian}{the Hessian matrix for the parameters at the estimated
    optimum. 
  }
  \item{convergence}{convergence code where 0 indicates successful
    convergence, 1 indicates the iteration limit was reached before
    convergence and 2 indicates the step factor was reduced below minimum
    before convergence was reached. 
  }
  \item{link}{character, the link function used.
  }
  \item{threshold}{character, the threshold structure used.
  }
  \item{call}{the mathed call.
  }
  \item{contrasts}{(where relevant) the contrasts used for the \code{formula} part of the
    model. 
  }
  \item{S.contrasts}{(where relevant) the contrasts used for the \code{scale} part of the
    model. 
  }
  \item{nom.contrasts}{(where relevant) the contrasts used for the \code{nominal} part of the
    model.
  }
  \item{terms}{the terms object for the \code{formula} part.
  }
  \item{S.terms}{(where relevant) the terms object for the \code{scale} part.
  }
  \item{nom.terms}{(where relevant) the terms object for the \code{nominal} part.
  }
  \item{xlevels}{(where relevant) a record of the levels of the factors used in fitting
    for the \code{formula} part.
  }
  \item{S.xlevels}{(where relevant) a record of the levels of the factors used in fitting
    for the \code{scale} part.
  }
  \item{nom.xlevels}{(where relevant) a record of the levels of the factors used in fitting
    for the \code{nominal} part.
  }
  \item{tJac}{the transpose of the Jacobian for the threshold structure.
  }
  \item{y.levels}{the levels of the response variable.
  }
  \item{info}{a table of basic model information for printing.
  }
  \item{model}{if requested (the default), the
    \code{\link{model.frame}} containing variables from \code{formula},
    \code{scale} and \code{nominal} parts. 
  }
  \item{na.action}{(where relevant) information returned by
    \code{\link{model.frame}} on the special handling of \code{NA}s. 
  }
  
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\author{
  Rune Haubo B Christensen
}
\examples{

data(wine)
fm1 <- clm(rating ~ temp * contact, data = wine)
fm1 ## print method
summary(fm1)
fm2 <- update(fm1, ~.-temp:contact)
anova(fm1, fm2)

dropterm(fm1, test = "Chi")
drop1(fm1, test = "Chi")
add1(fm1, ~.+judge, test = "Chi")
addterm(fm1, ~.+judge, test = "Chi")

fm2 <- step(fm1)
summary(fm2)
fm3 <- stepAIC(fm1)
summary(fm3)

coef(fm1)
vcov(fm1)
AIC(fm1)
extractAIC(fm1)
logLik(fm1)
fitted(fm1)

confint(fm1) ## type = "profile"
confint(fm1, type = "Wald")
pr1 <- profile(fm1)
confint(pr1)

## plotting the profiles:
par(mfrow = c(2, 2))
plot(pr1, root = TRUE) ## check for linearity 
par(mfrow = c(2, 2))
plot(pr1)
par(mfrow = c(2, 2))
plot(pr1, approx = TRUE)
par(mfrow = c(2, 2))
plot(pr1, Log = TRUE)
par(mfrow = c(2, 2))
plot(pr1, Log = TRUE, relative = FALSE)

## other link functions:
fm4.lgt <- update(fm1, link = "logit") ## default
fm4.prt <- update(fm1, link = "probit")
fm4.ll <- update(fm1, link = "loglog")
fm4.cll <- update(fm1, link = "cloglog")
fm4.cct <- update(fm1, link = "cauchit")
anova(fm4.lgt, fm4.prt, fm4.ll, fm4.cll, fm4.cct)

## structured thresholds:
fm5 <- update(fm1, threshold = "symmetric")
fm6 <- update(fm1, threshold = "equidistant")
anova(fm1, fm5, fm6)

## the slice methods:
slice.fm1 <- slice(fm1)
par(mfrow = c(3, 3))
plot(slice.fm1)
## see more at '?slice.clm'

## Example from MASS::polr:
data(housing, package = "MASS")
fm1 <- clm(Sat ~ Infl + Type + Cont, weights = Freq, data = housing)
summary(fm1)

## Another example:
data(soup, package = "ordinal")

fm.soup <- clm(SURENESS ~ PRODID, data = soup)
summary(fm.soup)

}
\keyword{models}

