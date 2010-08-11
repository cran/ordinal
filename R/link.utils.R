#################################
## pfun:

pgumbel <-
  function(q, location = 0, scale = 1, lower.tail = TRUE, max = TRUE)
### CDF for Gumbel max and min distributions 
### Currently only unit length location and scale are supported.
{
  if(max)  ## right skew, loglog link
    .C("pgumbel",
       q = as.double(q),
       length(q),
       as.double(location)[1],
       as.double(scale)[1],
       as.integer(lower.tail),
       NAOK = TRUE)$q
  else ## left skew, cloglog link
    .C("pgumbel2",
       q = as.double(q),
       length(q),
       as.double(location)[1],
       as.double(scale)[1],
       as.integer(lower.tail),
       NAOK = TRUE)$q
}

## pgumbel <- function(q, location = 0, scale = 1, lower.tail = TRUE)
## ### CDF for the gumbel distribution
## ### Currently only unit length location and scale are supported.
##     .C("pgumbel",
##        q = as.double(q),
##        length(q),
##        as.double(location)[1],
##        as.double(scale)[1],
##        as.integer(lower.tail),
##        NAOK = TRUE)$q
## 
## pgumbel2 <- function(q, location = 0, scale = 1, lower.tail = TRUE)
## ### CDF for the 'swapped' gumbel distribution
## ### Currently only unit length location and scale are supported.
##     .C("pgumbel2",
##        q = as.double(q),
##        length(q),
##        as.double(location)[1],
##        as.double(scale)[1],
##        as.integer(lower.tail),
##        NAOK = TRUE)$q

pgumbelR <- function(q, location = 0, scale = 1, lower.tail = TRUE)
### R equivalent of pgumbel()
{
    q <- (q - location)/scale
    p <- exp(-exp(-q))
    if (!lower.tail) 1 - p else p
}

pgumbel2R <- function(q, location = 0, scale = 1, lower.tail = TRUE)
{
    q <- (-q - location)/scale
    p <- exp(-exp(-q))
    if (!lower.tail) p else 1 - p
}

pAOR <- function(q, lambda, lower.tail = TRUE) {
    if(lambda < 1e-6)
        stop("'lambda' has to be positive. lambda = ", lambda, " was supplied")
    p <- 1 - (lambda * exp(q) + 1)^(-1/lambda)
    if(!lower.tail) 1 - p else p
}

pAO <- function(q, lambda, lower.tail = TRUE)
    .C("pAO",
       q = as.double(q),
       length(q),
       as.double(lambda[1]),
       as.integer(lower.tail),
       NAOK = TRUE)$q


plgammaR <- function(eta, lambda, lower.tail = TRUE) {
    q <- lambda
    v <- q^(-2) * exp(q * eta)
    if(q < 0)
        p <- 1 - pgamma(v, q^(-2))
    if(q > 0)
        p <- pgamma(v, q^(-2))
    if(isTRUE(all.equal(0, q, tolerance = 1e-6)))
        p <- pnorm(eta)
    if(!lower.tail) 1 - p else p
}

plgamma <- function(q, lambda, lower.tail = TRUE)
    .C("plgamma",
       q = as.double(q),
       length(q),
       as.double(lambda[1]),
       as.integer(lower.tail[1]),
       NAOK = TRUE)$q

#################################
## dfun:

dgumbel <-
  function(x, location = 0, scale = 1, log = FALSE, max = TRUE)
### PDF for the Gumbel max and mon distributions
{
  if(max)  ## right skew, loglog link
    .C("dgumbel",
       x = as.double(x),
       length(x),
       as.double(location)[1],
       as.double(scale)[1],
       as.integer(log),
       NAOK = TRUE)$x
  else ## left skew, cloglog link
    .C("dgumbel2",
       x = as.double(x),
       length(x),
       as.double(location)[1],
       as.double(scale)[1],
       as.integer(log),
       NAOK = TRUE)$x
} 

## Deprecated:
## dgumbel <- function(x, location = 0, scale = 1, log = FALSE)
## ### PDF for the gumbel distribution
## ### Currently only unit length location and scale are supported.
##     .C("dgumbel",
##        x = as.double(x),
##        length(x),
##        as.double(location)[1],
##        as.double(scale)[1],
##        as.integer(log),
##        NAOK = TRUE)$x
## 
## dgumbel2 <- function(x, location = 0, scale = 1, log = FALSE) {
## ### PDF for the 'swapped' gumbel distribution
## ### Currently only unit length location and scale are supported.
##   stopifnot(length(location) == 1 && ## test here?
##             length(scale) == 1 &&
##             length(log) == 1)
##   .C("dgumbel2",
##      x = as.double(x),
##      length(x),
##      as.double(location)[1],
##      as.double(scale)[1],
##      as.integer(log),
##      NAOK = TRUE)$x
## }

dgumbelR <- function(x, location = 0, scale = 1, log = FALSE)
### dgumbel in R
{
    q <- (x - location)/scale
    log.d <- -exp(-q) - q - log(scale)
    if (!log) exp(log.d) else log.d
}

dgumbel2R <- function(x, location = 0, scale = 1, log = FALSE)
{
    q <- (-x - location)/scale
    log.d <- -exp(-q) - q - log(scale)
    if (!log) exp(log.d) else log.d
}

dAOR <- function(eta, lambda, log = FALSE) {
### exp(eta) * (lambda * exp(eta) + 1)^(-1-1/lambda)
  stopifnot(length(lambda) == 1 &&
            length(log) == 1)
  if(lambda < 1e-6)
    stop("'lambda' has to be positive. lambda = ", lambda,
         " was supplied") 
  log.d <- eta - (1 + 1/lambda) * log(lambda * exp(eta) + 1)
  if(!log) exp(log.d) else log.d
}

dAO <- function(eta, lambda, log = FALSE) {
  stopifnot(length(lambda) == 1 &&
            length(log) == 1)
  .C("dAO",
     eta = as.double(eta),
     length(eta),
     as.double(lambda),
     as.integer(log),
     NAOK = TRUE)$eta
}

dlgammaR <- function(x, lambda, log = FALSE) {
    q <- lambda
    q.2 <- q^(-2)
    qx <- q * x
    log.d <- log(abs(q)) + q.2 * log(q.2) -
        lgamma(q.2) + q.2 * (qx - exp(qx))
    if (!log) exp(log.d) else log.d
}

dlgamma <- function(x, lambda, log = FALSE) {
  stopifnot(length(lambda) == 1 &&
            length(log) == 1)
  .C("dlgamma",
     x = as.double(x),
     length(x),
     as.double(lambda),
     as.integer(log),
     NAOK = TRUE)$x
}

#################################
## gfun:

ggumbel <- function(x, max = TRUE) {
### gradient of dgumbel(x) wrt. x
  if(max) ## right skew, loglog link
    .C("ggumbel",
       x = as.double(x),
       length(x),
       NAOK = TRUE)$x
  else ## left skew, cloglog link
    .C("ggumbel2",
       x = as.double(x),
       length(x),
       NAOK = TRUE)$x
}

## ggumbel <- function(x)
## ### gradient of dgumbel(x) wrt. x
##     .C("ggumbel",
##        x = as.double(x),
##        length(x),
##        NAOK = TRUE)$x
## 
## ggumbel2 <- function(x)
## ### gradient of dgumbel(x) wrt. x
##     .C("ggumbel2",
##        x = as.double(x),
##        length(x),
##        NAOK = TRUE)$x

ggumbelR <- function(x){
### ggumbel in R
  q <- exp(-x)
  ifelse(q == Inf, 0, {
    eq <- exp(-q)
    -eq*q + eq*q*q
  })
}

ggumbel2R <- function(x) -ggumbelR(-x)

glogis <- function(x)
### gradient of dlogis
    .C("glogis",
       x = as.double(x),
       length(x),
       NAOK = TRUE)$x

gnorm <- function(x)
### gradient of dnorm(x) wrt. x
    .C("gnorm",
       x = as.double(x),
       length(x),
       NAOK = TRUE)$x

gcauchy <- function(x)
### gradient of dcauchy(x) wrt. x
    .C("gcauchy",
       x = as.double(x),
       length(x),
       NAOK = TRUE)$x

glogisR <- function(x) {
### glogis in R
  res <- rep(0, length(x))
  isFinite <- !is.infinite(x)

  x <- x[isFinite]
  isNegative <- x < 0
  q <- exp(-abs(x))
  q <- 2*q^2*(1 + q)^-3 - q*(1 + q)^-2
  q[isNegative] <- -q[isNegative]
  res[isFinite] <- q
  res
}

gnormR <- function(x)
### gnorm in R
    -x * dnorm(x)

gcauchyR <- function(x)
### gcauchy(x) in R
    -2*x/pi*(1+x^2)^-2

gAOR <- function(eta, lambda) {
  stopifnot(length(lambda) == 1)
  lex <- lambda * exp(eta)
  dAO(eta, lambda) * (1 - (1 + 1/lambda) * lex/(1 + lex))
}

gAO <- function(eta, lambda) {
  stopifnot(length(lambda) == 1)
  .C("gAO",
     eta = as.double(eta),
     length(eta),
     as.double(lambda[1]),
     NAOK = TRUE)$eta
}

glgammaR <- function(x, lambda) {
  stopifnot(length(lambda) == 1)
  (1 - exp(lambda * x))/lambda * dlgamma(x, lambda)
}

glgammaR2 <- function(x, lambda) {
  stopifnot(length(lambda == 1))
  if(lambda == 0)
    return(gnorm(x))
  y <- dlgamma(x, lambda)
  y[!is.na(y) && y > 0] <- y * (1 - exp(lambda * x))
  return(y)
}

glgamma <- function(x, lambda) {
  stopifnot(length(lambda) == 1)
  .C("glgamma",
     x = as.double(x),
     length(x),
     as.double(lambda[1]),
     NAOK = TRUE)$x
}

##################################################################
## Random numbers:

rgumbel <- function(n, location = 0, scale = 1, max = TRUE) {
  if(max)
    location - scale * log(-log(runif(n)))
  else
    location + scale * log(-log(runif(n)))
}

##################################################################
## quantile functions:

qgumbel <-
  function(p, location = 0, scale = 1, lower.tail = TRUE, max = TRUE)
{ 
  if(max)  ## right skew, loglog link
    location - scale * log(-log(p))
  else ## left skew, cloglog link
    location + scale * log(-log(1 - p))
}


##################################################################
PFUN <- function(x, lambda = 1, link)
    .C("pfun",
       x = as.double(x),
       length(x),
       as.double(lambda),
       as.integer(link))$x

DFUN <- function(x, lambda = 1, link)
    .C("dfun",
       x = as.double(x),
       length(x),
       as.double(lambda),
       as.integer(link))$x

GFUN <- function(x, lambda = 1, link)
    .C("gfun",
       x = as.double(x),
       length(x),
       as.double(lambda),
       as.integer(link))$x
