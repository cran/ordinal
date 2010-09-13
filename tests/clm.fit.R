library(ordinal)
data(wine)

## clm.fit with nominal and scale effects:

## get simple model:
fm1 <- clm(rating ~ temp, scale=~temp, nominal=~ contact,
           data=wine, method="model")
str(fm1, give.attr=FALSE)

## construct some weights and offsets:
set.seed(1)
off1 <- runif(length(fm1$y))
set.seed(1)
off2 <- rnorm(length(fm1$y))
set.seed(1)
wet <- runif(length(fm1$y))

## Fit various models:
fit <- clm.fit(fm1$y, fm1$X, fm1$S, fm1$NOM, weights=wet)
Coef <- 
  c(-0.905224120279548, 1.31043498891987, 3.34235590523008,
    4.52389661722693,      -3.03954652971192, -1.56922389038976,
    -1.75662549320839, -1.16845464236365,      2.52988580848393,
    -0.0261457032829033) 
all.equal(fit$par, Coef)
str(fit)

fit <- clm.fit(fm1$y, fm1$X, fm1$S, fm1$NOM, offset=off1)
str(fit)

fit <- clm.fit(fm1$y, fm1$X, fm1$S, fm1$NOM, offset=off1,
               S.offset=off2) 
str(fit)

fit <- clm.fit(fm1$y, fm1$X, fm1$S)
str(fit)

fit <- clm.fit(fm1$y, fm1$X)
str(fit)

fit <- clm.fit(fm1$y)
str(fit)

## Remember: compare with corresponding .Rout file
