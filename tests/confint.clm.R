## get data:
library(ordinal)
data(soup, package = "ordinal")
tab26 <- with(soup, table("Product" = PROD, "Response" = SURENESS))
dat26 <- expand.grid(sureness = as.factor(1:6), prod = c("Ref", "Test"))
dat26$wghts <- c(t(tab26))

m1 <- clm(sureness ~ prod, scale = ~prod, data = dat26,
          weights = wghts, link = "logistic")
m2 <- update(m1, link = "probit")
m6 <- update(m1, link = "Aranda-Ordaz", lambda = 1e-3)
m7 <- update(m1, link = "Aranda-Ordaz")
m8 <- update(m1, link = "log-gamma", lambda = 1)
m9 <- update(m1, link = "log-gamma")

## profile
pr2 <- profile(m2, trace = 1)
confint(m2)
par(mfrow = c(2,2))
plot(pr2)

pr9 <- profile(m9, trace = 1)
par(mfrow = c(2,2))
plot(pr9)

summary(m7)
pr7 <- profile(m7, range = c(2, 7))
## Important warnings:
warnings()
## c(pr7)
par(mfrow = c(2,2))
plot(pr7)

pr6 <- profile(m6, range = c(1e-3, 9))
par(mfrow = c(2,2))
plot(pr6)

pr8 <- profile(m8)
par(mfrow = c(2,2))
plot(pr8)

par(mfrow = c(2,2))
plot(pr9, Log=TRUE, relative = TRUE)
par(mfrow = c(2,2))
plot(pr9, Log=TRUE, relative = TRUE, ylim = c(-4, 0))
par(mfrow = c(2,2))
plot(pr9, Log=TRUE, relative = FALSE)
par(mfrow = c(1,1))
plot(pr9, Log=TRUE, relative = FALSE, parm = 3)

## confint:
confint(pr2)
confint(pr9)
confint(pr6)
confint(pr8)


## Extend example from polr in package MASS:
## Fit model from polr example:
data(housing, package = "MASS")
fm1 <- clm(Sat ~ Infl + Type + Cont, weights = Freq,
           data = housing, method = "Newton")

pr1 <- profile(fm1, trace = 1)
confint(pr1)
par(mfrow=c(3,2))
plot(pr1)
