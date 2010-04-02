library(ordinal)
options(contrasts = c("contr.treatment", "contr.poly"))
data(soup, package = "ordinal")

## More manageable data set:
dat <- subset(soup, as.numeric(as.character(RESP)) <=  24)
dat$RESP <- dat$RESP[drop=TRUE]

m1 <- clmm(SURENESS ~ PROD, random = RESP, data = dat, link="probit",
           Hess = TRUE, method="ucminf", threshold = "symmetric")

m1
summary(m1)
logLik(m1)
vcov(m1)
extractAIC(m1)
anova(m1, update(m1, location = SURENESS ~ 1, Hess = FALSE))
anova(m1, update(m1, random = NULL))

## Use adaptive Gauss-Hermite quadrature rather than the Laplace
## approximation:
update(m1, Hess = FALSE, nAGQ = 3)


##################################################################
## Binomial example with data from the lme4-package:
data(cbpp, package = "lme4")
cbpp2 <- rbind(cbpp[,-(2:3)], cbpp[,-(2:3)])
cbpp2 <- within(cbpp2, {
    incidence <- as.factor(rep(0:1, each=nrow(cbpp)))
    freq <- with(cbpp, c(incidence, size - incidence))
})

## Fit with Laplace approximation:
fm1 <- clmm(incidence ~ period, random = herd, weights = freq,
            data = cbpp2, Hess = 1)
summary(fm1)

## Testing link functions:
(m1 <- update(fm1, link = "probit", Hess = 0))
update(fm1, link = "loglog", Hess = 0)
update(fm1, link = "cloglog", Hess = 0)
strt <- c(coef(m1)[1:4], log(coef(m1)[5]))
## FAILS:
fm1.try <-
    try(update(fm1, link = "cauchit", Hess = 0, start = strt,
               method = "nlminb",
               control = list(maxIter = 1000, maxLineIter = 1000,
               trace = 1)), silent = TRUE)
class(fm1.try) == "try-error" ## TRUE
(fm1 <- update(fm1, link = "Aranda-Ordaz", Hess = 0))
(fm1 <- update(fm1, link = "log-gamma", Hess = 0))
(fm1 <- update(fm1, link = "Aranda-Ordaz", Hess = 0, lambda = 1))
(fm1 <- update(fm1, link = "log-gamma", Hess = 0, lambda = 1))


## Fit with the adaptive Gauss-Hermite quadrature approximation:
fm2 <- clmm(incidence ~ period, random = herd, weights = freq,
            data = cbpp2, Hess = 1, nAGQ = 7)
summary(fm2)

#################################
### Setting sdFixed:
fm1 <- clmm(incidence ~ period, random = herd, weights = freq,
            data = cbpp2, Hess = 0, link = "Aranda-Ordaz",
            sdFixed = .5)
fm1


##################################################################
## Testing subset:
(m1 <- clmm(SURENESS ~ PROD, random = RESP, data = dat, link="probit",
           Hess = FALSE, method="ucminf", threshold = "symmetric",
           subset = RESP != "13"))


## Testing missing values:
dat$RESP[78] <- NA
(m1 <- clmm(SURENESS ~ PROD, random = RESP, data = dat, link="probit",
            threshold = "symmetric", subset = RESP != "13"))
stopifnot(nrow(m1$location) == 189)  # OK

dat$PROD[13] <- NA
(m1 <- clmm(SURENESS ~ PROD, random = RESP, data = dat, link="probit",
            threshold = "symmetric", subset = RESP != "13"))
stopifnot(nrow(m1$location) == 188) # OK

dat$SURENESS[7] <- NA
(m1 <- clmm(SURENESS ~ PROD, random = RESP, data = dat, link="probit",
            threshold = "symmetric", subset = RESP != "13"))
stopifnot(nrow(m1$location) == 187) # OK

## factor in call:
Resp <- as.numeric(as.character(dat$RESP))
m2 <- clmm(SURENESS ~ PROD, random = factor(Resp), data = dat,
           link = "probit", Hess = 0, threshold = "symmetric")
m2

## control arguments:
(fm1 <- clmm(incidence ~ period, random = herd, weights = freq,
            data = cbpp2, method = "nlminb", trace = 1))
(fm1 <- clmm(incidence ~ period, random = herd, weights = freq,
            data = cbpp2, method = "model.frame"))
(fm1 <- clmm(incidence ~ period, random = herd, weights = freq,
            data = cbpp2, method = "ucminf", trace = 1))
(fm1 <- clmm(incidence ~ period, random = herd, weights = freq,
             data = cbpp2,
             control = clmm.control(grtol = 1e-8, trace = 1),
             method = "ucminf"))

#################################
## Profile:
fm1 <- clmm(incidence ~ period, random = herd, weights = freq,
            data = cbpp2, Hess = 1)

summary(fm1)
pr.fm1 <- profile(fm1, alpha = .01, nSteps = 5)
attr(pr.fm1, "WaldCI")
c(pr.fm1)
confint(pr.fm1, level = .95)
logLik(fm1)
plot(pr.fm1, Log = FALSE, relative = TRUE)
par(mfrow = c(2,2))
plot(pr.fm1)
plot(pr.fm1, Log=TRUE, relative = TRUE)
plot(pr.fm1, Log=TRUE, relative = FALSE)


pr2.fm1 <- profile(fm1, nSteps = 5, range = c(.2, 1.3))
pr2.fm1
par(mfrow = c(1,1))
plot(pr2.fm1, Log = FALSE, relative = TRUE)

#################################

