library(ordinal)
options(contrasts = c("contr.treatment", "contr.poly"))
data(soup, package="ordinal")
## More manageable data set:
(tab26 <- with(soup, table("Product" = PROD, "Response" = SURENESS)))
dimnames(tab26)[[2]] <- c("Sure", "Not Sure", "Guess", "Guess", "Not Sure", "Sure")
dat26 <- expand.grid(sureness = as.factor(1:6), prod = c("Ref", "Test"))
dat26$wghts <- c(t(tab26))

## Tests based on bug caught by Simon Blomberg, 13. Dec. 2010.
mN1 <- clm(sureness ~ 1, nominal = ~ prod, data = dat26,
           weights = wghts, link = "logistic")
stopifnot(isTRUE(all.equal(predict(mN1),
                           predict(mN1, newdata = dat26))))

mN1 <- clm(sureness ~ prod, data = dat26,
           weights = wghts, link = "logistic", threshold = "symm")
stopifnot(isTRUE(all.equal(predict(mN1),
                           predict(mN1, newdata = dat26))))

mN1 <- clm(sureness ~ 1, nominal = ~ prod, data = dat26,
           weights = wghts, link = "logistic", threshold = "symm")
stopifnot(isTRUE(all.equal(predict(mN1),
                           predict(mN1, newdata = dat26))))

mN1 <- clm(sureness ~ prod, scale = ~ prod, data = dat26,
           weights = wghts, link = "logistic", threshold = "symm")
stopifnot(isTRUE(all.equal(predict(mN1),
                           predict(mN1, newdata = dat26))))

## Running clmm-instances:
dat <- subset(soup, as.numeric(as.character(RESP)) <=  24)
dat$RESP <- dat$RESP[drop=TRUE]
m1 <- clmm(SURENESS ~ 1, nominal = ~ PROD, random = RESP,
           data = dat, link="probit",
           Hess = TRUE, method="ucminf")
p1 <- predict(m1, newdata = dat)
p1 <- predict(m1)

m1 <- clmm(SURENESS ~ 1, nominal =~ PROD, random = RESP, data = dat,
           link="probit",
           Hess = TRUE, method="ucminf", threshold = "symmetric")
p1 <- predict(m1, newdata = dat)
p1 <- predict(m1)

m1 <- clmm(SURENESS ~ PROD, random = RESP, data = dat, link="probit",
           Hess = TRUE, method="ucminf", threshold = "symmetric")
p1 <- predict(m1, newdata = dat)
p1 <- predict(m1)

