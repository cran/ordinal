### library(ucminf)
### library(numDeriv)
### source("./ordinalDevel/R/clm.R")
### source("./ordinalDevel/R/clmmUnitUni.R")

library(ordinal)
options(contrasts = c("contr.treatment", "contr.poly"))
data(soup, package="ordinal")
## More manageable data set:
(tab26 <- with(soup, table("Product" = PROD, "Response" = SURENESS)))
dimnames(tab26)[[2]] <- c("Sure", "Not Sure", "Guess", "Guess", "Not Sure", "Sure")
dat26 <- expand.grid(sureness = as.factor(1:6), prod = c("Ref", "Test"))
dat26$wghts <- c(t(tab26))
m1 <- clm(sureness ~ prod, scale = ~prod, data = dat26,
          weights = wghts, link = "logistic")
## print, summary, vcov, logLik, AIC:
m1
summary(m1)
vcov(m1)

logLik(m1)
ll.m1 <- structure(-2687.74456343981, df = 7, class = "logLik")
stopifnot(all.equal(logLik(m1), ll.m1))

AIC(m1)

coef(m1)
cm1 <- c(-1.49125702755587, -0.45218462707814, -0.107208315524318, 0.163365282774162,
         0.88291347877514, 1.29587762626394, 0.147986162902775)
stopifnot(all.equal(as.vector(coef(m1)), cm1))

coef(summary(m1))
csm1 <- structure(c(-1.49125702755587, -0.45218462707814, -0.107208315524318,
0.163365282774162, 0.88291347877514, 1.29587762626394, 0.147986162902775,
0.0921506468161812, 0.0718240681909781, 0.069954084652323, 0.0702546879687391,
0.0795708692869622, 0.119032405993894, 0.065104213008022, -16.1828167145758,
-6.2957256316336, -1.53255261729392, 2.32532927691394, 11.0959385851501,
10.8867632762999, 2.27306584421104, 6.66732036748908e-59, 3.05965144996025e-10,
0.125386123756898, 0.0200543599621069, 1.31274723412040e-28,
1.33293711602276e-27, 0.0230222123418036), .Dim = c(7L, 4L), .Dimnames = list(
    c("1|2", "2|3", "3|4", "4|5", "5|6", "prodTest", "prodTest"
    ), c("Estimate", "Std. Error", "z value", "Pr(>|z|)")))
stopifnot(all.equal(coef(summary(m1)), csm1))

## link functions:
m2 <- update(m1, link = "probit")
m3 <- update(m1, link = "cloglog")

m4 <- update(m1, link = "loglog")
m5 <- update(m1, link = "cauchit", start = coef(m1))
m6 <- update(m1, link = "Aranda-Ordaz", lambda = 1)
m7 <- update(m1, link = "Aranda-Ordaz")
m8 <- update(m1, link = "log-gamma", lambda = 1)
m9 <- update(m1, link = "log-gamma")
## nominal effects:
mN1 <- clm(sureness ~ 1, nominal = ~ prod, data = dat26,
weights = wghts, link = "logistic")
anova(m1, mN1)
## optimizer / method:
update(m1, scale = ~ 1, method = "Newton")
update(m1, scale = ~ 1, method = "nlminb")
update(m1, scale = ~ 1, method = "optim")
update(m1, scale = ~ 1, method = "model.frame")
update(m1, location = ~.-prod, scale = ~ 1,
       nominal = ~ prod, method = "model.frame")
## threshold functions
mT1 <- update(m1, threshold = "symmetric")
mT2 <- update(m1, threshold = "equidistant")
anova(m1, mT1, mT2)
## Extend example from polr in package MASS:
## Fit model from polr example:
data(housing, package = "MASS")
fm1 <- clm(Sat ~ Infl + Type + Cont, weights = Freq, data = housing)
fm1
summary(fm1)
## With probit link:
summary(update(fm1, link = "probit"))
## Allow scale to depend on Cont-variable
summary(fm2 <- update(fm1, scale =~ Cont))
summary(fm3 <- update(fm1, location =~.-Cont, nominal =~ Cont))
anova(fm1, fm2, fm3)
## which seems to improve the fit

#################################
## Use of formula-objects in location, scale and nominal:
## Bug-report from Llu�s Marco Almagro <lluis.marco@upc.edu>
## 5 May 2010 17:58
f <- formula(sureness ~ prod)
fs <- formula( ~ prod)
m2 <- clm(f, scale = fs, data = dat26, weights = wghts)
summary(m2)

#################################
## Other ways to construct formulas:
y <- factor(sample(1:4,20,replace=TRUE))
x <- rnorm(20)
data <- data.frame(y=y,x=x)
rm(x, y)
fit <- clm(data$y ~ data$x)
fit <- clm(data[,1] ~ data[,2])
## This currently fails, but the union of variable names in location
## and scale could possibly be tested for:
## fit <- clm(data$y ~ data$x, ~data$x)

#################################
## Better handling of ill-defined variance-covariance matrix of the
## parameters in summary methods for clm and clmm objects:
dat26.2 <- data.frame(sureness = as.factor(1:12),
                      prod = rep(c("One", "Two", "Three"),each=4))
fm1 <- clm(sureness ~ prod, ~prod, data = dat26.2)
fm1
summary(fm1)
summary(fm1, corr = 1)

#################################
## Missing values:
## Bug-report from Jonathan Williams
## <Jonathan.Williams@dpag.ox.ac.uk>, 18 March 2010 12:42
data(soup, package = "ordinal")
soup$SURENESS[10] <- NA
c1a <- clm(ordered(SURENESS)~PROD, data=soup); summary(c1a)
c2a <- clm(ordered(SURENESS)~PROD, scale = ~PROD, data=soup); summary(c2a)
c3a <- clm(ordered(SURENESS)~1, scale = ~PROD, data=soup); summary(c3a)
data(soup, package = "ordinal")
soup$PROD[1] <- NA
c1a <- clm(ordered(SURENESS)~PROD, data=soup); summary(c1a)
c2a <- clm(ordered(SURENESS)~PROD, scale = ~PROD, data=soup); summary(c2a)
c3a <- clm(ordered(SURENESS)~1, scale = ~PROD, data=soup); summary(c3a)
soup$SURENESS[10] <- NA
c1a <- clm(ordered(SURENESS)~PROD, data=soup); summary(c1a)
c2a <- clm(ordered(SURENESS)~PROD, scale = ~PROD, data=soup); summary(c2a)
c3a <- clm(ordered(SURENESS)~1, scale = ~PROD, data=soup); summary(c3a)

## na.actions:
c4a <- clm(ordered(SURENESS)~PROD, scale = ~PROD, data=soup,
           na.action=na.omit)
summary(c4a)

tC1 <- try(clm(ordered(SURENESS)~PROD, scale = ~PROD, data=soup,
               na.action=na.fail), silent = TRUE)
stopifnot(class(tC1) == "try-error")

c4a <- clm(ordered(SURENESS)~PROD, scale = ~PROD, data=soup,
           na.action=na.exclude)
summary(c4a)

tC2 <- try(clm(ordered(SURENESS)~PROD, scale = ~PROD, data=soup,
               na.action=na.pass), silent = TRUE)
stopifnot(class(tC2) == "try-error")

## Subset:
c4a <- clm(ordered(SURENESS)~PROD, scale = ~PROD, data=soup,
           subset = 1:100)
c4a <- clm(ordered(SURENESS)~1, scale = ~PROD, data=soup,
           subset = 1:100)
c4a <- clm(ordered(SURENESS)~PROD, data=soup,
           subset = 1:100)
c4a <- clm(ordered(SURENESS)~1, data=soup,
           subset = 1:100)

## Offset:
data(soup, package = "ordinal")
set.seed(290980)
offs <- runif(nrow(soup))
c4a <- clm(ordered(SURENESS)~PROD + offset(offs),
           scale = ~PROD, data=soup, subset = 1:100)
summary(c4a)

c4a <- clm(ordered(SURENESS)~PROD + offset(offs),
           scale = ~PROD + offset(offs), data=soup, subset = 1:100)
summary(c4a)

c4a <- clm(ordered(SURENESS)~PROD,
           scale = ~PROD + offset(offs), data=soup, subset = 1:100)
summary(c4a)

## data as matrix:
dat26M <- as.matrix(dat26)
m1 <- clm(sureness ~ prod, scale = ~prod, data = dat26,
          weights = wghts, link = "logistic")
summary(m1)

## data in enclosing environment:
attach(soup)
m1 <- clm(SURENESS ~ PROD, scale = ~PROD)
summary(m1)
detach(soup)
