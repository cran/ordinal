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

## anova
anova(m1, update(m1, scale = ~.-prod))
mN1 <- clm(sureness ~ 1, nominal = ~prod, data = dat26,
           link = "logistic")
anova(m1, mN1)
anova(m1, update(m1, scale = ~.-prod), mN1)

## dropterm
dropterm(m1, test = "Chi")
mB1 <- clm(SURENESS ~ PROD + GENDER + SOUPTYPE,
           scale = ~ COLD, data = soup, link = "probit",
           Hess = FALSE)
dropterm(mB1, test = "Chi")       # or
dropterm(mB1, which = "location", test = "Chi")
dropterm(mB1, which = "scale", test = "Chi")

## addterm
addterm(mB1, scope = ~.^2, test = "Chi", which = "location")
addterm(mB1, scope = ~ . + AGEGROUP + SOUPFREQ,
        test = "Chi", which = "location")
addterm(mB1, scope = ~ . + GENDER + SOUPTYPE,
        test = "Chi", which = "scale")

## Fit model from polr example:
data(housing, package = "MASS")
fm1 <- clm(Sat ~ Infl + Type + Cont, weights = Freq, data = housing)
addterm(fm1, ~ Infl + Type + Cont, test= "Chisq", which = "scale")
dropterm(fm1, test = "Chisq")
fm2 <- update(fm1, scale =~ Cont)
fm3 <- update(fm1, location =~.-Cont, nominal =~ Cont)
anova(fm1, fm2, fm3)

