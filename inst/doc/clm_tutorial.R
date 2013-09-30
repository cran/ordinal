### R code from vignette source 'clm_tutorial.Rnw'

###################################################
### code chunk number 1: Initialize
###################################################

## Load common packages, functions and set settings:
library(ordinal)
library(xtable)
## 
RUN <- FALSE    #redo computations and write .RData files
## Change options:
op <- options() ## To be able to reset settings
options("digits" = 7)
options(help_type = "html")
## options("width" = 75)
options("SweaveHooks" = list(fig=function()
        par(mar=c(4,4,.5,0)+.5)))
options(continue=" ")



###################################################
### code chunk number 2: clm_tutorial.Rnw:152-155
###################################################
data(wine)
head(wine)
str(wine)


###################################################
### code chunk number 3: clm_tutorial.Rnw:176-190
###################################################
data(wine)
temp.contact.bottle <- with(wine, temp:contact:bottle)[drop=TRUE]
tab <- xtabs(as.numeric(rating) ~ temp.contact.bottle + judge,
             data=wine) 
class(tab) <- "matrix"
attr(tab, "call") <- NULL
mat <- cbind(rep(c("cold", "warm"), each = 4),
             rep(rep(c("no", "yes"), each=2), 2),
             1:8, tab)
colnames(mat) <-
  c("Temperature", "Contact", "Bottle", 1:9)
xtab <- xtable(mat)
print(xtab, only.contents=TRUE, include.rownames=FALSE,
      sanitize.text.function = function(x) x)


###################################################
### code chunk number 4: clm_tutorial.Rnw:232-234
###################################################
fm1 <- clm(rating ~ temp + contact, data=wine)
fm1


###################################################
### code chunk number 5: clm_tutorial.Rnw:237-238
###################################################
summary(fm1)


###################################################
### code chunk number 6: clm_tutorial.Rnw:257-258
###################################################
clm.control()$gradTol


###################################################
### code chunk number 7: clm_tutorial.Rnw:277-278
###################################################
exp(coef(fm1)[5])


###################################################
### code chunk number 8: clm_tutorial.Rnw:286-288
###################################################
fm2 <- clm(rating ~ temp, data=wine)
anova(fm2, fm1)


###################################################
### code chunk number 9: clm_tutorial.Rnw:294-295
###################################################
drop1(fm1, test = "Chi")


###################################################
### code chunk number 10: clm_tutorial.Rnw:300-302
###################################################
fm0 <- clm(rating ~ 1, data=wine)
add1(fm0, scope = ~ temp + contact, test = "Chi")


###################################################
### code chunk number 11: clm_tutorial.Rnw:308-309
###################################################
confint(fm1)


###################################################
### code chunk number 12: clm_tutorial.Rnw:315-316
###################################################
confint(fm1, type="Wald")


###################################################
### code chunk number 13: clm_tutorial.Rnw:323-324
###################################################
fm.cll <- clm(rating ~ contact + temp, data=wine, link="cloglog")


###################################################
### code chunk number 14: clm_tutorial.Rnw:355-357
###################################################
fm.nom <- clm(rating ~ temp, nominal=~contact, data=wine)
summary(fm.nom)


###################################################
### code chunk number 15: clm_tutorial.Rnw:381-382
###################################################
anova(fm1, fm.nom)


###################################################
### code chunk number 16: clm_tutorial.Rnw:393-394
###################################################
fm.nom2 <- clm(rating ~ temp + contact, nominal=~contact, data=wine)


###################################################
### code chunk number 17: clm_tutorial.Rnw:397-398
###################################################
summary(fm.nom2)


###################################################
### code chunk number 18: clm_tutorial.Rnw:475-477
###################################################
fm.sca <- clm(rating ~ temp + contact, scale=~temp, data=wine)
summary(fm.sca)


###################################################
### code chunk number 19: clm_tutorial.Rnw:487-488
###################################################
exp(fm.sca$zeta)


###################################################
### code chunk number 20: clm_tutorial.Rnw:506-509
###################################################
fm.equi <- clm(rating ~ temp + contact, data=wine,
               threshold="equidistant") 
summary(fm.equi)


###################################################
### code chunk number 21: clm_tutorial.Rnw:516-517
###################################################
c(with(fm.equi, tJac %*% alpha))


###################################################
### code chunk number 22: clm_tutorial.Rnw:522-523
###################################################
mean(diff(fm1$alpha))


###################################################
### code chunk number 23: clm_tutorial.Rnw:529-530
###################################################
anova(fm1, fm.equi)


###################################################
### code chunk number 24: clm_tutorial.Rnw:545-546
###################################################
predict(fm1, type = "class")


###################################################
### code chunk number 25: clm_tutorial.Rnw:553-556
###################################################
newData <- expand.grid(temp=levels(wine$temp),
                       contact=levels(wine$contact))
cbind(newData, predict(fm1, newdata=newData)$fit)


###################################################
### code chunk number 26: clm_tutorial.Rnw:563-564
###################################################
head(do.call("cbind", predict(fm1, se.fit=TRUE, interval=TRUE)))


###################################################
### code chunk number 27: clm_tutorial.Rnw:575-577
###################################################
fm.nom2 <- clm(rating ~ contact, nominal=~temp, data=wine)
summary(fm.nom2)


###################################################
### code chunk number 28: clm_tutorial.Rnw:595-596
###################################################
anova(fm1, fm.nom2)


###################################################
### code chunk number 29: clm_tutorial.Rnw:604-607
###################################################
data(soup)
fm.soup <- clm(SURENESS ~ PRODID * DAY, data=soup)
summary(fm.soup)


###################################################
### code chunk number 30: clm_tutorial.Rnw:612-613
###################################################
with(soup, table(DAY, PRODID))


###################################################
### code chunk number 31: clm_tutorial.Rnw:621-624
###################################################
mm <- model.matrix( ~ PRODID * DAY, data=soup)
ncol(mm)
qr(mm, LAPACK = FALSE)$rank


###################################################
### code chunk number 32: clm_tutorial.Rnw:656-657
###################################################
convergence(fm1)


###################################################
### code chunk number 33: clm_tutorial.Rnw:671-674
###################################################
slice.fm1 <- slice(fm1, lambda = 5)
par(mfrow = c(2, 3))
plot(slice.fm1)


###################################################
### code chunk number 34: slice11
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(slice.fm1, parm = 1)


###################################################
### code chunk number 35: slice12
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(slice.fm1, parm = 2)


###################################################
### code chunk number 36: slice13
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(slice.fm1, parm = 3)


###################################################
### code chunk number 37: slice14
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(slice.fm1, parm = 4)


###################################################
### code chunk number 38: slice15
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(slice.fm1, parm = 5)


###################################################
### code chunk number 39: slice16
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(slice.fm1, parm = 6)


###################################################
### code chunk number 40: slice2
###################################################
slice2.fm1 <- slice(fm1, lambda = 1e-5)
par(mfrow = c(2, 3))
plot(slice2.fm1)


###################################################
### code chunk number 41: slice21
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(slice2.fm1, parm = 1)


###################################################
### code chunk number 42: slice22
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(slice2.fm1, parm = 2)


###################################################
### code chunk number 43: slice23
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(slice2.fm1, parm = 3)


###################################################
### code chunk number 44: slice24
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(slice2.fm1, parm = 4)


###################################################
### code chunk number 45: slice25
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(slice2.fm1, parm = 5)


###################################################
### code chunk number 46: slice26
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(slice2.fm1, parm = 6)


###################################################
### code chunk number 47: profileLikelihood
###################################################
pr1 <- profile(fm1, alpha=1e-4)
plot(pr1)


###################################################
### code chunk number 48: prof1
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(pr1, which.par=1)


###################################################
### code chunk number 49: prof2
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(pr1, which.par=2)


###################################################
### code chunk number 50: misc (eval = FALSE)
###################################################
## 


