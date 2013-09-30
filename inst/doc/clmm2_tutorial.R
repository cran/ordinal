### R code from vignette source 'clmm2_tutorial.Rnw'

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
### code chunk number 2: clmm2_tutorial.Rnw:151-154
###################################################
data(wine)
head(wine)
str(wine)


###################################################
### code chunk number 3: clmm2_tutorial.Rnw:175-189
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
### code chunk number 4: clmm2_tutorial.Rnw:216-218
###################################################
fm1 <- clmm2(rating ~ temp + contact, random=judge, data=wine)
fm1


###################################################
### code chunk number 5: clmm2_tutorial.Rnw:225-228
###################################################
fm2 <- clmm2(rating ~ temp + contact, random=judge, data=wine,
            Hess=TRUE, nAGQ=10)
summary(fm2)


###################################################
### code chunk number 6: clmm2_tutorial.Rnw:264-265
###################################################
exp(coef(fm2)[5])


###################################################
### code chunk number 7: clmm2_tutorial.Rnw:273-275
###################################################
fm3 <- clmm2(rating ~ temp, random=judge, data=wine, nAGQ=10)
anova(fm3, fm2)


###################################################
### code chunk number 8: clmm2_tutorial.Rnw:281-283
###################################################
fm4 <- clm2(rating ~ temp + contact, data=wine)
anova(fm4, fm2)


###################################################
### code chunk number 9: clmm2_tutorial.Rnw:294-296
###################################################
pr2 <- profile(fm2, range=c(.1, 4), nSteps=30, trace=0)
confint(pr2)


###################################################
### code chunk number 10: profilePlot
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(pr2)


###################################################
### code chunk number 11: profileFig
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(pr2)


###################################################
### code chunk number 12: ranefPlot
###################################################
getOption("SweaveHooks")[["fig"]]()
ci <- fm2$ranef + qnorm(0.975) * sqrt(fm2$condVar) %o% c(-1, 1)
ord.re <- order(fm2$ranef)
ci <- ci[order(fm2$ranef),]
plot(1:9, fm2$ranef[ord.re], axes=FALSE, ylim=range(ci),
     xlab="Judge", ylab="Judge effect")
axis(1, at=1:9, labels = ord.re)
axis(2)
for(i in 1:9) segments(i, ci[i,1], i, ci[i, 2])
abline(h = 0, lty=2)


###################################################
### code chunk number 13: clmm2_tutorial.Rnw:347-348
###################################################
getOption("SweaveHooks")[["fig"]]()
ci <- fm2$ranef + qnorm(0.975) * sqrt(fm2$condVar) %o% c(-1, 1)
ord.re <- order(fm2$ranef)
ci <- ci[order(fm2$ranef),]
plot(1:9, fm2$ranef[ord.re], axes=FALSE, ylim=range(ci),
     xlab="Judge", ylab="Judge effect")
axis(1, at=1:9, labels = ord.re)
axis(2)
for(i in 1:9) segments(i, ci[i,1], i, ci[i, 2])
abline(h = 0, lty=2)


###################################################
### code chunk number 14: clmm2_tutorial.Rnw:360-361
###################################################
head(cbind(wine, fitted(fm2)))


###################################################
### code chunk number 15: clmm2_tutorial.Rnw:366-367
###################################################
head(cbind(wine, pred=predict(fm2, newdata=wine)))


###################################################
### code chunk number 16: clmm2_tutorial.Rnw:385-387
###################################################
plogis(fm2$Theta[3] - fm2$beta[2]) - 
  plogis(fm2$Theta[2] - fm2$beta[2])


###################################################
### code chunk number 17: clmm2_tutorial.Rnw:395-396
###################################################
qnorm(0.95) * c(-1, 1) * fm2$stDev


###################################################
### code chunk number 18: clmm2_tutorial.Rnw:401-409
###################################################
pred <-
  function(eta, theta, cat = 1:(length(theta)+1), inv.link = plogis)
{
  Theta <- c(-1e3, theta, 1e3)
  sapply(cat, function(j)
         inv.link(Theta[j+1] - eta) - inv.link(Theta[j] - eta) )
}
pred(qnorm(0.05) * fm2$stDev, fm2$Theta)


###################################################
### code chunk number 19: clmm2_tutorial.Rnw:415-433
###################################################
mat <- expand.grid(judge = qnorm(0.95) * c(-1, 0, 1) * fm2$stDev,
                   contact = c(0, fm2$beta[2]),
                   temp = c(0, fm2$beta[1]))
pred.mat <- pred(eta=rowSums(mat), theta=fm2$Theta)
lab <- paste("contact=", rep(levels(wine$contact), 2), ", ",
             "temp=", rep(levels(wine$temp), each=2), sep="")
par(mfrow=c(2, 2))
for(k in c(1, 4, 7, 10)) {
  plot(1:5, pred.mat[k,], lty=2, type = "l", ylim=c(0,1),
       xlab="Bitterness rating scale", axes=FALSE,
       ylab="Probability", main=lab[ceiling(k/3)], las=1)
  axis(1); axis(2)
  lines(1:5, pred.mat[k+1, ], lty=1)
  lines(1:5, pred.mat[k+2, ], lty=3)
  legend("topright",
         c("avg. judge", "5th %-tile judge", "95th %-tile judge"),
         lty=1:3, bty="n")
}


###################################################
### code chunk number 20: clmm2_tutorial.Rnw:438-448
###################################################
getOption("SweaveHooks")[["fig"]]()
k <- 1
plot(1:5, pred.mat[k,], lty=2, type = "l", ylim=c(0,1),
     xlab="Bitterness rating scale", axes=FALSE,
     ylab="Probability", main=lab[ceiling(k/3)], las=1)
axis(1); axis(2)
lines(1:5, pred.mat[k+1, ], lty=1)
lines(1:5, pred.mat[k+2, ], lty=3)
legend("topright",
       c("avg. judge", "5th %-tile judge", "95th %-tile judge"),
       lty=1:3, bty="n")


###################################################
### code chunk number 21: clmm2_tutorial.Rnw:450-460
###################################################
getOption("SweaveHooks")[["fig"]]()
k <- 4
plot(1:5, pred.mat[k,], lty=2, type = "l", ylim=c(0,1),
     xlab="Bitterness rating scale", axes=FALSE,
     ylab="Probability", main=lab[ceiling(k/3)], las=1)
axis(1); axis(2)
lines(1:5, pred.mat[k+1, ], lty=1)
lines(1:5, pred.mat[k+2, ], lty=3)
legend("topright",
       c("avg. judge", "5th %-tile judge", "95th %-tile judge"),
       lty=1:3, bty="n")


###################################################
### code chunk number 22: clmm2_tutorial.Rnw:462-472
###################################################
getOption("SweaveHooks")[["fig"]]()
k <- 7
plot(1:5, pred.mat[k,], lty=2, type = "l", ylim=c(0,1),
     xlab="Bitterness rating scale", axes=FALSE,
     ylab="Probability", main=lab[ceiling(k/3)], las=1)
axis(1); axis(2)
lines(1:5, pred.mat[k+1, ], lty=1)
lines(1:5, pred.mat[k+2, ], lty=3)
legend("topright",
       c("avg. judge", "5th %-tile judge", "95th %-tile judge"),
       lty=1:3, bty="n")


###################################################
### code chunk number 23: clmm2_tutorial.Rnw:474-484
###################################################
getOption("SweaveHooks")[["fig"]]()
k <- 10
plot(1:5, pred.mat[k,], lty=2, type = "l", ylim=c(0,1),
     xlab="Bitterness rating scale", axes=FALSE,
     ylab="Probability", main=lab[ceiling(k/3)], las=1)
axis(1); axis(2)
lines(1:5, pred.mat[k+1, ], lty=1)
lines(1:5, pred.mat[k+2, ], lty=3)
legend("topright",
       c("avg. judge", "5th %-tile judge", "95th %-tile judge"),
       lty=1:3, bty="n")


###################################################
### code chunk number 24: clmm2_tutorial.Rnw:494-495
###################################################
exp(2*qnorm(0.95) * fm2$stDev)


###################################################
### code chunk number 25: clmm2_tutorial.Rnw:501-502
###################################################
exp(2*qnorm(0.75) * fm2$stDev)


###################################################
### code chunk number 26: misc (eval = FALSE)
###################################################
## 


