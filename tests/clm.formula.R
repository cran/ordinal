library(ordinal)
data(wine)

#################################
## Appropriate evaluation of formulas:

## These fail and give appropriate error messages:
##  fm1 <- clm(rating ~ contact, scale=temp, data=wine)
##  fm1 <- clm(rating ~ contact, scale=~Temp, data=wine) 
##  fm1 <- clm(rating ~ contact, scale="temp", data=wine)
##  sca <- "temp"
##  fm1 <- clm(rating ~ contact, scale=sca, data=wine)
##  sca <- as.formula(sca)
##  sca <- as.formula(temp)
##  sca <- with(wine, as.formula(temp))

## These all work as intended with no warnings or errors:
fm1 <- clm(rating ~ contact, scale="~temp", data=wine)
fm1 <- clm(rating ~ contact, scale=~temp, data=wine)
sca <- "~temp"
fm1 <- clm(rating ~ contact, scale=sca, data=wine)
sca <- as.formula("~temp")
fm1 <- clm(rating ~ contact, scale=sca, data=wine)
fm1 <- clm(rating ~ contact, scale=as.formula(~temp), data=wine)
fm1 <- clm(rating ~ contact, scale=as.formula("~temp"), data=wine)

#################################
## can evaluate of if 'formula' is character:
f <- "rating ~ contact + temp"
clm(f, data=wine)
clm(as.formula(f), data=wine)

#################################

### finding variables in the environment of the formula:
data(wine)
makeform <- function() {
  f1 <- as.formula(rating ~ temp + contact)
  rating <- wine$rating
  temp <- wine$temp
  contact <- wine$contact
  f1
}
## 'makeform' makes are formula object in the environment of the
## function makeform:
f1 <- makeform()
f1 # print
class(f1)
## If we give the data, we can evaluate the model:
fm1 <- clm(f1, data=wine)
## We can also evaluate the model because the data are available in
## the environment associated with the formula:
fm1 <- clm(f1)
## For instance, the 'rating' variable is not found in the Global
## environment; we have to evaluate the 'name' of 'rating' in the
## appropriate environment:
(try(rating, silent=TRUE))
eval(as.name("rating"), envir=environment(f1))
## If instead we generate the formula in the Global environment where
## the variables are not found, we cannot evaluate the model:
f2 <- as.formula(rating ~ temp + contact)
(try(fm2 <- clm(f2), silent=TRUE))
environment(f2) <- environment(f1)
fm2 <- clm(f2)
#################################

