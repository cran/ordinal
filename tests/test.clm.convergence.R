library(ordinal)


## Testing that errors in chol() are caught soon enough:
cy <- with(wine, which(temp == "cold" & contact == "yes"))
wine2 <- subset(wine, subset=(!1:nrow(wine) %in% cy))
wine2[c(9, 15, 46), "rating"] <- NA
fm1 <- clm(rating ~ temp, scale=~contact, nominal=~contact,
           data=wine2)
fm1 <- try(clm(rating ~ temp, scale=~contact, nominal=~contact,
               data=wine2, control=list(gradTol=1e-12)), silent=TRUE)
fm2 <- try(clm(rating ~ temp, scale=~contact, nominal=~contact,
               data=wine2, control=list(gradTol=1e-15)), silent=TRUE)
## These gave errors in version 2014.11-12.
stopifnot(!inherits(fm1, "try-error"))
stopifnot(!inherits(fm2, "try-error"))
summary(fm1)
summary(fm2)
