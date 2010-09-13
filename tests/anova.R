library(ordinal)
data(wine)

fm1 <- clm(rating ~ temp, data=wine)
fmm1 <- clmm(rating ~ temp + (1|judge), data=wine)

## These now give identical printed results:
## Previously the printed model names were messed up when anova.clmm
## were called. 
anova(fm1, fmm1)
anova(fmm1, fm1)
