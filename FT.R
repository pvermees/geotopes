setwd('~/Documents/geotopes/R/')
library(IsoplotR)

FT1 <- read.data('FT1.csv',method='fissiontracks',format=1)
FT2 <- read.data('FT2.csv',method='fissiontracks',format=2)
FT3 <- read.data('FT3.csv',method='fissiontracks',format=3)

dens <- settings('mindens','zircon')

settings('mindens','apatite',3.23)
settings('etchfact','zircon',0.99)
settings('tracklength','apatite',15)

radialplot(FT1,transformation='arcsin')

radialplot(FT2,transformation='sqrt')

radialplot(FT2,exterr=TRUE)

radialplot(FT1,k='auto')

z <- set.zeta(x=FT1,tst=c(100,1))

age(FT1,exterr=TRUE)

weightedmean(FT3)

kde(FT2,log=TRUE)

cad(FT1)
