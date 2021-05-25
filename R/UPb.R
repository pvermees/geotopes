setwd('~/Documents/IsoplotR/')
library(IsoplotR)

UPb1 <- read.data('UPb1.csv',method='U-Pb',format=1)
UPb2 <- read.data('UPb6.csv',method='U-Pb',format=6,ierr=3)
UPb3 <- read.data('UPb8.csv',method='U-Pb',format=8,ierr=4)

# use the Steiger and Jaeger (1977) ratio with zero uncertainty
settings('iratio','U238U235',137.88,0)
# use the Schoene et al. (2006) decay constant and uncertainty
settings('lambda','U238',0.000154993,0.00000013)

settings('iratio','Pb207Pb206',1.106)
radialplot(UPb1,common.Pb=1)

settings('iratio','Pb206Pb204',9.307)
settings('iratio','Pb207Pb204',10.294)
weightedmean(UPb2,common.Pb=1)

settings('iratio','Pb208Pb206',3.1671)
settings('iratio','Pb208Pb207',2.8634)
weightedmean(UPb3,common.Pb=1)

settings(reset=TRUE)

kde(UPb2,common.Pb=2)

cad(UPb1,common.Pb=3)

d <- diseq(U48=list(x=1.05,option=2))
UPb4 <- read.data('diseq.csv',method='U-Pb',format=2,d=d)
concordia(UPb4,type=2)

d <- diseq(RaU=list(x=1.2,option=1))
UPb4 <- read.data('diseq.csv',method='U-Pb',format=2,d=d)
concordia(UPb4,type=2)

d <- diseq(ThU=list(x=4,option=3))
UPb4 <- read.data('diseq.csv',method='U-Pb',format=8,d=d)
radialplot(UPb4)

d <- diseq(U48=list(x=0,option=1),ThU=list(x=2,option=1),
           RaU=list(x=2,option=1),PaU=list(x=2,option=1))
UPb4 <- read.data('diseq.csv',method='U-Pb',format=2,d=d)
concordia(UPb4,type=2,xlim=c(0,5000),ylim=c(0.047,0.057))

concordia(UPb1,type=1) # Wetherill
concordia(UPb1,type=2) # Tera-Wasserburg
concordia(UPb3,type=3) # U-Th-Pb

concordia(UPb1,tlim=c(300,3000))

concordia(UPb1,type=2,xlim=c(10,30),ylim=c(0.05,0.06))

concordia(UPb1,alpha=0.14)

ThU <- UPb3$x[,'Th232U238']

concordia(UPb3,levels=ThU,ellipse.fill=c('blue','red'))

concordia(UPb3,levels=ThU,
          ellipse.fill=c(rgb(1,0,0,0.5),rgb(0,1,0,0.5)))

concordia(UPb3,levels=ThU,ellipse.fill='white',
          ellipse.stroke=heat.colors(n=10,alpha=0.8))

concordia(UPb3,levels=ThU,ellipse.fill=NA,
          ellipse.stroke=rev(rainbow(n=10)))

concordia(UPb3,levels=ThU,clabel='Th/U')

concordia(UPb1,ticks=c(230,240,250))

concordia(UPb1,ticks=10)

concordia(UPb1,show.age=1,omit=10)

concordia(UPb1,show.age=1,hide=10,exterr=TRUE)

concordia(UPb2,type=2,show.age=2)

concordia(UPb2,type=1,show.age=3)

settings('iratio','Pb206Pb204',9.307)
settings('iratio','Pb207Pb204',10.294)
concordia(UPb2,type=2,show.age=3,anchor=1)
settings(reset=TRUE)

concordia(UPb2,type=2,show.age=3,anchor=c(2,300))

concordia(UPb2,show.age=4,sigdig=3)

oldpar <- par(mfrow=c(1,2))
isochron(UPb2,type=1) # 204Pb/206Pb vs. 238U/206Pb
isochron(UPb2,type=2) # 204Pb/207Pb vs. 235U/207Pb
par(oldpar)

oldpar <- par(mfrow=c(2,2))
isochron(UPb3,type=1) # 208Pbc/206Pb
isochron(UPb3,type=2) # 208Pbc/207Pb
isochron(UPb3,type=3) # 206Pbc/208Pb
isochron(UPb3,type=4) # 207Pbc/208Pb
par(oldpar)

oldpar <- par(mfrow=c(2,2))
isochron(UPb2,type=1,model=3)
isochron(UPb2,type=2,model=2)
isochron(UPb3,type=2,model=1)
isochron(UPb3,type=3,model=3)
par(oldpar)

ThU <- UPb3$x[,'Th232U238']
fit <- isochron(UPb3,type=3,model=1,xlim=c(0,2000),
                ylim=c(0,0.6),alpha=0.05,exterr=TRUE,
                sigdig=3,show.numbers=TRUE,levels=ThU,
                ellipse.fill=c('white','red'),
                ellipse.stroke='blue',clabel='Th/U')

radialplot(UPb1,type=2)

radialplot(UPb1,type=4,cutoff.76=1200)

radialplot(UPb1,cutoff.disc=
             discfilter(option=5,cutoff=c(-2,7)))

dscf <- cutoff.disc=discfilter(option=2,before=FALSE,
                               cutoff=c(-5,15))
radialplot(UPb2,common.Pb=2,cutoff.disc=dscf)

settings('iratio','Pb208Pb206',2.089)
settings('iratio','Pb208Pb207',2.519)
dscf <- cutoff.disc=discfilter(option=3,before=FALSE,
                               cutoff=c(0,1))
radialplot(UPb3,common.Pb=1,cutoff.disc=dscf)
settings(reset=TRUE)

dscf <- cutoff.disc=discfilter(option=4,before=TRUE,
                               cutoff=c(-1,5))
radialplot(UPb3,common.Pb=3,cutoff.disc=dscf)

oldpar <- par(cex=0.9)
radialplot(UPb3,show.numbers=TRUE,exterr=TRUE,
           from=5,z0=30,to=150,pch=21,alpha=0.05,
           bg=c('white','red'),cex=2,
           levels=UPb3$x[,'Th232U238'],clabel='Th/U')
par(oldpar)

weightedmean(UPb1)

d <- diseq(U48=list(x=0,option=1),ThU=list(x=2,option=1),
           RaU=list(x=2,option=1),PaU=list(x=2,option=1))
UPb4 <- read.data('diseq.csv',method='U-Pb',format=2,d=d)
weightedmean(UPb4,type=2,common.Pb=3,
             cutoff.disc=discfilter(option=5))

weightedmean(UPb3,random.effects=TRUE,detect.outliers=TRUE,
             outlier.col='red',ranked=TRUE,from=-10,to=100,
             alpha=0.02,sigdig=1,levels=UPb3$x[,'Th232U238'],
             rect.col=rev(topo.colors(n=100)),clabel='Th/U')

d <- discfilter(option=4,before=FALSE,cutoff=c(-1,5))
kde(UPb2,common.Pb=2,type=2,cutoff.disc=d)

cad(UPb1,common.Pb=3,type=4,cutoff.76=1100)

age(UPb2,common.Pb=3)

age(UPb1,exterr=TRUE)

age(UPb3,common.Pb=3,discordance=
      discfilter(before=TRUE,option=6))

tt <- age(UPb1,sigdig=4)