library(IsoplotR)

PbPb <- read.data('PbPb3.csv',method='Pb-Pb',format=3)

settings('iratio','U238U235',137.88,0)
settings('lambda','U238',0.000154993,0.00000013)

isochron(PbPb,model=1)

isochron(PbPb,inverse=FALSE,growth=TRUE)

PbPb2 <- read.data('Kamber.csv',method='Pb-Pb',format=1)
isochron(PbPb2,growth=TRUE)

settings('alpha',0.01)
isochron(PbPb,exterr=TRUE,show.numbers=TRUE,
         xlim=c(0,0.03),ylim=c(0,0.75),
         oerr=3,sigdig=3,ellipse.fill=NA,
         ellipse.stroke="red")
settings('alpha',0.05)

age(PbPb,isochron=FALSE,common.Pb=2,i=2)

settings('iratio','Pb207Pb204',6.25)
settings('iratio','Pb206Pb204',10)
age(PbPb,isochron=FALSE,common.Pb=1,i=2)

age(PbPb,isochron=FALSE,common.Pb=2,i=2,exterr=TRUE)

age(PbPb,isochron=FALSE,common.Pb=2,i=2,exterr=TRUE,projerr=TRUE)

ThPb <- read.data('ThPb2.csv',method='Th-Pb',format=2)

settings('iratio','Pb208Pb206',30,0)
settings('lambda','Th232',0.00005,0.000003)

isochron(ThPb,inverse=FALSE,exterr=FALSE)

isochron(ThPb,inverse=TRUE,show.numbers=TRUE,
         ellipse.fill=NA,ellipse.stroke=rgb(0,0,1,0.5))

age(ThPb,isochron=FALSE,i2i=TRUE,i=18)

settings('iratio','Pb208Pb204',38.82)
age(ThPb,isochron=FALSE,i2i=FALSE,i=18)

radialplot(PbPb,show.numbers=TRUE,exterr=TRUE,
           transformation='linear',
           k=1,from=4555,to=4575,z0=4565)

settings('iratio','Pb206Pb204',9.15)
settings('iratio','Pb207Pb204',10.23)
weightedmean(PbPb,common.Pb=1,random.effects=TRUE,ranked=TRUE)

kde(PbPb,kde.col='orange',rug=FALSE,show.hist=FALSE,
    common.Pb=2,bw=20,from=4500,to=4650)

cad(ThPb,verticals=TRUE,pch='x')

radialplot(ThPb,i2i=TRUE)
