PbPb <- read.data('PbPb3.csv',method='Pb-Pb',format=3)

settings('iratio','U238U235',137.88,0)
settings('lambda','U238',0.000154993,0.00000013)

isochron(PbPb,model=1)

isochron(PbPb,inverse=TRUE)

isochron(PbPb,exterr=TRUE,show.numbers=TRUE,xlim=c(0,0.03),ylim=c(0,0.75),
         alpha=0.01,sigdig=3,ellipse.fill=NA,ellipse.stroke="red")

age(PbPb,common.Pb=2)

settings('iratio','Pb207Pb204',6.25)
settings('iratio','Pb206Pb204',10)
age(PbPb,common.Pb=1)

