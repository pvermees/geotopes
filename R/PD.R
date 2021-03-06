library(IsoplotR)

RbSr <- read.data('RbSr1.csv',method='Rb-Sr',format=1)
SmNd <- read.data('SmNd2.csv',method='Sm-Nd',format=2)
LuHf <- read.data('LuHf3.csv',method='Lu-Hf',format=3)
ReOs <- read.data('ReOs2.csv',method='Re-Os',format=2)

settings('lambda','Rb87',0.0000142,0)
settings('lambda','Sm147',0.000006540,0.000000049)
settings('lambda','Re187',0.000016689,0.000000038)
settings('lambda','Lu176',0.00001865,0.00000015)
settings('lambda',reset=TRUE)

isochron(RbSr,model=1)
isochron(SmNd,model=2)
isochron(LuHf,model=3)

isochron(ReOs,inverse=FALSE)
isochron(ReOs,inverse=TRUE)

age(RbSr,isochron=FALSE,exterr=TRUE,i2i=TRUE,projerr=TRUE,i=1)

radialplot(LuHf,show.numbers=TRUE,transformation='linear',
           cex=1.5,levels=LuHf$x[,'Luppm'],clabel='Lu [ppm]',
           bg=rev(topo.colors(n=100)),i2i=TRUE)

weightedmean(ReOs,detect.outliers=FALSE,
             ranked=TRUE,rect.col=NA,i2i=TRUE)

kde(RbSr,from=4494,to=4504,binwidth=1,adaptive=FALSE,i2i=TRUE)

settings('iratio','Nd143Nd144',0.508)
cad(SmNd,i2i=FALSE)
