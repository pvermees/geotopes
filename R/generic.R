library(IsoplotR)

data1 <- read.data('MountTom.csv',method='other')
data2 <- read.data('regression.csv',method='other')
data3 <- read.data('LudwigMean.csv',method='other')
data4 <- read.data('LudwigSpectrum.csv',method='other')
data5 <- read.data('LudwigKDE.csv',method='other')

fn <- system.file('MountTom.csv',package='IsoplotR')
data1 <- read.data(fn,method='other')

radialplot(data1)

radialplot(data1,show.numbers=TRUE)

radialplot(data1,transformation='sqrt')

radialplot(data1,k='auto')

radialplot(data1,k='min',from=10,to=200,z0=18.27,
           alpha=0.01,sigdig=3)

radialplot(data1,pch=23)

radialplot(data1,alpha=0.01)

radialplot(data1,sigdig=3)

radialplot(data1,bg='blue')

relerr <- data1[,2]/data1[,1]
radialplot(data1,levels=relerr,bg=c('white','red'))

radialplot(data1,levels=relerr,bg=c('white','red'),clabel='s[t]/t')

oldpar <- par(cex=0.5)
radialplot(data1,cex=2)
par(oldpar) # restore the old cex value

isochron(data2)

isochron(data2,model=3)

isochron(data2,show.numbers=TRUE)

isochron(data2,model=2,xlim=c(0,300),ylim=c(0,1500))

oldpar <- par(cex=1.1)
isochron(data2,alpha=0.01,sigdig=4)
par(oldpar)

isochron(data2,levels=data2[,5],clabel='rho',
         ellipse.fill=topo.colors(n=100),ellipse.stroke='red')

weightedmean(data3,random.effects=FALSE)

weightedmean(data3,detect.outliers=FALSE,outlier.col='red')

weightedmean(data3,ranked=TRUE)

weightedmean(data3,from=220,to=280,alpha=0.1,sigdig=1,
             rect.col='#FF000080',outlier.col=NA)

agespectrum(data4)

agespectrum(data4,plateau=FALSE)

agespectrum(data4,random.effects=TRUE)

agespectrum(data4,plateau.col='blue',non.plateau.col='yellow')

kde(data5)

kde(data5,from=230,to=280)

kde(data5,bw=.3,binwidth=2)

kde(data5,log=TRUE,bw=0.012,binwidth=0.008,from=230,to=280)

kde(data5,show.hist=FALSE)

kde(data5,adaptive=FALSE)

kde(data5,rug=FALSE)

cad(data5)

cad(data5,verticals=FALSE)

cad(data5,verticals=FALSE,pch=16)
