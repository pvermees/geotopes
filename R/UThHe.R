library(IsoplotR)

UThHe <- read.data('UThHe.csv',method='U-Th-He')
UThSmHe <- read.data('UThSmHe.csv',method='U-Th-He')

settings('lambda','Sm147',0.000006540,0.000000049)

oldpar <- par(mfrow=c(1,2))
helioplot(UThHe,logratio=TRUE)
helioplot(UThHe,logratio=FALSE)
par(oldpar)

helioplot(UThHe,show.numbers=TRUE)

helioplot(UThHe,show.central.comp=TRUE)

helioplot(UThHe,model=3)

helioplot(UThHe,oerr=1,sigdig=3)

helioplot(UThHe,logratio=TRUE,xlim=c(3,5),ylim=c(3,5))

helioplot(UThHe,logratio=FALSE,fact=c(150,2,3))

helioplot(UThSmHe,levels=UThSmHe[,'Sm'],clabel='Sm',
          ellipse.fill=c('#00005080','#00FF0080'),
          ellipse.stroke='grey30')

isochron(UThHe,show.ellipses = 1)

age(UThHe)

radialplot(UThSmHe,transformation='sqrt',k='min',
           clabel='Sm',levels=UThSmHe[,'Sm'],
           bg=c('blue','white','red'))

weightedmean(UThSmHe,random.effects=TRUE,ranked=TRUE,
             clabel='Sm',levels=UThSmHe[,'Sm'],
             rect.col=c('blue','white','red'))

kde(UThSmHe,log=TRUE,from=5,to=12)

cad(UThSmHe)
