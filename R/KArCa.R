library(IsoplotR)

ArAr <- read.data('ArAr1.csv',method='Ar-Ar',format=1)

settings('iratio','Ar40Ar36',295.5,0.5)

KCa <- read.data('KCa3.csv',method='K-Ca',format=3)

settings('lambda','K40',5.543e-4,0)

isochron(ArAr,model=3)

isochron(KCa,model=2)

isochron(ArAr,inverse=FALSE,exterr=FALSE)

isochron(KCa,inverse=TRUE)

isochron(KCa,inverse=TRUE,exterr=FALSE,show.numbers=TRUE,
         xlim=c(0,2),ylim=c(0,0.02),ellipse.fill=rgb(0.5,1,0.5,0.2))

age(ArAr,i2i=TRUE)

age(KCa,i2i=TRUE)

settings('iratio','Ca40Ca44',50,0)
age(KCa,i2i=FALSE)

agespectrum(ArAr,plateau=TRUE,i2i=FALSE)

agespectrum(ArAr,plateau=TRUE,exterr=TRUE)

agespectrum(ArAr,levels=ArAr$x[,'Ar40Ar36'],plateau=FALSE,
            plateau.col=topo.colors(n=100,alpha=0.5),
            exterr=TRUE,i2i=FALSE,random.effects=FALSE,
            clabel=expression(''^40*'Ar/'^36*'Ar'))

kde(KCa,log=TRUE,from=700,to=1000,bw=0.02,
    binwidth=0.02,adaptive=FALSE,i2i=TRUE)

radialplot(ArAr,from=60,to=64,z0=61,levels=ArAr$x[,'Ar40Ar36'],
           clabel=expression(''^40*'Ar/'^36*'Ar'))

settings('iratio','Ar40Ar36',295.5,0.5)
weightedmean(ArAr,i2i=FALSE,random.effects=TRUE,ranked=TRUE)

cad(KCa,verticals=FALSE)
