library(IsoplotR)

ThU1 <- read.data('ThU1.csv',method='Th-U',format=1)
ThU2 <- read.data('ThU2.csv',method='Th-U',format=2)
ThU3 <- read.data('ThU3.csv',method='Th-U',format=3)
ThU4 <- read.data('ThU4.csv',method='Th-U',format=4)

l34 <- settings('lambda','U234')[1]

settings('lambda','Th230',log(2)/75)

oldpar <- par(mfrow=c(2,2))
for (i in 1:4){ isochron(ThU1,type=i) }
par(oldpar)

isochron(ThU3,type=2)

isochron(ThU1,type=3,model=3)

isochron(ThU1,type=3,model=1,y0option=1)

evolution(ThU2)

evolution(ThU2,xlim=c(0.2,1.2),ylim=c(0.5,1.5))

evolution(ThU2,transform=TRUE)

evolution(ThU2,transform=TRUE,xlim=c(100,250),ylim=c(1,1.2))

evolution(ThU1,Th0i=1)

ThU2b <- read.data('ThU2.csv',method='Th-U',format=2,
                   Th02U48=c(1,0.02,2,0.1,1,0.1,0,0,0))
evolution(ThU2b,Th0i=2)

ThU2c <- read.data('ThU2.csv',method='Th-U',
                   format=2,Th02i=c(1,0.1))
evolution(ThU2c,Th0i=3)

evolution(ThU3)

evolution(ThU3b,Th0i=1)

ThU3b <- read.data('ThU3.csv',method='Th-U',
                   format=3,U8Th2=2)
evolution(ThU3b,Th0i=2)

evolution(ThU1,isochron=TRUE,model=1,exterr=TRUE)

evolution(ThU1,oerr=4,levels=ThU1$x[,'U238Th232'],
          clabel='U/Th',ellipse.fill=
          cm.colors(n=10,alpha=0.5),ellipse.stroke='orange')

age(ThU1,Th0i=1)

age(ThU2b,Th0i=2)

age(ThU2c,Th0i=3)

age(ThU3,Th0i=1)

radialplot(ThU1,show.numbers=TRUE,
           col='blue',Th0i=0,omit=1)

weightedmean(ThU2c,Th0i=2,ranked=TRUE)

kde(ThU3,Th0i=1,log=TRUE,from=120,
    to=200,bw=0.05,binwidth=0.05)

cad(ThU3b,Th0i=0)
