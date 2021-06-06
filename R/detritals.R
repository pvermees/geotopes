library(IsoplotR)

DZ <- read.data('DZ.csv',method='detritals')

kde(DZ,hide=c('T8','T13'))

kde(DZ,from=0,to=3000)

kde(DZ,log=TRUE,from=200,to=3000)

kde(DZ,samebandwidth=FALSE)

kde(DZ,log=TRUE,bw=0.05,binwidth=0.1)

kde(DZ,rug=TRUE)

kde(DZ,normalise=FALSE)

kde(DZ,adaptive=FALSE)

cad(DZ,col='rainbow')

mds(DZ,classical=TRUE)

mds(DZ,classical=FALSE,shepard=TRUE)

mds(DZ,nnlines=TRUE)

mds(DZ,pch=22,cex=0.9,pos=1,bg='yellow',col='red')
