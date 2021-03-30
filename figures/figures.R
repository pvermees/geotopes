library(IsoplotR)

setwd('/home/pvermees/Documents/geotopes/figures/')

pars <- function(mar=c(2.5,2.3,0.5,0.25),mgp=c(1.5,0.5,0),...){
    par(list(mar=mar,mgp=mgp,...))
}

cairo <- function(file,width,height,family="serif",pointsize=13,...){
    cairo_pdf(file=file,width=width,height=height,
              family=family,pointsize=pointsize,...)
}

cairo(file='spurious.pdf',width=8,height=2)
pars(mfrow=c(1,4),mar=c(2.5,2.5,1.6,0.2),mgp=c(1.8,0.3,0))
set.seed(2)
n <- 50
x <- 100+rnorm(n)
y <- 100+rnorm(n)
z <- 100+rnorm(n)*10
plot(x,y,pch=21,bg='white',axes=FALSE,xlab='',ylab='')
mtext(bquote('a) r'^2*'='*.(signif(cor(x,y),2))),line=0,adj=0,cex=0.9)
axis(1)
axis(2)
mtext('y',side=1,line=1.3,cex=1.0)
mtext('x',side=2,line=1.5,cex=1.0)
plot(x,z,pch=21,bg='white',axes=FALSE,xlab='',ylab='')
mtext(bquote('b) r'^2*'='*.(signif(cor(x,z),2))),line=0,adj=0,cex=0.9)
axis(1)
axis(2)
mtext('y',side=1,line=1.3,cex=1.0)
mtext('z',side=2,line=1.4,cex=1.0)
plot(y,z,pch=21,bg='white',axes=FALSE,xlab='',ylab='')
mtext(bquote('c) r'^2*'='*.(signif(cor(y,z),2))),line=0,adj=0,cex=0.9)
axis(1)
axis(2)
mtext('x',side=1,line=1.3,cex=1.0)
mtext('z',side=2,line=1.4,cex=1.0)
plot(x/z,y/z,pch=21,bg='white',axes=FALSE,xlab='',ylab='')
mtext(bquote('d) r'^2*'='*.(signif(cor(x/z,y/z),2))),line=0,adj=0,cex=0.9)
axis(1)
axis(2)
mtext('y/z',side=1,line=1.3,cex=1.0)
mtext('x/z',side=2,line=1.4,cex=1.0)
dev.off()

cairo(file='errorcorrelation.pdf',width=7,height=2)
ReOs <- read.data('Kendall2006.csv',method='Re-Os',ierr=1)
KCa <- examples$KCa
UPb <- examples$UPb
pars(mfrow=c(1,3),mar=c(2.5,2.75,0.5,0.5))
isochron(ReOs)
isochron(KCa)
concordia(UPb,hide=10)
dev.off()

cairo(file='inverrorcorrelation.pdf',width=7,height=2)
pars(mfrow=c(1,3),mar=c(2.5,2.75,0.5,0.5))
isochron(ReOs,inverse=TRUE)
isochron(KCa,inverse=TRUE)
concordia(UPb,hide=10,type=2,tlim=c(247,255),ylim=c(0.0508,0.0518))
dev.off()

dat <- cbind(c(10,20,28),c(1,1,1),c(20,30,42),c(1,1,1),c(0.9,0.9,-0.9))
colnames(dat) <- c('X','sX','Y','sY','rXY')
fit <- york(dat)
xyfitted <- IsoplotR:::get.york.xy(dat,a=fit$a[1],b=fit$b[1])

cairo(file='wtdmean.pdf',width=5,height=3)
pars()
X <- c(1005, 995, 1000, 1125)
sX <- c(10, 10, 10, 100)
avg <- mean(X)
wtdavg <- sum(X/sX^2)/sum(1/sX^2)
mswd <- sum(((X-wtdavg)^2)/(sX^2))
ns <- length(X)
plot(1:ns,X,xlim=c(1/2,ns+1/2),ylim=c(900,1335),
     xlab='i',ylab='t (Ma)',bty='n',xaxt='n',pch=16)
axis(side=1,at=1:ns,labels=1:ns)
arrows(1:ns,X-2*sX,1:ns,X+2*sX,angle=90,code=3,length=0.05)
lines(c(1/2,ns+1/2),rep(avg,2),lty=2)
lines(c(1/2,ns+1/2),rep(wtdavg,2))
dev.off()

cairo(file='chi2wtdmean.pdf',width=6,height=3)
pars(mar=c(3,3,1,1))
df <- length(X)-1
chi2stat <- mswd*df
x <- seq(from=0,to=20,length.out=100)
y <- dchisq(x,df=df)
plot(x,y,type='l',xlab=expression(chi^2),
     ylab=expression('f('*chi^2*')'),
     bty='n',ylim=c(0,0.3))
lines(x=rep(chi2stat,2),y=c(0,1),lty=2)
x <- seq(from=chi2stat,to=20,length.out=50)
y <- dchisq(x,df=df)
polygon(c(x,tail(x,1),x[1]),c(y,0,0),col='black')
text(x=chi2stat,y=0.2,pos=4,xpd=NA,
     labels=substitute(chi[stat]^2*'='*a,
                       list(a=signif(chi2stat,3)) ) )
p <- signif(pchisq(chi2stat,df=df,lower.tail=FALSE),3)
text(x=10,y=0.03,pos=4,xpd=NA,
     labels=substitute(p*'='*a,list(a=p)))
dev.off()

cairo(file='chi2.pdf',width=6,height=3)
pars(mar=c(3,3,1,1))
x <- seq(from=0.001,to=10,length.out=50)
y <- dchisq(x,df=1)
plot(x,y,type='l',xlab=expression(chi^2),
     ylab=expression('f('*chi^2*')'),
     bty='n',ylim=c(0,0.5))
lines(x=rep(fit$mswd,2),y=c(0,1),lty=2)
x <- seq(from=fit$mswd,to=10,length.out=50)
y <- dchisq(x,df=1)
polygon(c(x,tail(x,1),x[1]),c(y,0,0),col='black')
text(x=fit$mswd,y=0.5,pos=4,xpd=NA,
     labels=substitute(chi[stat]^2*'='*a,
                       list(a=signif(fit$mswd,3)) ) )
text(x=(10+fit$mswd)/3,y=y[1],
     pos=4,xpd=NA,
     labels=substitute(p*'='*a,list(a=signif(fit$p.value,2)) ) )
dev.off()

cairo(file='mswd.pdf',width=7,height=3.5)
pars(mar=c(3,3,1,1))
M <- 2.5
x <- seq(from=0,to=M,length.out=100)
y <- dchisq(x,df=1)
i <- which.max(y)
dof <- c(2,5,10,20,50)
plot(x,y,type='l',xlab='MSWD',ylab='f(MSWD)',bty='n',ylim=c(0,2.5))
text(x=0,y=M,labels='df=1',xpd=NA,pos=3)
xadj <- c(1,3,0,1,1)/2
yadj <- c(-1,-1,-1,-1,-1)/2
for (i in 1:length(dof)){
    df <- dof[i]
    x <- seq(from=0,to=M*df,length.out=100)
    y <- dchisq(x,df=df)
    lines(x=x/df,y=y*df)
    j <- which.max(y)
    text(x[j]/df,y[j]*df,labels=paste0('df=',df),adj=c(xadj[i],yadj[i]))
}
dev.off()

cairo(file='~/Documents/temp/isochronMSWD.pdf',height=6,width=2.5)
pars(mfrow=c(3,1),mar=c(3,3,3,0.5))
set.seed(36)
tt <- 1000
ns <- 10
x <- runif(ns,min=0,max=0.3)
SrRb <- IsoplotR:::get.RbSr.ratio(tt,st=0)[1]
y <- 0.7 + x*SrRb
sx <- 0.004
sy <- 0.0001
rxy <- 0.3
sxy <- rxy*sx*sy
E <- rbind(c(sx^2,sxy),c(sxy,sy^2))
dxy <- MASS::mvrnorm(ns,mu=c(0,0),Sigma=E)
X <- x + dxy[,1]
Y <- y + dxy[,2]
ydat <- cbind(X,sx,Y,sy,rxy)
colnames(ydat) <- c('Rb87/Sr86', 'err[Rb87/Sr86]',
                    'Sr87/Sr86', 'err[Sr87/Sr86]', 'rho')
RbSr <- list(x=ydat,format=1)
class(RbSr) <- c('RbSr','PD')
isochron(RbSr)
set.seed(3)
dxy <- MASS::mvrnorm(ns,mu=c(0,0),Sigma=E)
X <- x + dxy[,1]/10
Y <- y + dxy[,2]/10
ydat <- cbind(X,sx,Y,sy,rxy)
RbSr <- list(x=ydat,format=1)
class(RbSr) <- c('RbSr','PD')
isochron(RbSr)
set.seed(3)
dxy <- MASS::mvrnorm(ns,mu=c(0,0),Sigma=E)
X <- x + dxy[,1]*5
Y <- y + dxy[,2]*5
ydat <- cbind(X,sx,Y,sy,rxy)
RbSr <- list(x=ydat,format=1)
class(RbSr) <- c('RbSr','PD')
isochron(RbSr,model=1)
dev.off()

cairo(file='~/Documents/temp/model-1.pdf',height=3,width=6)
pars(mfrow=c(1,2),mar=c(3,3,3,0.5))
d <- c(99.37336,97.07454,101.75070,101.43449,97.02336,
       102.97518,101.93286,102.31634,99.57383,101.81529)
de <- c(8.46319657,1.86649249,-0.07321272,12.75171960,-12.72951972,
        -21.69428369,1.09848605,8.01051989,4.08519189,3.64344448)
dd <- d + de
d3 <- cbind(dd,2)
weightedmean(d3,random.effects=FALSE,detect.outliers=FALSE)
isochron(RbSr,model=1)
dev.off()

dat <- read.data('Ludwig.csv',method='other')
cairo(file='~/Documents/temp/frequency1.pdf',height=3,width=8)
pars(mfrow=c(1,2))
cad(dat[,1])
kde(dat[,1])
dev.off()

cairo(file='~/Documents/temp/frequency2.pdf',height=3,width=4)
pars()
radialplot(dat)
dev.off()

fn <- system.file('UPb8.csv',package='IsoplotR')
cairo(file='~/Documents/temp/logtrans1.pdf',height=3,width=8)
UPb <- read.data(fn,method='U-Pb',format=8)
pars(mfrow=c(1,2))
kde(UPb,type=5,from=-50,to=200,bw=20)
kde(UPb,log=TRUE,type=5,from=1,to=500,bw=0.5)
dev.off()

cairo(file='~/Documents/temp/logtrans2.pdf',height=3,width=4)
pars()
radialplot(UPb,type=5)
dev.off()

cairo(file='~/Documents/temp/3means.pdf',height=3,width=9)
set.seed(3)
pars(mfcol=c(2,3))
mu <- 1000
sigma <- 5
err <- 5
m <- mu-3*(sigma+err)
M <- mu+3*(sigma+err)
x <- seq(from=m,to=M,length.out=100)
ns <- 20
# column 1
plot(x=rep(mu,2),y=c(0,0.1),type='l',xlim=c(m,M),yaxt='n',bty='n',col='blue')
lines(x=c(m,M),y=rep(0,2))
legend('topleft',legend='1a',bty='n')
lines(x=x,y=dnorm(x,mean=mu,sd=err),col='red')
rx <- rnorm(ns,mean=mu,sd=sigma)
sx <- rep(err,ns)
dat <- cbind(x=rx,err=sx)
weightedmean(dat,random.effects=FALSE)
legend('topleft',legend='1b',bty='n')
# column 2
plot(x=x,y=dnorm(x,mean=mu,sd=sigma),type='l',
     xlim=c(m,M),yaxt='n',bty='n',col='blue')
legend('topleft',legend='2a',bty='n')
lines(x=x,y=dnorm(x,mean=mu,sd=sqrt(sigma^2+err^2)),col='red')
rx <- rnorm(ns,mean=mu,sd=sqrt(sigma^2+err^2))
sx <- rep(err,ns)
dat <- cbind(x=rx,err=sx)
weightedmean(dat,random.effects=TRUE)
legend('topleft',legend='2b',bty='n')
# column 3
lmu <- log(10)
lsigma <- 5/10
lerr <- 5/10
m <- lmu-3*(lsigma+lerr)
M <- lmu+1.5*(lsigma+lerr)
lx <- seq(from=m,to=M,length.out=100)
d1 <- dnorm(lx,mean=lmu,sd=lsigma)
d2 <- dnorm(lx,mean=lmu,sd=sqrt(lsigma^2+lerr^2))
x <- exp(lx)
dxin <- diff(lx)
dxout <- diff(x)
y1 <- c(dxin,utils::tail(dxin,1))*d1/c(dxout,utils::tail(dxout,1))
y2 <- c(dxin,utils::tail(dxin,1))*d2/c(dxout,utils::tail(dxout,1))
plot(x,y1,type='l',bty='n',col='blue')
lines(x,y2,col='red')
legend('topleft',legend='3a',bty='n')
rx <- exp(rnorm(ns,mean=lmu,sd=sqrt(lsigma^2+lerr^2)))
sx <- rx*rnorm(ns,mean=lerr,sd=lerr/3)
radialplot(cbind(rx,sx),transformation='log')
legend('topleft',legend='3b',bty='n')
if (FALSE){
    # column 4
    M <- lmu+3*(lsigma+lerr)
    lx <- seq(from=m,to=M,length.out=100)
    d1 <- dnorm(lx,mean=lmu,sd=lsigma)
    d2 <- dnorm(lx,mean=lmu,sd=sqrt(lsigma^2+lerr^2))
    x <- exp(lx)
    plot(x,d1,type='l',bty='n',log='x',col='blue')
    lines(x,d2,col='red')
    legend('topleft',legend='3c',bty='n')
    rx <- rnorm(ns,mean=lmu,sd=sqrt(lsigma^2+lerr^2))
    sx <- rnorm(ns,mean=lerr,sd=lerr/3)
    radialplot(cbind(rx,sx),transformation='log')
    legend('topleft',legend='3d',bty='n')
}
dev.off()

cairo(file='~/Documents/temp/3meansb.pdf',height=3,width=6)
pars()
radialplot(cbind(rx,sx),transformation='log')
legend('topleft',legend='3d',bty='n')
dev.off()

cairo(file='~/Documents/temp/concordia.pdf',height=2.5,width=7.5)
pars(mfrow=c(1,3))
fn <- system.file('UPb8.csv',package='IsoplotR')
UPb <- read.data(fn,method='U-Pb',format=8)
concordia(UPb,type=1)
concordia(UPb,type=2,ticks=c(25,50,100,500,2000,3000,4000))
concordia(UPb,type=3)
dev.off()

cairo(file='~/Documents/temp/concordia_age.pdf',height=2,width=8)
pars(mfrow=c(1,4))
conc <- concordia(examples$UPb,hide=10,show.age=1)
xconc <- IsoplotR:::age_to_Pb207U235_ratio(conc$age[1])[1]
yconc <- IsoplotR:::age_to_Pb206U238_ratio(conc$age[1])[1]
concordia(examples$UPb,hide=10,show.age=1,
          xlim=c(0.28,0.2825),ylim=c(0.0396,0.0399),
          ticks=seq(from=250,to=253,by=0.5))
points(xconc,yconc,pch=21,bg='white')
conc <- concordia(examples$UPb,hide=10,show.age=1,type=2,
                  xlim=c(24.9,25.4),ylim=c(0.0508,0.0518))
xconc <- IsoplotR:::age_to_U238Pb206_ratio(conc$age[1])[1]
yconc <- IsoplotR:::age_to_Pb207Pb206_ratio(conc$age[1])[1]
concordia(examples$UPb,hide=10,show.age=1,type=2,
          xlim=c(25.1,25.2),ylim=c(0.0512,0.05145),
          ticks=seq(from=250,to=253,by=0.5))
points(xconc,yconc,pch=21,bg='white')
dev.off()

UPb1 <- UPb
UPb1$x <- matrix(rep(UPb$x[1,,drop=FALSE],9),nrow=9,byrow=TRUE)
UPb1$x[,1] <- UPb1$x[,1] + rnorm(9,mean=0,sd=0.000001)
UPb1$x[,3] <- UPb1$x[,3] + rnorm(9,mean=0,sd=0.00001)
colnames(UPb1$x) <- colnames(UPb$x)

fn <- system.file('UPb1.csv',package='IsoplotR')
UPb <- read.data(fn,method='U-Pb',format=1)
conc <- concordia(UPb,show.age=1,hide=10,exterr=FALSE)
dev.off()
cairo(file='~/Documents/temp/concordia_MSWD.pdf',height=4.5,width=8)
set.seed(11)
pars(mfrow=c(2,3),mar=c(2,3,3,0.5))
xl <- c(0.277,0.287)
yl <- c(0.0394,0.0402)
xy <- MASS::mvrnorm(n=10,mu=conc$x,Sigma=conc$cov/1000)
UPb1 <- UPb
UPb1$x[,c(1,3)] <- xy
concordia(UPb1,show.age=1,exterr=FALSE,xlim=xl,ylim=yl)
xy <- MASS::mvrnorm(n=10,mu=conc$x,Sigma=conc$cov*20)
UPb2 <- UPb
UPb2$x[,c(1,3)] <- xy
concordia(UPb2,hide=10,show.age=1,exterr=FALSE,xlim=xl,ylim=yl)
xy <- MASS::mvrnorm(n=10,mu=conc$x,Sigma=conc$cov*100)
UPb3 <- examples$UPb
UPb3$x[,c(1,3)] <- xy
concordia(UPb3,hide=10,show.age=1,exterr=FALSE,xlim=xl,ylim=yl)
xl <- c(0.280,0.284)
yl <- c(0.0396,0.0399)
UPb4 <- UPb; UPb4$x[,1] <- UPb4$x[,1]-0.000615
concordia(UPb4,hide=10,show.age=1,exterr=FALSE,xlim=xl,ylim=yl)
UPb5 <- UPb; UPb5$x[,1] <- UPb5$x[,1]-0.0004
concordia(UPb5,hide=10,show.age=1,exterr=FALSE,xlim=xl,ylim=yl)
UPb6 <- UPb; UPb6$x[,1] <- UPb6$x[,1]+0.002
concordia(UPb6,hide=10,show.age=1,exterr=FALSE,xlim=xl,ylim=yl)
dev.off()

cairo(file='~/Documents/temp/UPbisochron13.pdf',height=3,width=3)
pars()
fn <- system.file('UPb5.csv',package='IsoplotR')
dat <- read.data(fn,method='U-Pb',format=5)
concordia(dat,type=2,show.age=2,xlim=c(0,22),ylim=c(0,0.8))
dev.off()

cairo(file='~/Documents/temp/UPbisochron46.pdf',height=3,width=6)
pars()
pars(mfrow=c(1,2))
isochron(dat,type=1,xlim=c(0,22),ylim=c(0,0.05))
isochron(dat,type=2,xlim=c(0,3.5),ylim=c(0,0.07))
dev.off()

cairo(file='~/Documents/temp/UPbisochron78.pdf',height=3,width=9)
pars()
pars(mfrow=c(1,3),mar=c(3,3,1,1))
fn <- system.file('UPb8.csv',package='IsoplotR')
dat <- read.data(fn,method='U-Pb',format=8)
fit <- isochron(dat,type=3,plot=FALSE)
ydat <- data2york(dat,option=6,tt=fit$age[1])
sn <- 19
Th2Pb6 <- dat$x[sn,'Th232U238']*dat$x[sn,'U238Pb206']
Pb8o6 <- ydat[sn,'Y']
l232 <- settings('lambda','Th232')[1]
Pb86 <- Pb8o6 + Th2Pb6*(exp(l232*fit$age[1])-1)
plot(x=c(0,Th2Pb6),y=c(Pb8o6,Pb86),type='p',
     xlim=c(0,1500),ylim=c(0,2.5),pch=16,cex=1.2)
arrows(x0=Th2Pb6,y0=Pb86,x1=0,y1=Pb8o6,code=2,length=0.1)
arrows(x0=Th2Pb6,y0=Pb86-2*ydat[sn,'sY'],
       x1=Th2Pb6,y1=Pb86+2*ydat[sn,'sY'],
       code=3,length=0.1,angle=90)
arrows(x0=0,y0=Pb8o6-2*ydat[sn,'sY'],
       x1=0,y1=Pb8o6+2*ydat[sn,'sY'],
       code=3,length=0.1,angle=90)
isochron(dat,type=1,xlim=c(0,600),ylim=c(0,2.5))
isochron(dat,type=2,xlim=c(0,100),ylim=c(0,3))
dev.off()

cairo(file='~/Documents/temp/Hourigan.pdf',height=6,width=6)
Hourigan <- read.csv('Hourigan.csv',header=TRUE)
UPb <- read.data(Hourigan[,1:5],method='U-Pb',format=2)
concordia(UPb,levels=Hourigan[,'ThU'],clabel='Th/U',
          ellipse.fill=c('white','red'),show.age=2)
dev.off()

cairo(file='~/Documents/temp/commonPbisochron13.pdf',height=4,width=7)
pars(mfrow=c(1,2),mar=c(2,2.5,4,0.5))
dat <- read.data('BH14.csv',method='U-Pb',format=2,ierr=3)
good <- dat$x[,'U238Pb206']>5
gdat <- dat
gdat$x <- dat$x[good,]
concordia(gdat,show.age=2,type=2)
concordia(gdat,show.age=2,type=2,common.Pb=2,xlim=c(97,108))
dev.off()

cairo(file='~/Documents/temp/commonPbisochron46.pdf',height=4,width=7)
pars(mfrow=c(1,2),mar=c(2,2.5,4,0.5))
fn <- system.file('UPb5.csv',package='IsoplotR')
dat <- read.data(fn,method='U-Pb',format=5)
concordia(dat,show.age=2,type=2)
concordia(dat,show.age=2,type=2,common.Pb=2,xlim=c(20.8,21.4))
dev.off()

cairo(file='~/Documents/temp/commonPbisochron46nominal.pdf',height=4,width=7)
pars(mfrow=c(1,2),mar=c(2,2.5,4,0.5))
settings('iratio','Pb206Pb204',20)
settings('iratio','Pb207Pb204',16)
fn <- system.file('UPb5.csv',package='IsoplotR')
dat <- read.data(fn,method='U-Pb',format=5)
concordia(dat,show.age=2,type=2,anchor=1)
concordia(dat,show.age=2,type=2,common.Pb=2,anchor=1,xlim=c(20.8,21.4))
dev.off()

cairo(file='~/Documents/temp/commonPbisochron13detrital.pdf',height=3,width=7)
pars(mfrow=c(1,2),mar=c(2,2.5,0.5,0.5))
settings('iratio','Pb207Pb206',0.8)
dat <- read.data('BH14.csv',method='U-Pb',format=2,ierr=3)
dat$x <- dat$x[1:2,]
dat$x[1,] <- c(50,0.75,0.1,0.01,0)
dat$x[2,] <- c(20,0.75,0.6,0.01,0)
concordia(dat,type=2,show.numbers=FALSE,
          xlim=c(0,110),ylim=c(0,0.8),ticks=c(100,200,1000,4000))
points(x=0,y=0.8,pch='.')
concordia(dat,type=2,common.Pb=1,xlim=c(50,85),
          ticks=c(70,80,90,100,110,120))
dev.off()

cairo(file='~/Documents/temp/SK.pdf',height=3,width=3)
pars()
sk <- IsoplotR:::stacey.kramers(c(500,1000,3000))
Pb76 <- sk[,'i74']/sk[,'i64']
x <- c(IsoplotR:::age_to_U238Pb206_ratio(500)[1],
       IsoplotR:::age_to_U238Pb206_ratio(1000)[1],
       IsoplotR:::age_to_U238Pb206_ratio(3000)[1])
y <- c(IsoplotR:::age_to_Pb207Pb206_ratio(500)[1],
       IsoplotR:::age_to_Pb207Pb206_ratio(1000)[1],
       IsoplotR:::age_to_Pb207Pb206_ratio(3000)[1])
concordia(type=2,xlim=c(0,13),ylim=c(0,1.15),
          tlim=c(100,5500),ticks=c(500,1000,2000,3000,4000))
lines(x=c(0,x[1]),y=c(Pb76[1],y[1]))
lines(x=c(0,x[2]),y=c(Pb76[2],y[2]))
lines(x=c(0,x[3]),y=c(Pb76[3],y[3]))
dev.off()

cairo(file='~/Documents/temp/diseq.pdf',height=4,width=9)
pars(mfrow=c(1,2),mar=c(3.2,3.2,2.8,0.5))
fn <- system.file("diseq.csv",package="IsoplotR")
UPb <- read.data(fn,method='U-Pb',format=2)
concordia(UPb,type=2,show.age=1,
          xlim=c(-500,5000),ylim=c(0.045,0.056),
          ticks=c(2,3,5,10,50,100,200,300,400))
d <- diseq(U48=list(x=0,option=1),ThU=list(x=2,option=1),
           RaU=list(x=2,option=1),PaU=list(x=2,option=1))
UPb <- read.data(fn,method='U-Pb',format=2,d=d)
concordia(UPb,type=2,show.age=1,
          xlim=c(-500,5000),ylim=c(0.047,0.056),
          ticks=c(2,3,5,10,50,100,200,300,400))
dev.off()

cairo(file='~/Documents/temp/PbPb.pdf',height=4,width=8)
pars(mfrow=c(1,2),mar=c(2.5,2.5,3,0))
fn <- system.file("PbPb1.csv",package="IsoplotR")
PbPb <- read.data(fn,method='Pb-Pb',format=1,ierr=1)
PbPb$x[,5] <- 0.95
isochron(PbPb,inverse=FALSE,show.numbers=TRUE)
isochron(PbPb,inverse=TRUE,show.numbers=TRUE)
dev.off()

cairo(file='~/Documents/temp/Kamber.pdf',height=4,width=7)
pars(mfrow=c(1,2),mar=c(2.5,2.5,4,0))
PbPb <- read.data('Kamber.csv',method='Pb-Pb',format=1,ierr=2)
isochron(PbPb,model=2,inverse=FALSE,
         growth=TRUE,xlim=c(10.9,15),sigdig=4)
isochron(PbPb,model=2,inverse=TRUE,growth=TRUE,sigdig=4,
         xlim=c(0.065,0.09),ylim=c(0.95,1.2))
dev.off()

cairo(file='~/Documents/temp/ThPb.pdf',height=4,width=7)
pars(mfrow=c(1,2),mar=c(2.5,2.5,4,0))
isochron(examples$ThPb,model=1,inverse=FALSE)
isochron(examples$ThPb,model=1,inverse=TRUE)
dev.off()

cairo(file='~/Documents/temp/ThPbSingleGrain.pdf',height=4,width=5)
pars(mfrow=c(2,2))
set.seed(3)
ns <- 10
buffer <- 100
x0 <- 1000
x <- runif(ns,min=buffer,max=x0-buffer)
a <- 1/settings('iratio','Pb208Pb204')[1]
b <- -a/x0
y <- a + b * x + rnorm(ns,0,0.003)
plot(x,y,type='n',xlim=c(0,x0+buffer),ylim=c(0,a),
     xlab=expression(''^232*'Th/'^208*'Pb'),
     ylab=expression(''^204*'Pb/'^208*'Pb'),
     bty='n')
lines(x=c(0,x0),y=c(a,0),lwd=2)
X <- rep(0,ns)
for (i in 1:ns){
    X[i] <- x[i] - y[i]/b
    lines(x=c(x[i],X[i]),y=c(y[i],0),col='grey50')
}
points(x,y,pch=21,bg='white')
legend('topright',legend='a)',bty='n')
a <- 1/a
b <- 1/x0
x <- x/y
y <- 1/y
plot(x,y,type='n',bty='n',
     xlim=c(0,max(x)),ylim=c(0,max(y)),
     xlab=expression(''^208*'Pb/'^204*'Pb'),
     ylab=expression(''^232*'Th/'^204*'Pb'))
xmax <- max(x)
b <- rep(0,ns)
for (i in 1:ns){
    b[i] <- (y[i]-a)/x[i]
    yi <- a + b[i]*xmax
    lines(x=c(0,xmax),y=c(a,yi),col='grey50')
}
points(x,y,pch=21,bg='white')
legend('topleft',legend='b)',bty='n')
kde(IsoplotR:::get.ThPb.age(1/X,0)[,1],log=TRUE,binwidth=0.15)
legend('topright',legend='c)',bty='n')
kde(IsoplotR:::get.ThPb.age(b,0)[,1],log=TRUE,binwidth=0.15)
legend('topleft',legend='d)',bty='n')
dev.off()

cairo(file='~/Documents/temp/ArAr.pdf',height=6,width=6)
pars(mfrow=c(2,2),mar=c(3,3,3,0))
isochron(examples$ArAr,inverse=FALSE)
isochron(examples$ArAr,inverse=TRUE)
radialplot(examples$ArAr,i2i=TRUE)
radialplot(examples$ArAr,i2i=FALSE)
dev.off()

cairo(file='~/Documents/temp/agespectrum.pdf',height=3.5,width=7)
pars(mfrow=c(1,2),mar=c(3.5,3,3,0.5))
weightedmean(examples$ArAr,random.effects=FALSE)
agespectrum(examples$ArAr,random.effects=FALSE)
dev.off()

cairo(file='~/Documents/temp/RbSr.pdf',height=2.5,width=9)
pars(mfrow=c(1,4),mar=c(3.5,3,3,0.5))
isochron(examples$RbSr,inverse=FALSE)
legend('topleft','a)',bty='n',cex=1.2)
isochron(examples$RbSr,inverse=TRUE)
legend('topright','b)',bty='n',cex=1.2)
kde(examples$RbSr,i2i=TRUE)
legend('topleft','c)',bty='n',cex=1.2)
cad(examples$RbSr,i2i=TRUE)
legend('topleft','d)',bty='n',cex=1.2)
dev.off()

cairo(file='~/Documents/temp/UThHe.pdf',height=3,width=9)
fn <- system.file('UThHe.csv',package='IsoplotR')
UThHe <- read.data(fn,method='U-Th-He')
pars(mfrow=c(1,3),mar=c(3.5,3,3,0.5))
isochron(UThHe,model=3)
helioplot(UThHe,logratio=FALSE)
helioplot(UThHe,logratio=TRUE)
dev.off()

if (FALSE){
    He <- c(0.21,0.18,0.69,0.38,0.77,0.50,0.72,0.99,0.38,0.780)
    U <- c(27,37,57,91,20,90,94,66,63,62)
    l38 <- settings('lambda','U238')[1]
    l35 <- settings('lambda','U235')[1]
    U85 <- settings('iratio','U238U235')[1]
    P <- (8*U85*l38 + 7*l35)*U/(U85+1)
    signif(mean(He/P),3)
    signif(1/mean(P/He),3)
}

cairo(file='~/Documents/temp/FTradial.pdf',height=4.5,width=4.5)
pars(mfrow=c(2,2),mar=c(2.5,2.5,3,0.5))
radialplot(examples$FT1,transformation='log')
legend('topleft',legend='a)',bty='n')
radialplot(examples$FT1,transformation='arcsin')
legend('topleft',legend='b)',bty='n')
radialplot(examples$FT2,transformation='log')
legend('topleft',legend='c)',bty='n')
radialplot(examples$FT2,transformation='sqrt')
legend('topleft',legend='d)',bty='n')
dev.off()
