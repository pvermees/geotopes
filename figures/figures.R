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
