library(IsoplotR)

setwd('/home/pvermees/Documents/geotopes/figures/')

pars <- function(mar=c(2.5,2.3,0.5,0.25),mgp=c(1.5,0.5,0),mfrow=c(1,1)){
    par(list(mar=mar,mgp=mgp,mfrow=mfrow))
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
fit <- IsoplotR::york(dat)
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
