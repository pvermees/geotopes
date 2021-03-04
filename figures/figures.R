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
