setwd('~/Documents/geotopes/R/')

blankcorr <- function(dat,bi,si){
    return(dat[si] - mean(dat[bi]))
}

meanUPbRatio <- function(D,P,bi,si){
    ratiovec <- blankcorr(D,bi,si)/blankcorr(P,bi,si)
    return(mean(ratiovec))
}

calibration <- function(D,P,tstd,lambda,bi,si){
    atomic <- exp(lambda*tstd)-1
    measured <- meanUPbRatio(D,P,bi,si)
    return(atomic/measured)
}

atomic <- function(samp,stand,num='Pb206',den='U238',tstd,lambda,tb,ts){
    bi <- samp$Time<tb
    si <- samp$Time>ts
    measured <- meanUPbRatio(D=samp[[num]],P=samp[[den]],bi=bi,si=si)
    cal <- calibration(D=std[[num]],P=std[[den]],
                       tstd=tstd,lambda=lambda,bi=bi,si=si)
    return(cal*measured)
}

tUPb <- function(samp,stand,num='Pb206',den='U238',tstd,lambda,tb,ts){
    a <- atomic(samp,stand,num,den,tstd,lambda,tb,ts)
    return(log(1+a)/lambda)
}

tb <- 20 # end of blank (in seconds)
ts <- 25 # start of signal (in seconds)
tstd <- 600.4 # Ma
l38 <- 1.55125e-4 # Myr-1
l35 <- 9.8485e-4 # Myr-1
l32 <- 0.495e-4 # Myr-1

smp <- read.csv("91500.csv",header=TRUE)
std <- read.csv("GJ1.csv",header=TRUE)
t0638 <- tUPb(samp=smp,stand=std,num='Pb206',den='U238',
              tstd=tstd,lambda=l38,tb=tb,ts=ts)
t0735 <- tUPb(samp=smp,stand=std,num='Pb207',den='U235',
              tstd=tstd,lambda=l35,tb=tb,ts=ts)
t0832 <- tUPb(samp=smp,stand=std,num='Pb208',den='Th232',
              tstd=tstd,lambda=l32,tb=tb,ts=ts)
